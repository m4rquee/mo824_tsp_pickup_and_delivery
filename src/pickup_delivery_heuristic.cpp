#include "pickup_delivery_utils.hpp"
#include <cstdlib>
#include <iostream>
#include <lemon/list_graph.h>
#include <lemon/preflow.h>
#include <string>

bool solve(Pickup_Delivery_Instance &P, double &LB, double &UB, DNodeVector &Sol) {
    // Generates the arborescence that will guide the route creation:
    MinCostArb arb_solver(P.g, P.weight);
    arb_solver.run(P.source);// root the arborescence in the source
    // As a spanning digraph rooted at the source this is itself a LB:
    LB = max(LB, arb_solver.arborescenceCost());
    cerr << LB;

    bool improved = arborescence_heuristic(P, LB, UB, Sol, arb_solver);
    cerr << ' ' << UB;

    if (LB != UB)// if can improve
        improved |= local_search(P, LB, UB, Sol);
    cerr << ' ' << UB;

    return improved;
}

int main(int argc, char *argv[]) {
    Digraph g;// graph declaration
    string digraph_filename, source_node_name, target_node_name;
    DNodeStringMap vname(g); // name of graph nodes
    DNodePosMap px(g), py(g);// xy-coodinates for each node
    DNodeColorMap vcolor(g); // color of nodes
    ArcStringMap aname(g);   // name for graph arcs
    ArcColorMap ecolor(g);   // color of edges
    ArcValueMap lpvar(g);    // used to obtain the contents of the LP variables
    ArcValueMap weight(g);   // edge weights
    vector<DNode> V;
    Digraph::NodeMap<DNode> del_pickup(g);// map a delivery to it's pickup
    DNodeBoolMap is_pickup(g, false);     // used to quickly check if a node is a pickup

    set_pdfreader("xdg-open");// the Linux will choose the default one
    if (argc < 3) {
        cout << endl
             << "Projeto de MO824: Rota com coleta e entrega de peso minimo;"
             << endl
             << "Usage: " << argv[0]
             << "  <pickup_delivery_digraph_filename> <epsilon>" << endl
             << endl;
        cout << "Example:" << endl
             << "\t" << argv[0] << " "
             << getpath(argv[0]) + "../instances/pickup_delivery_5.dig 2" << endl
             << endl
             << "\t" << argv[0] << " "
             << getpath(argv[0]) + "../instances/pickup_delivery_10.dig 0.5" << endl
             << endl;
        exit(0);
    }

    digraph_filename = argv[1];
    double LB = 0, UB = MY_INF;// considere MY_INF como infinito.
    if (argc >= 4)
        LB = atof(argv[3]);
    if (argc >= 5)
        UB = atof(argv[4]);
    DNodeVector pickup, delivery;
    DNode source, target;
    int npairs;

    if (!ReadPickupDeliveryDigraph(digraph_filename, g, vname, px, py, weight,
                                   source, target, npairs, pickup, delivery,
                                   del_pickup, is_pickup)) {
        cout << "Erro na leitura do grafo de entrada." << endl;
        exit(EXIT_FAILURE);
    }

    Pickup_Delivery_Instance P(g, vname, px, py, weight, source, target, npairs,
                               pickup, delivery, del_pickup, is_pickup, INFINITY);
    PrintInstanceInfo(P);
    P.eps_min = P.eps_max = atof(argv[2]);

    DNodeVector Solucao;

    bool melhorou = solve(P, LB, UB, Solucao);

    if (melhorou) {
        PrintSolution(P, Solucao, "Solucao obtida.");
        cout << "custo: " << UB << endl;
    }
    return 0;
}
