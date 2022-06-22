#include "pickup_delivery_utils.hpp"

Pickup_Delivery_Instance::Pickup_Delivery_Instance(
    Digraph &graph, DNodeStringMap &vvname, DNodePosMap &posx,
    DNodePosMap &posy, ArcValueMap &eweight, DNode &sourcenode,
    DNode &targetnode, int &vnpairs, DNodeVector &vpickup,
    DNodeVector &vdelivery, Digraph::NodeMap<DNode> &del_pickup,
    DNodeBoolMap &is_pickup, int &time_limit)
    : g(graph), vname(vvname), px(posx), py(posy), weight(eweight),
      nnodes(2 * vnpairs + 2), source(sourcenode), target(targetnode),
      npairs(vnpairs), pickup(vpickup), delivery(vdelivery),
      del_pickup(del_pickup), is_pickup(is_pickup), time_limit(time_limit) {
  // Store the out arcs of each node sorted by weight:
  ArcCmp arcCmp(weight);           // arc comparator based on this weight map
  min_arc_heap sorting_heap(arcCmp); // aux heap used for sorting
  // Sort the nodes out arcs for each node:
  for (DNodeIt n(g); n != INVALID; ++n) {
    ordered_arcs[n].reserve(npairs);
    for (OutArcIt a(g, n); a != INVALID; ++a) {
      sorting_heap.push(a);                   // add all out arcs to the heap
      weight_map[n][g.target(a)] = weight[a]; // add all weights to a map
    }
    while (!sorting_heap.empty()) { // sort the arcs by popping the heap
      ordered_arcs[n].push_back(sorting_heap.top());
      sorting_heap.pop();
    }
  }
}

void Pickup_Delivery_Instance::start_counter() {
  start = chrono::system_clock::now();
}

void PrintInstanceInfo(Pickup_Delivery_Instance &P) {
  cout << endl << endl;
  cout << "Pickup Delivery Graph Informations" << endl;
  cout << "\tTime limit = " << P.time_limit << "s" << endl;
  cout << "\tSource = " << P.vname[P.source] << endl;
  cout << "\tTarget = " << P.vname[P.target] << endl;
  for (int i = 0; i < P.npairs; i++) {
    cout << "\tPair pickup-->delivery: " << P.vname[P.pickup[i]] << " --> "
         << P.vname[P.delivery[i]] << endl;
  }
  cout << endl;
}

void PrintSolution(Pickup_Delivery_Instance &P, DNodeVector &Sol,
                   const string &msg) {
  // Imprime a solucao no terminal.
  cout << msg << endl << "\t";
  cout << P.vname[Sol[0]];
  for (int i = 1; i < P.nnodes; i++)
    cout << "-->" << P.vname[Sol[i]];
  cout << endl;
}

void graph_pruning(Digraph &g, const DNode &source, const DNode &target,
                   const int &npairs, const DNodeVector &pickup,
                   const DNodeVector &delivery) {
  cout << "Limpa arcos inválidos:" << endl;
  cout << "Arcos antes da limpeza: " << countArcs(g) << endl;
  // pickup[i] -> target (the last one before the target must be a delivery)
  // pickup[i] -> source (cannot go back to the source)
  for (int i = 0; i < npairs; i++)
    for (OutArcIt a(g, pickup[i]); a != INVALID; ++a)
      if (g.target(a) == target or g.target(a) == source)
        g.erase(a);
  // delivery[i] -> pickup[i] (a pickup is always before the corresponding delivery)
  // delivery[i] -> source (cannot go back to the source)
  // source -> delivery[i] (the second one after the source must be a pickup)
  for (int i = 0; i < npairs; i++) {
    for (OutArcIt a(g, delivery[i]); a != INVALID; ++a)
      if (g.target(a) == pickup[i] or g.target(a) == source)
        g.erase(a);
    for (InArcIt a(g, delivery[i]); a != INVALID; ++a)
      if (g.source(a) == source) {
        g.erase(a);
        break;
      }
  }
  // target -> v, for all nodes (target is the last node)
  for (OutArcIt a(g, target); a != INVALID; ++a)
    g.erase(a);
  // source -> target (the pickup/deliveries must be in the middle)
  for (OutArcIt a(g, source); a != INVALID; ++a)
    if (g.target(a) == target) {
      g.erase(a);
      break;
    }
  cout << "Arcos depois da limpeza: " << countArcs(g) << endl;
}

bool ReadPickupDeliveryDigraph(const string &filename, Digraph &g,
                               DNodeStringMap &vname, DNodePosMap &posx,
                               DNodePosMap &posy, ArcValueMap &weight,
                               DNode &source, DNode &target, int &npairs,
                               DNodeVector &pickup, DNodeVector &delivery,
                               Digraph::NodeMap<DNode> &del_pickup,
                               DNodeBoolMap &is_pickup) {
  ReadDigraph(filename, g, vname, posx, posy, weight);
  int n = countNodes(g);
  DNode DN[n];
  if ((n < 4) || (n % 2)) {
    cout << "Numero de vertices " << n
         << " no grafo nao eh par ou eh menor que 4." << endl;
    return false;
  }
  npairs = (n - 2) / 2;
  pickup.resize(npairs);
  delivery.resize(npairs);
  int i = 0;
  for (DNodeIt v(g); v != INVALID; ++v) {
    DN[i] = v;
    i++;
  }

  source = DN[0];
  target = DN[1];
  for (i = 0; i < npairs; i++) {
    pickup[i] = DN[2 * i + 2];
    delivery[i] = DN[2 * i + 3];
    is_pickup[pickup[i]] = true;
    del_pickup[delivery[i]] = pickup[i]; // map a delivery to it's pickup
  }

  // Remove invalid arcs (7n + 2):
  graph_pruning(g, source, target, npairs, pickup, delivery);
  return true;
}

double route_cost(Pickup_Delivery_Instance &P, const DNodeVector &Sol) {
  double cost = 0.0;
  for (int i = 1; i < P.nnodes; i++)
    if (P.weight_map[Sol[i - 1]].count(Sol[i]) == 0)
      return MY_INF; // the route is invalid
    else
      cost += P.weight_map[Sol[i - 1]][Sol[i]];
  return cost;
}

bool ViewPickupDeliverySolution(Pickup_Delivery_Instance &P, double &LB,
                                double &UB, DNodeVector &Sol,
                                const string &msg) {
  DigraphAttributes GA(P.g, P.vname, P.px, P.py);
  GA.SetDefaultDNodeAttrib(
      "color=LightGray style=filled width=0.2 height=0.2 fixedsize=true");
  for (ArcIt a(P.g); a != INVALID; ++a)
    GA.SetColor(a, "Invis");
  GA.SetColor(P.source, "Red"); // source and target are painted in White
  GA.SetColor(P.target, "Red");
  GA.SetShape(P.source, "star");
  GA.SetShape(P.target, "star");

  for (int i = 0; i < P.npairs; i++) { // distinguish the pickup
    GA.SetShape(P.pickup[i], "box");
  }

  if (P.npairs <= 16) { // se tiver poucos pares, dah para pintar os pares de mesma cor.
    for (int i = 0; i < P.npairs; i++) { // pinta
      GA.SetColor(P.pickup[i], ith_VisualDistinctColorName(i));
      GA.SetColor(P.delivery[i], ith_VisualDistinctColorName(i));
    }
  }
  for (int i = 1; i < P.nnodes; i++) {
    // pinta o arco Sol[i-1] -->  Sol[i]
    for (OutArcIt a(P.g, Sol[i - 1]); a != INVALID; ++a)
      if (P.g.target(a) == Sol[i]) {
        GA.SetColor(a, "Black");
        break;
      }
  }
  GA.SetLabel("Path from node " + P.vname[P.source] + " to node " +
              P.vname[P.target] + " of value " + DoubleToString(UB) +
              ". LB = " + DoubleToString(LB) + ". " + msg);
  GA.View();
  return true;
}

bool can_swap(Pickup_Delivery_Instance &P, DNodeVector &Sol, int i, int j) {
  bool i_is_pickup = P.is_pickup[Sol[i]];
  bool j_is_pickup = P.is_pickup[Sol[j]];
  if (i_is_pickup and j_is_pickup) {
    for (int k = i + 1; k < j; k++)
      if (!P.is_pickup[Sol[k]] and Sol[i] == P.del_pickup[Sol[k]])
        return false; // i is delivered at k
  } else if (i_is_pickup and !j_is_pickup) {
    if (i == 1) return false; // must start with a pickup
    if (j == P.nnodes - 2) return false; // must end with a delivery
    if (Sol[i] == P.del_pickup[Sol[j]]) return false; // i is delivered at j
    for (int k = i + 1; k < j; k++)
      if ((!P.is_pickup[Sol[k]] and Sol[i] == P.del_pickup[Sol[k]]) or
          Sol[k] == P.del_pickup[Sol[j]])
        return false; // i is delivered at k or k is delivered at j
  } else if (!i_is_pickup and j_is_pickup)
    return true; // no problem swapping
  else // (!i_is_pickup and !j_is_pickup):
    for (int k = i + 1; k < j; k++)
      if (Sol[k] == P.del_pickup[Sol[j]])
        return false; // k is delivered at j
  return true; // can swap the nodes
}

inline double swap_update(Pickup_Delivery_Instance &P, DNodeVector &Sol, int i,
                          int j) {
  double removed_arcs_weight =
      P.weight_map[Sol[i - 1]][Sol[i]] + P.weight_map[Sol[i]][Sol[i + 1]] +
      P.weight_map[Sol[j - 1]][Sol[j]] + P.weight_map[Sol[j]][Sol[j + 1]];
  swap(Sol[i], Sol[j]);
  double added_arcs_weight =
      P.weight_map[Sol[i - 1]][Sol[i]] + P.weight_map[Sol[i]][Sol[i + 1]] +
      P.weight_map[Sol[j - 1]][Sol[j]] + P.weight_map[Sol[j]][Sol[j + 1]];
  return added_arcs_weight - removed_arcs_weight;
}

bool _local_search(Pickup_Delivery_Instance &P, double &LB, double &UB,
                   DNodeVector &Sol) {
  bool improved = false;
  double best_cost, curr_cost;
  best_cost = curr_cost = route_cost(P, Sol);
  int n = P.nnodes - 2;
  // Search from end to begin because the heavier arcs are there:
  for (int i = n; i >= 1; i--)
    for (int j = i - 1; j >= 1; j--)
      if (can_swap(P, Sol, j, i)) { // the first index must be the smaller one
        curr_cost += swap_update(P, Sol, i, j); // evaluate and swap the nodes
        if (curr_cost < UB) { // found a better UB solution
          improved = true;
          UB = best_cost = curr_cost;
          NEW_UB_MESSAGE(Sol);
        } else if (curr_cost < best_cost) // just improved the solution
          best_cost = curr_cost;
        else { // reset
          swap(Sol[i], Sol[j]);
          curr_cost = best_cost;
        }
      }
  return improved;
}

bool local_search(Pickup_Delivery_Instance &P, double &LB, double &UB,
                  DNodeVector &Sol) {
  bool improved = false, aux;
  cout << "-----> Fazendo uma busca local." << endl;
  while ((aux = _local_search(P, LB, UB, Sol))) {
    improved |= aux;
    int elapsed = ELAPSED;
    if (elapsed >= P.time_limit) {
      cout << endl << "Tempo máximo de " << P.time_limit << "s atingido." << endl;
      break;
    }
  }
  cout << endl;
  return improved;
}

void _arborescence_transversal(Pickup_Delivery_Instance &P, MinCostArb &solver,
                               DNodeVector &Sol, DNode &currNode,
                               DNodeBoolMap &visited,
                               map<DNode, bool> &p_visited, int pos) {
  Sol[pos] = currNode;             // save the current node to the route
  if (pos == P.nnodes - 2) return; // if is in the last node before the target

  // Mark as visited:
  visited[currNode] = true;
  if (P.is_pickup[currNode]) // mark the currNode if it's a pickup
    p_visited[currNode] = true;

  // Min arc over all neighbours:
  DNode min_next;
  double min_cost = MY_INF;
  // Min arc over all neighbours, restricted to the arborescence arcs:
  DNode min_next_arb;
  double min_cost_arb = MY_INF;

  for (const Arc &a : P.ordered_arcs[currNode]) {
    DNode next = P.g.target(a);
    bool is_pickup = P.is_pickup[next];
    // Can only go to next if not visited, and it's a pickup or it's
    // corresponding pickup has already been visited:
    if (!visited[next] && (is_pickup || p_visited[P.del_pickup[next]])) {
      if (P.weight[a] < min_cost) { // min arc over all neighbours
        min_next = next;
        min_cost = P.weight[a];
      }
      if (P.weight[a] < min_cost_arb && solver.arborescence(a)) { // restricted min
        min_next_arb = next;
        min_cost_arb = P.weight[a];
        break; // the arcs are sorted, so the first one is the minimum
      }
    }
  }

  // Choose the min arc giving preference to the restricted one:
  DNode &nextNode = min_next;
  double eps = EPS_MIN + (EPS_MAX - EPS_MIN) * pos / (2 * P.npairs - 1); // linear from min to max
  if (min_cost_arb != MY_INF && min_cost_arb <= min_cost * eps)
    nextNode = min_next_arb;
  _arborescence_transversal(P, solver, Sol, nextNode, visited, p_visited, ++pos);
}

double arborescence_transversal(Pickup_Delivery_Instance &P, DNodeVector &Sol,
                                MinCostArb &solver) {
  map<DNode, bool> p_visited; // if each pickup has already been visited
  for (const auto &key : P.pickup) p_visited[key] = false; // init the map
  DNodeBoolMap visited(P.g, false); // if each node has already been visited

  // The target is always the last node:
  Sol[P.nnodes - 1] = P.target;
  visited[P.target] = true;

  // Transverse the graph from the source greedily guided by the arborescence:
  _arborescence_transversal(P, solver, Sol, P.source, visited, p_visited, 0);
  return route_cost(P, Sol); // Calculate the route cost
}

bool arborescence_heuristic(Pickup_Delivery_Instance &P, double &LB, double &UB,
                            DNodeVector &Sol, MinCostArb &solver) {
  DNodeVector arbSol(P.nnodes); // starts the solution vector
  double newUB = arborescence_transversal(P, arbSol, solver);
  if (newUB < UB) { // check if this is a better solution
    UB = newUB;
    Sol = arbSol;
    _NEW_UB_MESSAGE(Sol, "Novo UB pela heuristica da arborescência.");
    cout << "LB: " << LB << endl << endl;
    return true;
  }
  return false;
}
