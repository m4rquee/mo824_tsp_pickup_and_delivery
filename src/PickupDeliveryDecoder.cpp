#include "PickupDeliveryDecoder.hpp"

PickupDeliveryDecoder::PickupDeliveryDecoder(Pickup_Delivery_Instance &P)
    : P(P) {}

PickupDeliveryDecoder::~PickupDeliveryDecoder() = default;

double PickupDeliveryDecoder::decode(const std::vector<double> &chromosome) const {
  DNodeVector sol(P.nnodes);
  decode(chromosome, sol);
  return route_cost(P, sol);
}

void PickupDeliveryDecoder::decode(const vector<double> &chromosome, DNodeVector &Sol) const {
  map<DNode, bool> p_visited; // if each pickup has already been visited
  for (const auto &key : P.pickup) p_visited[key] = false; // init the map

  // The endpoints are fixed:
  Sol[0] = P.source;
  Sol[P.nnodes - 1] = P.target;

  vector<pair<double, unsigned>> pic_ranking(P.npairs);
  vector<pair<double, unsigned>> del_ranking(P.npairs);
  for (int i = 0; i < P.npairs; i++) {
    pic_ranking[i] = pair<double, unsigned>(chromosome[i], i);
    del_ranking[i] = pair<double, unsigned>(chromosome[P.npairs + i], i);
  }

  // This will then produce a permutation of [n] in pair::second:
  sort(pic_ranking.begin(), pic_ranking.end());
  sort(del_ranking.begin(), del_ranking.end());

  // Consume from the sorted rankings to produce a solution:
  int i, p, d;
  for (i = 1, p = 0, d = 0; i <= 2 * P.npairs and p < P.npairs; i++) {
    DNode &currDel = P.delivery[del_ranking[d].second];
    // Puts the node with higher rank first (if is a delivery, then it's
    // corresponding pickup must already been visited):
    if (del_ranking[d].first < pic_ranking[p].first and
        p_visited[P.del_pickup[currDel]]) {
      d++;
      Sol[i] = currDel;
    } else {
      Sol[i] = P.pickup[pic_ranking[p++].second];
      p_visited[Sol[i]] = true;
    }
  }
  // Copy the remaining deliveries:
  for (; d < P.npairs; i++, d++) Sol[i] = P.delivery[del_ranking[d].second];
}
