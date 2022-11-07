#pragma once
#include <utility>
#include <vector>
#include <algorithm>
#include "_main.hxx"
#include "vertices.hxx"
#include "edges.hxx"
#include "csr.hxx"
#include "slpa.hxx"

using std::tuple;
using std::vector;
using std::make_pair;
using std::swap;




// SLPA-MOVE-ITERATION
// -------------------

/**
 * Move each vertex to its best community.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 * @param vcom community set each vertex belongs to (updated)
 * @param x original graph
 * @param vtot total edge weight of each vertex
 * @param B belonging coefficient threshold
 * @returns number of changed vertices
 */
template <class G, class K, class V, size_t L, class FA, class FP>
K slpaMoveIteration(vector<K>& vcs, vector<V>& vcout, vector<Labelset<K, V, L>>& vcom, const G& x, const vector<V>& vtot, V B, FA fa, FP fp) {
  K a = K();
  x.forEachVertexKey([&](auto u) {
    if (!fa(u)) return;
    K d = vcom[u][0].first;
    slpaClearScan(vcs, vcout);
    slpaScanCommunities(vcs, vcout, x, u, vcom);
    slpaSortScan(vcs, vcout);
    vcom[u] = slpaChooseCommunity(x, u, vcom, vcs, vcout, B*vtot[u]);
    K c = vcom[u][0].first;
    if (c!=d) { ++a; fp(u); }
  });
  return a;
}




// SLPA-SEQ
// --------

template <size_t LABELS=SLPA_MAX_MEMBERSHIP, class G, class K, class FA, class FP>
SlpaResult<K> slpaSeq(const G& x, const vector<K>* q, const SlpaOptions& o, FA fa, FP fp) {
  using V = typename G::edge_value_type;
  const size_t L = LABELS;
  int l = 0;
  K S = x.span();
  K N = x.order();
  V B = V(1)/LABELS;
  vector<K> vcs;
  vector<V> vcout(S), vtot(S);
  vector<Labelset<K, V, L>> vcom(S);
  float t = measureDuration([&]() {
    slpaVertexWeights(vtot, x);
    slpaInitialize(vcom, x);
    for (l=0; l<o.maxIterations;) {
      K n = slpaMoveIteration(vcs, vcout, vcom, x, vtot, B, fa, fp); ++l;
      PRINTFD("slpaSeq(): l=%d, n=%d, N=%d, n/N=%f\n", l, n, N, float(n)/N);
      if (float(n)/N <= o.tolerance) break;
    }
  }, o.repeat);
  return {slpaBestCommunities(vcom), l, t};
}
template <size_t LABELS=SLPA_MAX_MEMBERSHIP, class G, class K, class FA>
inline SlpaResult<K> slpaSeq(const G& x, const vector<K>* q, const SlpaOptions& o, FA fa) {
  auto fp = [](auto u) {};
  return slpaSeq<LABELS>(x, q, o, fa, fp);
}
template <size_t LABELS=SLPA_MAX_MEMBERSHIP, class G, class K>
inline SlpaResult<K> slpaSeq(const G& x, const vector<K>* q, const SlpaOptions& o) {
  auto fa = [](auto u) { return true; };
  return slpaSeq<LABELS>(x, q, o, fa);
}




// SLPA-SEQ-STATIC
// ---------------

template <size_t LABELS=SLPA_MAX_MEMBERSHIP, class G, class K>
inline SlpaResult<K> slpaSeqStatic(const G& x, const vector<K>* q=nullptr, const SlpaOptions& o={}) {
  return slpaSeq<LABELS>(x, q, o);
}




// SLPA-SEQ-DYNAMIC-DELTA-SCREENING
// --------------------------------

template <size_t LABELS=SLPA_MAX_MEMBERSHIP, class G, class K, class V>
inline SlpaResult<K> slpaSeqDynamicDeltaScreening(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const SlpaOptions& o={}) {
  const size_t L = LABELS;
  K S = x.span();
  const vector<Labelset<K, V, L>>& vcom = *q;
  auto vaff = slpaAffectedVerticesDeltaScreening(x, deletions, insertions, vcom);
  auto fa   = [&](auto u) { return vaff[u]==true; };
  return slpaSeq<LABELS>(x, q, o, fa);
}




// SLPA-SEQ-DYNAMIC-FRONTIER
// -------------------------

template <size_t LABELS=SLPA_MAX_MEMBERSHIP, class G, class K, class V>
inline SlpaResult<K> slpaSeqDynamicFrontier(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const SlpaOptions& o={}) {
  const size_t L = LABELS;
  K S = x.span();
  const vector<Labelset<K, V, L>>& vcom = *q;
  auto vaff = slpaAffectedVerticesFrontier(x, deletions, insertions, vcom);
  auto fa = [&](auto u) { return vaff[u]==true; };
  auto fp = [&](auto u) { x.forEachEdgeKey(u, [&](auto v) { vaff[v] = true; }); };
  return slpaSeq<LABELS>(x, q, o, fa, fp);
}
