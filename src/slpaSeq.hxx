#pragma once
#include <utility>
#include <algorithm>
#include <random>
#include <vector>
#include "_main.hxx"
#include "vertices.hxx"
#include "edges.hxx"
#include "csr.hxx"
#include "slpa.hxx"

using std::tuple;
using std::vector;
using std::random_device;
using std::default_random_engine;
using std::uniform_int_distribution;
using std::make_pair;
using std::swap;
using std::min;




// SLPA-MOVE-ITERATION
// -------------------

/**
 * Move each vertex to its best community.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 * @param vcom community set each vertex belongs to (updated)
 * @param x original graph
 * @param l number of labels available
 * @param fr random number generator
 * @returns number of changed vertices
 */
template <bool STRICT=false, class G, class K, class V, size_t L, class FR, class FA, class FP>
K slpaMoveIteration(vector<K>& vcs, vector<V>& vcout, vector<Labelset<K, L>>& vcom, const G& x, K l, FR fr, FA fa, FP fp) {
  K a = K();
  x.forEachVertexKey([&](auto u) {
    if (!fa(u)) return;
    slpaClearScan(vcs, vcout);
    slpaScanCommunities(vcs, vcout, x, u, vcom, l, fr);
    vcom[u][l] = slpaChooseCommunity<STRICT>(vcs, vcout);
    if (vcom[u][l]!=vcom[u][l-1]) { ++a; fp(u); }
  });
  return a;
}




// SLPA-SEQ
// --------

template <size_t LABELS=SLPA_MEMORY, bool STRICT=false, class G, class K, class FA, class FP>
SlpaResult<K> slpaSeq(const G& x, const vector<K>* q, const SlpaOptions& o, FA fa, FP fp) {
  using V = typename G::edge_value_type;
  const size_t L = LABELS;
  int l = 0;
  K S = x.span();
  K N = x.order();
  vector<K> vcs;
  vector<V> vcout(S);
  vector<Labelset<K, L>> vcom(S);
  random_device dev;
  default_random_engine rnd(dev());
  uniform_int_distribution<K> dis(0, 65535);
  auto fr = [&]() { return dis(rnd); };
  float t = measureDuration([&]() {
    if (q) slpaInitializeFrom(vcom, x, *q);
    else   slpaInitialize(vcom, x);
    for (l=0; l < min(o.maxIterations, K(LABELS-1));) {
      K n = slpaMoveIteration<STRICT>(vcs, vcout, vcom, x, ++l, fr, fa, fp);
      PRINTFD("slpaSeq(): l=%d, n=%d, N=%d, n/N=%f\n", l, n, N, float(n)/N);
      if (float(n)/N <= o.tolerance) break;
    } ++l;
  }, o.repeat);
  slpaSortCommunities(vcom, l);
  return {slpaBestCommunities(vcom, l), l, t};
}
template <size_t LABELS=SLPA_MEMORY, bool STRICT=false, class G, class K, class FA>
inline SlpaResult<K> slpaSeq(const G& x, const vector<K>* q, const SlpaOptions& o, FA fa) {
  auto fp = [](auto u) {};
  return slpaSeq<LABELS, STRICT>(x, q, o, fa, fp);
}
template <size_t LABELS=SLPA_MEMORY, bool STRICT=false, class G, class K>
inline SlpaResult<K> slpaSeq(const G& x, const vector<K>* q, const SlpaOptions& o) {
  auto fa = [](auto u) { return true; };
  return slpaSeq<LABELS, STRICT>(x, q, o, fa);
}




// SLPA-SEQ-STATIC
// ---------------

template <size_t LABELS=SLPA_MEMORY, bool STRICT=false, class G, class K>
inline SlpaResult<K> slpaSeqStatic(const G& x, const vector<K>* q=nullptr, const SlpaOptions& o={}) {
  return slpaSeq<LABELS, STRICT>(x, q, o);
}




// SLPA-SEQ-DYNAMIC-DELTA-SCREENING
// --------------------------------

template <size_t LABELS=SLPA_MEMORY, bool STRICT=false, class G, class K, class V>
inline SlpaResult<K> slpaSeqDynamicDeltaScreening(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const SlpaOptions& o={}) {
  const size_t L = LABELS;
  K S = x.span();
  auto vaff = slpaAffectedVerticesDeltaScreening<STRICT>(x, deletions, insertions, *q);
  auto fa   = [&](auto u) { return vaff[u]==true; };
  return slpaSeq<LABELS, STRICT>(x, q, o, fa);
}




// SLPA-SEQ-DYNAMIC-FRONTIER
// -------------------------

template <size_t LABELS=SLPA_MEMORY, bool STRICT=false, class G, class K, class V>
inline SlpaResult<K> slpaSeqDynamicFrontier(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const SlpaOptions& o={}) {
  const size_t L = LABELS;
  K S = x.span();
  auto vaff = slpaAffectedVerticesFrontier(x, deletions, insertions, *q);
  auto fa = [&](auto u) { return vaff[u]==true; };
  auto fp = [&](auto u) { x.forEachEdgeKey(u, [&](auto v) { vaff[v] = true; }); };
  return slpaSeq<LABELS, STRICT>(x, q, o, fa, fp);
}
