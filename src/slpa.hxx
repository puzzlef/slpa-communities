#pragma once
#include <limits>
#include <utility>
#include <vector>
#include "_main.hxx"

using std::pair;
using std::tuple;
using std::vector;
using std::numeric_limits;
using std::make_pair;
using std::move;
using std::get;




// SLPA-OPTIONS
// ------------

// Maximum community membership memory per vertex.
#define SLPA_MEMORY 16


struct SlpaOptions {
  int   repeat;
  float tolerance;
  int   maxIterations;

  SlpaOptions(int repeat=1, float tolerance=0.05, int maxIterations=100) :
  repeat(repeat), tolerance(tolerance), maxIterations(maxIterations) {}
};




// SLPA-RESULT
// -----------

template <class K>
struct SlpaResult {
  vector<K> membership;
  int   iterations;
  float time;

  SlpaResult(vector<K>&& membership, int iterations=0, float time=0) :
  membership(membership), iterations(iterations), time(time) {}

  SlpaResult(vector<K>& membership, int iterations=0, float time=0) :
  membership(move(membership)), iterations(iterations), time(time) {}
};




// LABELSET
// --------

template <class K, size_t L>
using Labelset = array<K, L>;




// SLPA-INITIALIZE
// ---------------

/**
 * Initialize communities such that each vertex is its own community.
 * @param vcom community set each vertex belongs to (updated)
 * @param x original graph
 */
template <class G, class K, size_t L>
inline void slpaInitialize(vector<Labelset<K, L>>& vcom, const G& x) {
  x.forEachVertexKey([&](auto u) { vcom[u] = {u}; });
}


/**
 * Initialize communities from given initial communities.
 * @param vcom community set each vertex belongs to (updated)
 * @param x original graph
 * @param q initial community each vertex belongs to
 */
template <class G, class K, size_t L>
inline void splaInitializeFrom(vector<Labelset<K, L>>& vcom, const G& x, const vector<K>& q) {
  x.forEachVertexKey([&](auto u) { vcom[u] = {q[u]}; });
}




// SLPA-CHOOSE-COMMUNITY
// ---------------------

/**
 * Scan an edge community connected to a vertex.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 * @param u given vertex
 * @param v outgoing edge vertex
 * @param w outgoing edge weight
 * @param vcom community set each vertex belongs to
 * @param l number of labels available
 * @param fr random number generator
 */
template <bool SELF=false, class K, class V, size_t L, class FR>
inline void slpaScanCommunity(vector<K>& vcs, vector<V>& vcout, K u, K v, V w, const vector<Labelset<K, L>>& vcom, K l, FR fr) {
  if (!SELF && u==v) return;
  K c = vcom[v][fr() % l];
  if (!vcout[c]) vcs.push_back(c);
  vcout[c] += w;
}


/**
 * Scan communities connected to a vertex.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 * @param x original graph
 * @param u given vertex
 * @param vcom community set each vertex belongs to
 * @param l number of labels available
 * @param fr random number generator
 */
template <bool SELF=false, class G, class K, class V, size_t L, class FR>
inline void slpaScanCommunities(vector<K>& vcs, vector<V>& vcout, const G& x, K u, const vector<Labelset<K, L>>& vcom, K l, FR fr) {
  x.forEachEdge(u, [&](auto v, auto w) { slpaScanCommunity<SELF>(vcs, vcout, u, v, w, vcom, l, fr); });
}


/**
 * Clear communities scan data.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 */
template <class K, class V>
inline void slpaClearScan(vector<K>& vcs, vector<V>& vcout) {
  for (K c : vcs)
    vcout[c] = V();
  vcs.clear();
}


/**
 * Choose connected community with most weight.
 * @param x original graph
 * @param u given vertex
 * @param vcs communities vertex u is linked to
 * @param vcout total edge weight from vertex u to community C
 * @returns best community
 */
template <bool STRICT=false, class K, class V>
inline K slpaChooseCommunity(const vector<K>& vcs, const vector<V>& vcout) {
  K cmax = K();
  V wmax = V();
  for (K c : vcs) {
    // Do some basic randomization if multiple labels have max weight.
    if (vcout[c]>wmax || (!STRICT && vcout[c]==wmax && (c & 2))) { cmax = c; wmax = vcout[c]; }
  }
  return cmax;
}




// SLPA-BEST-COMMUNITIES
// ---------------------

template <class K, size_t L>
inline void slpaSortCommunities(vector<Labelset<K, L>>& vcom, K l) {
  K S = vcom.size();
  vector<K> a(S);
  for (K u=0; u<S; ++u) {
    auto ib = vcom[u].begin();
    auto ie = vcom[u].begin() + l;
    sort_values(ib, ie);
  }
}


template <class K, size_t L>
inline K slpaBestCommunity(const Labelset<K, L>& labs, K l) {
  return most_frequent(labs.begin(), labs.begin()+l);
}


template <class K, size_t L>
inline vector<K> slpaBestCommunities(const vector<Labelset<K, L>>& vcom, K l) {
  K S = vcom.size();
  vector<K> a(S);
  for (K u=0; u<S; ++u)
    a[u] = slpaBestCommunity(vcom[u], l);
  return a;
}




// SLPA-AFFECTED-VERTICES-DELTA-SCREENING
// --------------------------------------
// Using delta-screening approach.
// - All edge batches are undirected, and sorted by source vertex-id.
// - For edge additions across communities with source vertex `i` and highest modularity changing edge vertex `j*`,
//   `i`'s neighbors and `j*`'s community is marked as affected.
// - For edge deletions within the same community `i` and `j`,
//   `i`'s neighbors and `j`'s community is marked as affected.

/**
 * Find the vertices which should be processed upon a batch of edge insertions and deletions.
 * @param x original graph
 * @param deletions edge deletions for this batch update (undirected, sorted by source vertex id)
 * @param insertions edge insertions for this batch update (undirected, sorted by source vertex id)
 * @param vcom community set each vertex belongs to (sorted)
 * @param l number of labels available
 * @param fr random number generator
 * @returns flags for each vertex marking whether it is affected
 */
template <bool STRICT=false, class FLAG=bool, class G, class K, class V, size_t L, class FR>
auto slpaAffectedVerticesDeltaScreening(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<Labelset<K, L>>& vcom, K l, FR fr) {
  K S = x.span();
  vector<K> vcs; vector<V> vcout(S);
  vector<FLAG> vertices(S), neighbors(S), communities(S);
  for (const auto& [u, v] : deletions) {
    K cu = slpaBestCommunity(vcom[u], l);
    K cv = slpaBestCommunity(vcom[v], l);
    if (cu!=cv) continue;
    vertices[u]  = true;
    neighbors[u] = true;
    communities[cv] = true;
  }
  for (size_t i=0; i<insertions.size();) {
    K u  = get<0>(insertions[i]);
    K cu = slpaBestCommunity(vcom[u], l);
    slpaClearScan(vcs, vcout);
    for (; i<insertions.size() && get<0>(insertions[i])==u; ++i) {
      K v  = get<1>(insertions[i]);
      V w  = get<2>(insertions[i]);
      K cv = slpaBestCommunity(vcom[v], l);
      if (cu==cv) continue;
      slpaScanCommunity(vcs, vcout, u, v, w, vcom, l, fr);
    }
    K cl = slpaChooseCommunity<STRICT>(vcs, vcout);
    if (cl==cu) continue;
    vertices[u]  = true;
    neighbors[u] = true;
    communities[cl] = true;
  }
  x.forEachVertexKey([&](auto u) {
    K cu = slpaBestCommunity(vcom[u], l);
    if (neighbors[u]) x.forEachEdgeKey(u, [&](auto v) { vertices[v] = true; });
    if (communities[cu]) vertices[u] = true;
  });
  return vertices;
}




// SLPA-AFFECTED-VERTICES-FRONTIER
// -------------------------------
// Using frontier based approach.
// - All source and destination vertices are marked as affected for insertions and deletions.
// - For edge additions across communities with source vertex `i` and destination vertex `j`,
//   `i` is marked as affected.
// - For edge deletions within the same community `i` and `j`,
//   `i` is marked as affected.
// - Vertices whose communities change in local-moving phase have their neighbors marked as affected.

/**
 * Find the vertices which should be processed upon a batch of edge insertions and deletions.
 * @param x original graph
 * @param deletions edge deletions for this batch update (undirected, sorted by source vertex id)
 * @param insertions edge insertions for this batch update (undirected, sorted by source vertex id)
 * @param vcom community set each vertex belongs to (sorted)
 * @param l number of labels available
 * @returns flags for each vertex marking whether it is affected
 */
template <class FLAG=bool, class G, class K, class V, size_t L, class FR>
auto slpaAffectedVerticesFrontier(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<Labelset<K, L>>& vcom, K l) {
  K S = x.span();
  vector<FLAG> vertices(S);
  for (const auto& [u, v] : deletions) {
    K cu = slpaBestCommunity(vcom[u], l);
    K cv = slpaBestCommunity(vcom[v], l);
    if (cu!=cv) continue;
    vertices[u] = true;
  }
  for (const auto& [u, v, w] : insertions) {
    K cu = slpaBestCommunity(vcom[u], l);
    K cv = slpaBestCommunity(vcom[v], l);
    if (cu==cv) continue;
    vertices[u] = true;
  }
  return vertices;
}
