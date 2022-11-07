#include <utility>
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>
#include "src/main.hxx"

using namespace std;




// You can define datatype with -DTYPE=...
#ifndef TYPE
#define TYPE float
#endif




template <class G, class K, class V>
double getModularity(const G& x, const SlpaResult<K>& a, V M) {
  auto fc = [&](auto u) { return a.membership[u]; };
  return modularityBy(x, fc, M, V(1));
}


template <class G>
void runExperiment(const G& x, int repeat) {
  using K = typename G::key_type;
  using V = typename G::edge_value_type;
  vector<K> *init = nullptr;
  auto M = edgeWeight(x)/2;
  auto Q = modularity(x, M, 1.0f);
  printf("[%01.6f modularity] noop\n", Q);
  SlpaOptions o = {repeat};

  float tolerance = 0.05f;
  {
    // Find SLPA using a single thread (4 labels).
    auto ak = slpaSeqStatic<4, false>(x, init, {repeat, tolerance});
    printf("[%09.3f ms; %04d iters.; %01.9f modularity] slpaSeqStatic       {labels=%02d, tolerance=%.0e}\n", ak.time, ak.iterations, getModularity(x, ak, M), 4, tolerance);
    auto al = slpaSeqStatic<4, true>(x, init, {repeat, tolerance});
    printf("[%09.3f ms; %04d iters.; %01.9f modularity] slpaSeqStaticStrict {labels=%02d, tolerance=%.0e}\n", al.time, al.iterations, getModularity(x, al, M), 4, tolerance);
  }
  {
    // Find SLPA using a single thread (8 labels).
    auto ak = slpaSeqStatic<8, false>(x, init, {repeat, tolerance});
    printf("[%09.3f ms; %04d iters.; %01.9f modularity] slpaSeqStatic       {labels=%02d, tolerance=%.0e}\n", ak.time, ak.iterations, getModularity(x, ak, M), 8, tolerance);
    auto al = slpaSeqStatic<8, true>(x, init, {repeat, tolerance});
    printf("[%09.3f ms; %04d iters.; %01.9f modularity] slpaSeqStaticStrict {labels=%02d, tolerance=%.0e}\n", al.time, al.iterations, getModularity(x, al, M), 8, tolerance);
  }
  {
    // Find SLPA using a single thread (16 labels).
    auto ak = slpaSeqStatic<16, false>(x, init, {repeat, tolerance});
    printf("[%09.3f ms; %04d iters.; %01.9f modularity] slpaSeqStatic       {labels=%02d, tolerance=%.0e}\n", ak.time, ak.iterations, getModularity(x, ak, M), 16, tolerance);
    auto al = slpaSeqStatic<16, true>(x, init, {repeat, tolerance});
    printf("[%09.3f ms; %04d iters.; %01.9f modularity] slpaSeqStaticStrict {labels=%02d, tolerance=%.0e}\n", al.time, al.iterations, getModularity(x, al, M), 16, tolerance);
  }
  {
    // Find SLPA using a single thread (32 labels).
    auto ak = slpaSeqStatic<32, false>(x, init, {repeat, tolerance});
    printf("[%09.3f ms; %04d iters.; %01.9f modularity] slpaSeqStatic       {labels=%02d, tolerance=%.0e}\n", ak.time, ak.iterations, getModularity(x, ak, M), 32, tolerance);
    auto al = slpaSeqStatic<32, true>(x, init, {repeat, tolerance});
    printf("[%09.3f ms; %04d iters.; %01.9f modularity] slpaSeqStaticStrict {labels=%02d, tolerance=%.0e}\n", al.time, al.iterations, getModularity(x, al, M), 32, tolerance);
  }
  {
    // Find SLPA using a single thread (64 labels).
    auto ak = slpaSeqStatic<64, false>(x, init, {repeat, tolerance});
    printf("[%09.3f ms; %04d iters.; %01.9f modularity] slpaSeqStatic       {labels=%02d, tolerance=%.0e}\n", ak.time, ak.iterations, getModularity(x, ak, M), 64, tolerance);
    auto al = slpaSeqStatic<64, true>(x, init, {repeat, tolerance});
    printf("[%09.3f ms; %04d iters.; %01.9f modularity] slpaSeqStaticStrict {labels=%02d, tolerance=%.0e}\n", al.time, al.iterations, getModularity(x, al, M), 64, tolerance);
  }
}


int main(int argc, char **argv) {
  using K = int;
  using V = TYPE;
  char *file = argv[1];
  int repeat = argc>2? stoi(argv[2]) : 5;
  OutDiGraph<K, None, V> x;  // V w = 1;
  printf("Loading graph %s ...\n", file);
  readMtxW<true>(x, file); println(x);
  auto y = symmetricize(x); print(y); printf(" (symmetricize)\n");
  // auto fl = [](auto u) { return true; };
  // selfLoopU(y, w, fl); print(y); printf(" (selfLoopAllVertices)\n");
  runExperiment(y, repeat);
  printf("\n");
  return 0;
}
