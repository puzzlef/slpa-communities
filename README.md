Single-threaded CPU-based Speaker-listener Label Propagation Algorithm ([SLPA]) for [community detection].

I was trying out the **Speaker-listener Label Propagation Algorithm (SLPA)** which i had implemented yesterday (results out today). This is one of the algorithms given in literature review of *LabelRank*. Here each node has a *fixed memory* where it remembers all the popular labels that it **listened**. The size of this memory is fixed to the *total number of iterations to be performed + 1*. The algorithm is as follows.

1. Each vertex is initialized such that it remembers itself as *popular*.
2. Each neighbor speaks one of the *random labels* in its memory.
3. The vertex (listener) adds the most popular label to its memory.
4. Repeat from 2 until a *fixed number of iterations* is performed (labels - 1).
5. I allow *early convergence* if *at least n%* of vertices remember their *previous label*.
6. For each vertex, i pick the *most popular label* in its memory as its *community*.

Unlike **RAK** (also called **LPA**) and **COPRA**, this is a *randomized* *historical-label* based technique. A *random number generator (RNG)* is used by each neighbor of a vertex to pick a random label to **speak** from its memory. I originally used C++'s *default random engine* (which is probably mersenne twister) which seems slow. So i updated it to use a *xorshift32 random engine*. The listener vertex listens to all its neigbors, and based on edge weight, picks the *most popular label*. In **strict** mode, the listener pick the first *most* *popular label* as to record in its memory, and in **non-strict** mode, the listener pick randomly one of the most popular labels to record. It records this in the *next index* in its memory. The authors of this paper suggest to perform a *fixed number of iterations* (which is one less than the space we have reserved for storing the labels in the memory of each vertex). However, if the latest popular label (for a vertex) is same as the previous one in its memory, i consider it as a *sign of convergence*. So, when *n% of vertices* show this sign of convergence (**tolerance**), i perform an *early exit*. By default, a **tolerance** of `0.05` is used. As we want to find *disjoint communities*, i use the most popular label in the memory of each vertex as its final community.

[![](https://i.imgur.com/dq8i8HX.png)][sheetp]

[![](https://i.imgur.com/eTBXvF8.png)][sheetp]

<br>
<br>

The experiment is run with **labels** (memory space of each vertex) ranging from `4` to `64`.

[![](https://i.imgur.com/wNqbbqO.png)][sheetp]

[![](https://i.imgur.com/sqP2n0M.png)][sheetp]

[![](https://i.imgur.com/gGLZXwA.png)][sheetp]

<br>
<br>

It appears that *increasing* the number of **labels** *increases the time required for* *completion* (as expected), and also *increases the final modularity*. However, on graphs `soc-Slashdot*` modularity seems to decrease with increasing labels.

[![](https://i.imgur.com/9CYJIeF.png)][sheetp]

<br>
<br>

It is also observed that the **strict** variant is *slightly faster* than the **non-strict** one, but achieves *similar modularity*. On *road networks* it seems **strict** approach obtains *lower modularity*, and on other types of graphs, **strict** approach achieves *higher modularity*.

[![](https://i.imgur.com/P0WfY3A.png)][sheetp]

[![](https://i.imgur.com/reyxFuk.png)][sheetp]

<br>
<br>

However, its seems that this method of community detection does not provide a good return on investment (it takes longer but does not achieve great modularity). I do not perform any post-processing of communities as of now (such as splitting disconnected communities).

All outputs are saved in a [gist] and a small part of the output is listed here. Some [charts] are also included below, generated from [sheets]. The input data used for this experiment is available from the [SuiteSparse Matrix Collection]. This experiment was done with guidance from [Prof. Kishore Kothapalli] and [Prof. Dip Sankar Banerjee].


[SLPA]: https://arxiv.org/abs/1202.2465
[COPRA]: https://arxiv.org/abs/0910.5516
[RAK]: https://arxiv.org/abs/0709.2938
[community detection]: https://en.wikipedia.org/wiki/Community_search
[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu

<br>

```bash
$ g++ -std=c++17 -O3 main.cxx
$ ./a.out ~/data/web-Stanford.mtx
$ ./a.out ~/data/web-BerkStan.mtx
$ ...

# Loading graph /home/subhajit/data/web-Stanford.mtx ...
# order: 281903 size: 2312497 [directed] {}
# order: 281903 size: 3985272 [directed] {} (symmetricize)
# [-0.000497 modularity] noop
# [00337.151 ms; 0004 iters.; 0.474440932 modularity] slpaSeqStatic       {labels=04, tolerance=5e-02}
# [00304.351 ms; 0004 iters.; 0.479095578 modularity] slpaSeqStaticStrict {labels=04, tolerance=5e-02}
# [00901.554 ms; 0008 iters.; 0.752845287 modularity] slpaSeqStatic       {labels=08, tolerance=5e-02}
# [00866.540 ms; 0008 iters.; 0.679729104 modularity] slpaSeqStaticStrict {labels=08, tolerance=5e-02}
# [02229.548 ms; 0016 iters.; 0.826291323 modularity] slpaSeqStatic       {labels=16, tolerance=5e-02}
# [02193.593 ms; 0016 iters.; 0.789462328 modularity] slpaSeqStaticStrict {labels=16, tolerance=5e-02}
# [03706.551 ms; 0032 iters.; 0.847930729 modularity] slpaSeqStatic       {labels=32, tolerance=5e-02}
# [03655.065 ms; 0032 iters.; 0.817599833 modularity] slpaSeqStaticStrict {labels=32, tolerance=5e-02}
# [07928.522 ms; 0064 iters.; 0.857062638 modularity] slpaSeqStatic       {labels=64, tolerance=5e-02}
# [07843.098 ms; 0064 iters.; 0.826945186 modularity] slpaSeqStaticStrict {labels=64, tolerance=5e-02}
#
# ...
```

<br>
<br>


## References

- [Towards Linear Time Overlapping Community Detection in Social Networks; Jierui Xie et al. (2012)](https://link.springer.com/chapter/10.1007/978-3-642-30220-6_3)
- [Finding overlapping communities in networks by label propagation; Steve Gregory (2010)](https://iopscience.iop.org/article/10.1088/1367-2630/12/10/103018)
- [Near linear time algorithm to detect community structures in large-scale networks; Usha Nandini Raghavan et al. (2007)](https://arxiv.org/abs/0709.2938)
- [The University of Florida Sparse Matrix Collection; Timothy A. Davis et al. (2011)](https://doi.org/10.1145/2049662.2049663)
- [How to import VSCode keybindings into Visual Studio?](https://stackoverflow.com/a/62417446/1413259)
- [Configure X11 Forwarding with PuTTY and Xming](https://www.centlinux.com/2019/01/configure-x11-forwarding-putty-xming-windows.html)
- [Installing snap on CentOS](https://snapcraft.io/docs/installing-snap-on-centos)

<br>
<br>


[![](https://i.imgur.com/1dNrrfK.jpg)](https://www.youtube.com/watch?v=3X85rHyfg0k)<br>
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)
[![DOI](https://zenodo.org/badge/562802206.svg)](https://zenodo.org/badge/latestdoi/562802206)
![](https://ga-beacon.deno.dev/G-KD28SG54JQ:hbAybl6nQFOtmVxW4if3xw/github.com/puzzlef/slpa-communities)

[gist]: https://gist.github.com/wolfram77/8730879a0220d86091fc8e7e8e4f3b5d
[charts]: https://imgur.com/a/WV5T0cQ
[sheets]: https://docs.google.com/spreadsheets/d/1wfRvU5t_F53PP7el5gp_i5vjgcx2cNOxatXmA8tO9PE/edit?usp=sharing
[sheetp]: https://docs.google.com/spreadsheets/d/e/2PACX-1vRLFqylFOh1JvkRTIq1zGx7Vpy9zb6DhyV1WANf82NV2rcKp9JwkQBlkBaTXmtmwFqgqTFyFHnoQ5qU/pubhtml
