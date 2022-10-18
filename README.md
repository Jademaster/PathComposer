# PathComposer

Implements the algorithm described in [How to Compose Shortest Paths](https://arxiv.org/abs/2205.15306). This repository consists of:

* compositionpatterns.py
* semiring.py
* compositionproblem.py
* main.py

compositionpatterns.py generates the "composition symbols" as defined in the above paper. semiring.py defines classes for semirings and matrix semirings as well as necessary methods for them. compositionproblem.py defines classes for graphs and "composition problems" i.e. pairs of graphs which share a common set of vertices. Lastly, main.py gives an example of the use of the other three scripts. A random composition problem is generated and the shortest paths are computed using the compositional algorithm and compared to the networkx algorithm.
