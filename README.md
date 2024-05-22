This repository contains the code both for Tutte's algorithm and the first part of the Polygon Universal algorithm (see https://arxiv.org/abs/2103.06916) that creates the triangulation-respecting sketch
Both are implemented in Python using numpy, matplotlib, networkx and sagemaths for the generated graphs. PolygonUniversal uses pygeos to generate delaunay triangulations for the given polygon
Provided it in this repo are also a few illustrations created by both algorithms, as well as files containing 
generated 3-connected planar graphs of different sizes, made by using plantri (https://users.cecs.anu.edu.au/~bdm/plantri/) 

For Tutte's algorithm, execute in the following way:

Tutte.py [0-5] -> execute the algorithm on 6 cases we created ourselves (Various planar embeddings).
Tutte.py 6 file(graphs in Graph6format, one per line) -> draw a randomly chosen graph from the given file using the algorithm
Tutte.py 7 file(graphs in Graph6format, one per line) n_graphs:int -> Compute runtime of the algorithm over n_graphs graphs from file

For Polygon Universal, execute in the following way:

PolygonUniversal.py [0-6] debug[0-1]:
[0,4] debug[0,1]-> various cases
5: PolygonUniversal.py 5 debug[0,1] file(graphs in Graph6format, one per line) -> draw a randomly chosen graph from the given file using the algorithm
6: PolygonUniversal.py 6 debug[0,1] file(graphs in Graph6format, one per line) ngraphs:int -> Compute runtime of the algorithm over n_graphs graphs from file

The debug flag chooses whether to display the triangulation of the polygon and logs during the execution of the algorithm

Note that the algorithm computes the sketches, and then picks random point on each simplex of each vertex for the embedding, therefore result is not always a planar embedding
