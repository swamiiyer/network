#!/usr/bin/python
#
# This script prints as output the values of some of the basic properties of 
# the input network in GraphML/GML format.
#
# Usage: python network_stats.py <infile>

import networkx, numpy, sys

def main(args):
    if len(args) == 0:
        print "Usage: python network_stats.py <infile>"
        sys.exit()
    infile = args[0]
    if infile.endswith(".gml"):
        g = networkx.read_gml(infile)
    elif infile.endswith("graphml"):
        g = networkx.read_graphml(infile)
    else:
        print "Unknown file extension"
        sys.exit(1)
    n = len(g.nodes())
    m = len(g.edges())
    print "number of vertices = %d" %(n)
    print "number of edges = %d" %(m)
    print "connectance = %f" %(1.0 * m / (n * (n - 1) / 2))
    print "average degree = %f" %(numpy.average(g.degree().values()))
    components = networkx.connected_components(g)
    print "number of components = %d" %(len(components))
    print "fractional size of largest component = %f" \
        %(max([len(i) * 1. for i in components]) / \
              len(g.nodes()))
    print "average clustering coefficient = %f" \
        %(networkx.average_clustering(g))
    print "assortativity coefficient = %f" \
        %(networkx.degree_assortativity_coefficient(g))

if __name__ == "__main__":
    main(sys.argv[1:])
