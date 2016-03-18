#!/usr/bin/python

# This script converts a network in graphml format to one in gml format.

import sys, networkx

def main(argv):
    """
    Entry point.
    """
    
    if len(argv) != 1:
        print "Usage: python graphml2gml.py <infile>"
        sys.exit(0)
    infile = argv[0]
    g = networkx.read_graphml(infile)
    networkx.write_gml(networkx.from_edgelist(g.edges()), 
                       infile.replace(".graphml", ".gml"))

if __name__ == "__main__":
    main(sys.argv[1:])
