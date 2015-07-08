#!/usr/bin/python

# This script converts a network in gml format to one in graphml format.

import sys, networkx

def main(argv):
    """
    Entry point.
    """
    
    if len(argv) != 1:
        print "Usage: python gml2graphml.py <infile>"
        sys.exit(0)
    infile = argv[0]
    g = networkx.read_gml(infile)
    networkx.write_graphml(g, infile.replace(".gml", ".graphml"))

if __name__ == "__main__":
    main(sys.argv[1:])
