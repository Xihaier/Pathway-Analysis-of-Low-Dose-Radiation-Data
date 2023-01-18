#!/usr/bin/env python

# usage:  ./reformat_key.py  <kegg gene-pathway file>

import sys
from pprint import pprint

genes_in_path = {}

if  len( sys.argv ) >= 2:
    with open( sys.argv[1], "r" ) as f:
        for line in f:
            kgene, kpath = line.strip().split()
            gene = kgene.split( ':' )[1]    # remove prefixes
            path = kpath.split( ':' )[1] 
            if path not in genes_in_path:
                genes_in_path[ path ] = []
            genes_in_path[ path ].append( gene )
            
else:
    print( "specify KEGG gene-pathway file on command line" )

# print out table

pathways = genes_in_path.keys() 
pathways.sort()

for p in pathways:
    print( "\t".join( [ p ] + genes_in_path.get( p ) ) )
     


