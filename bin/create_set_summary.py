#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import gzip

def getinfo(fulltaxon_file, genomes_coords, output_file):
    coords={}
    with open(output_file, 'w') as OUT:
        OUT.write("#Id\tName\tTax\tStart\tEnd\tDescription\tSpecID\tSpecTax\tGenusId\tGenusTax\tFamId\tFamTax\n")
        with open(genomes_coords, 'r') as c:
            for line in c:
                line_str=line.strip()
                ary=line_str.split("\t")
                idf=str(ary[0]).strip(" ")
                st=str(ary[1])
                end=str(ary[2])
                coords[idf]=[st, end]
            # print (coords)
        with gzip.open(fulltaxon_file, 'r') as t:
            for line in t:
                line_str=line.decode('utf-8')
                ary=line_str.split("\t")
                id=str(ary[0])
            # print (id)
            #    print(coords[id][1])
                taxid=ary[1]
                desc=ary[2]
                for i in range(3,len(ary)):
                    info=ary[i].split(":")
                    if info[0]=="family":
                        fam=info[1]
                        famtax=info[2].strip("\n")
                    elif info[0]=="genus":
                        gen=info[1]
                        gentax=info[2].strip("\n")
                    elif info[0]=="species":
                        spec=info[1]
                        spectax=info[2].strip("\n")
                OUT.write("\t".join([id, desc, taxid, coords[id][0], coords[id][1], desc, spec, spectax,  gen, gentax, fam, famtax,"\n"])) #[desc, id, taxid, coords[id][0], coords[id][1], desc, spec, spectax,  gen, gentax, fam, famtax)
    OUT.close()
def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="This script Creates a summary tsv file required for the data visualization step on the CAPTVRED pipeline.")
    parser.add_argument("--fulltaxon_file", help="Taxonomizator output (GZ format)")
    parser.add_argument("--genomes_coords", help="Genome coordinates")
    parser.add_argument("--output_file", help="Output file")
    
    # Parse the arguments
    args = parser.parse_args()

    # Call the function to get summary table
    getinfo(args.fulltaxon_file, args.genomes_coords, args.output_file)


if __name__ == "__main__":
    main()