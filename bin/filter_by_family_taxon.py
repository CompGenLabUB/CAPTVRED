import re
import gzip
import argparse
from Bio import SeqIO

def get_ids_ofinterest(fam_list_file, fulltaxon_file):
    # Read the list of IDs to include
    with open(fam_list_file, 'r') as f:
      fams=set(line.strip() for line in f)
      print(fams)
      with gzip.open(fulltaxon_file, 'r') as t:
        ids_to_include=set()
        for line in t:
            line_str=line.decode('utf-8')
            ary=line_str.split("\t")
            if any(elem in ary for elem in fams):
                ids_to_include.add(ary[0])
        return ids_to_include



def split_fasta(fasta_file, ids_to_include, FOI_seqs, OTHER_seqs):
    try:
        with gzip.open(fasta_file, 'rt') as FASTA:
            # Open the output files
            #with gzip.open(FOI_seqs, 'wb') as FOI_handle, gzip.open(OTHER_seqs, 'wb') as OTHER_handle:
            with gzip.open(FOI_seqs, 'wt') as FOI_handle, gzip.open(OTHER_seqs, 'wt') as OTHER_handle:
            # Parse the FASTA file and write records to the appropriate output file
                for record in SeqIO.parse(FASTA, 'fasta'):
                    print(record.id)
                    my_id=record.id.split("|")[2]
                    if record.id in ids_to_include:
                        SeqIO.write(record, FOI_handle, 'fasta')
                    else:
                        SeqIO.write(record, OTHER_handle, 'fasta')
            # Close the output files
            FOI_handle.close()
            OTHER_handle.close()
    except IOError as e:
        print(f"Error: {e}")

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="This script splits a gzipped FASTA file based on a list of ID of interest (FOI and OTHER).")
    parser.add_argument("--fasta_file", help="Gzipped FASTA file")
    parser.add_argument("--fam_list_file", help="List of families of interest")
    parser.add_argument("--fulltaxon_file", help="Taxonomizator output (GZ format)")
    parser.add_argument("--FOI_seqs", help="Output gzipped FASTA file for IDs of interest")
    parser.add_argument("--OTHER_seqs", help="Output gzipped FASTA file for IDs not included in the list")

    # Parse the arguments
    args = parser.parse_args()

    # Call the function with parsed arguments
    ids_list=get_ids_ofinterest(args.fam_list_file, args.fulltaxon_file)
    print ("N ids=", len(ids_list))
    split_fasta(args.fasta_file, ids_list, args.FOI_seqs, args.OTHER_seqs)

if __name__ == "__main__":
    main()