#!/usr/bin/env python3
import argparse

from binning import setup_binning
from binning import run_binning
# Needs to take the following arguments:
# threads,contigs,reads,individual_binners,minimum_contig_length
# 
def main():
    
    parser = argparse.ArgumentParser(description="TooManyBinners spoil the broth?")
    
    parser.add_argument("-t", "--threads", help="Threads", required=True)
    parser.add_argument("-fw", "--forward-reads", help="Forward read path", required=True)
    parser.add_argument("-rev", "--reverse-reads", help="reverse read path", required=True)
    parser.add_argument("-contigs", "--contig-path", help="contig file path", required=True)
    parser.add_argument("-b", "--individual-binners", help="Pick individual binners with commas, choices are: Semibin2,Maxbin2,Metabat2", required=True) # KEEP UPDATING THIS
    parser.add_argument("-o", "--output-directory", help="Output directory_path", required=True)
    
    args = parser.parse_args()
    contig_abundance_gen,the_binner = setup_binning(args)

    run_binning(args.output_directory, the_binner, ",".join(args.individual_binners))

# first needs to generate coverage depths file


# Finally running of ensemble binner



main()