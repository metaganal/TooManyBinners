TooManyBinners binning pipeline

Pipeline tool built in python used to generate metagenomic bins with a variety of different binning tools for MAG generation. The aim of the tool is to simplify the MAG generation process by minimising the configurations required and the individual steps needed. 

Currently can generate bins from the following binners:
- Vamb
- SemiBin2
- CONCOCT
- MaxBin2
- Metabat2

Able to be used from contigs or reads (reads generated through Metaspades, alignment information for binning is generated through Bowtie2).

Takes the following arguments:

"-t", "--threads", help="Threads", required=True
"-fw", "--forward-reads", help="Forward read path", required=True
"-rev", "--reverse-reads", help="reverse read path", required=True
"-contigs", "--contig-path", help="contig file path, if not provided will auto assemble"
"--minimum-contig-length", help="minium size of contigs for binning. Default is 2000."
"-b", "--individual-binners", help="Pick individual binners with commas, choices are: Semibin2,Maxbin2,Metabat2,Vamb,CONCOCT"
"-o", "--output-directory", help="Output directory_path", required=True
"-c", "--custom-kmer-lengths", help="Custom kmer lengths for metaspades assembly (will default to auto if not)".
"-us", "--using-scaffolds", help="Using metaspades assembly scaffolds?"

