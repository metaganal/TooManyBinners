# TooManyBinners binning pipeline


Pipeline tool built in python used to generate metagenomic bins with a variety of different binning tools for MAG generation. The aim of the tool is to simplify the MAG generation process by minimising the configurations required and the individual steps needed. 

It can generate bin sets from either just one binning tool or all of them, which can be used for evaluating the best binning tool for the dataset, or used downstream for ensemble binning - this tool will also soon include an ensemble binner.

To be used as singularity image (singularity definition file included in the repo).


## Currently can generate bins from the following binners:

- Vamb
- SemiBin2
- CONCOCT
- MaxBin2
- Metabat2

Able to be used from contigs or reads (reads generated through Metaspades, alignment information for binning is generated through Bowtie2).


## Takes the following arguments:

Required:
- "-t", "--threads". Number of threads to use.
- "-fw", "--forward-reads". Path to forward reads.
- "-rev", "--reverse-reads", Path to reverse reads.
- "-b", "--individual-binners", Individual binners with commas, choices are currently: Semibin2,Maxbin2,Metabat2,Vamb,CONCOCT.
- "-o", "--output-directory", Output directory path.

Optional:
- "-contigs", "--contig-path", path to contigs if starting from contigs rather than reads.
- "--minimum-contig-length", Minimum size of contigs used for binning, default is 2000. Lower size can potentially increase bins generated but also increase contamination.
- "-c", "--custom-kmer-lengths", Custom kmer lengths for metaspades assembly (will default to auto if not). Option for metaspades assembly which can impact assembly results, tweaking this can improve results but automatic does fine.
- "-us", "--using-scaffolds", Whether to use the scaffolds or contigs produced by assembly when running metaspades assembly option.


## Example:

```
singularity run --cleanenv TooManyBinners.sif -fw samplereads_1_P1.fastq.gz \
-rev samplereads_1_P2.fastq.gz \
-b Semibin2,Maxbin2,Metabat2,Vamb,CONCOCT \
-o test_output/ -t 10 \
-us True
```


## Output:

Output currently is a final bins directory consisting of the output bin sets generated from each binning tool.
