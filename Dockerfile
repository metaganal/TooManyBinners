FROM alpine
LABEL TooManyBinners_version=0.1.1


RUN apt-get update && apt-get install -y wget git bash gcc gfortran g++ make file python3 nano curl


RUN wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh" -O mamba_file.sh
RUN bash mamba_file.sh -b -f -p /opt/mamba
RUN PATH=/opt/mamba/bin:$PATH
RUN mamba init

RUN mamba install -c conda-forge psutil
RUN mamba create -n prebinning -c conda-forge -c bioconda bowtie2 samtools spades
RUN mamba create -n CONCOCTMetabat2MaxBin2SemiBin2 -c conda-forge -c bioconda concoct maxbin2 metabat2 semibin scikit-learn=1.1

RUN mamba create -n Vamb4 pip
RUN bash -c "source /opt/mamba/bin/activate Vamb4 && pip install vamb"
RUN wget https://raw.githubusercontent.com/RasmussenLab/vamb/master/src/create_fasta.py -O /opt/create_fasta.py
    git clone https://www.github.com/lobrien20/TooManyBinners -b tests
    mv TooManyBinners /opt/
