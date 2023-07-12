import subprocess
import os
from utils import run_and_log_a_subprocess

class MetaSpadesAssemblyRunner:
    
    def __init__(self, output_directory, using_scaffolds, read_fwd_path, read_rev_path, threads):
        self.output_directory = output_directory
        self.using_scaffolds = using_scaffolds
        self.read_fwd_path = read_fwd_path
        self.read_rev_path = read_rev_path
        self.threads = threads

    
    def run_assembly(self):
        metaspades_result_directory = f"{self.output_directory}/metaspades_assembly_directory/"
        metaspades_assembly_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/prebinning', 'spades.py', '-1', self.read_fwd_path, '-2', self.read_rev_path, '-t', self.threads, '--meta', '-o', metaspades_result_directory]
        run_and_log_a_subprocess(self.output_directory, metaspades_assembly_args, "metaspades_assembly_auto_kmer")

        if self.using_scaffolds == True:

            print("Using metaspades generated scaffolds as contigs")
            return f"{metaspades_result_directory}/scaffolds.fasta"

        else:

            return f"{metaspades_result_directory}/contigs.fasta"


    def run_assembly_with_custom_kmer_lengths(self, custom_kmer_lengths):
        metaspades_result_directory = f"{self.output_directory}/metaspades_assembly_directory/"
        metaspades_assembly_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/prebinning', 'spades.py', '-1', self.read_fwd_path, '-2', self.read_rev_path, '-t', self.threads, '-k', custom_kmer_lengths, '--meta', '-o', metaspades_result_directory]
        run_and_log_a_subprocess(self.output_directory, metaspades_assembly_args, "metaspades_assembly_auto_kmer")

        if self.using_scaffolds == True:

            print("Using metaspades generated scaffolds as contigs")
            return f"{metaspades_result_directory}/scaffolds.fasta"

        else:

            return f"{metaspades_result_directory}/contigs.fasta"



    def check_if_assembly_step_has_already_been_completed(self):
        metaspades_result_directory = f"{self.output_directory}/metaspades_assembly_directory/"
        if os.path.isdir(metaspades_result_directory):
            if not os.path.exists(f"{metaspades_result_directory}/contigs.fasta") or not os.path.exists(f"{metaspades_result_directory}/scaffolds.fasta"):
                print("Found metaspades directory but on contigs or scaffolds. Error, exiting.")
                exit()
            if self.using_scaffolds == True:
                return f"{metaspades_result_directory}/scaffolds.fasta"
            else:
                return f"{metaspades_result_directory}/contigs.fasta"
        
        return False
        
class ContigAbundances:
    
    
    def __init__(self, output_directory, contig_file_path, read_fwd_path, read_rev_path, threads):
        
        self.output_directory = output_directory
        self.contig_file_path = contig_file_path
        self.read_fwd_path = read_fwd_path
        self.read_rev_path = read_rev_path
        self.threads = threads
        
        self.indice_basename = self.build_indices()
        self.aligned_sam_file = self.align_reads_to_contigs()
        self.contig_abundance_file = self.convert_and_sort_aligned_sam_file()
        
        
        
    def build_indices(self):
        
        indice_basename = f"{self.output_directory}/contig_indices"
        bowtie_indice_build_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/prebinning', 'bowtie2-build', '--threads', self.threads, self.contig_file_path, indice_basename]
        run_and_log_a_subprocess(self.output_directory, bowtie_indice_build_args, "bowtie_build_contig_indices")
        return indice_basename
        
    def align_reads_to_contigs(self):
        
        sam_output_file = f"{self.output_directory}/aligned_reads_to_contigs.sam"
        
        bowtie2_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/prebinning', 'bowtie2', '-x', self.indice_basename, '-1', self.read_fwd_path, 
                        '-2', self.read_rev_path, '-p', self.threads, '--very-sensitive', 
                        '--no-unal', '-S', sam_output_file]
        
        run_and_log_a_subprocess(self.output_directory, bowtie2_args, 'bowtie_sam_log')
        
        return sam_output_file
    
        
    def convert_and_sort_aligned_sam_file(self):
        bam_output_file = f"{self.output_directory}/aligned_reads_to_contigs.bam"
        sorted_bam_output_file = f"{self.output_directory}/aligned_reads_to_contigs.bam"
        samtools_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/prebinning', 'samtools', 'view', '-@', self.threads, '-Sb', self.aligned_sam_file, '-o', bam_output_file] 
        samtools_sort_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/prebinning', 'samtools', 'sort', '-@', self.threads, '-O', 'bam', '-o', bam_output_file, sorted_bam_output_file]
        samtools_index_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/prebinning', 'samtools', 'index', '-@', self.threads, sorted_bam_output_file]
        
        run_and_log_a_subprocess(self.output_directory, samtools_args, 'samtools_conversion')
        run_and_log_a_subprocess(self.output_directory, samtools_sort_args, 'samtools_sort')
        run_and_log_a_subprocess(self.output_directory, samtools_index_args, 'samtools_index')
        
        return bam_output_file
    
        
class BinSet:
     
    def __init__(self, directory_of_bins):
        self.directory_of_bins = directory_of_bins
         
         

class Binner:
    contigs_path = ""
    fwd_read_path = ""
    rev_read_path = ""
    abundance_information_path = ""
    log_directory_path = ""
    threads = "1"
    min_contig_length = "2000"
    


    @classmethod
    def add_read_contig_and_abundance_paths_to_base_binner_class(cls, contigs_path, fwd_read_path, rev_read_path, abundance_file_path):
        
        cls.contigs_path = contigs_path
        cls.fwd_read_path = fwd_read_path
        cls.rev_read_path = rev_read_path
        cls.abundance_information_path = abundance_file_path

    @classmethod
    def add_or_change_threads(cls, new_threads):
        cls.threads = new_threads

    @classmethod

    def add_or_change_log_directory(cls, log_directory_path):
        
        cls.log_directory_path = log_directory_path

    @classmethod
    def add_or_change_min_contig_length(cls, min_contig_length):
        
        cls.min_contig_length = min_contig_length
        
        
    @classmethod
    def calculate_read_depth(cls, depth_output_path):
        
        get_contig_depth_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/CONCOCTMetabat2MaxBin2SemiBin2', 'jgi_summarize_bam_contig_depths', '--outputDepth', depth_output_path, cls.abundance_information_path]
        run_and_log_a_subprocess(cls.log_directory_path, get_contig_depth_args, "metabat2_contig_read_depth_gen")
        cls.read_depths_path = depth_output_path
       
        
    def run_maxbin2(self, output_directory):
        
        maxbin2_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/CONCOCTMetabat2MaxBin2SemiBin2', 'run_MaxBin.pl', '-contig', self.contigs_path,
                        '-thread', self.threads, '-abund', self.abundance_information_path, '-out', f"{output_directory}/sample_result_"]
        
        run_and_log_a_subprocess(output_directory, maxbin2_args, "maxbin2_binning")

        

    def run_metabat2(self, output_directory):
        
        metabat2_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/CONCOCTMetabat2MaxBin2SemiBin2','metabat2', '-m', self.min_contig_length, '-t', self.threads, '--unbinned',
                        '--seed', '0', '-i', self.contigs_path, '-a', self.read_depths_path, '-o', output_directory]
        
        run_and_log_a_subprocess(output_directory, metabat2_args, "metabat2_binning")
        
        
    def run_semibin2(self, output_directory):
        
        semibin_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/CONCOCTMetabat2MaxBin2SemiBin2', 'SemiBin', 'single_easy_bin', '-i', self.contigs_path, '-b', self.abundance_information_path, '-o', output_directory, '-t', self.threads,
                        '--write-pre-reclustering-bins', '--training-type', 'self']
        
        run_and_log_a_subprocess(output_directory, semibin_args, "semibin_binning")

    def run_concoct(self, output_directory):
        contig_bed_file_path = f"{output_directory}/contigs_10k.bed"
        contig_10k_file_path = f"{output_directory}/contigs_10k.fa"
        step_1_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/CONCOCTMetabat2MaxBin2SemiBin2', 'cut_up_fasta.py', self.contigs_path, '-c', '10000', '-o', '0', '--merge_last', '-b', contig_bed_file_path]
        run_and_log_a_subprocess(output_directory, step_1_args, "concoct_cutup_fasta", alternate_stdout_path=contig_10k_file_path)
        
        coverage_table_path = f"{output_directory}/coverage_table.tsv"
        step_2_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/CONCOCTMetabat2MaxBin2SemiBin2', 'concoct_coverage_table.py', contig_bed_file_path, self.abundance_information_path]
        run_and_log_a_subprocess(output_directory, step_2_args, "concoct_gen_cov_table", alternate_stdout_path=coverage_table_path)
        
        step_3_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/CONCOCTMetabat2MaxBin2SemiBin2', 'concoct', '--composition_file', self.contigs_path, '--length_threshold', self.min_contig_length, '--coverage_file', coverage_table_path, '-b', output_directory]
        run_and_log_a_subprocess(output_directory, step_3_args, "concoct_step_3_run_concoct")
        
        clustered_file_path = f"{output_directory}/clustering_gt{self.min_contig_length}.csv"
        clustering_merged_path = f"{output_directory}/clustering_merged.csv"
        step_4_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/CONCOCTMetabat2MaxBin2SemiBin2', 'merge_cutup_clustering.py', clustered_file_path]
        run_and_log_a_subprocess(output_directory, step_4_args, "concoct_gen_cov_table", alternate_stdout_path=clustering_merged_path)
        final_concoct_bins_path = f"{output_directory}/fasta_bins/"
        os.mkdir(final_concoct_bins_path)
        
        step_1_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/CONCOCTMetabat2MaxBin2SemiBin2', 'extract_fasta_bins.py', self.contigs_path, clustering_merged_path, '--output_path', final_concoct_bins_path]
        run_and_log_a_subprocess(output_directory, step_3_args, "concoct_final_step_extract_fasta_bins")
        
        
    def run_vamb(self, output_directory):
            # vamb --outdir path/to/outdir --fasta /path/to/catalogue.fna.gz --bamfiles /path/to/bam/*.bam -o C
        vamb_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/Vamb4', 'vamb', '--outdir', output_directory, '--fasta', self.contigs_path, '--bamfiles', self.abundance_information_path, '-m', self.min_contig_length,
                         '-p', self.threads]
        run_and_log_a_subprocess(output_directory, vamb_args, "vamb_binning")



def setup_binning(args):
    contig_path = generate_contigs(args)

    contig_abundance_gen = ContigAbundances(output_directory=args.output_directory, contig_file_path=contig_path, read_fwd_path=args.forward_reads, read_rev_path=args.reverse_reads, threads=args.threads)

    Binner.add_read_contig_and_abundance_paths_to_base_binner_class(contigs_path=contig_path, fwd_read_path=args.forward_reads, rev_read_path=args.reverse_reads, abundance_file_path=contig_abundance_gen.contig_abundance_file)
    log_directory = f"{args.output_directory}/log_directory/"
    os.mkdir(log_directory)
    Binner.add_or_change_log_directory(log_directory)
    the_binner = Binner()
    the_binner.calculate_read_depth(f"{args.output_directory}/read_depths.txt")
    
    return contig_abundance_gen,the_binner

def generate_contigs(args):
    if args.contig_path:
        return args.contig_path
    else:
        if not args.using_scaffolds:
            using_scaffolds = False
        else:
            using_scaffolds = True
    
        assembler = MetaSpadesAssemblyRunner(args.output_directory, using_scaffolds, args.forward_reads, args.reverse_reads, args.threads)
        contig_path = assembler.check_if_assembly_step_has_already_been_completed()
        
        if contig_path == False:
            
            if not args.custom_kmer_lengths:
            
                contig_path = assembler.run_assembly()
            
            else:
                
                contig_path = assembler.run_assembly_with_custom_kmer_lengths(args.custom_kmer_lengths)
                
    return contig_path


def run_binning(output_directory, the_binner, binner_option_list):
    
    bin_methods_dict = {'Semibin2' : the_binner.run_semibin2(f"{output_directory}/Semibin2/"), 'Metabat2' : the_binner.run_metabat2(f"{output_directory}/Metabat2/"),
                        'Maxbin2' : the_binner.run_maxbin2(f"{output_directory}/Maxbin2/"), "Vamb" : the_binner.run_vamb(f"{output_directory}/Vamb/"), "CONCOCT" : the_binner.run_concoct(f"{output_directory}/CONCOCT/")}
    
    for binner in binner_option_list:
        
        try:
            
            bin_methods_dict[binner]()
        
        except KeyError:
            
            print(f"Could not find individual binner chosen: {binner} - please check your binner options. Exiting.")
            exit()