import shutil
import subprocess
import os
from .utils import run_and_log_a_subprocess

class MetaSpadesAssemblyRunner:
    
    def __init__(self, output_directory, using_scaffolds, read_fwd_path, read_rev_path, threads, sample_name, log_directory, memory):
        self.output_directory = output_directory
        self.using_scaffolds = using_scaffolds
        self.read_fwd_path = read_fwd_path
        self.read_rev_path = read_rev_path
        self.threads = threads
        self.sample_name = sample_name
        self.log_directory = log_directory
        self.memory = memory

    
    def run_assembly(self):
        
        metaspades_result_directory = f"{self.output_directory}/metaspades_assembly_directory/"
        metaspades_assembly_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/prebinning', 'spades.py', '-1', self.read_fwd_path, '-2', self.read_rev_path, '-t', self.threads, '--memory', self.memory, '--meta', '-o', metaspades_result_directory]

        run_and_log_a_subprocess(self.log_directory, metaspades_assembly_args, "metaspades_assembly_auto_kmer")

        if self.using_scaffolds == True:

            print("Using metaspades generated scaffolds as contigs")
            os.rename(f"{metaspades_result_directory}/scaffolds.fasta", f"{metaspades_result_directory}/{self.sample_name}_scaffolds.fasta")
            return f"{metaspades_result_directory}/{self.sample_name}_scaffolds.fasta"

        else:
            os.rename(f"{metaspades_result_directory}/contigs.fasta", f"{metaspades_result_directory}/{self.sample_name}_contigs.fasta")
            return f"{metaspades_result_directory}/{self.sample_name}_contigs.fasta"


    def run_assembly_with_custom_kmer_lengths(self, custom_kmer_lengths):
        metaspades_result_directory = f"{self.output_directory}/metaspades_assembly_directory/"
        metaspades_assembly_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/prebinning', 'spades.py', '-1', self.read_fwd_path, '-2', self.read_rev_path, '-t', self.threads, '-k', custom_kmer_lengths, '--memory', self.memory, '--meta', '-o', metaspades_result_directory]

        run_and_log_a_subprocess(self.log_directory, metaspades_assembly_args, "metaspades_assembly_auto_kmer")


        if self.using_scaffolds == True:

            print("Using metaspades generated scaffolds as contigs")
            return f"{metaspades_result_directory}/{self.sample_name}_scaffolds.fasta"

        else:

            return f"{metaspades_result_directory}/contigs.fasta"



    def check_if_assembly_step_has_already_been_completed(self):
        metaspades_result_directory = f"{self.output_directory}/metaspades_assembly_directory/"
        if os.path.isdir(metaspades_result_directory):
            if not os.path.exists(f"{metaspades_result_directory}/{self.sample_name}_contigs.fasta") and not os.path.exists(f"{metaspades_result_directory}/{self.sample_name}_scaffolds.fasta"):
                print("Found metaspades directory but on contigs or scaffolds. Error, exiting.")
                exit()
            if self.using_scaffolds == True:
                return f"{metaspades_result_directory}/{self.sample_name}_scaffolds.fasta"
            else:
                return f"{metaspades_result_directory}/{self.sample_name}_contigs.fasta"
        
        return False
        
class ContigAbundances:
    
    
    def __init__(self, output_directory, contig_file_path, read_fwd_path, read_rev_path, threads, log_directory):
        
        self.output_directory = output_directory
        self.log_directory = log_directory
        self.contig_file_path = contig_file_path
        self.read_fwd_path = read_fwd_path
        self.read_rev_path = read_rev_path
        self.threads = threads
        
        self.indice_basename = self.build_indices()
        self.aligned_sam_file = self.align_reads_to_contigs()
        self.contig_abundance_file = self.convert_and_sort_aligned_sam_file()
        
        
        
    def build_indices(self):
        
        indice_basename = f"{self.output_directory}/contig_indices"

        if len([file for file in os.listdir(self.output_directory) if ".bt2" in file]) == 6:
            return indice_basename
        bowtie_indice_build_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/prebinning', 'bowtie2-build', '--threads', self.threads, self.contig_file_path, indice_basename]
        run_and_log_a_subprocess(self.log_directory, bowtie_indice_build_args, "bowtie_build_contig_indices")

        return indice_basename
        
    def align_reads_to_contigs(self):
        sam_output_file = f"{self.output_directory}/aligned_reads_to_contigs.sam"

        if os.path.exists(sam_output_file):
            return sam_output_file

        bowtie2_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/prebinning', 'bowtie2', '-x', self.indice_basename, '-1', self.read_fwd_path, 
                        '-2', self.read_rev_path, '-p', self.threads, '--very-sensitive', 
                        '--no-unal', '-S', sam_output_file]
        
        run_and_log_a_subprocess(self.log_directory, bowtie2_args, 'bowtie_sam_log')
        
        return sam_output_file
    
        
    def convert_and_sort_aligned_sam_file(self):
        
        bam_output_file = f"{self.output_directory}/aligned_reads_to_contigs.bam"

        sorted_bam_output_file = f"{self.output_directory}/sorted_reads_to_contigs.bam"
        
        if os.path.exists(sorted_bam_output_file):
            return sorted_bam_output_file

        
        samtools_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/prebinning', 'samtools', 'view', '-@', self.threads, '-Sb', self.aligned_sam_file, '-o', bam_output_file] 
        samtools_sort_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/prebinning', 'samtools', 'sort', '-@', self.threads, '-O', 'bam', '-o', sorted_bam_output_file, bam_output_file]
        samtools_index_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/prebinning', 'samtools', 'index', '-@', self.threads, sorted_bam_output_file]
        
        run_and_log_a_subprocess(self.log_directory, samtools_args, 'samtools_conversion')
        run_and_log_a_subprocess(self.log_directory, samtools_sort_args, 'samtools_sort')
        run_and_log_a_subprocess(self.log_directory, samtools_index_args, 'samtools_index')
        
        return sorted_bam_output_file
    
        
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
    sample_name = ""
    final_bins_directory = ""



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
    def add_or_change_final_bins_directory(cls, final_bins_directory):
        
        cls.final_bins_directory = final_bins_directory
        

    @classmethod
    def add_or_change_min_contig_length(cls, min_contig_length):
        
        cls.min_contig_length = min_contig_length
        
    @classmethod
    def add_or_change_sample_name(cls, sample_name):
        cls.sample_name = sample_name

    @classmethod
    def calculate_read_depth(cls, depth_output_path):

        if os.path.exists(depth_output_path):
            cls.read_depths_path = depth_output_path
            return

        get_contig_depth_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/CONCOCTMetabat2MaxBin2SemiBin2', 'jgi_summarize_bam_contig_depths', '--outputDepth', depth_output_path, cls.abundance_information_path]
        run_and_log_a_subprocess(cls.log_directory_path, get_contig_depth_args, "metabat2_contig_read_depth_gen")
        cls.read_depths_path = depth_output_path
       
    def move_bin_results_to_main_bin_dir(self, bin_directory, binner_name):
        for result in os.listdir(bin_directory):
            result_path = f"{bin_directory}/{result}"
            
            if os.path.isdir(result_path):
            
                self.move_bin_results_to_main_bin_dir(result_path, binner_name)
        
            if ".fa" or ".fasta" in result_path:
                shutil.copy(result_path, f"{self.final_bins_directory}/{binner_name}_{result}")



    def run_maxbin2(self, output_directory):

        if os.path.exists(f"{output_directory}"):
            return
        os.mkdir(output_directory)
        maxbin2_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/CONCOCTMetabat2MaxBin2SemiBin2', 'run_MaxBin.pl', '-contig', self.contigs_path, '-min_contig_length', self.min_contig_length,
                        '-thread', self.threads, '-abund', self.read_depths_path, '-out', f"{output_directory}/sample_result_"]
        
        run_and_log_a_subprocess(self.log_directory_path, maxbin2_args, "maxbin2_binning")
        self.move_bin_results_to_main_bin_dir(output_directory, "maxbin2")

        

    def run_metabat2(self, output_directory):
        if os.path.exists(f"{output_directory}"):
            print("Already found output metabat2 directory.")
            return
        
        metabat2_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/CONCOCTMetabat2MaxBin2SemiBin2','metabat2', '-m', self.min_contig_length, '-t', self.threads, '--unbinned',
                        '--seed', '0', '-i', self.contigs_path, '-a', self.read_depths_path, '-o', output_directory]
        
        run_and_log_a_subprocess(self.log_directory_path, metabat2_args, "metabat2_binning")
        self.move_bin_results_to_main_bin_dir(output_directory, "metabat2")

        
        
    def run_semibin2(self, output_directory):

        if os.path.exists(f"{output_directory}/output_recluster_bins/"):
            return

        semibin_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/CONCOCTMetabat2MaxBin2SemiBin2', 'SemiBin', 'single_easy_bin', '-i', self.contigs_path, '-b', self.abundance_information_path, '-o', output_directory, '-t', self.threads,
                        '--write-pre-reclustering-bins', '--training-type', 'self']
        
        run_and_log_a_subprocess(self.log_directory_path, semibin_args, "semibin_binning")
        self.move_bin_results_to_main_bin_dir(f"{output_directory}/output_recluster_bins/", "semibin2")


    def run_concoct(self, output_directory):
        if os.path.exists(f"{output_directory}"):
            return
        contig_bed_file_path = f"{output_directory}/contigs_10k.bed"
        contig_10k_file_path = f"{output_directory}/contigs_10k.fa"

        os.mkdir(output_directory)
        step_1_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/CONCOCTMetabat2MaxBin2SemiBin2', 'cut_up_fasta.py', self.contigs_path, '-c', '10000', '-o', '0', '--merge_last', '-b', contig_bed_file_path]
        run_and_log_a_subprocess(self.log_directory_path, step_1_args, "concoct_cutup_fasta", alternate_stdout_path=contig_10k_file_path)
        
        coverage_table_path = f"{output_directory}/coverage_table.tsv"
        step_2_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/CONCOCTMetabat2MaxBin2SemiBin2', 'concoct_coverage_table.py', contig_bed_file_path, self.abundance_information_path]
        run_and_log_a_subprocess(self.log_directory_path, step_2_args, "concoct_gen_cov_table", alternate_stdout_path=coverage_table_path)
        
        step_3_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/CONCOCTMetabat2MaxBin2SemiBin2', 'concoct', '--composition_file', contig_10k_file_path, '--length_threshold', self.min_contig_length, '--coverage_file', coverage_table_path, '-b', output_directory]
        run_and_log_a_subprocess(self.log_directory_path, step_3_args, "concoct_step_3_run_concoct")

        
        clustered_file_path = f"{output_directory}/clustering_gt{self.min_contig_length}.csv"
        clustering_merged_path = f"{output_directory}/clustering_merged.csv"
        step_4_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/CONCOCTMetabat2MaxBin2SemiBin2', 'merge_cutup_clustering.py', clustered_file_path]

        run_and_log_a_subprocess(self.log_directory_path, step_4_args, "concoct_cutup_clustering", alternate_stdout_path=clustering_merged_path)
        final_concoct_bins_path = f"{output_directory}/fasta_bins"
        os.mkdir(final_concoct_bins_path)
        
        step_5_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/CONCOCTMetabat2MaxBin2SemiBin2', 'extract_fasta_bins.py', self.contigs_path, clustering_merged_path, '--output_path', final_concoct_bins_path]
        run_and_log_a_subprocess(self.log_directory_path, step_5_args, "concoct_final_step_extract_fasta_bins")
        self.move_bin_results_to_main_bin_dir(final_concoct_bins_path, "concoct")

        
        
    def run_vamb(self, output_directory):
        if os.path.exists(output_directory):
            return 
            # vamb --outdir path/to/outdir --fasta /path/to/catalogue.fna.gz --bamfiles /path/to/bam/*.bam -o C

        os.mkdir(output_directory)
        vamb_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/Vamb4', 'vamb', '--outdir', f"{output_directory}/results/", '--fasta', self.contigs_path, '--bamfiles', self.abundance_information_path, '-m', self.min_contig_length,

                         '-p', self.threads]
        run_and_log_a_subprocess(self.log_directory_path, vamb_args, "vamb_binning")
        # note this filters out any genomes with less than 100k bases in size!
        create_fasta_args = ['mamba', 'run', '--prefix', '/opt/mamba/envs/Vamb4', 'python3', '/opt/create_fasta.py', self.contigs_path, f"{output_directory}/results/vae_clusters.tsv", '100000', f"{output_directory}/bins"]
        run_and_log_a_subprocess(self.log_directory_path, create_fasta_args, "vamb_binning_step2")
        self.move_bin_results_to_main_bin_dir(f"{output_directory}/bins", "vamb")
        




def setup_binning(args, sample_name):
    log_directory = f"{args.output_directory}/log_directory/"

    if not os.path.exists(log_directory):
        os.mkdir(log_directory)
    
    contig_path = generate_contigs(args, sample_name, log_directory)


    contig_abundance_gen = ContigAbundances(output_directory=args.output_directory, contig_file_path=contig_path, read_fwd_path=args.forward_reads, read_rev_path=args.reverse_reads, threads=args.threads, log_directory=log_directory)

    Binner.add_read_contig_and_abundance_paths_to_base_binner_class(contigs_path=contig_path, fwd_read_path=args.forward_reads, rev_read_path=args.reverse_reads, abundance_file_path=contig_abundance_gen.contig_abundance_file)
    Binner.add_or_change_sample_name = sample_name
    if args.minimum_contig_length:
        Binner.add_or_change_min_contig_length(args.min_contig_length)
    Binner.add_or_change_log_directory(log_directory)
    the_binner = Binner()
    the_binner.calculate_read_depth(f"{args.output_directory}/{sample_name}_read_depths.txt")
    
    return contig_abundance_gen,the_binner



def generate_contigs(args, sample_name, log_directory):
    if args.contig_path:
        return args.contig_path
    else:
        if not args.using_scaffolds:
            using_scaffolds = False
        else:
            using_scaffolds = True
    
        assembler = MetaSpadesAssemblyRunner(args.output_directory, using_scaffolds, args.forward_reads, args.reverse_reads, args.threads, sample_name, log_directory, args.memory)
        contig_path = assembler.check_if_assembly_step_has_already_been_completed()
        
        if contig_path == False:
            
            if not args.custom_kmer_lengths:
            
                contig_path = assembler.run_assembly()
            
            else:
                
                contig_path = assembler.run_assembly_with_custom_kmer_lengths(args.custom_kmer_lengths)
                
    return contig_path



def run_binning(output_directory, the_binner, binner_option_list):
    the_binner.add_or_change_final_bins_directory(f"{output_directory}/final_bins/")
    if not os.path.isdir(f"{output_directory}/final_bins/"):
        os.mkdir(f"{output_directory}/final_bins/")
    bin_methods_dict = {'Semibin2' : the_binner.run_semibin2(f"{output_directory}/Semibin2/"), 'Metabat2' : the_binner.run_metabat2(f"{output_directory}/Metabat2/"),
                        'Maxbin2' : the_binner.run_maxbin2(f"{output_directory}/Maxbin2/"), "Vamb" : the_binner.run_vamb(f"{output_directory}/Vamb/"), "CONCOCT" : the_binner.run_concoct(f"{output_directory}/CONCOCT/")}
    
    for binner in binner_option_list:
        
        try:
            
            bin_methods_dict[binner]()
        
        except KeyError:
            
            print(f"Could not find individual binner chosen: {binner} - please check your binner options. Exiting.")
            exit()



