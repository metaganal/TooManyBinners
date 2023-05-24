import subprocess
import os
from utils import run_and_log_a_subprocess

class ContigAbundances:
    
    
    def __init__(self, output_directory, contig_file_path, read_fwd_path, read_rev_path, threads):
        
        self.output_directory = output_directory
        self.contig_file_path = contig_file_path
        self.read_fwd_path = read_fwd_path
        self.red_rev_path = read_rev_path
        self.threads = threads
        
        self.indice_basename = self.build_indices()
        self.aligned_sam_file = self.align_reads_to_contigs()
        self.contig_abundance_file = self.convert_and_sort_aligned_sam_file()
        
        
        
    def build_indices(self):
        
        indice_basename = f"{self.output_directory}/contig_indices"
        bowtie_indice_build_args = ['conda', 'run', '--prefix', '/opt/miniconda3/envs/ContigAbundDepthsEnv', 'bowtie2-build', '--threads', self.threads, self.contig_file_path, indice_basename]
        run_and_log_a_subprocess(self.output_directory, bowtie_indice_build_args, "bowtie_build_contig_indices")
        return indice_basename
        
    def align_reads_to_contigs(self):
        
        sam_output_file = f"{self.output_directory}/aligned_reads_to_contigs.sam"
        
        bowtie2_args = ['conda', 'run', '--prefix', '/opt/miniconda3/envs/ContigAbundDepthsEnv', 'bowtie2', '-x', self.indice_basename, '-1', self.read_fwd_path, 
                        '-2', self.read_fwd_path, '-p', self.threads, '--very-sensitive', 
                        '--no-unal', '-S', sam_output_file]
        
        run_and_log_a_subprocess(self.output_directory, bowtie2_args, 'bowtie_sam_log')
        
        return sam_output_file
    
        
    def convert_and_sort_aligned_sam_file(self):
        bam_output_file = f"{self.output_directory}/aligned_reads_to_contigs.bam"
        sorted_bam_output_file = f"{self.output_directory}/aligned_reads_to_contigs.bam"
        samtools_args = ['conda', 'run', '--prefix', '/opt/miniconda3/envs/ContigAbundDepthsEnv', 'samtools', 'view', '-@', self.threads, '-Sb', self.aligned_sam_file, '-o', bam_output_file] 
        samtools_sort_args = ['conda', 'run', '--prefix', '/opt/miniconda3/envs/ContigAbundDepthsEnv', 'samtools', 'sort', '-@', self.threads, '-O', 'bam', '-o', bam_output_file, sorted_bam_output_file]
        samtools_index_args = ['conda', 'run', '--prefix', '/opt/miniconda3/envs/ContigAbundDepthsEnv', 'samtools', 'index', '-@', self.threads, sorted_bam_output_file]
        
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
    threads = 1
    min_contig_length = 1000
    


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
        
        get_contig_depth_args = ['conda', 'run' '--prefix', '/opt/miniconda3/envs/Metabat2' 'jgi_summarize_bam_contig_depths', '--outputDepth', depth_output_path, cls.abundance_information_path]
        run_and_log_a_subprocess(cls.log_directory_path, get_contig_depth_args, "metabat2_contig_read_depth_gen")
        cls.read_depths_path = depth_output_path
       
        
    def run_maxbin2(self, output_directory):
        
        maxbin2_args = ['conda', 'run' '--prefix', '/opt/miniconda3/envs/Maxbin2', 'run_MaxBin.pl', '-contig', self.contigs_path,
                        '-threads', self.threads, '-abund', self.abundance_information_path, '-out', output_directory]
        
        run_and_log_a_subprocess(output_directory, maxbin2_args, "maxbin2_binning")

        


    def run_metabat2(self, output_directory):
        
        metabat2_args = ['conda', 'run' '--prefix', '/opt/miniconda3/envs/Metabat2','metabat2', '-m', self.min_contig_length, '-t', self.threads, '--unbinned',
                        '--seed', '0', '-i', self.contigs_path, '-a', self.read_depths_path, '-o', output_directory]
        
        run_and_log_a_subprocess(output_directory, metabat2_args, "metabat2_binning")
        
        
    def run_semibin2(self, output_directory):
        
        semibin_args = ['conda', 'run', '--prefix', '/opt/miniconda3/envs/Semibin2', 'SemiBin', 'single_easy_bin', '-i', self.contigs_path, '-b', self.read_depths_path, '-o', output_directory, '-t', self.threads,
                        '--write-pre-reclustering-bins', '--training-type', 'self']
        
        run_and_log_a_subprocess(output_directory, semibin_args, "semibin_binning")





def setup_binning(args):
    contig_abundance_gen = ContigAbundances(output_directory=args.output_directory, contig_file_path=args.contig_path, read_fwd_path=args.forward_reads, read_rev_path=args.reverse_reads, threads=args.threads)

# Then running of individual binners and verify each option
    Binner.add_read_contig_and_abundance_paths_to_base_binner_class(contigs_path=args.contig_path, fwd_read_path=args.forward_reads, rev_read_path=args.reverse_reads, 
                                                                    abundance_file_path=contig_abundance_gen.contig_abundance_file)
    
    the_binner = Binner()
    the_binner.calculate_read_depth(f"{args.output_directory}/read_depths.txt")
    
    return contig_abundance_gen,the_binner


def run_binning(output_directory, the_binner, binner_option_list):
    
    bin_methods_dict = {'Semibin2' : the_binner.run_semibin2(f"{output_directory}/Semibin2/"), 'Metabat2' : the_binner.run_metabat2(f"{output_directory}/Metabat2/"), 
                        'Maxbin2' : the_binner.run_maxbin2(f"{output_directory}/Maxbin2/")}
    
    for binner in binner_option_list:
        
        try:
            
            bin_methods_dict[binner]()
        
        except KeyError:
            
            print(f"Could not find individual binner chosen: {binner} - please check your binner options. Exiting.")
            exit()