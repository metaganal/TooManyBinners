import os
import subprocess
from classes import Binner
from classes import ContigAbundances

def run_and_log_a_subprocess(process_log_dir, process_args, name_of_process, name_of_sample=None, alternate_stdout_path=False):

    try:
        
        os.mkdir(process_log_dir)
    
    except FileExistsError:
        
        pass
    


    if name_of_sample == None:
        
        current_process_num = len(os.listdir(process_log_dir)) + 1

        process_log_stdout = f"{process_log_dir}/process_{current_process_num}_stdoutput.txt"
        process_log_stderr = f"{process_log_dir}/process_{current_process_num}_stderr.txt"
    
    else:
        
        process_log_stdout = f"{process_log_dir}/process_{name_of_sample}_stdoutput.txt"
        process_log_stderr = f"{process_log_dir}/process_{name_of_sample}_stderr.txt"
    
    
    result = subprocess.run(process_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    if alternate_stdout_path != False:
        process_log_stdout = alternate_stdout_path

    with open(process_log_stdout, "w") as f:
            
        f.write(result.stdout)
    
    with open(process_log_stderr, "w") as f:
        
        f.write(result.stderr)

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