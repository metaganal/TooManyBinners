import os
import subprocess

def run_and_log_a_subprocess(process_log_dir, process_args, name_of_sample=None, alternate_stdout_path=False):

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

