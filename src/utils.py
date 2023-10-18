import os
import subprocess

def run_and_log_a_subprocess(process_log_dir, process_args, process_name=None, alternate_stdout_path=False):

    try:
        
        os.mkdir(process_log_dir)
    
    except FileExistsError:
        
        pass
    


    if process_name == None:
        
        current_process_num = len(os.listdir(process_log_dir)) + 1

        process_log_stdout = f"{process_log_dir}/process_{current_process_num}_stdoutput.txt"
        process_log_stderr = f"{process_log_dir}/process_{current_process_num}_stderr.txt"
    
    else:
        
        process_log_stdout = f"{process_log_dir}/process_{process_name}_stdoutput.txt"
        process_log_stderr = f"{process_log_dir}/process_{process_name}_stderr.txt"
    
    
    result = subprocess.run(process_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    if alternate_stdout_path != False:
        process_log_stdout = alternate_stdout_path

    with open(process_log_stdout, "w") as f:
            
        f.write(result.stdout)
    
    with open(process_log_stderr, "w") as f:
        
        f.write(result.stderr)

    write_to_checkpoint_file(process_log_dir, process_name)


def write_to_checkpoint_file(process_log_dir, checkpoint_string):
    with open(f"{process_log_dir}/checkpoints.txt", "w") as f:
        f.write(f"{checkpoint_string},")


def get_name_of_sample(rev_path, fwd_path):
    sample_name = ""
    for i in range(0, len(rev_path), 1):
        if rev_path[i] != fwd_path[i]:
            break
        sample_name = sample_name + rev_path[i]
        
    return sample_name

