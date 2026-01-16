import pandas as pd

def generate_download_script(jobs_file, n_jobs, file_dir):

    # Read the formatted jobs CSV file to get ALICE grid directories
    jobs_data = pd.read_csv(jobs_file, header=None, names=[
        'id', 
        'dir', 
        'throughput', 
        'total', 
        'done', 
        'active', 
        'wait', 
        'error', 
        'status'
    ]).head(n_jobs)

    with open('download.sh', 'w') as script_file:

        # Write header for the bash script
        script_file.write('#!/bin/bash\n\n')
        script_file.write('echo "Downloading input data files..."\n')

        # Write download commands for each job
        for _, row in jobs_data.iterrows():
            source = row['dir']
            output_dir = file_dir + source.split('/')[-2]
            script_file.write(f'alien.py cp {source}*AO2D.root file:{output_dir}\n')

    n_jobs = len(jobs_data)

    print(f'Download script "download.sh" generated for {n_jobs} jobs.')
    print(f'With an average 100 MB per job, this amounts to approximately {n_jobs * 100} MB of data, or {(n_jobs * 100) / 1024:.2f} GB.')
    print(f'At a speed of .5 MB/s, this will take around {(n_jobs * 100) / 0.5 / 60:.2f} minutes to download.')

MC_name = 'DQ'
n_jobs = int(1e12)
# n_jobs = 100
# MC_name = 'HF'
# n_jobs = 10
# MC_name = 'genpurp'

jobs_file = f'input_multi/{MC_name}/jobs_formatted.csv'
# output_dir = '/home/kaareendrup/analysis/localMuons_table/input_data/' + MC_name + '/'
output_dir = '/media/kaareendrup/ec65dbb9-11dd-4abf-ae3e-f029466fe958/analysis/input_data/' + MC_name + '/'

generate_download_script(jobs_file, n_jobs, output_dir)