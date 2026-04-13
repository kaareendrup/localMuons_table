
import os
import json
import math

def generate_input_multi(
    data_dir, 
    config_dir,
    config_name, 
    outputdirector_name,
    command_base,
    n_jobs, 
    n_files,
    MC = True
):

    with open(data_dir+'/input_data.txt', 'r') as f:
        lines = f.readlines()

    if n_files != -1:
        assert n_files <= len(lines), "Not enough files!"
    else:
        n_files = len(lines)
    files_per_job = math.ceil(n_files/n_jobs)

    print('Jobs:', n_jobs)
    print('Total files:', n_files)
    print('Files per job:', files_per_job)

    if MC:
        config_file = f"{config_dir}/configuration_dqEfficiency_withAssoc_{config_name}.json"
    else:
        config_file = f"{config_dir}/configuration_dqTableReader_withAssoc_{config_name}.json"
    output_director_file = f"{config_dir}/{outputdirector_name}.json"
    command_base = command_base.replace("--configuration json://", f"--configuration json://{config_file}")

    # Set outdir
    out_dir = f'{data_dir}/multi/{config_name}'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    command_file = f'{out_dir}/commands.txt'
    run_file = f'{out_dir}/run.sh'

    with open(command_file, 'w') as cf, open(run_file, 'w') as rf:

        rf.write('#!/bin/bash'+'\n')

        for j in range(n_jobs):

            out_file = f'{out_dir}/input_data_multi_{j}.txt'

            # Write ouputdirector for different AOD name
            with open(output_director_file, 'r') as jf:
                data = json.load(jf)

            data['OutputDirector']['resfile'] += str(j)

            output_director_file_j = f'{out_dir}/outputDirector_withAssoc_{config_name}_{j}.json'
            with open(output_director_file_j, 'w') as jf:
                json.dump(data, jf, indent=4)

            # Write unique command
            out_command = f"{command_base} --aod-file @{out_file} --aod-writer-json {output_director_file_j}"
            cf.write(out_command+'\n')

            # Write commands to bash file
            rf.write('\n')
            # rf.write(f'cd {j} || exit\n')
            rf.write(out_command+'\n')
            # rf.write('cd ..\n')

            # Write output to file
            with open(f'{out_file}', 'w') as f:
                for i in range(files_per_job):

                    f.write(lines[j*files_per_job+i])

                    if j*files_per_job+i >= n_files-1:
                        break

            print(f'Wrote {i+1} lines to {out_file}')
        
            if j*files_per_job+i >= n_files-1:
                print('All files assigned!')
                break

    print("To run, run:")
    print(f"bash {out_dir}/run.sh")

config_dir = "/home/kaareendrup/analysis/localMuons_table/config"
input_dir = "/home/kaareendrup/analysis/input_data"

############## ACTUAL CONFIGURABLES ##############
# data_name = "c3"
data_name = "f4d"
# data_name = "DQ_data"

muon_type = "standalone"
# muon_type = "global"

data_type = "reco"
# data_type = "gen"

# n_files = -1
n_files = 740
# n_files = 400
# n_files = 200
# n_files = 16
# n_files = 1

# n_files_per_job = 8
n_files_per_job = 4

############## AUTOMATIC ##############
data_dir = f"{input_dir}/{data_name}_{muon_type}"
config_name = f"{muon_type}_{data_type}"

if data_name == "DQ_data":
    MC = False
    outputdirector_name = "outputDirector_data"
    command_base = "o2-analysis-dq-table-reader-with-assoc -b --configuration json:// | o2-analysis-dq-model-converter-event-extended -b --configuration json://"
else:
    MC = True
    outputdirector_name = "outputDirector_MC"
    command_base = "o2-analysis-dq-efficiency-with-assoc -b --configuration json://"

assert n_files % n_files_per_job == 0, "N files must be divisible by files per job!"
n_jobs = int(n_files / n_files_per_job)

generate_input_multi(data_dir, config_dir, config_name, outputdirector_name,command_base, n_jobs, n_files, MC)
