import os

def make_input_list(input_dir):
    
    input_files = []

    runs = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
    for run in runs:
        run_dir = os.path.join(input_dir, run+'/AOD/')

        subruns = [d for d in os.listdir(run_dir) if os.path.isdir(os.path.join(run_dir, d))]

        for subrun in subruns:
            input_files.append(os.path.join(run_dir, subrun, 'AO2D.root'))

    with open(f'{input_dir}/input_files.txt', 'w') as f:
        for file in input_files:
            f.write(os.path.join(input_dir, file) + '\n')

run_dir = '/home/kaareendrup/analysis/localMuons_table/input_data/DQ_data'
make_input_list(run_dir)