import os

def make_input_list(input_dir):
    
    input_files = []

    runs = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
    if len(runs) > 1:
        print('Found nested runs')

        for run in runs:

            if run == 'multi':
                continue

            run_dir = os.path.join(input_dir, run+'/AOD/')

            subruns = [d for d in os.listdir(run_dir) if os.path.isdir(os.path.join(run_dir, d))]

            for subrun in subruns:
                file_path = os.path.join(run_dir, subrun, 'AO2D.root')
                if os.path.isfile(file_path):
                    input_files.append(file_path)
                # input_files.append(os.path.join(run_dir, subrun, 'AO2D.root'))
    else:
        print('No nested run structure')
        input_files = [f for f in os.listdir(input_dir) if f.endswith('.root')]
        input_files = [os.path.join(input_dir, f) for f in input_files]

    with open(f'{input_dir}/input_data.txt', 'w') as f:
        for file in input_files:
            f.write(os.path.join(input_dir, file) + '\n')

    print(f'Saved {len(input_files)} files to {input_dir}/input_data.txt')

# run_dir = '/home/kaareendrup/analysis/localMuons_table/input_data/DQ_data'
run_dir = '/home/kaareendrup/analysis/input_data/c3_global'
# run_dir = '/home/kaareendrup/analysis/input_data/c3_standalone'
# run_dir = '/home/kaareendrup/analysis/input_data/f4d_global'
# run_dir = '/home/kaareendrup/analysis/input_data/f4d_standalone'
# run_dir = '/home/kaareendrup/analysis/input_data/DQ_data_global'
# run_dir = '/home/kaareendrup/analysis/input_data/DQ_data_standalone'
make_input_list(run_dir)