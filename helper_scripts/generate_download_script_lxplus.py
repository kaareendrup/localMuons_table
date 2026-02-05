
def generate_download_script(run_dir, out_dir, begin, end):

    with open('file_list.txt', 'w') as f:
        for i in range(begin, end + 1):
            file_path = f"{i}/reducedAod{i}.root"
            f.write(file_path + '\n')

    print(f"File list generated from {begin} to {end} in 'file_list.txt'.")
    print("To download the files, use the following command:")
    print(f"rsync -av --progress --no-relative --files-from=file_list.txt -e ssh kiversen@lxplus.cern.ch:{run_dir} {out_dir}")

run_dir = '/afs/cern.ch/user/k/kiversen/eos/O2DQ_withTableMaker_MC_genpurp/rundirs/k4h_global'
out_dir = '/media/kaareendrup/ec65dbb9-11dd-4abf-ae3e-f029466fe958/analysis/input_data/k4h'
begin = 20
end = 24

generate_download_script(run_dir, out_dir, begin, end)