def make_command(n_files, data_name, tracktype, recotype, MC, PbPb, config_dir):

    datatype = "PbPb" if PbPb else ""

    if MC:
        base_command = f"o2-analysis-dq-efficiency-with-assoc -b --configuration json://{config_dir}/configuration_dqEfficiency_withAssoc_{tracktype}_{recotype}.json"
    elif PbPb:
        base_command = f"o2-analysis-dq-table-reader -b --configuration json://{config_dir}/configuration_dqTableReader_{tracktype}{datatype}.json | o2-analysis-dq-model-converter-event-extended -b --configuration json://{config_dir}/configuration_dqTableReader_{tracktype}{datatype}.json"
    else:
        base_command = f"o2-analysis-dq-table-reader-with-assoc -b --configuration json://{config_dir}/configuration_dqTableReader_withAssoc_{tracktype}{datatype}.json | o2-analysis-dq-model-converter-event-extended -b --configuration json://{config_dir}/configuration_dqTableReader_withAssoc_{tracktype}{datatype}.json"

    if PbPb:
        input_file = f"/home/kaareendrup/analysis/input_data/{data_name}/input_data_{n_files}.txt"
    else:
        input_file = f"/home/kaareendrup/analysis/input_data/{data_name}_{tracktype}/input_data_{n_files}.txt"

    command = f"{base_command} --aod-file @{input_file} --aod-writer-json {config_dir}/outputDirector_{recotype}.json"

    print("To run, run:")
    print(command)

config_dir = "/home/kaareendrup/analysis/localMuons_table/config"

n_files = 10

# data_name = "c3"
# data_name = "f4d"
# data_name = "DQ_data"
data_name = "DQ_PbPb"

# tracktype = "global"
tracktype = "standalone"

recotype = "reco"
# recotype = "gen"

MC = data_name in ["c3", "f4d"]
PbPb = "PbPb" in data_name

make_command(n_files, data_name, tracktype, recotype, MC, PbPb, config_dir)