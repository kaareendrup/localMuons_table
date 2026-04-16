import json

def import_bins(config_file):
    with open(config_file) as f:
        config = json.load(f)

    with open(f'results/{config["hepdata_name"]}.json') as f:
        hepdata = json.load(f)

    bins = []
    for point in hepdata["values"]:
        pt_low  = float(point["x"][0]["low"])
        pt_high = float(point["x"][0]["high"])

        if pt_low >= config["cuts_JPsi"]["pT_JPsi_min"] and pt_low <= config["cuts_JPsi"]["pT_JPsi_max"]:
            bins.append(pt_low)

    # Export to config json
    config["hists"]["pT_bins"] = [bin for bin in bins]
    with open(config_file, 'w') as f:
        json.dump(config, f, indent=4)

    print(f"Imported {len(bins)} bins from HEPData and updated {config_file}.")

config_file = 'localMuons_table/config/config_analysis.json'
import_bins(config_file)