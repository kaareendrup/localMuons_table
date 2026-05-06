
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TString.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>

TTree* get_tree(TKey *key, TFile *file, TString treeName) {
    
    // Get the tree stored in the dataframe identified by the key
    const char* dirName = key->GetName();
    TDirectoryFile *dir = (TDirectoryFile*)file->Get(dirName);
    TTree *tree = (TTree*)dir->Get(treeName);
    return tree;
}

void count_input(TString data_name, int n_files, bool MC = false) {

    // Load input data file list
    TString input_file_list = TString::Format("input_data/%s/input_data_%d.txt", data_name.Data(), n_files);

    std::ifstream infile(input_file_list.Data());
    std::string line;

    // Initialize event counting variable
    int event_count = 0;
    int muon_count = 0;
    int gen_muon_count = 0;
    Int_t pdgCode;
    int file_count = 0;

    // Loop over input files
    while (std::getline(infile, line)) {

        std::cout << TString::Format("Processing file %d of %d", file_count, n_files) << std::endl;
        TFile *file = TFile::Open(line.c_str());
        auto *keys = file->GetListOfKeys();
        TTree *event_tree, *muon_tree, *track_tree;

        // Loop over dataframes to count events and muons
        for (int i = 0; i < keys->GetEntries()-1; ++i) {
            event_tree = get_tree((TKey*)keys->At(i), file, "O2reducedevent");
            event_count += event_tree->GetEntries();
            muon_tree = get_tree((TKey*)keys->At(i), file, "O2reducedmuon");
            muon_count += muon_tree->GetEntries();

            delete event_tree;
            delete muon_tree;

            if (!MC) continue;
            track_tree = get_tree((TKey*)keys->At(i), file, "O2reducedmctrack");
            track_tree->SetBranchAddress("fPdgCode", &pdgCode);

            for (Long64_t j = 0; j < track_tree->GetEntries(); ++j) {
                track_tree->GetEntry(j);
                if (abs(pdgCode) == 13) {
                    gen_muon_count++;
                }
            }

            delete track_tree;
        }

        file->Close();
        file_count++;
    }
    std::cout << std::endl;

    // Print the total number of unique events
    std::cout << "Input data stats for: " << data_name << std::endl;
    std::cout << "  Total events: " << event_count << std::endl;
    std::cout << "  Total reconstructed muon tracks (global and standalone): " << muon_count << std::endl;
    if (MC) {
        std::cout << "  Total generated muons: " << gen_muon_count << std::endl << std::endl;
    }
}

void count_analysis(TString data_name, int n_files) {

    // Initialize event counting variable
    int event_count = 0;
    int muon_count = 0;

    for (int f = 0; f < n_files; ++f) {
        std::cout << "Processing file " << f << " of " << n_files << std::endl;
        TString data_file = "results/" + data_name + "/muonAOD" + f + ".root";

        // Load the dataframe keys
        TFile *file = TFile::Open(data_file);
        auto *keys = file->GetListOfKeys();

        // Loop over dataframes
        for (int i = 0; i < keys->GetEntries()-1; ++i) {

            // std::cout << "Processing " << i << "/" << keys->GetEntries()-1 << std::endl;

            TTree *tree = get_tree((TKey*)keys->At(i), file, "O2dqmuontable");

            std::set<ULong64_t> unique_event_indices;
            ULong64_t fEventIdx;
            tree->SetBranchAddress("fEventIdx", &fEventIdx);
            
            Long64_t n = tree->GetEntries();
            for (Long64_t j = 0; j < n; ++j) {
                tree->GetEntry(j);
                unique_event_indices.insert(fEventIdx);
            }

            event_count += unique_event_indices.size();
            muon_count += n;

            delete tree;
        }
    }
    std::cout << std::endl;
        
    // Print the total number of unique events
    std::cout << "Analysis data stats for: " << data_name << std::endl;
    std::cout << "  Total events: " << event_count << std::endl;
    std::cout << "  Total muons: " << muon_count << std::endl << std::endl;
}


void count_events() {

    // TString MC_name = "c3_global";
    // count_input(MC_name, 200, MC = true);
    // TString MC_gen_name = "c3_global/gen";
    // count_analysis(MC_gen_name, 25);
    // TString MC_reco_name = "c3_global/reco";
    // count_analysis(MC_reco_name, 25);
    
    // TString MC_name = "c3_standalone";
    // count_input(MC_name, 200), MC = true);
    // TString MC_gen_name = "c3_standalone/gen";
    // count_analysis(MC_gen_name, 25);
    // TString MC_reco_name = "c3_standalone/reco";
    // count_analysis(MC_reco_name, 25);

    TString data_name = "DQ_data_standalone";
    count_input(data_name, 400);
    TString data_reco_name = "DQ_data_standalone/reco";
    count_analysis(data_reco_name, 50);
}