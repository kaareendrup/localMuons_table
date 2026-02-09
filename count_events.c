
TTree* get_tree(TKey *key, TFile *file, TString treeName) {
    
    // Get the tree stored in the dataframe identified by the key
    const char* dirName = key->GetName();
    TDirectoryFile *dir = (TDirectoryFile*)file->Get(dirName);
    TTree *tree = (TTree*)dir->Get(treeName);
    return tree;
}

void count_input(TString MC_name) {

    // Load input data file list
    TString input_file_list = "input_data/" + MC_name + "/input_data_local.txt";
    std::ifstream infile(input_file_list.Data());
    std::string line;

    // Initialize event counting variable
    int event_count = 0;
    int muon_count = 0;

    // Loop over input files
    while (std::getline(infile, line)) {

        TFile *file = TFile::Open(line.c_str());
        auto *keys = file->GetListOfKeys();
        TTree *event_tree, *muon_tree;

        // Loop over dataframes to count events and muons
        for (int i = 0; i < keys->GetEntries()-1; ++i) {
            event_tree = get_tree((TKey*)keys->At(i), file, "O2reducedevent");
            event_count += event_tree->GetEntries();
            muon_tree = get_tree((TKey*)keys->At(i), file, "O2reducedmuon");
            muon_count += muon_tree->GetEntries();

            delete event_tree;
            delete muon_tree;
        }
    }

    // Print the total number of unique events
    std::cout << "Input data stats for: " << MC_name << std::endl;
    std::cout << "  Total events: " << event_count << std::endl;
    std::cout << "  Total muons: " << muon_count << std::endl << std::endl;
}

void count_analysis(TString MC_name) {
    TString data_file = "results/" + MC_name + "/muonAOD.root";
    
    // Load the dataframe keys
    TFile *file = TFile::Open(data_file);
    auto *keys = file->GetListOfKeys();
    
    // Initialize event counting variable
    int event_count = 0;
    int muon_count = 0;
    
    // Loop over dataframes
    for (int i = 0; i < keys->GetEntries()-1; ++i) {

        std::cout << "Processing " << i << "/" << keys->GetEntries()-1 << std::endl;

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
    
    // Print the total number of unique events
    std::cout << "Analysis data stats for: " << MC_name << std::endl;
    std::cout << "  Total events: " << event_count << std::endl;
    std::cout << "  Total muons: " << muon_count << std::endl << std::endl;
}

    
void count_events() {
    // TString MC_name = "DQ";
    // TString MC_name = "HF";
    // TString MC_name = "genpurp";
    TString MC_name = "k4h_baseline/reco";
    // TString MC_name = "k4h_standalone/reco";
    // TString MC_name = "k4h_standalone/gen";

    // count_input(MC_name);
    count_analysis(MC_name);
}