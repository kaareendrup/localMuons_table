#include "setALICEStyle.c"

// Set style
std::cout << std::fixed << std::setprecision(1);
SetALICEStyle();

TTree* get_tree(TKey *key, TFile *file) {

    // Get the tree stored in the dataframe identified by the key
    const char* dirName = key->GetName();
    TDirectoryFile *dir = (TDirectoryFile*)file->Get(dirName);
    TTree *tree = (TTree*)dir->Get("O2dqmuontable");

    return tree;
}

void insert_unique(std::set<Long64_t> *uniqueValues, TKey *key, TFile *file, const char* branchName) {

    // Get dataframe tree and set mother/grandmother PDG branch address
    Long64_t PDG;
    TTree *tree = get_tree(key, file);
    tree->SetBranchAddress(branchName, &PDG);

    // Loop over entries and insert unique values (absolute value for non-null PDG)
    for (Long64_t j = 0; j < tree->GetEntries(); ++j) {
        tree->GetEntry(j);
        if (PDG != -9999){
            uniqueValues->insert(std::abs(PDG));
        } else {
            uniqueValues->insert(PDG);
        }
    }

    delete tree;
}

void plot_PID_hist(){

    // Load the tree
    // TString MC_name = "DQ";
    TString MC_name = "DQ_gen";
    // TString MC_name = "HF";
    // TString MC_name = "genpurp";
    TString data_file = "results/" + MC_name + "/muonAOD.root";
    
    // Load dataframe keys 
    TFile *file = TFile::Open(data_file);
    auto *keys = file->GetListOfKeys();
    
    // Initialize counting variables
    std::set<Long64_t> uniqueValues_motherPDG, uniqueValues_grandmotherPDG;
    std::map<Long64_t, int> motherPDG_counts, grandmotherPDG_counts;
    
    // First loop over dataframes to find unique values
    for (int i = 0; i < keys->GetEntries()-1; ++i) {
        insert_unique(&uniqueValues_motherPDG, (TKey*)keys->At(i), file, "fMotherPDG");
        insert_unique(&uniqueValues_grandmotherPDG, (TKey*)keys->At(i), file, "fGrandmotherPDG");
    }
    
    // Initialize counts to zero
    for (const auto& val : uniqueValues_motherPDG) {
        motherPDG_counts[val] = 0;
    }
    for (const auto& val : uniqueValues_grandmotherPDG) {
        grandmotherPDG_counts[val] = 0;
    }
    
    // Intialize varliables for PDG and dataframe tree
    Long64_t motherPDG, grandmotherPDG;
    TTree *tree;
    
    // Second pass to count occurrences
    for (int i = 0; i < keys->GetEntries()-1; ++i) {

        // Get the tree and set branch addresses
        tree = get_tree((TKey*)keys->At(i), file);
        tree->SetBranchAddress("fMotherPDG", &motherPDG);
        tree->SetBranchAddress("fGrandmotherPDG", &grandmotherPDG);

        // Loop over entries and increment counts for the corresponding PDG values
        for (Long64_t j = 0; j < tree->GetEntries(); ++j) {
            tree->GetEntry(j);
            if (motherPDG != -9999){
                motherPDG_counts[std::abs(motherPDG)]++;
            } else {
                motherPDG_counts[motherPDG]++;
            }
            if (grandmotherPDG != -9999){
                grandmotherPDG_counts[std::abs(grandmotherPDG)]++;
            } else {
                grandmotherPDG_counts[grandmotherPDG]++;
            }
        }

        delete tree;
    }

    // Plot with both log and linear y-axis
    for (bool y_log : {false, true}) {

        // Draw histograms
        TCanvas *c1 = new TCanvas("c1", "PID Histograms", 1200, 600);
        c1->Divide(1,2);
        TH1F *motherPID = new TH1F("h1","Mother PDG Values",uniqueValues_motherPDG.size(),0,uniqueValues_motherPDG.size());
        TH1F *grandmotherPID = new TH1F("h2","Grandmother PDG Values",uniqueValues_grandmotherPDG.size(),0,uniqueValues_grandmotherPDG.size());

        for (int i = 0; i < uniqueValues_motherPDG.size(); ++i) {
            auto& pair = *std::next(motherPDG_counts.begin(), i);
            TString label = TString::Format("%.0lld", pair.first);

            motherPID->SetBinContent(i+1, pair.second);
            motherPID->GetXaxis()->SetBinLabel(i+1, label);
        }
        for (int i = 0; i < uniqueValues_grandmotherPDG.size(); ++i) {
            auto& pair = *std::next(grandmotherPDG_counts.begin(), i);
            TString label = TString::Format("%.0lld", pair.first);

            grandmotherPID->SetBinContent(i+1, pair.second);
            grandmotherPID->GetXaxis()->SetBinLabel(i+1, label);
        }

        // Plot mother PID
        c1->cd(1);
        motherPID->GetXaxis()->SetTitle("Mother PDG");
        motherPID->GetYaxis()->SetTitle("Counts");
        motherPID->SetFillColor(kBlue-10);
        motherPID->Draw();
        drawLabel(MC_name, 0.82);

        if (y_log) {
            gPad->SetLogy();
        }

        // Plot grandmother PID
        c1->cd(2);
        grandmotherPID->GetXaxis()->SetTitle("Grandmother PDG");
        grandmotherPID->GetYaxis()->SetTitle("Counts");
        grandmotherPID->SetFillColor(kRed-10);
        grandmotherPID->Draw();
        drawLabel(MC_name, 0.82);

        if (y_log) {
            gPad->SetLogy();
        }

        // Style
        increasePadMargins(c1, 2);
        motherPID->GetYaxis()->SetTitleOffset(.6);
        grandmotherPID->GetYaxis()->SetTitleOffset(.6);
        grandmotherPID->GetXaxis()->SetTitleOffset(1.7);
        gPad->SetBottomMargin(0.2);

        if (y_log) {
            c1->SaveAs("results/" + MC_name + "/PID_histograms_log.png");
        } else {
            c1->SaveAs("results/" + MC_name + "/PID_histograms.png");
        }

        delete c1;
        delete motherPID;
        delete grandmotherPID;
    }
}