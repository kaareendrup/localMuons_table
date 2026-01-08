#include "setALICEStyle.c"

std::cout << std::fixed << std::setprecision(1);
SetALICEStyle();

TTree* get_tree(TKey *key, TFile *file) {

    // Get the tree stored in the dataframe identified by the key
    const char* dirName = key->GetName();
    TDirectoryFile *dir = (TDirectoryFile*)file->Get(dirName);
    TTree *tree = (TTree*)dir->Get("O2dqmuontable");

    return tree;
}

struct RangeMap {
    int min;
    int max;
    TString name;
};

static const std::vector<RangeMap> particleRanges = {
    {13, 13, "mu-"},
    {211, 211, "pi+"},
    {321, 321, "K+"},
    {411, 435, "c mesons"},
    {441, 441, "Eta_c(1S)"},
    {443, 445, "J/Psi"},
    {511, 557, "b mesons"},
    {2212, 2212, "p"},
    {3112, 3334, "s baryons"},
    {4101, 4103, "c diquarks"},
    {4112, 4444, "c baryons"},
    {5101, 5103, "b diquarks"},
    {5112, 5554, "b baryons"},
    {100443, 100443, "Psi(2S)"},
    {9999, 9999, "undefined"},
};

static const std::vector<RangeMap> particleRangesSimple = {
    {13, 13, "mu-"},
    {211, 211, "pi+"},
    {321, 321, "strange"},
    {411, 435, "charm"},
    {441, 445, "J/Psi"},
    {511, 557, "beauty"},
    {2212, 2212, "p"},
    {3112, 3334, "strange"},
    {4112, 4444, "charm"},
    {5112, 5554, "beauty"},
    {100443, 100443, "Psi(2S)"},
    {9999, 9999, "undefined"},
};

TString particleName(int code, std::vector<RangeMap> Ranges) {
    for (const auto& r : Ranges) {
        if (code >= r.min && code <= r.max)
            return r.name;
    }
    return TString::Format("PDG %d", code);
}

void insert_unique(std::set<TString> *uniqueValues, TKey *key, TFile *file, const char* branchName, const bool simple) {

    // Get dataframe tree and set mother/grandmother PDG branch address
    Long64_t PDG;
    TTree *tree = get_tree(key, file);
    tree->SetBranchAddress(branchName, &PDG);

    // Loop over entries and insert unique values converted to categories
    for (Long64_t j = 0; j < tree->GetEntries(); ++j) {
        tree->GetEntry(j);
        if (simple){
            uniqueValues->insert(particleName(std::abs(PDG), particleRangesSimple));
        }
        else{
            uniqueValues->insert(particleName(std::abs(PDG), particleRanges));
        }
    }
}

void plot_PID_hist_collapsed(){

    // TString MC_name = "DQ";
    TString MC_name = "HF";
    // TString MC_name = "genpurp";
    TString data_file = "results/" + MC_name + "/muonAOD.root";
    bool simple = false;
    // bool simple = true;
    
    // Load dataframe keys 
    TFile *file = TFile::Open(data_file);
    auto *keys = file->GetListOfKeys();

    // Initialize counting variables
    std::set<TString> uniqueValues_motherPDG, uniqueValues_grandmotherPDG;
    std::map<TString, int> motherPDG_counts, grandmotherPDG_counts;
    
    // First loop over dataframes to find unique values
    for (int i = 0; i < keys->GetEntries()-1; ++i) {
        insert_unique(&uniqueValues_motherPDG, (TKey*)keys->At(i), file, "fMotherPDG", simple);
        insert_unique(&uniqueValues_grandmotherPDG, (TKey*)keys->At(i), file, "fGrandmotherPDG", simple);
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

        // Loop over entries and increment counts for the corresponding PDG categories
        for (Long64_t j = 0; j < tree->GetEntries(); ++j) {
            tree->GetEntry(j);
            if (simple){
                motherPDG_counts[particleName(std::abs(motherPDG), particleRangesSimple)]++;
                grandmotherPDG_counts[particleName(std::abs(grandmotherPDG), particleRangesSimple)]++;
                
            }
            else{
                motherPDG_counts[particleName(std::abs(motherPDG), particleRanges)]++;
                grandmotherPDG_counts[particleName(std::abs(grandmotherPDG), particleRanges)]++;
            }
        }
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
            TString label = pair.first;
            
            motherPID->SetBinContent(i+1, pair.second);
            motherPID->GetXaxis()->SetBinLabel(i+1, label);
        }

        for (int i = 0; i < uniqueValues_grandmotherPDG.size(); ++i) {
            auto& pair = *std::next(grandmotherPDG_counts.begin(), i);
            TString label = pair.first;

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
        motherPID->GetXaxis()->SetLabelSize(.07);
        grandmotherPID->GetYaxis()->SetTitleOffset(.6);
        grandmotherPID->GetXaxis()->SetLabelSize(.07);

        if (y_log) {
            c1->SaveAs("results/" + MC_name + "/PID_histograms_collapsed_log.png");
        } else {
            c1->SaveAs("results/" + MC_name + "/PID_histograms_collapsed.png");
        }

        delete c1;
        delete motherPID;
        delete grandmotherPID;
    }
}