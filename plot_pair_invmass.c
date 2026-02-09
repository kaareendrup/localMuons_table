#include "setALICEStyle.c"

std::cout << std::fixed << std::setprecision(1);
SetALICEStyle();

TTree* get_tree(TKey *key, TFile *file) {

    const char* dirName = key->GetName();
    TDirectoryFile *dir = (TDirectoryFile*)file->Get(dirName);
    TTree *tree = (TTree*)dir->Get("O2dqmuontable");

    return tree;
}

void plot_pair_invmass(){

    TString data_name = "DQ_data";
    // TString data_name = "DQ";
    // TString data_name = "HF";
    // TString data_name = "genpurp";
    // TString data_name = "k4h_standalone";
    
    TString type = "reco";
    // TString type = "gen";

    TString data_file = "results/" + data_name + "/" + type + "/muonAOD.root";
    
    // float range_min = 2.5;
    // float range_max = 3.5;
    float range_min = 0;
    float range_max = 5;

    // Load the dataframe keys
    TFile *file = TFile::Open(data_file);
    auto *keys = file->GetListOfKeys();

    // Initialize inv mass variables
    std::vector<Double_t> candidate_inv_masses;

    // Loop over dataframes
    for (int i = 0; i < keys->GetEntries()-1; ++i) {
        TTree *tree = get_tree((TKey*)keys->At(i), file);

        // Group muons by event index
        std::map<ULong64_t, std::vector<Long64_t>> muon_groups;
        ULong64_t fEventIdx;
        tree->SetBranchAddress("fEventIdx", &fEventIdx);
        
        // First pass: build groups of muons from the same event
        Long64_t n = tree->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            tree->GetEntry(i);
            muon_groups[fEventIdx].push_back(i);
        }

        // Prepare to read muon kinematics
        float fPt, fPhi, fEta;
        tree->SetBranchAddress("fPtassoc", &fPt);
        tree->SetBranchAddress("fPhiassoc", &fPhi);
        tree->SetBranchAddress("fEtaassoc", &fEta);

        // Second pass: process each group
        for (auto& event : muon_groups) {
            ULong64_t eventID = event.first;
            auto& muon_entries = event.second;

            if (muon_entries.size() < 2) continue; // Needs at least 2 muons to form a pair
            // if (muon_entries.size() < 3) continue; // Needs at least 3 muons to form a pair + associate

            std::vector<ROOT::Math::PtEtaPhiMVector> muon_vectors;

            // Read muon kinematics and build 4-vectors
            for (auto entry : muon_entries) {
                tree->GetEntry(entry);
                ROOT::Math::PtEtaPhiMVector muon_vec(fPt, fEta, fPhi, 0.105658); // Muon mass ~105.658 MeV/c^2
                muon_vectors.push_back(muon_vec);
            }

            // Store the indexes of the best candidate muon pair for potential use later
            float inv_mass_cand = 1e9;
            int idx_cand_1, idx_cand_2;

            // Find the muon pair with invariant mass closest to J/Psi mass
            for (size_t j = 0; j < muon_vectors.size(); ++j) {
                for (size_t k = j + 1; k < muon_vectors.size(); ++k) {
                    auto inv_mass = (muon_vectors[j] + muon_vectors[k]).M();

                    // Replace candidate if closer to J/Psi mass
                    if (std::abs(inv_mass - 3.0969) < std::abs(inv_mass_cand - 3.0969)) { // J/Psi mass ~3.0969 GeV/c^2
                        inv_mass_cand = inv_mass;
                        idx_cand_1 = j;
                        idx_cand_2 = k;
                    }
                }
            }

            // Store the best candidate invariant mass
            candidate_inv_masses.push_back(inv_mass_cand);
        }

        tree->ResetBranchAddresses();
        delete tree;
    }

    // Plot histogram of candidate invariant masses
    TCanvas *c1 = new TCanvas("c1", "Invariant Mass of Muon Pairs", 800, 600);
    TH1F *invMassHist = new TH1F("h1","Invariant Mass of Muon Pairs;Invariant Mass (GeV/c^{2});Counts",100,range_min,range_max);
    invMassHist->FillN(candidate_inv_masses.size(), candidate_inv_masses.data(), nullptr);
    invMassHist->Draw();
    invMassHist->Sumw2();

    // Style
    invMassHist->SetLineWidth(3); 
    invMassHist->SetLineColor(kBlack); 
    invMassHist->SetMarkerStyle(kFullCircle); 
    invMassHist->SetMarkerSize(.5); 
    invMassHist->SetMarkerColor(kRed);
    increaseMargins(c1);

    drawLabel(data_name, type);
    TString out_name = TString::Format("results/%s/%s/muon_pair_invariant_mass_%.1f_%.1f", data_name.Data(), type.Data(), range_min, range_max);
    out_name.ReplaceAll(".", "_");
    c1->SaveAs(out_name + ".png");
}