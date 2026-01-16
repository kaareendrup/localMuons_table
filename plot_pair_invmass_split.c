#include "setALICEStyle.c"

std::cout << std::fixed << std::setprecision(1);
SetALICEStyle();

TTree* get_tree(TKey *key, TFile *file) {

    const char* dirName = key->GetName();
    TDirectoryFile *dir = (TDirectoryFile*)file->Get(dirName);
    TTree *tree = (TTree*)dir->Get("O2dqmuontable");

    return tree;
}

void plot_pair_invmass_split(){

    TString MC_name = "DQ";
    // TString MC_name = "HF";
    // TString MC_name = "genpurp";
    TString data_file = "results/" + MC_name + "/muonAOD.root";
    
    float range_min = 0;
    float range_max = 6;
    int n_bins = 100;

    TString motherLabel = "fMother";
    // TString motherLabel = "fGrandmother";

    // Load the dataframe keys
    TFile *file = TFile::Open(data_file);
    auto *keys = file->GetListOfKeys();

    // Initialize inv mass variables
    std::vector<Double_t> all_inv_masses, JPsi_inv_masses, Psi2S_inv_masses, other_inv_masses;

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

        // Prepare to read MC truth info
        Long64_t fMotherID, fMotherPDG;
        tree->SetBranchAddress(motherLabel + "ID", &fMotherID);
        tree->SetBranchAddress(motherLabel + "PDG", &fMotherPDG);

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

            std::vector<Long64_t> muon_motherIDs, muon_motherPDGs;
            std::vector<ROOT::Math::PtEtaPhiMVector> muon_vectors;

            // Read muon kinematics and build 4-vectors
            for (auto entry : muon_entries) {
                tree->GetEntry(entry);
                ROOT::Math::PtEtaPhiMVector muon_vec(fPt, fEta, fPhi, 0.105658); // Muon mass ~105.658 MeV/c^2
                muon_vectors.push_back(muon_vec);

                muon_motherIDs.push_back(fMotherID);
                muon_motherPDGs.push_back(fMotherPDG);
            }

            // Pair all muons and store invariant masses
            for (size_t j = 0; j < muon_vectors.size(); ++j) {
                for (size_t k = j + 1; k < muon_vectors.size(); ++k) {

                    auto inv_mass = (muon_vectors[j] + muon_vectors[k]).M();
                    all_inv_masses.push_back(inv_mass);

                    // Categorize based on common mother ID and PDG
                    if (muon_motherIDs[j] == muon_motherIDs[k]) {
                        if (muon_motherPDGs[j] == 443 && muon_motherPDGs[k] == 443) {
                            JPsi_inv_masses.push_back(inv_mass);
                        }
                        else if (muon_motherPDGs[j] == 100443 && muon_motherPDGs[k] == 100443) {
                            Psi2S_inv_masses.push_back(inv_mass);
                        }
                        else {
                            other_inv_masses.push_back(inv_mass);
                        }
                    } else {
                        other_inv_masses.push_back(inv_mass);
                    }
                }
            }
        }
    }

    // Plot histogram of candidate invariant masses
    TCanvas *c1 = new TCanvas("c1", "Invariant Mass of Muon Pairs", 800, 600);
    TH1F *invMassHist = new TH1F("h1","Invariant Mass of Muon Pairs;Invariant Mass (GeV/c^{2});Counts",n_bins,range_min,range_max);
    invMassHist->FillN(all_inv_masses.size(), all_inv_masses.data(), nullptr);
    invMassHist->SetLineColor(kBlack); 
    invMassHist->SetLineWidth(3);
    invMassHist->Sumw2();
    invMassHist->Draw("HIST");
    
    // Add secondary histograms
    TH1F *jpsiHist = new TH1F("h2","",n_bins,range_min,range_max);
    jpsiHist->FillN(JPsi_inv_masses.size(), JPsi_inv_masses.data(), nullptr);
    jpsiHist->SetLineColor(kGreen+1);
    jpsiHist->SetLineWidth(2);
    jpsiHist->Sumw2();
    jpsiHist->Draw("HIST SAME");
    jpsiHist->Draw("E SAME");
    
    TH1F *psi2sHist = new TH1F("h3","",n_bins,range_min,range_max);
    psi2sHist->FillN(Psi2S_inv_masses.size(), Psi2S_inv_masses.data(), nullptr);
    psi2sHist->SetLineColor(kGreen+2);
    psi2sHist->SetLineWidth(2);
    psi2sHist->Sumw2();
    psi2sHist->Draw("HIST SAME");
    psi2sHist->Draw("E SAME");
    
    TH1F *otherHist = new TH1F("h4","",n_bins,range_min,range_max);
    otherHist->FillN(other_inv_masses.size(), other_inv_masses.data(), nullptr);
    otherHist->SetLineColor(kRed);
    otherHist->SetLineWidth(2);
    otherHist->Sumw2();
    otherHist->Draw("HIST SAME");
    otherHist->Draw("E SAME");
    
    // Style
    invMassHist->SetLineWidth(3); 
    increaseMargins(c1);

    // Legend
    TLegend *legend = new TLegend(0.7,0.7,0.88,0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(invMassHist, "All", "l");
    legend->AddEntry(jpsiHist, "J/#psi", "l");
    legend->AddEntry(psi2sHist, "#psi(2S)", "l");
    legend->AddEntry(otherHist, "Other", "l");
    legend->Draw();

    drawLabel(MC_name);
    TString out_name = TString::Format("results/%s/muon_pair_invariant_mass_%.1f_%.1f_%s", MC_name.Data(), range_min, range_max, motherLabel.Data());
    out_name.ReplaceAll(".", "_");
    c1->SaveAs(out_name + ".png");
}