#include "setALICEStyle.c"

std::cout << std::fixed << std::setprecision(1);
SetALICEStyle();

void inspect_mothers(TString MC_name, TString type) {

    TString data_file = "results/" + MC_name + "/" + type + "/muonAOD.root";
    TFile *file = TFile::Open(data_file, "READ");

    // Loop over all keys in the file (directories etc.)
    TIter nextKey(file->GetListOfKeys());
    TKey* key;
    int dirCount = 0;

    std::vector<Long64_t> motherPDGs;
    std::vector<Long64_t> grandmotherPDGs;
    std::vector<Long64_t> nMothers;

    Long64_t totalMuons = 0;

    while ((key = (TKey*) nextKey())) {
        TObject* obj = key->ReadObj();

        // --- Handle only TDirectoryFile objects ---
        if (obj->InheritsFrom("TDirectory")) {
            TDirectory* dir = (TDirectory*) obj;
            printf("Reading tracks from dir %d of %d: %s\n", dirCount, file->GetListOfKeys()->GetEntries(), dir->GetName());

            TTree* MCMuons = (TTree*)dir->Get("O2dqmuontable");

            Long64_t fMotherPDG, fGrandmotherPDG, fNMothers;
            MCMuons->SetBranchAddress("fMotherPDG", &fMotherPDG);
            MCMuons->SetBranchAddress("fGrandmotherPDG", &fGrandmotherPDG);
            MCMuons->SetBranchAddress("fNMothers", &fNMothers);

            for (Long64_t i = 0; i < MCMuons->GetEntries(); ++i) {
                MCMuons->GetEntry(i);
                totalMuons++;

                if ((std::abs(fMotherPDG) == 443) || (std::abs(fMotherPDG) == 100443) ||
                    (std::abs(fMotherPDG) >= 411 && std::abs(fMotherPDG) <= 445) || 
                    (std::abs(fMotherPDG) >= 4101 && std::abs(fMotherPDG) <= 4444)) {
                    
                    if ((std::abs(fGrandmotherPDG) >= 511 && std::abs(fGrandmotherPDG) <= 557) || 
                        (std::abs(fGrandmotherPDG) >= 5101 && std::abs(fGrandmotherPDG) <= 5554)) {
                        motherPDGs.push_back(std::abs(fMotherPDG));
                        grandmotherPDGs.push_back(std::abs(fGrandmotherPDG));
                        nMothers.push_back(fNMothers);
                    }
                }
            }

            MCMuons->ResetBranchAddresses();
            delete MCMuons;
        }

        dirCount++;
    }

    file->Close();

    // Print summary
    std::cout << "Found " << motherPDGs.size() << " out of " << totalMuons << " muons from charm hadrons with beauty grandmothers." << std::endl;
    std::cout << "Unique mother PDGs:" << std::endl;
    std::set<Long64_t> uniqueMothers(motherPDGs.begin(), motherPDGs.end());
    for (const auto& pdg : uniqueMothers) {
        std::cout << "  PDG: " << pdg << " Count: " 
                  << std::count(motherPDGs.begin(), motherPDGs.end(), pdg) << std::endl;
    }
    std::cout << "Unique grandmother PDGs:" << std::endl;
    std::set<Long64_t> uniqueGrandmothers(grandmotherPDGs.begin(), grandmotherPDGs.end());
    for (const auto& pdg : uniqueGrandmothers) {
        std::cout << "  PDG: " << pdg << " Count: " 
                  << std::count(grandmotherPDGs.begin(), grandmotherPDGs.end(), pdg) << std::endl;
    }
    std::cout << "Distribution of number of mothers:" << std::endl;
    std::set<Long64_t> uniqueNMothers(nMothers.begin(), nMothers.end());
    for (const auto& n : uniqueNMothers) {
        std::cout << "  nMothers: " << n << " Count: " 
                  << std::count(nMothers.begin(), nMothers.end(), n) << std::endl;
    }

    // Get counts of unique combinations of mother and grandmother PDGs
    std::map<std::pair<Long64_t, Long64_t>, int> combinationCounts;
    for (size_t i = 0; i < motherPDGs.size(); ++i) {
        combinationCounts[{motherPDGs[i], grandmotherPDGs[i]}]++;
    }

    // Plot 2D histogram of mother vs grandmother PDG
    TCanvas *c1 = new TCanvas("c1", "Mother vs Grandmother PDG", 800, 600);
    TH2F *hist = new TH2F("h1", "Mother vs Grandmother PDG;Mother PDG;Grandmother PDG", 
        uniqueMothers.size(),0,uniqueMothers.size(), 
        uniqueGrandmothers.size(),0,uniqueGrandmothers.size());

    // Fill histogram with counts
    for (const auto& pair : combinationCounts) {
        Long64_t motherPDG = pair.first.first;
        Long64_t grandmotherPDG = pair.first.second;
        int count = pair.second;
        int xBin = std::distance(uniqueMothers.begin(), uniqueMothers.find(motherPDG)) + 1;
        int yBin = std::distance(uniqueGrandmothers.begin(), uniqueGrandmothers.find(grandmotherPDG)) + 1;
        hist->SetBinContent(xBin, yBin, count);
    }

    // Set bin labels
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        hist->GetXaxis()->SetBinLabel(i, std::to_string(*std::next(uniqueMothers.begin(), i-1)).c_str());
    }
    for (int i = 1; i <= hist->GetNbinsY(); ++i) {
        hist->GetYaxis()->SetBinLabel(i, std::to_string(*std::next(uniqueGrandmothers.begin(), i-1)).c_str());
    }

    hist->Draw("COLZ");
    increaseMargins(c1);
    c1->SetRightMargin(0.15);

    c1->SaveAs(("results/" + MC_name + "/" + type + "/mother_vs_grandmother_PDG.png").Data());
}

void inspect_muon_mothers() {

    TString MC_name = "DQ";
    // TString MC_name = "DQ_gen";
    // TString MC_name = "HF";
    // TString MC_name = "genpurp";

    TString type = "reco";
    // TString type = "gen";

    inspect_mothers(MC_name, type);
}