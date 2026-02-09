#include "setALICEStyle.c"

std::cout << std::fixed << std::setprecision(1);
SetALICEStyle();

void analysis_JPsi_MC() {

    TString MC_name = "DQ";
    // TString MC_name = "DQ_gen";
    // TString MC_name = "HF";
    // TString MC_name = "genpurp";

    TString type = "reco";
    // TString type = "gen";

    TString data_file = "results/" + MC_name + "/" + type + "/muonAOD.root";
    
    TString motherLabel = "fMother";
    // TString motherLabel = "fGrandmother";

    int n_bins = 20;
    
    // Load the dataframe keys
    TFile *file = TFile::Open(data_file);
    TIter nextKey(file->GetListOfKeys());
    TKey* key;
    
    int triggers_JPsi = 0;
    int triggers_Psi2S = 0;
    std::vector<double> deltaEta_JPsi, deltaPhi_JPsi;
    std::vector<double> deltaEta_Psi2S, deltaPhi_Psi2S;

    // Loop over dataframes
    while ((key = (TKey*) nextKey())) {

        // Load directory and tree
        TObject* obj = key->ReadObj();
        if (!(obj->InheritsFrom("TDirectory"))) continue;

        TDirectory* dir = (TDirectory*) obj;
        TTree *tree = (TTree*)dir->Get("O2dqmuontable");

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
        Long64_t fGlobalIndexMCtrack, fMotherID, fMotherPDG;
        ULong64_t fGlobalIndexAssoc;
        tree->SetBranchAddress("fGlobalIndexassoc", &fGlobalIndexAssoc);
        tree->SetBranchAddress("fGlobalIndexMCtrack", &fGlobalIndexMCtrack);
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

            std::vector<std::pair<int, int>> pair_idxs;
            std::set<int> pair_components;

            // Search for pairs with same mother 
            for (int j = 0; j < muon_vectors.size(); ++j) {
                for (int k = j + 1; k < muon_vectors.size(); ++k) {
                    if (muon_motherIDs[j] == muon_motherIDs[k] && muon_motherIDs[j] != -9999) {

                        // Store pair indices and components
                        pair_idxs.push_back(std::make_pair(j, k));
                        pair_components.insert(j);
                        pair_components.insert(k);

                        // Count triggers based on PDG
                        if (muon_motherPDGs[j] == 443 && muon_motherPDGs[k] == 443) {
                            triggers_JPsi++;
                        } else if (muon_motherPDGs[j] == 100443 && muon_motherPDGs[k] == 100443) {
                            triggers_Psi2S++;
                        }
                    }
                }
            }

            // Match pairs to associates
            for (std::pair<int, int> pair : pair_idxs) {
                for (int a = 0; a < muon_vectors.size(); ++a) {

                    // Ensure the associate is not part of a pair
                    if (!(pair_components.count(a))) {
                        double deltaEta = (muon_vectors[pair.first] + muon_vectors[pair.second]).Eta() - muon_vectors[a].Eta();
                        double deltaPhi = (muon_vectors[pair.first] + muon_vectors[pair.second]).Phi() - muon_vectors[a].Phi();
                        if (deltaPhi > 3.0 / 2.0 * M_PI) {
                            deltaPhi -= 2.0 * M_PI;
                        }
                        if (deltaPhi < -0.5 * M_PI) {
                            deltaPhi += 2.0 * M_PI;
                        }

                        // Store deltaEta and deltaPhi values based on PID
                        if (muon_motherPDGs[pair.first] == 443 && muon_motherPDGs[pair.second] == 443) {
                            deltaEta_JPsi.push_back(deltaEta);
                            deltaPhi_JPsi.push_back(deltaPhi);
                        } else if (muon_motherPDGs[pair.first] == 100443 && muon_motherPDGs[pair.second] == 100443) {
                            deltaEta_Psi2S.push_back(deltaEta);
                            deltaPhi_Psi2S.push_back(deltaPhi);
                        }
                    }
                }
            }
        }

        tree->ResetBranchAddresses();
        delete tree;
    }

    // Plot histograms of deltaEta and deltaPhi for JPsi and Psi2S
    TCanvas *c1 = new TCanvas("c1", "Delta Eta and Delta Phi", 1300, 600);
    c1->Divide(2,1);
    c1->cd(1);
    // TH1F *deltaEtaHist_JPsi = new TH1F("h1","Delta Eta JPsi;#Delta#eta;Counts",n_bins,-5,3);
    TH1F *deltaEtaHist_JPsi = new TH1F("h1","Delta Eta JPsi;#Delta#eta;#frac{1}{N_{trig}} dN/d#Delta#eta",n_bins,-5,3);
    deltaEtaHist_JPsi->FillN(deltaEta_JPsi.size(), deltaEta_JPsi.data(), nullptr);
    deltaEtaHist_JPsi->SetLineColor(kGreen+1);
    deltaEtaHist_JPsi->SetLineWidth(2);
    deltaEtaHist_JPsi->Sumw2();
    deltaEtaHist_JPsi->Scale(1.0 / (triggers_JPsi*deltaEtaHist_JPsi->GetBinWidth(1)));
    deltaEtaHist_JPsi->Draw();
    

    TH1F *deltaEtaHist_Psi2S = new TH1F("h2","Delta Eta Psi2S;#Delta#eta;#frac{1}{N_{trig}} dN/d#Delta#eta",n_bins,-5,3);
    deltaEtaHist_Psi2S->FillN(deltaEta_Psi2S.size(), deltaEta_Psi2S.data(), nullptr);
    deltaEtaHist_Psi2S->SetLineColor(kGreen+2);
    deltaEtaHist_Psi2S->SetLineWidth(2);
    deltaEtaHist_Psi2S->Sumw2();
    deltaEtaHist_Psi2S->Scale(1.0 / (triggers_Psi2S*deltaEtaHist_Psi2S->GetBinWidth(1)));
    deltaEtaHist_Psi2S->Draw("SAME");

    c1->cd(2);
    TH1F *deltaPhiHist_JPsi = new TH1F("h3","Delta Phi JPsi;#Delta#varphi (rad);#frac{1}{N_{trig}} dN/d#Delta#varphi (rad^{-1})",n_bins,-0.5*M_PI, 3.0/2.0*M_PI);
    deltaPhiHist_JPsi->FillN(deltaPhi_JPsi.size(), deltaPhi_JPsi.data(), nullptr);
    deltaPhiHist_JPsi->SetLineColor(kGreen+1);
    deltaPhiHist_JPsi->SetLineWidth(2);
    deltaPhiHist_JPsi->Sumw2();
    deltaPhiHist_JPsi->Scale(1.0 / (triggers_JPsi*deltaPhiHist_JPsi->GetBinWidth(1)));
    deltaPhiHist_JPsi->SetMinimum(0);
    deltaPhiHist_JPsi->Draw();
    
    TH1F *deltaPhiHist_Psi2S = new TH1F("h4","Delta Phi Psi2S;#Delta#varphi (rad);#frac{1}{N_{trig}} dN/d#Delta#varphi (rad^{-1})",n_bins,-0.5*M_PI, 3.0/2.0*M_PI);
    deltaPhiHist_Psi2S->FillN(deltaPhi_Psi2S.size(), deltaPhi_Psi2S.data(), nullptr);
    deltaPhiHist_Psi2S->SetLineColor(kGreen+2);
    deltaPhiHist_Psi2S->SetLineWidth(2);
    deltaPhiHist_Psi2S->Sumw2();
    deltaPhiHist_Psi2S->Scale(1.0 / (triggers_Psi2S*deltaPhiHist_Psi2S->GetBinWidth(1)));
    deltaPhiHist_Psi2S->SetMinimum(0);
    deltaPhiHist_Psi2S->Draw("SAME");
    
    // Scale ymax
    deltaEtaHist_JPsi->SetMaximum(std::max(deltaEtaHist_JPsi->GetMaximum(), deltaEtaHist_Psi2S->GetMaximum()) * 1.2);
    deltaEtaHist_Psi2S->SetMaximum(std::max(deltaEtaHist_JPsi->GetMaximum(), deltaEtaHist_Psi2S->GetMaximum()) * 1.2);
    deltaPhiHist_JPsi->SetMaximum(std::max(deltaPhiHist_JPsi->GetMaximum(), deltaPhiHist_Psi2S->GetMaximum()) * 1.2);
    deltaPhiHist_Psi2S->SetMaximum(std::max(deltaPhiHist_JPsi->GetMaximum(), deltaPhiHist_Psi2S->GetMaximum()) * 1.2);

    //Adjust label positions
    deltaEtaHist_JPsi->GetYaxis()->SetTitleOffset(1.6);
    deltaPhiHist_JPsi->GetYaxis()->SetTitleOffset(1.8);

    // Adjust margins
    increasePadMargins(c1, 2);
    for (int i = 1; i <= 2; ++i) {
        c1->cd(i);
        gPad->SetLeftMargin(0.23); 
        gPad->SetRightMargin(0.05); 
        drawLabel(MC_name, 0.28, 0.85);
    }

    // Legend
    TLegend *legend = new TLegend(0.7,0.3,0.88,0.48);
    legend->SetBorderSize(0);
    legend->AddEntry(deltaEtaHist_JPsi, "J/#Psi", "l");
    legend->AddEntry(deltaEtaHist_Psi2S, "#Psi(2S)", "l");
    legend->Draw();


    TString out_name = TString::Format("results/%s/%s/deltaEtaDeltaPhi_JPsi_Psi2S_%s", MC_name.Data(), type.Data(), motherLabel.Data());
    out_name.ReplaceAll(".", "_");
    c1->SaveAs(out_name + ".png");
}