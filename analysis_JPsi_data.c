#include "setALICEStyle.c"

std::cout << std::fixed << std::setprecision(1);
SetALICEStyle();

struct triggerWithAssocs {
    ROOT::Math::PtEtaPhiMVector JPsiTrack;
    std::vector<ROOT::Math::PtEtaPhiMVector> associates;

    triggerWithAssocs(ROOT::Math::PtEtaPhiMVector track) : JPsiTrack(track) {};
};

void analysis_JPsi_data() {

    TString data_name = "DQ_data";
    TString type = "reco";
    
    TString data_file = "results/" + data_name + "/" + type + "/muonAOD.root";
    
    int n_bins = 20;
    float signal_range_min = 2.9;
    float signal_range_max = 3.3;
    
    float signal_width = signal_range_max - signal_range_min;
    float background_range_min = signal_range_min - (signal_width / 2.0);
    float background_range_max = signal_range_max + (signal_width / 2.0);
    
    // Load the dataframe keys
    TFile *file = TFile::Open(data_file);

    // Loop over all keys in the file (directories etc.)
    TIter nextKey(file->GetListOfKeys());
    TKey* key;

    std::vector<triggerWithAssocs> JPsiCandidates, backgroundCandidates;

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

        // Prepare to label
        ULong64_t fGlobalIndexAssoc;
        tree->SetBranchAddress("fGlobalIndexassoc", &fGlobalIndexAssoc);

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
            }

            // Store the indexes of the best candidate muon pair
            ROOT::Math::PtEtaPhiMVector JPsiCandidate(-9999,0,0,-1); 
            std::cout << "Inital candidate mass: " << JPsiCandidate.M() << std::endl;
            int idx_cand_1, idx_cand_2;

            // Find the muon pair with invariant mass closest to J/Psi mass
            for (size_t j = 0; j < muon_vectors.size(); ++j) {
                for (size_t k = j + 1; k < muon_vectors.size(); ++k) {
                    auto track = muon_vectors[j] + muon_vectors[k];

                    // Replace candidate if closer to J/Psi mass
                    if (std::abs(track.M() - 3.0969) < std::abs(JPsiCandidate.M() - 3.0969)) { // J/Psi mass ~3.0969 GeV/c^2
                        std::cout << "  New candidate mass: " << track.M() << std::endl;
                        JPsiCandidate = track;
                        idx_cand_1 = j;
                        idx_cand_2 = k;
                    }
                }
            }

            triggerWithAssocs candidate(JPsiCandidate);

            // Loop over muons again to find assocs
            for (size_t j = 0; j < muon_vectors.size(); ++j) {
                if (j == idx_cand_1 || j == idx_cand_2) continue; // Skip the trigger muons
                auto assocTrack = muon_vectors[j];
                candidate.associates.push_back(assocTrack);
            }
            
            JPsiCandidates.push_back(candidate);
        }

        tree->ResetBranchAddresses();
        delete tree;
    }

    std::vector<double> deltaEta_JPsi, deltaPhi_JPsi;
    std::vector<double> deltaEta_background, deltaPhi_background;
    int triggers_JPsi = 0, triggers_background = 0;

    for (auto& candidate : JPsiCandidates) {
        for (auto& assoc : candidate.associates) {
            double deltaPhi = candidate.JPsiTrack.Phi() - assoc.Phi();
            if (deltaPhi > 3.0 / 2.0 * M_PI) {
                deltaPhi -= 2.0 * M_PI;
            }
            if (deltaPhi < -0.5 * M_PI) {
                deltaPhi += 2.0 * M_PI;
            }

            if (candidate.JPsiTrack.M() > signal_range_min && candidate.JPsiTrack.M() < signal_range_max) {
                deltaEta_JPsi.push_back(candidate.JPsiTrack.Eta() - assoc.Eta());
                deltaPhi_JPsi.push_back(deltaPhi);
                triggers_JPsi++;
            } else if (candidate.JPsiTrack.M() > background_range_min && candidate.JPsiTrack.M() < signal_range_min || 
                       candidate.JPsiTrack.M() > signal_range_max && candidate.JPsiTrack.M() < background_range_max) {
                deltaEta_background.push_back(candidate.JPsiTrack.Eta() - assoc.Eta());
                deltaPhi_background.push_back(deltaPhi);
                triggers_background++;
            }
        }
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

    TH1F *deltaEtaHist_background = new TH1F("h2","Delta Eta Psi2S;#Delta#eta;#frac{1}{N_{trig}} dN/d#Delta#eta",n_bins,-5,3);
    deltaEtaHist_background->FillN(deltaEta_background.size(), deltaEta_background.data(), nullptr);
    deltaEtaHist_background->SetLineColor(kGreen+2);
    deltaEtaHist_background->SetLineWidth(2);
    deltaEtaHist_background->Sumw2();
    deltaEtaHist_background->Scale(1.0 / (triggers_background*deltaEtaHist_background->GetBinWidth(1)));
    deltaEtaHist_background->Draw("SAME");

    c1->cd(2);
    TH1F *deltaPhiHist_JPsi = new TH1F("h3","Delta Phi JPsi;#Delta#varphi (rad);#frac{1}{N_{trig}} dN/d#Delta#varphi (rad^{-1})",n_bins,-0.5*M_PI, 3.0/2.0*M_PI);
    deltaPhiHist_JPsi->FillN(deltaPhi_JPsi.size(), deltaPhi_JPsi.data(), nullptr);
    deltaPhiHist_JPsi->SetLineColor(kGreen+1);
    deltaPhiHist_JPsi->SetLineWidth(2);
    deltaPhiHist_JPsi->Sumw2();
    deltaPhiHist_JPsi->Scale(1.0 / (triggers_JPsi*deltaPhiHist_JPsi->GetBinWidth(1)));
    deltaPhiHist_JPsi->SetMinimum(0);
    deltaPhiHist_JPsi->Draw();
    
    TH1F *deltaPhiHist_background = new TH1F("h4","Delta Phi Psi2S;#Delta#varphi (rad);#frac{1}{N_{trig}} dN/d#Delta#varphi (rad^{-1})",n_bins,-0.5*M_PI, 3.0/2.0*M_PI);
    deltaPhiHist_background->FillN(deltaPhi_background.size(), deltaPhi_background.data(), nullptr);
    deltaPhiHist_background->SetLineColor(kGreen+2);
    deltaPhiHist_background->SetLineWidth(2);
    deltaPhiHist_background->Sumw2();
    deltaPhiHist_background->Scale(1.0 / (triggers_background*deltaPhiHist_background->GetBinWidth(1)));
    deltaPhiHist_background->SetMinimum(0);
    deltaPhiHist_background->Draw("SAME");
    
    // Scale ymax
    double max_JPsi = std::max(deltaEtaHist_JPsi->GetMaximum(), deltaPhiHist_JPsi->GetMaximum());
    double max_background = std::max(deltaEtaHist_background->GetMaximum(), deltaPhiHist_background->GetMaximum());
    double global_max = std::max(max_JPsi, max_background);
    deltaEtaHist_JPsi->SetMaximum(1.5 * global_max);
    deltaPhiHist_JPsi->SetMaximum(1.5 * global_max);
    deltaEtaHist_background->SetMaximum(1.5 * global_max);
    deltaPhiHist_background->SetMaximum(1.5 * global_max);

    //Adjust label positions
    deltaEtaHist_JPsi->GetYaxis()->SetTitleOffset(1.6);
    deltaPhiHist_JPsi->GetYaxis()->SetTitleOffset(1.8);

    // Adjust margins
    increasePadMargins(c1, 2);
    for (int i = 1; i <= 2; ++i) {
        c1->cd(i);
        gPad->SetLeftMargin(0.23); 
        gPad->SetRightMargin(0.05); 
        // drawLabel(data_name, 0.28, 0.85);
    }

    // Legend
    TLegend *legend = new TLegend(0.7,0.3,0.88,0.48);
    legend->SetBorderSize(0);
    legend->AddEntry(deltaEtaHist_JPsi, "J/#Psi", "l");
    legend->AddEntry(deltaEtaHist_background, "Background", "l");
    legend->Draw();


    TString out_name = TString::Format("results/%s/%s/deltaEtaDeltaPhi", data_name.Data(), type.Data());
    out_name.ReplaceAll(".", "_");
    c1->SaveAs(out_name + ".png");
}