#include "setALICEStyle.c"

std::cout << std::fixed << std::setprecision(1);
SetALICEStyle();

struct triggerWithAssocs {
    ROOT::Math::PtEtaPhiMVector JPsiTrack;
    std::vector<ROOT::Math::PtEtaPhiMVector> associates;

    triggerWithAssocs(ROOT::Math::PtEtaPhiMVector track) : JPsiTrack(track) {};
};

void analysis_JPsi_data() {

    // Data
    TString data_name = "DQ";
    // TString data_name = "DQ_data";
    // TString data_name = "DQ_data_global";
    TString type = "reco";
    
    TString data_file = "results/" + data_name + "/" + type + "/muonAOD.root";
    
    // Cuts
    // float pT_trigger_min = 1.0;
    float pT_trigger_min = 0.0;
    float pT_trigger_max = 20.0;
    // float pT_assoc_min = 3.0;
    float pT_assoc_min = 1.0;
    float pT_assoc_max = 20.0;

    // Signal region
    int n_bins = 20;
    float signal_range_min = 2.7;
    float signal_range_max = 3.4;
    
    float signal_width = signal_range_max - signal_range_min;
    float background_range_min = signal_range_min - (signal_width / 2.0);
    float background_range_max = signal_range_max + (signal_width / 2.0);
    
    // Load the dataframe keys
    TFile *file = TFile::Open(data_file);

    // Loop over all keys in the file (directories etc.)
    TIter nextKey(file->GetListOfKeys());
    TKey* key;

    std::vector<triggerWithAssocs> JPsiCandidates, backgroundCandidates;
    int dirCount = 0;

    // Loop over dataframes
    while ((key = (TKey*) nextKey())) {

        // Load directory and tree
        TObject* obj = key->ReadObj();
        if (!(obj->InheritsFrom("TDirectory"))) continue;
        
        TDirectory* dir = (TDirectory*) obj;
        TTree *tree = (TTree*)dir->Get("O2dqmuontable");

        printf("Reading tracks from dir %d of %d: %s\n", dirCount, file->GetListOfKeys()->GetEntries(), dir->GetName());

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
            // std::cout << "Inital candidate mass: " << JPsiCandidate.M() << std::endl;
            int idx_cand_1, idx_cand_2;

            // Find the muon pair with invariant mass closest to J/Psi mass
            for (size_t j = 0; j < muon_vectors.size(); ++j) {
                for (size_t k = j + 1; k < muon_vectors.size(); ++k) {
                    auto track = muon_vectors[j] + muon_vectors[k];

                    // Replace candidate if closer to J/Psi mass
                    if (std::abs(track.M() - 3.0969) < std::abs(JPsiCandidate.M() - 3.0969)) { // J/Psi mass ~3.0969 GeV/c^2
                        // std::cout << "  New candidate mass: " << track.M() << std::endl;
                        JPsiCandidate = track;
                        idx_cand_1 = j;
                        idx_cand_2 = k;
                    }
                }
            }

            // Create best candiddate from vector
            triggerWithAssocs candidate(JPsiCandidate);

            // Loop over muons again to find assocs
            for (size_t j = 0; j < muon_vectors.size(); ++j) {
                if (j == idx_cand_1 || j == idx_cand_2) continue; // Skip the trigger muons
                auto assocTrack = muon_vectors[j];
                candidate.associates.push_back(assocTrack);
            }
            
            // Store candidate and associates
            JPsiCandidates.push_back(candidate);
        }

        tree->ResetBranchAddresses();
        delete tree;

        dirCount++;
    }

    // Extract values for histograms
    std::vector<double> deltaEta_JPsi, deltaPhi_JPsi;
    std::vector<double> deltaEta_background, deltaPhi_background;
    int triggers_JPsi = 0, triggers_background = 0, triggers_total = 0;
    std::vector<double> candidate_inv_masses;
    std::vector<double> trigger_pT, assoc_pT;

    for (auto& candidate : JPsiCandidates) {

        if ((candidate.JPsiTrack.Pt() < pT_trigger_min) || (candidate.JPsiTrack.Pt() > pT_trigger_max)) continue; 
        candidate_inv_masses.push_back(candidate.JPsiTrack.M());
        trigger_pT.push_back(candidate.JPsiTrack.Pt());

        for (auto& assoc : candidate.associates) {
            if ((assoc.Pt() < pT_assoc_min) || (assoc.Pt() > pT_assoc_max)) continue;
            assoc_pT.push_back(assoc.Pt());

            // Calculate deltaPhi
            double deltaPhi = candidate.JPsiTrack.Phi() - assoc.Phi();
            if (deltaPhi > 3.0 / 2.0 * M_PI) {
                deltaPhi -= 2.0 * M_PI;
            }
            if (deltaPhi < -0.5 * M_PI) {
                deltaPhi += 2.0 * M_PI;
            }

            // Add to the appropriate histogram based on candidate mass
            if (candidate.JPsiTrack.M() > signal_range_min && candidate.JPsiTrack.M() < signal_range_max) {
                deltaEta_JPsi.push_back(candidate.JPsiTrack.Eta() - assoc.Eta());
                deltaPhi_JPsi.push_back(deltaPhi);
            } else if (candidate.JPsiTrack.M() > background_range_min && candidate.JPsiTrack.M() < signal_range_min || 
                       candidate.JPsiTrack.M() > signal_range_max && candidate.JPsiTrack.M() < background_range_max) {
                deltaEta_background.push_back(candidate.JPsiTrack.Eta() - assoc.Eta());
                deltaPhi_background.push_back(deltaPhi);
            }
        }

        // Count triggers
        if (candidate.JPsiTrack.M() > signal_range_min && candidate.JPsiTrack.M() < signal_range_max) {
            triggers_JPsi++;
        } else if (candidate.JPsiTrack.M() > background_range_min && candidate.JPsiTrack.M() < signal_range_min || 
                    candidate.JPsiTrack.M() > signal_range_max && candidate.JPsiTrack.M() < background_range_max) {
            triggers_background++;
        }
        triggers_total++;

    }

    /////////// Plot invariant mass distribution of candidates ///////////
    TCanvas *c0 = new TCanvas("c0", "Invariant Mass of Muon Pairs", 800, 600);
    TH1F *invMassHist = new TH1F("h0","Invariant Mass of Muon Pairs;Invariant Mass (GeV/c^{2});Counts",100,1.0,5.0);
    invMassHist->FillN(candidate_inv_masses.size(), candidate_inv_masses.data(), nullptr);
    invMassHist->Draw();

    double yMinSingle = invMassHist->GetMinimum();
    double yMaxSingle = invMassHist->GetMaximum();

    TBox* box1s = new TBox(background_range_min, yMinSingle, signal_range_min, yMaxSingle);
    box1s->SetFillColorAlpha(kGray, 1); // semi-transparent
    box1s->SetFillStyle(3004);
    box1s->Draw("same");

    TBox* box2s = new TBox(signal_range_max, yMinSingle, background_range_max, yMaxSingle);
    box2s->SetFillColorAlpha(kGray, 1); // semi-transparent
    box2s->SetFillStyle(3004);
    box2s->Draw("same");    

    // Adjust margins
    increasePadMargins(c0, 1);

    drawLabel(data_name, "", 0.4, 0.89);

    TString out_name_invmass = TString::Format("results/%s/%s/JPsi_invariant_mass", data_name.Data(), type.Data());
    out_name_invmass.ReplaceAll(".", "_");
    c0->SaveAs(out_name_invmass + ".png");

    /////////// Plot pT distributions of triggers and associates ///////////
    TCanvas *c3 = new TCanvas("c3", "pT of triggers and associates", 1300, 600);
    c3->Divide(2,1);

    c3->cd(1);
    TH1F *triggerpTHist = new TH1F("h7","pT of Trigger J/Psi;pT (GeV/c);Counts",20,0,20);
    triggerpTHist->FillN(trigger_pT.size(), trigger_pT.data(), nullptr);
    triggerpTHist->SetLineWidth(2);
    triggerpTHist->SetLineColor(kRed);
    triggerpTHist->Draw();
    gPad->SetLogy();
    
    c3->cd(2);
    TH1F *assocpTHist = new TH1F("h8","pT of Associated Muons;pT (GeV/c);Counts",20,0,20);
    assocpTHist->FillN(assoc_pT.size(), assoc_pT.data(), nullptr);
    assocpTHist->SetLineWidth(2);
    assocpTHist->SetLineColor(kBlue);
    assocpTHist->Draw();
    gPad->SetLogy();

    increaseMargins(c3);
    increasePadMargins(c3, 2);
    triggerpTHist->SetMaximum(1000 * triggerpTHist->GetMaximum());
    assocpTHist->SetMaximum(1000 * assocpTHist->GetMaximum());
    triggerpTHist->SetMinimum(1);
    assocpTHist->SetMinimum(1);

    // Legend
    TLegend *legendpT = new TLegend(0.6,0.8,0.7,0.9);
    legendpT->SetBorderSize(0);
    legendpT->SetFillStyle(0);
    legendpT->AddEntry(triggerpTHist, "Trigger", "l");
    legendpT->AddEntry(assocpTHist, "Assoc", "l");
    legendpT->Draw();

    drawLabel(data_name, "", 0.55, 0.89);

    TString out_name_pT = TString::Format("results/%s/%s/JPsi_assoc_pT", data_name.Data(), type.Data());
    out_name_pT.ReplaceAll(".", "_");
    c3->SaveAs(out_name_pT + ".png");

    /////////// Plot histograms of deltaEta and deltaPhi for JPsi and Psi2S ///////////
    TCanvas *c1 = new TCanvas("c1", "Delta Eta and Delta Phi", 1300, 600);
    c1->Divide(2,1);
    c1->cd(1);

    // TH1F *deltaEtaHist_JPsi = new TH1F("h1","Delta Eta JPsi;#Delta#eta;Counts",n_bins,-5,3);
    // TH1F *deltaEtaHist_JPsi = new TH1F("h1","Delta Eta JPsi;#Delta#eta;#frac{1}{N_{trig}} dN/d#Delta#eta",n_bins,-5,3);
    TH1F *deltaEtaHist_JPsi = new TH1F("h1","Delta Eta JPsi;#Delta#eta;Counts",n_bins,-5,3);
    deltaEtaHist_JPsi->FillN(deltaEta_JPsi.size(), deltaEta_JPsi.data(), nullptr);
    deltaEtaHist_JPsi->SetLineColor(kGreen+1);
    deltaEtaHist_JPsi->SetLineWidth(2);
    deltaEtaHist_JPsi->Sumw2();
    // deltaEtaHist_JPsi->Scale(1.0 / (triggers_JPsi*deltaEtaHist_JPsi->GetBinWidth(1)));
    deltaEtaHist_JPsi->Draw();

    // TH1F *deltaEtaHist_background = new TH1F("h2","Delta Eta Psi2S;#Delta#eta;#frac{1}{N_{trig}} dN/d#Delta#eta",n_bins,-5,3);
    TH1F *deltaEtaHist_background = new TH1F("h2","Delta Eta Psi2S;#Delta#eta;Counts",n_bins,-5,3);
    deltaEtaHist_background->FillN(deltaEta_background.size(), deltaEta_background.data(), nullptr);
    deltaEtaHist_background->SetLineColor(kBlack);
    deltaEtaHist_background->SetLineWidth(2);
    deltaEtaHist_background->Sumw2();
    // deltaEtaHist_background->Scale(1.0 / (triggers_background*deltaEtaHist_background->GetBinWidth(1)));
    deltaEtaHist_background->Draw("SAME");

    c1->cd(2);
    // TH1F *deltaPhiHist_JPsi = new TH1F("h3","Delta Phi JPsi;#Delta#varphi (rad);#frac{1}{N_{trig}} dN/d#Delta#varphi (rad^{-1})",n_bins,-0.5*M_PI, 3.0/2.0*M_PI);
    TH1F *deltaPhiHist_JPsi = new TH1F("h3","Delta Phi JPsi;#Delta#varphi (rad);Counts",n_bins,-0.5*M_PI, 3.0/2.0*M_PI);
    deltaPhiHist_JPsi->FillN(deltaPhi_JPsi.size(), deltaPhi_JPsi.data(), nullptr);
    deltaPhiHist_JPsi->SetLineColor(kGreen+1);
    deltaPhiHist_JPsi->SetLineWidth(2);
    deltaPhiHist_JPsi->Sumw2();
    // deltaPhiHist_JPsi->Scale(1.0 / (triggers_JPsi*deltaPhiHist_JPsi->GetBinWidth(1)));
    deltaPhiHist_JPsi->SetMinimum(0);
    deltaPhiHist_JPsi->Draw();
    
    // TH1F *deltaPhiHist_background = new TH1F("h4","Delta Phi Psi2S;#Delta#varphi (rad);#frac{1}{N_{trig}} dN/d#Delta#varphi (rad^{-1})",n_bins,-0.5*M_PI, 3.0/2.0*M_PI);
    TH1F *deltaPhiHist_background = new TH1F("h4","Delta Phi Psi2S;#Delta#varphi (rad);Counts",n_bins,-0.5*M_PI, 3.0/2.0*M_PI);
    deltaPhiHist_background->FillN(deltaPhi_background.size(), deltaPhi_background.data(), nullptr);
    deltaPhiHist_background->SetLineColor(kBlack);
    deltaPhiHist_background->SetLineWidth(2);
    deltaPhiHist_background->Sumw2();
    // deltaPhiHist_background->Scale(1.0 / (triggers_background*deltaPhiHist_background->GetBinWidth(1)));
    deltaPhiHist_background->SetMinimum(0);
    deltaPhiHist_background->Draw("SAME");
    
    // Scale ymax
    double max_Eta = std::max(deltaEtaHist_JPsi->GetMaximum(), deltaEtaHist_background->GetMaximum());
    double max_Phi = std::max(deltaPhiHist_JPsi->GetMaximum(), deltaPhiHist_background->GetMaximum());
    deltaEtaHist_JPsi->SetMaximum(1.5 * max_Eta);
    deltaPhiHist_JPsi->SetMaximum(1.5 * max_Phi);
    deltaEtaHist_background->SetMaximum(1.5 * max_Eta);
    deltaPhiHist_background->SetMaximum(1.5 * max_Phi);
    deltaEtaHist_JPsi->SetMinimum(0);
    deltaPhiHist_JPsi->SetMinimum(0);
    deltaEtaHist_background->SetMinimum(0);
    deltaPhiHist_background->SetMinimum(0);

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
    TLegend *legendCorr = new TLegend(0.65,0.7,0.9,0.85);
    legendCorr->SetBorderSize(0);
    legendCorr->AddEntry(deltaEtaHist_JPsi, "Signal", "l");
    legendCorr->AddEntry(deltaEtaHist_background, "Background", "l");
    legendCorr->Draw();
    drawLabel(data_name, "", 0.55, 0.89);

    TString out_name = TString::Format("results/%s/%s/deltaEtaDeltaPhi", data_name.Data(), type.Data());
    out_name.ReplaceAll(".", "_");
    c1->SaveAs(out_name + ".png");

    /////////// Make subtraction plots ///////////
    TCanvas *c2 = new TCanvas("c2", "Subtracted Delta Eta and Delta Phi", 1300, 600);
    c2->Divide(2,1);
    c2->cd(1);
    TH1F *deltaEtaHist_subtracted = (TH1F*)deltaEtaHist_JPsi->Clone("h5");
    deltaEtaHist_subtracted->SetTitle("Subtracted Delta Eta;#Delta#eta;#frac{1}{N_{trig}} dN/d#Delta#eta");
    deltaEtaHist_subtracted->Add(deltaEtaHist_background, -1);
    deltaEtaHist_subtracted->Sumw2();
    deltaEtaHist_subtracted->SetLineColor(kBlue+1);
    deltaEtaHist_subtracted->Scale(1.0 / (triggers_total * deltaEtaHist_subtracted->GetBinWidth(1)));
    deltaEtaHist_subtracted->Draw();

    c2->cd(2);
    TH1F *deltaPhiHist_subtracted = (TH1F*)deltaPhiHist_JPsi->Clone("h6");
    deltaPhiHist_subtracted->SetTitle("Subtracted Delta Phi;#Delta#varphi (rad);#frac{1}{N_{trig}} dN/d#Delta#varphi (rad^{-1})");
    deltaPhiHist_subtracted->Add(deltaPhiHist_background, -1);
    deltaPhiHist_subtracted->Sumw2();
    deltaPhiHist_subtracted->SetLineColor(kBlue+1);
    deltaPhiHist_subtracted->Scale(1.0 / (triggers_total * deltaPhiHist_subtracted->GetBinWidth(1)));
    deltaPhiHist_subtracted->Draw();

    // Scale ymax
    deltaEtaHist_subtracted->SetMaximum(1.5 * deltaEtaHist_subtracted->GetMaximum());
    deltaPhiHist_subtracted->SetMaximum(1.5 * deltaPhiHist_subtracted->GetMaximum());
    deltaEtaHist_subtracted->SetMinimum(0);
    deltaPhiHist_subtracted->SetMinimum(0);

    // Adjust margins
    increasePadMargins(c2, 2);
    for (int i = 1; i <= 2; ++i) {
        c2->cd(i);
        gPad->SetLeftMargin(0.23); 
        gPad->SetRightMargin(0.05); 
        // drawLabel(data_name, 0.28, 0.85);
    }
    
    drawLabel(data_name, "", 0.55, 0.89);

    c2->SaveAs(out_name + "_subtracted.png");
}