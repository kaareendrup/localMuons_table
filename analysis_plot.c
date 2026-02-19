#include "setALICEStyle.c"

std::cout << std::fixed << std::setprecision(1);
SetALICEStyle();

std::map<std::string, int> get_trigger_counts(TTree* tree) {
    // Load trigger counts from the "Triggers" tree into a map for easy access when normalizing histograms
    std::map<std::string, int> trigger_counts;
    std::string *category = nullptr;
    int count;

    tree->SetBranchAddress("category", &category);
    tree->SetBranchAddress("count", &count);

    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {

        tree->GetEntry(i);
        std::string category_str = *category;
        trigger_counts[category_str] = count;
    }

    return trigger_counts;
}

void drawHist(TH1F* hist, TString title, int line_color, float scale_factor, bool same = false) {
    // Draw any histogram with consistent styling
    hist->SetLineColor(line_color);
    hist->SetLineWidth(2);
    hist->SetTitle(title);
    hist->Sumw2();
    hist->Scale(scale_factor);
    hist->SetMinimum(0);
    if (same) {
        hist->Draw("SAME");
    } else {
        hist->Draw();
    }
}

void setMax(std::vector<TH1F*> hists) {
    // Adjust y-axis maximum to be 1.5 times the largest maximum among the provided histograms
    double max_val = 0;
    for (auto hist : hists) {
        if (hist->GetMaximum() > max_val) {
            max_val = hist->GetMaximum();
        }
    }
    for (auto hist : hists) {
        hist->SetMaximum(1.5 * max_val);
    }
}

void analysis_plot_data(TFile* file, TString data_name, TString type) {

    // Load trigger counts for normalization
    TTree *triggersTree = nullptr;
    file->GetObject("Triggers", triggersTree);
    std::map<std::string, int> trigger_counts = get_trigger_counts(triggersTree);

    // Load signal ranges and pT sigments from the "SignalRanges" tree
    TTree *signalRangesTree = nullptr;
    file->GetObject("MetaData", signalRangesTree);
    std::vector<float> *signal_ranges = nullptr;
    std::vector<float> *segments_pt = nullptr;
    std::vector<float> *pTCuts = nullptr;
    std::vector<float> *etaCuts = nullptr;
    signalRangesTree->SetBranchAddress("signal_range", &signal_ranges);
    signalRangesTree->SetBranchAddress("segments_pt", &segments_pt);
    signalRangesTree->SetBranchAddress("pTCuts", &pTCuts);
    signalRangesTree->SetBranchAddress("etaCuts", &etaCuts);
    signalRangesTree->GetEntry(0); 

    /////////// Plot invariant mass distribution of candidates ///////////
    TCanvas *c0 = new TCanvas("c0", "Invariant Mass of Muon Pairs", 800, 600);
    TH1F* invMassHist = (TH1F*) file->Get("All_invMass");
    invMassHist->Draw();

    double yMinSingle = invMassHist->GetMinimum();
    double yMaxSingle = invMassHist->GetMaximum();

    // Draw shaded boxes for background regions
    TBox* box1s = new TBox(signal_ranges->at(0), yMinSingle, signal_ranges->at(1), yMaxSingle);
    box1s->SetFillColorAlpha(kGray, 1); // semi-transparent
    box1s->SetFillStyle(3004);
    box1s->Draw("same");

    TBox* box2s = new TBox(signal_ranges->at(2), yMinSingle, signal_ranges->at(3), yMaxSingle);
    box2s->SetFillColorAlpha(kGray, 1); // semi-transparent
    box2s->SetFillStyle(3004);
    box2s->Draw("same");    

    // Style and save
    increasePadMargins(c0, 1);
    drawLabel_cuts(data_name, type, pTCuts, etaCuts, 0.55, 0.89);

    TString out_name_invmass = TString::Format("results/%s/%s/JPsi_invariant_mass", data_name.Data(), type.Data());
    out_name_invmass.ReplaceAll(".", "_");
    c0->SaveAs(out_name_invmass + ".png");

    /////////// Plot histograms of deltaEta and deltaPhi ///////////
    TCanvas *c1 = new TCanvas("c1", "Delta Eta and Delta Phi", 1300, 600);
    c1->Divide(2,1);
    c1->cd(1);

    // Load and draw signal and background histograms for deltaEta
    TH1F* deltaEtaHistSig = (TH1F*) file->Get("All_signal_deltaEta");
    TH1F* deltaEtaHistBkg = (TH1F*) file->Get("All_background_deltaEta");

    drawHist(deltaEtaHistSig, "Delta Eta;#Delta#eta;#frac{1}{N_{trig}} dN/d#Delta#eta", kGreen+1, 1.0);
    drawHist(deltaEtaHistBkg, "Delta Eta;#Delta#eta;#frac{1}{N_{trig}} dN/d#Delta#eta", kBlack, 1.0, true);
    setMax({deltaEtaHistSig, deltaEtaHistBkg});

    // Load and draw signal and background histograms for deltaPhi
    c1->cd(2);
    TH1F* deltaPhiHistSig = (TH1F*) file->Get("All_signal_deltaPhi");
    TH1F* deltaPhiHistBkg = (TH1F*) file->Get("All_background_deltaPhi");

    drawHist(deltaPhiHistSig, "Delta Phi;#Delta#varphi (rad);#frac{1}{N_{trig}} dN/d#Delta#varphi (rad^{-1})", kGreen+1, 1.0);
    drawHist(deltaPhiHistBkg, "Delta Phi;#Delta#varphi (rad);#frac{1}{N_{trig}} dN/d#Delta#varphi (rad^{-1})", kBlack, 1.0, true);
    setMax({deltaPhiHistSig, deltaPhiHistBkg});

    //Adjust label positions
    deltaEtaHistSig->GetYaxis()->SetTitleOffset(1.6);
    deltaPhiHistSig->GetYaxis()->SetTitleOffset(1.8);

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
    legendCorr->AddEntry(deltaEtaHistSig, "Signal", "l");
    legendCorr->AddEntry(deltaEtaHistBkg, "Background", "l");
    legendCorr->Draw();
    c1->cd(1);
    drawLabel_cuts(data_name, type, pTCuts, etaCuts, 0.55, 0.89);
    
    TString out_name = TString::Format("results/%s/%s/deltaEtaDeltaPhi", data_name.Data(), type.Data());
    out_name.ReplaceAll(".", "_");
    c1->SaveAs(out_name + ".png");

    /////////// Make subtraction plots ///////////
    TCanvas *c2 = new TCanvas("c2", "Subtracted Delta Eta and Delta Phi", 1300, 600);
    c2->Divide(2,1);

    // Subtract background from signal for deltaEta
    c2->cd(1);
    TH1F *deltaEtaHist_subtracted = (TH1F*)deltaEtaHistSig->Clone("h5");
    deltaEtaHist_subtracted->Add(deltaEtaHistBkg, -1);
    drawHist(deltaEtaHist_subtracted, "Subtracted Delta Eta;#Delta#eta;#frac{1}{N_{trig}} dN/d#Delta#eta", kBlue+1, 1.0 / (trigger_counts["All"] * deltaEtaHist_subtracted->GetBinWidth(1)));
    setMax({deltaEtaHist_subtracted});

    // Subtract background from signal for deltaPhi
    c2->cd(2);
    TH1F *deltaPhiHist_subtracted = (TH1F*)deltaPhiHistSig->Clone("h6");
    deltaPhiHist_subtracted->Add(deltaPhiHistBkg, -1);
    drawHist(deltaPhiHist_subtracted, "Subtracted Delta Phi;#Delta#varphi (rad);#frac{1}{N_{trig}} dN/d#Delta#varphi (rad^{-1})", kBlue+1, 1.0 / (trigger_counts["All"] * deltaPhiHist_subtracted->GetBinWidth(1)));
    setMax({deltaPhiHist_subtracted});

    // Adjust margins
    increasePadMargins(c2, 2);
    for (int i = 1; i <= 2; ++i) {
        c2->cd(i);
        gPad->SetLeftMargin(0.23); 
        gPad->SetRightMargin(0.05); 
        // drawLabel(data_name, 0.28, 0.85);
    }
    c2->cd(1);
    drawLabel_cuts(data_name, type, pTCuts, etaCuts, 0.55, 0.89);
    c2->SaveAs(out_name + "_subtracted.png");

    /////////// Make signal/background and subtraction plots for pT bins ///////////
    TCanvas *c3 = new TCanvas("c3", "Signal/background Delta Phi by pT", 1300, 600);
    c3->Divide(3,3);
    TCanvas *c4 = new TCanvas("c4", "Subtracted Delta Phi by pT", 1300, 600);
    c4->Divide(3,3);

    // Enable title
    gStyle->SetOptTitle(1);

    // Loop over pT segments
    for (int p = 0; p < 9; ++p) {

        // Load the correct histograms for this pT bin
        float ptmin = segments_pt->at(p);
        float ptmax = segments_pt->at(p + 1);
        TH1F* deltaPhiHistSig_pT = (TH1F*) file->Get(TString::Format("All_signal_deltaPhi_pT_%.1f_%.1f", ptmin, ptmax));
        TH1F* deltaPhiHistBkg_pT = (TH1F*) file->Get(TString::Format("All_background_deltaPhi_pT_%.1f_%.1f", ptmin, ptmax));
        
        if (!deltaPhiHistSig_pT || !deltaPhiHistBkg_pT) {std::cerr << "Could not retrieve histograms for pT range " << TString::Format("All_signal_deltaPhi_pT_%.1f_%.1f", ptmin, ptmax) << "\n"; continue;}
        
        // Draw signal and background for this pT bin
        c3->cd(p + 1);
        drawHist(deltaPhiHistSig_pT, TString::Format("%.1f < p_{T} < %.1f GeV/c;#Delta#varphi (rad);#frac{1}{N_{trig}} dN/d#Delta#varphi (rad^{-1})", ptmin, ptmax), kGreen+1, 1.0);
        drawHist(deltaPhiHistBkg_pT, TString::Format("%.1f < p_{T} < %.1f GeV/c;#Delta#varphi (rad);#frac{1}{N_{trig}} dN/d#Delta#varphi (rad^{-1})", ptmin, ptmax), kBlack, 1.0, true);
        setMax({deltaPhiHistSig_pT, deltaPhiHistBkg_pT});

        // Draw subtracted histogram for this pT bin
        c4->cd(p + 1);
        TH1F *deltaPhiHist_subtracted_pT = (TH1F*)deltaPhiHistSig_pT->Clone(TString::Format("h6_pT_%g_%g", ptmin, ptmax));
        deltaPhiHist_subtracted_pT->Add(deltaPhiHistBkg_pT, -1);
        drawHist(deltaPhiHist_subtracted_pT, TString::Format("%.1f < p_{T} < %.1f GeV/c;#Delta#varphi (rad);#frac{1}{N_{trig}} dN/d#Delta#varphi (rad^{-1})", ptmin, ptmax), kBlue+1, 1.0 / (trigger_counts["All"] * deltaPhiHist_subtracted_pT->GetBinWidth(1)));
        setMax({deltaPhiHist_subtracted_pT});

    }
    
    c3->SaveAs(out_name + "_deltaPhi_by_pT.png");
    c4->SaveAs(out_name + "_deltaPhi_by_pT_subtracted.png");
    gStyle->SetOptTitle(0);
}

void analysis_plot_MC(TFile* file, TString data_name, TString type) {

    // Load signal ranges and pT sigments from the "SignalRanges" tree
    TTree *signalRangesTree = nullptr;
    file->GetObject("MetaData", signalRangesTree);
    std::vector<float> *pTCuts = nullptr;
    std::vector<float> *etaCuts = nullptr;
    signalRangesTree->SetBranchAddress("pTCuts", &pTCuts);
    signalRangesTree->SetBranchAddress("etaCuts", &etaCuts);
    signalRangesTree->GetEntry(0); 

    TTree *triggersTree = nullptr;
    file->GetObject("Triggers", triggersTree);
    std::map<std::string, int> trigger_counts = get_trigger_counts(triggersTree);
    
    // Import histograms split by charm/non-charm and J-Psi/Psi(2S)
    TH1F *deltaEtaHist_JPsi_charm = (TH1F*) file->Get("JPsi_charm_deltaEta");
    TH1F *deltaPhiHist_JPsi_charm = (TH1F*) file->Get("JPsi_charm_deltaPhi");
    TH1F *deltaEtaHist_JPsi_noncharm = (TH1F*) file->Get("JPsi_noncharm_deltaEta");
    TH1F *deltaPhiHist_JPsi_noncharm = (TH1F*) file->Get("JPsi_noncharm_deltaPhi");
    
    TH1F *deltaEtaHist_Psi2S_charm = (TH1F*) file->Get("Psi2S_charm_deltaEta");
    TH1F *deltaPhiHist_Psi2S_charm = (TH1F*) file->Get("Psi2S_charm_deltaPhi");
    TH1F *deltaEtaHist_Psi2S_noncharm = (TH1F*) file->Get("Psi2S_noncharm_deltaEta");
    TH1F *deltaPhiHist_Psi2S_noncharm = (TH1F*) file->Get("Psi2S_noncharm_deltaPhi");   
    
    // Clone histograms to create combined plots
    TH1F *deltaEtaHist_JPsi = (TH1F*)deltaEtaHist_JPsi_charm->Clone("h1");
    TH1F *deltaEtaHist_Psi2S = (TH1F*)deltaEtaHist_Psi2S_charm->Clone("h2");
    TH1F *deltaPhiHist_JPsi = (TH1F*)deltaPhiHist_JPsi_charm->Clone("h3");
    TH1F *deltaPhiHist_Psi2S = (TH1F*)deltaPhiHist_Psi2S_charm->Clone("h4");

    /////// PLot split by charm/non-charm and plot for J/Psi only ///////
    TCanvas *c1 = new TCanvas("c1", "Delta Eta and Delta Phi", 1300, 600);
    c1->Divide(2,1);
    
    // dEta
    c1->cd(1);
    drawHist(deltaEtaHist_JPsi_charm, "Delta Eta JPsi;#Delta#eta;#frac{1}{N_{trig}} dN/d#Delta#eta", kRed+1, 1.0 / (trigger_counts["JPsi"] * deltaEtaHist_JPsi_charm->GetBinWidth(1)));
    drawHist(deltaEtaHist_JPsi_noncharm, "Delta Eta JPsi;#Delta#eta;#frac{1}{N_{trig}} dN/d#Delta#eta", kBlack, 1.0 / (trigger_counts["JPsi"] * deltaEtaHist_JPsi_noncharm->GetBinWidth(1)), true);
    setMax({deltaEtaHist_JPsi_charm, deltaEtaHist_JPsi_noncharm});

    // dPhi
    c1->cd(2);
    drawHist(deltaPhiHist_JPsi_charm, "Delta Phi JPsi;#Delta#varphi (rad);#frac{1}{N_{trig}} dN/d#Delta#varphi (rad^{-1})", kRed+1, 1.0 / (trigger_counts["JPsi"] * deltaPhiHist_JPsi_charm->GetBinWidth(1)));
    drawHist(deltaPhiHist_JPsi_noncharm, "Delta Phi JPsi;#Delta#varphi (rad);#frac{1}{N_{trig}} dN/d#Delta#varphi (rad^{-1})", kBlack, 1.0 / (trigger_counts["JPsi"] * deltaPhiHist_JPsi_noncharm->GetBinWidth(1)), true);
    setMax({deltaPhiHist_JPsi_charm, deltaPhiHist_JPsi_noncharm});

    // Legend
    TLegend *legend1 = new TLegend(0.65,0.7,0.9,0.85);
    legend1->SetBorderSize(0);
    legend1->AddEntry(deltaEtaHist_JPsi_charm, "J/#psi - #mu from charm", "l");
    legend1->AddEntry(deltaEtaHist_JPsi_noncharm, "J/#psi - #mu from non-charm", "l");
    legend1->Draw();
    c1->cd(1);
    drawLabel_cuts(data_name, type, pTCuts, etaCuts, 0.55, 0.89);

    increasePadMargins(c1, 2);

    c1->SaveAs(TString::Format("results/%s/%s/JPsi_charm_noncharm.png", data_name.Data(), type.Data()));    

    /////// Plot split by J-Pps/Psi(2S) combining charm and non-charm ///////
    TCanvas *c2 = new TCanvas("c2", "Delta Eta and Delta Phi", 1300, 600);
    c2->Divide(2,1);

    // dEta, adding non-charm to charm to get total
    c2->cd(1);
    deltaEtaHist_JPsi->Add(deltaEtaHist_JPsi_noncharm, 1);
    drawHist(deltaEtaHist_JPsi, "Delta Eta JPsi;#Delta#eta;#frac{1}{N_{trig}} dN/d#Delta#eta", kGreen+1, 1.0 / ((trigger_counts["JPsi"]+trigger_counts["Psi2S"]) * deltaEtaHist_JPsi->GetBinWidth(1)));

    deltaEtaHist_Psi2S->Add(deltaEtaHist_Psi2S_noncharm, 1);
    drawHist(deltaEtaHist_Psi2S, "Delta Eta Psi2S;#Delta#eta;#frac{1}{N_{trig}} dN/d#Delta#eta", kGreen+2, 1.0 / ((trigger_counts["JPsi"]+trigger_counts["Psi2S"]) * deltaEtaHist_Psi2S->GetBinWidth(1)), true);
    setMax({deltaEtaHist_JPsi, deltaEtaHist_Psi2S});

    // dPhi, adding non-charm to charm to get total
    c2->cd(2);
    deltaPhiHist_JPsi->Add(deltaPhiHist_JPsi_noncharm, 1);
    drawHist(deltaPhiHist_JPsi, "Delta Phi JPsi;#Delta#varphi (rad);#frac{1}{N_{trig}} dN/d#Delta#varphi (rad^{-1})", kGreen+1, 1.0 / ((trigger_counts["JPsi"]+trigger_counts["Psi2S"]) * deltaPhiHist_JPsi->GetBinWidth(1)));

    deltaPhiHist_Psi2S->Add(deltaPhiHist_Psi2S_noncharm, 1);
    drawHist(deltaPhiHist_Psi2S, "Delta Phi Psi2S;#Delta#varphi (rad);#frac{1}{N_{trig}} dN/d#Delta#varphi (rad^{-1})", kGreen+2, 1.0 / ((trigger_counts["JPsi"]+trigger_counts["Psi2S"]) * deltaPhiHist_Psi2S->GetBinWidth(1)), true);
    setMax({deltaPhiHist_JPsi, deltaPhiHist_Psi2S});

    // Legend
    TLegend *legend2 = new TLegend(0.65,0.7,0.85,0.85);
    legend2->SetBorderSize(0);
    legend2->AddEntry(deltaEtaHist_JPsi, "J/#psi", "l");
    legend2->AddEntry(deltaEtaHist_Psi2S, "#psi(2S)", "l");
    legend2->Draw();
    c2->cd(1);
    drawLabel_cuts(data_name, type, pTCuts, etaCuts, 0.55, 0.89);

    increasePadMargins(c2, 2);

    c2->SaveAs(TString::Format("results/%s/%s/JPsi_all.png", data_name.Data(), type.Data()));
}

void analysis_plot() {

    TString data_name = "DQ";
    // TString data_name = "DQ_data";

    TString type = "gen";
    // TString type = "reco";

    TFile* file = TFile::Open(TString::Format("results/%s/%s/analysis.root", data_name.Data(), type.Data()), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open file\n";
        return;
    }

    if (data_name == "DQ_data") {
        analysis_plot_data(file, data_name, type);
    } else {
        analysis_plot_MC(file, data_name, type);
    }

}