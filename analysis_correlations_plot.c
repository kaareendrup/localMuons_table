#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TString.h>
#include <TMath.h>
#include <Math/Vector4D.h>
#include <TH1F.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "utils/plots.c"
#include "utils/analysis.c"

void analysis_plot_data(json config, TFile* file, TString data_name, TString type) {

    // Load trigger counts for normalization
    TTree *triggersTree = nullptr;
    file->GetObject("Correlations", triggersTree);
    std::map<std::string, int> trigger_counts = get_trigger_counts(triggersTree);

    float signal_range_min = config["signal_range"]["min"];
    float signal_range_max = config["signal_range"]["max"];
    float background_range_min = config["background_range"]["min"];
    float background_range_max = config["background_range"]["max"];

    std::vector<double> pT_bins = config["hists"]["pT_bins"].get<std::vector<double>>();

    /////////// Plot invariant mass distribution of candidates ///////////
    TCanvas *c0 = new TCanvas("c0", "Invariant Mass of Muon Pairs", 800, 600);
    TH1F* invMassHist = (TH1F*) file->Get("All_-1_invMass");
    invMassHist->Draw();

    double yMinSingle = invMassHist->GetMinimum();
    double yMaxSingle = invMassHist->GetMaximum();

    // Draw shaded boxes for background regions
    TBox* box1s = new TBox(background_range_min, yMinSingle, signal_range_min, yMaxSingle);
    box1s->SetFillColorAlpha(kGray, 1); // semi-transparent
    box1s->SetFillStyle(3004);
    box1s->Draw("same");

    TBox* box2s = new TBox(signal_range_max, yMinSingle, background_range_max, yMaxSingle);
    box2s->SetFillColorAlpha(kGray, 1); // semi-transparent
    box2s->SetFillStyle(3004);
    box2s->Draw("same");

    // Style and save
    increasePadMargins(c0, 1);
    drawLabel_cuts(data_name, type, &config, 0.55, 0.89);

    TString out_name_invmass = TString::Format("results/%s/%s/JPsi_invariant_mass", data_name.Data(), type.Data());
    out_name_invmass.ReplaceAll(".", "_");
    c0->SaveAs(out_name_invmass + ".png");

    /////////// Plot histograms of deltaEta and deltaPhi ///////////
    TCanvas *c1 = new TCanvas("c1", "Delta Eta and Delta Phi", 1300, 600);
    c1->Divide(2,1);
    c1->cd(1);

    // Load and draw signal and background histograms for deltaEta
    TH1F* deltaEtaHistSig = (TH1F*) file->Get("All_-1_signal_deltaEta");
    TH1F* deltaEtaHistBkg = (TH1F*) file->Get("All_-1_background_deltaEta");

    drawHist(deltaEtaHistSig, "Delta Eta;#Delta#eta;#frac{1}{N_{trig}} dN/d#Delta#eta", kGreen+1, 1.0);
    drawHist(deltaEtaHistBkg, "Delta Eta;#Delta#eta;#frac{1}{N_{trig}} dN/d#Delta#eta", kBlack, 1.0, true);
    setMax({deltaEtaHistSig, deltaEtaHistBkg});

    // Load and draw signal and background histograms for deltaPhi
    c1->cd(2);
    TH1F* deltaPhiHistSig = (TH1F*) file->Get("All_-1_signal_deltaPhi");
    TH1F* deltaPhiHistBkg = (TH1F*) file->Get("All_-1_background_deltaPhi");

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
    legendCorr->AddEntry(deltaEtaHistSig, "Signal", "l");
    legendCorr->AddEntry(deltaEtaHistBkg, "Background", "l");
    legendCorr->Draw();
    c1->cd(1);
    drawLabel_cuts(data_name, type, &config, 0.55, 0.89);
    
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
    drawHist(deltaEtaHist_subtracted, "Subtracted Delta Eta;#Delta#eta;#frac{1}{N_{trig}} dN/d#Delta#eta", kBlue+1, 1.0 / (trigger_counts["All_-1"] * deltaEtaHist_subtracted->GetBinWidth(1)));
    setMax({deltaEtaHist_subtracted});

    // Subtract background from signal for deltaPhi
    c2->cd(2);
    TH1F *deltaPhiHist_subtracted = (TH1F*)deltaPhiHistSig->Clone("h6");
    deltaPhiHist_subtracted->Add(deltaPhiHistBkg, -1);
    drawHist(deltaPhiHist_subtracted, "Subtracted Delta Phi;#Delta#varphi (rad);#frac{1}{N_{trig}} dN/d#Delta#varphi (rad^{-1})", kBlue+1, 1.0 / (trigger_counts["All_-1"] * deltaPhiHist_subtracted->GetBinWidth(1)));
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
    drawLabel_cuts(data_name, type, &config, 0.55, 0.89);
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
        float ptmin = pT_bins[p];
        float ptmax = pT_bins[p + 1];
        TH1F* deltaPhiHistSig_pT = (TH1F*) file->Get(TString::Format("All_-1_signal_deltaPhi_pT_%.1f_%.1f", ptmin, ptmax));
        TH1F* deltaPhiHistBkg_pT = (TH1F*) file->Get(TString::Format("All_-1_background_deltaPhi_pT_%.1f_%.1f", ptmin, ptmax));
        
        if (!deltaPhiHistSig_pT || !deltaPhiHistBkg_pT) {std::cerr << "Could not retrieve histograms for pT range " << TString::Format("All_-1_signal_deltaPhi_pT_%.1f_%.1f", ptmin, ptmax) << "\n"; continue;}
        
        // Draw signal and background for this pT bin
        c3->cd(p + 1);
        drawHist(deltaPhiHistSig_pT, TString::Format("%.1f < p_{T} < %.1f GeV/c;#Delta#varphi (rad);#frac{1}{N_{trig}} dN/d#Delta#varphi (rad^{-1})", ptmin, ptmax), kGreen+1, 1.0);
        drawHist(deltaPhiHistBkg_pT, TString::Format("%.1f < p_{T} < %.1f GeV/c;#Delta#varphi (rad);#frac{1}{N_{trig}} dN/d#Delta#varphi (rad^{-1})", ptmin, ptmax), kBlack, 1.0, true);
        setMax({deltaPhiHistSig_pT, deltaPhiHistBkg_pT});

        // Draw subtracted histogram for this pT bin
        c4->cd(p + 1);
        TH1F *deltaPhiHist_subtracted_pT = (TH1F*)deltaPhiHistSig_pT->Clone(TString::Format("h6_pT_%g_%g", ptmin, ptmax));
        deltaPhiHist_subtracted_pT->Add(deltaPhiHistBkg_pT, -1);
        drawHist(deltaPhiHist_subtracted_pT, TString::Format("%.1f < p_{T} < %.1f GeV/c;#Delta#varphi (rad);#frac{1}{N_{trig}} dN/d#Delta#varphi (rad^{-1})", ptmin, ptmax), kBlue+1, 1.0 / (trigger_counts["All_-1"] * deltaPhiHist_subtracted_pT->GetBinWidth(1)));
        setMax({deltaPhiHist_subtracted_pT});

    }
    
    c3->SaveAs(out_name + "_deltaPhi_by_pT.png");
    c4->SaveAs(out_name + "_deltaPhi_by_pT_subtracted.png");
    gStyle->SetOptTitle(0);
}

void analysis_correlations_plot() {

    std::cout << std::fixed << std::setprecision(1);
    SetALICEStyle();

    ////////////////////////////////////////////////////////////////////
    ////            Load configuration, setup up filenames          ////
    ////////////////////////////////////////////////////////////////////
    std::ifstream jsonFile("localMuons_table/config/config_analysis.json");
    json config;
    jsonFile >> config;

    // Data
    std::string dataset_name = config["data_name"];
    std::string muon_type = config["muon_type"];
    TString data_name = TString::Format("%s_%s", dataset_name.c_str(), muon_type.c_str());

    // TString type = "gen";
    TString type = "reco";

    TFile* file = TFile::Open(TString::Format("results/%s/%s/analysis.root", data_name.Data(), type.Data()), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open file\n";
        return;
    }

    analysis_plot_data(config, file, data_name, type);

}