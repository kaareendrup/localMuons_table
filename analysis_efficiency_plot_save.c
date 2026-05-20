
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TString.h>
#include <TMath.h>
#include <Math/Vector4D.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TGraphErrors.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

#include "utils/plots.c"
#include <nlohmann/json.hpp>
using json = nlohmann::json;

void analysis_efficiency_plot_save() {

    std::cout << std::fixed << std::setprecision(1);
    SetALICEStyle();

    std::ifstream jsonFile("localMuons_table/config/config_analysis.json");
    json config;
    jsonFile >> config;

    // Data
    std::string data_name = config["data_name"];
    std::string muon_type = config["muon_type"];
    TString MC_name = TString::Format("%s_%s", data_name.c_str(), muon_type.c_str());

    TString in_file = TString::Format("results/%s/efficiency.root", MC_name.Data());
    TFile* file = TFile::Open(in_file, "READ");

    // Efficiency plots
    TH1F *muonEffHist = (TH1F*)file->Get("muonEffHist");
    TH1F *muonEffTrueHist = (TH1F*)file->Get("muonEffTrueHist");
    TH1F *JPsiEffHist = (TH1F*)file->Get("JPsiEffHist");
    TH1F *JPsiEffTrueHist = (TH1F*)file->Get("JPsiEffTrueHist");

    TCanvas *c1 = new TCanvas("c1", "Muon and J/Psi eff", 1300, 600);
    c1->Divide(2,1);

    gStyle->SetOptTitle(1);
    c1->cd(1);
    
    drawHist(muonEffHist, "#mu Efficiency;p_{T} (GeV/c);A #times #epsilon_{#mu}", kBlue, 1.0, false, true);
    drawLabel_cuts(MC_name, "", &config, 0.85, 0.59);

    c1->cd(2);
    drawHist(JPsiEffHist, "J/#Psi Efficiency;p_{T} (GeV/c);A #times #epsilon_{J/#Psi}", kBlue, 1.0, false, true);
    drawLabel_cuts(MC_name, "", &config, 0.85, 0.59);
    // drawLabel_cuts(MC_name, "", &config, 0.2, 0.85, true, false);
    gPad->SetLeftMargin(0.15);

    c1->SaveAs(TString::Format("results/%s/efficiency_plots.png", MC_name.Data()));

    // Efficiency plots (true pT)
    TCanvas *c2 = new TCanvas("c2", "Muon and J/Psi eff", 1300, 600);
    c2->Divide(2,1);

    c2->cd(1);
    drawHist(muonEffTrueHist, "#mu Efficiency;p_{T} (GeV/c);A #times #epsilon_{#mu}", kBlue, 1.0, false, true);
    drawLabel_cuts(MC_name, "", &config, 0.85, 0.59);

    c2->cd(2);
    drawHist(JPsiEffTrueHist, "J/#Psi Efficiency;p_{T} (GeV/c);A #times #epsilon_{J/#Psi}", kBlue, 1.0, false, true);
    drawLabel_cuts(MC_name, "", &config, 0.85, 0.59);

    c2->SaveAs(TString::Format("results/%s/efficiency_plots_true.png", MC_name.Data()));

    // pT spectrum plots
    TH1F *pTMuonRecoHist = (TH1F*)file->Get("pTMuonRecoHist");
    TH1F *pTMuonRecoTrueHist = (TH1F*)file->Get("pTMuonRecoTrueHist");
    TH1F *pTJPsiRecoHist = (TH1F*)file->Get("pTJPsiRecoHist");
    TH1F *pTJPsiRecoTrueHist = (TH1F*)file->Get("pTJPsiRecoTrueHist");
    TH1F *pTMuonGenHist = (TH1F*)file->Get("pTMuonGenHist");
    TH1F *pTJPsiGenHist = (TH1F*)file->Get("pTJPsiGenHist");
    
    pTMuonRecoTrueHist->Scale(1, "width");
    pTMuonRecoHist->Scale(1, "width");
    pTJPsiRecoTrueHist->Scale(1, "width");
    pTJPsiRecoHist->Scale(1, "width");
    pTMuonGenHist->Scale(1, "width");
    pTJPsiGenHist->Scale(1, "width");
    
    TCanvas *c3 = new TCanvas("c3", "Muon and J/Psi p_{T}", 1300, 600);
    c3->Divide(2,1);
    
    c3->cd(1);
    drawHist(pTMuonRecoTrueHist, "#mu p_{T};p_{T} (GeV/c);Counts", kRed, 1.0);
    drawHist(pTMuonRecoHist, "#mu p_{T};p_{T} (GeV/c);Counts", kBlue, 1.0, true);
    TLegend *legend1 = new TLegend(0.7, 0.7, 0.85, 0.85);
    legend1->AddEntry(pTMuonRecoTrueHist, "true p_{T}", "l");
    legend1->AddEntry(pTMuonRecoHist, "reco p_{T}", "l");
    legend1->Draw();
    drawLabel_cuts(MC_name, "", &config, 0.85, 0.59);
    
    c3->cd(2);
    drawHist(pTJPsiRecoTrueHist, "J/#Psi p_{T};p_{T} (GeV/c);Counts", kRed, 1.0);
    drawHist(pTJPsiRecoHist, "J/#Psi p_{T};p_{T} (GeV/c);Counts", kBlue, 1.0, true);
    TLegend *legend2 = new TLegend(0.7, 0.7, 0.85, 0.85);
    legend2->AddEntry(pTJPsiRecoTrueHist, "true p_{T}", "l");
    legend2->AddEntry(pTJPsiRecoHist, "reco p_{T}", "l");
    legend2->Draw();
    drawLabel_cuts(MC_name, "", &config, 0.85, 0.59);
    
    c3->SaveAs(TString::Format("results/%s/pT_plots_true_reco.png", MC_name.Data()));
    
    TCanvas *c4 = new TCanvas("c4", "Muon and J/Psi p_{T}", 1300, 600);
    c4->Divide(2,1);
    
    c4->cd(1);
    drawHist(pTMuonGenHist, "#mu p_{T};p_{T} (GeV/c);Counts", kRed, 1.0);
    drawHist(pTMuonRecoHist, "#mu p_{T};p_{T} (GeV/c);Counts", kBlue, 1.0, true);
    TLegend *legend3 = new TLegend(0.7, 0.7, 0.85, 0.85);
    legend3->AddEntry(pTMuonGenHist, "gen p_{T}", "l");
    legend3->AddEntry(pTMuonRecoHist, "reco p_{T}", "l");
    legend3->Draw();
    drawLabel_cuts(MC_name, "", &config, 0.85, 0.59);

    c4->cd(2);
    drawHist(pTJPsiGenHist, "J/#Psi p_{T};p_{T} (GeV/c);Counts", kRed, 1.0);
    drawHist(pTJPsiRecoHist, "J/#Psi p_{T};p_{T} (GeV/c);Counts", kBlue, 1.0, true);
    TLegend *legend4 = new TLegend(0.7, 0.7, 0.85, 0.85);
    legend4->AddEntry(pTJPsiGenHist, "gen p_{T}", "l");
    legend4->AddEntry(pTJPsiRecoHist, "reco p_{T}", "l");
    legend4->Draw();
    drawLabel_cuts(MC_name, "", &config, 0.85, 0.59);

    c4->SaveAs(TString::Format("results/%s/pT_plots_gen_reco.png", MC_name.Data()));
    
    // Resolution plots
    TString MC_name_1 = TString::Format("%s_%s", data_name.c_str(), "standalone");
    TString MC_name_2 = TString::Format("%s_%s", data_name.c_str(), "global");

    TString in_file_1 = TString::Format("results/%s/efficiency.root", MC_name_1.Data());
    TFile* file_1 = TFile::Open(in_file_1, "READ");
    TString in_file_2 = TString::Format("results/%s/efficiency.root", MC_name_2.Data());
    TFile* file_2 = TFile::Open(in_file_2, "READ");

    TH1F *pTMuonRecoTrueRecoHist_1 = (TH1F*)file_1->Get("pTMuonRecoTrueRecoHist");
    TH1F *pTMuonRecoTrueRecoHist_2 = (TH1F*)file_2->Get("pTMuonRecoTrueRecoHist");

    // Check files and histograms
    if (!file_1 || file_1->IsZombie() || !pTMuonRecoTrueRecoHist_1) {
        std::cerr << "Error: Could not open file or find histogram for " << MC_name_1.Data() << std::endl;
        return;
    }
    if (!file_2 || file_2->IsZombie() || !pTMuonRecoTrueRecoHist_2) {
        std::cerr << "Error: Could not open file or find histogram for " << MC_name_2.Data() << std::endl;
        return;
    }

    TCanvas *c5 = new TCanvas("c5", "Muon pT resolution", 700, 600);
    drawHist(pTMuonRecoTrueRecoHist_1, "Muon p_{T} resolution;|p_{T}^{true} - p_{T}^{reco}|/p_{T}^{true};Counts", kBlue, 1.2);
    drawHist(pTMuonRecoTrueRecoHist_2, "Muon p_{T} resolution;|p_{T}^{true} - p_{T}^{reco}|/p_{T}^{true};Counts", kRed, 1.2, true);

    // Add legend
    TLegend *legend_res = new TLegend(0.22, 0.7, 0.55, 0.85);
    legend_res->AddEntry(pTMuonRecoTrueRecoHist_1, "Standalone Muons", "l");
    legend_res->AddEntry(pTMuonRecoTrueRecoHist_2, "Global Muons", "l");
    legend_res->Draw();
    setMax({pTMuonRecoTrueRecoHist_1, pTMuonRecoTrueRecoHist_2});

    pTMuonRecoTrueRecoHist_1->SetMinimum(1e3);
    pTMuonRecoTrueRecoHist_2->SetMinimum(1e3);
    gPad->SetLogy();

    drawLabel_cuts("DQ", "", &config, 0.85, 0.8);

    c5->SaveAs(TString::Format("results/%s/pt_resolution_overlay.png", MC_name_1.Data()));

    TGraphErrors* gr_1 = (TGraphErrors*)file_1->Get("gr_ptResolution");
    TGraphErrors* gr_2 = (TGraphErrors*)file_2->Get("gr_ptResolution");

    TCanvas *c6 = new TCanvas("c6", "Muon pT resolution vs pT", 700, 600);
    drawGraph(gr_1, "Muon p_{T} resolution vs p_{T};p_{T} (GeV/c);#bar{|p_{T}^{true} - p_{T}^{reco}|/p_{T}^{true}}", kBlue);
    drawGraph(gr_2, "Muon p_{T} resolution vs p_{T};p_{T} (GeV/c);#bar{|p_{T}^{true} - p_{T}^{reco}|/p_{T}^{true}}", kRed, true);

    // Add legend
    TLegend *legend_res_pt = new TLegend(0.42, 0.7, 0.75, 0.85);
    legend_res_pt->AddEntry(gr_1, "Standalone Muons", "l");
    legend_res_pt->AddEntry(gr_2, "Global Muons", "l");
    legend_res_pt->Draw();
    setMax({gr_1->GetHistogram(), gr_2->GetHistogram()});
    gr_1->GetXaxis()->SetLimits(1, 20);
    gr_2->GetXaxis()->SetLimits(1, 20);
    increaseMargins(c6);
    c6->SetLeftMargin(0.20);
    c6->SetTopMargin(0.10);

    TGraphErrors* gr_3 = (TGraphErrors*)file_1->Get("gr_ptMean");
    TGraphErrors* gr_4 = (TGraphErrors*)file_2->Get("gr_ptMean");

    c6->SaveAs(TString::Format("results/%s/pt_resolution_vs_pt_overlay.png", MC_name_1.Data()));

    TCanvas *c7 = new TCanvas("c7", "Muon pT diff mean vs pT", 700, 600);
    drawGraph(gr_3, "Muon p_{T} resolution mean vs p_{T};p_{T} (GeV/c);#bar{|p_{T}^{true} - p_{T}^{reco}|/p_{T}^{true}}", kBlue);
    drawGraph(gr_4, "Muon p_{T} resolution mean vs p_{T};p_{T} (GeV/c);#bar{|p_{T}^{true} - p_{T}^{reco}|/p_{T}^{true}}", kRed, true);

    // Add legend
    TLegend *legend_res_pt_mean = new TLegend(0.22, 0.7, 0.55, 0.85);
    legend_res_pt_mean->AddEntry(gr_3, "Standalone Muons", "l");
    legend_res_pt_mean->AddEntry(gr_4, "Global Muons", "l");
    legend_res_pt_mean->Draw();
    setMax({gr_3->GetHistogram(), gr_4->GetHistogram()});
    gr_3->GetXaxis()->SetLimits(1, 20);
    gr_4->GetXaxis()->SetLimits(1, 20);
    increaseMargins(c7);
    c7->SetLeftMargin(0.20);
    c7->SetTopMargin(0.10);

    drawLabel_cuts("DQ", "", &config, 0.9, 0.8);

    c7->SaveAs(TString::Format("results/%s/mean_pt_resolution_vs_pt_overlay.png", MC_name_1.Data()));
}