
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

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "utils/plots.c"

void analysis_efficiency_plot_save_overlay() {

    std::cout << std::fixed << std::setprecision(1);
    SetALICEStyle();

    std::ifstream jsonFile("localMuons_table/config/config_analysis.json");
    json config;
    jsonFile >> config;

    float pT_JPsi_min = config["cuts_JPsi"]["pT_JPsi_min"];
    float pT_JPsi_max = config["cuts_JPsi"]["pT_JPsi_max"];

    float eta_JPsi_min = config["cuts_JPsi"]["eta_JPsi_min"];
    float eta_JPsi_max = config["cuts_JPsi"]["eta_JPsi_max"];
    std::vector<float> pT_JPsi_cuts = {pT_JPsi_min, pT_JPsi_max, 0, 20};
    std::vector<float> eta_JPsi_cuts = {eta_JPsi_min, eta_JPsi_max, -4.0, 4.0};

    // Data
    std::string data_name_1 = "c3";
    std::string data_name_2 = "f4d";
    std::string muon_type = config["muon_type"];
    TString MC_name_1 = TString::Format("%s_%s", data_name_1.c_str(), muon_type.c_str());
    TString MC_name_2 = TString::Format("%s_%s", data_name_2.c_str(), muon_type.c_str());

    TString in_file_1 = TString::Format("results/%s/efficiency.root", MC_name_1.Data());
    TFile* file_1 = TFile::Open(in_file_1, "READ");
    TString in_file_2 = TString::Format("results/%s/efficiency.root", MC_name_2.Data());
    TFile* file_2 = TFile::Open(in_file_2, "READ");

    // Efficiency plots
    TH1F *JPsiEffHist_1 = (TH1F*)file_1->Get("JPsiEffHist");
    TH1F *JPsiEffHist_2 = (TH1F*)file_2->Get("JPsiEffHist");

    TCanvas *c1 = new TCanvas("c1", "J/Psi eff", 700, 600);
    gStyle->SetOptTitle(1);
    // JPsiEffHist_1->GetXaxis()->SetRangeUser(0, 10);
    
    drawHist(JPsiEffHist_1, "J/#Psi Efficiency;p_{T} (GeV/c);Eff_{J/#Psi}", kBlue, 1.2);
    drawHist(JPsiEffHist_2, "J/#Psi Efficiency;p_{T} (GeV/c);Eff_{J/#Psi}", kRed, 1.2, true);

    // Add legend
    TLegend *legend = new TLegend(0.5, 0.7, 0.85, 0.85);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(JPsiEffHist_1, "DQ MC", "l");
    legend->AddEntry(JPsiEffHist_2, "General purpose MC", "l");
    legend->Draw();
    setMax({JPsiEffHist_1, JPsiEffHist_2});

    drawLabel_cuts(MC_name_1, "", &pT_JPsi_cuts, &eta_JPsi_cuts, 0.49, 0.79);
    c1->SaveAs(TString::Format("results/%s/efficiency_plots_overlay.png", MC_name_1.Data()));

}