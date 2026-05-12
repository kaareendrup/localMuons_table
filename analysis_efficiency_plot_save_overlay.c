
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
    c1->Divide(1,2);

    c1->cd(1);
    gPad->SetPad(0,0.3,1,1);    // Top 70%
    // JPsiEffHist_1->GetXaxis()->SetRangeUser(0, 10);
    
    drawHist(JPsiEffHist_1, "J/#Psi Efficiency;p_{T} (GeV/c);Eff_{J/#Psi}", kBlue, 1.2);
    drawHist(JPsiEffHist_2, "J/#Psi Efficiency;p_{T} (GeV/c);Eff_{J/#Psi}", kRed, 1.2, true);

    // Add legend
    TLegend *legend = new TLegend(0.42, 0.7, 0.75, 0.85);
    legend->AddEntry(JPsiEffHist_1, "DQ MC", "l");
    legend->AddEntry(JPsiEffHist_2, "General purpose MC", "l");
    legend->Draw();
    setMax({JPsiEffHist_1, JPsiEffHist_2});

    drawLabel_cuts(MC_name_1, "", &config, 0.4, 0.79);
    gPad->SetBottomMargin(0); // Remove bottom margin for top pad

    c1->cd(2);
    gStyle->SetOptTitle(0);
    TH1F *ratio_hist = createRatioPlot(JPsiEffHist_1, JPsiEffHist_2);
    ratio_hist->SetMinimum(0);
    ratio_hist->SetMaximum(2);

    c1->SaveAs(TString::Format("results/%s/efficiency_plots_overlay.png", MC_name_1.Data()));
}