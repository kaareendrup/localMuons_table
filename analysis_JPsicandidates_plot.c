
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

void analysis_JPsicandidates_plot() {
    // This function loads the histograms created in 
    // analysis_JPsicandidates, applies scaling if needed, 
    // and creates plots for the invariant mass spectra and pT distributions.

    ////////////////////////////////////////////////////////////////////
    ////            Load configuration, setup up filenames          ////
    ////////////////////////////////////////////////////////////////////
    std::cout << std::fixed << std::setprecision(1);
    SetALICEStyle();

    std::ifstream jsonFile("localMuons_table/config/config_analysis.json");
    json config;
    jsonFile >> config;

    std::string data_name = config["data_name"];
    std::string muon_type = config["muon_type"];
    float JPsi_branching_ratio = config["branching_ratio"];

    TString MC_name = TString::Format("%s_%s", data_name.c_str(), muon_type.c_str());

    ////////////////////////////////////////////////////////////////////
    ////                Load histograms from file                   ////
    ////////////////////////////////////////////////////////////////////

    // Load and plot invariant mass spectra by pT bins
    createInvMassHist("reco", config, MC_name);
    createInvMassHist("gen", config, MC_name);
    
    // Load and plot pT spectra, extracting signal-background
    TH1F* pT_reco = createPTHist("reco", config, MC_name);
    TH1F* pT_gen = createPTHist("gen", config, MC_name);

    // Import metadata
    int nEvents = getNEvents(MC_name);

    // Get efficiency for scaling
    std::vector<double> efficiency = get_efficiency(config, MC_name);

    ////////////////////////////////////////////////////////////////////
    ////                Plot corrected J/Psi spectra                ////
    ////////////////////////////////////////////////////////////////////
    
    // Create uncorrected and efficiency-corrected pT spectra, and plot them together with the generated distribution
    TCanvas *c3 = new TCanvas("c3", "pT bin counts", 900, 400);
    c3->Divide(2,1);

    // Plot uncorrected distribution
    c3->cd(1);
    pT_reco->Draw();
    TLegend *leg3s = new TLegend(0.4,0.6,0.9,0.9);
    leg3s->AddEntry(pT_reco, "Reconstructed (not corrected)", "l");
    leg3s->Draw();

    // Plot generated distribution
    c3->cd(2);
    pT_gen->SetLineColor(kBlue);
    pT_gen->Draw();

    // Plot efficiency-corrected distribution
    TH1F *pT_reco_scale = (TH1F*)pT_reco->Clone("pTscale");
    scale_histogram(pT_reco_scale, efficiency, nEvents);
    pT_reco_scale->Scale(JPsi_branching_ratio);
    pT_reco_scale->Draw("same");
    
    setMax({pT_gen, pT_reco_scale});
    TLegend *leg3 = new TLegend(0.4,0.7,0.9,0.9);
    leg3->AddEntry(pT_gen, "Generated", "l");
    leg3->AddEntry(pT_reco_scale, "Reconstructed\n (corrected)", "l");
    leg3->Draw();
    c3->SaveAs(TString::Format("results/%s/pTspectrascaled.png", MC_name.Data()));

    // Plot only the corrected distribution with ratio to generated
    TCanvas *c4 = new TCanvas("c4", "pT bin counts", 600, 600);
    c4->Divide(1,2);
  
    c4->cd(1);
    gPad->SetPad(0,0.3,1,1);    // Top 70%
    pT_gen->Draw();
    pT_reco_scale->Draw("same");
    gPad->SetBottomMargin(0); // Remove bottom margin for top pad
    
    TLegend *leg4 = new TLegend(0.4,0.7,0.9,0.9);
    leg4->AddEntry(pT_gen, "Generated", "l");
    leg4->AddEntry(pT_reco_scale, "Reconstructed\n (corrected)", "l");
    leg4->Draw();

    // Plot ratio
    c4->cd(2);
    gPad->SetPad(0,0,1,0.3);    // Bottom 30%
    gPad->SetTopMargin(0); // Remove top margin for bottom pad
    TH1F *ratio_hist = (TH1F*)pT_gen->Clone("ratio_hist");
    ratio_hist->Divide(pT_reco_scale);
    ratio_hist->Draw();
    gPad->SetGridy();
    ratio_hist->SetMinimum(.95);
    ratio_hist->SetMaximum(1.08);
    ratio_hist->GetYaxis()->SetTitle("Ratio");
    ratio_hist->GetYaxis()->SetTitleOffset(0.5);
    ratio_hist->GetYaxis()->SetTitleSize(0.09);
    ratio_hist->GetYaxis()->SetLabelSize(0.08);
    ratio_hist->GetYaxis()->SetNdivisions(505);
    ratio_hist->GetXaxis()->SetTitleSize(0.09);
    ratio_hist->GetXaxis()->SetLabelSize(0.08);

    gPad->SetBottomMargin(0.3);
    c4->SetBottomMargin(0.8);
    c4->SaveAs(TString::Format("results/%s/pTspectraratio.png", MC_name.Data()));
}