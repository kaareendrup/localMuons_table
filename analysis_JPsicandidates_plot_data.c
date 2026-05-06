
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
#include "utils/efficiency.c"


void analysis_JPsicandidates_plot_data() {
    // This function loads the histograms created in 
    // analysis_JPsicandidates, applies scaling if needed, 
    // and creates plots for the invariant mass spectra and pT distributions.

    ////////////////////////////////////////////////////////////////////
    ////            Load configuration, setup up filenames          ////
    ////////////////////////////////////////////////////////////////////
    std::cout << std::fixed << std::setprecision(1);
    SetALICEStyle();

    // Load config
    std::ifstream jsonFile("localMuons_table/config/config_analysis.json");
    json config;
    jsonFile >> config;

    std::string data_name = config["data_name"];
    std::string muon_type = config["muon_type"];
    std::string eff_source = config["data_eff_source"];
    std::vector<double> pT_bins = config["hists"]["pT_bins"].get<std::vector<double>>();
    float cuts_pT_JPsi_min = config["cuts_JPsi"]["pT_JPsi_min"];
    float cuts_pT_JPsi_max = config["cuts_JPsi"]["pT_JPsi_max"];

    TString data = TString::Format("%s_%s", data_name.c_str(), muon_type.c_str());
    TString MC_name = TString::Format("%s_%s", eff_source.c_str(), muon_type.c_str());
    
    ////////////////////////////////////////////////////////////////////
    ////                Load metadata and hepdata                   ////
    ////////////////////////////////////////////////////////////////////
    int nEvents = getNEvents(data);
    
    // Import efficiency
    std::vector<double> efficiency = get_efficiency(config, MC_name);

    // Import hepdata points
    std::string hepdata_name = config["hepdata_name"];
    double crossSection = config["hepdata_crosssection"];

    TString hepdata_file = TString::Format("results/%s.json", hepdata_name.c_str());
    std::ifstream hepdata_json(hepdata_file);
    json hepdata;
    hepdata_json >> hepdata;

    std::vector<double> hepdata_values;
    std::vector<double> hepdata_errors;
    for (const auto& point : hepdata["values"]) {
        double pt_low  = std::stod(point["x"][0]["low"].get<std::string>());
        double pt_high = std::stod(point["x"][0]["high"].get<std::string>());

        const auto& y = point["y"][0];
        double value = std::stod(y["value"].get<std::string>());

        double err_stat = 0.0;
        double err_sys  = 0.0;

        for (const auto& err : y["errors"]) {
            std::string label = err["label"];
            double e = std::stod(err["symerror"].get<std::string>());

            if (label == "stat") err_stat = e;
            if (label == "sys")  err_sys  = e;
        }

        if (pt_low >= cuts_pT_JPsi_min && pt_high <= cuts_pT_JPsi_max) { // Only include points within the plotted range
            hepdata_values.push_back(value);
            hepdata_errors.push_back(sqrt(pow(err_stat, 2) + pow(err_sys, 2)));
        }
    }
    
    ////////////////////////////////////////////////////////////////////
    ////                Load histograms from file                   ////
    ////////////////////////////////////////////////////////////////////
    createInvMassHist("reco", config, data);
    TH1F* pT_reco = createPTHist("reco", config, data);
    
    ////////////////////////////////////////////////////////////////////
    ////        Plot corrected and uncorrected J/Psi spectra        ////
    ////////////////////////////////////////////////////////////////////
    TCanvas *c3 = new TCanvas("c3", "pT bin counts", 900, 400);
    c3->Divide(2,1);
    
    // Create uncorrected histogram
    c3->cd(1);
    pT_reco->Draw();
    TLegend *leg3s = new TLegend(0.4,0.6,0.9,0.9);
    leg3s->AddEntry(pT_reco, "Reconstructed (not corrected)", "l");
    leg3s->Draw();

    // Create scaled histogram
    c3->cd(2);

    TH1F *pT_reco_scale = (TH1F*)pT_reco->Clone("pTscale");
    scale_histogram(pT_reco_scale, efficiency, nEvents);
    pT_reco_scale->Draw("same");
    
    pT_reco_scale->GetYaxis()->SetTitle("#frac{1}{N_{events}} d^{2}N/(dp_{T} dy) (GeV/c)^{-1}");
    setMax({pT_reco_scale});
    TLegend *leg3 = new TLegend(0.4,0.7,0.9,0.9);
    leg3->AddEntry(pT_reco_scale, "Reconstructed\n (corrected)", "l");
    leg3->Draw();
    c3->SaveAs(TString::Format("results/%s/pTspectrascaled.png", data.Data()));

    // Create a canvas with two pads: top for the histogram, bottom for the ratio
    TCanvas *c4 = new TCanvas("c4", "pT bin counts", 700, 600);
    c4->Divide(1,2);
  
    c4->cd(1);
    gPad->SetPad(0,0.3,1,1);    // Top 70%
    pT_reco_scale->Draw();
    
    // Add hepdata points with error bars
    TH1F *pTHepData = new TH1F("pTHepData", "pTHepData", pT_bins.size()-1, pT_bins.data());
    for (size_t i = 0; i < hepdata_values.size(); ++i) {
        pTHepData->SetBinContent(i+1, hepdata_values[i]);
        pTHepData->SetBinError(i+1, hepdata_errors[i]);
    }
    pTHepData->Scale(crossSection*1e-6); // Scale by cross section if needed
    // pTHepData->SetMarkerStyle(20);
    pTHepData->SetMarkerColor(kBlack);
    pTHepData->SetLineColor(kBlack);
    pTHepData->Draw("E1 SAME");
    
    pT_reco_scale->SetMinimum(1e-5);
    
    TLegend *leg4 = new TLegend(0.1,0.1,0.6,0.25);
    leg4->AddEntry(pT_reco_scale, "Reconstructed\n (corrected)", "l");
    leg4->AddEntry(pTHepData, "ALICE 2017", "lep");
    leg4->Draw();

    setMax({pT_reco_scale, pTHepData});
    gPad->SetLogy();
    increaseMargins(c4);
    drawLabel_cuts(data, "", &config, 0.85, 0.8, false, true, 0.05);
    gPad->SetBottomMargin(0); // Remove bottom margin for top pad

    c4->cd(2);
    TH1F *ratio_hist = createRatioPlot(pT_reco_scale, pTHepData);
    ratio_hist->SetMinimum(.3);
    ratio_hist->SetMaximum(1.7);

    c4->SaveAs(TString::Format("results/%s/pTspectracompare.png", data.Data()));
}