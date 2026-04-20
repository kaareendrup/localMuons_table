
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1F.h>

#include <vector>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

std::vector<double> get_efficiency(json config, TString MC_name) {

    // Data
    TString in_file = TString::Format("results/%s/efficiency.root", MC_name.Data());
    TFile* file = TFile::Open(in_file, "READ");

    // Efficiency plots
    TH1F *JPsiEffHist = (TH1F*)file->Get("JPsiEffHist");

    int nbins = config["hists"]["pT_bins"].size();
    std::vector<double> efficiency(nbins);
    for (int i = 0; i < nbins; ++i) {
        efficiency[i] = JPsiEffHist->GetBinContent(i+1);
    }
    return efficiency;
}
