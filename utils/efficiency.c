
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1F.h>

#include <vector>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

std::vector<double> get_efficiency(json config, TString MC_name, bool muons = false) {

    // Data
    TString in_file = TString::Format("results/%s/efficiency.root", MC_name.Data());
    TFile* file = TFile::Open(in_file, "READ");

    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open efficiency file: " << in_file << "\n";
        return {};
    }

    // Efficiency plots
    TH1F *EffHist = nullptr;
    if (!muons) {
        EffHist = (TH1F*)file->Get("JPsiEffHist");
    } else {
        EffHist = (TH1F*)file->Get("muonEffHist");
    }

    int nbins = config["hists"]["pT_bins"].size();
    std::vector<double> efficiency(nbins);
    for (int i = 0; i < nbins; ++i) {
        efficiency[i] = EffHist->GetBinContent(i+1);
    }
    return efficiency;
}
