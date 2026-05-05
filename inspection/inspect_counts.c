
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TString.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

void SetALICEStyle() {
    gStyle->SetOptStat(0);          // No stat box
    gStyle->SetOptTitle(0);         // No default title

    gStyle->SetPadTickX(1);         // Ticks on both sides
    gStyle->SetPadTickY(1);

    gStyle->SetFrameLineWidth(2);
    gStyle->SetLineWidth(2);

    gStyle->SetLabelSize(0.05, "XY");
    gStyle->SetTitleSize(0.06, "XY");
    gStyle->SetTitleOffset(1.2, "X");
    gStyle->SetTitleOffset(1.3, "Y");

    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFont(42);

    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "XY");
    gStyle->SetTitleFont(42, "XY");

    gStyle->SetLegendBorderSize(0);     // No border
    gStyle->SetLegendFillStyle(0);     // Transparent background
}

template <typename T>
void process_file(TString reco_file, TString treename, TString event_branch, int &n_muons_reco, int &n_events, std::map<int, int>& muons_per_event, T value_type = 0) {

    TFile *recoFile = TFile::Open(reco_file);
    TIter nextRecoKey(recoFile->GetListOfKeys());
    TKey* recoKey;
    
    int dirCount = 0;

    // Loop over dataframes
    while ((recoKey = (TKey*) nextRecoKey())) {

        // Load directory and tree
        TObject* obj = recoKey->ReadObj();
        if (!(obj->InheritsFrom("TDirectory"))) continue;

        TDirectory* recoDir = (TDirectory*) obj;
        TTree *muonRecoTree = (TTree*)recoDir->Get(treename);

        std::cout << TString::Format("Reading tracks from dir %d of %d: %s\r", dirCount, recoFile->GetListOfKeys()->GetEntries(), recoDir->GetName()) << std::flush;

        std::map<T, int> event_groups;
        T fEventIdx;
        muonRecoTree->SetBranchAddress(event_branch, &fEventIdx);

        // First pass: build groups of muons from the same event
        Long64_t n = muonRecoTree->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            muonRecoTree->GetEntry(i);
            event_groups[fEventIdx]++;
            n_muons_reco++;
        }

        dirCount++;
        for (const auto& group : event_groups) {
            muons_per_event[group.second]++;
            n_events++;
        }
    }

    delete recoKey;
    recoFile->Close();
    delete recoFile;
}

TH1F* make_hist(std::map<int, int>& muons_per_event, TString name) {

    TH1F* hist = new TH1F(name, "Muons per Event;Number of Muons;Number of Events", 15, -0.5, 14.5);
    for (const auto& pair : muons_per_event) {
        hist->Fill(pair.first, pair.second);
    }
    return hist;
}

void print_stats(int &n_muons_reco, int &n_events, std::map<int, int>& muons_per_event) {

    std::cout << std::endl;
    std::cout << "Total muons reco: " << n_muons_reco << std::endl;
    std::cout << "Total events: " << n_events << std::endl;
    std::cout << "Muons per event:" << std::endl;
    int tot_threeormore = 0;
    for (const auto& pair : muons_per_event) {
        std::cout << "  " << pair.first << ": " << pair.second << std::endl;
        if (pair.first >= 3) tot_threeormore += pair.second;
    }
    std::cout << "Total events with 3 or more muons: " << tot_threeormore << std::endl;
}

TH1F* inspect_input(TString data_name, int n_files=0) {

    // Load input data file list
    TString input_file_list;
    if (n_files > 0) {
        input_file_list = TString::Format("input_data/%s/input_data_%d.txt", data_name.Data(), n_files);
    } else {
        input_file_list = TString::Format("input_data/%s/input_data.txt", data_name.Data());
    }

    std::ifstream infile(input_file_list.Data());
    std::string line;

    int n_muons_reco = 0;
    int n_events = 0;
    std::map<int, int> muons_per_event;
    int file_count = 0;

    // Loop over input files
    while (std::getline(infile, line)) {

        std::cout << TString::Format("Processing file %d of %d", file_count, n_files) << std::flush;

        process_file(line.c_str(), "O2reducedmuon", "fIndexReducedEvents", n_muons_reco, n_events, muons_per_event, 0);

        file_count++;
    }

    print_stats(n_muons_reco, n_events, muons_per_event);
    return make_hist(muons_per_event, "inputhist");
}

TH1F* inspect_analysis(TString data_name, int n_files, std::vector<int> exceptions = {}) {
    
    int n_muons_reco = 0;
    int n_events = 0;
    std::map<int, int> muons_per_event;

    for (int i = 0; i < n_files; ++i) {
        if (std::find(std::begin(exceptions), std::end(exceptions), i) != std::end(exceptions)) {
            continue;
        }
        
        std::cout << "Processing file " << i << " of " << n_files << std::flush;

        TString reco_file = TString::Format("results/%s/reco/muonAOD%d.root", data_name.Data(), i);
        process_file(reco_file, "O2dqmuontable", "fEventIdx", n_muons_reco, n_events, muons_per_event, 0uLL);

    }

    print_stats(n_muons_reco, n_events, muons_per_event);
    return make_hist(muons_per_event, "analysishist");
}

void inspect_triggers(TString data_name, bool cuts = true) {

    int trigger_count = 0;
    int assoc_count = 0;

    TString cut_suffix = cuts ? "" : "_nocuts";

    TString data_file = TString::Format("results/%s/reco/eventmuons%s.root", data_name.Data(), cut_suffix.Data());
    TFile *file = TFile::Open(data_file);

    TTree *tree = nullptr;
    file->GetObject("Triggers", tree);

    std::vector<double> *eta_assocs = nullptr;
    tree->SetBranchAddress("eta_assocs",  &eta_assocs);

    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {

        tree->GetEntry(i);
        std::cout << "Processing entry " << i+1 << " of " << tree->GetEntries() << "\r" << std::flush;

        trigger_count++;
        assoc_count += eta_assocs->size();
    }
    std::cout << std::endl;

    std::cout << "Cuts: " << (cuts ? "ON" : "OFF") << std::endl;
    std::cout << "Total triggers: " << trigger_count << std::endl;
    std::cout << "Total associated particles: " << assoc_count << std::endl;
}

void inspect_counts() {

    SetALICEStyle();

    TString data_name = "DQ_data_standalone";
    TH1F* hist_input = inspect_input(data_name, 0);
    TH1F* hist_analysis = inspect_analysis(data_name, 172);
    inspect_triggers(data_name, false);
    inspect_triggers(data_name);

    // Draw hists
    TCanvas *c1 = new TCanvas("c1", "Muons per Event", 800, 600);
    hist_input->SetLineColor(kBlue);
    hist_input->SetMarkerColor(kBlue);
    hist_input->Draw("HIST");
    hist_analysis->SetLineColor(kRed);
    hist_analysis->SetMarkerColor(kRed);
    hist_analysis->Draw("HIST SAME");
    gPad->SetLogy();
    // c1->BuildLegend();
    c1->SaveAs(TString::Format("results/%s/muons_per_event.png", data_name.Data()));
}