
#include "setALICEStyle.c"

SetALICEStyle();

TTree* get_tree(TKey *key, TFile *file) {

    const char* dirName = key->GetName();
    TDirectoryFile *dir = (TDirectoryFile*)file->Get(dirName);
    TTree *tree = (TTree*)dir->Get("O2dqmuontable");

    return tree;
}

void plot_source_matrix(TFile *file, TString MC_name) {

    auto *keys = file->GetListOfKeys();
    const TString sources[6] = {"fIsPrimary", "fIsProducedInTransport", "fIsProducedByGenerator", "fIsFromBackgroundEvent", "fIsHEPMCFinalState", "fIsPowhegDY"};
    bool flags[6];

    TH2I *sourceMatrix = new TH2I("sourceMatrix", "Source Matrix; ; ", 6, 0, 6, 6, 0, 6);

    // Loop over dataframes
    for (int i = 0; i < keys->GetEntries()-1; ++i) {
        TTree *tree = get_tree((TKey*)keys->At(i), file);

        // Set branch addresses for source flags
        for (int i =0; i < 6; ++i){
            tree->SetBranchAddress(sources[i], &flags[i]);
        }

        // Loop over entries
        Long64_t n = tree->GetEntries();
        for (Long64_t j = 0; j < n; ++j) {
            tree->GetEntry(j);

            // Fill source matrix
            for (int m = 0; m < 6; ++m){
                if (!flags[m]) continue;

                for (int n = m; n < 6; ++n){
                    if (flags[m] && flags[n]){
                        sourceMatrix->Fill(m, n);
                    }
                }
            }
        }
    }

    for (int i = 0; i < 6; ++i){
        sourceMatrix->GetXaxis()->SetBinLabel(i+1, sources[i]);
        sourceMatrix->GetYaxis()->SetBinLabel(i+1, sources[i]);
    }

    TCanvas *c1 = new TCanvas("c1", "Source Matrix", 900, 600);
    gStyle->SetOptStat(0);
    sourceMatrix->Draw("COLZ TEXT");

    c1->SetLeftMargin(0.25);
    c1->SetBottomMargin(0.15);
    c1->SetRightMargin(0.15);
    c1->SetTopMargin(0.15);

    c1->SaveAs("results/" + MC_name + "/source_matrix.png");
    delete c1;
}

void plot_source_pT(TFile *file, TString MC_name) {

    auto *keys = file->GetListOfKeys();
    const TString sources[2] = {"fIsPrimary", "fIsProducedInTransport"};
    bool flags[2];
    std::vector<Double_t> pT_Primary, pT_Transport, pT_Other;
    float pT;
    int n_entries = 0;

    // Loop over dataframes
    for (int i = 0; i < keys->GetEntries()-1; ++i) {
        printf("Processing dataframe %d / %d\n", i+1, keys->GetEntries()-1);
        TTree *tree = get_tree((TKey*)keys->At(i), file);

        // Set branch addresses for source flags and pT
        tree->SetBranchAddress("fPtassoc", &pT);
        tree->SetBranchAddress("fIsPrimary", &flags[0]);
        tree->SetBranchAddress("fIsProducedInTransport", &flags[1]);

        // Loop over entries
        for (Long64_t j = 0; j < tree->GetEntries(); ++j) {
            tree->GetEntry(j); 
            n_entries++;

            if (flags[0]) {
                pT_Primary.push_back(pT);
            }
            else if (flags[1]) {
                pT_Transport.push_back(pT);
            }
            else {
                pT_Other.push_back(pT);
            }
        }
    }

    printf("Total entries processed: %d\n", n_entries);

    TCanvas *c1 = new TCanvas("c1", "Source pT Distributions", 1200, 600);
    TH1F *pTPrimaryHist = new TH1F("h1","Primary Muon pT Distribution;p_{T} (GeV/c);Counts",100,0,10);
    pTPrimaryHist->FillN(pT_Primary.size(), pT_Primary.data(), nullptr);
    pTPrimaryHist->SetLineColor(kBlue+1);
    pTPrimaryHist->SetLineWidth(2);
    pTPrimaryHist->Draw();

    TH1F *pTTransportHist = new TH1F("h2","Transport Muon pT Distribution;p_{T} (GeV/c);Counts",100,0,10);
    pTTransportHist->FillN(pT_Transport.size(), pT_Transport.data(), nullptr);
    pTTransportHist->SetLineColor(kRed+1);
    pTTransportHist->SetLineWidth(2);
    pTTransportHist->Draw("SAME");

    TH1F *pTOtherHist = new TH1F("h3","Other Muon pT Distribution;p_{T} (GeV/c);Counts",100,0,10);
    pTOtherHist->FillN(pT_Other.size(), pT_Other.data(), nullptr);
    pTOtherHist->SetLineColor(kBlack);
    pTOtherHist->SetLineWidth(2);
    pTOtherHist->Draw("SAME");

    increaseMargins(c1);
    drawLabel(MC_name);
    TLegend *legend = new TLegend(0.6,0.7,0.88,0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(pTPrimaryHist, "Primary Muons", "l");
    legend->AddEntry(pTTransportHist, "Transport Muons", "l");
    legend->AddEntry(pTOtherHist, "Other Muons", "l");
    legend->Draw();

    // Adjust ymax 
    double maxY = std::max(pTPrimaryHist->GetMaximum(), pTTransportHist->GetMaximum());
    pTPrimaryHist->SetMaximum(1.2 * maxY);
    pTTransportHist->SetMaximum(1.2 * maxY);

    // Set log Y axis
    gPad->SetLogy();

    c1->SaveAs("results/" + MC_name + "/source_pT_distributions.png");
    delete c1;
}

void inspect_sources() {

    // TString MC_name = "DQ";
    // TString MC_name = "HF";
    TString MC_name = "genpurp";
    TString data_file = "results/" + MC_name + "/muonAOD.root";
    
    // Load the dataframe keys
    TFile *file = TFile::Open(data_file);

    plot_source_matrix(file, MC_name);
    plot_source_pT(file, MC_name);
}