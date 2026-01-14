
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
}

void inspect_sources() {

    // TString MC_name = "DQ";
    // TString MC_name = "HF";
    TString MC_name = "genpurp";
    TString data_file = "results/" + MC_name + "/muonAOD.root";
    
    // Load the dataframe keys
    TFile *file = TFile::Open(data_file);

    plot_source_matrix(file, MC_name);
}