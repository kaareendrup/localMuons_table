#include "setALICEStyle.c"

std::cout << std::fixed << std::setprecision(1);
SetALICEStyle();

TTree* get_tree(TKey *key, TFile *file) {

    const char* dirName = key->GetName();
    TDirectoryFile *dir = (TDirectoryFile*)file->Get(dirName);
    TTree *tree = (TTree*)dir->Get("O2dqmuontable");

    return tree;
}

void plot_pT_split(){

    TString MC_name = "DQ";
    // TString MC_name = "HF";
    // TString MC_name = "genpurp";
    TString data_file = "results/" + MC_name + "/muonAOD.root";
    
    float range_min = 0;
    float range_max = 10;

    // Load the dataframe keys
    TFile *file = TFile::Open(data_file);
    auto *keys = file->GetListOfKeys();

    // Initialize inv mass variables
    std::vector<Double_t> all_pT, JPsi_pT, Psi2S_pT, charm_pT, b_pT, decay_pT, secondary_pT, other_pT;

    // Loop over dataframes
    for (int i = 0; i < keys->GetEntries()-1; ++i) {
        TTree *tree = get_tree((TKey*)keys->At(i), file);
                
        // Prepare to read muon info
        Long64_t fMotherPDG;
        float fPt;
        bool fIsSecondary;
        tree->SetBranchAddress("fMotherPDG", &fMotherPDG);
        tree->SetBranchAddress("fPtassoc", &fPt);
        tree->SetBranchAddress("fIsSecondaryMuon", &fIsSecondary);
        
        // First pass: build groups of muons from the same event
        Long64_t n = tree->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            tree->GetEntry(i);

            all_pT.push_back(fPt);

            if (std::abs(fMotherPDG) == 443) {
                JPsi_pT.push_back(fPt);
            }
            else if (std::abs(fMotherPDG) == 100443) {
                Psi2S_pT.push_back(fPt);
            }
            else if ((std::abs(fMotherPDG) >= 411 && std::abs(fMotherPDG) <= 441) || 
                     (std::abs(fMotherPDG) >= 4101 && std::abs(fMotherPDG) <= 4444)) {
                charm_pT.push_back(fPt);
            }
            else if ((std::abs(fMotherPDG) >= 511 && std::abs(fMotherPDG) <= 557) || 
                     (std::abs(fMotherPDG) >= 5101 && std::abs(fMotherPDG) <= 5554)) {
                b_pT.push_back(fPt);
            }
            else if (std::abs(fMotherPDG) == 211 || std::abs(fMotherPDG) == 321) {
                decay_pT.push_back(fPt);
            }else if (fIsSecondary) {
                secondary_pT.push_back(fPt);
            } 
            else {
                other_pT.push_back(fPt);
            }
        }
    }

    // Plot histogram of candidate invariant masses
    TCanvas *c1 = new TCanvas("c1", "pT of muons", 800, 600);
    TH1F *pTHist = new TH1F("h1","pT of Muons;pT (GeV/c);Counts",100,range_min,range_max);
    pTHist->FillN(all_pT.size(), all_pT.data(), nullptr);
    pTHist->SetLineWidth(3); 
    pTHist->SetLineColor(kBlack); 
    pTHist->Draw();
    // Add secondary histograms
    TH1F *jpsiHist = new TH1F("h2","",100,range_min,range_max);
    jpsiHist->FillN(JPsi_pT.size(), JPsi_pT.data(), nullptr);
    jpsiHist->SetLineColor(kGreen+1);
    jpsiHist->SetLineWidth(2);
    jpsiHist->Draw("SAME");

    TH1F *psi2sHist = new TH1F("h3","",100,range_min,range_max);
    psi2sHist->FillN(Psi2S_pT.size(), Psi2S_pT.data(), nullptr);
    psi2sHist->SetLineColor(kGreen+2);
    psi2sHist->SetLineWidth(2);
    psi2sHist->Draw("SAME");

    TH1F *charmHist = new TH1F("h4","",100,range_min,range_max);
    charmHist->FillN(charm_pT.size(), charm_pT.data(), nullptr);
    charmHist->SetLineColor(kMagenta);
    charmHist->SetLineWidth(2);
    charmHist->Draw("SAME");

    TH1F *bHist = new TH1F("h5","",100,range_min,range_max);
    bHist->FillN(b_pT.size(), b_pT.data(), nullptr);
    bHist->SetLineColor(kRed);
    bHist->SetLineWidth(2);
    bHist->Draw("SAME");

    TH1F *decayHist = new TH1F("h6","",100,range_min,range_max);
    decayHist->FillN(decay_pT.size(), decay_pT.data(), nullptr);
    decayHist->SetLineColor(kBlue);
    decayHist->SetLineWidth(2);
    decayHist->Draw("SAME");

    TH1F *secondaryHist = new TH1F("h7","",100,range_min,range_max);
    secondaryHist->FillN(secondary_pT.size(), secondary_pT.data(), nullptr);
    secondaryHist->SetLineColor(kYellow);
    secondaryHist->SetLineWidth(2);
    secondaryHist->Draw("SAME");

    TH1F *otherHist = new TH1F("h8","",100,range_min,range_max);
    otherHist->FillN(other_pT.size(), other_pT.data(), nullptr);
    otherHist->SetLineColor(kCyan);
    otherHist->SetLineWidth(2);
    otherHist->Draw("SAME");
    
    // Style
    gPad->SetLogy();
    increaseMargins(c1);

    // Legend
    TLegend *legend = new TLegend(0.6,0.6,0.88,0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(pTHist, "All", "l");
    legend->AddEntry(jpsiHist, "J/#psi", "l");
    legend->AddEntry(psi2sHist, "#psi(2S)", "l");
    legend->AddEntry(charmHist, "Charm", "l");
    legend->AddEntry(bHist, "Beauty", "l");
    legend->AddEntry(decayHist, "Decay muons", "l");
    legend->AddEntry(secondaryHist, "Secondary muons", "l");
    legend->AddEntry(otherHist, "Other", "l");
    legend->Draw();

    drawLabel(MC_name, 0.18, 0.89);
    TString out_name = TString::Format("results/%s/muon_pT_split", MC_name.Data());
    out_name.ReplaceAll(".", "_");
    c1->SaveAs(out_name + ".png");
}