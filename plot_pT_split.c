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

    // TString MC_name = "DQ";
    // TString MC_name = "HF";
    // TString MC_name = "genpurp";
    // TString MC_name = "k4h_baseline";
    // TString MC_name = "k4h_cuts";
    TString MC_name = "k4h_standalone";
    
    // TString type = "reco";
    TString type = "gen";

    // bool showtype = false;
    bool showtype = true;

    TString data_file = "results/" + MC_name + "/" + type + "/muonAOD.root";
    int n_bins = 20;
    float range_min = 0;
    float range_max = 10;

    // Load the dataframe keys
    TFile *file = TFile::Open(data_file);
    auto *keys = file->GetListOfKeys();
    TTree *tree;

    // Initialize inv mass variables
    // std::vector<Double_t> all_pT, JPsi_pT, Psi2S_pT, charm_pT, b_pT, K_pT, Pi_pT, LM_pT, noMC_pT, other_pT;
    std::vector<Double_t> all_pT, JPsi_pT, Psi2S_pT, charm_pT, b_pT, LM_pT, noMC_pT, other_pT;

    // Loop over dataframes
    for (int i = 0; i < keys->GetEntries()-1; ++i) {

        std::cout << "Processing " << i << "/" << keys->GetEntries()-1 << std::endl;

        tree = get_tree((TKey*)keys->At(i), file);

        // Prepare to read muon info
        Long64_t fMotherPDG;
        float fPt;
        bool fIsProducedInTransport;
        // tree->SetBranchAddress("fMotherPDG", &fMotherPDG);
        tree->SetBranchAddress("fGrandmotherPDG", &fMotherPDG);
        tree->SetBranchAddress("fPtassoc", &fPt);
        tree->SetBranchAddress("fIsProducedInTransport", &fIsProducedInTransport);
        
        // First pass: build groups of muons from the same event
        Long64_t n = tree->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            tree->GetEntry(i);

            all_pT.push_back(fPt);

            if (std::abs(fMotherPDG) == 443) {
                JPsi_pT.push_back(fPt);
            } else if (std::abs(fMotherPDG) == 100443) {
                Psi2S_pT.push_back(fPt);
            } else if ((std::abs(fMotherPDG) >= 411 && std::abs(fMotherPDG) <= 445) || 
                       (std::abs(fMotherPDG) >= 4101 && std::abs(fMotherPDG) <= 4444)) {
                charm_pT.push_back(fPt);
            } else if ((std::abs(fMotherPDG) >= 511 && std::abs(fMotherPDG) <= 557) || 
                       (std::abs(fMotherPDG) >= 5101 && std::abs(fMotherPDG) <= 5554)) {
                b_pT.push_back(fPt);
            } else if (std::abs(fMotherPDG) == 321 || std::abs(fMotherPDG) == 311 || 
                       std::abs(fMotherPDG) == 130 || std::abs(fMotherPDG) == 323 ) {
                // K_pT.push_back(fPt);
                LM_pT.push_back(fPt);
            } else if (std::abs(fMotherPDG) == 211 || std::abs(fMotherPDG) == 113 || 
                       std::abs(fMotherPDG) == 111 || std::abs(fMotherPDG) == 213 ) {
                // Pi_pT.push_back(fPt);
                LM_pT.push_back(fPt);
            } else if (std::abs(fMotherPDG) == 221 || std::abs(fMotherPDG) == 331 || 
                       std::abs(fMotherPDG) == 223 || std::abs(fMotherPDG) == 333 ) {
                LM_pT.push_back(fPt);
            } else if (fMotherPDG == -9999) {
                noMC_pT.push_back(fPt);
            } else {
                other_pT.push_back(fPt);
            }
        }

        tree->ResetBranchAddresses();
        delete tree;
    }

    // Plot histogram of pT distributions
    TCanvas *c1 = new TCanvas("c1", "pT of muons", 800, 600);
    TH1F *pTHist = new TH1F("h1","pT of Muons;pT (GeV/c);Counts",n_bins,range_min,range_max);
    pTHist->FillN(all_pT.size(), all_pT.data(), nullptr);
    pTHist->SetLineWidth(3); 
    pTHist->SetLineColor(kBlack); 
    pTHist->Draw();

    // Add secondary histograms
    TH1F *jpsiHist = new TH1F("h2","",n_bins,range_min,range_max);
    jpsiHist->FillN(JPsi_pT.size(), JPsi_pT.data(), nullptr);
    jpsiHist->SetLineColor(kMagenta);
    jpsiHist->SetLineWidth(2);
    jpsiHist->Sumw2();
    jpsiHist->Draw("HIST SAME");
    jpsiHist->Draw("E SAME");

    TH1F *psi2sHist = new TH1F("h3","",n_bins,range_min,range_max);
    psi2sHist->FillN(Psi2S_pT.size(), Psi2S_pT.data(), nullptr);
    psi2sHist->SetLineColor(kMagenta+1);
    psi2sHist->SetLineWidth(2);
    psi2sHist->Sumw2();
    psi2sHist->Draw("HIST SAME");
    psi2sHist->Draw("E SAME");

    TH1F *charmHist = new TH1F("h4","",n_bins,range_min,range_max);
    charmHist->FillN(charm_pT.size(), charm_pT.data(), nullptr);
    charmHist->SetLineColor(kGreen);
    charmHist->SetLineWidth(2);
    charmHist->Sumw2();
    charmHist->Draw("HIST SAME");
    charmHist->Draw("E SAME");

    TH1F *bHist = new TH1F("h5","",n_bins,range_min,range_max);
    bHist->FillN(b_pT.size(), b_pT.data(), nullptr);
    bHist->SetLineColor(kRed);
    bHist->SetLineWidth(2);
    bHist->Sumw2();
    bHist->Draw("HIST SAME");
    bHist->Draw("E SAME");

    // TH1F *KHist = new TH1F("h6","",n_bins,range_min,range_max);
    // KHist->FillN(K_pT.size(), K_pT.data(), nullptr);
    // KHist->SetLineColor(kBlue);
    // KHist->SetLineWidth(2);
    // KHist->Sumw2();
    // KHist->Draw("HIST SAME");
    // KHist->Draw("E SAME");

    // TH1F *PiHist = new TH1F("h7","",n_bins,range_min,range_max);
    // PiHist->FillN(Pi_pT.size(), Pi_pT.data(), nullptr);
    // PiHist->SetLineColor(kOrange);
    // PiHist->SetLineWidth(2);
    // PiHist->Sumw2();
    // PiHist->Draw("HIST SAME");
    // PiHist->Draw("E SAME");

    TH1F *LMHist = new TH1F("h8","",n_bins,range_min,range_max);
    LMHist->FillN(LM_pT.size(), LM_pT.data(), nullptr);
    // LMHist->SetLineColor(kOrange+2);
    LMHist->SetLineColor(kBlue);
    LMHist->SetLineWidth(2);
    LMHist->Sumw2();
    LMHist->Draw("HIST SAME");
    LMHist->Draw("E SAME");

    TH1F *noMCHist = new TH1F("h9","",n_bins,range_min,range_max);
    noMCHist->FillN(noMC_pT.size(), noMC_pT.data(), nullptr);
    noMCHist->SetLineColor(kGray);
    noMCHist->SetLineWidth(2);
    noMCHist->Sumw2();
    noMCHist->Draw("HIST SAME");
    noMCHist->Draw("E SAME");

    TH1F *otherHist = new TH1F("h10","",n_bins,range_min,range_max);
    otherHist->FillN(other_pT.size(), other_pT.data(), nullptr);
    otherHist->SetLineColor(kCyan);
    otherHist->SetLineWidth(2);
    otherHist->Sumw2();
    otherHist->Draw("HIST SAME");
    otherHist->Draw("E SAME");
    
    // Style
    gPad->SetLogy();
    increaseMargins(c1);
    pTHist->SetMaximum(8 * pTHist->GetMaximum());
    pTHist->SetMinimum(1);

    // Legend
    TLegend *legend = new TLegend(0.57,0.6,0.94,0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(pTHist, "All #mu", "l");
    legend->AddEntry(jpsiHist, "#mu from J/#psi", "l");
    legend->AddEntry(psi2sHist, "#mu from #psi(2S)", "l");
    legend->AddEntry(charmHist, "#mu from Charm", "l");
    legend->AddEntry(bHist, "#mu from Beauty", "l");
    // legend->AddEntry(KHist, "#mu from Kaons", "l");
    // legend->AddEntry(PiHist, "#mu from Pions", "l");
    // legend->AddEntry(LMHist, "#mu from other light mesons (#eta,#omega,#phi)", "l");
    legend->AddEntry(LMHist, "Decay #mu", "l");
    legend->AddEntry(otherHist, "Other #mu", "l");
    legend->AddEntry(noMCHist, "No mother", "l");
    legend->Draw();

    // Print info
    std::cout << "Total muons: " << all_pT.size() << std::endl;
    std::cout << "Muons from J/psi: " << JPsi_pT.size() << std::endl;
    std::cout << "Muons from Psi(2S): " << Psi2S_pT.size() << std::endl;
    std::cout << "Muons from Charm: " << charm_pT.size() << std::endl;
    std ::cout << "Muons from Beauty: " << b_pT.size() << std::endl;
    std::cout << "Muons from Light Mesons: " << LM_pT.size() << std::endl;
    std::cout << "Muons with no MC info: " << noMC_pT.size() << std::endl;
    std::cout << "Other muons: " << other_pT.size() << std::endl;

    if (showtype) {
        drawLabel(MC_name, type, 0.55, 0.89);
    } else {
        drawLabel(MC_name, "", 0.55, 0.89);
    }

    TString out_name = TString::Format("results/%s/%s/muon_pT_split", MC_name.Data(), type.Data());
    out_name.ReplaceAll(".", "_");
    c1->SaveAs(out_name + ".png");
}