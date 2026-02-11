#include "setALICEStyle.c"

std::cout << std::fixed << std::setprecision(1);
SetALICEStyle();

void plot_pT_data(){

    TString MC_name = "DQ_data";    
    TString type = "reco";
    
    TString data_file = "results/" + MC_name + "/" + type + "/muonAOD.root";
    int n_bins = 20;
    float range_min = 0;
    float range_max = 10;

    // Load the dataframe keys
    TFile *file = TFile::Open(data_file);
    TIter nextKey(file->GetListOfKeys());
    TKey* key;

    // Initialize inv mass variables
    std::vector<Double_t> pT;
    int dirCount = 0;

    // Loop over dataframes
    while ((key = (TKey*) nextKey())) {
        
        TObject*obj = key->ReadObj();
        if (!(obj->InheritsFrom("TDirectory"))) continue;
        TDirectory* dir = (TDirectory*) obj;
        TTree *tree = (TTree*)dir->Get("O2dqmuontable");

        printf("Reading tracks from dir %d of %d: %s\n", dirCount, file->GetListOfKeys()->GetEntries(), dir->GetName());

        // Prepare to read muon info
        float fPt;
        tree->SetBranchAddress("fPtassoc", &fPt);
        
        // First pass: build groups of muons from the same event
        Long64_t n = tree->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            tree->GetEntry(i);
            pT.push_back(fPt);
        }

        tree->ResetBranchAddresses();
        delete tree;

        dirCount++;
    }

    // Plot histogram of pT distributions
    TCanvas *c1 = new TCanvas("c1", "pT of muons", 800, 600);
    TH1F *pTHist = new TH1F("h1","pT of Muons;pT (GeV/c);Counts",n_bins,range_min,range_max);
    pTHist->FillN(pT.size(), pT.data(), nullptr);
    pTHist->Sumw2();
    pTHist->SetLineWidth(2); 
    pTHist->SetLineColor(kBlue); 
    pTHist->Draw();
    
    // Style
    gPad->SetLogy();
    increaseMargins(c1);
    pTHist->SetMaximum(1000 * pTHist->GetMaximum());
    pTHist->SetMinimum(1);

    // Legend
    TLegend *legend = new TLegend(0.6,0.8,0.7,0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(pTHist, "All #mu", "l");
    legend->Draw();

    // Print info
    std::cout << "Total muons: " << pT.size() << std::endl;
    drawLabel(MC_name, "", 0.55, 0.89);

    TString out_name = TString::Format("results/%s/%s/muon_pT", MC_name.Data(), type.Data());
    out_name.ReplaceAll(".", "_");
    c1->SaveAs(out_name + ".png");
}