
double getRapidity(double pT, double eta) {
    double mJPsi = 3.096916; // J/Psi mass in GeV/c^2
    return log((sqrt(pow(mJPsi, 2) + (pow(pT, 2) * pow(cosh(eta), 2))) + pT * sinh(eta)) / (sqrt(pow(mJPsi, 2) + pow(pT, 2))));
}

double getDeltaY(double pT, double eta_min, double eta_max) {

    double y_min = getRapidity(pT, eta_min);
    double y_max = getRapidity(pT, eta_max);
    // std::cout << "pT = " << pT << " GeV/c, y_min = " << y_min << ", y_max = " << y_max << std::endl;
    return y_max - y_min;
}

void fillHist(int idx, const TString& name, std::map<int, std::unique_ptr<TH1F>>& hists, double value, int n_bins, double x_min, double x_max) {

    // If histogram doesn't exist yet, create it
    if (hists.find(idx) == hists.end()) {
        hists[idx] = std::make_unique<TH1F>(name, name, n_bins, x_min, x_max);
        hists[idx]->SetDirectory(nullptr); 
    }

    hists[idx]->Fill(value);
}
