void prac_pi0_kinematics(TString file, Int_t numBins) {

    // draw some of the nice kinematics functions: Q2, t1, the missing mass that shows the pion missing mass Mx2_1, phi2
    // from what I understand of the source code Mx2_i where i is the non-electron final state particle that is INCLUDED ie not missing so the proton is included for Mx2_1 ie the pion is missing

    TFile *f = TFile::Open(file);
    TTree *t = (TTree*)f->Get("PhysicsEvents");

    TH1F *hQ2 = new TH1F("hQ2", "hQ2", numBins, 1, 9);
    t->Draw("Q2>>hQ2", "", "goff");
    TH1F *ht1 = new TH1F("ht1", "ht1", numBins, 0, 2);
    t->Draw("-1*t1>>ht1", "", "goff");
    TH1F *hMx2_1 = new TH1F("hMx2_1", "hMx2_1", numBins, -0.5, 1);
    t->Draw("Mx2_1>>hMx2_1", "", "goff");
    TH1F *hphi2 = new TH1F("hphi2", "hphi2", numBins, 0, 2*TMath::Pi());
    t->Draw("phi2>>hphi2", "", "goff");

    TCanvas *c1 = new TCanvas("c1", "c1", 6000, 4500);
    c1->Divide(2,2);

    c1->cd(1);
    hQ2->SetTitle("Fall 2018 DV#pi^{0}P Q^{2}");
    hQ2->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
    hQ2->SetMarkerColor(1);
    hQ2->SetMarkerStyle(21);
    hQ2->Draw("PE");

    c1->cd(2);
    ht1->SetTitle("Fall 2018 DV#pi^{0}P -t");
    ht1->GetXaxis()->SetTitle("-t (GeV^{2})");
    ht1->SetMarkerColor(1);
    ht1->SetMarkerStyle(21);
    ht1->Draw("PE");

    c1->cd(3);
    hMx2_1->SetTitle("Fall 2018 DV#pi^{0}P #pi^{0} Mx^{2}");
    hMx2_1->GetXaxis()->SetTitle("Mx^{2} (GeV^{2})");
    hMx2_1->SetMarkerColor(1);
    hMx2_1->SetMarkerStyle(21);
    hMx2_1->Draw("PE");

    c1->cd(4);
    hphi2->SetTitle("Fall 2018 DV#pi^{0}P #phi");
    hphi2->GetXaxis()->SetTitle("#phi (rad)");
    hphi2->SetMarkerColor(1);
    hphi2->SetMarkerStyle(21);
    hphi2->Draw("PE");

    c1->SaveAs("prac_pi0-kinematics.png");
}

void prac_pi0_asym(TString file, Int_t numBins, Double_t beamPol = 1.0){

    TFile *f = TFile::Open(file);
    TTree *t = (TTree*)f->Get("PhysicsEvents");

    TH1F *hPhi_pos = new TH1F("hPhi_pos", "hPhi_pos", numBins, 0, 2*TMath::Pi());
    t->Draw("phi2>>hPhi_pos", "helicity==+1", "goff");
    TH1F *hPhi_neg = new TH1F("hPhi_neg", "hPhi_neg", numBins, 0, 2*TMath::Pi());
    t->Draw("phi2>>hPhi_neg", "helicity==-1", "goff");

    TH1F *hAdd = new TH1F("hAdd", "hAdd", numBins, 0, 2*TMath::Pi());
    TH1F *hSub = new TH1F("hSub", "hSub", numBins, 0, 2*TMath::Pi());
    hAdd->Add(hPhi_pos, hPhi_neg);
    hSub->Add(hPhi_pos, hPhi_neg, 1, -1);

    TH1F *hAsymPi0 = new TH1F("hAsymPi0", "hAsymPi0", numBins, 0, 2*TMath::Pi());
    hAsymPi0->Divide(hSub, hAdd);
    hAsymPi0->Scale(1./beamPol);

    for (int i=1; i<=numBins; i++) {
        hPhi_pos->SetBinError(i, TMath::Sqrt(hPhi_pos->GetBinContent(i)));
        hPhi_neg->SetBinError(i, TMath::Sqrt(hPhi_neg->GetBinContent(i)));
        hAsymPi0->SetBinError(i, 1/TMath::Sqrt((hAdd->GetBinContent(i))*beamPol));
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 2000, 1500);
    hPhi_pos->SetTitle("DV#pi^{0}P N+ and N- versus #phi");
    hPhi_pos->GetXaxis()->SetTitle("#phi (rad)");
    hPhi_pos->SetMarkerColor(2);
    hPhi_pos->SetMarkerStyle(21);
    hPhi_pos->SetStats(0);

    hPhi_neg->SetMarkerColor(4);
    hPhi_neg->SetMarkerStyle(21);
    hPhi_neg->SetStats(0);

    if (hPhi_pos->GetBinContent(hPhi_pos->GetMaximumBin()) > hPhi_neg->GetBinContent(hPhi_neg->GetMaximumBin())) {
        hPhi_pos->GetYaxis()->SetRangeUser(0.0, hPhi_pos->GetBinContent(hPhi_pos->GetMaximumBin())+2000.);
    }
    else {
        hPhi_pos->GetYaxis()->SetRangeUser(0.0, hPhi_neg->GetBinContent(hPhi_neg->GetMaximumBin())+2000.);
    }

    hPhi_pos->Draw("PE");
    hPhi_neg->Draw("PE sames");
    //gStyle->SetErrorX(0.0); // this sets the error bars to zero on the x axis, not sure if we actually want this or not

    TLegend *leg = new TLegend(0.8,0.7,0.9,0.9);
    leg->AddEntry(hPhi_pos,"N+","ep");
    leg->AddEntry(hPhi_neg,"N-","ep");
    leg->Draw();

    c1->SaveAs("prac_pi0_asym-pos_vs_neg.png");

    TCanvas *c2 = new TCanvas("c2", "c2", 2000, 1500);
    hAsymPi0->SetTitle("DV#pi^{0}P A_{LU} versus #phi");
    hAsymPi0->GetXaxis()->SetTitle("#phi (rad)");
    hAsymPi0->GetYaxis()->SetTitle("#frac{N^{+}-N^{-}}{N^{+}+N^{-}}");
    hAsymPi0->SetMarkerColor(1);
    hAsymPi0->SetMarkerStyle(21);

    TF1 *fit = new TF1("fit", "[0]*sin(x)/(1+[1]*cos(x))+[2]", 0, 2*TMath::Pi());
    fit->SetParameter(0,1);
    fit->SetParameter(1,1);
    fit->SetParameter(2,0);
    fit->SetLineColor(1);

    hAsymPi0->Draw("PE");
    hAsymPi0->Fit("fit");
    gStyle->SetOptFit(0001);

    TLatex *l = new TLatex();
    l->SetTextSize(0.04);
    l->DrawLatexNDC(0.64, 0.53, "([p0]*sin(x))/(1+[p1]*cos(x))+[p2]");

    c2->Update();
    c2->SaveAs("prac_pi0_asym-asym.png");
}