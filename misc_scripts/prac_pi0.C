// extracing kinematics
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
    hQ2->SetTitle("DV#pi^{0}P Q^{2}");
    hQ2->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
    hQ2->SetMarkerColor(1);
    hQ2->SetMarkerStyle(21);
    hQ2->Draw("PE");

    c1->cd(2);
    ht1->SetTitle("DV#pi^{0}P -t");
    ht1->GetXaxis()->SetTitle("-t (GeV^{2})");
    ht1->SetMarkerColor(1);
    ht1->SetMarkerStyle(21);
    ht1->Draw("PE");

    c1->cd(3);
    hMx2_1->SetTitle("DV#pi^{0}P #pi^{0} Mx^{2}");
    hMx2_1->GetXaxis()->SetTitle("Mx^{2} (GeV^{2})");
    hMx2_1->SetMarkerColor(1);
    hMx2_1->SetMarkerStyle(21);
    hMx2_1->Draw("PE");

    c1->cd(4);
    hphi2->SetTitle("DV#pi^{0}P #phi");
    hphi2->GetXaxis()->SetTitle("#phi (rad)");
    hphi2->SetMarkerColor(1);
    hphi2->SetMarkerStyle(21);
    hphi2->Draw("PE");

    c1->SaveAs("prac_pi0-kinematics.png");
}

// extracing asymmetry
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
    c2->SaveAs("prac_pi0-asym.png");
}

// comparing data with reconstructed MC
void prac_pi0_recon_comp(TString dataFile, TString reconFile, Int_t numBins){
    // opening files
    TFile *f_data = TFile::Open(dataFile);
    TFile *f_recon = TFile::Open(reconFile);

    // extracing trees
    TTree *t_data = (TTree*)f_data->Get("PhysicsEvents");
    TTree *t_recon = (TTree*)f_recon->Get("PhysicsEvents");

    // determining number of events in each file
    Long64_t numData = t_data->GetEntries();
    Long64_t numRecon = t_recon->GetEntries();

    // projecting tree branches onto histograms
    TH1F *hData_Q2 = new TH1F("hData_Q2", "hData_Q2", numBins, 1, 9);
    t_data->Draw("Q2>>hData_Q2","","goff");
    TH1F *hRecon_Q2 = new TH1F("hRecon_Q2", "hRecon_Q2", numBins, 1, 9);
    t_recon->Draw("Q2>>hRecon_Q2","","goff");

    TH1F *hData_t1 = new TH1F("hData_t1", "hData_t1", numBins, 0, 2);
    t_data->Draw("-1*t1>>hData_t1","","goff");
    TH1F *hRecon_t1 = new TH1F("hRecon_t1", "hRecon_t1", numBins, 0, 2);
    t_recon->Draw("-1*t1>>hRecon_t1","","goff");

    TH1F *hData_Mx2 = new TH1F("hData_Mx2", "hData_Mx2", numBins, -0.5, 1);
    t_data->Draw("Mx2>>hData_Mx2","","goff");
    TH1F *hRecon_Mx2 = new TH1F("hRecon_Mx2", "hRecon_Mx2", numBins, -0.5, 1);
    t_recon->Draw("Mx2>>hRecon_Mx2","","goff");

    TH1F *hData_phi2 = new TH1F("hData_phi2", "hData_phi2", numBins, 0, 2*TMath::Pi());
    t_data->Draw("phi2>>hData_phi2","","goff");
    TH1F *hRecon_phi2 = new TH1F("hRecon_phi2", "hRecon_phi2", numBins, 0, 2*TMath::Pi());
    t_recon->Draw("phi2>>hRecon_phi2","","goff");

    TCanvas *c1 = new TCanvas("c1", "c1", 4000, 3000);
    c1->Divide(2,2);

    // Q2
    c1->cd(1);
    hData_Q2->Scale(1./(double)numData);
    hData_Q2->SetMarkerColor(2);
    hData_Q2->SetMarkerStyle(24);
    hData_Q2->SetStats(0);
    hData_Q2->SetTitle("DV#pi^{0}P Q^{2}");
    hData_Q2->GetXaxis()->SetTitle("Q^{2} (Gev^{2})");

    hRecon_Q2->Scale(1/(double)numRecon);
    hRecon_Q2->SetMarkerColor(4);
    hRecon_Q2->SetMarkerStyle(24);

    if (hData_Q2->GetBinContent(hData_Q2->GetMaximumBin()) > hRecon_Q2->GetBinContent(hRecon_Q2->GetMaximumBin())) {
        hData_Q2->GetYaxis()->SetRangeUser(0.0, hData_Q2->GetBinContent(hData_Q2->GetMaximumBin())+0.01);
    }
    else {
        hData_Q2->GetYaxis()->SetRangeUser(0.0, hRecon_Q2->GetBinContent(hRecon_Q2->GetMaximumBin())+0.01);
    }

    hData_Q2->Draw("PE");
    hRecon_Q2->Draw("PE same");

    TLegend *leg1 = new TLegend(0.8,0.4,0.9,0.6);
    leg1->AddEntry(hData_Q2,"data","ep");
    leg1->AddEntry(hRecon_Q2,"mc","ep");
    leg1->Draw();

    // t1
    c1->cd(2);
    hData_t1->Scale(1./(double)numData);
    hData_t1->SetMarkerColor(2);
    hData_t1->SetMarkerStyle(24);
    hData_t1->SetStats(0);
    hData_t1->SetTitle("DV#pi^{0}P -t");
    hData_t1->GetXaxis()->SetTitle("-t (GeV^{2})");

    hRecon_t1->Scale(1./(double)numRecon);
    hRecon_t1->SetMarkerColor(4);
    hRecon_t1->SetMarkerStyle(24);

    if (hData_t1->GetBinContent(hData_t1->GetMaximumBin()) > hRecon_t1->GetBinContent(hRecon_t1->GetMaximumBin())) {
        hData_t1->GetYaxis()->SetRangeUser(0.0, hData_t1->GetBinContent(hData_t1->GetMaximumBin())+0.005);
    }
    else {
        hData_t1->GetYaxis()->SetRangeUser(0.0, hRecon_t1->GetBinContent(hRecon_t1->GetMaximumBin())+0.005);
    }

    hData_t1->Draw("PE");
    hRecon_t1->Draw("PE same");

    TLegend *leg4 = new TLegend(0.8,0.5,0.9,0.7);
    leg4->AddEntry(hData_t1,"data","ep");
    leg4->AddEntry(hRecon_t1,"mc","ep");
    leg4->Draw();

    // Mx2
    c1->cd(3);
    hData_Mx2->Scale(1./(double)numData);
    hData_Mx2->SetMarkerColor(2);
    hData_Mx2->SetMarkerStyle(24);
    hData_Mx2->SetStats(0);
    hData_Mx2->SetTitle("DV#pi^{0}P Mx^{2}");
    hData_Mx2->GetXaxis()->SetTitle("Mx^{2} (GeV^{2})");

    hRecon_Mx2->Scale(1./(double)numRecon);
    hRecon_Mx2->SetMarkerColor(4);
    hRecon_Mx2->SetMarkerStyle(24);

    if (hData_Mx2->GetBinContent(hData_Mx2->GetMaximumBin()) > hRecon_Mx2->GetBinContent(hRecon_Mx2->GetMaximumBin())) {
        hData_Mx2->GetYaxis()->SetRangeUser(0.0, hData_Mx2->GetBinContent(hData_Mx2->GetMaximumBin())+0.05);
    }
    else {
        hData_Mx2->GetYaxis()->SetRangeUser(0.0, hRecon_Mx2->GetBinContent(hRecon_Mx2->GetMaximumBin())+0.05);
    }

    hData_Mx2->Draw("PE");
    hRecon_Mx2->Draw("PE same");

    TLegend *leg2 = new TLegend(0.8,0.4,0.9,0.6);
    leg2->AddEntry(hData_Mx2,"data","ep");
    leg2->AddEntry(hRecon_Mx2,"mc","ep");
    leg2->Draw();

    // phi2
    c1->cd(4);
    hData_phi2->Scale(1./(double)numData);
    hData_phi2->SetMarkerColor(2);
    hData_phi2->SetMarkerStyle(24);
    hData_phi2->SetStats(0);
    hData_phi2->SetTitle("DV#pi^{0}P #phi");
    hData_phi2->GetXaxis()->SetTitle("#phi (rad)");

    hRecon_phi2->Scale(1./(double)numRecon);
    hRecon_phi2->SetMarkerColor(4);
    hRecon_phi2->SetMarkerStyle(24);

    if (hData_phi2->GetBinContent(hData_phi2->GetMaximumBin()) > hRecon_phi2->GetBinContent(hRecon_phi2->GetMaximumBin())) {
        hData_phi2->GetYaxis()->SetRangeUser(0.0, hData_phi2->GetBinContent(hData_phi2->GetMaximumBin())+0.01);
    }
    else {
        hData_phi2->GetYaxis()->SetRangeUser(0.0, hRecon_phi2->GetBinContent(hRecon_phi2->GetMaximumBin())+0.01);
    }

    hData_phi2->Draw();
    hRecon_phi2->Draw("SAME");

    TLegend *leg3 = new TLegend(0.45,0.7,0.55,0.9);
    leg3->AddEntry(hData_phi2,"data","ep");
    leg3->AddEntry(hRecon_phi2,"mc","ep");
    leg3->Draw();

    c1->Update();
    c1->SaveAs("prac_pi0-recon_comp.png");
}

// comparing generated and reconstructed MC to extract bin by bin efficiency
void prac_pi0_efficiency(TString dataFile, TString reconFile, TString genFile, Int_t numBins){
    TFile *f_data = TFile::Open(dataFile);  
    TFile *f_recon = TFile::Open(reconFile);
    TFile *f_gen = TFile::Open(genFile);

    TTree *t_data = (TTree*)f_data->Get("PhysicsEvents");
    TTree *t_recon = (TTree*)f_recon->Get("PhysicsEvents");
    TTree *t_gen = (TTree*)f_gen->Get("PhysicsEvents");

    TH1F *hData_phi2 = new TH1F("hData_phi2", "hData_phi2", numBins, 0, 2*TMath::Pi());
    t_data->Draw("phi2>>hData_phi2","","goff");
    TH1F *hRecon_phi2 = new TH1F("hRecon_phi2", "hRecon_phi2", numBins, 0, 2*TMath::Pi());
    t_recon->Draw("phi2>>hRecon_phi2","","goff");
    TH1F *hGen_phi2 = new TH1F("hGen_phi2", "hGen_phi2", numBins, 0, 2*TMath::Pi());
    t_gen->Draw("phi2>>hGen_phi2","","goff");
    
    TH1F *hRatio = new TH1F("hRatio", "hRatio", numBins, 0, 2*TMath::Pi());
    for (int i=1; i<=12; i++) {
        hRatio->SetBinContent(i, (double)hRecon_phi2->GetBinContent(i)/(double)hGen_phi2->GetBinContent(i));
    }

    TH1F *hData_phi2_scaled = new TH1F("hData_phi2_scaled", "hData_phi2_scaled", numBins, 0, 2*TMath::Pi());
    for (int i=1; i<=12; i++) {
        hData_phi2_scaled->SetBinContent(i, hData_phi2->GetBinContent(i)*(1/hRatio->GetBinContent(i)));
    }

    // summary plot
    TCanvas *c1 = new TCanvas("c1", "c1", 6000, 4500);
    c1->Divide(2,3);

    c1->cd(1);
    hRecon_phi2->SetMarkerColor(1);
    hRecon_phi2->SetMarkerStyle(21);
    hRecon_phi2->SetStats(0);
    hRecon_phi2->SetTitle("#phi Distribution for DV#pi^{0}P Reconstructed MC");
    hRecon_phi2->GetXaxis()->SetTitle("#phi (rad)");
    hRecon_phi2->Draw("PE");

    c1->cd(2);
    hGen_phi2->SetMarkerColor(1);
    hGen_phi2->SetMarkerStyle(21);
    hGen_phi2->SetStats(0);
    hGen_phi2->SetTitle("#phi Distribution for DV#pi^{0}P Generated MC");
    hGen_phi2->GetXaxis()->SetTitle("#phi (rad)");
    hGen_phi2->Draw("PE");

    c1->cd(3);
    hRatio->SetMarkerColor(1);
    hRatio->SetMarkerStyle(21);
    hRatio->SetStats(0);
    hRatio->SetTitle("DV#pi^{0}P Reconstructed to Generated MC Event Efficiency in #phi");
    hRatio->GetXaxis()->SetTitle("#phi (rad)");
    hRatio->Draw("P");

    c1->cd(5);
    hData_phi2->SetMarkerColor(1);
    hData_phi2->SetMarkerStyle(21);
    hData_phi2->SetStats(0);
    hData_phi2->SetTitle("#phi Distribution for DV#pi^{0}P Raw Data");
    hData_phi2->GetXaxis()->SetTitle("#phi (rad)");
    hData_phi2->Draw("PE");

    c1->cd(6);
    hData_phi2_scaled->SetMarkerColor(1);
    hData_phi2_scaled->SetMarkerStyle(21);
    hData_phi2_scaled->SetStats(0);
    hData_phi2_scaled->SetTitle("#phi Distribution for DV#pi^{0}P Scaled Data");
    hData_phi2_scaled->GetXaxis()->SetTitle("#phi (rad)");
    hData_phi2_scaled->Draw("PE");

    c1->SaveAs("prac_pi0-efficiency_1.png");

    TCanvas *c2 = new TCanvas("c2", "c2", 2000, 1500);
    hGen_phi2->Scale(1./(double)hGen_phi2->Integral());
    hGen_phi2->SetMarkerColor(4);
    hGen_phi2->SetMarkerStyle(21);
    hGen_phi2->SetStats(0);
    hGen_phi2->SetTitle("#phi Distribution for DV#pi^{0}P Normalized Generated MC and Scaled Data");
    hGen_phi2->GetXaxis()->SetTitle("#phi (rad)");

    hData_phi2_scaled->Scale(1./(double)hData_phi2_scaled->Integral());
    hData_phi2_scaled->SetMarkerColor(2);
    hData_phi2_scaled->SetMarkerStyle(21);
    hData_phi2_scaled->SetStats(0);
    
    hGen_phi2->Draw("PE");
    hData_phi2_scaled->Draw("PE same");

    if (hData_phi2_scaled->GetBinContent(hData_phi2_scaled->GetMaximumBin()) > hGen_phi2->GetBinContent(hGen_phi2->GetMaximumBin())) {
        hGen_phi2->GetYaxis()->SetRangeUser(0.0, hData_phi2_scaled->GetBinContent(hData_phi2_scaled->GetMaximumBin())+0.005);
        if (hData_phi2_scaled->GetBinContent(hData_phi2_scaled->GetMinimumBin()) < hGen_phi2->GetBinContent(hGen_phi2->GetMinimumBin())){
            hGen_phi2->GetYaxis()->SetRangeUser(hData_phi2_scaled->GetBinContent(hData_phi2_scaled->GetMinimumBin())-0.005, hData_phi2_scaled->GetBinContent(hData_phi2_scaled->GetMaximumBin())+0.005);
        }
        else {
            hGen_phi2->GetYaxis()->SetRangeUser(hGen_phi2->GetBinContent(hGen_phi2->GetMinimumBin())-0.005, hData_phi2_scaled->GetBinContent(hData_phi2_scaled->GetMaximumBin())+0.005);
        }
    }
    else {
        hGen_phi2->GetYaxis()->SetRangeUser(0.0, hGen_phi2->GetBinContent(hGen_phi2->GetMaximumBin())+0.005);
        if (hData_phi2_scaled->GetBinContent(hData_phi2_scaled->GetMinimumBin()) < hGen_phi2->GetBinContent(hGen_phi2->GetMinimumBin())){
            hGen_phi2->GetYaxis()->SetRangeUser(hData_phi2_scaled->GetBinContent(hData_phi2_scaled->GetMinimumBin())-0.005, hGen_phi2->GetBinContent(hGen_phi2->GetMaximumBin())+0.005);
        }
        else {
            hGen_phi2->GetYaxis()->SetRangeUser(hGen_phi2->GetBinContent(hGen_phi2->GetMinimumBin())-0.005, hGen_phi2->GetBinContent(hGen_phi2->GetMaximumBin())+0.005);
        }
    }

    TLegend *leg1 = new TLegend(0.45,0.1,0.55,0.3);
    leg1->AddEntry(hGen_phi2,"mc","ep");
    leg1->AddEntry(hData_phi2_scaled,"data","ep");
    leg1->Draw();

    c2->SaveAs("prac_pi0-efficiency_2.png");
}

// dvcsPionFile refers to running the pi0 data through the processing_dvcs analysis script
// dvcsFile refers to running the actual data through the processing_dvcs analysis script (ie the number of detected DVCS events)
void prac_pi0_contamination(TString dataFile, TString reconFile, TString dvcsPionFile, TString dvcsFile, Double_t xB, Double_t q2, Double_t t, Int_t numBins) {
    TFile *f_data = TFile::Open(dataFile);
    TFile *f_recon = TFile::Open(reconFile);
    TFile *f_dvcsPion = TFile::Open(dvcsPionFile);
    TFile *f_dvcs = TFile::Open(dvcsFile);

    TTree *t_data = (TTree*)f_data->Get("PhysicsEvents");
    TTree *t_recon = (TTree*)f_recon->Get("PhysicsEvents");
    TTree *t_dvcsPion = (TTree*)f_dvcsPion->Get("PhysicsEvents");
    TTree *t_dvcs = (TTree*)f_dvcs->Get("PhysicsEvents");

    Double_t xB_min, xB_max, Q2_min, Q2_max, t_min, t_max, tpos, tpos_min, tpos_max;
    xB_min = xB - 0.01;
    xB_max = xB + 0.01;
    Q2_min = q2 - 0.1;
    Q2_max = q2 + 0.1;
    t_min = t - 0.005;
    t_max = t + 0.005;
    tpos = -1*t;
    tpos_min = -1*t_max;
    tpos_max = -1*t_min;

    TString cuts = Form("x > %0.3f && x < %0.3f && Q2 > %0.3f && Q2 < %0.3f && t1 > %0.3f && t1 < %0.3f", xB_min, xB_max, Q2_min, Q2_max, t_min, t_max);

    TH1F *hData = new TH1F("hData", "hData", numBins, 0, 2*TMath::Pi());
    //t_data->Draw("phi2>>hData", "x > 0.125 && x < 0.127 && Q2 > 1.69 && Q2 < 1.89 && t1 > -0.720 && t1 < -0.520", "goff");
    t_data->Draw("phi2>>hData", cuts, "goff");
    TH1F *hRecon = new TH1F("hRecon", "hRecon", numBins, 0, 2*TMath::Pi());
    //t_recon->Draw("phi2>>hRecon", "x > 0.125 && x < 0.127 && Q2 > 1.69 && Q2 < 1.89 && t1 > -0.720 && t1 < -0.520", "goff");
    t_recon->Draw("phi2>>hRecon", cuts, "goff");
    TH1F *hDVCSPion = new TH1F("hDVCSPion", "hDVCSPion", numBins, 0, 2*TMath::Pi());
    //t_dvcsPion->Draw("phi2>>hDVCSPion", "x > 0.125 && x < 0.127 && Q2 > 1.69 && Q2 < 1.89 && t1 > -0.720 && t1 < -0.520", "goff");
    t_dvcsPion->Draw("phi2>>hDVCSPion", cuts, "goff");    
    TH1F *hDVCS = new TH1F("hDVCS", "hDVCS", numBins, 0, 2*TMath::Pi());
    //t_dvcs->Draw("phi2>>hDVCS", "x > 0.125 && x < 0.127 && Q2 > 1.69 && Q2 < 1.89 && t1 > -0.720 && t1 < -0.520", "goff");
    t_dvcs->Draw("phi2>>hDVCS", cuts, "goff");

    TH1F *hDiv1 = new TH1F("hDiv1", "hDiv1", numBins, 0, 2*TMath::Pi());
    hDiv1->Divide(hDVCSPion, hRecon);

    TH1F *hDiv2 = new TH1F("hDiv2", "hDiv2", numBins, 0, 2*TMath::Pi());
    hDiv2->Divide(hDiv1, hDVCS);

    TH1F *hCont = new TH1F("hCont", "hCont", numBins, 0, 2*TMath::Pi());
    hCont->Multiply(hData, hDiv2);

    for (int i=1; i<=numBins; i++) {
        if (hCont->GetBinContent(i) == 0) {
            hCont->SetBinContent(i, 0.0);
        }
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 2000, 1500);
    hCont->SetTitle("DV#pi^{0}P Contamination Constants");
    hCont->GetXaxis()->SetTitle("#phi (rad)");
    hCont->GetYaxis()->SetTitle("c_{i}");
    hCont->SetMarkerColor(1);
    hCont->SetMarkerStyle(21);
    hCont->SetStats(0);
    hCont->GetYaxis()->SetRangeUser(-0.01, hCont->GetBinContent(hCont->GetMaximumBin())+0.05);
    hCont->Draw("P");

    TLatex *l = new TLatex();
    l->SetTextSize(0.035);
    l->DrawLatexNDC(0.625, 0.53, "([p0]*sin(x))/(1+[p1]*cos(x))+[p2]");
    l->DrawLatexNDC(0.7, 0.900, Form("x_{B}: %.3f", xB));
    l->DrawLatexNDC(0.7, 0.850, Form("Q^{2}: %.3f", q2));
    l->DrawLatexNDC(0.7, 0.800, Form("-t: %.3f", tpos));

    c1->SaveAs("prac_pi0-contamination.png");
}

// extracts this for specific xB, q2, t kinematics so it can be directly applied to the DVCS asymmetry for pi0 subtraction
void prac_pi0_contamination_return(TString dataFile, TString reconFile, TString dvcsFile, Double_t xB, Double_t q2, Double_t t, Int_t numBins) {
    // kinematic settings
    Double_t xB_min, xB_max, Q2_min, Q2_max, t_min, t_max, tpos, tpos_min, tpos_max;
    xB_min = xB - 0.02;
    xB_max = xB + 0.02;
    Q2_min = q2 - 0.2;
    Q2_max = q2 + 0.2;
    t_min = t - 0.1;
    t_max = t + 0.1;
    tpos = -1*t;
    tpos_min = -1*t_max;
    tpos_max = -1*t_min;

    TFile *f_data = TFile::Open(dataFile);
    TFile *f_recon = TFile::Open(reconFile);
    TFile *f_dvcs = TFile::Open(dvcsFile);

    TTree *t_data = (TTree*)f_data->Get("PhysicsEvents");
    TTree *t_recon = (TTree*)f_recon->Get("PhysicsEvents");
    TTree *t_dvcs = (TTree*)f_dvcs->Get("PhysicsEvents");
}