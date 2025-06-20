// extracting bin-by-bin scaling ratio by comparing reconstructed and generated mc
void efficiency(TString proc, TString gen, TString data, Int_t numBins) {
    TFile *f_proc = TFile::Open(proc);
    TFile *f_gen = TFile::Open(gen);
    TFile *f_data = TFile::Open(data);

    TTree *t_proc = (TTree*)f_proc->Get("PhysicsEvents");
    TTree *t_gen = (TTree*)f_gen->Get("PhysicsEvents");
    TTree *t_data = (TTree*)f_data->Get("PhysicsEvents");

    TH1F *hProc = new TH1F("hProc", "hProc", numBins, 0, 2*TMath::Pi());
    t_proc->Draw("phi2>>hProc","","goff");
    TH1F *hGen = new TH1F("hGen", "hGen", numBins, 0, 2*TMath::Pi());
    t_gen->Draw("phi2>>hGen","","goff");    
    TH1F *hData = new TH1F("hData", "hData", numBins, 0, 2*TMath::Pi());
    t_data->Draw("phi2>>hData","","goff");

    TH1F *hRatio = new TH1F("hRatio", "hRatio", numBins, 0, 2*TMath::Pi());
    for (int i=1; i<=numBins; i++) {
        hRatio->SetBinContent(i, (double)hProc->GetBinContent(i)/(double)hGen->GetBinContent(i));
    }

    TH1F *hData_scaled = new TH1F("hData_scaled", "hData_scaled", numBins, 0, 2*TMath::Pi());
    for (int i=1; i<=numBins; i++) {
        hData_scaled->SetBinContent(i, hData->GetBinContent(i)*(1/hRatio->GetBinContent(i)));
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 6000, 4500);
    c1->Divide(2,3);

    c1->cd(1);
    hProc->SetMarkerColor(1);
    hProc->SetMarkerStyle(21);
    hProc->SetStats(0);
    hProc->SetTitle("#phi Distributions for Processed DVCS MC (Fall 2018)");
    hProc->GetXaxis()->SetTitle("#phi (rad)");
    hProc->Draw("PE");

    c1->cd(2);
    hGen->SetMarkerColor(1);
    hGen->SetMarkerStyle(21);
    hGen->SetStats(0);
    hGen->SetTitle("#phi Distributions for Generated DVCS MC (Fall 2018)");
    hGen->GetXaxis()->SetTitle("#phi (rad)");
    hGen->Draw("PE");

    c1->cd(3);
    hRatio->SetMarkerColor(1);
    hRatio->SetMarkerStyle(21);
    hRatio->SetStats(0);
    hRatio->SetTitle("Processed to Generated DVCS MC Event Efficiency in #phi (Fall 2018)");
    hRatio->GetXaxis()->SetTitle("#phi (rad)");
    hRatio->Draw("P");

    c1->cd(5);
    hData->SetMarkerColor(1);
    hData->SetMarkerStyle(21);
    hData->SetStats(0);
    hData->SetTitle("#phi Distributions for Raw DVCS Data (Fall 2018)");
    hData->GetXaxis()->SetTitle("#phi (rad)");
    hData->Draw("PE");

    c1->cd(6);
    hData_scaled->SetMarkerColor(1);
    hData_scaled->SetMarkerStyle(21);
    hData_scaled->SetStats(0);
    hData_scaled->SetTitle("#phi Distributions for Scaled DVCS Data (Fall 2018)");
    hData_scaled->GetXaxis()->SetTitle("#phi (rad)");
    hData_scaled->Draw("PE");

    c1->SaveAs("prac_16-06-2025_efficiency.png");

    TCanvas *c2 = new TCanvas("c2", "c2", 2000, 1500);
    hGen->Scale(1./(double)hGen->Integral());
    hGen->SetMarkerColor(4);
    hGen->SetMarkerStyle(21);
    hGen->SetStats(0);
    hGen->SetTitle("#phi Distributions for Normalized Generated DVCS MC and Scaled DVCS Data (Fall 2018)");
    hGen->GetXaxis()->SetTitle("#phi (rad)");
    hData_scaled->Scale(1./(double)hData_scaled->Integral());
    hData_scaled->SetMarkerColor(2);
    hData_scaled->SetMarkerStyle(21);
    
    hGen->Draw("PE");
    hData_scaled->Draw("SAME PE");

    TLegend *leg1 = new TLegend(0.45,0.7,0.55,0.9);
    leg1->AddEntry(hGen,"mc","ep");
    leg1->AddEntry(hData_scaled,"data","ep");
    leg1->Draw();

    c2->SaveAs("prac_16-06-2025_comp.png");
}

// extracting bin-by-bin scaling ratio by comparing reconstructed and generated mc; distinguishing helicity for easier application to asymmetries
void efficiency_pos_neg(TString proc, TString gen, TString data, Int_t numBins) {
    TFile *f_proc = TFile::Open(proc);
    TFile *f_gen = TFile::Open(gen);
    TFile *f_data = TFile::Open(data);

    TTree *t_proc = (TTree*)f_proc->Get("PhysicsEvents");
    TTree *t_gen = (TTree*)f_gen->Get("PhysicsEvents");
    TTree *t_data = (TTree*)f_data->Get("PhysicsEvents");

    TH1F *hProc_pos = new TH1F("hProc_pos", "hProc_pos", numBins, 0, 2*TMath::Pi());
    TH1F *hProc_neg = new TH1F("hProc_neg", "hProc_neg", numBins, 0, 2*TMath::Pi());
    t_proc->Draw("phi2>>hProc_pos", "helicity==+1", "goff");
    t_proc->Draw("phi2>>hProc_neg", "helicity==-1", "goff");

    TH1F *hGen_pos = new TH1F("hGen_pos", "hGen_pos", numBins, 0, 2*TMath::Pi());
    TH1F *hGen_neg = new TH1F("hGen_neg", "hGen_neg", numBins, 0, 2*TMath::Pi());
    t_gen->Draw("phi2>>hGen_pos", "helicity==+1", "goff");
    t_gen->Draw("phi2>>hGen_neg", "helicity==-1", "goff");

    TH1F *hData_pos = new TH1F("hData_pos", "hData_pos", numBins, 0, 2*TMath::Pi());
    TH1F *hData_neg = new TH1F("hData_neg", "hData_neg", numBins, 0, 2*TMath::Pi());
    t_data->Draw("phi2>>hData_pos", "helicity==+1", "goff");
    t_data->Draw("phi2>>hData_neg", "helicity==-1", "goff");

    TH1F *hRatio_pos = new TH1F("hRatio_pos", "hRatio_pos", numBins, 0, 2*TMath::Pi());
    TH1F *hRatio_neg = new TH1F("hRatio_neg", "hRatio_neg", numBins, 0, 2*TMath::Pi());
    for (int i=1; i<=numBins; i++) {
        hRatio_pos->SetBinContent(i, (double)hProc_pos->GetBinContent(i)/(double)hGen_pos->GetBinContent(i));
        hRatio_neg->SetBinContent(i, (double)hProc_neg->GetBinContent(i)/(double)hGen_neg->GetBinContent(i));
    }

    TH1F *hData_scaled_pos = new TH1F("hData_scaled_pos", "hData_scaled_pos", numBins, 0, 2*TMath::Pi());
    TH1F *hData_scaled_neg = new TH1F("hData_scaled_neg", "hData_scaled_neg", numBins, 0, 2*TMath::Pi());
    for (int i=1; i<=numBins; i++) {
        hData_scaled_pos->SetBinContent(i, hData_pos->GetBinContent(i)*(1/hRatio_pos->GetBinContent(i)));
        hData_scaled_neg->SetBinContent(i, hData_neg->GetBinContent(i)*(1/hRatio_neg->GetBinContent(i)));
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 6000, 4500);
    c1->Divide(2,3);

    c1->cd(1);
    hProc_pos->SetMarkerColor(2);
    hProc_pos->SetMarkerStyle(21);
    hProc_pos->SetStats(0);
    hProc_pos->SetTitle("#phi Distributions for Processed DVCS MC (Fall 2018)");
    hProc_pos->GetXaxis()->SetTitle("#phi (rad)");

    hProc_neg->SetMarkerColor(4);
    hProc_neg->SetMarkerStyle(21);
    hProc_neg->SetStats(0);

    if (hProc_pos->GetBinContent(hProc_pos->GetMaximumBin()) > hProc_neg->GetBinContent(hProc_neg->GetMaximumBin())) {
        hProc_pos->GetYaxis()->SetRangeUser(0.0, hProc_pos->GetBinContent(hProc_pos->GetMaximumBin())+1000.);
    }
    else {
        hProc_pos->GetYaxis()->SetRangeUser(0.0, hProc_neg->GetBinContent(hProc_neg->GetMaximumBin())+1000.);
    }

    hProc_pos->Draw("PE");
    hProc_neg->Draw("PE sames");

    TLegend *leg1 = new TLegend(0.8,0.7,0.9,0.9);
    leg1->AddEntry(hProc_pos,"N+","ep");
    leg1->AddEntry(hProc_neg,"N-","ep");
    leg1->Draw();
    c1->Update();

    c1->cd(2);
    hGen_pos->SetMarkerColor(2);
    hGen_pos->SetMarkerStyle(21);
    hGen_pos->SetStats(0);
    hGen_pos->SetTitle("#phi Distributions for Processed DVCS MC (Fall 2018)");
    hGen_pos->GetXaxis()->SetTitle("#phi (rad)");

    hGen_neg->SetMarkerColor(4);
    hGen_neg->SetMarkerStyle(21);
    hGen_neg->SetStats(0);

    if (hGen_pos->GetBinContent(hGen_pos->GetMaximumBin()) > hGen_neg->GetBinContent(hGen_neg->GetMaximumBin())) {
        hGen_pos->GetYaxis()->SetRangeUser(0.0, hGen_pos->GetBinContent(hGen_pos->GetMaximumBin())+2000.);
    }
    else {
        hGen_pos->GetYaxis()->SetRangeUser(0.0, hGen_neg->GetBinContent(hGen_neg->GetMaximumBin())+2000.);
    }

    hGen_pos->Draw("PE");
    hGen_neg->Draw("PE sames");

    TLegend *leg2 = new TLegend(0.8,0.7,0.9,0.9);
    leg2->AddEntry(hGen_pos,"N+","ep");
    leg2->AddEntry(hGen_neg,"N-","ep");
    leg2->Draw();
    c1->Update();

    c1->cd(3);
    hRatio_pos->SetMarkerColor(2);
    hRatio_pos->SetMarkerStyle(21);
    hRatio_pos->SetStats(0);
    hRatio_pos->SetTitle("#phi Distributions for Processed DVCS MC (Fall 2018)");
    hRatio_pos->GetXaxis()->SetTitle("#phi (rad)");

    hRatio_neg->SetMarkerColor(4);
    hRatio_neg->SetMarkerStyle(21);
    hRatio_neg->SetStats(0);

    if (hRatio_pos->GetBinContent(hRatio_neg->GetMaximumBin()) > hRatio_neg->GetBinContent(hRatio_neg->GetMaximumBin())) {
        hRatio_pos->GetYaxis()->SetRangeUser(0.0, hRatio_pos->GetBinContent(hRatio_pos->GetMaximumBin())+2000.);
    }
    else {
        hRatio_pos->GetYaxis()->SetRangeUser(0.0, hRatio_neg->GetBinContent(hRatio_neg->GetMaximumBin())+2000.);
    }

    hRatio_pos->Draw("PE");
    hRatio_neg->Draw("PE sames");

    TLegend *leg3 = new TLegend(0.8,0.7,0.9,0.9);
    leg3->AddEntry(hRatio_pos,"N+","ep");
    leg3->AddEntry(hRatio_neg,"N-","ep");
    leg3->Draw();
    c1->Update();

    /*hRatio->SetMarkerColor(1);
    hRatio->SetMarkerStyle(21);
    hRatio->SetStats(0);
    hRatio->SetTitle("Processed to Generated DVCS MC Event Efficiency in #phi (Fall 2018)");
    hRatio->GetXaxis()->SetTitle("#phi (rad)");
    hRatio->Draw("P");

    c1->cd(5);
    hData->SetMarkerColor(1);
    hData->SetMarkerStyle(21);
    hData->SetStats(0);
    hData->SetTitle("#phi Distributions for Raw DVCS Data (Fall 2018)");
    hData->GetXaxis()->SetTitle("#phi (rad)");
    hData->Draw("PE");

    c1->cd(6);
    hData_scaled->SetMarkerColor(1);
    hData_scaled->SetMarkerStyle(21);
    hData_scaled->SetStats(0);
    hData_scaled->SetTitle("#phi Distributions for Scaled DVCS Data (Fall 2018)");
    hData_scaled->GetXaxis()->SetTitle("#phi (rad)");
    hData_scaled->Draw("PE");

    c1->SaveAs("prac_16-06-2025_efficiency.png");

    TCanvas *c2 = new TCanvas("c2", "c2", 2000, 1500);
    // for this one we could maybe have two canvases and have one on each
    hGen->Scale(1./(double)hGen->Integral());
    hGen->SetMarkerColor(4);
    hGen->SetMarkerStyle(21);
    hGen->SetStats(0);
    hGen->SetTitle("#phi Distributions for Normalized Generated DVCS MC and Scaled DVCS Data (Fall 2018)");
    hGen->GetXaxis()->SetTitle("#phi (rad)");
    hData_scaled->Scale(1./(double)hData_phi2_scaled->Integral());
    hData_scaled->SetMarkerColor(2);
    hData_scaled->SetMarkerStyle(21);
    
    hGen->Draw("PE");
    hData_scaled->Draw("SAME PE");

    TLegend *leg1 = new TLegend(0.45,0.7,0.55,0.9);
    leg1->AddEntry(hGen,"mc","ep");
    leg1->AddEntry(hData_phi2_scaled,"data","ep");
    leg1->Draw();

    c2->SaveAs("prac_16-06-2025_comp.png");*/
}