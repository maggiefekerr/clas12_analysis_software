#include <TSystem.h>
#include <cstdio>
#include <iostream>
#include <string>
#include <sstream>

// xB: 0.126, Q2: 1.759, t: -0.670, beamE: 10.604, numBins: 9, beamPol: 0.89
// eventually want to reimplement the KM15 code to read out to a .txt file as well for consistency which is why I defined
// a TString km15File here

void prac_comp_data_models(TString dataFile, TString vggFile, TString km15File, Double_t xB, Double_t q2, Double_t t, Int_t numBins, Double_t E_beam=10.604, Double_t beamPol=0.89){
    // data
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

    TFile *f = TFile::Open(dataFile);
    TH1F *hPhi_pos = new TH1F("hPhi_pos", "hPhi_pos", numBins, 0, 2*TMath::Pi());
    TH1F *hPhi_neg = new TH1F("hPhi_neg", "hPhi_neg", numBins, 0, 2*TMath::Pi());

    TTreeReader r("PhysicsEvents", f);
    TTreeReaderValue<double> phi2(r, "phi2");
    TTreeReaderValue<int> helicity(r, "helicity");
    TTreeReaderValue<double> x(r, "x");
    TTreeReaderValue<double> Q2(r, "Q2");
    TTreeReaderValue<double> t1(r, "t1");

    while (r.Next()) {
        if ((xB_min < *x < xB_max) && (Q2_min < *Q2 < Q2_max) && (tpos_min < *t1*-1 < tpos_max)){
            if (*helicity == +1) {
                hPhi_pos->AddBinContent(hPhi_pos->FindBin(*phi2));
            }
            if (*helicity == -1) {
                hPhi_neg->AddBinContent(hPhi_pos->FindBin(*phi2));
            }
        }
    }

    TH1F *hAdd = new TH1F("hAdd", "hAdd", numBins, 0, 2*TMath::Pi());
    hAdd->Add(hPhi_pos, hPhi_neg);
    TH1F *hSub = new TH1F("hSub", "hSub", numBins, 0, 2*TMath::Pi());
    hSub->Add(hPhi_pos, hPhi_neg, 1, -1);
    TH1F *hAsymData = new TH1F("hAsymData", "hAsymData", numBins, 0, 2*TMath::Pi());
    hAsymData->Divide(hSub, hAdd);
    hAsymData->Scale(1./beamPol);

    // models
    Int_t numBinDiv = 10;

    // vgg
    TH1F *hAsymVGG = new TH1F("hAsymVGG", "hAsymVGG", numBins*numBinDiv, 0, 2*TMath::Pi());

    std::ostringstream oss1;
    oss1 << "python -u prac_dvcsgens.py vgg_model " << vggFile << " " << xB << " " << q2 << " " << tpos << " " << numBins << " " << numBinDiv << " " << E_beam;
    std::string asyRun1 = oss1.str();
    FILE* pipe1 = gSystem->OpenPipe(asyRun1.c_str(), "r");
    if (!pipe1) {
        std::cerr << "Failed to run script" << std::endl;
        return;
    }

    ifstream if_vgg(vggFile, ifstream::in);
    Double_t phiVal1, xsPosVal, xsNegVal;
    for(int i=1; i <=(numBins*numBinDiv); i++) {
        if_vgg >> phiVal1 >> xsPosVal >> xsNegVal;
        hAsymVGG->SetBinContent(i, -1*(xsPosVal-xsNegVal)/(xsPosVal+xsNegVal));
    }
    if_vgg.close();

    // km15
    TH1F *hAsymKM15 = new TH1F("hAsymKM15", "hAsymKM15", numBins*numBinDiv, 0, 2*TMath::Pi());
    
    std::ostringstream oss2;
    oss2 << "python -u prac_gepard.py km15_model " << km15File << " " << xB << " " << q2 << " " << tpos << " " << numBins << " " << numBinDiv << " " << E_beam << " " << "ALU";
    std::string asyRun2 = oss2.str();
    FILE* pipe2 = gSystem->OpenPipe(asyRun2.c_str(), "r");
    if (!pipe2) {
        std::cerr << "Failed to run script" << std::endl;
        return;
    }

    ifstream if_km15(km15File, ifstream::in);
    Double_t phiVal2, asyVal;
    for(int i=1; i <=(numBins*numBinDiv); i++) {
        if_km15 >> phiVal2 >> asyVal;
        //cout << phiVal2 << " " << asyVal << endl;
        hAsymKM15->SetBinContent(i, asyVal);
    }
    if_km15.close();

    cout << "here and alive" << endl;

    // plot and compare
    TCanvas *c1 = new TCanvas("c1", "c1", 2000, 1500);
    hAsymData->SetMarkerColor(1);
    hAsymData->SetLineColor(1);
    hAsymData->SetMarkerStyle(21);
    hAsymData->SetStats(1);

    hAsymVGG->SetMarkerColor(2);
    hAsymVGG->SetLineColor(2);
    hAsymVGG->SetMarkerStyle(21);
    hAsymVGG->SetTitle("A_{LU} from Fall 2018 RG-A vs VGG and KM15 Models");
    hAsymVGG->GetXaxis()->SetTitle("#phi (rad)");
    hAsymVGG->GetYaxis()->SetTitle("#frac{N^{+}-N^{-}}{N^{+}+N^{-}}");
    hAsymVGG->SetStats(0);

    hAsymKM15->SetMarkerColor(4);
    hAsymKM15->SetLineColor(4);
    hAsymKM15->SetMarkerStyle(21);
    hAsymKM15->SetStats(0);

    for (Int_t i=1; i<=numBins; i++) {
        hAsymData->SetBinError(i, 1/TMath::Sqrt((hAdd->GetBinContent(i))*beamPol));
    }

    TF1 *fit = new TF1("fit", "[0]*sin(x)/(1+[1]*cos(x))+[2]", 0, 2*TMath::Pi());
    fit->SetParameter(0,1);
    fit->SetParameter(1,1);
    fit->SetParameter(2,0);
    fit->SetLineColor(1);

    hAsymVGG->Draw("P");
    hAsymKM15->Draw("P sames");
    hAsymData->Draw("P sames");
    hAsymData->Fit("fit");
    gStyle->SetOptFit(0001);

    TLegend *leg1 = new TLegend(0.1,0.1,0.2,0.3);
    leg1->AddEntry(hAsymData,"data","ep");
    leg1->AddEntry(hAsymVGG,"vgg","p");
    leg1->AddEntry(hAsymKM15,"km15","p");
    leg1->Draw();

    TLatex *l = new TLatex();
    l->SetTextSize(0.035);
    l->DrawLatexNDC(0.625, 0.53, "([p0]*sin(x))/(1+[p1]*cos(x))+[p2]");
    l->DrawLatexNDC(0.225, 0.250, Form("x_{B}: %.3f", xB));
    l->DrawLatexNDC(0.225, 0.200, Form("Q^{2}: %.3f", q2));
    l->DrawLatexNDC(0.225, 0.150, Form("-t: %.3f", tpos));

    c1->Update();
    c1->SaveAs("prac_comp-data-models.png");
}