#include <TSystem.h>
#include <cstdio>
#include <iostream>
#include <string>
#include <sstream>

// xB: 0.126, Q2: 1.759, t: -0.670, beamE: 10.604, numBins: 9, beamPol: 0.89

void prac_comp_data_dvcsgen(TString dataFile, TString vggFile, Double_t xB, Double_t q2, Double_t t, Int_t numBins, Double_t E_beam=10.604, Double_t beamPol=0.89){ // lowercase q2 is just so ROOT doesn't crash out

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

    // reading through events and handling them
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

    // making model asymmetry
    Int_t numBinDiv = 10;
    TH1F *hAsymVGG = new TH1F("hAsymVGG", "hAsymVGG", numBins*numBinDiv, 0, 2*TMath::Pi());

    std::ostringstream oss1;
    oss1 << "python -u prac_dvcsgens.py vgg_model " << vggFile << " " << xB << " " << q2 << " " << tpos << " " << numBins << " " << numBinDiv << " " << E_beam;
    std::string posRun = oss1.str();
    FILE* pipe1 = gSystem->OpenPipe(posRun.c_str(), "r");
    if (!pipe1) {
        std::cerr << "Failed to run positive cross section script" << std::endl;
        return;
    }

    ifstream if_vgg(vggFile, ifstream::in);
    Float_t phiVal, xsPosVal, xsNegVal;
    for(int i=1; i <=(numBins*numBinDiv); i++) {
        if_vgg >> phiVal >> xsPosVal >> xsNegVal;
        cout << phiVal*TMath::DegToRad() << " " << hAsymVGG->GetBinCenter(i) << endl;
        hAsymVGG->SetBinContent(i, (xsPosVal-xsNegVal)/(xsPosVal+xsNegVal));
    }
    if_vgg.close();
    
    /*TFile *f_vgg = TFile::Open("vgg_xs_phi-xsPos-xsNeg.txt", "r");
    for (int i=1; i<=(numBins*numBinDiv); i++) {
        
        
    }*/
    // making model asymmetry
    /*Int_t numBinDiv = 10;
    TH1F *hAsymVGG = new TH1F("hAsymVGG", "hAsymVGG", numBins*numBinDiv, 0, 2*TMath::Pi());
    for (int i=1; i<=numBins; i++) {
        for (int j=0; j<numBinDiv; j++) {
            Double_t phi_pos = hAsymData->GetBinCenter(i) - 0.5*hAsymData->GetBinWidth(i) + j*hAsymData->GetBinWidth(i)/numBinDiv;
            // retrieving positively polarized cross section
            std::ostringstream oss1;
            oss1 << "python -u prac_dvcsgens.py vgg_model_xs_pos " << xB << " " << q2 << " " << tpos << " " << phi_pos << " " << E_beam;
            std::string posRun = oss1.str();
            FILE* pipe1 = gSystem->OpenPipe(posRun.c_str(), "r");
            if (!pipe1) {
                std::cerr << "Failed to run positive cross section script" << std::endl;
                return;
            }
            char buffer1[128];
            std::string result1;
            while (fgets(buffer1, sizeof(buffer1), pipe1) != nullptr) {
                result1 += buffer1;
            }
            gSystem->ClosePipe(pipe1);
            Double_t pos_xs = std::stod(result1);

            // retrieving negatively polarized cross section
            std::ostringstream oss2;
            oss2 << "python -u prac_dvcsgen.py vgg_model_xs_neg " << xB << " " << q2 << " " << tpos << " " << phi_pos << " " << E_beam;
            std::string negRun = oss2.str();
            FILE* pipe2 = gSystem->OpenPipe(negRun.c_str(), "r");
            if (!pipe2) {
                std::cerr << "Failed to run negative cross section script" << std::endl;
                return;
            }
            char buffer2[128];
            std::string result2;
            while (fgets(buffer2, sizeof(buffer2), pipe2) != nullptr) {
                result2 += buffer2;
            }
            gSystem->ClosePipe(pipe2);
            Double_t neg_xs = std::stod(result2);

            hAsymVGG->SetBinContent((i-1)*numBinDiv + (j+1), (pos_xs-neg_xs)/(pos_xs+neg_xs));
        }
    }

    // plot and compare
    TCanvas *c1 = new TCanvas("c1", "c1", 2000, 1500);
    hAsymData->SetMarkerColor(4);
    hAsymData->SetLineColor(4);
    hAsymData->SetMarkerStyle(21);
    hAsymData->SetStats(1);

    hAsymVGG->SetMarkerColor(2);
    hAsymVGG->SetLineColor(2);
    hAsymVGG->SetMarkerStyle(21);
    hAsymVGG->SetTitle("A_{LU} from Fall 2018 RG-A and VGG Model");
    hAsymVGG->GetXaxis()->SetTitle("#phi (rad)");
    hAsymVGG->GetYaxis()->SetTitle("#frac{N^{+}-N^{-}}{N^{+}+N^{-}}");
    hAsymVGG->SetStats(0);

    for (Int_t i=1; i<=numBins; i++) {
        hAsymData->SetBinError(i, 1/TMath::Sqrt((hAdd->GetBinContent(i))*beamPol));
    }

    TF1 *fit = new TF1("fit", "[0]*sin(x)/(1+[1]*cos(x))+[2]", 0, 2*TMath::Pi());
    fit->SetParameter(0,1);
    fit->SetParameter(1,1);
    fit->SetParameter(2,0);
    fit->SetLineColor(4);

    hAsymVGG->Draw("P");
    hAsymData->Draw("P sames");
    hAsymData->Fit("fit");
    gStyle->SetOptFit(0001);

    TLegend *leg1 = new TLegend(0.1,0.1,0.2,0.3);
    leg1->AddEntry(hAsymData,"data","ep");
    leg1->AddEntry(hAsymVGG,"vgg","p");
    leg1->Draw();

    TLatex *l = new TLatex();
    l->SetTextSize(0.035);
    l->DrawLatexNDC(0.625, 0.53, "([p0]*sin(x))/(1+[p1]*cos(x))+[p2]");
    l->DrawLatexNDC(0.225, 0.250, Form("x_{B}: %.3f", xB));
    l->DrawLatexNDC(0.225, 0.200, Form("Q^{2}: %.3f", q2));
    l->DrawLatexNDC(0.225, 0.150, Form("-t: %.3f", tpos));

    c1->Update();
    c1->SaveAs("prac_comp-data-dvcsgen_comp.png");*/
}