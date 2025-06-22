// hoping for this to sort of a one-stop code for a lot of the asymmetry calculations
// measured asymmetry
// models
// measured asymmetry with both of the models (individually and collectively)
// pi0 contamination
// pi0 subtraction
// comparing pure ep->ep\gamma asymmetry to both of the models
#include <TSystem.h>
#include <cstdio>
#include <iostream>
#include <string>
#include <sstream>

void prac_ALU_measured_asymmetry(TH1F* hDVCS_measAsym, TF1* fit_measAsym) {
    TCanvas *c1 = new TCanvas("c1", "c1", 2000, 1500);

    hDVCS_measAsym->Draw();
}

void prac_ALU_pi0_asymmetry(TH1F* hPi0_Asym, TF1* fit_pi0Asym) {

}

void prac_ALU_dvcs_asymmetry(TH1F* hDVCS_measAsym, TH1F* hPi0_Asym, TH1F* hDVCS_Asym, TF1* fit_measAsym, TF1* fit_pi0, TF1* fit_dvcsAsym){
    TCanvas *c1 = new TCanvas("c1", "c1", 2000, 1500);
    hDVCS_measAsym->SetMarkerColor(2);
    hDVCS_measAsym->SetLineColor(2);
    hDVCS_measAsym->SetMarkerStyle(24);
    hDVCS_measAsym->SetStats(0);
    hDVCS_measAsym->SetTitle("A_{LU}");
    hDVCS_measAsym->GetXaxis()->SetTitle("#phi (rad)");
    hDVCS_measAsym->GetYaxis()->SetTitle("A_{LU}");

    hPi0_Asym->SetMarkerColor(4);
    hPi0_Asym->SetLineColor(4);
    hPi0_Asym->SetMarkerStyle(26);
    hPi0_Asym->SetStats(0);

    hDVCS_Asym->SetMarkerColor(1);
    hDVCS_Asym->SetLineColor(1);
    hDVCS_Asym->SetMarkerStyle(25);
    hDVCS_Asym->SetStats(1);  

    hDVCS_measAsym->Draw("P");
    hPi0_Asym->Draw("P sames");
    hDVCS_Asym->Draw("P sames");
    hDVCS_measAsym->Fit("fit_measAsym");
    hPi0_Asym->Fit("fit_pi0");
    hDVCS_Asym->Fit("fit_dvcsAsym");
    gStyle->SetOptFit(0001);
}

void prac_ALU(TString dvcsFile, TString pi0File, TString vggFile, TString km15File, TString pi0ReconFile, TString pi0DVCSFile, Double_t xB, Double_t q2, Double_t t, Int_t numBins, Double_t beamE=10.604, Double_t beamPol=1.0) {
    // opening files
    TFile *f_dvcs = TFile::Open(dvcsFile);
    TFile *f_pi0 = TFile::Open(pi0File);
    TFile *f_pi0Recon = TFile::Open(pi0ReconFile);
    TFile *f_pi0DVCS = TFile::Open(pi0DVCSFile);

    // defining trees
    TTree *t_dvcs = (TTree*)f_dvcs->Get("PhysicsEvents");
    TTree *t_pi0 = (TTree*)f_pi0->Get("PhysicsEvents");
    TTree *t_pi0Recon = (TTree*)f_pi0Recon->Get("PhysicsEvents");
    TTree *t_pi0DVCS = (TTree*)f_pi0DVCS->Get("PhysicsEvents");

    // kinematic ranges
    Double_t xB_min, xB_max, Q2_min, Q2_max, t_min, t_max, tpos, tpos_min, tpos_max;
    xB_min = xB - 0.1;
    xB_max = xB + 0.1;
    Q2_min = q2 - 0.2;
    Q2_max = q2 + 0.2;
    t_min = t - 0.1;
    t_max = t + 0.1;
    tpos = -1*t;
    tpos_min = -1*t_max;
    tpos_max = -1*t_min;
    TString cuts = Form("x > %0.3f && x < %0.3f && Q2 > %0.3f && Q2 < %0.3f && t1 > %0.3f && t1 < %0.3f", xB_min, xB_max, Q2_min, Q2_max, t_min, t_max);

    // measured asymmetry (without pi0 subtraction)
    TH1F* hDVCS_phiPos = new TH1F("hDVCS_phiPos", "hDVCS_phiPos", numBins, 0, 2*TMath::Pi());
    t_dvcs->Draw("phi2>>hDVCS_phiPos", "helicity == +1 && " + cuts, "goff");
    TH1F* hDVCS_phiNeg = new TH1F("hDVCS_phiNeg", "hDVCS_phiNeg", numBins, 0, 2*TMath::Pi());
    t_dvcs->Draw("phi2>>hDVCS_phiNeg", "helicity == -1 && " + cuts, "goff");

    TH1F *hAdd1 = new TH1F("hAdd1", "hAdd1", numBins, 0, 2*TMath::Pi());
    hAdd1->Add(hDVCS_phiPos, hDVCS_phiNeg);
    TH1F *hSub1 = new TH1F("hSub1", "hSub1", numBins, 0, 2*TMath::Pi());
    hSub1->Add(hDVCS_phiPos, hDVCS_phiNeg, 1, -1);
    TH1F *hDVCS_measAsym = new TH1F("hDVCS_measAsym", "hDVCS_measAsym", numBins, 0, 2*TMath::Pi());
    hDVCS_measAsym->Divide(hSub1, hAdd1);
    hDVCS_measAsym->Scale(1./beamPol);
    for (Int_t i=1; i<=numBins; i++) {
        hDVCS_measAsym->SetBinError(i, 1/TMath::Sqrt((hAdd1->GetBinContent(i))*beamPol));
    }
    TF1 *fit_measAsym = new TF1("fit_measAsym", "[0]*sin(x)/(1+[1]*cos(x))+[2]", 0, 2*TMath::Pi());
    fit_measAsym->SetParameter(0,1);
    fit_measAsym->SetParameter(1,1);
    fit_measAsym->SetParameter(2,0);
    //prac_ALU_measured_asymmetry(hDVCS_measAsym, fit_measAsym);

    // pion asymmetry
    TH1F* hPi0_phiPos = new TH1F("hPi0_phiPos", "hPi0_phiPos", numBins, 0, 2*TMath::Pi());
    t_pi0->Draw("phi2>>hPi0_phiPos", "helicity == +1 && " + cuts, "goff");
    TH1F* hPi0_phiNeg = new TH1F("hPi0_phiNeg", "hPi0_phiNeg", numBins, 0, 2*TMath::Pi());
    t_pi0->Draw("phi2>>hPi0_phiNeg", "helicity == -1 && " + cuts, "goff");

    TH1F *hAdd2 = new TH1F("hAdd2", "hAdd2", numBins, 0, 2*TMath::Pi());
    hAdd2->Add(hPi0_phiPos, hPi0_phiNeg);
    TH1F *hSub2 = new TH1F("hSub2", "hSub2", numBins, 0, 2*TMath::Pi());
    hSub2->Add(hPi0_phiPos, hPi0_phiNeg, 1, -1);
    TH1F *hPi0_Asym = new TH1F("hPi0_Asym", "hPi0_Asym", numBins, 0, 2*TMath::Pi());
    hPi0_Asym->Divide(hSub2, hAdd2);
    hPi0_Asym->Scale(1./beamPol);
    for (Int_t i=1; i<=numBins; i++) {
        hPi0_Asym->SetBinError(i, 1/TMath::Sqrt((hAdd2->GetBinContent(i))*beamPol));
    }
    TF1 *fit_pi0Asym = new TF1("fit_pi0Asym", "[0]*sin(x)/(1+[1]*cos(x))+[2]", 0, 2*TMath::Pi());
    fit_pi0Asym->SetParameter(0,1);
    fit_pi0Asym->SetParameter(1,1);
    fit_pi0Asym->SetParameter(2,0);
    //prac_ALU_pi0_asymmetry(hPi0_Asym, fit_pi0Asym);

    // models
    Int_t numBinDiv = 10;

    // vgg model
    TH1F *hDVCS_VGGAsym = new TH1F("hDVCS_VGGAsym", "hDVCS_VGGAsym", numBins*numBinDiv, 0, 2*TMath::Pi());
    std::ostringstream oss1;
    oss1 << "python -u prac_dvcsgens.py vgg_model " << vggFile << " " << xB << " " << q2 << " " << tpos << " " << numBins << " " << numBinDiv << " " << beamE;
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
        hDVCS_VGGAsym->SetBinContent(i, -1*(xsPosVal-xsNegVal)/(xsPosVal+xsNegVal));
    }
    if_vgg.close();

    // km15 model
    TH1F *hDVCS_KM15AsymK = new TH1F("hDVCS_KM15AsymK", "hDVCS_KM15AsymK", numBins*numBinDiv, 0, 2*TMath::Pi());
    std::ostringstream oss2;
    oss2 << "python -u prac_gepard.py km15_model " << km15File << " " << xB << " " << q2 << " " << tpos << " " << numBins << " " << numBinDiv << " " << beamE << " " << "ALU";
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
        hDVCS_KM15AsymK->SetBinContent(i, asyVal);
    }
    if_km15.close();

    // contamination
    TH1F *hPi0_phi = new TH1F("hPi0_phi", "hPi0_phi", numBins, 0, 2*TMath::Pi());
    t_pi0->Draw("phi2>>hPi0_phi", cuts, "goff");
    TH1F *hPi0Recon_phi = new TH1F("hPi0Recon_phi", "hPi0Recon_phi", numBins, 0, 2*TMath::Pi());
    t_pi0Recon->Draw("phi2>>hPi0Recon_phi", cuts, "goff");
    TH1F *hPi0DVCS_phi = new TH1F("hPi0DVCS_phi", "hPi0DVCS_phi", numBins, 0, 2*TMath::Pi());
    t_pi0DVCS->Draw("phi2>>hPi0DVCS_phi", cuts, "goff");    
    TH1F *hDVCS_phi = new TH1F("hDVCS_phi", "hDVCS_phi", numBins, 0, 2*TMath::Pi());
    t_dvcs->Draw("phi2>>hDVCS_phi", cuts, "goff");
    TH1F *hDiv1 = new TH1F("hDiv1", "hDiv1", numBins, 0, 2*TMath::Pi());
    hDiv1->Divide(hPi0DVCS_phi, hPi0Recon_phi);
    TH1F *hDiv2 = new TH1F("hDiv2", "hDiv2", numBins, 0, 2*TMath::Pi());
    hDiv2->Divide(hDiv1, hDVCS_phi);
    TH1F *hContamination = new TH1F("hContamination", "hContamination", numBins, 0, 2*TMath::Pi());
    hContamination->Multiply(hPi0_phi, hDiv2);
    for (int i=1; i<=numBins; i++) {
        if (hContamination->GetBinContent(i) == 0) {
            hContamination->SetBinContent(i, 0.0);
        }
    }

    // dvcs asymmetry
    TH1F *hDVCS_Asym = new TH1F("hDVCS_Asym", "hDVCS_Asym", numBins, 0, 2*TMath::Pi());
    for (int i=1; i<=numBins; i++) {
        hDVCS_Asym->SetBinContent(i, (1.0/(1.0 - hContamination->GetBinContent(i)))*(hDVCS_measAsym->GetBinContent(i) - (hContamination->GetBinContent(i))*(hPi0_Asym->GetBinContent(i))));
    }
    hDVCS_Asym->Scale(1./beamPol);
    TF1 *fit_dvcsAsym = new TF1("fit_dvcsAsym", "[0]*sin(x)/(1+[1]*cos(x))+[2]", 0, 2*TMath::Pi());
    fit_dvcsAsym->SetParameter(0,1);
    fit_dvcsAsym->SetParameter(1,1);
    fit_dvcsAsym->SetParameter(2,0);

    prac_ALU_dvcs_asymmetry(hDVCS_measAsym, hPi0_Asym, hDVCS_Asym, fit_measAsym, fit_pi0, fit_dvcsAsym);
}
