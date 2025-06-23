#include <TSystem.h>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

void km15_aut(TString km15File, Int_t numBins, Double_t beamE=10.604){
    // currently there are 28 binnings, so just hard-coding that to make it easier
    std::ifstream if_input_km15(km15File);
    for (int i=0; i<4; i++) {
        TCanvas *c1 = new TCanvas("c1", "c1", 6000, 4500);
        c1->Divide(3,3,0,0);
        for (int j=0; j<7; j++) {
            c1->cd(j+3);
            TH1F *hKM15 = new TH1F("hKM15", "", numBins, 0, 360);
            Double_t tpos, xB, Q2;
            if_input_km15 >> tpos >> xB >> Q2;

            TString histFile = Form("./output_-t_%.3f_xB_%.3f_Q2_%.3f.txt", tpos, xB, Q2);
            std::ostringstream oss;
            oss << "python -u km15_aut.py km15_aut " << histFile << " " << tpos << " " << xB << " " << Q2 << " " << numBins << " " << beamE;
            std::string run = oss.str();
            FILE* pipe = gSystem->OpenPipe(run.c_str(), "r");
            if (!pipe) {
                std::cerr << "Failed to run script" << std::endl;
                return;
            }

            ifstream if_output_km15(histFile, ifstream::in);
            Double_t phi, asy;
            for(int k=1; k <=numBins; k++) {
                if_output_km15 >> phi >> asy;
                hKM15->SetBinContent(k, asy);
            }
            if_output_km15.close();

            hKM15->GetYaxis()->SetRangeUser(-2.0,1.0);
            hKM15->SetMarkerColor(1);
            hKM15->SetMarkerStyle(21);
            hKM15->SetMarkerSize(1);
            hKM15->SetStats(0);
            hKM15->Draw("P");

            Double_t textSize[2] = {0.05, 0.06};
            Double_t xPos[2] = {0.15, 0.10};
            TLatex *l = new TLatex();
            if ((j == 4) || (j == 5) || (j == 6)) {
                l->SetTextSize(0.0575);
            }
            else {
                l->SetTextSize(0.065);
            }
            if (j == 4) {
                l->DrawLatexNDC(0.15, 0.300, Form("-t = %.3f GeV^{2}", tpos));
                l->DrawLatexNDC(0.15, 0.225, Form("x_{B} = %.3f", xB));
                l->DrawLatexNDC(0.15, 0.150, Form("Q^{2} = %.3f GeV^{2}", Q2));
            }
            else if (j == 1) {
                l->DrawLatexNDC(0.15, 0.200, Form("-t = %.3f GeV^{2}", tpos));
                l->DrawLatexNDC(0.15, 0.125, Form("x_{B} = %.3f", xB));
                l->DrawLatexNDC(0.15, 0.050, Form("Q^{2} = %.3f GeV^{2}", Q2));
            }
            else if ((j == 5) || (j == 6)){
                l->DrawLatexNDC(0.05, 0.300, Form("-t = %.3f GeV^{2}", tpos));
                l->DrawLatexNDC(0.05, 0.225, Form("x_{B} = %.3f", xB));
                l->DrawLatexNDC(0.05, 0.150, Form("Q^{2} = %.3f GeV^{2}", Q2));
            }
            else {
                l->DrawLatexNDC(0.05, 0.200, Form("-t = %.3f GeV^{2}", tpos));
                l->DrawLatexNDC(0.05, 0.125, Form("x_{B} = %.3f", xB));
                l->DrawLatexNDC(0.05, 0.050, Form("Q^{2} = %.3f GeV^{2}", Q2));
            }

            if ((j == 1) || (j == 4)) {
                hKM15->GetYaxis()->SetTitle("TTSA");
            }
            if ((j == 4) || (j == 5) || (j == 6)) {
                hKM15->GetXaxis()->SetTitle("#phi [#circ]");
            }
            c1->Update();
        }
        c1->Update();
        c1->SaveAs(Form("./km15_aut_output/output_%d.png", i));   
    }
    if_input_km15.close();
}