#include <TSystem.h>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

void vgg_aut(TString paramFile, Int_t numBins, Double_t beamE=10.604){
    std::ifstream if_input_vgg(paramFile);
    for (int i=0; i<4; i++) {
        TCanvas *c1 = new TCanvas("c1", "c1", 6000, 4500);
        c1->Divide(3,3,0,0);
        for (int j=0; j<7; j++) {
            c1->cd(j+3);
            TH1F *hVGG = new TH1F("hVGG", "", numBins, 0, 360);
            Double_t tpos, xB, Q2;
            if_input_vgg >> tpos >> xB >> Q2;

            TString histFile = Form("./vgg_aut_output/output_-t_%.3f_xB_%.3f_Q2_%.3f.txt", tpos, xB, Q2);
            std::ostringstream oss;
            oss << "python -u vgg_aut.py vgg_aut " << histFile << " " << tpos << " " << xB << " " << Q2 << " " << numBins << " " << beamE;
            std::string run = oss.str();
            FILE* pipe = gSystem->OpenPipe(run.c_str(), "r");
            if (!pipe) {
                std::cerr << "Failed to run script" << std::endl;
                return;
            }

            ifstream if_output_vgg(histFile, ifstream::in);
            Double_t phi, asy;
            for(int k=1; k <=numBins; k++) {
                if_output_vgg >> phi >> asy;
                hVGG->SetBinContent(k, asy);
            }
            if_output_vgg.close();

            hVGG->GetYaxis()->SetRangeUser(-2.0,1.0);
            hVGG->SetMarkerColor(1);
            hVGG->SetMarkerStyle(21);
            hKMhVGG15->SetMarkerSize(1);
            hVGG->SetStats(0);
            hVGG->Draw("P");

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
                hVGG->GetYaxis()->SetTitle("A_{UT}");
            }
            if ((j == 4) || (j == 5) || (j == 6)) {
                hVGG->GetXaxis()->SetTitle("#phi [#circ]");
            }
            c1->Update();
        }
        c1->Update();
        c1->SaveAs(Form("./vgg_aut_output/output_%d.png", i));   
    }
    if_input_vgg.close();
}