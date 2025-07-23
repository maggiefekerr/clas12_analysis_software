#include <TSystem.h>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

// this code is calculating ratios to show the difference in model prediction for asy and xs with all GPDs 
// contributing vs only H -> part of process for high level analysis
// data outputs to ratios directory
// python file used to generate data is binnedRatios_km15_asy_xs_LU.py

/*void binnedRatios_gpdCont() {
    // reading in meanVals file
    Double_t Q2mean[64], xBmean[64], tmean[64];
    std::ifstream if_meanVals("meanVals.txt");
    for (int i=0; i<64; i++) {
        Double_t a, b, c;
        if_meanVals >> a >> b >> c;
        Q2mean[i] = a;
        xBmean[i] = b;
        tmean[i] = c;
    }
    if_meanVals.close();

    // Broken points: 27, 39, 50, 52, 57, 58, 59, 60, 62

    // model
    Double_t amp[64];
    Double_t sig[64];
    TGraph *gKM15[64];
    TCanvas *c1[64];
    for (int i=0; i<64; i++) {
        c1[i] = new TCanvas(Form("c1_%i", i), Form("c1_%i", i), 1000, 750);
        gKM15[i] = new TGraph();
        std::ifstream if_KM15(Form("./km15_alu_output/continuous/fall2018_Q2_%.3f_xB_%.3f_t_%.3f.txt", Q2mean[i], xBmean[i], tmean[i]));
        Double_t phi, asy;
        // 27, 39, 50, 52, 57, 58, 59, 60, 62 - broken points
        if ((i == 26) || (i == 38) || (i == 49) || (i == 51) || (i == 56) || (i == 57) || (i == 58) || (i == 59) || (i == 61)) {
            c1[i]->Close();
            continue;
        }
        else {
            while (if_KM15 >> phi >> asy) {
                gKM15[i]->AddPoint(phi, asy);
            }
        }
        if_KM15.close();

        TF1 *fit1 = new TF1(Form("fit1_%i", i), "[0]*sin(x)/(1+[1]*cos(x))+[2]", 0, 2*TMath::Pi());
        fit1->SetParameter(0,1);
        fit1->SetParameter(1,1);
        fit1->SetParameter(2,0);

        gKM15[i]->Fit(Form("fit1_%i", i));
        amp[i] = fit1->GetParameter(0);
        sig[i] = fit1->GetParError(0);
        c1[i]->Close();
    }

    // model onlyH
    Double_t amp_onlyH[64];
    Double_t sig_onlyH[64];
    TGraph *gKM15_onlyH[64];
    TCanvas *c2[64];
    for (int i=0; i<64; i++) {
        c2[i] = new TCanvas(Form("c2_%i", i), Form("c2_%i", i), 1000, 750);
        gKM15_onlyH[i] = new TGraph();
        std::ifstream if_KM15_onlyH(Form("./km15_alu_output/onlyH_continuous/fall2018_Q2_%.3f_xB_%.3f_t_%.3f.txt", Q2mean[i], xBmean[i], tmean[i]));
        Double_t phi, asy;
        // 27, 39, 50, 52, 57, 58, 59, 60, 62 - broken points
        if ((i == 26) || (i == 38) || (i == 49) || (i == 51) || (i == 56) || (i == 57) || (i == 58) || (i == 59) || (i == 61)) {
            c2[i]->Close();
            continue;
        }
        else {
            while (if_KM15_onlyH >> phi >> asy) {
                gKM15_onlyH[i]->AddPoint(phi, asy);
            }
        }
        if_KM15_onlyH.close();

        TF1 *fit2 = new TF1(Form("fit2_%i", i), "[0]*sin(x)/(1+[1]*cos(x))+[2]", 0, 2*TMath::Pi());
        fit2->SetParameter(0,1);
        fit2->SetParameter(1,1);
        fit2->SetParameter(2,0);

        gKM15_onlyH[i]->Fit(Form("fit2_%i", i));
        amp_onlyH[i] = fit2->GetParameter(0);
        sig_onlyH[i] = fit2->GetParError(0);
        c2[i]->Close();
    }

    // plot ratios
    TGraph *g_inQ2 = new TGraph();
    TGraph *g_inxB = new TGraph();
    TGraph *g_int = new TGraph();

    for (int i=0; i<64; i++) {
        // 27, 39, 50, 52, 57, 58, 59, 60, 62 - broken points
        if ((i == 26) || (i == 38) || (i == 49) || (i == 51) || (i == 56) || (i == 57) || (i == 58) || (i == 59) || (i == 61)) {
            continue;
        }
        else {
            g_inQ2->AddPoint(Q2mean[i], amp_onlyH[i]/amp[i]);
            g_inxB->AddPoint(xBmean[i], amp_onlyH[i]/amp[i]);
            g_int->AddPoint(-1*tmean[i], amp_onlyH[i]/amp[i]);

            if (amp_onlyH[i]/amp[i] > 1.5) {
                cout << "Ratio " << amp_onlyH[i]/amp[i] << ": Q2 = " << Q2mean[i] << " GeV^2, xB = " << xBmean[i] << ", -t = " << -1*tmean[i] << " GeV^2" << endl; 
            }
        }
    }

    TCanvas *c3 = new TCanvas("c3", "c3", 4000, 3500);
    g_inQ2->SetTitle("Ratio of A_{LU} from Only #it{H} to Full Model");
    g_inQ2->GetXaxis()->SetTitle("Q^2 [GeV^2]");
    g_inQ2->GetYaxis()->SetTitle("A_{LU}^{#it{H}} / A_{LU}^{all}");
    g_inQ2->SetMarkerStyle(21);
    g_inQ2->SetMarkerSize(1);
    g_inQ2->GetYaxis()->SetRangeUser(0.0, 5.0);
    g_inQ2->Draw("AP");
    c3->SaveAs("ratioHtoAll_inQ2.png");
    c3->Close();

    TCanvas *c4 = new TCanvas("c4", "c4", 4000, 3500);
    g_inxB->SetTitle("Ratio of A_{LU} from Only #it{H} to Full Model");
    g_inxB->GetXaxis()->SetTitle("x_{B}");
    g_inxB->GetYaxis()->SetTitle("A_{LU}^{#it{H}} / A_{LU}^{all}");
    g_inxB->SetMarkerStyle(21);
    g_inxB->SetMarkerSize(1);
    g_inxB->GetYaxis()->SetRangeUser(0.0, 5.0);
    g_inxB->Draw("AP");
    c4->SaveAs("ratioHtoAll_inxB.png");
    c4->Close();

    TCanvas *c5 = new TCanvas("c5", "c5", 4000, 3500);
    g_int->SetTitle("Ratio of A_{LU} from Only #it{H} to Full Model");
    g_int->GetXaxis()->SetTitle("-t [GeV^2]");
    g_int->GetYaxis()->SetTitle("A_{LU}^{#it{H}} / A_{LU}^{all}");
    g_int->SetMarkerStyle(21);
    g_int->SetMarkerSize(1);
    g_int->GetYaxis()->SetRangeUser(0.0, 5.0);
    g_int->Draw("AP");
    c5->SaveAs("ratioHtoAll_int.png");
    c5->Close();
}*/

//sqrt(1 - y - (y**2*eps2)/4)
// K = sqrt(K2(Q2, xB, t, y, eps2))

bool kinematic_check(Double_t Q2, Double_t xB, Double_t t) {
    // setting needed kinematics
    Double_t Mp = 0.938272;
    Double_t eps2 = TMath::Power(2.*xB*Mp/TMath::Sqrt(Q2), 2.0);
    Double_t W = TMath::Sqrt((Q2/xB) - Q2 + (Mp*Mp));
    Double_t s = 2.*Mp*10.6 + (Mp*Mp);
    Double_t y = ((W*W) + Q2 - (Mp*Mp))/(s - (Mp*Mp));
    Double_t tmin = -1*Q2*(2.*(1. - xB)*(1. - TMath::Sqrt(1.+eps2)) + eps2)/(4.*xB*(1.-xB) + eps2);
    
    // calculate value of 1 - y - (y**2*eps2)/4 -> should be non-negative as it is in a sqrt
    Double_t test1 = 1 - y - (y*y*eps2)/4.;
    Double_t test2 = (TMath::Sqrt(1.+eps2) + (4.*xB *(1.-xB)+eps2)/(4.*(1.-xB)) * (t - tmin)/Q2) *
                    (-1*(t/Q2) * (1.-xB) * (1.-y-(y*y*eps2)/4.) * (1. - (tmin/t)));

    cout << Q2 << " " << xB << " " << t << endl;            
    if (test1 < 0) {
        cout << "test 1 negative" << endl;
    }
    if  (test2 < 0) {
        cout << "test 2 negative" << endl;
    }
    if ((test1 > 0) && (test2 > 0)) {
        cout << "all clear" << endl;
        return true;
    }
    else {
        return false;
    }
}

void make_model_inputs() {
    Double_t t[4] = {0.1, 0.4, 0.7, 1.0};
    for (int k=0; k<4; k++) {
        std::ofstream of_inputs(Form("./ratios/data/input_kinematics_-t_%.2f.txt", t[k]));
        for (int i=0; i<61; i++) {
            Double_t Q2 = 1.0 + 0.1*i;

            for (int j=0; j<65; j++) {
                Double_t xB = 0.06 + 0.01*j;

                if (kinematic_check(Q2, xB, -1*t[k])) {
                    of_inputs << Q2 << " " << xB << endl;
                }
            }     
        }
        of_inputs.close();
    }
}

void get_model_data_asy(TString directory, Int_t tVal){ // current options are "all" and "only_h"
    Double_t t[4] = {0.1, 0.4, 0.7, 1.0};
    std::ifstream if_inputs(Form("./ratios/data/input_kinematics_-t_%.2f.txt", t[tVal]));
    Double_t Q2, xB;
    while (if_inputs >> Q2 >> xB) {
        TString fileName = Form("alu_Q2-%.2f_xB-%.2f_-t-%.2f", Q2, xB, t[tVal]);
        std::ostringstream oss;
        oss << "python -u binnedRatios_km15_asy_LU.py asy " << directory << " " << fileName << " " << Q2 << " " << xB << " " << t[tVal] << " 10.6 100";
        std::string run = oss.str();
        FILE* pipe = gSystem->OpenPipe(run.c_str(), "r");
        if (!pipe) {
            std::cerr << "Failed to run script" << std::endl;
            return;
        }
    }
}

void get_model_data_xs(TString directory, Int_t tVal){ // current options are "all" and "only_h"
    Double_t t[4] = {0.1, 0.4, 0.7, 1.0};
    std::ifstream if_inputs(Form("./ratios/data/input_kinematics_-t_%.2f.txt", t[tVal]));
    Double_t Q2, xB;
    while (if_inputs >> Q2 >> xB) {
        TString fileName = Form("xslu_Q2-%.2f_xB-%.2f_-t-%.2f", Q2, xB, t[tVal]);
        std::ostringstream oss;
        oss << "python -u binnedRatios_km15_xs_LU.py xs " << directory << " " << fileName << " " << Q2 << " " << xB << " " << t[tVal] << " 10.6 100";
        std::string run = oss.str();
        FILE* pipe = gSystem->OpenPipe(run.c_str(), "r");
        if (!pipe) {
            std::cerr << "Failed to run script" << std::endl;
            return;
        }
    }
}

// plots full, only H; fits; writes out fit to file
void plot_model_data_asy() {
    Double_t t[4] = {0.1, 0.4, 0.7, 1.0};

    for (int i=0; i<4; i++) {
        std::ifstream if_inputs(Form("./ratios/data/input_kinematics_-t_%.2f.txt", t[i]));
        std::ofstream of_outputs(Form("./ratios/data/alu_output_kinematics_&_amps_-t_%.2f.txt", t[i]));
        Double_t Q2, xB;

        while (if_inputs >> Q2 >> xB) {
            // dealing with all GPDs
            std::ifstream if_all(Form("./ratios/data/all/asy/-t_%.1f/alu_Q2-%.2f_xB-%.2f_-t-%.2f.txt", t[i], Q2, xB, t[i]));
            TGraph *gAll = new TGraph();
            Double_t phiA, asyA;
            while (if_all >> phiA >> asyA) {
                gAll->AddPoint(phiA, asyA);
            }
            TF1 *fitA = new TF1("fitA", "[0]*sin(x)/(1+[1]*cos(x))+[2]", 0, 2*TMath::Pi());
            fitA->SetParameter(0,1);
            fitA->SetParameter(1,1);
            fitA->SetParameter(2,0);
            gAll->Fit("fitA");
            if_all.close();

            // dealing with only H
            std::ifstream if_onlyH(Form("./ratios/data/only_h/asy/-t_%.1f/alu_Q2-%.2f_xB-%.2f_-t-%.2f.txt", t[i], Q2, xB, t[i]));
            TGraph *gOnlyH = new TGraph();
            Double_t phiH, asyH;
            while (if_onlyH >> phiH >> asyH) {
                gOnlyH->AddPoint(phiH, asyH);
            }
            TF1 *fitH = new TF1("fitH", "[0]*sin(x)/(1+[1]*cos(x))+[2]", 0, 2*TMath::Pi());
            fitH->SetParameter(0,1);
            fitH->SetParameter(1,1);
            fitH->SetParameter(2,0);
            gOnlyH->Fit("fitH");
            if_onlyH.close();

            Double_t ampA = fitA->GetParameter(0);
            Double_t ampH = fitH->GetParameter(0);

            of_outputs << Q2 << " " << xB << " " << " " << ampA << " " << ampH << " " << ampH/ampA << endl;
        }

        if_inputs.close();
        of_outputs.close();
    }
}

// plots full, only H; fits; writes out fit to file
void plot_model_data_xs(){
    Double_t t[4] = {0.1, 0.4, 0.7, 1.0};

    for (int i=0; i<4; i++) {
        std::ifstream if_inputs(Form("./ratios/data/input_kinematics_-t_%.2f.txt", t[i]));
        std::ofstream of_outputs(Form("./ratios/data/xslu_output_kinematics_&_amps_-t_%.2f.txt", t[i]));
        Double_t Q2, xB;

        while (if_inputs >> Q2 >> xB) {
            // dealing with all GPDs
            std::ifstream if_all(Form("./ratios/data/all/xs/-t_%.1f/xslu_Q2-%.2f_xB-%.2f_-t-%.2f.txt", t[i], Q2, xB, t[i]));
            TGraph *gAll = new TGraph();
            Double_t phiA, xsA;
            while (if_all >> phiA >> xsA) {
                gAll->AddPoint(phiA, xsA);
            }
            TF1 *fitA = new TF1("fitA", "[0]+[1]*cos(x)+[2]*cos(2*x)", 0, 2*TMath::Pi());
            fitA->SetParameter(0,0);
            fitA->SetParameter(1,-0.5);
            fitA->SetParameter(2,0);
            gAll->Fit("fitA");
            if_all.close();

            // dealing with only H
            std::ifstream if_onlyH(Form("./ratios/data/only_h/xs/-t_%.1f/xslu_Q2-%.2f_xB-%.2f_-t-%.2f.txt", t[i], Q2, xB, t[i]));
            TGraph *gOnlyH = new TGraph();
            Double_t phiH, xsH;
            while (if_onlyH >> phiH >> xsH) {
                gOnlyH->AddPoint(phiH, xsH);
            }
            TF1 *fitH = new TF1("fitH", "[0]+[1]*cos(x)+[2]*cos(2*x)", 0, 2*TMath::Pi());
            fitH->SetParameter(0,0);
            fitH->SetParameter(1,-0.5);
            fitH->SetParameter(2,0);
            gOnlyH->Fit("fitH");
            if_onlyH.close();

            Double_t intA = fitA->Integral(0, 2*(TMath::Pi()));
            Double_t intH = fitH->Integral(0, 2*(TMath::Pi()));

            of_outputs << Q2 << " " << xB << " " << " " << intA << " " << intH << " " << intH/intA << endl;
        }

        if_inputs.close();
        of_outputs.close();
    }
}

void plot_ratios_asy() {
    Double_t t[4] = {0.1, 0.4, 0.7, 1.0};
    TCanvas *c1[4];
    TH2F *hRatio[4];

    for (int i=0; i<4; i++) {
        c1[i] = new TCanvas(Form("c1_%i", i), Form("c1_%i", i), 2000, 1500);
        hRatio[i] = new TH2F(Form("hRatio_%i", i), "", 65, 0.055, 0.705, 61, 0.95, 7.05);
        
        std::ifstream if_asyms(Form("./ratios/data/alu_output_kinematics_&_amps_-t_%.2f.txt", t[i]));
        Double_t Q2, xB, asyA, asyH, ratio;

        while (if_asyms >> Q2 >> xB >> asyA >> asyH >> ratio) {
            hRatio[i]->SetBinContent(hRatio[i]->GetXaxis()->FindBin(xB), hRatio[i]->GetYaxis()->FindBin(Q2), ratio);
        }

        hRatio[i]->SetTitle(Form("A^{#it{H}}_{LU} / A^{all}_{LU}, -t = %.2f GeV^{2}", t[i]));
        hRatio[i]->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");
        hRatio[i]->GetXaxis()->SetTitle("x_{B}");
        hRatio[i]->GetZaxis()->SetTitle("A^{#it{H}}_{LU} / A^{all}_{LU}");
        hRatio[i]->GetZaxis()->SetRangeUser(0.0, 5.0);
        hRatio[i]->SetStats(0);

        hRatio[i]->Draw("COLZ");

        c1[i]->SaveAs(Form("./ratios/ratio_asym_-t_%.2f.png", t[i]));
        c1[i]->Close();
        if_asyms.close();
    }
}

void plot_ratios_xs() {
    Double_t t[4] = {0.1, 0.4, 0.7, 1.0};
    TCanvas *c1[4];
    TH2F *hRatio[4];

    for (int i=0; i<4; i++) {
        c1[i] = new TCanvas(Form("c1_%i", i), Form("c1_%i", i), 2000, 1500);
        hRatio[i] = new TH2F(Form("hRatio_%i", i), "", 65, 0.055, 0.705, 61, 0.95, 7.05);
        
        std::ifstream if_xs(Form("./ratios/data/xslu_output_kinematics_&_amps_-t_%.2f.txt", t[i]));
        Double_t Q2, xB, intA, intH, ratio;

        while (if_xs >> Q2 >> xB >> intA >> intH >> ratio) {
            hRatio[i]->SetBinContent(hRatio[i]->GetXaxis()->FindBin(xB), hRatio[i]->GetYaxis()->FindBin(Q2), ratio);
        }

        hRatio[i]->SetTitle(Form("#int_{0}^{2#pi}d#sigma^{#it{H}} / #int_{0}^{2#pi}d#sigma^{all}, -t = %.2f GeV^{2}", t[i]));
        hRatio[i]->GetYaxis()->SetTitle("Q^2 [GeV^{2}]");
        hRatio[i]->GetXaxis()->SetTitle("x_{B}");
        //hRatio[i]->GetZaxis()->SetTitle("#int_{0}^{2#pi}d#sigma^{#it{H}} / #int_{0}^{2#pi}d#sigma^{all}");
        hRatio[i]->GetZaxis()->SetRangeUser(0.0, 2.0);
        hRatio[i]->SetStats(0);
        hRatio[i]->Draw("COLZ");

        c1[i]->SaveAs(Form("./ratios/ratio_xs_-t_%.2f.png", t[i]));
        c1[i]->Close();
        if_xs.close();
    }
}