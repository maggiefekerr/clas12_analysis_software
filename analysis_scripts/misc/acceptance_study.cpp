#include <TMath.h>
#include <TRandom3.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TAxis.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <utility> // For std::pair
#include <TH1D.h>     // For histograms
#include <TLegend.h>  // For legends
#include <TMinuit.h>
#include <Math/Integrator.h>
#include <Math/IntegratorMultiDim.h>
#include <functional>

// Global structure to hold deviations. 
// Each pair corresponds to an exclusion case, holding two vectors for deviations in B and C respectively.
std::vector<std::pair<std::vector<double>, std::vector<double>>> deviationsForCases(6);
std::vector<std::pair<std::vector<double>, std::vector<double>>> deviationsForCasesMLM(6); 
std::vector<double> phiVecGlobal;
int binsToExcludeGlobal;

// Fit function: a trigonometric polynomial
double fitFunction(double x, double *par) {
    return par[0] * (1 + par[1] * cos(x) + par[2] * cos(2*x));
}

#include <cmath> // For sin, cos, log

void negLogLikelihood(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
    double B = par[0];
    double C = par[1];

    double a = binsToExcludeGlobal * 2 * M_PI / 24;
    double b = (24 - binsToExcludeGlobal) * 2 * M_PI / 24;

    // Calculate normalization factor alpha
    double integral = (b - a) + B * (sin(b) - sin(a)) + C / 2 * (sin(2*b) - sin(2*a));
    double alpha = 1 / integral;

    double sum = 0;
    for (size_t phi = 0; phi < phiVecGlobal.size(); ++phi) {
        if (phiVecGlobal[phi] >= a && phiVecGlobal[phi] < b) {
            double probDensity = alpha * (1 + B * cos(phiVecGlobal[phi]) + C * cos(2*phiVecGlobal[phi]));
            if (probDensity > 0) {
                sum += log(probDensity);
            } 
        }
    }

    f = -sum;
}


void generateData(double B, double C, long unsigned int N, std::vector<double>& phiVec) {
    TRandom3 rand(0); // Seed for reproducibility
    double A = 1.0; // Assuming A is constant
    long unsigned int accepted = 0;
    while (accepted < N) {
        double phi = rand.Uniform(0, 2*TMath::Pi());
        // double phi = rand.Uniform(0 .5, 2*TMath::Pi()-0.5);
        double value = A * (1 + B * TMath::Cos(phi) + C * TMath::Cos(2 * phi));
        double randVal = rand.Uniform(0, A * (1 + fabs(B) + fabs(C))); // Max possible value

        if (randVal < value) {
            phiVec.push_back(phi);
            accepted++;
        }
    }
}

void setupCanvas(TCanvas* masterCanvas) {
    masterCanvas->Divide(3, 2, 0, 0); // Set gaps to 0

    for (int i = 1; i <= 6; ++i) {
        masterCanvas->cd(i);
        TPad *pad = (TPad*)gPad;
        
        pad->SetTopMargin((i <= 3) ? 0.02 : 0.00); // Smaller margin for top row
        pad->SetBottomMargin((i > 3) ? 0.15 : 0.00); // Larger bottom margin only for bottom row
        pad->SetLeftMargin(((i-1) % 3 == 0) ? 0.15 : 0.00); // Larger left margin only for first column
        pad->SetRightMargin(((i-1) % 3 == 2) ? 0.01 : 0.00); // Uniform right margin
    }
}


void plotForExclusion(const std::vector<double>& phiVec, double B, double C, int canvasIndex, int binsToExclude, TCanvas* masterCanvas) {
    double binWidth = (2 * TMath::Pi()) / 24;
    TGraphErrors *graphIncluded = new TGraphErrors();
    TGraphErrors *graphExcluded = new TGraphErrors();
    graphExcluded->SetMarkerStyle(4); // Open circles for excluded points

    // Initialize bin counts
    std::vector<int> binCountsIncluded(24, 0), binCountsExcluded(24, 0);

    for (double phi : phiVec) {
        int binIndex = static_cast<int>(phi / binWidth);
        if (binIndex >= binsToExclude && binIndex < 24 - binsToExclude) {
            binCountsIncluded[binIndex]++;
        } else {
            binCountsExcluded[binIndex]++;
        }
    }

    // Set points for included and excluded graphs
    for (int i = 0; i < 24; ++i) {
        double binCenter = binWidth * i + binWidth / 2;
        if (binCountsIncluded[i] > 0) {
            graphIncluded->SetPoint(graphIncluded->GetN(), binCenter, binCountsIncluded[i]);
            graphIncluded->SetPointError(graphIncluded->GetN() - 1, 0, TMath::Sqrt(binCountsIncluded[i]));
        }
        if (binCountsExcluded[i] > 0) {
            graphExcluded->SetPoint(graphExcluded->GetN(), binCenter, binCountsExcluded[i]);
            graphExcluded->SetPointError(graphExcluded->GetN() - 1, 0, TMath::Sqrt(binCountsExcluded[i]));
        }
    }

    double exclusionPercentage = (binsToExclude / 24.0) * 100.0; // Update formula to accurately represent the exclusion percentage

    masterCanvas->cd(canvasIndex);
    graphIncluded->SetMarkerStyle(20); // Filled circles for included points
    graphIncluded->Draw("AP");
    graphExcluded->Draw("P SAME");
    graphExcluded->GetXaxis()->SetLimits(0, TMath::TwoPi());
    graphExcluded->GetYaxis()->SetLimits(13000, 17000);

    // Find the maximum value of the distribution for the initial guess of [0]
	double maxY = 0;
	for (int i = 0; i < graphIncluded->GetN(); ++i) {
	    double x, y;
	    graphIncluded->GetPoint(i, x, y);
	    if (y > maxY) maxY = y;
	}

    /****** CHI2 MINIMIZATION PORTION *******/
	// Perform the fit within the limited range
	TF1 *fitFuncLimited = new TF1("fitFuncLimited", "[0]*(1 + [1]*cos(x) + [2]*cos(2*x))", binWidth * binsToExclude, 2*TMath::Pi() - binWidth * binsToExclude);
	// Set initial parameters
	fitFuncLimited->SetParameter(0, maxY); // Initial guess for [0] is the max Y value
	fitFuncLimited->SetParameter(1, B);    // Initial guess for [1]
	fitFuncLimited->SetParameter(2, C);    // Initial guess for [2]
	// Set parameter limits for [1] and [2]
	fitFuncLimited->SetParLimits(1, -1, 1); // Constrain [1] to be between -1 and 1
	fitFuncLimited->SetParLimits(2, -1, 1); // Constrain [2] to be between -1 and 1
	// Fit the graph with these settings
	graphIncluded->Fit(fitFuncLimited, "Q R");
	// Extract fitted parameters and their errors
    double fittedB = fitFuncLimited->GetParameter(1);
    double fittedC = fitFuncLimited->GetParameter(2);
    double errB = fitFuncLimited->GetParError(1);
    double errC = fitFuncLimited->GetParError(2);
    // Calculate deviations in sigma
    double deviationSigmaB = (fittedB - B) / errB;
    double deviationSigmaC = (fittedC - C) / errC;
    // Store deviations
    deviationsForCases[canvasIndex-1].first.push_back(deviationSigmaB); // Store deviation for B
    deviationsForCases[canvasIndex-1].second.push_back(deviationSigmaC); // Store deviation for C
	// Define a new TF1 that uses the parameters from the fit but extends across the full range
	TF1 *fitFuncFullRange = new TF1("fitFuncFullRange", "[0]*(1 + [1]*cos(x) + [2]*cos(2*x))", 0, 2*TMath::Pi());
	fitFuncFullRange->SetParameters(fitFuncLimited->GetParameters());
	fitFuncFullRange->SetLineColor(kRed);

    /****** MLM MINIMIZATION PORTION *******/
    double arglist[10]; arglist[0] = 50000; arglist[1] = 0.01;
    int ierflg = 0;
    TMinuit minuit(2);
    minuit.SetPrintLevel(-1);
    minuit.SetErrorDef(0.5);
    minuit.SetFCN(negLogLikelihood);
    // minuit.DefineParameter(0, "A", maxY, 0.001, 0, 0);
    minuit.DefineParameter(0, "B", B, 0.0001, 0, 0);
    minuit.DefineParameter(1, "C", C, 0.0001, 0, 0);
    minuit.mnexcm("MIGRAD", arglist, 2, ierflg);

    double fittedA, errA; minuit.GetParameter(0,fittedA,errA);
    // std::cout << fittedB << " " << fittedC << std::endl;
    minuit.GetParameter(0, fittedB, errB);
    minuit.GetParameter(1, fittedC, errC);
    // std::cout << fittedB << " " << fittedC << std::endl << std::endl;
    // Calculate deviations in sigma
    deviationSigmaB = (fittedB - B) / errB;
    deviationSigmaC = (fittedC - C) / errC;
    // Store deviations
    deviationsForCasesMLM[canvasIndex-1].first.push_back(deviationSigmaB);
    deviationsForCasesMLM[canvasIndex-1].second.push_back(deviationSigmaC); 

    // Assuming 'fittedA' holds the A value from the chi2 fit,
    // and 'fittedB', 'fittedC' are from the MLM fit:
    TF1* fitFuncMLMChi2 = new TF1("fitFuncMLMChi2", "[0]*(1 + [1]*cos(x) + [2]*cos(2*x))", 0, 2*TMath::Pi());

    // Set the parameters with A from chi2 fit and B, C from MLM fit
    fitFuncMLMChi2->SetParameter(0, fitFuncLimited->GetParameter(0)); // Use A from chi2 fit
    fitFuncMLMChi2->SetParameter(1, fittedB); // Use B from MLM fit
    fitFuncMLMChi2->SetParameter(2, fittedC); // Use C from MLM fit
    fitFuncMLMChi2->SetLineColor(kBlue); // Set a different color for distinction

    /****** PLOTTING PORTION *******/
	masterCanvas->cd(canvasIndex);
	graphIncluded->Draw("AP");
	graphIncluded->GetYaxis()->SetRangeUser(0.4*fitFuncLimited->GetParameter(0), 1.6 *fitFuncLimited->GetParameter(0));
	graphExcluded->Draw("P SAME");
	fitFuncFullRange->Draw("SAME");
    fitFuncMLMChi2->Draw("SAME");

    // Customize axis limits and labels
    graphIncluded->GetXaxis()->SetTitle("#phi");
    graphIncluded->GetYaxis()->SetTitle("Counts");
    graphIncluded->GetXaxis()->SetTitleSize(0.05);
    graphIncluded->GetYaxis()->SetTitleSize(0.05);
    if (canvasIndex<5) {
    	graphIncluded->GetXaxis()->SetLimits(0, TMath::TwoPi());
    } else {
    	graphIncluded->GetXaxis()->SetLimits(0.001, TMath::TwoPi());
    }

    // Calculate the start position for TPaveText based on the column, directly here
    double textStartX = (canvasIndex % 3 == 1) ? 0.15 : 0.0; // Example: 0.15 for left column, adjust 0.12 for middle/right columns
    TPaveText *pt = new TPaveText(textStartX, 0.65, textStartX + 0.575, 1.0, "NDC");
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->SetTextSize(0.04); // Adjust the text size if needed
    pt->AddText(Form("Exclusion: %.1f%%", exclusionPercentage));
    // For chi^2 lines in red
    TText *chi2Text1 = pt->AddText(Form("#chi^{2}: A_{UU}^{cos#phi} = %.3f #pm %.3f", fitFuncLimited->GetParameter(1),fitFuncLimited->GetParError(1)));
    chi2Text1->SetTextColor(kRed);
    TText *chi2Text2 = pt->AddText(Form("#chi^{2}: A_{UU}^{cos2#phi} = %.3f #pm %.3f", fitFuncLimited->GetParameter(2),fitFuncLimited->GetParError(2)));
    chi2Text2->SetTextColor(kRed);
    double chi2 = fitFuncLimited->GetChisquare();
    double ndf = fitFuncLimited->GetNDF();
    TText *chi2NdfText = pt->AddText(Form("#chi^{2}/ndf = %.3f", chi2 / ndf));
    chi2NdfText->SetTextColor(kRed);
    // For MLM lines in blue
    TText *mlmText1 = pt->AddText(Form("MLM: A_{UU}^{cos#phi} = %.3f #pm %.3f", fittedB, errB));
    mlmText1->SetTextColor(kBlue);
    TText *mlmText2 = pt->AddText(Form("MLM: A_{UU}^{cos2#phi} = %.3f #pm %.3f", fittedC, errC));
    mlmText2->SetTextColor(kBlue);
    pt->Draw();

    // After drawing the graphs and fit function, adjust axis titles visibility if necessary
    bool isLeftColumn = (canvasIndex % 3 == 1); // Adjust based on your actual layout
    bool isBottomRow = (canvasIndex > 3); // Adjust based on your actual layout

    if (!isLeftColumn) {
        graphIncluded->GetYaxis()->SetTitle("");
    }
    if (!isBottomRow) {
        graphIncluded->GetXaxis()->SetTitle("");
    }
}

void acceptanceStudy(double B, double C, int iterations) {
    for (int loop = 0; loop < iterations; ++loop) {
        std::cout << "Starting loop " << loop + 1 << " of " << iterations << ". " << std::endl;
        phiVecGlobal.clear();

        std::vector<double> phiVec;
        generateData(B, C, 1e5, phiVec); // Generate fresh data for each iteration
        phiVecGlobal = phiVec;

        // Create the example plot only during the first iteration
        TCanvas *masterCanvas = nullptr;
        masterCanvas = new TCanvas("masterCanvas", "Phi Distribution Fits", 1200, 800);
        setupCanvas(masterCanvas);

        int exclusionSteps[] = {0, 1, 2, 3, 4, 5};
        for (int i = 0; i < sizeof(exclusionSteps)/sizeof(exclusionSteps[0]); ++i) {
            int binsToExclude = exclusionSteps[i];
            binsToExcludeGlobal = binsToExclude;
            plotForExclusion(phiVec, B, C, i + 1, binsToExclude, masterCanvas);
        }

        if (loop == 0) {
            std::ostringstream filename;
            filename << "output/example_B=" << B << "_C=" << C << ".png";
            masterCanvas->SaveAs(filename.str().c_str());
        }
        delete masterCanvas;
    }
}

void plotDeviationsDistributionsB(double B, double C) {
    TCanvas* canvas = new TCanvas("canvas", "Deviations Distribution", 1200, 800);
    canvas->Divide(3, 2, 0.0, 0.0); // No gaps between pads

    int exclusionSteps[] = {0, 1, 2, 3, 4, 5}; 
    for (int i = 0; i < 6; ++i) {
        canvas->cd(i + 1);
        TPad *pad = (TPad*)gPad;

        // Set margins to remove space between subplots. Adjust values as needed.
        float leftMargin = (i % 3 == 0) ? 0.15 : 0.0;
        float rightMargin = 0.0;
        float topMargin = (i < 3) ? 0.05 : 0.0;
        float bottomMargin = (i >= 3) ? 0.15 : 0.0;
        pad->SetMargin(leftMargin, rightMargin, bottomMargin, topMargin);

        // The rest of the code to draw histograms and legends goes here
        double exclusionPercentage = (exclusionSteps[i] / 24.0) * 100.0;
        TH1D* hischi2 = new TH1D(Form("histB_%d", i), Form(";#Delta#sigma;Counts",""), 40, -3, 3);
        TH1D* histMLM = new TH1D(Form("histC_%d", i), "", 40, -3, 3); // No need for title, shared with hischi2

        // deviationsForCases[i].first is the vector for B deviations
        for (double deviation : deviationsForCases[i].first) {
            hischi2->Fill(deviation);
        }
        for (double deviation : deviationsForCasesMLM[i].first) {
            histMLM->Fill(deviation);
        }

        hischi2->SetLineColor(kRed);
        histMLM->SetLineColor(kBlue);
        hischi2->SetStats(false); // Hide stats box
        histMLM->SetStats(false);

        // Inside your loop, before drawing the histograms
        if (i >= 4) { // Adjust for pads 5 and 6
            hischi2->GetXaxis()->SetRangeUser(-2.95, 3);
            histMLM->GetXaxis()->SetRangeUser(-2.95, 3);
        }

        hischi2->Draw();
        histMLM->Draw("SAME");
        // Find the maximum bin content in hischi2 and scale it
		double maxValB = hischi2->GetMaximum();
		hischi2->SetMaximum(maxValB * 1.25); // Set Y-axis max to 1.25 times the max bin content

        // Adjust TPaveText position based on column
        double textStartX = 0.175; // Default for left column
        if (i % 3 == 1) textStartX = 0.10; // Adjust for middle column
        if (i % 3 == 2) textStartX = 0.10; // Adjust for right column
        
        TPaveText* pt = new TPaveText(textStartX, 0.725, textStartX + 0.55, 0.95, "NDC");
        pt->SetBorderSize(1); // Enable border
        pt->SetFillColor(0); // Set fill color to white or transparent
        pt->SetFillStyle(1001); // Solid fill
        pt->SetShadowColor(1); // Enable shadow, usually black
        pt->SetTextAlign(12);
        pt->SetTextSize(0.04);
        pt->AddText(Form("Exclusion: %.1f%%", exclusionPercentage));
        // For chi^2 lines in red
        TText *chi2Text1 = pt->AddText(Form("#chi^{2}, A_{UU}^{cos#phi}: #mu=%.2f, #sigma=%.2f", hischi2->GetMean(), hischi2->GetStdDev()));
        chi2Text1->SetTextColor(kRed);
        // For MLM lines in blue
        TText *mlmText1 = pt->AddText(Form("MLM, A_{UU}^{cos#phi}: #mu=%.2f, #sigma=%.2f", histMLM->GetMean(), histMLM->GetStdDev()));
        mlmText1->SetTextColor(kBlue);
        pt->Draw();
    }

    std::ostringstream filename;
    filename << "output/systematic_study_cosphi_B=" << B << "_C=" << C << ".png";
    canvas->SaveAs(filename.str().c_str());

    delete canvas; // Clean up
}

void plotDeviationsDistributionsC(double B, double C) {
    TCanvas* canvas = new TCanvas("canvas", "Deviations Distribution", 1200, 800);
    canvas->Divide(3, 2, 0.0, 0.0); // No gaps between pads

    int exclusionSteps[] = {0, 1, 2, 3, 4, 5}; 
    for (int i = 0; i < 6; ++i) {
        canvas->cd(i + 1);
        TPad *pad = (TPad*)gPad;

        // Set margins to remove space between subplots. Adjust values as needed.
        float leftMargin = (i % 3 == 0) ? 0.15 : 0.0;
        float rightMargin = 0.0;
        float topMargin = (i < 3) ? 0.05 : 0.0;
        float bottomMargin = (i >= 3) ? 0.15 : 0.0;
        pad->SetMargin(leftMargin, rightMargin, bottomMargin, topMargin);

        // The rest of the code to draw histograms and legends goes here
        double exclusionPercentage = (exclusionSteps[i] / 24.0) * 100.0;
        TH1D* hischi2 = new TH1D(Form("histB_%d", i), Form(";#Delta#sigma;Counts",""), 40, -3, 3);
        TH1D* histMLM = new TH1D(Form("histC_%d", i), "", 40, -3, 3); // No need for title, shared with hischi2

        // deviationsForCases[i].first is the vector for B deviations
        for (double deviation : deviationsForCases[i].second) {
            hischi2->Fill(deviation);
        }
        for (double deviation : deviationsForCasesMLM[i].second) {
            histMLM->Fill(deviation);
        }

        hischi2->SetLineColor(kRed);
        histMLM->SetLineColor(kBlue);
        hischi2->SetStats(false); // Hide stats box
        histMLM->SetStats(false);

        // Inside your loop, before drawing the histograms
        if (i >= 4) { // Adjust for pads 5 and 6
            hischi2->GetXaxis()->SetRangeUser(-2.95, 3);
            histMLM->GetXaxis()->SetRangeUser(-2.95, 3);
        }

        hischi2->Draw();
        histMLM->Draw("SAME");
        // Find the maximum bin content in hischi2 and scale it
        double maxValB = hischi2->GetMaximum();
        hischi2->SetMaximum(maxValB * 1.25); // Set Y-axis max to 1.25 times the max bin content

        // Adjust TPaveText position based on column
        double textStartX = 0.175; // Default for left column
        if (i % 3 == 1) textStartX = 0.10; // Adjust for middle column
        if (i % 3 == 2) textStartX = 0.10; // Adjust for right column
        
        TPaveText* pt = new TPaveText(textStartX, 0.725, textStartX + 0.55, 0.95, "NDC");
        pt->SetBorderSize(1); // Enable border
        pt->SetFillColor(0); // Set fill color to white or transparent
        pt->SetFillStyle(1001); // Solid fill
        pt->SetShadowColor(1); // Enable shadow, usually black
        pt->SetTextAlign(12);
        pt->SetTextSize(0.04);
        pt->AddText(Form("Exclusion: %.1f%%", exclusionPercentage));
        // For chi^2 lines in red
        TText *chi2Text1 = pt->AddText(Form("#chi^{2}, A_{UU}^{cos2#phi}: #mu=%.2f, #sigma=%.2f", hischi2->GetMean(), hischi2->GetStdDev()));
        chi2Text1->SetTextColor(kRed);
        // For MLM lines in blue
        TText *mlmText1 = pt->AddText(Form("MLM, A_{UU}^{cos2#phi}: #mu=%.2f, #sigma=%.2f", histMLM->GetMean(), histMLM->GetStdDev()));
        mlmText1->SetTextColor(kBlue);
        pt->Draw();
    }

    std::ostringstream filename;
    filename << "output/systematic_study_cos2phi_B=" << B << "_C=" << C << ".png";
    canvas->SaveAs(filename.str().c_str());

    delete canvas; // Clean up
}



int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <B> <C>" << std::endl;
        return 1;
    }
    double B = atof(argv[1]);
    double C = atof(argv[2]);

    // Run the acceptance study n times
    acceptanceStudy(B, C, 10000);

    // In your main function or at the end of acceptanceStudy
	plotDeviationsDistributionsB(B, C);
    plotDeviationsDistributionsC(B, C);

}
