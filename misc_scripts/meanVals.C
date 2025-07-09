// extracting mean values from the binnings described in the BSA result for RG-A
// this will be rough means as the data is just with the loose exclusivity cuts and not with the exclusive cuts, pi0 subtraction, and radiative correction

// Q2 < 1.4; xB: < 0.13, [0.13, 0.21], > 0.21 
// Q2 = [1.4, 1.8]; xB: < 0.13, [0.13, 0.21], > 0.21 
// Q2 = [1.8, 2.4]; xB: < 0.16, [0.16, 0.26], > 0.26 
// Q2 = [2.4, 3.2]; xB: < 0.21, [0.21, 0.33], > 0.33 
// Q2 = [3.2, 5.0]; xB: < 0.33, > 0.33
// Q2 > 5.0; xB: < 0.56, > 0.56

void meanVals(TString fileName, TString runName) {
    TString outputName = Form("%s_meanVals_Q2_xB_t.txt", runName.Data());
    std::ofstream of_meanVals(outputName);
    Double_t mp = 0.93827;
    Double_t Q2List[7] = {0.0001, 1.4, 1.8, 2.4, 3.2, 5.0, 10.0};
    Int_t xB_numBins[6] = {3, 3, 3, 3, 2, 2};
    Double_t xBList[6][4] = {   {0.0001, 0.13, 0.21, 1.0},
                                {0.0001, 0.13, 0.21, 1.0},
                                {0.0001, 0.16, 0.26, 1.0},
                                {0.0001, 0.21, 0.33, 1.0},
                                {0.0001, 0.33, 1.0, 0.0}, // setting last point to zero, will handle this in for loop
                                {0.0001, 0.56, 1.0, 0.0}  // setting last point to zero, will handle this in for loop
                            };
    TFile *file = new TFile(fileName); 
    TTree *t = (TTree*)file->Get("PhysicsEvents");

    Int_t num = 0;
    
    for (int i=0; i<5; i++) {
        Double_t Q2min = Q2List[i];
        Double_t Q2max = Q2List[i+1];
        for(int j=0; j<xB_numBins[i]; j++) {
            Double_t xBmin = xBList[i][j];
            Double_t xBmax = xBList[i][j+1];
        
            Double_t eps = 2.*xBmin*Q2min/mp;
            Double_t tposMin = (Q2min*2.*(1.-xBmin)*(1.-TMath::Sqrt(1.+eps*eps))+(eps*eps))/(4.*xBmin*(1.-xBmin)+(eps*eps));
            Double_t tposMax = (Q2min*2.*(1.-xBmin)*(1.+TMath::Sqrt(1.+eps*eps))+(eps*eps))/(4.*xBmin*(1.-xBmin)+(eps*eps));
            Double_t tDiff = tposMax - tposMin;
            cout << eps << " " << tposMin << " " << tposMax << " " << tDiff << endl;
            Double_t tposList[5] = {tposMin, tposMin + tDiff/4., tposMin + tDiff/2., tposMin + (3./4.)*tDiff, tposMax}; // these are all positive values
            for (int k=0; k<4; k++){
                Double_t tposmin = tposList[k];
                Double_t tposmax = tposList[k+1];
                cout << Q2min << " " << Q2max << " " << xBmin << " " << xBmax << " " << tposmin << " " << tposmax << endl;
                /*TH1F *hQ2 = new TH1F("hQ2", "hQ2", 100, Q2min, Q2max);
                TH1F *hxB = new TH1F("hQ2", "hQ2", 100, xBmin, xBmax);
                TH1F *ht = new TH1F("hQ2", "hQ2", 100, tposmin, tposmax);

                TString cuts = Form("x > %0.3f && x < %0.3f && Q2 > %0.3f && Q2 < %0.3f && t1 > %0.3f && t1 < %0.3f", xBmin, xBmax, Q2min, Q2max, -1*tposmax, -1*tposmin);

                t->Draw("Q2>>hQ2", cuts, "goff");
                t->Draw("x>>hxB", cuts, "goff");
                t->Draw("t1>>ht", cuts, "goff");

                of_meanVals << hQ2->GetMean() << " " << hxB->GetMean() << " " << ht->GetMean() << endl;*/
                num += 1;
            }
        }
    }
    cout << num << endl;
}