// using the binning from G. Christiaens thesis (24 bins) to extract the correct binning

// getting the mean Q2, xB, t values for each of the three runs
void meanVals_forRun(TString fileName, TString runName) {
    TString outputName = Form("%s_meanVals_Q2_xB_t.txt", runName.Data());
    std::ofstream of_meanVals(outputName);

    Double_t minxBList[24] = {0.0, 0.16, 0.26, 0.0, 0.16, 0.26, 0.0, 0.16, 0.26, 0.0, 0.16, 0.26, 0.0, 0.16, 0.26, 0.0, 0.16, 0.26, 0.0, 0.16, 0.26, 0.0, 0.16, 0.26};
    Double_t maxxBList[24] ={0.16, 0.26, 1.0, 0.16, 0.26, 1.0, 0.16, 0.26, 1.0, 0.16, 0.26, 1.0, 0.16, 0.26, 1.0, 0.16, 0.26, 1.0, 0.16, 0.26, 1.0, 0.16, 0.26, 1.0};
    Double_t mintList[24] = {-0.15, -0.22, -0.40, -0.25, -0.40, -0.70, -0.45, -0.80, -1.15, -2.5, -2.5, -2.5, -0.15, -0.22, -0.40, -0.25, -0.40, -0.70, -0.45, -0.80, -1.15, -2.5, -2.5, -2.5};
    Double_t maxtList[24] = {0.0, 0.0, 0.0, -0.15, -0.22, -0.40, -0.25, -0.40, -0.70, -0.45, -0.80, -1.15, 0.0, 0.0, 0.0, -0.15, -0.22, -0.40, -0.25, -0.40, -0.70, -0.45, -0.80, -1.15};
    Double_t minQ2List[24] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.75, 2.40, 3.25, 1.75, 2.40, 3.25, 1.75, 2.40, 3.25, 1.75, 2.40, 3.25};
    Double_t maxQ2List[24] = {1.75, 2.40, 3.25, 1.75, 2.40, 3.25, 1.75, 2.40, 3.25, 1.75, 2.40, 3.25, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,10.0, 10.0, 10.0, 10.0, 10.0, 10.0};

    /*for (int i=0; i<24; i++){
        cout << minxBList[i] << " " << maxxBList[i] << " " << mintList[i] << " " << maxtList[i] << " " << minQ2List[i] << " " << maxQ2List[i] << endl;
    }*/

    TFile *file = new TFile(fileName); 
    TTree *t = (TTree*)file->Get("PhysicsEvents");

    for (int i=0; i<24; i++) {
        TH1F *hQ2 = new TH1F("hQ2", "hQ2", 100, minQ2List[i], maxQ2List[i]);
        TH1F *hxB = new TH1F("hxB", "hxB", 100, minxBList[i], maxxBList[i]);
        TH1F *ht = new TH1F("ht", "ht", 100, mintList[i], maxtList[i]);

        TString cuts = Form("x > %0.3f && x < %0.3f && Q2 > %0.3f && Q2 < %0.3f && t1 > %0.3f && t1 < %0.3f", minxBList[i], maxxBList[i], minQ2List[i], maxQ2List[i], mintList[i], maxtList[i]);

        t->Draw("Q2>>hQ2", cuts, "goff");
        t->Draw("x>>hxB", cuts, "goff");
        t->Draw("t1>>ht", cuts, "goff");

        cout << hQ2->GetMean() << " " << hxB->GetMean() << " " << ht->GetMean() << endl;
        of_meanVals << hQ2->GetMean() << " " << hxB->GetMean() << " " << ht->GetMean() << endl;
    }
}

// getting weighted average across the runs
void meanVals(){

}