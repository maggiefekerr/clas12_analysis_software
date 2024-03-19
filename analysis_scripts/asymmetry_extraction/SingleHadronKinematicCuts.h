#pragma once
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <string>
#include "common_vars.h"
#include "BaseKinematicCuts.h" // Include BaseKinematicCuts

class SingleHadronKinematicCuts : public BaseKinematicCuts {
public:
    SingleHadronKinematicCuts(TTreeReader& reader);
    bool applyCuts(int currentFits, bool isMC) override;

private:
    TTreeReaderValue<int> runnum;
    TTreeReaderValue<double> e_theta;
    TTreeReaderValue<double> e_phi;
    TTreeReaderValue<double> p_p;
    TTreeReaderValue<double> p_theta;
    TTreeReaderValue<double> p_phi;
    TTreeReaderValue<double> Q2;
    TTreeReaderValue<double> W;
    TTreeReaderValue<double> Mx;
    TTreeReaderValue<double> x;
    TTreeReaderValue<double> y;
    TTreeReaderValue<double> z;
    TTreeReaderValue<double> pT;
    TTreeReaderValue<double> xF;
    TTreeReaderValue<double> phi;
    TTreeReaderValue<double> target_pol;
};