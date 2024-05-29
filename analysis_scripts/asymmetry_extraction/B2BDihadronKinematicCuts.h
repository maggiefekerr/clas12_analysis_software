#pragma once
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <string>
#include "common_vars.h"
#include "BaseKinematicCuts.h" 

class B2BDihadronKinematicCuts : public BaseKinematicCuts { // Inherit from BaseKinematicCuts
public:
    B2BDihadronKinematicCuts(TTreeReader& reader);
    bool applyCuts(int currentFits, bool isMC) override; // Override applyCuts method

private:
    TTreeReaderValue<int> runnum;
    TTreeReaderValue<double> Q2;
    TTreeReaderValue<double> W;
    TTreeReaderValue<double> p1_p;
    TTreeReaderValue<double> x;
    TTreeReaderValue<double> y;
    TTreeReaderValue<double> z1;
    TTreeReaderValue<double> xF1;
    TTreeReaderValue<double> xF2;
    TTreeReaderValue<double> Mx;
    TTreeReaderValue<double> Mx1;
    TTreeReaderValue<double> Mx2;
    TTreeReaderValue<double> target_pol;
};