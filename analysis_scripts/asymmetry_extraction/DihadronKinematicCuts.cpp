#include "DihadronKinematicCuts.h"
#include "common_vars.h"
#include "BaseKinematicCuts.h" // Include BaseKinematicCuts
#include <string>
#include <cmath>

using std::string;

DihadronKinematicCuts::DihadronKinematicCuts(TTreeReader& reader)
    : BaseKinematicCuts(reader), // Initialize BaseKinematicCuts
      runnum(reader, "runnum"), Q2(reader, "Q2"), W(reader, "W"), 
      x(reader, "x"), y(reader, "y"), z1(reader, "z1"), z2(reader, "z2"), Mh23(reader, "Mh23"), 
      xF1(reader, "xF1"), xF2(reader, "xF2"), Mx(reader, "Mx"), target_pol(reader, "target_pol") {}

bool DihadronKinematicCuts::applyCuts(int currentFits, bool isMC) {
    bool goodEvent = false;
    string property = binNames[currentFits];

    if (property == "Mh") {
        goodEvent = *Q2 > 1 && *W > 2 && *y < 0.75;
    } else if (property == "epippimX") {
        goodEvent = *Q2 > 1 && *W > 2 && *y < 0.75 && *z1 > 0.2 && *z2 > 0.2;
    } else if (property == "xepippimX") {
        goodEvent = *Q2 > 1 && *W > 2 && *y < 0.75 && *z1 > 0.2 && *z2 > 0.2 && 
            *xF1 > 0 && *xF2 > 0;
    } else if (property == "ekpkmX") {
        goodEvent = *Q2 > 1 && *W > 2 && *y < 0.75 && *z1 > 0.2 && *z2 > 0.2;
    } else if (property == "xekpkmX") {
        goodEvent = *Q2 > 1 && *W > 2 && *y < 0.75 && *z1 > 0.2 && *z2 > 0.2 && 
            *xF1 > 0 && *xF2 > 0;
    } else {
      std::cout << "Property, " << property << ", not detected." << std::endl;
    }
    
    if (isMC || (*runnum < 16042 || *runnum > 17811)) {
      return goodEvent;
    } else {
      // return goodEvent && *target_pol!=0;
        return goodEvent;
    }
    return false;
}