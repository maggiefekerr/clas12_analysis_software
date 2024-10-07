#ifndef PLOT_DVCS_DATA_MC_COMPARISON_H
#define PLOT_DVCS_DATA_MC_COMPARISON_H

#include <string>
#include <TTreeReader.h>

void plot_dvcs_data_mc_comparison(const std::string& output_dir, int xB_bin, TTreeReader& data_reader, TTreeReader& mc_gen_reader, TTreeReader& mc_rec_reader);

#endif