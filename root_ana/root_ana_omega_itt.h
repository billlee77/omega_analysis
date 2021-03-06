#ifndef ROOT_ANA_HEEP_itt_H
#define ROOT_ANA_HEEP_itt_H


#include <iostream>
#include <string>
#include <algorithm>

#include "TTree.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TString.h"

#include "read_setting.h"
#include "analysis.h"
#include "analysis_omega.h"
#include "print_data.h"


void Initialization();

ReadFile::kinematics_profile kin_setting;
//ReadFile::efficiency_profile eff_setting;
ReadFile::center_profile cen_setting;

ReadFile::sim_profile sim_norm_pro;


ReadFile::eff_pro* eff_setting;

ReadFile::eff_pro* eff_setting_dummy;
ReadFile::eff_pro* eff_setting_target;



Int_t num_runs;
Int_t set_runs;

TString Dummy_Name(int);
TString Target_Name(int);
TString List_Name(int);
TString Sim_Name(int);

TString Target_Name(TString, int);
TString Standard_Name_Str(TString, TString, Double_t, Double_t, Double_t);
TString Standard_Name_Str(TString, TString, TString,  Double_t);

TFile **output_file;
TFile **output_file_dummy;
TFile **output_file_target;
TFile **output_file_sim;

int Return_Setting(int, int, int);
int Return_Setting_dummy(int, int, int);
int Return_Setting_target(int, int, int);

void Fill_Setting_lst();



std::vector<TString> Q2_vec, Q2_vec_mod;
std::vector<TString>::iterator Q2_vec_mod_itt;

std::vector< std::vector<int> > setting_run_list;
std::vector< std::vector<int> > setting_run_list_dummy;
std::vector< std::vector<int> > setting_run_list_target;

std::vector<TString>  lists_vec, setting_vec;

std::vector<TString>  setting_vec_target, setting_vec_dummy;

std::vector<TString>::iterator lists_vec_itt;

std::vector<TString>  lists_vec_target, lists_vec_dummy;



std::vector<TString> sim_vec;

TString list_name;
TString sim_ana_name;

// Analysis **ana;
// Analysis **ana_target;
// Analysis **ana_dummy;


// Analysis_heep **ana;
// Analysis_heep **ana_target;
// Analysis_heep **ana_dummy;
// Analysis_heep **ana_sim;
// 

Analysis_omega_itt **ana_sim;


ReadFile* rf;
ReadFile* sim_rf;


File_Output *fo; 


TString target_name;
bool dummy_tar;

TString list_dir;

std::vector<int> run_list_vec;
std::vector<int> run_list_vec_target, run_list_vec_dummy;

void Cout_out_lists();

int setting_num;

/// Loading Functions
void Sim_Load();


/// Analysis Functions]
void Sim_Analysis();


void Cleanup();


TString sim_file_inc;


#endif
