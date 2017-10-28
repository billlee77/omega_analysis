#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <iostream>
#include <vector>
#include <string>


#include <fstream>

#include <stdio.h>      
#include <stdlib.h> 

#include <math.h>       
#include <cmath>
#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TLine.h"
#include "TPad.h"
#include "TGraph.h"
#include "TGraphPolar.h"

#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TAttLine.h"
#include "TF1.h"
#include "TText.h"
#include "TPaveText.h"
#include "TPaletteAxis.h"
#include "TPaveStats.h"
#include "TLatex.h"

#include "TList.h"
#include "TFile.h"
#include "TChain.h"
#include "TCutG.h"

//#include "root_ana_pl.h"
#include "read_setting.h"
#include "print_data.h"

#include "cut.h"

class Analysis {

	public:

		Analysis();
		Analysis(ReadFile::efficiency_profile eff_struc);

		~Analysis();

		void Init();
		void Run_by_Run_Analysis(Int_t run_num);
		void Yield_Out();

		ReadFile::kinematics_profile kin_ana;
		ReadFile::efficiency_profile eff_ana;
		ReadFile::center_profile 	 cen_ana;
		ReadFile::sim_profile 	 	 sim_ana;

		bool is_run_dummy;

		void Load_Data_Tree();

		void Para_Run_Def(Int_t num);
		void Para_Run_Clean();

		void Print_out_data();
		void Missing_Mass_Plot();
		void Overall_Dia_Plot();

		TFile* file_out_ana;
		TDirectory* out_dir;

		TH1F* coin_all; 
		TH1F* coin_real;
		TH1F* coin_rand;

		TTree* data_tree_in;

		TTree* run_tree;	
		TTree* tree_out;
		TTree* Ceate_File();

		int phi_bin_num;

		float Get_kset() {return kset;}

		void Para_Init();

		float Get_Coin_Center() {return coin_center;}


		void Set_Coin_Min(Double_t min) { coin_min = min; }
		void Set_Coin_Max(Double_t max) { coin_max = max; }
		void Set_Coin_Bin(Int_t bin) { coin_bin = bin; }
	

		void Set_Expected_MM(Double_t mm_set) { expected_mm = mm_set; }

		void Set_MM_Offset(Double_t* array);

     	TH1F* mm_offset_setting;
     	TH1F* mm_setting;	


//		TString sim_file;

		

// 		void Set_Eff(TString* name) {eff_file = name;};
// 		void Set_Off(TString* name) {off_file = name;};
// 
// 

	protected:

		TString all_coin_cut; 
		TString real_coin_cut;
		TString rand_coin_cut;

		float charge_tot, acccc_temp, errrr_temp; 

		void SetYield(float yield_setting) { yield = yield_setting; }

		Double_t Get_Total_Eff() { return tot_eff; } 
		Double_t Get_Total_Eff_Err() {return tot_eff_err;}

		Double_t Get_Charge() { return charge; }
 		Double_t Get_Charge_Err() {return charge_err; }





 		TString eff_file;
 		TString off_file;

		TChain* chain;
		
		TString Get_Current_Run_File(){return data_file;};
//		void Set_Coin_Bin(TString dir_str) {  data_file_dir = dir_str; }


	
		/*--------------------------------------------------*/
		// Data file in and directory
		TFile* file_in;
		TString data_file;
		TString data_file_dir; 

		Int_t Get_Run_Num()  {return run_num; }
		Float_t Get_E_Beam() {return ebeam_run; }
		Float_t Get_Q2()     {return q2_run; }

		std::ofstream yield_file_out;

		TList* list;
		
		TString sim_ana_name;

	private:

		Int_t run_num;

		Int_t num_runs, cen_runs;

		float charge, tot_eff, charge_err, tot_eff_err;

		Float_t ebeam_run, q2_run;

		Float_t t_missmass, t_W, t_t, t_Q2, t_th_pq, t_phi_pq, t_PmPar, t_PmPer, t_PmOop;
	 	Float_t t_hsdelta, t_hsytar, t_hsyptar, t_hsxptar, t_ssdelta, t_ssytar, t_ssyptar, t_ssxptar; 
		
		
		/*--------------------------------------------------*/
		// File out and firectory
		TString out_file_name;
		TString root_out_dir_name;
		TString list_file;

		float yield, yield_err, phi, tb;

		float  tmnset, tmxset;

		TCanvas *c1;
		TCanvas *c2;
		TCanvas *c3;
		TCanvas *c4;

		TCanvas *phi_no_sub;
		TCanvas *phi_sub;

		std::vector<std::vector<double> > yield_err_sum;

		TH1F* real_event[6];

		Double_t coin_center;

		TString miss_mass_offset_str;
		
		/*--------------------------------------------------*/
		/// Define real and random cut
		
		float real_err, rand_err, real_event_error;



		TString target;	

		TH1F* cointime_setting; 


		TMultiGraph* diamond_setting;    
    	TMultiGraph* diamond_setting_cut;

		/// Cut specifically for Fp2 test code
		TString without_diamond_cut;
		TString yield_t_bin_cut;

		/// Parameter Initialization

		// Double_t coin_min;
		// Double_t coin_max;

		Double_t* miss_mass_offset;

		Int_t t_bin_set;
		float t_max, t_min, t_width;
		std::vector<float> yield_vec, yield_err_vec, phi_vec, tb_vec;

		float pstp, pmn, pmx;  

		float Q_2;
		int kset;

		Double_t coin_min;
		Double_t coin_max;
		Int_t    coin_bin;

		void Define_Cuts();

		Double_t expected_mm;
		
};

#endif 
