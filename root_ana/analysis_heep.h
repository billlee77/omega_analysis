#ifndef __ANALYSIS_HEEP_H__
#define __ANALYSIS_HEEP_H__

#include "analysis.h"
//#include "cut.h"


class Analysis_heep: public  Analysis {


	public:

		Analysis_heep();
		Analysis_heep(ReadFile::efficiency_profile eff_struc);

		void Heep_anna(Int_t run_itt);

		~Analysis_heep();

		void Yield_Out();

		void Heep_sim_anna(Int_t);

	protected:

		void Yield_output();

		void Acceptance_check(std::vector<TString>, TString, TString); 

		void Apt_Setting_Check(TString, TString, TString);

		Float_t Get_Scale_Factor(){return data_scale_factor; }
		void Set_Scale_Factor(Float_t factor){data_scale_factor = factor; }

		Double_t yield_setting; 

		Double_t yield_setting_err; 

		void Correct_Offset_Tree();
		
		TTree* Create_File(); 

		void Proton_Obsorption_Study();

		void Proton_Obsorption_Study(TTree*);
		void hsbeta_Omega_Cut_Study();

	private:

		Int_t run_num;

		/*--------------------------------------------------*/
		/// Define real and random cut

		Double_t coin_center;

// 		TString all_coin_cut; 
// 		TString real_coin_cut;
// 		TString rand_coin_cut;

		TString target;	

		void Define_Cuts();


// 		void Para_Init();
// 
// 		Double_t* miss_mass_offset;
// 		float acccc_temp, errrr_temp; 
// 
// 		Double_t coin_min;
// 		Double_t coin_max;
// 		Int_t    coin_bin;
// 
// 		Int_t t_bin_set;
// 		float t_max, t_min, t_width;
//   
// 		std::vector<float> yield_vec, yield_err_vec, phi_vec, tb_vec;
// 
// 		float pstp, pmn, pmx;  
// 
// 		float Q_2;
// 
// 		int kset;

		Int_t Calculate_Yield(TH2F* hist_target, Float_t hsbeta_limit);

		Float_t Calculate_miss_x (Float_t x_boundary, Float_t y_boudary, Float_t y_pos);

		float beta_cut;                  
		float gradi;
		float global_centered_hsbeta_cut;

		float yield; 
		float yield_err; 
		

		void Print_out_data();


		Double_t charge_setting; 
		
		
		void Check_Plot_Out(TString);
		void Check_Plot_Out(TString plot_var_1, TString plot_var_2);

		void Sim_Setting_Check_Plot_Out(TString, Int_t, Double_t, Double_t);


		void Setting_Check_Plot_Out(TString, Int_t, Double_t, Double_t);
		void Setting_Check_Plot_Out(TString);

		void Setting_Beta_Cointime_Check();

		void Setting_Beta_Cointime_Check_Setting_Tree(TTree*);

		//void GetBetaCenter(); 






		TString sim_selection;

		TString acceptence_cut; 
		TString acceptence_missmass_cut; 

		void DrawGoodBox(float x_center, float y_center);

		void DrawRandomBox(float x_center, float y_center);

		Float_t Calculate_tail_x (Float_t xx, Float_t yy);

		Float_t coin_center_off, beta_center_off;

		Float_t dummy_correction;
		
		Float_t data_scale_factor;

		std::vector<TString> Get_Acceptence_Vector();



		Int_t Calculate_Yield_Direct(); 

		TCutG* coin_cut;
		TCutG* rand_early_cut;
		TCutG* rand_late_cut;

		TCutG* coin_cor_cut;
		TCutG* rand_early_cor_cut;
		TCutG* rand_late_cor_cut;

		Float_t rand_early_l;
		Float_t rand_early_r;
		Float_t rand_late_l;
		Float_t rand_late_r;
	
		Float_t cointime_cut_l;
		Float_t cointime_cut_r;

		Double_t top_limit;
		Double_t bot_limit;
		
		Float_t random_real_ratio;

		TString Diamond_Cut;
		TString Cointime_Cut;
		TString Rand_Cut;

		Double_t hsbeta_center;

		TString hms_accept_cut;
        TString sos_accept_cut;
        TString heep_pid_cut;

		TTree* data_tree_in_cl;
		TBranch* mm_cor_branch;

};

#endif 
