#ifndef __ANALYSIS_OMEGA_H__
#define __ANALYSIS_OMEGA_H__

#include "analysis.h"
#include "analysis_heep.h"
//#include "cut.h"


class Analysis_omega: public Analysis_heep {


	public:

		Analysis_omega();
		Analysis_omega(ReadFile::efficiency_profile eff_struc);

		Analysis_omega(TString);

		void Omega_anna(Int_t run_itt);

		~Analysis_omega();

 		void Yield_Out();

		Int_t st_itt;

 		void Omega_sim_anna(Int_t);

		TString weight_app;

 	protected:
// 
// 
// 		void Acceptance_check(std::vector<TString>, TString, TString); 
// 
// 		void Apt_Setting_Check(TString, TString, TString);
// 
// 		Float_t Get_Scale_Factor(){return data_scale_factor; }
// 		void Set_Scale_Factor(Float_t factor){data_scale_factor = factor; }
// 
 	private:
// 
 		Int_t run_num;

// 		void Yield_Out();

// 
// 		/*--------------------------------------------------*/
// 		/// Define real and random cut
// 
 		Double_t coin_center;

 		TString target;	
// 
 		void Define_Cuts();

		void Para_Init();

 		Int_t Calculate_Yield(TH2F* hist_target, Float_t hsbeta_limit);

		Int_t Calculate_Yield_Direct();
	
		


		TCutG* coin_cut;
		TCutG* rand_early_cut;
		TCutG* rand_late_cut;

		TCutG* coin_cor_cut;
		TCutG* rand_early_cor_cut;
		TCutG* rand_late_cor_cut;

//
 		Float_t Calculate_miss_x (Float_t x_boundary, Float_t y_boudary, Float_t y_pos);
 
		float beta_cut;                  
		float gradi;
		float global_centered_hsbeta_cut;
 
 		float yield; 
 		float yield_err; 

 		void Print_out_data();
 
 		Double_t yield_setting; 
 
 		float yield_setting_err; 
 		
 		Double_t charge_setting; 
		
 		void Check_Plot_Out(TString);
 		void Check_Plot_Out(TString plot_var_1, TString plot_var_2);
// 
// 		TString sim_selection;
// 
// 		TString acceptence_cut; 
// 
// 
		void DrawGoodBox(float x_center, float y_center);

		void DrawRandomBox(float x_center, float y_center);
// 
 		Float_t Calculate_tail_x (Float_t xx, Float_t yy);
// 
 		Float_t coin_center_off, beta_center_off;
// 
 		Float_t dummy_correction;
// 		
// 		Float_t data_scale_factor;
// 
 		std::vector<TString> Get_Acceptence_Vector();
// 
// 
 		void Yield_output();


		TGraphErrors* graph_yield_check;
		TGraphErrors* mm_peak_check;

		void Check_MissingMass_Peak();

		void Correct_Offset_Tree();

		TList* list;

		TH2F* coin_beta_check_offset_setting;


		Float_t rand_early_l;
		Float_t rand_early_r;
		Float_t rand_late_l;
		Float_t rand_late_r;
	
		Float_t cointime_cut_l;
		Float_t cointime_cut_r;

		Double_t top_limit;
		Double_t bot_limit;

		Float_t random_real_ratio;

		Double_t GetBetaCenter();

		Double_t hsbeta_center;

		TMultiGraph* diamond_setting;

		void Overall_Dia_Plot();

		TTree* data_tree_in_cl;
		TBranch* mm_cor_branch;
	
		TCutG* diamond_cut_160; 

		TCutG* diamond_cut_245; 

		TCutG* diamond_cut;
	
		TString Diamond_Cut;
		TString Cointime_Cut;
		TString Rand_Cut;


	
		void Overall_Dia_Plot_After(TTree*);   /// Setting Diamond plot
		void Setting_u_Phi_plot(TTree*);       /// Setting Plot u-phi bullseye plot
 		void Setting_Beta_Cointime_Check();
		void Setting_Beta_Cointime_Check_Setting_Tree(TTree*);
 		void Setting_Check_Plot_Out(TTree*, TString, Int_t, Double_t, Double_t);




		Int_t   kset;

 		Int_t   u_bin_set;
 		Float_t t_max, t_min, t_width;
   
 		Int_t   phi_bin_num;

		Float_t phi_stp, phi_min, phi_max, phi_offset;  

		std::vector<float> yield_vec, yield_err_vec, phi_vec, tb_vec;


		
		/// Only plot missmass in u-phi bins
		void Yield_u_Phi_Bin(TTree*);

		/// Only a given kinematics variable in u-phi bins
		void Yield_u_Phi_Bin(TTree* obj_tree, TString plot_var);

//		void Yield_t_Phi_Bin();


		Int_t u_bin_num;
		Float_t* u_upper_limit;
		Float_t* u_lower_limit;

		void Sim_Setting_Check_Plot_Out(TString, Int_t, Double_t, Double_t);

 		void Sim_Setting_Check_U(TString, Int_t, Double_t, Double_t);

		void Sim_Diamond_Plot();

		TH1F* Dummy_Scale(TH1F*); 

   		void U_Distribution_Check(TTree*, TString, Int_t, Double_t, Double_t);
   		void U_Distribution_Check_Sim(TTree*, TString, Int_t, Double_t, Double_t);

		Int_t event;

		Double_t normfac;

		Bool_t is_sim_data;

		Float_t sim_Q2;
		Float_t q2_setting;

// 		Double_t sim_Q2;
// 		Double_t q2_setting;

		
		TString omega_pid_cut;
		TString hms_accept_cut;
		TString sos_accept_cut;
		TString acceptence_cut; 
		TString missmass_cut; 

		TString full_sim_cut;
	
		TString omega_u_cut;


		/// For analyzing the simulation data 
		TString sim_selection;
		TString sim_selection_no_dia;
		TString sim_scale_fac_str;


		TTree* sim_tree_acc;
		TTree* sim_tree_acc_dia;


		TString weight_str;
		TString u_plot_var;
		TString t_plot_var;

		void Sim_Set_Histo_Err(TH1F*, TH1F*);

		void Overall_2_D_Plot(TTree*, TString, TString);

		void Sim_2_D_Plot(TString, TString);

		
        void Yield_u_Phi_Bin_Sim(TTree* obj_tree);
        void Yield_u_Phi_Bin_Sim(TTree* obj_tree, TString plot_var);






};





class Analysis_omega_itt: public Analysis_omega {
 
 	public:
 
 		Analysis_omega_itt(); 
 		Analysis_omega_itt(TString);
 		~Analysis_omega_itt();

		


};

#endif 


