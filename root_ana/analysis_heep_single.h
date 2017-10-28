#ifndef __ANALYSIS_HEEP_SINGLE_H__
#define __ANALYSIS_HEEP_SINGLE_H__

//#include "analysis.h"
#include "analysis_heep.h"
//#include "cut.h"


//class Analysis_heep_single: public Analysis_heep {

class Analysis_heep_single: public Analysis_heep {


	public:

		Analysis_heep_single();

		Analysis_heep_single(ReadFile::efficiency_profile eff_struc);

 		void Heep_anna_single(Int_t run_itt);
 
		void Heep_sim_anna_single(Int_t);

		void Yield_Out(); 

 		~Analysis_heep_single();


		std::vector<TString> check_list_vec;



	private:

 		Int_t run_num;

		void Load_Data_Tree();

		void Define_Cuts();

		TString single_cut;

		void Check_Plot_Out(TString);

		float ssdelta, ssyptar, ssxptar, w_data, q2_data;
		float sszbeam, ssp, ssenergy;
		float evtype;

		float ssxtar_data, ssytar_data;		

		float Get_SOS_Momentum(float q2_tmp); 

		float Get_SOS_Angle(float q2_tmp);

		void Check_Out(TTree*);
		void Write_his();

		TH1F * ssdelta_h;    
		TH1F * ssxptar_h;    
		TH1F * ssyptar_h;    
		                    
		TH1F * w_h;          
		TH1F * w_h_org;          
		TH1F * q2_h;         
		                    
		TH1F * ssenergy_h;   
		TH1F * sszbeam_h;    
		
		TH1F * ssp_h;        
		                    
		TH1F * missingmass_h;
		                    
		TH1F * E_m_h;        
		TH1F * P_m_h;       
		                    
		TH1F * P_mx_h;       
		TH1F * P_my_h;       
		TH1F * P_mz_h;       
		                    
		TH2F * w_p;          
		TH2F * delta_m;  

		TH1F * test_h;

		TH1F * ssxtar_h;
		TH1F * ssytar_h;


		TTree* data_tree_in_cut;

		float ssxptar_offset, ssyptar_offset;
		
		float Get_Real_SOS_P(float);

		TString sim_selection;

		TString acceptence_cut; 

//		void Acceptance_check(std::vector<TString>, TString, TString); 

//		void Apt_Setting_Check(TString, TString, TString);

		std::vector<TString> Get_Acceptence_Vector();

//		float data_scale_factor;

		void Calculate_Yield(TH1F*); 

		float yield;
		float yield_err;

		float dummy_correction;

		float w_cut_up, w_cut_down;

		void Set_Hist_Err(TH1F*);

//		void Yield_output();

};

#endif 
