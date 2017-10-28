#ifndef read_setting_H
#define read_setting_H

#include <string>

#include "TString.h"
#include "TNtuple.h"


class ReadFile {

	public:

		ReadFile();
		ReadFile(TString);
		ReadFile(TString eff_file_name, TString offset_file_name);
		
		void Read_init();

		void Setting_file_loading();
		void Setting_file_loading_name();

		void Eff_file_loading(TString list_file_in);

		void Kin_Pro_Array_Load();
		void Eff_Pro_Array_Load();


		void Cen_Pro_Array_Load();

		void Sim_Pro_Array_Load();
		void Sim_Pro_Array_Omega_Load();

		void Calculate_t_Width();



		/*--------------------------------------------------*/
		// Data list file in and directory

		TString list_dir;
		TString list_file;
 
		struct kinematics_profile {

			Double_t* polset; 
			Double_t* Q2set;  
			Double_t* epsset; 
			                 
			Double_t* thpqset;
			Double_t* tmnset; 
			Double_t* tmxset; 
			                 
			Int_t* NBtset; 
			Double_t* kset;

		} kin_pro;


// 		struct kinematics_profile get() {
//   			return {0};
// 		}


		typedef struct efficiency_profile {

			Int_t* run_num;  

			Double_t* Q2;
			Double_t* ebeam;  

			Double_t* charge;
			Double_t* charge_err;  

			Double_t* tot_eff;     
			Double_t* tot_eff_err; 

			Double_t* epsset;     
			Double_t* thpqset; 	  	  

		} eff_pro;


		eff_pro eff_pro1;



		struct center_profile {

			Int_t*  center_run_num;
			Double_t*  center_mean;     
			Double_t*  center_sigma;    

		} cen_pro;



		struct sim_profile {

			Double_t* q2_setting;  
			Double_t* ebeam;  

			Double_t* epsilon;  
			Double_t* hms_theta;  

			Double_t* normfac;     
			Int_t* event_num; 	  	  

		} sim_pro;


		~ReadFile();

		Double_t* Return_array_D(TString variable_name, TTree* target_tree);

		Int_t* Return_array_I(TString variable_name, TTree* target_tree);



		Int_t t_bin_set  ;
		Double_t t_min   ; 
		Double_t t_max   ; 
		Double_t t_width ; 


		Int_t Get_Set_Runs() { return set_runs; };

		Int_t Get_Num_Runs() { return num_runs; };

		Int_t Get_Cen_Runs() { return cen_runs; };

		bool Get_Dummy() { return dummy_tar; };

		TString Get_List_File() { return list_file; };

		bool dummy_tar;

		void Ntuple_Reset();

		inline bool File_exist(TString);
		inline bool File_empty(TString);

// 		void Set_Eff_File(TString name) { eff_file_name = name; };
// 		void Set_Off_File(TString name) { offset_file_name = name; };

	private:

		std::string list_name;

		TNtuple* setting_ntp;
		TNtuple* eff_ntp;
		TNtuple* center_ntp; 
		TNtuple* sim_ntp; 

		Int_t set_runs;
		Int_t num_runs;
		Int_t cen_runs;

		TString eff_file_name;
		TString off_file_name;
		TString sim_file_name;


};


#endif
