/*--------------------------------------------------*/
// Program for analysing the omega data 
// Author: wl
// Email: wenliang.billlee@gmail.com

#include <iostream>
#include "root_ana_omega_itt.h"

//#include "cut.h"

using namespace std;

int main() {



 	cout << endl << "/*--------------------------------------------------*/ " << endl;
 	cout         << "                   Analysis Starts                     " << endl << endl;
//	cout         << "/*--------------------------------------------------*/ " << endl << endl;
	


 	sim_ana_name = "omega";
 	Initialization();
	Sim_Load();
	Sim_Analysis();
	Cleanup();

//	cout << endl << "/*--------------------------------------------------*/ "  << endl;
	cout << endl << "                   Analysis ends                       "  << endl;
	cout         << "/*--------------------------------------------------*/ "  << endl << endl;

	return 0;

}









/*--------------------------------------------------*/
/*--------------------------------------------------*/
// Initialization 
//

void Initialization() {

	
	TString tmp1;
	TString tmp2;
//	TString sim_file_inc;


//	sim_ana_name = "rho";
//	sim_ana_name = "omega";
//	sim_ana_name = "xphsp";


	tmp1     = "list.settings.omega";
	tmp2     = "offset.dat";
	sim_file_inc = "list_" + sim_ana_name + "_sim_normfac.dat";


	rf     = new ReadFile(tmp1, tmp2);
	sim_rf = new ReadFile(sim_file_inc); 

	fo = new File_Output();

	kin_setting  = rf->kin_pro;
//	eff_setting  = rf->eff_pro;
	cen_setting  = rf->cen_pro;

	sim_norm_pro = sim_rf->sim_pro;

	
	

//	num_runs = rf->Get_Num_Runs();

	set_runs = rf->Get_Set_Runs();

	dummy_tar = rf->Get_Dummy(); 

//  	output_file = new TFile*[set_runs];
//  	ana 		= new Analysis*[set_runs];
//  	eff_setting = new ReadFile::eff_pro[set_runs];

	/*--------------------------------------------------*/
	/// Mod
	
 	output_file_dummy  = new TFile*[set_runs];
 	output_file_target = new TFile*[set_runs];
 	output_file_sim    = new TFile*[set_runs];

 	ana_sim            = new Analysis_omega_itt*[set_runs];
// 
//  	ana_dummy 	       = new Analysis*[set_runs];
//  	ana_target         = new Analysis*[set_runs];


 	eff_setting_dummy  = new ReadFile::eff_pro[set_runs];
 	eff_setting_target = new ReadFile::eff_pro[set_runs];


	cout << "Set Runs !!!  : " << set_runs << endl;

	list_dir 			= "lists/";

	ofstream u_bin_int;
	u_bin_int.open(fo->Get_Dir() + "u_bin_interval", ios::out);

}




/*--------------------------------------------------*/
/*--------------------------------------------------*/

TString Sim_Name(int num_itt) {

//  	TString dir_str;
// 
// 	if (int(sim_norm_pro.hms_theta[num_itt]*1000) >= int(0)) {
// 
// 		dir_str.Form("yields." + sim_ana_name + "_%i_%i_+%4.4i.sim", int(round(sim_norm_pro.q2_setting[num_itt]*100)), int(sim_norm_pro.epsilon[num_itt]*100.),int(sim_norm_pro.hms_theta[num_itt]*1000));	
// 
// 	} else {
// 
// 		dir_str.Form("yields." + sim_ana_name + "_%i_%i_%4.4i.sim", int(round(sim_norm_pro.q2_setting[num_itt]*100)), int(sim_norm_pro.epsilon[num_itt]*100.), int(sim_norm_pro.hms_theta[num_itt]*1000));	
// 
// 	}
// 
// 	return dir_str;

 	TString dir_str_1;

	dir_str_1 = Standard_Name_Str(sim_ana_name, "sim", sim_norm_pro.q2_setting[num_itt], sim_norm_pro.epsilon[num_itt], sim_norm_pro.hms_theta[num_itt]);

	return dir_str_1;

}





/*--------------------------------------------------*/
/*--------------------------------------------------*/

TString List_Name(int num_itt) {

	TString dir_str;

// 
//     cout << kin_setting.Q2set[num_itt]*100 << "   " << kin_setting.Q2set[num_itt] << "    " 
// 		 << kin_setting.epsset[num_itt]*100 << "   " << int(round(kin_setting.epsset[num_itt]*100)) << endl;
 
	dir_str.Form("%i_%i", int(round(kin_setting.Q2set[num_itt]*100)), int(round(kin_setting.epsset[num_itt]*100)));	

	// cout << "gyyyy !" << dir_str << endl;

	return dir_str;

}



/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Sim_Load() {

	for (Int_t i = 0; i < set_runs; i++) {

		if (kin_setting.polset[i] > 0) {

			output_file_sim[i] = fo->Create_Out_File(Sim_Name(i));

//			list_name = List_Name(i) + "omega";
			list_name = List_Name(i) + "_" + sim_ana_name;
	
			sim_vec.push_back(list_name);
				
			ana_sim[i] = new Analysis_omega_itt(sim_ana_name);

			// cout << "sim file :: "  << sim_vec.size()  << "   "  << list_name << endl;

			cout << Sim_Name(i)  << endl;

		}

	}

//	exit(0);

}



/*--------------------------------------------------*/
/*--------------------------------------------------*/


void Sim_Analysis() {

 	for(int i = 0; i < sim_vec.size(); i++) {

		cout << "Sim file:   "  << sim_vec.size() << "    " << sim_vec[i] << endl;

//		output_file_sim[i]->ls();


//		ana_sim[i]->sim_file = sim_file_inc;

		ana_sim[i]->file_out_ana = output_file_sim[i];

		ana_sim[i]->is_run_dummy = true;		
		ana_sim[i]->Omega_sim_anna(i);


//		exit(0);

//		delete ana_sim[i];

  	}


}


/*--------------------------------------------------*/

void Cleanup() {


	delete rf;
	delete sim_rf;
	delete fo;


	delete output_file_dummy; 
    delete output_file_target;
    delete output_file_sim;   
                      
    delete ana_sim;           

	delete eff_setting_dummy;
	delete eff_setting_target;

	sim_vec.clear();

}


/*--------------------------------------------------*/

TString Standard_Name_Str(TString particle_type, TString data_type, Double_t q2, Double_t epsilon, Double_t theta) {

	TString dir_str;

	if (int(theta*1000) >= int(0)) {

		dir_str.Form("yields." + particle_type + "_%i_%i_+%4.4i." + data_type, int(round(q2*100)), int(round(epsilon*100)), int(round(theta*1000)));	

	} else {

		dir_str.Form("yields." + particle_type + "_%i_%i_%4.4i." + data_type, int(round(q2*100)), int(round(epsilon*100)), int(round(theta*1000)));	

	}

	return dir_str;

}



/*--------------------------------------------------*/
TString Standard_Name_Str(TString particle_type, TString data_type, TString q2_eps_str, Double_t theta) {

	TString dir_str;

	if (int(theta*1000) >= int(0)) {

		dir_str.Form("yields." + particle_type + "_" + q2_eps_str + "_+%4.4i." + data_type, int(round(theta*1000)));	

	} else {

		dir_str.Form("yields." + particle_type + "_" + q2_eps_str + "_%4.4i." + data_type, int(round(theta*1000)));	

	}

	return dir_str;

}
