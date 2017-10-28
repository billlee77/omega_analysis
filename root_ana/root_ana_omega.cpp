/*--------------------------------------------------*/
// Program for analysing the omega data 
// Author: wl
// Email: wenliang.billlee@gmail.com

#include <iostream>
#include "root_ana_omega.h"

//#include "cut.h"

using namespace std;

int main() {

	sim_ana_name = "omega";
	Initialization();
 
  	Dummy_Load();
  	Dummy_Setup();
  	Dummy_Analysis();
 
   
       
      	Target_Load();
      	Target_Setup();
      	Target_Analysis();
  
      cout << " -------------------Target Analysis Done----------------- " << endl;        
   
         Sim_Load();
         Sim_Analysis();
         Cleanup();
     
         sim_ana_name = "rho";
         Initialization();
         Sim_Load();
         Sim_Analysis();
         Cleanup();
         
         sim_ana_name = "xphsp";
         Initialization();
         Sim_Load();
         Sim_Analysis();
         Cleanup();
         
         sim_ana_name = "eta";
         Initialization();
         Sim_Load();
         Sim_Analysis();
         Cleanup();
         
         sim_ana_name = "eta_prime";
         Initialization();
         Sim_Load();
         Sim_Analysis();
         Cleanup();
     
 // 	cout << endl << "/*--------------------------------------------------*/ " << endl;
 // 	cout         << "                   Analysis Starts                     " << endl << endl;
 // //	cout         << "/*--------------------------------------------------*/ " << endl << endl;
 
 	cout << endl << "/*--------------------------------------------------*/ "  << endl;
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

 	ana_dummy 	       = new Analysis_omega*[set_runs];
 	ana_target         = new Analysis_omega*[set_runs];
 	ana_sim            = new Analysis_omega*[set_runs];
// 
//  	ana_dummy 	       = new Analysis*[set_runs];
//  	ana_target         = new Analysis*[set_runs];


 	eff_setting_dummy  = new ReadFile::eff_pro[set_runs];
 	eff_setting_target = new ReadFile::eff_pro[set_runs];


	cout << "Set Runs !!!  : " << set_runs << endl;

	list_dir 			= "lists/";


	ofstream u_bin_int;
	u_bin_int.open(fo->Get_Dir() + "u_bin_interval", ios::out);

//	u_bin_int.open(file_out_name(0, file_out_name.First("\/"))+ "/u_bin_interval", ostream::out);
//	cout << fo->Get_Dir() + "u_bin_interval" << endl;
//	u_bin_int << " aaaaaaaaaaaaaaaaa " << endl;

	u_bin_int.close();


//	exit(0);

	
}




/*--------------------------------------------------*/
/*--------------------------------------------------*/

TString Dummy_Name(int num_itt) {

	TString dir_str;

	if (int(kin_setting.thpqset[num_itt]*1000) >= int(0)) {

		dir_str.Form("yields.omega_%i_%i_+%4.4i.dummy", int(round(kin_setting.Q2set[num_itt]*100)), int(round(kin_setting.epsset[num_itt]*100)), int(kin_setting.thpqset[num_itt]*1000));	

	} else {

		dir_str.Form("yields.omega_%i_%i_%4.4i.dummy", int(round(kin_setting.Q2set[num_itt]*100)), int(round(kin_setting.epsset[num_itt]*100)), int(kin_setting.thpqset[num_itt]*1000));	

	}

	return dir_str;

}


/*--------------------------------------------------*/
/*--------------------------------------------------*/

TString Target_Name(int num_itt) {

 	TString dir_str;

	if (int(kin_setting.thpqset[num_itt]*1000) >=  int(0)) {

		dir_str.Form("yields.omega_%i_%i_+%4.4i.target", int(round(kin_setting.Q2set[num_itt]*100)), int(round(kin_setting.epsset[num_itt]*100)), int(kin_setting.thpqset[num_itt]*1000));	

	} else {

		dir_str.Form("yields.omega_%i_%i_%4.4i.target", int(round(kin_setting.Q2set[num_itt]*100)), int(round(kin_setting.epsset[num_itt]*100)), int(kin_setting.thpqset[num_itt]*1000));	

	}

	return dir_str;

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

TString Target_Name(TString data_type, int num_itt) {

// 	TString dir_str;
// 
// 	if (int(kin_setting.thpqset[num_itt]*1000) >= int(0)) {
// 
// 		dir_str.Form("yields.omega_%i_%i_+%4.4i." + data_type, int(round(kin_setting.Q2set[num_itt]*100)), int(round(kin_setting.epsset[num_itt]*100)), int(kin_setting.thpqset[num_itt]*1000));	
// 
// 	} else {
// 
// 		dir_str.Form("yields.omega_%i_%i_%4.4i." + data_type, int(round(kin_setting.Q2set[num_itt]*100)), int(round(kin_setting.epsset[num_itt]*100)), int(kin_setting.thpqset[num_itt]*1000));	
// 
// 	}
// 
// 	return dir_str; 

	TString dir_str_1;

	dir_str_1 = Standard_Name_Str("omega", data_type, kin_setting.Q2set[num_itt], kin_setting.epsset[num_itt], kin_setting.thpqset[num_itt]);
	
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
/// Return Setting Number for Dummy Runs
int Return_Setting_dummy(int list_itt, int run_itt, int set_itt) {

	Int_t iii = 0;

	
 	TString run_setting_name; 

// 	if (int(eff_setting_dummy[set_itt].thpqset[run_itt]*1000) >= int(0)) {
// 
// 		run_setting_name.Form( "yields.omega_" + Q2_vec_mod[list_itt] + "_+%4.4i.dummy", int(eff_setting_dummy[set_itt].thpqset[run_itt]*1000));
// 
// 	} else {
// 
// 		run_setting_name.Form( "yields.omega_" + Q2_vec_mod[list_itt] + "_%4.4i.dummy", int(eff_setting_dummy[set_itt].thpqset[run_itt]*1000));
// 
// 	}

	cout << "****!!!!!!~~~~~~~~     ::      " << setting_vec_dummy[set_itt] << "   aaa   " 
		 << run_setting_name  << "     " << eff_setting_dummy[set_itt].thpqset[run_itt] << endl;

	run_setting_name = Standard_Name_Str("omega", "dummy", Q2_vec_mod[list_itt], eff_setting_dummy[set_itt].thpqset[run_itt]);
	
	if (setting_vec_dummy[set_itt] == run_setting_name) {

	// 	cout << "AAAA" << setting_vec_dummy[set_itt] << "  " << list_itt << "   " << run_itt  << "  " << set_itt << endl;
		cout << setting_vec_dummy[set_itt] << "  " << list_itt << "      " << run_itt  
			 <<  "     " <<  set_itt << endl;

		setting_num = set_itt;
	}

//	cout << "  rr  " << run_setting_name << endl;

	return setting_num;

}


/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Return Setting Number for Target Runs

int Return_Setting_target(int list_itt, int run_itt, int set_itt) {

 	Int_t iii = 0;

 	TString run_setting_name; 

// 	if (int(eff_setting_target[set_itt].thpqset[run_itt]*1000) >= int(0)) {
// 
// 		run_setting_name.Form( "yields.omega_" + Q2_vec_mod[list_itt] + "_+%4.4i.target", int(eff_setting_target[set_itt].thpqset[run_itt]*1000));
// 
// 	} else {
// 
// 		run_setting_name.Form( "yields.omega_" + Q2_vec_mod[list_itt] + "_%4.4i.target", int(eff_setting_target[set_itt].thpqset[run_itt]*1000));
// 
// 	}

	run_setting_name = Standard_Name_Str("omega", "target", Q2_vec_mod[list_itt], eff_setting_target[set_itt].thpqset[run_itt]);
	
	if (setting_vec_target[set_itt] == run_setting_name) {
	// 	cout << "AAAA" << setting_vec_target[set_itt] << "  " << list_itt << "   " << run_itt  << "  " << set_itt << endl;
		cout << setting_vec_target[set_itt] << "  " << list_itt << "   " << run_itt  << "  " << set_itt << endl;
		setting_num = set_itt;
	}

	return setting_num;

}



// /*--------------------------------------------------*/
// /*--------------------------------------------------*/
// 
// int Return_Setting(int list_itt, int run_itt, int set_itt) {
// 
// 	Int_t iii = 0;
// 	
//  	TString run_setting_name; 
// 
// 	if (int(eff_setting[set_itt].thpqset[run_itt]*1000) == 0) {
// 
// 		run_setting_name.Form( "yields.omega_" + Q2_vec[list_itt] + "_0000.dummy", int(eff_setting[set_itt].thpqset[run_itt]*1000)); 
// 
// 	} else {
// 
// 		run_setting_name.Form( "yields.omega_" + Q2_vec[list_itt] + "_%4.4i.dummy", int(eff_setting[set_itt].thpqset[run_itt]*1000));
// 
// 	}
// 
// 
// 	if (setting_vec[set_itt] == run_setting_name) {
// 		cout << setting_vec[set_itt] << "  " << list_itt << "   " << run_itt  << "  " << set_itt << endl;
// 		setting_num = set_itt;
// 	}
// 
// 	return setting_num;
// 
// }
// 


/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Cout_out_lists() {

	cout << endl;
	cout << endl;

	/*--------------------------------------------------*/
	/// cout check on the lists 

 	for (int a=0; a < Q2_vec.size(); a++) {
		cout << "Q2 value:  " << a << "    " << Q2_vec[a] << endl;
	}

	cout << endl;

 	for (int a=0; a < Q2_vec_mod.size(); a++) {
		cout << "Q2 value mod:  " << a << "    " << Q2_vec_mod[a] << endl;
	}

	cout << endl;

 	for (int a=0; a < lists_vec.size(); a++) {
		cout << "Run list:  " << a << "    " << lists_vec[a] << endl;
	}

	cout << endl;

 	for (int a=0; a < lists_vec_target.size(); a++) {
		cout << "Target run list:  " << a << "    " << lists_vec_target[a] << endl;
	}

	cout << endl;

 	for (int a=0; a < lists_vec_dummy.size(); a++) {
		cout << "Dummy run list:  " << a << "    " << lists_vec_dummy[a] << endl;
	}




 	for (int a=0; a < setting_vec.size(); a++) {
		cout << "Setting list:  " << a << "    " << setting_vec[a] << endl;
	}

	cout << endl;

 	for (int a=0; a < setting_vec_target.size(); a++) {
		cout << "Target Setting list:  " << a << "    " << setting_vec_target[a] << endl;
	}

	cout << endl;

 	for (int a=0; a < setting_vec_dummy.size(); a++) {
		cout << "Dummy Setting list:  " << a << "    " << setting_vec_dummy[a] << endl;
	}

	cout << endl;

 	for (int a=0; a < run_list_vec.size(); a++) {
		cout << "Run list vec:  " << a << "    " << run_list_vec[a] << endl;
	}

	cout << endl;


	for (int a=0; a < run_list_vec_target.size(); a++) {
		cout << "Run list target vec:  " << a << "    " << run_list_vec_target[a] << endl;
	}

	cout << endl;

	for (int a=0; a < run_list_vec_dummy.size(); a++) {
		cout << "Run list dummy vec:  " << a << "    " << run_list_vec_dummy[a] << endl;
	}

	cout << endl;

}






/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Loading Functions

/*--------------------------------------------------*/
/// Loading efficiency fill the untility arrays etc for the Dummy Runs

void Dummy_Load() {

	Int_t ii = 0;

	for (Int_t i = 0; i < set_runs; i++) {

		if (kin_setting.polset[i] > 0) {

//			output_file_dummy[ii] = fo->Create_Out_File(Dummy_Name(i));
			output_file_dummy[ii] = fo->Create_Out_File(Target_Name("dummy", i));

//			list_name = list_dir + "list.dummy_" + List_Name(i) + "_plus";
			list_name = list_dir + "list.dummy_" + List_Name(i) + "_omega";

			cout << list_name << endl;

//			exit(0);


			lists_vec_itt = find (lists_vec.begin(), lists_vec.end(), list_name);

			if (!rf->File_exist(list_name) && rf->File_empty(list_name)) {

				cout << "List file for Setting " << list_name << " doesn't exist !! " << endl;

			} else {


				if(lists_vec_itt == lists_vec.end()) {
	
					lists_vec.push_back(list_name);
					rf->Eff_file_loading(list_name);
	
					num_runs = rf->Get_Num_Runs();
		
					/// Mod
	 				run_list_vec_dummy.push_back(num_runs);
	
				}
	
				Q2_vec_mod_itt = find(Q2_vec_mod.begin(), Q2_vec_mod.end(), List_Name(i));
	
				if(Q2_vec_mod_itt == Q2_vec_mod.end()) {
					Q2_vec_mod.push_back(List_Name(i));
				}
	
				eff_setting_dummy[ii] = rf->eff_pro1;
	
	//			ana_dummy[ii] = new Analysis(rf->eff_pro1);
				ana_dummy[ii] = new Analysis_omega(rf->eff_pro1);
	
//				setting_vec_dummy.push_back(Dummy_Name(i));	
	// 			cout  << Dummy_Name(i) << endl;

				setting_vec_dummy.push_back(Target_Name("dummy", i));	
//	 			cout  << "dummy " << Target_Name("dummy", i) << endl;
	
//				exit(0);

				rf->Ntuple_Reset();
	
				ii++;
	
			}
				
		}

	}


//	exit(0);

	setting_run_list_dummy.resize(ii);

}

//*--------------------------------------------------*/
/// Loading efficiency fill the untility arrays etc for the Target Runs

void Target_Load() {

	Int_t iiii = 0;

	for (Int_t i = 0; i < set_runs; i++) {

		if (kin_setting.polset[i] > 0) {

			// output_file_target[iiii] = fo->Create_Out_File(Target_Name(i));
			output_file_target[iiii] = fo->Create_Out_File(Target_Name("target", i));

			// list_name = list_dir + "list." + List_Name(i) + "_plus";
			list_name = list_dir + "list." + List_Name(i) + "_omega";

			// cout << list_name << endl;

			if (!rf->File_exist(list_name) && rf->File_empty(list_name)) {

				cout << "List file for Setting " << list_name << " doesn't exist !! " << endl;

			} else {

				lists_vec_itt = find (lists_vec.begin(), lists_vec.end(), list_name);

				if(lists_vec_itt == lists_vec.end()) {

 					lists_vec.push_back(list_name);

					cout << "!!!!!!!  " << list_name << endl;

 					rf->Eff_file_loading(list_name);

					num_runs = rf->Get_Num_Runs();

					/// Mod
 					run_list_vec_target.push_back(num_runs);

				}

				Q2_vec_mod_itt = find(Q2_vec_mod.begin(), Q2_vec_mod.end(), List_Name(i));

				if(Q2_vec_mod_itt == Q2_vec_mod.end()) {
					Q2_vec_mod.push_back(List_Name(i));
				}

				eff_setting_target[iiii] = rf->eff_pro1;

				ana_target[iiii] = new Analysis_omega(rf->eff_pro1);
//				ana_target[iiii] = new Analysis(rf->eff_pro1);

// 				setting_vec_target.push_back(Target_Name(i));
// 				cout << Target_Name(i) << endl;


				setting_vec_target.push_back(Target_Name("target", i));
	 			cout << "target " << Target_Name("target", i) << endl;

				

				rf->Ntuple_Reset();

				iiii++;

			}

		}
		
	
	}

	setting_run_list_target.resize(iiii);

//	exit(0);

}




/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Setting-up Functions

void Dummy_Setup() {

	cout << endl;
	cout << "/*--------------------------------------------------*/" << endl;
	cout << "             Setting up the dummy runs                " << endl << endl;
//	cout << "/*--------------------------------------------------*/" << endl;

	/*--------------------------------------------------*/
	/// Setting up the setting_run_list for the Dummy runs

	for (int i=0;  i < run_list_vec_dummy.size(); i++) {

	    for (int ii = 0; ii < run_list_vec_dummy[i]; ii++) { 

			int setting_num_tmp;

			for (int iii = 0; iii < setting_vec_dummy.size(); iii++) {
				setting_num_tmp = Return_Setting_dummy(i, ii, iii);
			}

			cout << "Dummy Seeeeeettting num tmp    " << setting_num_tmp << endl;

			setting_run_list_dummy[setting_num_tmp].push_back(ii);

   	 	}	
	}
}

/*--------------------------------------------------*/
/// Setting up the setting_run_list for the Target runs

void Target_Setup() {

	cout << endl;
	cout << "/*--------------------------------------------------*/" << endl;
	cout << "             Setting up the target runs                " << endl << endl;



	for (int i=0;  i < run_list_vec_target.size(); i++) {

		cout << "Is this working ? " << i << "     " << run_list_vec_target[i] << endl;

// 		if(run_list_vec_target.size()==0)
// 			break;


	    for (int ii = 0; ii < run_list_vec_target[i]; ii++) { 

			cout << i << "    " <<  ii << endl;

			int setting_num_tmp;

			for (int iii = 0; iii < setting_vec_target.size(); iii++) {
				setting_num_tmp = Return_Setting_target(i, ii, iii);

			}

			// cout << "Target Seeeeeettting num tmp    " << setting_num_tmp << endl;

			setting_run_list_target[setting_num_tmp].push_back(ii);

   	 	}	

	}

}




/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Analyzing Functions

/*--------------------------------------------------*/
/// Analyzing data run by run for the Dummy Runs

void Dummy_Analysis() {

	for(int i = 0; i < setting_run_list_dummy.size(); i++) {

		if (setting_run_list_dummy[i].size() != 0 ) {

 			ana_dummy[i]->file_out_ana =  output_file_dummy[i];
 
 			for (int ii = 0; ii < setting_run_list_dummy[i].size(); ii++) {
 
// 				cout << "run " << setting_run_list_dummy[i][ii] << endl;
 
				ana_dummy[i]->is_run_dummy = true;

				ana_dummy[i]->Omega_anna(setting_run_list_dummy[i][ii]);
				
 				cout << "~  " << setting_run_list_dummy.size() << "    " << setting_run_list_dummy[i].size() 
 					 << "     " << i << "    " << ii << "     " <<  setting_run_list_dummy[i][ii] << endl;
 
 			}
 
// 			cout << "******************" << endl;
 			ana_dummy[i]->Yield_Out();

			delete ana_dummy[i];

		}


	}

}

/*--------------------------------------------------*/
/// Analyzing data run by run for the Target Runs

void Target_Analysis() {

	for(int i = 0; i < setting_run_list_target.size(); i++) {

		if (setting_run_list_target[i].size() != 0 ) {

			cout << "kissy kissy" << endl;

			ana_target[i]->file_out_ana =  output_file_target[i];

			for (int ii = 0; ii < setting_run_list_target[i].size(); ii++) {
	
//				cout << "run " << setting_run_list_target[i][ii] << endl;

				ana_target[i]->is_run_dummy = false;
//				ana_target[i]->Run_by_Run_Analysis(setting_run_list_target[i][ii]);
				ana_target[i]->Omega_anna(setting_run_list_target[i][ii]);

			}

			cout << "******************" << endl;
 			ana_target[i]->Yield_Out();

			delete ana_target[i];

		}

	}

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
				
			ana_sim[i] = new Analysis_omega(sim_ana_name);

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



// 		cout << "Are we here????  ++++++++++++++++ " << endl;
// 		exit(0);

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
                      
    delete ana_dummy; 	      
    delete ana_target; 
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
