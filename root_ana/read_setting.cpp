
#include <iostream>
#include <fstream>
 
// #include <stdio.h>   
// #include <stdlib.h>  
 
#include "read_setting.h"

#include "general_utility.h"



/*--------------------------------------------------*/
/*--------------------------------------------------*/

using namespace std;



ReadFile::ReadFile() {

	static int tsttst = 0;

	if (tsttst == 0) {

		cout << " /*--------------------------------------------------*/" << endl;
		cout << "             loading the scalar information            " << endl << endl;
//		cout << " /*--------------------------------------------------*/" << endl;
		tsttst++;
	}


	Read_init();
	Setting_file_loading();

	Kin_Pro_Array_Load();
	// eff_pro_array_load();
	Cen_Pro_Array_Load();

	Calculate_t_Width();



}


// 
// /*--------------------------------------------------*/
// 
ReadFile::ReadFile(TString eff_file_name_tmp, TString off_file_name_tmp) {

	cout << "12312312 " << off_file_name_tmp << endl;

	static int tsttst = 0;

	if (tsttst == 0) {

		cout << " /*--------------------------------------------------*/" << endl;
		cout << "             loading the scalar information            " << endl << endl;
//		cout << " /*--------------------------------------------------*/" << endl;
		tsttst++;
	}


	cout << eff_file_name_tmp << endl;

	eff_file_name = eff_file_name_tmp;
	off_file_name = off_file_name_tmp;

	Read_init();
	Setting_file_loading_name();

	Kin_Pro_Array_Load();

	// eff_pro_array_load();

	Cen_Pro_Array_Load();

	// Sim_Pro_Array_Load();

	Calculate_t_Width();



}



// /*--------------------------------------------------*/
// 
ReadFile::ReadFile(TString sim_file_name_tmp) {

	///*--------------------------------------------------*/
	// Read Simulated data information file
	/*--------------------------------------------------*/

	sim_file_name = sim_file_name_tmp;

	if( sim_file_name.Contains("omega") || sim_file_name.Contains("Omega") 
		|| sim_file_name.Contains("xphsp") || sim_file_name.Contains("rho") 
		|| sim_file_name.Contains("eta") ) {

		Sim_Pro_Array_Omega_Load();

	} else {

		Sim_Pro_Array_Load();

	}

//	cout << "Angle Check "<< sim_pro.hms_theta[2] << endl;
//	exit(0);

}




 

void ReadFile::Read_init() {


	/*--------------------------------------------------*/
	/// Loading the offset settings

	setting_ntp = new TNtuple("setting","setting", "polset:Q2set:epsset:thpqset:tmnset:tmxset:NBtset:kset");
	eff_ntp = new TNtuple("eff","eff", "run_num:Q2:ebeam:charge:charge_err:tot_eff:tot_eff_err:epsset:thpqset");
	center_ntp = new TNtuple("center","center", "center_run_num:center_mean:center_sigma");

	// sim_ntp = new TNtuple("sim","sim", "q2:ebeam:normfac:event_num");

	dummy_tar = false;

//	string temp_str;
//	temp_str = list_file;

// 	if (temp_str.find("dummy") != string::npos) {
// 		dummy_tar = true;
// 	}
	

}

 




/*--------------------------------------------------*/
/*--------------------------------------------------*/

void ReadFile::Setting_file_loading() {

	/*--------------------------------------------------*/	
	/// Knematics setting loading 
	// setting_ntp->ReadFile("list.settings.fpi2");
	setting_ntp->ReadFile("list.settings.heep");

	///*--------------------------------------------------*/
	// Read center file
	/*--------------------------------------------------*/
	// center_ntp->ReadFile("fit_piplus/cointime_pl.dat");
	center_ntp->ReadFile("offset.dat");


}


void ReadFile::Setting_file_loading_name() {

	/*--------------------------------------------------*/	
	/// Knematics setting loading 
	// setting_ntp->ReadFile("list.settings.fpi2");
	setting_ntp->ReadFile(eff_file_name);

	///*--------------------------------------------------*/
	// Read center file
	/*--------------------------------------------------*/
	// center_ntp->ReadFile("fit_piplus/cointime_pl.dat");
	center_ntp->ReadFile(off_file_name);


	///*--------------------------------------------------*/
	// Read Simulated data information file
	/*--------------------------------------------------*/
	// sim_ntp->ReadFile(sim_file_name);



}





void ReadFile:: Eff_file_loading(TString list_file_in) {


	/*--------------------------------------------------*/	
	/// Directory
	// list_dir 			= "lists/";  

	/*--------------------------------------------------*/
	/// Input file name
	// list_file = "list.dummy_245_27_plus";


// 	string list_file11;
// 	list_file11 = list_file_in;



	/// Check file existed and empty
 	if (File_exist(list_file_in) && !File_empty(list_file_in)) {

		eff_ntp->ReadFile(list_file_in);
		Eff_Pro_Array_Load();


		if (eff_ntp->ReadFile(list_file_in)==0) {
		    // cout << "Bare "<< eff_ntp->ReadFile(list_file_in) << "   +1 " << eff_ntp->ReadFile(list_file_in)+1 << endl;
			num_runs = 0;
		}


 	} else{

		num_runs = 0;

		//kin_pro = get();

	}



	/*--------------------------------------------------*/
	/// Efficiency loading 
 	



}


/*--------------------------------------------------*/
/*--------------------------------------------------*/


void ReadFile::Kin_Pro_Array_Load() {

	kin_pro.polset     = Return_array_D("polset", setting_ntp);
	kin_pro.Q2set      = Return_array_D("Q2set",  setting_ntp);
	kin_pro.epsset     = Return_array_D("epsset", setting_ntp);
	                                               
	kin_pro.thpqset    = Return_array_D("thpqset", setting_ntp);
	kin_pro.tmnset     = Return_array_D("tmnset",  setting_ntp);
	kin_pro.tmxset     = Return_array_D("tmxset",  setting_ntp);
	                                               
	Double_t* NBtset_d = Return_array_D("NBtset", setting_ntp);
	kin_pro.kset       = Return_array_D("kset",   setting_ntp);	
	
//	set_runs = SizeOfArray(NBtset_d);
	set_runs = SizeOfArray(kin_pro.kset);




//	cout << "******************* REad Check :: " << set_runs << endl;

	kin_pro.NBtset = Array_D_to_I (NBtset_d, set_runs);

	
	// cout << " Read: "<< int(kin_pro.Q2set[0]*100) << endl;
	// cout << " Read: "<< kin_pro.polset[0] << endl;


}



void ReadFile::Eff_Pro_Array_Load() {


	
	Double_t* run_num_d  = Return_array_D("run_num",         eff_ntp);

	eff_pro1.Q2          = Return_array_D("Q2",              eff_ntp);
	eff_pro1.ebeam       = Return_array_D("ebeam",           eff_ntp);

	eff_pro1.charge      = Return_array_D("charge",          eff_ntp);
	eff_pro1.charge_err  = Return_array_D("charge_err/100",  eff_ntp);
	       
	eff_pro1.tot_eff     = Return_array_D("tot_eff",         eff_ntp);
	eff_pro1.tot_eff_err = Return_array_D("tot_eff_err/100", eff_ntp);
	       
	eff_pro1.epsset      = Return_array_D("epsset",          eff_ntp);
	eff_pro1.thpqset     = Return_array_D("thpqset",    	 eff_ntp);



	num_runs = SizeOfArray(run_num_d);

	eff_pro1.run_num = Array_D_to_I (run_num_d, num_runs);


//	cout << "READ !!!!!!!!! : " << run_num_d[0] << endl;


//	cout << "size: " <<  SizeOfArray(run_num_d) << endl;

	

	
//	cout << " Read: "<< kin_pro.polset[0] << endl;
//	cout << " Read: "<< kin_pro.thpqset[0] << endl;

}



void ReadFile::Cen_Pro_Array_Load() {

	Double_t* center_run_num_d = Return_array_D("center_run_num", center_ntp);
 	cen_pro.center_mean      = Return_array_D("center_mean",    center_ntp);
 	cen_pro.center_sigma     = Return_array_D("center_sigma",   center_ntp);


	cen_runs = SizeOfArray(center_run_num_d);
                                           
	cen_pro.center_run_num = Array_D_to_I (center_run_num_d, cen_runs);
//	cout << "center size: " << SizeOfArray(center_run_num_d) << endl; 

}





void ReadFile::Sim_Pro_Array_Load() {

	sim_ntp = new TNtuple("sim","sim", "ebeam:q2:normfac:event_num");
	sim_ntp->ReadFile(sim_file_name);


	sim_pro.ebeam 		     = Return_array_D("ebeam",     sim_ntp);
	sim_pro.q2_setting 		 = Return_array_D("q2",        sim_ntp);
 	sim_pro.normfac    		 = Return_array_D("normfac",   sim_ntp);
	Double_t* event_num_temp = Return_array_D("event_num", sim_ntp);


// 	cout << sim_pro.q2_setting[0] << endl;
// 	cout << sim_pro.q2_setting[1] << endl;
// 	cout << sim_pro.q2_setting[2] << endl;
// 	cout << sim_pro.q2_setting[3] << endl;


	Int_t num_tem = SizeOfArray(event_num_temp);
                                           
	sim_pro.event_num = Array_D_to_I(event_num_temp, num_tem);
//	cout << "center size: " << SizeOfArray(center_run_num_d) << endl; 

}

/*--------------------------------------------------*/
/*--------------------------------------------------*/

void ReadFile::Sim_Pro_Array_Omega_Load() {

	sim_ntp = new TNtuple("sim","sim", "q2:epsilon:hms_theta:normfac:event_num");
	sim_ntp->ReadFile(sim_file_name);


	sim_pro.q2_setting 		 = Return_array_D("q2",        sim_ntp);
	sim_pro.epsilon 		 = Return_array_D("epsilon",   sim_ntp);
	sim_pro.hms_theta 		 = Return_array_D("hms_theta", sim_ntp);
 	sim_pro.normfac    		 = Return_array_D("normfac",   sim_ntp);
	Double_t* event_num_temp = Return_array_D("event_num", sim_ntp);


// 	cout << sim_pro.q2_setting[0] << endl;
// 	cout << sim_pro.q2_setting[1] << endl;
// 	cout << sim_pro.q2_setting[2] << endl;
// 	cout << sim_pro.q2_setting[3] << endl;


	Int_t num_tem = SizeOfArray(event_num_temp);
                                           
	sim_pro.event_num = Array_D_to_I(event_num_temp, num_tem);
//	cout << "center size: " << SizeOfArray(center_run_num_d) << endl; 

}




/*--------------------------------------------------*/
/*--------------------------------------------------*/



Double_t* ReadFile::Return_array_D(TString variable_name, TTree* target_tree) {

	target_tree->Draw(variable_name, "", "goff");

	TArrayD* array_tmp = new TArrayD(target_tree->GetSelectedRows(), target_tree->GetV1());
	Double_t* array = array_tmp->GetArray();

	return array;
	delete array_tmp;

	/*--------------------------------------------------*/
	/// Create object and return address
	// Double_t array = *array_tmp->GetArray();
	// return &array;

}



/*--------------------------------------------------*/
/*--------------------------------------------------*/

void ReadFile::Calculate_t_Width() {

	t_bin_set = kin_pro.NBtset[0];
	t_min   = kin_pro.tmnset[0]; 
	t_max   = kin_pro.tmxset[0]; 
	t_width = ( t_max - t_min ) / t_bin_set; 

}



/*--------------------------------------------------*/
/*--------------------------------------------------*/



// 
// 
// Int_t* ReadFile::Return_array_I(TString variable_name, TTree* target_tree) {
// 
// 	target_tree->Draw(variable_name+"/I", "", "goff");
// 
// 	TArrayI* array_tmp = new TArrayI(target_tree->GetSelectedRows(), target_tree->GetV1());
// 	Int_t* array = array_tmp->GetArray();
// 
// 	return array;
// 	delete array_tmp;
// 
// 	/*--------------------------------------------------*/
// 	/// Create object and return address
// 	// Double_t array = *array_tmp->GetArray();
// 	// return &array;
// 
// }


/*--------------------------------------------------*/
/*--------------------------------------------------*/
// 
// Double_t SizeOfArray(Double_t* tar_arr) {
// 
// 	int size_of_array = 0;
// 
// 	while( int (*tar_arr) > 0 ) {
// 	    tar_arr++;
// 		size_of_array++;
// 	}
// 
// 	return size_of_array;
// }
// 



/*--------------------------------------------------*/
/*--------------------------------------------------*/

void ReadFile::Ntuple_Reset() {

	eff_ntp->Reset();

}






/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Check if the file exists

inline bool ReadFile::File_exist (TString name) {

	string name_str;
	name_str = name;

    return ( access( name_str.c_str(), F_OK ) != -1 );
}





/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Check if the file is empty

inline bool ReadFile::File_empty (TString name) {

	string name_str;
	name_str = name;

	ifstream infile;

	infile.open(name_str.c_str(), ifstream::in);
 
    return infile.peek() == ifstream::traits_type::eof();

}





ReadFile::~ReadFile() {


}

