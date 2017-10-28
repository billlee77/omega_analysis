#include <iostream>

#include "analysis_omega.h"
//#include "cut.h"

using namespace std;

/*--------------------------------------------------*/

Analysis_omega::Analysis_omega() {

}

/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Only for experimental data

Analysis_omega::Analysis_omega(ReadFile::efficiency_profile eff_struc) {
	
	eff_file = "list.settings.omega";
	off_file = "offset.dat";
// 
// 	cout << off_file << "sdafasfs" << endl;

	eff_ana = eff_struc;
	Init();


	graph_yield_check = new TGraphErrors();
	mm_peak_check = new TGraphErrors();

	yield_setting     = 0.0;
	yield_setting_err = 0.0;

	chain = new TChain();
	list = new TList();

	coin_beta_check_offset_setting = new TH2F("coin_beta_check_offset", "coin_beta_check", 200, -15, 15,  200, 0, 2);

	diamond_setting = new TMultiGraph();

	st_itt = 0;

}



/*--------------------------------------------------*/
/// Only for simulation data

Analysis_omega::Analysis_omega(TString sim_data_type) {

	sim_ana_name = sim_data_type;
	weight_app = "";

}



/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Analysis for Heep 
void Analysis_omega::Omega_anna(Int_t run_itt) {
 
	is_sim_data = false;

 	cout << "Analyzing Omega data !" << endl;
 
 	cout << " sdaf " << run_itt << endl;
 
  	run_tree = Create_File();
  
  	Para_Run_Def(run_itt);

	Para_Init();

 	Load_Data_Tree();

 	out_dir->cd();

  	coin_center = Analysis::Get_Coin_Center();

 	Define_Cuts();

	beta_cut                   = 0.1;
	gradi                      = 1./7.;	
	global_centered_hsbeta_cut = -0.1;

	
	rand_early_l = 9.425;
	rand_early_r = 3.275;
	rand_late_l  = 5.25;
	rand_late_r  = 13.45;  


	cointime_cut_l = 1.1;  
	cointime_cut_r = 0.95;  

	top_limit = 0.7;
	bot_limit = 0.6;

	random_real_ratio = (fabs(rand_late_l - rand_late_r) + fabs(rand_early_l-rand_early_r))/(cointime_cut_l + cointime_cut_r);

	// Target thickness difference between the dummy target and Loop 1 hydrogen target
	dummy_correction = 0.142;



	hsbeta_center = GetBetaCenter();


// 	/*--------------------------------------------------*/
// 	Plot to check: missmass, W, t, EXP/SIM, Q2, th_pq,
// 				   phi_pq, PmPar, PmPer, PmOop, hsdelta, hsytar,
// 					   hsyptar, hsxptar, ssdelta, ssytar, ssyptar, ssxptar

	Calculate_Yield_Direct();

	Check_Plot_Out("missmass");
	Check_Plot_Out("Em");
	Check_Plot_Out("Pm");
	Check_Plot_Out("W");

	Check_Plot_Out("hsdelta");
	Check_Plot_Out("t");
	Check_Plot_Out("Q2");
	Check_Plot_Out("th_pq");
	Check_Plot_Out("phi_pq");
	Check_Plot_Out("Pmpar");
	Check_Plot_Out("Pmper");
	Check_Plot_Out("Pmoop");
	Check_Plot_Out("hsytar");
	Check_Plot_Out("hsyptar");
	Check_Plot_Out("hsxptar");
	Check_Plot_Out("ssdelta");
	Check_Plot_Out("ssytar");
	Check_Plot_Out("ssyptar");
	Check_Plot_Out("ssxptar");
	Check_Plot_Out("cointime");
 
// 	Check_Plot_Out("hsbeta", "cointime");


	Check_MissingMass_Peak();

//	t_phi_check->Write("t_phi_check");

//	coin_beta_check->Write("coin_beta_check");

//	delete t_phi_check;
//	delete coin_beta_check;




// delete coin_beta_check_offset;

//	File_Clean();

//  data_tree_in->SetName("111");

//	Correct_Offset_Tree();

//
//
//
//
//
//
//


//	TH2F* coin_beta_check_offset = new TH2F("coin_beta_check_offset", "coin_beta_check", 150, -15, 15,  150, 0, 2);



	///*--------------------------------------------------*/
	/// Plot Run by Run hsbeta vs cointime

	TH2F* coin_beta_check_offset = new TH2F(*coin_beta_check_offset_setting);

	TString offset_str;

	if (coin_center < 0){
		offset_str.Form("hsbeta:cointime+%f", fabs(coin_center));
	} else { 
		offset_str.Form("hsbeta:cointime-%f", coin_center);
	} 

	data_tree_in->Draw(offset_str + ">>coin_beta_check_offset", all_coin_cut , "goff");

	coin_beta_check_offset->Write("coin_beta_check_offset");

	coin_beta_check_offset_setting->Add(coin_beta_check_offset);


	///*--------------------------------------------------*/
	/// Add the cointime corrected tree to the list

	Correct_Offset_Tree();

	Setting_Beta_Cointime_Check_Setting_Tree(data_tree_in_cl);

	list->Add(data_tree_in_cl);

 	delete data_tree_in;

 	file_in->Close();
 	delete file_in;

// 	chain->Add(Get_Current_Run_File()+"/h9500");

}






/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Analysis_omega::Para_Init() {


// 	q2_setting = kin_ana.Q2set[0];
// 	q2_setting = kin_ana.Q2set[0];
// 
// 	kset       = kin_ana.kset[0];
// 	t_min	   = kin_ana.tmnset[0];
// 	t_max      = kin_ana.tmxset[0];
 

//	u_bin_set = kin_ana.NBtset[0];
//	t_width   = ( t_max - t_min ) / u_bin_set;

	// phi_bin_num = 8;
	// phi_bin_num = 6;
// 	phi_bin_num = 4;
// 
// 	pstp = 360./phi_bin_num;
// 	pmn  = pstp/2.;
// 	pmx  = 360. - pstp/2.;
 
//	phi_bin_num = 6;


	if (is_sim_data == true) {
		q2_setting = sim_Q2;	
	} else {
		q2_setting = Get_Q2();		
	}

	cout << endl << "bbbbbbbbbbbbbbbbbbbb   " << is_sim_data << "    " << q2_setting << endl;


	///*--------------------------------------------------*/
	/// How many u bins and phi bins
	
	///*--------------------------------------------------*/
	/// 3 u bins 8 phi bins

	// phi_bin_num = 9;
	// phi_bin_num = 8;

//	phi_bin_num = 7;

//	u_bin_num = 4;

//	u_lower_limit= new Float_t[u_bin_num]; 
//	u_upper_limit= new Float_t[u_bin_num]; 
// 
//    ///*--------------------------------------------------*/
//    /// 1s u Binning
//    /// Data: 14 Nov 2016
//    /// Unequal u binning
//    /// stop at 0
//    //
//      	if( q2_setting == float(1.6)) {
//      
//      		u_lower_limit[0] = 0.0;
//      		u_upper_limit[0] = 0.10;
//      
//      		u_lower_limit[1] = 0.10;
//      		u_upper_limit[1] = 0.17;
//      
//      		u_lower_limit[2] = 0.17;
//      		u_upper_limit[2] = 0.32;
//      
//      	} else if ( q2_setting == float(2.45)) {
//      
//      		u_lower_limit[0] = 0.0; 
//      		u_upper_limit[0] = 0.19;
//      
//      		u_lower_limit[1] = 0.19;
//      		u_upper_limit[1] = 0.30;
//      
//      		u_lower_limit[2] = 0.30;
//      		u_upper_limit[2] = 0.50;
//      
//      	}
    


//   	if( q2_setting == float(1.6)) {
//   
//   		u_lower_limit[0] = -0.02;
//   		u_upper_limit[0] = 0.12;
//   
//   		u_lower_limit[1] = 0.08;
//   		u_upper_limit[1] = 0.19;
//   
//   		u_lower_limit[2] = 0.15;
//   		u_upper_limit[2] = 0.34;
//   
//   	} else if ( q2_setting == float(2.45)) {
//   
//   		u_lower_limit[0] = -0.02; 
//   		u_upper_limit[0] =  0.21;
//   
//   		u_lower_limit[1] = 0.17;
//   		u_upper_limit[1] = 0.32;
//   
//   		u_lower_limit[2] = 0.28;
//   		u_upper_limit[2] = 0.52;
//   
//   	}
// 



                                        
// 
//   	if( q2_setting == float(1.6)) {
//   
//   		u_lower_limit[0] = 0.0;
//   		u_upper_limit[0] = 0.12;
//   
//   		u_lower_limit[1] = 0.12;
//   		u_upper_limit[1] = 0.15;
//   
//   		u_lower_limit[2] = 0.15;
//   		u_upper_limit[2] = 0.32;
//   
//   	} else if ( q2_setting == float(2.45)) {
//   
//   		u_lower_limit[0] = 0.0; 
//   		u_upper_limit[0] = 0.24;
//   
//   		u_lower_limit[1] = 0.24;
//   		u_upper_limit[1] = 0.28;
//   
//   		u_lower_limit[2] = 0.28;
//   		u_upper_limit[2] = 0.50;
//   
//   	}
//  
// 
// 






// 	///*--------------------------------------------------*/
// 	/// 2 u bins 10 phi bins
// 	/// Data: 21 Feb 2017
// 
// 	phi_bin_num = 10;
// 	u_bin_num = 2;
// 
// 	u_lower_limit= new Float_t[u_bin_num]; 
// 	u_upper_limit= new Float_t[u_bin_num]; 
// 
//   	if( q2_setting == float(1.6)) {
//   
//   		u_lower_limit[0] = 0.0;
//   		u_upper_limit[0] = 0.14;
//   
//   		u_lower_limit[1] = 0.14;
//   		u_upper_limit[1] = 0.32;
// 
//   	} else if ( q2_setting == float(2.45)) {
//   
//   		u_lower_limit[0] = 0.0; 
//   		u_upper_limit[0] = 0.25;
//   
//   		u_lower_limit[1] = 0.25;
//   		u_upper_limit[1] = 0.50;
// 
//   	}
 


  	///*--------------------------------------------------*/
  	/// 4 u bins 7 phi bins
  	/// Data: 21 Feb 2017
  
  	phi_bin_num = 7;
  	u_bin_num = 4;
  
  	u_lower_limit= new Float_t[u_bin_num]; 
  	u_upper_limit= new Float_t[u_bin_num]; 
  
    	if( q2_setting == float(1.6)) {
    
    		u_lower_limit[0] = 0.0;
    		u_upper_limit[0] = 0.1;
    
    		u_lower_limit[1] = 0.08;
    		u_upper_limit[1] = 0.16;
    
    		u_lower_limit[2] = 0.14;
    		u_upper_limit[2] = 0.22;
 
    		u_lower_limit[3] = 0.20;
    		u_upper_limit[3] = 0.35;
 
 
    
    	} else if ( q2_setting == float(2.45)) {
    
    		u_lower_limit[0] = 0.0; 
    		u_upper_limit[0] = 0.20;
    
    		u_lower_limit[1] = 0.18;
    		u_upper_limit[1] = 0.26;
    
    		u_lower_limit[2] = 0.24;
    		u_upper_limit[2] = 0.32;
 
    		u_lower_limit[3] = 0.30;
    		u_upper_limit[3] = 0.55;
    
    	}
 


	/*--------------------------------------------------*/
	///

	phi_stp = 360./phi_bin_num;
	phi_min  = 0;
	phi_max  = 360;

	phi_offset = 0;

//	phi_offset = 25.7;

//	phi_offset = 22.5;
//	phi_offset = 10;


	TString file_out_name;

	file_out_name = file_out_ana->GetName();

	cout << file_out_name  << "   " << file_out_name.First("\/") << endl;
	cout << file_out_name(0, file_out_name.First("\/")) << endl;

	ofstream u_bin_int;
	u_bin_int.open(file_out_name(0, file_out_name.First("\/"))+ "/u_bin_interval", ios::app);


	static Int_t lo_q_entry = 0;
	static Int_t hi_q_entry = 0;


  	if( q2_setting == float(1.6)) {
  
		if (lo_q_entry == 0) {

			u_bin_int << q2_setting <<  "    " << u_bin_num  << "    " << phi_bin_num << endl;
			u_bin_int << "     " << u_lower_limit[0];

			for (Int_t itt=0; itt < u_bin_num; itt++) {
				u_bin_int <<  "    " << u_upper_limit[itt];
			}
			u_bin_int << endl;
				
			lo_q_entry++;

		} 
  
  	} else if ( q2_setting == float(2.45)) {
  

		if (hi_q_entry == 0) {

			u_bin_int << q2_setting <<  "    " << u_bin_num  << "    " << phi_bin_num << endl;
			u_bin_int << "     " << u_lower_limit[0];

			for (Int_t itt=0; itt < u_bin_num; itt++) {
				u_bin_int <<  "    " << u_upper_limit[itt];
			}
			u_bin_int << endl;
				
			hi_q_entry++;

		} 

  	}

 	
// 	u_bin_int << "111111111111111111111 " << endl;



	u_bin_int.close();
//	exit(0);

	


	// cout << Get_Q2() << "   "<< u_lower_limit[0] << "   " << u_upper_limit[0] << endl;

	if (st_itt == 0 ) {

		yield_vec.resize(u_bin_num*phi_bin_num);
		yield_err_vec.resize(u_bin_num*phi_bin_num);
		phi_vec.resize(u_bin_num*phi_bin_num);
		tb_vec.resize(u_bin_num*phi_bin_num);

		cout<< "1111111111 " << st_itt << endl;
		
		st_itt = 2;
		
	} 

}


/*--------------------------------------------------*/

/*--------------------------------------------------*/
/// Dufine cut

void Analysis_omega::Define_Cuts() {

    /// Define four point for low and high Q2
    

	Double_t xx, yy;

	diamond_cut_245 = new TCutG("diamond_cut_245",4);

	diamond_cut_245->SetVarX("Q2");
	diamond_cut_245->SetVarY("W");

  	diamond_cut_245->SetPoint(0, 1.842, 2.388);
   	diamond_cut_245->SetPoint(1, 2.346,  2.279);
   	diamond_cut_245->SetPoint(2, 3.13, 1.997);
   	diamond_cut_245->SetPoint(3, 2.54,  2.142);

   	diamond_cut_245->GetPoint(0, xx, yy);
   	diamond_cut_245->SetPoint(4, xx, yy);

//	cout << xx << "   " << yy << endl;

	diamond_cut_160 = new TCutG("diamond_cut_160",4);

	diamond_cut_160->SetVarX("Q2");
	diamond_cut_160->SetVarY("W");

  	diamond_cut_160->SetPoint(0, 1.17, 2.355 );
   	diamond_cut_160->SetPoint(1, 1.59, 2.262 );
   	diamond_cut_160->SetPoint(2, 2.08,  2.048);
   	diamond_cut_160->SetPoint(3, 1.62,  2.16);
   	diamond_cut_160->GetPoint(0, xx, yy);
   	diamond_cut_160->SetPoint(4, xx, yy);
 
//	cout << xx << "   " << yy << endl;			
//	cout << "Q2 Setting " << q2_setting << endl;

	/// Select which diamond cut you want based on the Q2 setting 

	cout << endl << "aaaaaaaaaaaaaa   " << q2_setting << endl;
	
    if (q2_setting == float(1.6)) {
		
		diamond_cut = (TCutG*) diamond_cut_160->Clone();
		cout << "Q2 Setting " << q2_setting << endl;

	} else if (q2_setting == float(2.45)) {

		diamond_cut = (TCutG*) diamond_cut_245->Clone();
		cout << "Q2 Setting " << q2_setting << endl;

	} 

//	cout << "aaaaaaaaaaaaaaaaaaaaaaa "<< endl;
	
//	diamond_cut->SetName("diamond_cut");
     

/*--------------------------------------------------*/
/// Calculate the 4 boundaries of the box
///
/// Line #1
	Double_t x_l1_p1, y_l1_p1;
	Double_t x_l1_p2, y_l1_p2;
	Double_t gradient_l1;
	Double_t x_l1_inter, y_l1_inter;
	

	diamond_cut->GetPoint(0, x_l1_p1, y_l1_p1);
	diamond_cut->GetPoint(1, x_l1_p2, y_l1_p2);
	
	gradient_l1 = (y_l1_p1 - y_l1_p2) / (x_l1_p1 - x_l1_p2) ;

	cout << "aaaaaaaaaaaaaa " << endl;



//	cout << x_1  << "  " << y_1 << "  " << x_2 << "  " << y_2 << endl;
	y_l1_inter = y_l1_p1 - x_l1_p1 * gradient_l1 ;
	x_l1_inter = x_l1_p1 - y_l1_p1 / gradient_l1 ;

 	cout << "dx: " << (x_l1_p1 - x_l1_p2) << "    dy: " << (y_l1_p1 - y_l1_p2) <<  "    g: "  << gradient_l1 << endl;
 	cout << "x intercept   " << x_l1_inter << "       y intercept    " << y_l1_inter << endl;
// 
// 	cout << x_1  << "  " << y_1 << "  " << x_2 << "  " << y_2 << endl;
// 	y_inter = y_2 - x_2 * gradient ;
// 	x_inter = x_2 - y_2 / gradient ;
// 
// 
// 	cout << "dx: " << (x_1 - x_2) << "    dy: " << (y_1 - y_2) <<  "    g: "  << gradient << endl;
// 	cout << "x intercept check  " << x_inter << "       y intercept    " << y_inter << endl;
 


/*--------------------------------------------------*/	
/// Line #2
	Double_t x_l2_p1, y_l2_p1;
	Double_t x_l2_p2, y_l2_p2;
	Double_t gradient_l2;
	Double_t x_l2_inter, y_l2_inter;
	
	diamond_cut->GetPoint(1, x_l2_p1, y_l2_p1);
	diamond_cut->GetPoint(2, x_l2_p2, y_l2_p2);
	
	gradient_l2 = (y_l2_p1 - y_l2_p2) / (x_l2_p1 - x_l2_p2) ;

//	cout << x_1  << "  " << y_1 << "  " << x_2 << "  " << y_2 << endl;
	y_l2_inter = y_l2_p1 - x_l2_p1 * gradient_l2 ;
	x_l2_inter = x_l2_p1 - y_l2_p1 / gradient_l2 ;


 	cout << "dx: " << (x_l2_p1 - x_l2_p2) << "    dy: " << (y_l2_p1 - y_l2_p2) <<  "    g: "  << gradient_l2 << endl;
 	cout << "x intercept   " << x_l2_inter << "       y intercept    " << y_l2_inter << endl;


/*--------------------------------------------------*/	
/// Line #3
	Double_t x_l3_p1, y_l3_p1;
	Double_t x_l3_p2, y_l3_p2;
	Double_t gradient_l3;
	Double_t x_l3_inter, y_l3_inter;
	
	diamond_cut->GetPoint(2, x_l3_p1, y_l3_p1);
	diamond_cut->GetPoint(3, x_l3_p2, y_l3_p2);
	
	gradient_l3 = (y_l3_p1 - y_l3_p2) / (x_l3_p1 - x_l3_p2) ;

//	cout << x_1  << "  " << y_1 << "  " << x_2 << "  " << y_2 << endl;
	y_l3_inter = y_l3_p1 - x_l3_p1 * gradient_l3 ;
	x_l3_inter = x_l3_p1 - y_l3_p1 / gradient_l3 ;


 	cout << "dx: " << (x_l3_p1 - x_l3_p2) << "    dy: " << (y_l3_p1 - y_l3_p2) <<  "    g: "  << gradient_l3 << endl;
 	cout << "x intercept   " << x_l3_inter << "       y intercept    " << y_l3_inter << endl;

/*--------------------------------------------------*/	
/// Line #4
	Double_t x_l4_p1, y_l4_p1;
	Double_t x_l4_p2, y_l4_p2;
	Double_t gradient_l4;
	Double_t x_l4_inter, y_l4_inter;
	
	diamond_cut->GetPoint(3, x_l4_p1, y_l4_p1);
	diamond_cut->GetPoint(4, x_l4_p2, y_l4_p2);
	
	gradient_l4 = (y_l4_p1 - y_l4_p2) / (x_l4_p1 - x_l4_p2) ;

//	cout << x_1  << "  " << y_1 << "  " << x_2 << "  " << y_2 << endl;
	y_l4_inter = y_l4_p1 - x_l4_p1 * gradient_l4 ;
	x_l4_inter = x_l4_p1 - y_l4_p1 / gradient_l4 ;


 	cout << "dx: " << (x_l4_p1 - x_l4_p2) << "    dy: " << (y_l4_p1 - y_l4_p2) <<  "    g: "  << gradient_l4 << endl;
 	cout << "x intercept   " << x_l4_inter << "       y intercept    " << y_l4_inter << endl;


	TString line1;
	line1.Form("W < (%f * Q2 + %f)", gradient_l1, y_l1_inter);

	TString line2;
	line2.Form("W < (%f * Q2 + %f)", gradient_l2, y_l2_inter);

	TString line3;
	line3.Form("W > (%f * Q2 + %f)", gradient_l3, y_l3_inter);

	TString line4;
	line4.Form("W > (%f * Q2 + %f)", gradient_l4, y_l4_inter);



	/// Combining 4 boundary lines into a cut (String).

	Diamond_Cut = line1 + " && " + line2 + " && " + line3 + " && " + line4;

	cout << Diamond_Cut << endl;

//	exit(0);




// 	hms_accept_cut = HMS_cut();
// 	sos_accept_cut = SOS_cut();
// 	omega_pid_cut  = Omega_PID_cut();

// 	hms_accept_cut = Omega_HMS_Cut_Garth();
// 	sos_accept_cut = Omega_SOS_Cut_Garth();
// 	omega_pid_cut  = Omega_PID_Cut_Garth();
 

	hms_accept_cut = Omega_HMS_Cut();
	sos_accept_cut = Omega_SOS_Cut();
	omega_pid_cut  = Omega_PID_Cut_Garth();
	missmass_cut   = Omega_Missmass_Cut();

	///*--------------------------------------------------*/
	/// Definition of u and t 

	if(is_sim_data) {

		u_plot_var = "minus_u";
        t_plot_var = "t";
		omega_u_cut    = Omega_U_Sim_Cut();

	} else {

		/*--------------------------------------------------*/
		/// Definition of u from Garth Huber:
		//
		// 		-u = 2*0.938272*0.938272 - 0.78265*0.78265 - t + W**2


		/*--------------------------------------------------*/
		/// Definition of u from Garth Huber:
		//
		// 		-t = Q2 - 2*0.938272*0.938272 - 0.78265*0.78265 - t + W**2
		// 		-u = Pm**2-(0.938272-Em)**2
		//

		u_plot_var = "(Pm**2-(0.938272-Em)**2)";
        t_plot_var = "-2*0.938272*(0.938272-sqrt(hsp**2 + 0.938272 * 0.938272))";

		omega_u_cut    = Omega_U_Cut();

	}

// 	delete diamond_cut_245;
// 	delete diamond_cut_160;


	acceptence_cut = hms_accept_cut + " && " + sos_accept_cut;

//	all_coin_cut = acceptence_cut + " && " + omega_pid_cut + " && " + omega_u_cut +  " && " + missmass_cut;

	all_coin_cut = acceptence_cut + " && " + omega_pid_cut  +  " && " + missmass_cut;


//	full_sim_cut = acceptence_cut + " && " + omega_u_cut + " && " + missmass_cut;

	full_sim_cut = acceptence_cut + " && " + missmass_cut;

	cout << endl << all_coin_cut << endl;

//	exit(0);



	weight_str = "Weight" + weight_app;

//	weight_str = " 1e-09 " + weight_app;

	cout << weight_str << endl;

//	exit(0);


//	weight_str = "Weight";


}



// 
// 
// 
// 
// 
/*--------------------------------------------------*/
/// Calculate the efficiency

Int_t Analysis_omega::Calculate_Yield(TH2F* hist_target, Float_t hsbeta_limit) {

// 	float x_mean = hist_target->GetMean(1);

 	float x_mean = coin_center;

//	cout << Get_Run_Num()<< "     " << coin_center << endl; 
//	exit(0);

 	float y_mean = hist_target->GetMean(2);

	coin_center_off = x_mean;
	beta_center_off = y_mean;

// 
//	float x_mean = coin_center;
//	float y_mean = hist_target->GetMean(2);

	float y_max  = hist_target->GetYaxis()->GetXmax();
	float y_min  = hist_target->GetYaxis()->GetXmin();

	/*--------------------------------------------------*/
	Int_t good_events_beta_cut = hist_target->Integral(hist_target->GetXaxis()->FindBin(x_mean-cointime_cut_l), 
			hist_target->GetXaxis()->FindBin(x_mean + cointime_cut_r), 
			hist_target->GetYaxis()->FindBin(hsbeta_limit), 
			hist_target->GetYaxis()->FindBin(y_max)); 

	Int_t random_late_beta_cut = hist_target->Integral(hist_target->GetXaxis()->FindBin(x_mean + rand_late_l), 
			hist_target->GetXaxis()->FindBin(x_mean + rand_late_r), 
			hist_target->GetYaxis()->FindBin(hsbeta_limit), 
			hist_target->GetYaxis()->FindBin(y_max));

	Int_t random_early_beta_cut = hist_target->Integral(hist_target->GetXaxis()->FindBin(x_mean - rand_early_l), 
			hist_target->GetXaxis()->FindBin(x_mean - rand_early_r), 
			hist_target->GetYaxis()->FindBin(hsbeta_limit), 
			hist_target->GetYaxis()->FindBin(y_max)); 


	Int_t tail_count_good  = 0;
	Int_t tail_count_early = 0;
	Int_t tail_count_late  = 0;

	Int_t blank  = 0;
	
	Float_t top_good;
	Float_t down_good;

	Float_t top_early;
	Float_t down_early;

	Float_t top_late;
	Float_t down_late;



	for (int j = 0; j < hist_target->GetYaxis()->GetNbins(); j++) {
		for (int i = 0; i < hist_target->GetXaxis()->GetNbins(); i++) {

				top_good  = Calculate_miss_x (x_mean+cointime_cut_r, global_centered_hsbeta_cut, hist_target->GetYaxis()->GetBinCenter(j));

				down_good = Calculate_miss_x (x_mean-cointime_cut_l, global_centered_hsbeta_cut, hist_target->GetYaxis()->GetBinCenter(j));

				top_early  = Calculate_miss_x (x_mean-rand_early_r, global_centered_hsbeta_cut, hist_target->GetYaxis()->GetBinCenter(j));
				down_early = Calculate_miss_x (x_mean-rand_early_l, global_centered_hsbeta_cut, hist_target->GetYaxis()->GetBinCenter(j));

				top_late  = Calculate_miss_x (x_mean+rand_late_r, global_centered_hsbeta_cut, hist_target->GetYaxis()->GetBinCenter(j));

				down_late = Calculate_miss_x (x_mean+rand_late_l, global_centered_hsbeta_cut, hist_target->GetYaxis()->GetBinCenter(j));


				if ( hist_target->GetYaxis()->GetBinCenter(j) > -0.6 && 
					 hist_target->GetYaxis()->GetBinCenter(j) < global_centered_hsbeta_cut) {

					if ( hist_target->GetXaxis()->GetBinCenter(i) < top_good &&  
						 hist_target->GetXaxis()->GetBinCenter(i) > down_good  ) {

//   						 cout << "top and down check "  <<  hist_target->GetYaxis()->GetBinCenter(j) <<  "  "
// 							 << hist_target->GetXaxis()->GetBinCenter(i) << "  " << top_good << "  " << down_good
// 							 << "  "  << hist_target->GetBinContent(i, j) << endl;			
// 
						 tail_count_good = tail_count_good + hist_target->GetBinContent(i, j);

					} else if (hist_target->GetXaxis()->GetBinCenter(i) < top_early &&  
						 hist_target->GetXaxis()->GetBinCenter(i) > down_early  ) {

//   						 cout << "top and down check "  <<  hist_target->GetYaxis()->GetBinCenter(j) <<  "  "
// 							 << hist_target->GetXaxis()->GetBinCenter(i) << "  " << top_early << "  " << down_early
// 							 << "  "  << hist_target->GetBinContent(i, j) << endl;			
// 
						 tail_count_early = tail_count_early + hist_target->GetBinContent(i, j);

					} else if (hist_target->GetXaxis()->GetBinCenter(i) < top_late &&  
						 hist_target->GetXaxis()->GetBinCenter(i) > down_late  ) {

//   						 cout << "top and down check "  <<  hist_target->GetYaxis()->GetBinCenter(j) <<  "  "
// 							 << hist_target->GetXaxis()->GetBinCenter(i) << "  " << top_late << "  " << down_late
// 							 << "  "  << hist_target->GetBinContent(i, j) << endl;			

						 tail_count_late = tail_count_late+ hist_target->GetBinContent(i, j);

					} else {
						blank = blank + hist_target->GetBinContent(i, j);		
					}

			}
		}
    }


	Int_t good_events_total  = good_events_beta_cut + tail_count_good;
	Int_t random_early_total = random_late_beta_cut + tail_count_late;
	Int_t random_late_total  = random_early_beta_cut + tail_count_early;

	/*--------------------------------------------------*/

	if (is_run_dummy) {

		cout <<" dummy Correction " << endl;  
		good_events_total  = dummy_correction * good_events_total;
        random_early_total = dummy_correction * random_early_total;
        random_late_total  = dummy_correction * random_late_total;

	} 

	yield = (float) good_events_total - ((float)random_early_total + (float)random_late_total)/random_real_ratio;;

	cout << "BB yield           " << yield << endl;

	float real_err = 0.0;
	float rand_err = 0.0;

// 	real_err = sqrt((float)good_events_total);
// 	rand_err = sqrt((float)random_early_total + (float)random_late_total)/random_real_ratio;

	if (is_run_dummy) {
		rand_err = sqrt((float)random_early_total + (float)random_late_total) * dummy_correction / random_real_ratio;
	} else { 
		rand_err = sqrt((float)random_early_total + (float)random_late_total) / random_real_ratio;
	}

 	if (is_run_dummy) {
		real_err = sqrt((float)good_events_total) * dummy_correction;
	} else { 
		real_err = sqrt((float)good_events_total);
	}

	yield_err = sqrt( pow(real_err,2) + pow(rand_err, 2) );

//	yield_err = 0.0;

/*--------------------------------------------------*/
/*--------------------------------------------------*/
// 	yield_setting = yield_setting + yield/Get_Total_Eff();
// 	charge_setting = charge_setting + Get_Charge();


//	acccc_temp_setting = acccc_temp_setting + acccc_temp;

	yield_setting     = yield_setting + yield;
// 	charge_setting = charge_setting + Get_Charge()*Get_Total_Eff();

	yield_setting_err = yield_setting_err + pow(yield_err,2);



/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Yield Calculation Example
// // 	yield_setting_err = yield_setting_err/pow(yield_setting, 2)  + errrr_temp/pow(acccc_temp, 2);
// // 
// // 	yield_setting     = yield_setting / acccc_temp ;
// //  yield_setting_err =  yield_setting * sqrt(yield_setting_err);
// // 
// // 	yield_setting     = yield_setting / 1000.;
// // 	yield_setting_err = yield_setting_err / 1000.;
/*--------------------------------------------------*/
/*--------------------------------------------------*/


/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Calculate Normalized yield and error run by run

	Double_t run_acc_err, run_acc;

	run_acc = Get_Charge()*Get_Total_Eff();

	run_acc_err = pow(Get_Charge_Err(), 2) + pow(Get_Total_Eff_Err(), 2); 

	yield_err = pow(yield_err,2)/pow(yield, 2) + run_acc_err;

	yield = yield / run_acc;
	
	yield_err = yield * sqrt(yield_err);

	yield = yield/1000;

	yield_err = yield_err/1000;

	graph_yield_check->SetPoint(graph_yield_check->GetN(), Get_Run_Num(), yield);

	graph_yield_check->SetPointError(graph_yield_check->GetN()-1, 0, yield_err);




	cout << "Yield:   " << yield  
//		 << "   Yield Error:  " << yield_err << endl 
// 		 << random_early_total << "    " << pow(sqrt(random_early_total)/4,2) << endl 
// 		 << random_late_total << "    " << pow(sqrt(random_late_total)/4,2) << endl
		 << "     Efficiency " << Get_Total_Eff()
		 << "     Charge " << Get_Charge() 
		 << "     acccc " << acccc_temp << endl << endl;


	run_tree->Fill();

	Print_out_data();

	return 0;

}




/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// 





Int_t Analysis_omega::Calculate_Yield_Direct() {

 	float x_mean = coin_center;
// 	float y_mean = hist_target->GetMean(2);
 	float y_mean = hsbeta_center;


	cout << hsbeta_center << endl;

	// exit(0);	

	float y_center = y_mean;

/*--------------------------------------------------*/


	TString coin_str;

	if (coin_center < 0){
		coin_str.Form("cointime+%f", fabs(coin_center));
	} else { 
		coin_str.Form("cointime-%f", coin_center);
	} 

	TString hsbeta_str;
//	hsbeta_str.Form("hsbeta-%f", y_mean);

	hsbeta_str = "hsbeta";

	coin_cut = new TCutG("coin_cut",6);

	coin_cut->SetVarX(coin_str);
	coin_cut->SetVarY(hsbeta_str);

   	coin_cut->SetPoint(0, cointime_cut_r,  y_center + top_limit);
   	coin_cut->SetPoint(1, cointime_cut_r,  y_center + global_centered_hsbeta_cut);
   	coin_cut->SetPoint(2, Calculate_tail_x (cointime_cut_r,  y_center), y_center-bot_limit);
   	coin_cut->SetPoint(3, Calculate_tail_x (-cointime_cut_l, y_center), y_center-bot_limit);
   	coin_cut->SetPoint(4, -cointime_cut_l, y_center + global_centered_hsbeta_cut);
   	coin_cut->SetPoint(5, -cointime_cut_l, y_center + top_limit);
   	coin_cut->SetPoint(6, cointime_cut_r,  y_center + top_limit);

/*--------------------------------------------------*/

	rand_early_cut = new TCutG("rand_early_cut",6);

	rand_early_cut->SetVarX(coin_str);
	rand_early_cut->SetVarY(hsbeta_str);

	rand_early_cut->SetPoint(0, -rand_early_r, y_center + top_limit);
	rand_early_cut->SetPoint(1, -rand_early_r, y_center + global_centered_hsbeta_cut);
   	rand_early_cut->SetPoint(2, Calculate_tail_x (-rand_early_r, y_center), y_center-bot_limit);
   	rand_early_cut->SetPoint(3, Calculate_tail_x (-rand_early_l, y_center), y_center-bot_limit);
	rand_early_cut->SetPoint(4, -rand_early_l, y_center + global_centered_hsbeta_cut);
	rand_early_cut->SetPoint(5, -rand_early_l, y_center + top_limit);
   	rand_early_cut->SetPoint(6, -rand_early_r, y_center + top_limit);

/*--------------------------------------------------*/

	rand_late_cut  = new TCutG("rand_late_cut", 6);

	rand_late_cut->SetVarX(coin_str);
	rand_late_cut->SetVarY(hsbeta_str);

	rand_late_cut->SetPoint(0, rand_late_r, y_center + top_limit);
	rand_late_cut->SetPoint(1, rand_late_r, y_center + global_centered_hsbeta_cut);
   	rand_late_cut->SetPoint(2, Calculate_tail_x(rand_late_r, y_center), y_center-bot_limit);
   	rand_late_cut->SetPoint(3, Calculate_tail_x(rand_late_l, y_center), y_center-bot_limit);
 	rand_late_cut->SetPoint(4, rand_late_l, y_center + global_centered_hsbeta_cut);
 	rand_late_cut->SetPoint(5, rand_late_l, y_center + top_limit);
 	rand_late_cut->SetPoint(6, rand_late_r, y_center + top_limit);



/*--------------------------------------------------*/

	TCutG* coin_cut_corrected = (TCutG*)coin_cut->Clone();
	coin_cut_corrected->SetName("ct_corrected_cut");
	coin_cut_corrected->SetVarX("ct_corrected");
	coin_cut_corrected->SetVarY("hsbeta");

/*--------------------------------------------------*/	
/// Line #1
	Double_t x_l1_p1, y_l1_p1;
	Double_t x_l1_p2, y_l1_p2;
	Double_t gradient_l1;
	Double_t x_l1_inter, y_l1_inter;

	coin_cut_corrected->GetPoint(0, x_l1_p1, y_l1_p1);
	coin_cut_corrected->GetPoint(1, x_l1_p2, y_l1_p2);
	
	gradient_l1 = (y_l1_p1 - y_l1_p2) / (x_l1_p1 - x_l1_p2) ;

//	cout << x_1  << "  " << y_1 << "  " << x_2 << "  " << y_2 << endl;
	y_l1_inter = y_l1_p1 - x_l1_p1 * gradient_l1 ;
	x_l1_inter = x_l1_p1 - y_l1_p1 / gradient_l1 ;

 	cout << "dx: " << (x_l1_p1 - x_l1_p2) << "    dy: " << (y_l1_p1 - y_l1_p2) <<  "    g: "  << gradient_l1 << endl;
 	cout << "x intercept   " << x_l1_inter << "       y intercept    " << y_l1_inter << endl;

// 
// 	cout << x_1  << "  " << y_1 << "  " << x_2 << "  " << y_2 << endl;
// 	y_inter = y_2 - x_2 * gradient ;
// 	x_inter = x_2 - y_2 / gradient ;
// 
// 
// 	cout << "dx: " << (x_1 - x_2) << "    dy: " << (y_1 - y_2) <<  "    g: "  << gradient << endl;
// 	cout << "x intercept check  " << x_inter << "       y intercept    " << y_inter << endl;


/*--------------------------------------------------*/	
/// Line #2
	Double_t x_l2_p1, y_l2_p1;
	Double_t x_l2_p2, y_l2_p2;
	Double_t gradient_l2;
	Double_t x_l2_inter, y_l2_inter;
	
	coin_cut_corrected->GetPoint(1, x_l2_p1, y_l2_p1);
	coin_cut_corrected->GetPoint(2, x_l2_p2, y_l2_p2);
	
	gradient_l2 = (y_l2_p1 - y_l2_p2) / (x_l2_p1 - x_l2_p2) ;

//	cout << x_1  << "  " << y_1 << "  " << x_2 << "  " << y_2 << endl;
	y_l2_inter =  y_l2_p1 - x_l2_p1 * gradient_l2;
	x_l2_inter = x_l2_p1 - y_l2_p1 / gradient_l2 ;



	cout << "!!! " << x_l2_p1 << "    " << y_l2_p1 << "     " << x_l2_p2 << "    " << y_l2_p2  << endl;
 	cout << "dx: " << (x_l2_p1 - x_l2_p2) << "    dy: " << (y_l2_p1 - y_l2_p2) <<  "    g: "  << gradient_l2 << endl;
 	cout << "x intercept   " << x_l2_inter << "       y intercept    " << y_l2_inter << endl;

	// exit(0);


/*--------------------------------------------------*/	
/// Line #3
	Double_t x_l3_p1, y_l3_p1;
	Double_t x_l3_p2, y_l3_p2;
	Double_t gradient_l3;
	Double_t x_l3_inter, y_l3_inter;
	
	coin_cut_corrected->GetPoint(2, x_l3_p1, y_l3_p1);
	coin_cut_corrected->GetPoint(3, x_l3_p2, y_l3_p2);
	
	gradient_l3 = (y_l3_p1 - y_l3_p2) / (x_l3_p1 - x_l3_p2) ;

//	cout << x_1  << "  " << y_1 << "  " << x_2 << "  " << y_2 << endl;
	y_l3_inter = y_l3_p1 - x_l3_p1 * gradient_l3 ;
	x_l3_inter = x_l3_p1 - y_l3_p1 / gradient_l3 ;


 	cout << "dx: " << (x_l3_p1 - x_l3_p2) << "    dy: " << (y_l3_p1 - y_l3_p2) <<  "    g: "  << gradient_l3 << endl;
 	cout << "x intercept   " << x_l3_inter << "       y intercept    " << y_l3_inter << endl;

/*--------------------------------------------------*/	
/// Line #4
	Double_t x_l4_p1, y_l4_p1;
	Double_t x_l4_p2, y_l4_p2;
	Double_t gradient_l4;
	Double_t x_l4_inter, y_l4_inter;
	
	coin_cut_corrected->GetPoint(3, x_l4_p1, y_l4_p1);
	coin_cut_corrected->GetPoint(4, x_l4_p2, y_l4_p2);
	
	gradient_l4 = (y_l4_p1 - y_l4_p2) / (x_l4_p1 - x_l4_p2) ;

//	cout << x_1  << "  " << y_1 << "  " << x_2 << "  " << y_2 << endl;
	y_l4_inter = x_l4_p1 * gradient_l4 - y_l4_p1;
	x_l4_inter = x_l4_p1 - y_l4_p1 / gradient_l4;


 	cout << "dx: " << (x_l4_p1 - x_l4_p2) << "    dy: " << (y_l4_p1 - y_l4_p2) <<  "    g: "  << gradient_l4 << endl;
 	cout << "x intercept   " << x_l4_inter << "       y intercept    " << y_l4_inter << endl;


/*--------------------------------------------------*/	
/// Line #5
	Double_t x_l5_p1, y_l5_p1;
	Double_t x_l5_p2, y_l5_p2;
	Double_t gradient_l5;
	Double_t x_l5_inter, y_l5_inter;
	
	coin_cut_corrected->GetPoint(4, x_l5_p1, y_l5_p1);
	coin_cut_corrected->GetPoint(5, x_l5_p2, y_l5_p2);
	
	gradient_l5 = (y_l5_p1 - y_l5_p2) / (x_l5_p1 - x_l5_p2) ;

//	cout << x_1  << "  " << y_1 << "  " << x_2 << "  " << y_2 << endl;
	y_l4_inter = y_l5_p1 - x_l5_p1 * gradient_l4 ;
	x_l4_inter = x_l5_p1 - y_l5_p1 / gradient_l4 ;


 	cout << "dx: " << (x_l5_p1 - x_l5_p2) << "    dy: " << (y_l5_p1 - y_l5_p2) <<  "    g: "  << gradient_l5 << endl;
 	cout << "x intercept   " << x_l5_inter << "       y intercept    " << y_l5_inter << endl;


/*--------------------------------------------------*/	
/// Line #6
	Double_t x_l6_p1, y_l6_p1;
	Double_t x_l6_p2, y_l6_p2;
	Double_t gradient_l6;
	Double_t x_l6_inter, y_l6_inter;
	
	coin_cut_corrected->GetPoint(5, x_l6_p1, y_l6_p1);
	coin_cut_corrected->GetPoint(6, x_l6_p2, y_l6_p2);
	
	gradient_l6 = (y_l6_p1 - y_l6_p2) / (x_l6_p1 - x_l6_p2) ;

//	cout << x_1  << "  " << y_1 << "  " << x_2 << "  " << y_2 << endl;
	y_l6_inter = y_l6_p1 - x_l6_p1 * gradient_l6 ;
	x_l6_inter = x_l6_p1 - y_l6_p1 / gradient_l6 ;


 	cout << "dx: " << (x_l6_p1 - x_l6_p2) << "    dy: " << (y_l6_p1 - y_l6_p2) <<  "    g: "  << gradient_l6 << endl;
 	cout << "x intercept   " << x_l6_inter << "       y intercept    " << y_l6_inter << endl;


	TString box1_1;
	box1_1.Form("ct_corrected <  %f", cointime_cut_r);

	TString box1_2;
	box1_2.Form("hsbeta > %f", y_center + global_centered_hsbeta_cut);

	TString box1_3;
	box1_3.Form("ct_corrected > %f", -cointime_cut_l);

	TString box1_4;
	box1_4.Form("hsbeta < %f", y_center + top_limit);





	TString box2_1;
	box2_1.Form("hsbeta <= %f", y_center + global_centered_hsbeta_cut);

	TString box2_2;
	box2_2.Form("hsbeta >= (%f * ct_corrected + %f)", gradient_l2, y_l2_inter);

	TString box2_3;
	box2_3.Form("hsbeta > %f", y_center-bot_limit);

	TString box2_4;
	box2_4.Form("hsbeta <= (%f * ct_corrected + %f)", gradient_l4, y_l4_inter);




// 	TString line2;
// 	line2.Form("hsbeta > (%f * ct_corrected + %f)", gradient_l2, y_l2_inter);
// 
// 	TString line3;
// 	line3.Form("hsbeta > %f", y_l3_inter);
// 
// 	TString line4;
// 	line4.Form("hsbeta < (%f * ct_corrected + %f)", gradient_l4, y_l4_inter);
// 
// 	TString line5;
// 	line5.Form("ct_corrected > %f", -cointime_cut_l);
// 
// 	TString line6;
// 	line6.Form("hsbeta < %f", y_l6_inter);

//	Cointime_Cut = line1 + " && " + line2 + " && " + line3 + " && " + line4 + " && " + line5 + " && "+ line6;
//	Cointime_Cut = line1 + " && " + line3 + " && " + line5 + " && "+ line6;

//	Cointime_Cut = line1 + " && " + line4 ;


	TString box1;

	box1 = box1_1 + " && " + box1_2  + " && "+ box1_3 + " && "+ box1_4;

	TString box2;

	// box2 = box2_1 + " && " + box2_2  + " && "+ box2_3 + " && "+ box2_4;


	box2 = box2_1 + " && " + box2_2+ " && "+ box2_3 + " && " + box2_4;


//	Cointime_Cut = "(" + box1 + ")" ;








	/*--------------------------------------------------*/
	/*--------------------------------------------------*/

	coin_cor_cut = new TCutG("coin_cor_cut", 4);

	coin_cor_cut->SetVarX(coin_str);
	coin_cor_cut->SetVarY(hsbeta_str);

   	coin_cor_cut->SetPoint(0,  cointime_cut_r,  0.1);
   	coin_cor_cut->SetPoint(1,  cointime_cut_r, -0.1);
   	coin_cor_cut->SetPoint(2, -cointime_cut_l, -0.1);
   	coin_cor_cut->SetPoint(3, -cointime_cut_l,  0.1);
   	coin_cor_cut->SetPoint(4,  cointime_cut_r,  0.1);

	/*--------------------------------------------------*/

	rand_early_cor_cut = new TCutG("rand_early_cor_cut", 4);
                  
	rand_early_cor_cut->SetVarX(coin_str);
	rand_early_cor_cut->SetVarY(hsbeta_str);
                  
   	rand_early_cor_cut->SetPoint(0, -rand_early_r,  0.1);
   	rand_early_cor_cut->SetPoint(1, -rand_early_r, -0.1);
   	rand_early_cor_cut->SetPoint(2, -rand_early_l, -0.1);
   	rand_early_cor_cut->SetPoint(3, -rand_early_l,  0.1);
   	rand_early_cor_cut->SetPoint(4, -rand_early_r,  0.1);

	/*--------------------------------------------------*/

	rand_late_cor_cut = new TCutG("rand_late_cor_cut", 4);

	rand_late_cor_cut->SetVarX(coin_str);
	rand_late_cor_cut->SetVarY(hsbeta_str);

   	rand_late_cor_cut->SetPoint(0, rand_late_r,  0.1);
   	rand_late_cor_cut->SetPoint(1, rand_late_r, -0.1);
   	rand_late_cor_cut->SetPoint(2, rand_late_l, -0.1);
   	rand_late_cor_cut->SetPoint(3, rand_late_l,  0.1);
   	rand_late_cor_cut->SetPoint(4, rand_late_r,  0.1);




	TString coin_cor_box;
	// TString rand_early_cor_box;
	// TString rand_late_cor_box;
	
	TString box_coin_1;
	box_coin_1.Form("ct_corrected < %f", cointime_cut_r);

	TString box_coin_2;
	box_coin_2.Form("hsbeta > %f", -0.1);

	TString box_coin_3;
	box_coin_3.Form("ct_corrected > %f", -cointime_cut_l);

	TString box_coin_4;
	box_coin_4.Form("hsbeta < %f", 0.1);
	
	coin_cor_box = "(" + box_coin_1 + " && " + box_coin_2 + " && " + box_coin_3 + " && " + box_coin_4 + ")";




	///*--------------------------------------------------*/
	// Cointime_Cut = "(" + box1 + " || " + box2 + ")";

	Cointime_Cut = "(" + box1 + " || " + box2 +  "||" + coin_cor_box + ")";

	Rand_Cut = "!" + Cointime_Cut;
















///*--------------------------------------------------*/
// 	TString rand_box_1;
// 	rand_box_1.Form("ct_corrected > %f", cointime_cut_r + 0.1);
// 
// 	TString rand_box_2;
// 	rand_box_2.Form("hsbeta > %f", y_center-bot_limit);
// 
// 	TString rand_box_3;
// 	rand_box_3.Form("ct_corrected < %f", -cointime_cut_l - 0.1);
// 
// 	TString rand_box_4;
// 	rand_box_4.Form("hsbeta < %f", y_center + top_limit);
// 
// 	TString r_box1;
// 	TString r_box2;
// 
// 	r_box1 = rand_box_1 + " && " + rand_box_2 + " && " + rand_box_4;
// 
// 	r_box2 = rand_box_2 + " && " + rand_box_3 + " && " + rand_box_4;
// 
// 	Rand_Cut = r_box1 + " || " + r_box2;
///*--------------------------------------------------*/
 




//	cout << box2_1 << "   " << box2_2 << "   " << box2_3 << "   " << box2_4 << endl;
//	cout << Cointime_Cut << endl;
//	exit(0);
//	Cointime_Cut = box1;
//	cout << box1_1 << "   " << box1_2 << "   " << box1_4 << "   " << box1_4 << endl;
//	cout << "asdasdas " << Cointime_Cut << endl;


/*--------------------------------------------------*/

	TString offset_str;

	if (coin_center < 0){
//		offset_str.Form("hsbeta-%f:cointime+%f", y_mean, abs(coin_center));
		offset_str.Form("hsbeta:cointime+%f", fabs(coin_center));
	} else { 
//		offset_str.Form("hsbeta-%f:cointime-%f", y_mean, coin_center);
		offset_str.Form("hsbeta:cointime+%f", fabs(coin_center));
	} 

/*--------------------------------------------------*/

	TCanvas* v1 = new TCanvas();
	
// 	v1->Update();
//
//
//	data_tree_in->Draw(offset_str, "(coin_cut || rand_early_cut || rand_late_cut) && " + all_coin_cut);
	data_tree_in->Draw(offset_str, "(coin_cut || rand_early_cut || rand_late_cut) &&" + all_coin_cut + " && " + Diamond_Cut);
 

//	data_tree_in->Draw(offset_str, "(coin_cut || rand_early_cut || rand_late_cut || coin_cor_cut || rand_early_cor_cut || rand_late_cor_cut) &&" + all_coin_cut + " && " + Diamond_Cut);


 	coin_cut->Draw("Lsame");
 	rand_early_cut->Draw("Lsame");
 	rand_late_cut->Draw("Lsame");

	v1->Update();
	v1->Write("testtest");	
 
// 	data_tree_in->Draw(offset_str, "coin_cut && " + all_coin_cut + " && " + Diamond_Cut, "P");
//  	TGraph* coin_gr = (TGraph*)gPad->GetPrimitive("Graph");
// 	TGraph* coin_gr_clone = (TGraph*)coin_gr->Clone();
// 	Int_t good_events_total = coin_gr->GetN(); 
// 
// 	data_tree_in->Draw(offset_str, "rand_early_cut && " + all_coin_cut + " && " + Diamond_Cut, "P");
//  	TGraph* rand_early_gr = (TGraph*)gPad->GetPrimitive("Graph");
// 	TGraph* rand_early_clone = (TGraph*)rand_early_gr->Clone();
// 	Int_t random_early_total = rand_early_gr->GetN(); 
// 
// 	data_tree_in->Draw(offset_str, "rand_late_cut && " + all_coin_cut + " && " + Diamond_Cut, "P");
//  	TGraph* rand_late_gr = (TGraph*)gPad->GetPrimitive("Graph");
// 	TGraph* rand_late_clone = (TGraph*)rand_late_gr->Clone();
// 	Int_t random_late_total = rand_late_gr->GetN(); 


	data_tree_in->Draw(offset_str, "coin_cut && " + all_coin_cut + " && " + Diamond_Cut, "P");
	Int_t good_events_total = data_tree_in->GetSelectedRows(); 

	data_tree_in->Draw(offset_str, "rand_early_cut && " + all_coin_cut + " && " + Diamond_Cut, "P");
	Int_t random_early_total = data_tree_in->GetSelectedRows(); 

	data_tree_in->Draw(offset_str, "rand_late_cut && " + all_coin_cut + " && " + Diamond_Cut, "P");
	Int_t random_late_total = data_tree_in->GetSelectedRows();


	/*--------------------------------------------------*/
	/// This part is to account the strange events that has hsbeta=0

	/// Coin window
	TString strange_hsbeta_coin_str;
	strange_hsbeta_coin_str.Form("abs(hsbeta) < %f &&  (" + coin_str + ") < %f && (" + coin_str + ") > %f", 0.1 , cointime_cut_r, -cointime_cut_l);


	cout << coin_str << "   " << cointime_cut_r + 0.1  << "   " << cointime_cut_l - 0.1 << endl;

	data_tree_in->Draw(offset_str,  strange_hsbeta_coin_str + " && " + all_coin_cut + " && " + Diamond_Cut);

 	TGraph* strange_hsbeta_coin_gr = (TGraph*)gPad->GetPrimitive("Graph");
	TGraph* strange_hsbeta_coin_gr_clone = (TGraph*)strange_hsbeta_coin_gr->Clone();
	Int_t strange_hsbeta_coin_total = strange_hsbeta_coin_gr->GetN(); 

// 	/// Early events
// 
// 	TString strange_hsbeta_rand_early_str;
// 	strange_hsbeta_rand_early_str.Form("abs(hsbeta) < %f &&  ct_corrected < %f && ct_corrected > %f", 0.1 , rand_early_l + 0.1 , rand_early_r - 0.1);
// 
// 	data_tree_in->Draw(offset_str, strange_hsbeta_rand_early_str + " && " + all_coin_cut + " && " + Diamond_Cut);
// 
//  	TGraph* strange_hsbeta_rand_early_gr = (TGraph*)gPad->GetPrimitive("Graph");
// 	TGraph* strange_hsbeta_rand_early_gr_clone = (TGraph*)strange_hsbeta_rand_early_gr->Clone();
// 	Int_t strange_hsbeta_rand_early_total = strange_hsbeta_rand_early_gr->GetN(); 
// 
// 
// 	/// Late events
// 
// 	TString strange_hsbeta_rand_late_str;
// 	strange_hsbeta_rand_late_str.Form("abs(hsbeta) < %f &&  ct_corrected < %f && ct_corrected > %f", 0.1 , rand_late_l + 0.1 , rand_late_r - 0.1);
// 
// 	data_tree_in->Draw(offset_str, strange_hsbeta_rand_late_str + " && " + all_coin_cut + " && " + Diamond_Cut);
// 
//  	TGraph* strange_hsbeta_rand_late_gr = (TGraph*)gPad->GetPrimitive("Graph");
// 	TGraph* strange_hsbeta_rand_late_gr_clone = (TGraph*)strange_hsbeta_rand_late_gr->Clone();
// 	Int_t strange_hsbeta_rand_late_total = strange_hsbeta_rand_late_gr->GetN(); 
// 


	cout << good_events_total << "    " << random_early_total << "    " << random_late_total << endl;

	cout << strange_hsbeta_coin_total << endl;

//	cout << strange_hsbeta_coin_total << "    " << strange_hsbeta_rand_early_total << "    " 
//         << strange_hsbeta_rand_late_total << endl;


//	exit(0);


// 	/*--------------------------------------------------*/
// 
// 
//  	v1->Update();
// 
// 	coin_gr_clone->SetMarkerColor(2);
// 	coin_gr_clone->Draw("AP");
// 		
// 	v1->Update();
// 
// 	coin_gr_clone->GetXaxis()->SetLimits(-15.0, 15.0);
// //	coin_gr_clone->GetYaxis()->SetLimits( -1.0,  1.0);
// 
// 	coin_gr_clone->SetMaximum(2.0);
// 	coin_gr_clone->SetMinimum(0.0);
// 
// 	v1->Update();
// 
// 	rand_early_clone->SetMarkerColor(6);
// 	rand_early_clone->Draw("Psame");
// 
// 	rand_late_clone->SetMarkerColor(4);
// 	rand_late_clone->Draw("Psame");
// 	
//  	coin_cut->SetLineColor(2); 
//  	coin_cut->Draw("Lsame"); 
// 
//  	rand_early_cut->SetLineColor(6); 
//  	rand_early_cut->Draw("Lsame");
// 	
//  	rand_late_cut->SetLineColor(4); 
//  	rand_late_cut->Draw("Lsame");
// 
// 	v1->Update();
// 	v1->Write("hebeta_cointime");

	cout << "/*--------------------------------------------------*/" << endl;
	cout << good_events_total << endl; 
	cout << random_early_total << endl; 
	cout << random_late_total << endl; 
	cout << "/*--------------------------------------------------*/" << endl;



	if (is_run_dummy) {

		cout <<" dummy Correction " << endl;  
		good_events_total  = dummy_correction * good_events_total;
        random_early_total = dummy_correction * random_early_total;
        random_late_total  = dummy_correction * random_late_total;

	} 

	yield = (float) good_events_total - ((float)random_early_total + (float)random_late_total)/random_real_ratio;;

	cout << "BB yield           " << yield << endl;

	float real_err = 0.0;
	float rand_err = 0.0;

	if (is_run_dummy) {
		rand_err = sqrt((float)random_early_total + (float)random_late_total) * dummy_correction / random_real_ratio;
	} else { 
		rand_err = sqrt((float)random_early_total + (float)random_late_total) / random_real_ratio;
	}

 	if (is_run_dummy) {
		real_err = sqrt((float)good_events_total) * dummy_correction;
	} else { 
		real_err = sqrt((float)good_events_total);
	}

	yield_err = sqrt( pow(real_err,2) + pow(rand_err, 2) );

//	yield_err = 0.0;

/*--------------------------------------------------*/
/*--------------------------------------------------*/
// 	yield_setting = yield_setting + yield/Get_Total_Eff();
// 	charge_setting = charge_setting + Get_Charge();


//	acccc_temp_setting = acccc_temp_setting + acccc_temp;

	yield_setting     = yield_setting + yield;
// 	charge_setting = charge_setting + Get_Charge()*Get_Total_Eff();

	yield_setting_err = yield_setting_err + pow(yield_err,2);



/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Yield Calculation Example
// // 	yield_setting_err = yield_setting_err/pow(yield_setting, 2)  + errrr_temp/pow(acccc_temp, 2);
// // 
// // 	yield_setting     = yield_setting / acccc_temp ;
// //  yield_setting_err =  yield_setting * sqrt(yield_setting_err);
// // 
// // 	yield_setting     = yield_setting / 1000.;
// // 	yield_setting_err = yield_setting_err / 1000.;
/*--------------------------------------------------*/
/*--------------------------------------------------*/


/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Calculate Normalized yield and error run by run

	Double_t run_acc_err, run_acc;

	run_acc = Get_Charge()*Get_Total_Eff();

	run_acc_err = pow(Get_Charge_Err(), 2) + pow(Get_Total_Eff_Err(), 2); 

	yield_err = pow(yield_err,2)/pow(yield, 2) + run_acc_err;

	yield = yield / run_acc;
	
	yield_err = yield * sqrt(yield_err);

	yield = yield/1000;

	yield_err = yield_err/1000;

	graph_yield_check->SetPoint(graph_yield_check->GetN(), Get_Run_Num(), yield);

	graph_yield_check->SetPointError(graph_yield_check->GetN()-1, 0, yield_err);

	cout << "Yield:   " << yield  
//		 << "   Yield Error:  " << yield_err << endl 
// 		 << random_early_total << "    " << pow(sqrt(random_early_total)/4,2) << endl 
// 		 << random_late_total << "    " << pow(sqrt(random_late_total)/4,2) << endl
		 << "     Efficiency " << Get_Total_Eff()
		 << "     Charge " << Get_Charge() 
		 << "     acccc " << acccc_temp << endl << endl;


	run_tree->Fill();

	Print_out_data();


	v1->cd();




	data_tree_in->Draw("W:Q2", "coin_cut  &&  " + all_coin_cut, "P");
 
	TGraph *w_q2_tmp = (TGraph*)gPad->GetPrimitive("Graph");	
	TGraph *w_q2_gr = (TGraph*)w_q2_tmp->Clone();	

	v1->Write("W_Q2");

// 	TString ttl;
// 	
// 	ttl.Form("%i", Get_Run_Num());
// 
// 
// 	w_q2_gr->SetName(ttl);
// 	w_q2_gr->SetTitle(ttl);

 
	diamond_setting->Add(w_q2_gr);



//	data_tree_in_cl = data_tree_in->CopyTree("(coin_cut || rand_early_cut || rand_late_cut) && " + all_coin_cut);

	data_tree_in_cl = data_tree_in->CopyTree("(coin_cut || rand_early_cut || rand_late_cut || coin_cor_cut || rand_early_cor_cut || rand_late_cor_cut) && " + all_coin_cut);


//	data_tree_in_cl = data_tree_in->CopyTree(all_coin_cut);



//	data_tree_in_cl = data_tree_in->CopyTree("( coin_cut || coin_cor_cut ) && " + all_coin_cut);

	delete v1;

}


/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Analysis_omega::Print_out_data() {

	out_dir->cd();

	run_tree->Write("yield");
	run_tree->SetName("tree_out");

	yield = 0.0;

	yield_err = 0.0;





}


/*--------------------------------------------------*/
/*--------------------------------------------------*/
void Analysis_omega::Yield_Out() {

	TTree* setting_tree_all = TTree::MergeTrees(list);

	cout << "Check Error: " << yield_setting_err/pow(yield_setting, 2) << "      " 
		 <<  errrr_temp/pow(acccc_temp, 2) << endl;

	yield_setting_err = yield_setting_err/pow(yield_setting, 2)  + errrr_temp/pow(acccc_temp, 2);

	yield_setting     = yield_setting / acccc_temp ;
  	yield_setting_err =  yield_setting * sqrt(yield_setting_err);

	yield_setting     = yield_setting / 1000.;
	yield_setting_err = yield_setting_err / 1000.;

	cout << "Total yield:       " << yield_setting     << endl;
	cout << "Total yield Error: " << yield_setting_err << endl;
	cout << "Total charge:      " << charge_setting    << endl;
	cout << "Effective charge:  " << acccc_temp        << endl;

//	cout << "Real yield Unit :  " << yield_setting / charge_setting /1000 << endl;
//	cout << "Real yield Unit :  " << yield_setting / charge_setting /1000 * 2.4 << endl;
	cout << "Real yield Unit :  " << yield_setting << endl;

	Float_t scale_factor;

	scale_factor = 0.0;

	scale_factor = 1.0/(acccc_temp*1000);

	cout  << " WWWWWWWWWWWW !!!! : " << dummy_correction << endl;

 	if (is_run_dummy)
 		scale_factor =  dummy_correction * scale_factor;

	cout  << endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!! asdasdasdasda "<<  endl << endl;

	Set_Scale_Factor(scale_factor);

	file_out_ana->cd();
	

	Proton_Obsorption_Study(setting_tree_all);


	SetYield(yield_setting);



  	Setting_Check_Plot_Out(setting_tree_all, "missmass", 100,  0.4,   1.4  );
  	Setting_Check_Plot_Out(setting_tree_all, "Em",       100,  0.4,   1.4  );
  	Setting_Check_Plot_Out(setting_tree_all, "Pm",       100,  0.1,   1.1  );

  	Setting_Check_Plot_Out(setting_tree_all, "hsdelta",  80, -10,    10   );
 	Setting_Check_Plot_Out(setting_tree_all, "W",        80,  1,     3    );
 	Setting_Check_Plot_Out(setting_tree_all, "Q2",       80,  1,     4    );
 	Setting_Check_Plot_Out(setting_tree_all, "th_pq",    80, -0.02,  0.2  );
 	Setting_Check_Plot_Out(setting_tree_all, "phi_pq",   80,  0,     6.5  );
 	Setting_Check_Plot_Out(setting_tree_all, "Pmpar",    80,  0,     1.2  );
 	Setting_Check_Plot_Out(setting_tree_all, "Pmper",    80, -0.2,   0.8  );
 	Setting_Check_Plot_Out(setting_tree_all, "Pmoop",    80, -0.4,   0.4  );
 	Setting_Check_Plot_Out(setting_tree_all, "hsytar",   80, -2.5,   2.5  );
 	Setting_Check_Plot_Out(setting_tree_all, "hsyptar",  80, -0.04,  0.04 );
 	Setting_Check_Plot_Out(setting_tree_all, "hsxptar",  80, -0.1,   0.1  );
 	Setting_Check_Plot_Out(setting_tree_all, "ssdelta",  80, -20,    20   );
 	Setting_Check_Plot_Out(setting_tree_all, "ssytar",   80, -2.5,   2    );
 	Setting_Check_Plot_Out(setting_tree_all, "ssyptar",  80, -0.07,  0.07 );
	Setting_Check_Plot_Out(setting_tree_all, "ssxptar",  80, -0.05,  0.05 );

	Setting_Check_Plot_Out(setting_tree_all, "haero_su", 160, -5.0,  35.0 );
 
  	Setting_Check_Plot_Out(setting_tree_all, "t",        100,  3,     6   );
   	Setting_Check_Plot_Out(setting_tree_all, "u",        100,  -0.5,  1   );
//  
// // 	Setting_Check_Plot_Out("Em", 600, -0.1,  0.7  );
// // 	Setting_Check_Plot_Out("Pm",       800, -0.1,  0.7  );
// 
// //	Setting_Check_Plot_Out("W", "Q2");
// 
	Setting_Beta_Cointime_Check();
 
 	Setting_u_Phi_plot(setting_tree_all);

	Setting_Beta_Cointime_Check_Setting_Tree(setting_tree_all);



	cout << " aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa " << endl;

//	exit(0);

	TCanvas* cc1 = new TCanvas();

	cc1->cd();

	graph_yield_check->SetMarkerStyle(5);
	graph_yield_check->SetMarkerSize(1.4);
	graph_yield_check->Draw("AP");

	cc1->Update();
	cc1->Write("yield_check");

	cout << " aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa " << endl;

	cc1->Clear();
	cc1->Update();

	mm_peak_check->SetMarkerStyle(5);
	mm_peak_check->SetMarkerSize(1.4);
	mm_peak_check->Draw("AP");



	cc1->Update();
	cc1->Write("mm_peak_check");

//	Acceptance_check(Get_Acceptence_Vector(), Heep_PID_Cut(), "data");

	Overall_Dia_Plot();

	Overall_Dia_Plot_After(setting_tree_all);

	Overall_2_D_Plot(setting_tree_all, "Q2" , "u");


	cout << " BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB " << endl;


 	TDirectory* t_phi_dir;
 	t_phi_dir = file_out_ana->mkdir("u_phi_bin");
	t_phi_dir->cd();

  	Yield_u_Phi_Bin(setting_tree_all);
	Yield_u_Phi_Bin(setting_tree_all, "Em");
	Yield_u_Phi_Bin(setting_tree_all, "Pm");

 	Yield_u_Phi_Bin(setting_tree_all, "hsdelta");
 	Yield_u_Phi_Bin(setting_tree_all, "hsxptar");
 	Yield_u_Phi_Bin(setting_tree_all, "hsyptar");
 

 	Yield_u_Phi_Bin(setting_tree_all, "phi_pq");
 	Yield_u_Phi_Bin(setting_tree_all, "u");


	U_Distribution_Check(setting_tree_all, "missmass", 100,  0.4, 1.4 );
	U_Distribution_Check(setting_tree_all, "Em",       100,  0.4, 1.4 );
	U_Distribution_Check(setting_tree_all, "Pm",       100,  0.1, 1.1 );

	U_Distribution_Check(setting_tree_all, "W",         80,    1, 3   );
	U_Distribution_Check(setting_tree_all, "Q2",        80,    1, 4   );
	U_Distribution_Check(setting_tree_all, "u",         80, -0.5, 1   );
	U_Distribution_Check(setting_tree_all, "t",         80,    3, 6   );

	U_Distribution_Check(setting_tree_all, "phi_pq",       80,    0, 6.4 );
 

//	exit(0);

 	Yield_output();


	cout << " aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa " << endl;


// 	cout  << endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!! asdasdasdasda "<<  endl << endl;
// 
// 	delete t_phi_dir;
//     delete cc1;
 	delete setting_tree_all;

	cout << " BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB " << endl;

}


/*--------------------------------------------------*/
/*--------------------------------------------------*/

Float_t Analysis_omega::Calculate_miss_x (Float_t x_boundary, Float_t y_boudary, Float_t y_pos) {

	Float_t x_pos;

	x_pos = x_boundary - (y_boudary - ( y_pos )) / gradi;

	return x_pos;

}

/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Analysis_omega::Check_Plot_Out(TString plot_name) {
	
 	out_dir->cd();
	
//	data_tree_in->Draw(plot_name, all_coin_cut);

//	data_tree_in->Draw(plot_name);

//	data_tree_in->Draw(plot_name, real_coin_cut);

	TCanvas* c5 = new TCanvas(); 

	c5->cd();

	data_tree_in->Draw(plot_name, "coin_cut && " + all_coin_cut + " && " + Diamond_Cut);
	TH1F *htemp_real_tmp = (TH1F*)gPad->GetPrimitive("htemp");
	TH1F *htemp_real = (TH1F*) htemp_real_tmp->Clone();
	
	htemp_real->SetOption("E");

// 	
// 	htemp_real->SetName(plot_name);
// 	htemp_real->Write(plot_name);
// 
// 	Double_t min  = htemp_real->GetXaxis()->GetXmin();
// 	Double_t max  = htemp_real->GetXaxis()->GetXmax();
// 	Int_t bin_num = htemp_real->GetXaxis()->GetNbins();
// 
// 	// cout << min << "   " << max << "   " << bin_num << endl;
// 
// 	data_tree_in->Draw(plot_name, rand_coin_cut);


	TCanvas* c6 = new TCanvas(); 
	c6->cd();


	
   	TH1F *htemp_rand = (TH1F*)htemp_real->Clone();
  	htemp_rand->Reset();
   	htemp_rand->SetName("htemp_rand");
	
 	data_tree_in->Draw(plot_name + ">> htemp_rand", "(rand_early_cut || rand_late_cut) && " + all_coin_cut, "goff");

 // 	TH1F *htemp_rand = (TH1F*)gPad->GetPrimitive("htemp");
// // 
// //	htemp_real->SetName(plot_name);
// //	htemp_real->Write(plot_name);

	htemp_rand->SetOption("E");
	
 	htemp_real->Add(htemp_rand, -1/random_real_ratio);
 
	htemp_real->Write(plot_name);
	
	delete htemp_real;
	delete htemp_rand;

	delete c5;
	delete c6;

}

// 
// 
// 
void Analysis_omega::Check_Plot_Out(TString plot_var_1, TString plot_var_2) {
	
 	out_dir->cd();

	TCanvas* v3 = new TCanvas();

	v3->cd();	
	data_tree_in->Draw(plot_var_1 + ":" + plot_var_2, all_coin_cut, "P");

 	TGraph *gtemp = (TGraph*)gPad->GetPrimitive("Graph");	
 	gtemp->Draw("P");

	v3->Write(plot_var_1 + "_vs_" + plot_var_2);

	delete gtemp;
	delete v3;

}
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 




/*--------------------------------------------------*/
/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Analysing Heep Simulation Data


void Analysis_omega::Omega_sim_anna(Int_t run_itt) {

	is_sim_data = true;

	chain = new TChain("chain");

 	TString sim_file;

 	sim_file = "list_" + sim_ana_name + "_sim_normfac.dat"; 

	ReadFile* rff_sim = new ReadFile(sim_file);

	sim_ana = rff_sim->sim_pro;

	sim_Q2 = Float_t (sim_ana.q2_setting[run_itt]);

	cout << "      " << sim_Q2 << "     " <<  endl;

	Para_Init();

	Define_Cuts();

// 	if( sim_ana.q2_setting[run_itt] == Float_t(1.60)) {
// 
// 	/// t from 0.01 to 0.45
// 
// 		u_lower_limit[0] = 0.01;
// 		u_upper_limit[0] = 0.12;
// 
// 		u_lower_limit[1] = 0.12;
// 		u_upper_limit[1] = 0.2;
// 
// 		u_lower_limit[2] = 0.2;
// 		u_upper_limit[2] = 0.4;
// 
// 	} else if (sim_ana.q2_setting[run_itt]== Float_t(2.45)) {
// 
// 		u_lower_limit[0] = 0.01; 
// 		u_upper_limit[0] = 0.212;
// 
// 		u_lower_limit[1] = 0.212;
// 		u_upper_limit[1] = 0.33;
// 
// 		u_lower_limit[2] = 0.33;
// 		u_upper_limit[2] = 0.60;
// 
// 	}
// 





// 
// 
// 
// 
// 
// 
// /*--------------------------------------------------*/
// /*--------------------------------------------------*/
// 
// 
// 	Double_t xx, yy;
// 
// 	diamond_cut_245 = new TCutG("diamond_cut_245",4);
// 
// 	diamond_cut_245->SetVarX("Q2");
// 	diamond_cut_245->SetVarY("W");
// 
//   	diamond_cut_245->SetPoint(0, 1.842, 2.388);
//    	diamond_cut_245->SetPoint(1, 2.346,  2.279);
//    	diamond_cut_245->SetPoint(2, 3.13, 1.997);
//    	diamond_cut_245->SetPoint(3, 2.54,  2.142);
// 
//    	diamond_cut_245->GetPoint(0, xx, yy);
//    	diamond_cut_245->SetPoint(4, xx, yy);
//  
// 	cout << xx << "   " << yy << endl;
// 
// 	diamond_cut_160 = new TCutG("diamond_cut_160",4);
// 
// 	diamond_cut_160->SetVarX("Q2");
// 	diamond_cut_160->SetVarY("W");
// 
//   	diamond_cut_160->SetPoint(0, 1.17, 2.355 );
//    	diamond_cut_160->SetPoint(1, 1.59, 2.262 );
//    	diamond_cut_160->SetPoint(2, 2.08,  2.048);
//    	diamond_cut_160->SetPoint(3, 1.62,  2.16);
//    	diamond_cut_160->GetPoint(0, xx, yy);
//    	diamond_cut_160->SetPoint(4, xx, yy);
//  
// 	cout << xx << "   " << yy << endl;
// 
// 
// //   	diamond_cut_160->SetPoint(4, 1.16, 2.36);
//  
// 
// 
// //	exit(0);
// 	 
//      
//     if ( sim_ana.q2_setting[run_itt] == Float_t(1.60)) {
// 		
// 		diamond_cut = (TCutG*) diamond_cut_160->Clone();
// 
// 		cout << Get_Q2() << endl;
// 
// 	} else if ( sim_ana.q2_setting[run_itt] == Float_t(2.45)) {
// 
// 		diamond_cut = (TCutG*) diamond_cut_245->Clone();
// 		cout << Get_Q2() << endl;
// 
// 	} 
// 
// 
// //	exit(0);
// 	
// 	diamond_cut->SetName("diamond_cut");
//      
// 
// /*--------------------------------------------------*/	
// /// Line #1
// 	Double_t x_l1_p1, y_l1_p1;
// 	Double_t x_l1_p2, y_l1_p2;
// 	Double_t gradient_l1;
// 	Double_t x_l1_inter, y_l1_inter;
// 	
// 
// 	diamond_cut->GetPoint(0, x_l1_p1, y_l1_p1);
// 	diamond_cut->GetPoint(1, x_l1_p2, y_l1_p2);
// 	
// 	gradient_l1 = (y_l1_p1 - y_l1_p2) / (x_l1_p1 - x_l1_p2) ;
// 
// //	cout << x_1  << "  " << y_1 << "  " << x_2 << "  " << y_2 << endl;
// 	y_l1_inter = y_l1_p1 - x_l1_p1 * gradient_l1 ;
// 	x_l1_inter = x_l1_p1 - y_l1_p1 / gradient_l1 ;
// 
//  	cout << "dx: " << (x_l1_p1 - x_l1_p2) << "    dy: " << (y_l1_p1 - y_l1_p2) <<  "    g: "  << gradient_l1 << endl;
//  	cout << "x intercept   " << x_l1_inter << "       y intercept    " << y_l1_inter << endl;
// 
// 
// /*--------------------------------------------------*/	
// /// Line #2
// 	Double_t x_l2_p1, y_l2_p1;
// 	Double_t x_l2_p2, y_l2_p2;
// 	Double_t gradient_l2;
// 	Double_t x_l2_inter, y_l2_inter;
// 	
// 	diamond_cut->GetPoint(1, x_l2_p1, y_l2_p1);
// 	diamond_cut->GetPoint(2, x_l2_p2, y_l2_p2);
// 	
// 	gradient_l2 = (y_l2_p1 - y_l2_p2) / (x_l2_p1 - x_l2_p2) ;
// 
// //	cout << x_1  << "  " << y_1 << "  " << x_2 << "  " << y_2 << endl;
// 	y_l2_inter = y_l2_p1 - x_l2_p1 * gradient_l2 ;
// 	x_l2_inter = x_l2_p1 - y_l2_p1 / gradient_l2 ;
// 
// 
//  	cout << "dx: " << (x_l2_p1 - x_l2_p2) << "    dy: " << (y_l2_p1 - y_l2_p2) <<  "    g: "  << gradient_l2 << endl;
//  	cout << "x intercept   " << x_l2_inter << "       y intercept    " << y_l2_inter << endl;
// 
// 
// /*--------------------------------------------------*/	
// /// Line #3
// 	Double_t x_l3_p1, y_l3_p1;
// 	Double_t x_l3_p2, y_l3_p2;
// 	Double_t gradient_l3;
// 	Double_t x_l3_inter, y_l3_inter;
// 	
// 	diamond_cut->GetPoint(2, x_l3_p1, y_l3_p1);
// 	diamond_cut->GetPoint(3, x_l3_p2, y_l3_p2);
// 	
// 	gradient_l3 = (y_l3_p1 - y_l3_p2) / (x_l3_p1 - x_l3_p2) ;
// 
// //	cout << x_1  << "  " << y_1 << "  " << x_2 << "  " << y_2 << endl;
// 	y_l3_inter = y_l3_p1 - x_l3_p1 * gradient_l3 ;
// 	x_l3_inter = x_l3_p1 - y_l3_p1 / gradient_l3 ;
// 
// 
//  	cout << "dx: " << (x_l3_p1 - x_l3_p2) << "    dy: " << (y_l3_p1 - y_l3_p2) <<  "    g: "  << gradient_l3 << endl;
//  	cout << "x intercept   " << x_l3_inter << "       y intercept    " << y_l3_inter << endl;
// 
// /*--------------------------------------------------*/	
// /// Line #4
// 	Double_t x_l4_p1, y_l4_p1;
// 	Double_t x_l4_p2, y_l4_p2;
// 	Double_t gradient_l4;
// 	Double_t x_l4_inter, y_l4_inter;
// 	
// 	diamond_cut->GetPoint(3, x_l4_p1, y_l4_p1);
// 	diamond_cut->GetPoint(4, x_l4_p2, y_l4_p2);
// 	
// 	gradient_l4 = (y_l4_p1 - y_l4_p2) / (x_l4_p1 - x_l4_p2) ;
// 
// //	cout << x_1  << "  " << y_1 << "  " << x_2 << "  " << y_2 << endl;
// 	y_l4_inter = y_l4_p1 - x_l4_p1 * gradient_l4 ;
// 	x_l4_inter = x_l4_p1 - y_l4_p1 / gradient_l4 ;
// 
// 
//  	cout << "dx: " << (x_l4_p1 - x_l4_p2) << "    dy: " << (y_l4_p1 - y_l4_p2) <<  "    g: "  << gradient_l4 << endl;
//  	cout << "x intercept   " << x_l4_inter << "       y intercept    " << y_l4_inter << endl;
// 
// 
// 	TString line1;
// 	line1.Form("W < (%f * Q2 + %f)", gradient_l1, y_l1_inter);
// 
// 	TString line2;
// 	line2.Form("W < (%f * Q2 + %f)", gradient_l2, y_l2_inter);
// 
// 	TString line3;
// 	line3.Form("W > (%f * Q2 + %f)", gradient_l3, y_l3_inter);
// 
// 	TString line4;
// 	line4.Form("W > (%f * Q2 + %f)", gradient_l4, y_l4_inter);
// 
// 
// 
// 	Diamond_Cut = line1 + " && " + line2 + " && " + line3 + " && " + line4;

	cout << Diamond_Cut << endl;

//	exit(0);


/*--------------------------------------------------*/
/*--------------------------------------------------*/

	Int_t q2  = Int_t( Float_t(q2_setting)*100 );
	Double_t eps = sim_ana.epsilon[run_itt]*100.;
	Int_t hms_theta = sim_ana.hms_theta[run_itt];
	
	TString sim_in_file_name;

	if(hms_theta > 0) {

		sim_in_file_name.Form(sim_ana_name+ "_%i_%.0lf_+%i" + weight_app, q2, eps, hms_theta);

	} else if (hms_theta == 0 ) {

		sim_in_file_name.Form(sim_ana_name + "_%i_%.0lf_00" + weight_app, q2, eps);

	} else {

		sim_in_file_name.Form(sim_ana_name + "_%i_%.0lf_%i" + weight_app, q2, eps, hms_theta);

	}


 	cout << sim_in_file_name << endl;
	
 	event = sim_ana.event_num[run_itt];
 	normfac = sim_ana.normfac[run_itt];

	sim_scale_fac_str.Form("%f*" + weight_str + "*", normfac/event);

	cout << sim_scale_fac_str << endl;

//	exit(0);

//	Int_t event = 5000;
//	Double_t normfac = 5000;

//	TString run_number("heep_ebeam_4p21_q2_2p406");
	TString run_number;

	run_number = sim_in_file_name;

	data_file_dir = "sim_data/";

//	data_file_dir = "sim_data_no_rad/";

	data_file = data_file_dir + run_number + ".root";


	cout << " !!!!!!!!!!!!!!!  " << data_file << endl;
	cout << " ???NAME  "	<< file_out_ana->GetName() << endl;


//	if (!rff_sim->File_exist(data_file)) return;
	if (!rff_sim->File_exist(data_file)) exit(0);

	file_in = TFile::Open(data_file);



	cout << endl << endl << " HMS Cut: " << hms_accept_cut << endl << endl;


	data_tree_in = (TTree*)file_in->Get("h666");;

	TH1F* mm = new TH1F("mm", run_number + " mm", 100, 0.4, 1.4);
	TH1F* mm_org = new TH1F("mm_org", run_number + " mm", 100, 0.4, 1.4);
	TH1F* mm_try = new TH1F("mm_try", run_number + " mm", 100, 0.4, 1.4);



//	sim_selection_no_dia.Form("%f*Weight*", normfac/event);
//	sim_selection_no_dia = sim_selection_no_dia + " ( " + acceptence_cut + " ) ";
	sim_selection_no_dia = sim_scale_fac_str + " ( " + full_sim_cut+ " ) ";



//	sim_selection.Form("%f*Weight*", normfac/event);
//	sim_selection = sim_selection + " ( " + acceptence_cut + " && " + Diamond_Cut + " ) ";

	sim_selection = sim_scale_fac_str + " ( " + full_sim_cut + " && " + Diamond_Cut + " ) ";

//	sim_selection = sim_scale_fac_str + " ( " + full_sim_cut + " && " + Diamond_Cut + " && "  + weight_str + " > 0.0 " + " ) ";

	cout << "!!!!!!!!!!!!!!!!!!!!!!!! " << endl;

	cout << sim_selection << endl;

//	exit(0);


//	data_tree_in->Draw("missmass>>mm", sim_selection + "(" + acceptence_cut + ")" , "goff");
	data_tree_in->Draw("missmass>>mm", sim_selection, "goff");

//	exit(0);

	data_tree_in->Draw("missmass>>mm_org", acceptence_cut, "goff");
	data_tree_in->Draw("missmass>>mm_try", weight_str + " * (" + acceptence_cut +")", "goff");

	cout << endl << "Sum yield: " <<  mm->Integral(0, -1) << endl << endl;;

	

		





	TCanvas* c1 = new TCanvas();

	c1->cd();
	mm->Draw();

	c1->Update();

	c1->Print( run_number + ".eps");

	TCanvas* c2 = new TCanvas();
	c2->Update();
	mm_org->Draw();
	c2->Update();

	TCanvas* c3 = new TCanvas();
	c3->Update();
	mm_try->Draw("E");
	c3->Update();


/*--------------------------------------------------*/
/*--------------------------------------------------*/



	cout << "Now Calculating the yield: " << endl;
	
	// cout << "Get Bins: " << mm->GetNbinsX() << endl;

	float yield_raw       = 0.0;		
	float yield_weighted  = 0.0;		

	float yield_err   = 0.0;		
	float yield_err_2 = 0.0;		

	for(int i=0; i <= mm->GetNbinsX();i++) {

//		cout << i << "   " << mm->GetBinContent(i) << endl;
		
		yield_raw      = yield_raw + mm_org->GetBinContent(i) ;
		yield_weighted = yield_weighted + mm->GetBinContent(i) ;
		
		
	}

	yield_err = sqrt(yield_err_2);

	cout << "Final raw yield : "   << yield_raw << endl;
	cout << "Weighed raw yield : " << yield_weighted << endl;



	float average_weight = yield_weighted / yield_raw;
	float yield_weighted_error = sqrt(yield_raw) * average_weight;

	cout << "Average Weight: " << average_weight << endl;
	cout << "Garth's Error:  " << yield_weighted_error << endl;


//	cout << mm_try->Integral(0, -1) << "    " <<  sqrt(mm_try->Integral(0, -1)) * normfac/event << endl;


	float upppp = mm_try->Integral(0, -1)* normfac ; 


//	float pp = upppp / event;

	cout << "New yield error : " <<   sqrt(upppp)/event << endl;


	chain->Add(data_file + "/h666");

	file_out_ana->cd();


	sim_tree_acc = chain->CopyTree(full_sim_cut);
    sim_tree_acc_dia =  chain->CopyTree(full_sim_cut + "&&" + Diamond_Cut );








	Sim_Setting_Check_Plot_Out("missmass", 100,  0.4,  1.4 );

	Sim_Setting_Check_Plot_Out("Em",       100,  0.4,  1.4 );
	Sim_Setting_Check_Plot_Out("Pm",       100,  0.1,  1.1 );

 	Sim_Setting_Check_Plot_Out("W",        80,  1,    3    );
 	Sim_Setting_Check_Plot_Out("Q2",       80,  1,    4    );

 	Sim_Setting_Check_Plot_Out("th_pq",    80, -0.02, 0.2  );


 	Sim_Setting_Check_Plot_Out("phi_pq",   80,  0,    6.5  );

//	exit(0);

	Sim_Setting_Check_Plot_Out("Pmpar",    80,  0,    1.2  );
	Sim_Setting_Check_Plot_Out("Pmper",    80, -0.2,  0.8  );
	Sim_Setting_Check_Plot_Out("Pmoop",    80, -0.04, 0.04 );


//  	Sim_Setting_Check_Plot_Out("hsztar",   80, -2.5,  2.5  );




 	Sim_Setting_Check_Plot_Out("hsdelta",  80, -10,   10   );
  	Sim_Setting_Check_Plot_Out("hsytar",   80, -2.5,  2.5  );
 	Sim_Setting_Check_Plot_Out("hsyptar",  80, -0.04, 0.04 );
 	Sim_Setting_Check_Plot_Out("hsxptar",  80, -0.1, 0.1 );

	Sim_Setting_Check_Plot_Out("ssdelta",  80, -20,   20   );
	Sim_Setting_Check_Plot_Out("ssytar",   80, -2.5,  2    );
 	Sim_Setting_Check_Plot_Out("ssyptar",  80, -0.07, 0.07 );
 	Sim_Setting_Check_Plot_Out("ssxptar",  80, -0.05, 0.05 );

 	Sim_Setting_Check_Plot_Out("t",        100, 3,    6    );
 	Sim_Setting_Check_Plot_Out("u",        80,  -0.5, 1    );

// 	Sim_Setting_Check_Plot_Out("tt",        100, 3,    6    );

// 	Sim_Setting_Check_U("u",  80, 0,  1);


// 	file_out_ana->mkdir("u_phi_bin");
// 	file_out_ana->cd("u_phi_bin");

	Sim_Diamond_Plot();

	Sim_2_D_Plot("Q2", "u");



	TDirectory* u_phi_dir;
 	u_phi_dir = file_out_ana->mkdir("u_phi_bin");
	u_phi_dir->cd();



	


	

//	exit(0);


	U_Distribution_Check_Sim(chain, "missmass", 100,  0.4, 1.4 );
	U_Distribution_Check_Sim(chain, "Em",       100,  0.4, 1.4 );
	U_Distribution_Check_Sim(chain, "Pm",       100,  0.1, 1.1 );

	U_Distribution_Check_Sim(chain, "W",         80,    1, 3   );
	U_Distribution_Check_Sim(chain, "Q2",        80,    1, 4   );
	U_Distribution_Check_Sim(chain, "u",         80, -0.5, 1   );
	U_Distribution_Check_Sim(chain, "t",         80,    3, 6   );


	U_Distribution_Check_Sim(chain, "phi_pq",    80,    0, 6.4 );


	Yield_u_Phi_Bin_Sim(chain);
	Yield_u_Phi_Bin_Sim(chain, "Em");
	Yield_u_Phi_Bin_Sim(chain, "Pm");

	Yield_u_Phi_Bin_Sim(chain, "hsdelta");
	Yield_u_Phi_Bin_Sim(chain, "hsxptar");
	Yield_u_Phi_Bin_Sim(chain, "hsyptar");

	Yield_u_Phi_Bin_Sim(chain, "u");
	Yield_u_Phi_Bin_Sim(chain, "phi_pq");


// 	Yield_u_Phi_Bin_Sim(chain, "Q2");
// 	Yield_u_Phi_Bin_Sim(chain, "W");
// 
// 	Yield_u_Phi_Bin_Sim(chain, "u");
// 





 
	delete chain;



	file_in->Close();
	delete file_in;
	
	file_out_ana->Close();
	delete file_out_ana;

}





/*--------------------------------------------------*/
/*--------------------------------------------------*/




void Analysis_omega::Sim_Setting_Check_Plot_Out(TString plot_name , Int_t bin_num, Double_t lower, Double_t upper) {

	TString plot_name_var;

	if (plot_name == "u") {

		plot_name_var = u_plot_var;

	} else if (plot_name == "th_pq" ) {

		plot_name_var = "thetapq";

	} else if (plot_name == "phi_pq" ) {

		plot_name_var = "phipq";

	} else if (plot_name == "Pmpar" ) {

		plot_name_var = "pmpar";

	} else if (plot_name == "Pmper" ) {

		plot_name_var = "pmper";

	} else if (plot_name == "Pmoop" ) {

		plot_name_var = "pmoop";

	} else {
	
		plot_name_var = plot_name ;

	} 


	TH1F* mm_hist = new TH1F(plot_name, plot_name, bin_num, lower, upper);

//	cout << "No Diamond Cut:   " << sim_selection_no_dia << endl;

 	chain->Draw(plot_name_var + ">>" + plot_name, sim_selection_no_dia);

//	exit(0);

	TH1F* mm_hist_raw = new TH1F(plot_name + "_raw", plot_name + "_raw", bin_num, lower, upper);

 	chain->Draw(plot_name_var + ">>" + plot_name+ "_raw", full_sim_cut);




	Sim_Set_Histo_Err( mm_hist, mm_hist_raw);

// 	for (int i = 1;  i <= mm_hist->GetNbinsX(); i++) {
// 		
// 		Double_t raw_count  = mm_hist_raw->GetBinContent(i);
// 		Double_t real_count = mm_hist->GetBinContent(i);
// 
// 		Double_t bin_err;
// 
// 		if (real_count != 0.0) {
// 			Float_t average_weight = raw_count / real_count;
// 			bin_err = sqrt(raw_count) * average_weight;
// 		} else { 
// 			bin_err = 0.0;
// 		}
// 
// 		mm_hist->SetBinError(i, bin_err);
// 
//  	}

 	mm_hist->Write(plot_name);


	TH1F* mm_hist_dia = new TH1F(plot_name + "_dia", plot_name, bin_num, lower, upper);

 	chain->Draw(plot_name_var + ">>" + plot_name + "_dia", sim_selection);

	TH1F* mm_hist_dia_raw = (TH1F*) mm_hist_dia->Clone(plot_name+ "_dia_raw");

 	chain->Draw(plot_name_var + ">>" + plot_name+ "_dia_raw", full_sim_cut + " && " + Diamond_Cut);


	Sim_Set_Histo_Err( mm_hist_dia, mm_hist_dia_raw);

// 	for (int i = 1;  i <= mm_hist_dia->GetNbinsX(); i++) {
// 		
// 		Double_t raw_count  = mm_hist_dia_raw->GetBinContent(i);
// 		Double_t real_count = mm_hist_dia->GetBinContent(i);
// 
// 		Double_t bin_err;
// 
// 		if (real_count != 0.0) {
// 			Float_t average_weight = raw_count / real_count;
// 			bin_err = sqrt(raw_count) * average_weight;
// 		} else { 
// 			bin_err = 0.0;
// 		}
// 
// 		mm_hist_dia->SetBinError(i, bin_err);
// 
//  	}

	mm_hist_dia->Write(plot_name + "_dia");


	delete mm_hist;
	delete mm_hist_raw;

	delete mm_hist_dia;
	delete mm_hist_dia_raw;


}




/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Analysis_omega::Setting_Beta_Cointime_Check() {
	
// 	TCanvas* beta_coin_can = new TCanvas();
// 	
// 	TH2F* setting_coin_beta_check = new TH2F("setting_coin_beta_check", "coin_beta_check", 150, -20, 10,  150, 0, 2);
// 
// 	chain->Draw("hsbeta:cointime>>setting_coin_beta_check", all_coin_cut, "goff");
// 
// //	setting_coin_beta_check->Write("coin_beta_check");
// 
// 	beta_coin_can->cd();
// 
// 	setting_coin_beta_check->Draw();
// 
// 	float x_cen = setting_coin_beta_check->GetMean(1);
// 	float y_cen = setting_coin_beta_check->GetMean(2) + global_centered_hsbeta_cut;
// 
//  	DrawGoodBox(x_cen, y_cen);
//  	DrawRandomBox(x_cen, y_cen);
// 
// 	beta_coin_can->Write("coin_beta_check");
// 
// 	delete beta_coin_can;


/*--------------------------------------------------*/
/*--------------------------------------------------*/

// 	TTree *setting_tree = TTree::MergeTrees(list);
// 	
// 	TCanvas* beta_coin_can = new TCanvas();
// // 	
// // 	TH2F* setting_coin_beta_check = new TH2F("setting_coin_beta_check", "coin_beta_check", 150, -20, 10,  150, 0, 2);
// // 
// // 	chain->Draw("hsbeta:cointime>>setting_coin_beta_check", all_coin_cut, "goff");
// // 
// // //	setting_coin_beta_check->Write("coin_beta_check");
// // 
// // 	beta_coin_can->cd();
// // 
// // //	setting_coin_beta_check->Draw();
// // 
// // 	float x_cen = setting_coin_beta_check->GetMean(1);
// // 	float y_cen = setting_coin_beta_check->GetMean(2) + global_centered_hsbeta_cut;
// // 
// // 	setting_coin_beta_check->Draw("scat=0.5");
// // 
// //  	DrawGoodBox(x_cen, y_cen);
// //  	DrawRandomBox(x_cen, y_cen);
// 
// 
// 	float x_cen;
// 	float y_cen;
// 
// 	TH2F* setting_coin_beta_check_cor = new TH2F("setting_coin_beta_check_cor", "", 150, -10, 10,  150, 0, 2);
// 
// 	setting_tree->Draw("hsbeta:ct_corrected>>setting_coin_beta_check_cor", all_coin_cut, "goff");
// 
//  	beta_coin_can->Update();
// 
// 	x_cen = setting_coin_beta_check_cor->GetMean(1);
// 	y_cen = setting_coin_beta_check_cor->GetMean(2) + global_centered_hsbeta_cut;
// 
// 	setting_coin_beta_check_cor->SetMarkerColor(2);
// 
// 	//setting_coin_beta_check_cor->Draw("scat=0.5 same");
// 	setting_coin_beta_check_cor->Draw();
// 
//   	DrawGoodBox(x_cen, y_cen);
//   	DrawRandomBox(x_cen, y_cen);
// 
// 
//     setting_coin_beta_check_cor-> GetYaxis()->SetTitle("hsbeta");
//     setting_coin_beta_check_cor-> GetYaxis()->CenterTitle();
// 
//     setting_coin_beta_check_cor-> GetXaxis()->SetTitle("Coincidence");
//     setting_coin_beta_check_cor-> GetXaxis()->CenterTitle();
// 
// 
// 
// // 	beta_coin_can->Modified();
// 
//  	beta_coin_can->Update();
//  	beta_coin_can->Write("coin_beta_check");
// 
// 	delete beta_coin_can;
// //	delete setting_tree;
// 
/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Dosenot work
//
// 	TH2F *setting_tree = new TH2F(); 
// 	TH2F* coin_temp = (TH2F*) list->FindObject("coin_beta_check");
// 
// //	TH2F* setting_tree = (TH2F*) coin_temp->Clone();
// 
// //	setting_tree->Reset();
// 
// 	list->ls();
// 
// 	coin_temp->GetEntries();


	TCanvas* beta_coin_can = new TCanvas();
	
  	float x_cen;
 	float y_cen;
 
// 	x_cen = coin_beta_check_offset_setting->GetMean(1);
 	x_cen = 0;

 	y_cen = coin_beta_check_offset_setting->GetMean(2) + global_centered_hsbeta_cut;

    coin_beta_check_offset_setting->SetMarkerColor(2);
    coin_beta_check_offset_setting->Draw();
 
   	DrawGoodBox(x_cen, y_cen);
   	DrawRandomBox(x_cen, y_cen);

    coin_beta_check_offset_setting->GetYaxis()->SetTitle("hsbeta");
    coin_beta_check_offset_setting->GetYaxis()->CenterTitle();

    coin_beta_check_offset_setting->GetXaxis()->SetTitle("Coincidence");
    coin_beta_check_offset_setting->GetXaxis()->CenterTitle();

 	beta_coin_can->Update();
 	beta_coin_can->Write("coin_beta_check");


	TCutG *cutg = new TCutG("mycut",1);
	cutg->SetPoint(0,-15,-0.2);
	cutg->SetPoint(1, 15,-0.2);



	TH1D* h11 = coin_beta_check_offset_setting->ProjectionX("cointime_pro", 0, -1, "mycut");

	Double_t h_max = h11->GetMaximum() * 0.75;

	beta_coin_can->cd();
	h11->Draw();

	TLine *line1 = new TLine();
 	line1->SetLineStyle(2);
 	line1->SetLineWidth(1.4);

	line1->DrawLine(rand_late_l,   0, rand_late_l,   h_max);
	line1->DrawLine(rand_late_r,   0, rand_late_r,   h_max);
	line1->DrawLine(-rand_early_l, 0, -rand_early_l, h_max);
	line1->DrawLine(-rand_early_r, 0, -rand_early_r, h_max);

	TLine *line2 = new TLine();
 	line2->SetLineStyle(2);
 	line2->SetLineColor(2);
 	line2->SetLineWidth(1.4);
	line2->DrawLine( cointime_cut_r, 0,  cointime_cut_r, h_max);
	line2->DrawLine(-cointime_cut_l, 0, -cointime_cut_l, h_max);


	beta_coin_can->Update();

	beta_coin_can->Write("cointime_projection");

	delete line1;
	delete line2;

	delete beta_coin_can;




}









/*--------------------------------------------------*/

void Analysis_omega::Setting_Beta_Cointime_Check_Setting_Tree(TTree* obj_tree) {
	
	TCanvas* beta_coin_can = new TCanvas();

//	TTree *setting_tree = TTree::MergeTrees(list);
//  	setting_tree->Draw("hsbeta:ct_corrected", Cointime_Cut + " || " + Rand_Cut);
 

	Int_t coin_prime, coin_cor;
	Int_t rand_prime, rand_cor;
	Float_t cor_prime_ratio;

 	obj_tree->Draw("hsbeta:ct_corrected", Cointime_Cut + "&&" + Diamond_Cut + " && hsbeta > 0.1", "goff");
	coin_prime = obj_tree->GetSelectedRows();

 	obj_tree->Draw("hsbeta:ct_corrected", Cointime_Cut + "&&" + Diamond_Cut + "&& hsbeta < 0.1", "goff");
	coin_cor = obj_tree->GetSelectedRows();

 	obj_tree->Draw("hsbeta:ct_corrected", Rand_Cut + "&&" + Diamond_Cut + "&& hsbeta > 0.1", "goff");
	rand_prime = obj_tree->GetSelectedRows();

	obj_tree->Draw("hsbeta:ct_corrected", Rand_Cut + "&&" + Diamond_Cut + "&& hsbeta < 0.1", "goff");
	rand_cor = obj_tree->GetSelectedRows();

	cor_prime_ratio = Float_t(coin_cor - rand_cor/random_real_ratio) / Float_t(coin_prime - rand_prime/random_real_ratio);

// 	cout << coin_prime << "  " << rand_prime << "  " << coin_cor   << "  " << rand_cor << "  " << cor_prime_ratio << endl;
// 
// 	cout << "Rand subtracted Prime Events: "    << Float_t(coin_prime - rand_prime/random_real_ratio) << endl;
// 	cout << "Rand subtracted hsbeta=0 Events: " << Float_t(coin_cor - rand_cor/random_real_ratio)     << endl;
// 	cout << "(hsbeta=0)/Prime Ratio: " << cor_prime_ratio     << endl;

//	exit(0);


	TH1F *frame = beta_coin_can->DrawFrame(-15, -0.3, 15, 1.7);
	
 	obj_tree->Draw("hsbeta:ct_corrected", "(" + Cointime_Cut + " || "+ Rand_Cut + ")" + "&&" + Diamond_Cut + " && hsbeta < 1.5 && hsbeta > -1.5", "psame");

//	TGraph *gtemp = (TGraph*)gPad->GetPrimitive("Graph");	
//	gtemp->SetMinimum(-0.5);
	
 	coin_cut->SetLineColor(2); 
 	coin_cor_cut->SetLineColor(2);
 
 	coin_cut->Draw("Lsame"); 
	coin_cor_cut->Draw("Lsame"); 

 	rand_early_cut->SetLineColor(6); 
 	rand_early_cor_cut->SetLineColor(6); 

 	rand_early_cut->Draw("Lsame");
 	rand_early_cor_cut->Draw("Lsame"); 
	
 	rand_late_cut->SetLineColor(4); 
 	rand_late_cor_cut->SetLineColor(4); 

 	rand_late_cut->Draw("Lsame");
 	rand_late_cor_cut->Draw("Lsame");


	TString ratiostring;
	ratiostring.Form("(hsbeta=0) / Prime Ratio: %f", cor_prime_ratio);

	TString primestring;
	primestring.Form("Rand Sub Prime : %f", Float_t(coin_prime - rand_prime/random_real_ratio));

	TString randstring;
	randstring.Form("Rand Sub (hsbeta=0): %f", Float_t(coin_cor - rand_cor/random_real_ratio));


	TText* ratio_str = new TText();
	ratio_str->SetNDC();
	ratio_str->SetTextSize(0.035);
	ratio_str->DrawText(.5,.75, ratiostring);
	ratio_str->DrawText(.5,.8, primestring);
	ratio_str->DrawText(.5,.85, randstring);



	beta_coin_can->Write("hsbeta_cointime_setting_tree");	

	delete beta_coin_can;

//	delete setting_tree;

}









// 
// /*--------------------------------------------------*/

void Analysis_omega::Setting_Check_Plot_Out(TTree* obj_tree, TString plot_name , Int_t bin_num, Double_t lower, Double_t upper) {

// 
// // 	TString selection_string;
// // 
// // 	selection_string.Form("%f*(",1/(acccc_temp*1000));
// // 	
// // //	selection_string = selection_string  + all_coin_cut + ")" ;
// // 	selection_string = selection_string  + real_coin_cut + ")" ;
// // 
// // 	TH1F* mm_hist_real = new TH1F(plot_name+"_real", plot_name, bin_num, lower, upper);
// // 
// //  	chain->Draw(plot_name + ">>" + plot_name+"_real", selection_string);
// // 
// // // 	chain->Draw(plot_name + ">>" + plot_name+"_real");
// // 
// // 
// // // 
// // // 	TString selection_string1;
// // // 
// // // 	selection_string1.Form("%f*(",1/(acccc_temp*1000));
// // // 	
// // // 	selection_string1 = selection_string1  + rand_coin_cut + ")" ;
// // // 
// // // 	TH1F* mm_hist_rand = new TH1F(plot_name+"_rand", plot_name, bin_num, lower, upper);
// // // 
// // //  	chain->Draw(plot_name + ">>" + plot_name+"_rand", selection_string1);
// // // 
// // // 	mm_hist_real->Add(mm_hist_rand, -0.333333333);
// // // 
// // //	plot_name.ToLower();
// // 
// // 	mm_hist_real->Write(plot_name);
// // 
// // 	delete mm_hist_real;
// 
// 
// 
// 	TTree *setting_tree = TTree::MergeTrees(list);
// 
// 
//  	TH1F* mm_hist_real = new TH1F(plot_name+"_real", plot_name, bin_num, lower, upper);
//  
// 	chain->Draw(plot_name + ">>" + plot_name+"_real", all_coin_cut);
//  
// 
//  	TH1F* mm_hist_rand = new TH1F(plot_name+"_rand", plot_name, bin_num, lower, upper);
//  
//   	chain->Draw(plot_name + ">>" + plot_name+"_rand", all_coin_cut);
// 
// 	TH1F *mm_hist_real_sub = (TH1F*) mm_hist_real->Clone();
//  	mm_hist_real_sub->Add(mm_hist_rand, -(1/random_real_ratio));
// 
// 
// 
// 	Double_t real_bin;
// 	Double_t rand_bin;	
// 	
// 	Double_t real_bin_err;
// 	Double_t rand_bin_err;
// 
// 	Double_t sub_err_2;
// 
// 	Double_t norm_yield;
// 	Double_t norm_yield_err;
// 
// 	for (int i = 1;  i <= mm_hist_real_sub->GetNbinsX(); i++) {
// 		
// //		if ( mm_hist_real_sub->GetBinContent(i) != 0) {
// 
// // // 	yield_setting_err = yield_setting_err/pow(yield_setting, 2)  + errrr_temp/pow(acccc_temp, 2);
// // // 
// // // 	yield_setting     = yield_setting / acccc_temp ;
// // //   	yield_setting_err =  yield_setting * sqrt(yield_setting_err);
// // // 
// // // 	yield_setting     = yield_setting / 1000.;
// // // 	yield_setting_err = yield_setting_err / 1000.;
// 
// 			real_bin = mm_hist_real->GetBinContent(i);
// 			rand_bin = mm_hist_rand->GetBinContent(i);
// 	
// 			Double_t real_bin_err = sqrt(real_bin);
// 			Double_t rand_bin_err = sqrt(rand_bin) / random_real_ratio;
// 	
// 			sub_err_2 = pow(real_bin_err,2) + pow(rand_bin_err,2);
// 
// 
// 			/*--------------------------------------------------*/
// 			/*--------------------------------------------------*/
// 
// 
// 			norm_yield = mm_hist_real_sub->GetBinContent(i);
// 
// 			if (norm_yield == 0) { 
// 
// 
// 			} else {
// 
// 				norm_yield_err = sub_err_2/pow(norm_yield, 2)  + errrr_temp/pow(acccc_temp, 2);
// 
// 			}
// 
// 
// 
// 			norm_yield = norm_yield/acccc_temp;
// 				
// 			norm_yield_err =  norm_yield * sqrt(norm_yield_err);
// 
// 
// 
// 
// 			norm_yield = norm_yield / 1000;
// 			
// 			norm_yield_err = norm_yield_err / 1000;
// 
// 			cout  << i << "   "   <<  mm_hist_real_sub->GetBinCenter(i) << "    "<< norm_yield << endl;
// 			
//  			mm_hist_real_sub->SetBinContent(i, norm_yield);
//  			mm_hist_real_sub->SetBinError(i, norm_yield_err);
// 
// //		}	
// 		 		
// //
// //
// // //		cout << i << endl;		
// // 
// // 		if ( mm_hiyst_real->GetBinContent(i) != 0)
// // 			mm_hist_real->SetBinError(i, 0.1);
// // 
//  	}
// // 
// // 
// //
// //    Double_t scale_fac = 1/(acccc_temp*1000);
// 
// //	mm_hist_real_sub->Scale(scale_fac);
// 
// 	mm_hist_real_sub->SetOption("E");
// 
// 	mm_hist_real_sub->Write(plot_name);
// // 
// 
// 
// //    exit(0);
// 
// 	delete setting_tree;


	TString plot_name_var;

	if (plot_name == "u") {

//		plot_name_var = "(0.938272-sqrt(hsp**2 + 0.938272 * 0.938272))**2 - Pm**2";
//		plot_name_var = "Pm**2-(0.938272-Em)**2";
		plot_name_var = u_plot_var;

	} else if (plot_name == "t") { 

//		plot_name_var = "Q2 - 2*0.938272*0.938272 - 0.78265*0.78265 - t + W**2";
//		plot_name_var = "-2*0.938272*(0.938272-sqrt(hsp**2 + 0.938272 * 0.938272))";
		plot_name_var = t_plot_var;
//
//		plot_name_var = "t";

	} else {

		plot_name_var = plot_name;

	}


// 	TTree *setting_tree = TTree::MergeTrees(list);
// 
// 	if (plot_name == "Pmpar") {
// 
// 		plot_name == "-Pmpar";
// 
// 	} else if (plot_name == "pmper") {
// 
// 		plot_name = "Pmper";
//  }



 	TH1F* mm_hist_real = new TH1F(plot_name+"_real", plot_name, bin_num, lower, upper);
 
//  	chain->Draw(plot_name + ">>" + plot_name+"_real", real_coin_cut);
  	obj_tree->Draw(plot_name_var + ">>" + plot_name+"_real", Cointime_Cut, "goff");
 

 	TH1F* mm_hist_rand = new TH1F(plot_name+"_rand", plot_name, bin_num, lower, upper);
 
//  	chain->Draw(plot_name + ">>" + plot_name+"_rand", rand_coin_cut);
  	obj_tree->Draw(plot_name_var + ">>" + plot_name+"_rand", Rand_Cut, "goff");


	TH1F *mm_hist_real_sub = (TH1F*) mm_hist_real->Clone();
 	mm_hist_real_sub->Add(mm_hist_rand, -1/random_real_ratio);


	Double_t real_bin;
	Double_t rand_bin;	
	
	Double_t real_bin_err;
	Double_t rand_bin_err;

	Double_t sub_err_2;

	Double_t norm_yield;
	Double_t norm_yield_err;

	for (int i = 1;  i <= mm_hist_real_sub->GetNbinsX(); i++) {
		
//		if ( mm_hist_real_sub->GetBinContent(i) != 0) {

// // 	yield_setting_err = yield_setting_err/pow(yield_setting, 2)  + errrr_temp/pow(acccc_temp, 2);
// // 
// // 	yield_setting     = yield_setting / acccc_temp ;
// //   	yield_setting_err =  yield_setting * sqrt(yield_setting_err);
// // 
// // 	yield_setting     = yield_setting / 1000.;
// // 	yield_setting_err = yield_setting_err / 1000.;

			real_bin = mm_hist_real->GetBinContent(i);
			rand_bin = mm_hist_rand->GetBinContent(i);
	
			Double_t real_bin_err = sqrt(real_bin);
			Double_t rand_bin_err = sqrt(rand_bin) * 1/random_real_ratio;
	
			sub_err_2 = pow(real_bin_err,2) + pow(rand_bin_err,2);




			/*--------------------------------------------------*/
			/*--------------------------------------------------*/


			norm_yield = mm_hist_real_sub->GetBinContent(i);

			if (norm_yield == 0) { 

				norm_yield_err = 0;

			} else {

				norm_yield_err = sub_err_2/pow(norm_yield, 2)  + errrr_temp/pow(acccc_temp, 2);

			}

			norm_yield = norm_yield/acccc_temp;
				
			norm_yield_err =  norm_yield * sqrt(norm_yield_err);

			norm_yield = norm_yield / 1000;
			
			norm_yield_err = norm_yield_err / 1000;

//			cout  << i << "   "   <<  mm_hist_real_sub->GetBinCenter(i) << "    "<< norm_yield << endl;
	

			if (is_run_dummy) { 
				norm_yield = norm_yield * dummy_correction;
				norm_yield_err = norm_yield * dummy_correction;
			}


		
 			mm_hist_real_sub->SetBinContent(i, norm_yield);
 			mm_hist_real_sub->SetBinError(i, norm_yield_err);


			


//		}	
		 		
//
//
// //		cout << i << endl;		
// 
// 		if ( mm_hiyst_real->GetBinContent(i) != 0)
// 			mm_hist_real->SetBinError(i, 0.1);
// 
 	}
// 
// 
//
//    Double_t scale_fac = 1/(acccc_temp*1000);

//	mm_hist_real_sub->Scale(scale_fac);

	mm_hist_real_sub->SetOption("E");

	mm_hist_real_sub->Write(plot_name);




// 
// /*--------------------------------------------------*/
// /*--------------------------------------------------*/
// /*--------------------------------------------------*/

// 	TTree *setting_tree = TTree::MergeTrees(list);
// 
// 	Double_t real_bin;
// 	Double_t rand_bin;	
// 	
// 	Double_t real_bin_err;
// 	Double_t rand_bin_err;
// 
// 	Double_t sub_err_2;
// 
// 	Double_t norm_yield;
// 	Double_t norm_yield_err;
// 

 	TH1F* mm_hist_real_dia = new TH1F(plot_name+"_real_dia", plot_name, bin_num, lower, upper);
 
//  	chain->Draw(plot_name + ">>" + plot_name+"_real", real_coin_cut);
  	obj_tree->Draw(plot_name_var + ">>" + plot_name+"_real_dia", Diamond_Cut + " && "+ Cointime_Cut, "goff");
 

 	TH1F* mm_hist_rand_dia = new TH1F(plot_name+"_rand_dia", plot_name, bin_num, lower, upper);
 
//  	chain->Draw(plot_name + ">>" + plot_name+"_rand", rand_coin_cut);
  	obj_tree->Draw(plot_name_var + ">>" + plot_name+"_rand_dia", Diamond_Cut + " && " + Rand_Cut, "goff");


	TH1F *mm_hist_real_dia_sub = (TH1F*) mm_hist_real_dia->Clone();
 	mm_hist_real_dia_sub->Add(mm_hist_rand_dia, -1/random_real_ratio);

	for (int i = 1;  i <= mm_hist_real_dia_sub->GetNbinsX(); i++) {
		
//		if ( mm_hist_real_sub->GetBinContent(i) != 0) {

// // 	yield_setting_err = yield_setting_err/pow(yield_setting, 2)  + errrr_temp/pow(acccc_temp, 2);
// // 
// // 	yield_setting     = yield_setting / acccc_temp ;
// //   	yield_setting_err =  yield_setting * sqrt(yield_setting_err);
// // 
// // 	yield_setting     = yield_setting / 1000.;
// // 	yield_setting_err = yield_setting_err / 1000.;

			real_bin = mm_hist_real_dia->GetBinContent(i);
			rand_bin = mm_hist_rand_dia->GetBinContent(i);
	
			Double_t real_bin_err = sqrt(real_bin);
			Double_t rand_bin_err = sqrt(rand_bin) * 1/random_real_ratio;
	
			sub_err_2 = pow(real_bin_err,2) + pow(rand_bin_err,2);




			/*--------------------------------------------------*/
			/*--------------------------------------------------*/


			norm_yield = mm_hist_real_dia_sub->GetBinContent(i);

			if (norm_yield == 0) { 

				norm_yield_err = 0;

			} else {

				norm_yield_err = sub_err_2/pow(norm_yield, 2)  + errrr_temp/pow(acccc_temp, 2);

			}

			norm_yield = norm_yield/acccc_temp;
				
			norm_yield_err =  norm_yield * sqrt(norm_yield_err);

			norm_yield = norm_yield / 1000;
			
			norm_yield_err = norm_yield_err / 1000;

//			cout  << i << "   "   <<  mm_hist_real_sub->GetBinCenter(i) << "    "<< norm_yield << endl;


			if (is_run_dummy) { 
				norm_yield = norm_yield * dummy_correction;
				norm_yield_err = norm_yield * dummy_correction;
			}
			
 			mm_hist_real_dia_sub->SetBinContent(i, norm_yield);
 			mm_hist_real_dia_sub->SetBinError(i, norm_yield_err);

//		}	
		 		
//
//
// //		cout << i << endl;		
// 
// 		if ( mm_hiyst_real->GetBinContent(i) != 0)
// 			mm_hist_real->SetBinError(i, 0.1);
// 
 	}
// 
// 
//
//    Double_t scale_fac = 1/(acccc_temp*1000);

//	mm_hist_real_sub->Scale(scale_fac);

	mm_hist_real_dia_sub->SetOption("E");

	mm_hist_real_dia_sub->Write(plot_name + "_dia");


	TCanvas* c1 = new TCanvas(); 

	c1->Update();

	mm_hist_real_sub->Draw();

	c1->Update();

 	mm_hist_real_dia_sub->SetLineColor(2);
 	mm_hist_real_dia_sub->Draw("same");
// 
//	mm_hist_real_dia_sub->Draw();

	c1->Update();


	c1->Write(plot_name+ "_same");
// 
// 

// delete c1;




	delete c1;

	delete mm_hist_real;
	delete mm_hist_rand;
	delete mm_hist_real_sub;
	
	delete mm_hist_real_dia;
	delete mm_hist_rand_dia;
	delete mm_hist_real_dia_sub;

//	delete setting_tree;


}



// 
// 

// /*--------------------------------------------------*/
// 
// void Analysis_omega::Setting_Check_Plot_Out(TString plot_name) {
// 	
// 	chain->Draw(plot_name, all_coin_cut);
//  	TH1F *htemp = (TH1F*)gPad->GetPrimitive("htemp");
// 	htemp->Write(plot_name);
// 
// }
// 
// 
// 
// void Analysis_omega::Setting_Check_Plot_Out(TString plot_var_1, TString plot_var_2) {
// 	
// 	TCanvas* v3 = new TCanvas();
// 	v3->cd();	
// 
// 	chain->Draw(plot_var_1 + ":" + plot_var_2, all_coin_cut, "P");
// 
//  	TGraph *gtemp = (TGraph*)gPad->GetPrimitive("Graph");	
//  	gtemp->Draw("P");
// 
// 	v3->Write(plot_var_1 + "_vs_" + plot_var_2);
// 
// //	cout << "aasdasd a" << endl;
// 
// 	delete gtemp;
// 	delete v3;
// 
// }
// 
// 



/*--------------------------------------------------*/
/// Plot boxes


void Analysis_omega::DrawGoodBox(float x_center, float y_center) {

//	cout << x_center - 1 << endl; 

	TLine *line1 = new TLine();
 	line1->SetLineWidth(1.4);
 	line1->SetLineColor(2);
	line1->DrawLine(x_center-cointime_cut_l, y_center, x_center - cointime_cut_l, y_center + top_limit);
	line1->DrawLine(x_center+cointime_cut_r, y_center, x_center + cointime_cut_r, y_center + top_limit);
	line1->DrawLine(x_center-cointime_cut_l, y_center, x_center + cointime_cut_r, y_center);

	line1->DrawLine(x_center-cointime_cut_l, y_center, Calculate_tail_x (x_center - cointime_cut_l, y_center ), y_center - bot_limit);
	line1->DrawLine(x_center+cointime_cut_r, y_center, Calculate_tail_x (x_center + cointime_cut_r, y_center ), y_center - bot_limit);


}


// 
void Analysis_omega::DrawRandomBox(float x_center, float y_center) {

	TLine *line1 = new TLine();
 	line1->SetLineStyle(2);
 	line1->SetLineWidth(1.4);
	line1->DrawLine(x_center+rand_late_l, y_center, x_center + rand_late_l, y_center + top_limit);
	line1->DrawLine(x_center+rand_late_r, y_center, x_center + rand_late_r, y_center + top_limit);
	line1->DrawLine(x_center+rand_late_l, y_center, x_center + rand_late_r, y_center);

	line1->DrawLine(x_center+rand_late_l, y_center, Calculate_tail_x (x_center + rand_late_l, y_center ), y_center - bot_limit);
	line1->DrawLine(x_center+rand_late_r, y_center, Calculate_tail_x (x_center + rand_late_r, y_center ), y_center - bot_limit);


	TLine *line2 = new TLine();
 	line2->SetLineStyle(2);
 	line2->SetLineWidth(1.4);
	line2->DrawLine(x_center-rand_early_l, y_center, x_center - rand_early_l, y_center + top_limit);
	line2->DrawLine(x_center-rand_early_r, y_center, x_center - rand_early_r, y_center + top_limit);
	line2->DrawLine(x_center-rand_early_l, y_center, x_center - rand_early_r, y_center);

	line2->DrawLine(x_center-rand_early_r, y_center, Calculate_tail_x (x_center - rand_early_r, y_center), y_center - bot_limit);
	line2->DrawLine(x_center-rand_early_l, y_center, Calculate_tail_x (x_center - rand_early_l, y_center), y_center - bot_limit);

}
// 
// 
// 

Float_t Analysis_omega::Calculate_tail_x (Float_t xx, Float_t yy) {
	
	Float_t xx_tail;

	xx_tail = xx - (yy - (yy - bot_limit)) / gradi;

//	cout << "Calculation Check: "  << xx_tail << "    " << xx << "    " << float (yy - (-0.6)) * gradi << endl;

	return xx_tail;


}

// 
// 
// 

/*--------------------------------------------------*/

void Analysis_omega::Yield_output() {

	static int dummy_count = 0; 
	static int sim_count = 0; 
	static int target_count = 0; 

	
	TString file_name_str;

	file_name_str = file_out_ana->GetName();

//	cout << file_name_str << "    " << file_name_str.Contains("asdummy") << endl;


	if (file_name_str.Contains("dummy") )  {

		if (dummy_count ==0) {
			yield_file_out.open("dummy_yield.txt",std::fstream::out );
			dummy_count++;
		} else {
			yield_file_out.open("dummy_yield.txt",std::fstream::app );
		}


	} else if (file_name_str.Contains("sim")) {

		if (sim_count ==0) {
			yield_file_out.open("sim_yield.txt", std::fstream::out );		
			sim_count++;
		} else {
			yield_file_out.open("sim_yield.txt",std::fstream::app );
		}

	} else {

		if (target_count ==0) {
			yield_file_out.open("target_yield.txt", std::fstream::out );		
			target_count++;
		} else {
			yield_file_out.open("target_yield.txt",std::fstream::app );
		}



	}

	yield_file_out << yield_setting  << "  "<< yield_setting_err << endl;

}

// 
// 
// 
// 
// 
// /*--------------------------------------------------*/
// 
// 
// 
// 
// void Analysis_heep::Acceptance_check(vector<TString> list_vec, TString cut_str, TString data_type) {
// 
// 
// 	vector<TString>::iterator it;
// 
// 	TDirectory* acceptence_dir;
// 
// 	acceptence_dir = file_out_ana->mkdir("acceptence_check");
// 
// 	acceptence_dir->cd();
// 
//   	for (it = list_vec.begin(); it != list_vec.end(); ++it) {
// 
// 		TString tmp;
// 		tmp = *it;
// 
// 		Apt_Setting_Check(tmp, cut_str, data_type);
// 
// //    	cout << "    " << tmp << endl;
// 
// 	}
// 
// //    	cout << ' ' << *it << endl;
// // 
// // 	exit(0);
// 
// }
// 
// 
// /*--------------------------------------------------*/
// /*--------------------------------------------------*/
// 
// void Analysis_heep::Apt_Setting_Check(TString plot_name, TString cut_str, TString data_type) {
// 	
// // 	ssdelta_h     = new TH1F("ssdelta", "ssdelta", 200, -20, 20);
// // 	ssxptar_h     = new TH1F("ssxptar", "ssxptar", 200, -0.05, 0.05);
// // 	ssyptar_h     = new TH1F("ssyptar", "ssyptar", 200, -0.1, 0.1);
// // 
// // 	w_h           = new TH1F("w_h",  "w_h",  200, 0.5, 1.5);
// // 	q2_h          = new TH1F("q2_h", "q2_h", 200, 4, 6);
// // 
// //	sszbeam_h     = new TH1F("sszbeam",  "sszbeam" , 200, -5, 10);
// //
// //
// 
// // w_h           = new TH1F("w_h",  "w_h",  200, 0.5, 1.5)
// 
// 
// 
// 	
// 	TString plot_op_str;
// 	
// 
// 	if (plot_name == "ssxptar") {
// 	
// 		plot_op_str = plot_name + ">> htemp(200, -0.05, 0.05)";	
// 
// 	} else if (plot_name == "ssyptar") {
// 
// 		plot_op_str = plot_name + ">> htemp(200, -0.1, 0.1)";	
//  
// 	} else if (plot_name == "Q2") {
// 
// 		plot_op_str = plot_name + ">> htemp(200, 3, 8)";	
// 
// 	} else if (plot_name == "sszbeam") {
// 
// 		plot_op_str = plot_name + ">> htemp(200, -5, 10)";	
// 
// 	} else if (plot_name == "ssytar") {
// 
// 		plot_op_str = plot_name + ">> htemp(200, -5, 5)";	
// 
// 	} else if (plot_name == "ssdelta") {
// 
// 		plot_op_str = plot_name + ">> htemp(200, -40, 40)";	
// 
// 	} else if (plot_name == "W") {
// 
// 		plot_op_str = plot_name;
// 
// 		if( data_type == "data") {
// 
// 			plot_op_str = "w";
// 
// 		}
// 
//  		plot_op_str = plot_op_str + " >> htemp(200, 0.5, 1.5)";
// 
// 	}  else if (plot_name == "ssxfp") {
// 
// 		plot_op_str = plot_name + " >> htemp(200, -20, 15)";
// 
// 	}  else if (plot_name == "ssxpfp") {
// 
// 		plot_op_str = plot_name + " >> htemp(200, -0.15, 0.14)";
// 
// 	} else if (plot_name == "missmass") {
// 
// 		plot_op_str = plot_name + " >> htemp(200, -0.15, 0.15)";
// 
// 	} else if (plot_name == "hsdelta") {
// 
// 		plot_op_str = plot_name + " >> htemp(200, -15, 15)";
// 
// 	} else if (plot_name == "hsxptar") {
// 
// 		plot_op_str = plot_name + " >> htemp(200, -0.1, 0.1)";
// 
// 	} else if (plot_name == "hsyptar") {
// 
// 		plot_op_str = plot_name + " >> htemp(200, -0.1, 0.1)";
// 
// 	} else {
// 
// 		plot_op_str = plot_name;
// 
// 	}
// 
// 
// 
// 	// cout << " asdasda    " << cut_str << endl;
// 
// 
// 
// 
//  	chain->Draw(plot_op_str, cut_str);
//  
//    	TH1F* htempp = (TH1F*)gPad->GetPrimitive("htemp");
// // 
// // 
// 	if (data_type == "data") {
// 	
// 		htempp ->Scale(Get_Scale_Factor());
// 
// 	} 
// // 
// // 
// //
// 
// 	cout << plot_name << endl;
// 
//  	htempp->Write(plot_name);
// 
// 	delete htempp;
// 
// 
// }
// 
/*--------------------------------------------------*/

vector<TString> Analysis_omega::Get_Acceptence_Vector() {

	static const TString list_arr[] = { "missmass", "Q2","ssdelta", "ssxfp", "ssxpfp", "ssytar", "ssxptar", "ssyptar", "hsdelta", "hsxptar", "hsyptar"};

	vector<TString> vec (list_arr, list_arr + sizeof(list_arr) / sizeof(list_arr[0]) );

	return vec;

}


/*--------------------------------------------------*/


void Analysis_omega::Check_MissingMass_Peak() {

	TCanvas *c2 = new TCanvas();

	c2->cd();
	
	data_tree_in->Draw("missmass>>htemp(100, 0.4, 1.4)", real_coin_cut);
 	TH1F *htemp_real = (TH1F*)gPad->GetPrimitive("htemp");
	
//	Double_t missing_mass_max = htemp_real->GetMaximum(); 


	htemp_real->GetXaxis()->SetRangeUser(0.78, 0.82);

	c2->Update();



//	Double_t missing_mass_max = 0.783; 

	Double_t missing_mass_max = htemp_real->GetBinCenter(htemp_real->GetMaximumBin()); 


	htemp_real->GetXaxis()->SetRangeUser(0.4, 1.4);
	

	TF1* fit_fun = new TF1("my_function", "gaus");
	TF1* fit_fun_improve = new TF1("my_function_imp", "gaus");


	cout << missing_mass_max << endl;
	cout << missing_mass_max << endl;
	cout << missing_mass_max << endl;
	cout << missing_mass_max << endl;

	fit_fun->SetRange(missing_mass_max-0.025, missing_mass_max+0.02);

	htemp_real->Fit("my_function", "R");


 	fit_fun_improve = new TF1("my_function_imp", "gaus");
 	
 	fit_fun_improve->SetRange(fit_fun->GetParameter("Mean")-0.02, fit_fun->GetParameter("Mean")+0.016);

	htemp_real->Fit("my_function_imp", "R");

	fit_fun_improve->Draw("same");

	c2->Update();

	c2->Write("mm_peak_check");

	mm_peak_check->SetPoint(mm_peak_check->GetN(), Get_Run_Num(), fit_fun_improve->GetParameter("Mean"));	
	mm_peak_check->SetPointError(mm_peak_check->GetN()-1, 0, fit_fun_improve->GetParError(2));	




 	delete htemp_real;
 	delete fit_fun;
 	delete fit_fun_improve;
 
	delete c2;


}


/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Analysis_omega::Correct_Offset_Tree() { 

	Float_t missmass_cor;
	Float_t missmass_org;

	missmass_cor = 1;


    mm_cor_branch = data_tree_in_cl->Branch("ct_corrected", &missmass_cor);
    
	data_tree_in_cl->SetBranchAddress("cointime", &missmass_org);

	data_tree_in_cl->SetBranchStatus("*", 0);
	data_tree_in_cl->SetBranchStatus("ct_corrected", 1);


	//TTree* data_tree_clone = (TTree* ) data_tree_in->Clone(0);

	Int_t nentries = data_tree_in_cl->GetEntries();

    for (Int_t i=0; i < nentries; i++) { 

		data_tree_in_cl->SetBranchStatus("cointime", 1);

		data_tree_in_cl->GetEvent(i);

		if (Get_Coin_Center() < 0) { 
			missmass_cor  = missmass_org + fabs(Get_Coin_Center()); 
		} else {
			missmass_cor  = missmass_org - Get_Coin_Center(); 
		}


//  		cout << "Org  " << nentries << "   " << Get_Coin_Center() << "    " << missmass_org << endl;
//  		cout << "cor  " << missmass_cor << endl;

		data_tree_in_cl->SetBranchStatus("cointime", 0);

		mm_cor_branch->Fill();

	}


//      tree->GetEvent(i);
// 
//      //  compute variables for new branches and fill these branches
//   ...
//    bnew1->Fill();
//    bnew2->Fill();

	data_tree_in_cl->SetBranchStatus("*", 1);

//	return data_tree_in; 


}








/*--------------------------------------------------*/
/*--------------------------------------------------*/

Double_t Analysis_omega::GetBetaCenter() { 


	TString coin_cut_str;

	coin_cut_str.Form("cointime > %f && cointime <%f", coin_center - cointime_cut_l, coin_center + cointime_cut_r);


	TCanvas* v2 = new TCanvas();  

	TH1F* hsbeta_check = new TH1F("hsbeta_check", "hsbeta",  400, 0, 2);

	data_tree_in->Draw("hsbeta >> hsbeta_check", all_coin_cut + " && " + coin_cut_str, "goff");

	v2->cd();

	hsbeta_check->GetXaxis()->SetRangeUser(0.85, 1.05);

	hsbeta_check->Draw("hist");

	Double_t hsbeta_max = hsbeta_check->GetMaximum();	

	v2->Update();

	v2->Write("hsbeta_check");

	Double_t hsbeta_center_tmp = hsbeta_check->GetMean();

	// cout << hsbeta_center_tmp << " asdasdas " << endl;

	return hsbeta_center_tmp;

}







	
/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Output the overall diamond plot for the setting

void Analysis_omega::Overall_Dia_Plot() {

	TCanvas* dia_canva = new TCanvas();

	dia_canva->Update();

	diamond_setting->SetName("diamond_gr");
		
	diamond_setting->Draw("AP");
//	diamond_setting_cut->Draw("P");

	// cout << "  " << diamond_setting_cut->GetListOfGraphs()->GetSize() << endl; 
	// diamond_setting_cut->GetListOfGraphs()->ls();
	// diamond_setting_cut->GetListOfGraphs()->Dump();
	// diamond_setting_cut->GetListOfGraphs()->Print();

	TList * graph_list = diamond_setting->GetListOfGraphs();
 	TIter next(graph_list);
 	TObject* object = 0;
	Int_t event_num = 0;

 	while ((object = next()))
 	{
		event_num = event_num + ((TGraph*)object)->GetN();
 	}


// 	TList * graph_list_cut = diamond_setting_cut->GetListOfGraphs();
//  	TIter next_cut(graph_list_cut);
//  	TObject* object_cut = 0;
// 	Int_t event_num_cut = 0;
// 
//  	while ((object_cut = next_cut()))
//  	{
// 		event_num_cut = event_num_cut + ((TGraph*)object_cut)->GetN();
//  	}
// 





//	cout << "  fsafdsf  " << event_num << endl;

 	TText* number_points = new TText();
 	number_points->SetTextSize(0.035);

 	TString points_str; 
//  	points_str.Form("Total Events: %i", event_num + event_num_cut);
//  	number_points->DrawTextNDC(0.7, 0.85, points_str);
 
 	points_str.Form("Black Events: %i", event_num);
 	number_points->DrawTextNDC(0.7, 0.75, points_str);

//  	points_str.Form("Red Events: %i", event_num_cut);
//  	number_points->DrawTextNDC(0.7, 0.65, points_str);
// 



	diamond_cut->Draw("Lsame");


 	dia_canva->Update();
 	dia_canva->Write("diamand");

	delete number_points;
	delete dia_canva;


}

/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Output the overall diamond plot for the setting

void Analysis_omega::Overall_Dia_Plot_After(TTree* obj_tree) {


	TMultiGraph* mg = new TMultiGraph();

	TCanvas* dia_canva_1 = new TCanvas();

	dia_canva_1->Update();

	// TTree *setting_tree = TTree::MergeTrees(list);
 	// setting_tree->Draw("W:Q2", "", "P");

 	obj_tree->Draw("W:Q2", "", "P");


 	TGraph *diamond_setting_before_tmp = (TGraph*)gPad->GetPrimitive("Graph");	
	TGraph *diamond_setting_before = (TGraph*)diamond_setting_before_tmp->Clone();
	mg->Add(diamond_setting_before);

	TCanvas* dia_canva_2 = new TCanvas();

	dia_canva_2->Update();


 	// setting_tree->Draw("W:Q2", Diamond_Cut, "Psame");
 
	obj_tree->Draw("W:Q2", Diamond_Cut, "P");

 	TGraph *diamond_setting_after_tmp = (TGraph*)gPad->GetPrimitive("Graph");	
	TGraph *diamond_setting_after = (TGraph*)diamond_setting_after_tmp->Clone();
 	diamond_setting_after->SetMarkerColor(2);
// 	diamond_setting_after->SetDrawOption("P");

 	mg->Add(diamond_setting_after);


	TCanvas* dia_canva_3 = new TCanvas();


// //	TGraph *w_q2_gr = (TGraph*)w_q2_tmp->Clone();	
// 
// //	diamond_setting_b->Draw("AP");
// 
// 	Int_t event_num = diamond_setting_before->GetN();
// 
// //  	TText* number_points = new TText();
// //  	number_points->SetTextSize(0.035);
// // 
// //  	TString points_str; 
// // //  	points_str.Form("Total Events: %i", event_num + event_num_cut);
// // //  	number_points->DrawTextNDC(0.7, 0.85, points_str);
// //  
// //  	points_str.Form("Black Events: %i", event_num);
// //  	number_points->DrawTextNDC(0.7, 0.75, points_str);
// // 
// 
// 
 	mg->Draw("AP");
// 
//
//	diamond_setting_after->Draw();

 	diamond_cut->SetLineColor(2); 
 	diamond_cut->Draw("Lsame");
// 
  	dia_canva_3->Update();
  	dia_canva_3->Write("diamond_after_cut");
// 
// //	delete number_points;
// 	delete dia_canva;
// 
//
//
//

	delete diamond_setting_before_tmp;
	delete diamond_setting_after_tmp;
	
	delete mg;

	delete dia_canva_1;
	delete dia_canva_2;
	delete dia_canva_3;

// 	delete setting_tree;

}







/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Output the overall 2D plot for the setting

void Analysis_omega::Overall_2_D_Plot(TTree* obj_tree, TString var_1, TString var_2 ) {

	TString plot_name_str = var_1 + "_" + var_2;

	if ( var_1 == "u") {

		var_1 = u_plot_var;

	} else if (var_1 == "t") { 

		var_1 = t_plot_var;

	} 

	if ( var_2 == "u") {

		var_2 = u_plot_var;

	} else if (var_2 == "t") { 

		var_2 = t_plot_var;

	}


	TString plot_var = var_1 + ":"+ var_2;

	TMultiGraph* mg = new TMultiGraph();

	TCanvas* dia_canva_1 = new TCanvas();

	dia_canva_1->Update();

	// TTree *setting_tree = TTree::MergeTrees(list);
 	// setting_tree->Draw("W:Q2", "", "P");

 	obj_tree->Draw(plot_var, "", "P");




 	TGraph *diamond_setting_before_tmp = (TGraph*)gPad->GetPrimitive("Graph");	
	TGraph *diamond_setting_before = (TGraph*)diamond_setting_before_tmp->Clone();
	mg->Add(diamond_setting_before);


	TCanvas* dia_canva_2 = new TCanvas();

	dia_canva_2->Update();


 	// setting_tree->Draw("W:Q2", Diamond_Cut, "Psame");
 
	obj_tree->Draw(plot_var, Diamond_Cut, "P");


   	TGraph *diamond_setting_after_tmp = (TGraph*)gPad->GetPrimitive("Graph");	
	TGraph *diamond_setting_after = (TGraph*)diamond_setting_after_tmp->Clone();
 	diamond_setting_after->SetMarkerColor(2);
 	mg->Add(diamond_setting_after);

 	mg->Draw("P");
 
  	dia_canva_2->Update();
  	dia_canva_2->Write(plot_name_str);


	delete diamond_setting_before_tmp;
	delete diamond_setting_after_tmp;

	delete mg;

	delete dia_canva_1;
	delete dia_canva_2;


}







// /*--------------------------------------------------*/
// /*--------------------------------------------------*/
// /// Output the overall 2D plot for the setting
// 
// void Analysis_omega::Overall_2_D_Plot(TTree* obj_tree, TString var_1, TString var_2 ) {
// 
// 	TString plot_name_str = var_1 + "_" + var_2;
// 
// 	if ( var_1 == "u") {
// 
// 		var_1 = u_plot_var;
// 
// 	} else if (var_1 == "t") { 
// 
// 		var_1 = t_plot_var;
// 
// 	} 
// 
// 	if ( var_2 == "u") {
// 
// 		var_2 = u_plot_var;
// 
// 	} else if (var_2 == "t") { 
// 
// 		var_2 = t_plot_var;
// 
// 	}
// 
// 	TString plot_var = var_1 + ":"+ var_2;
// 
// 	TCanvas* canva_1 = new TCanvas();
// 
//  	obj_tree->Draw("Q2");
// 
// 	canva_1->Write("Q2_test");
// 
// 	delete canva_1;
// 
// 
// }
// 






/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Plotting the Polar plot in u and phi. Bull's eye plot! 

void Analysis_omega::Setting_u_Phi_plot(TTree* obj_tree) {

	TCanvas* u_phi_can = new TCanvas();

//	TTree * setting_tree = TTree::MergeTrees(list);

	obj_tree->Draw("phi_pq:" + u_plot_var, Cointime_Cut + " && " + Diamond_Cut);

 	TGraph* g_u_phi = (TGraph*)gPad->GetPrimitive("Graph");	

	TCanvas* u_phi_can_pol = new TCanvas("u_phi_can_pol", "", 700,700);

	TGraphPolar* g_u_phi_polar = new TGraphPolar(g_u_phi->GetN(), g_u_phi->GetY(), g_u_phi->GetX());

	g_u_phi_polar->SetTitle("");	


	u_phi_can_pol->cd();

	g_u_phi_polar->Draw("AOP");

	u_phi_can_pol->Update();
	
 	g_u_phi_polar->GetPolargram()->SetNdivPolar(4);
 	g_u_phi_polar->GetPolargram()->SetNdivRadial(4);
 	g_u_phi_polar->GetPolargram()->SetRangeRadial(0, 0.8); 

	u_phi_can_pol->Update();

 	u_phi_can_pol->Write("u_phi_check");


	delete u_phi_can_pol;
	delete u_phi_can;
//	delete setting_tree;

}




/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Missing mass u-phi bin

void Analysis_omega::Yield_u_Phi_Bin(TTree* obj_tree) {

	Double_t real_err = 0.0;
	Double_t rand_err = 0.0;

	Double_t real_bin, real_bin_err, norm_yield    ;
	Double_t rand_bin, rand_bin_err, norm_yield_err;	

	Double_t sub_err_2;


	TH1F* h_missmass      = new TH1F("h_miss",      "missmass_real", 70,  0.4,  1.1);
	TH1F* h_missmass_rand = new TH1F("h_miss_rand", "missmass_rand", 70,  0.4,  1.1);

	TH1F* u_check_h       = new TH1F("u_check", "u_check", 30, 0, 0.7);

	TH1F* event_phi_real  = new TH1F("phi_real", "phi_real", phi_bin_num, 0, 360);
	TH1F* event_phi_rand  = new TH1F("phi_rand", "phi_rand", phi_bin_num, 0, 360);

	TCanvas * u_phi  = new TCanvas("u_phi" ,"u_phi" , 1600, 300);
	u_phi->Divide(u_bin_num, 1);

	TCanvas * u_phi_miss  = new TCanvas("u_phi_miss", "u_phi_miss", 2000, 1000 );
	u_phi_miss->Divide(phi_bin_num, u_bin_num);


	for(int u_bin_set_num = 1; u_bin_set_num <= u_bin_num; u_bin_set_num++) {

		/// This is for computing the missmass in u-bin 

 		float yield_cut_u_l = u_lower_limit[u_bin_set_num-1];
 		float yield_cut_u_h = u_upper_limit[u_bin_set_num-1];

		TString yield_u_bin_cut;


		//// No Lower limits for the first u bin

		if(u_bin_set_num == 1) {
			yield_u_bin_cut.Form( u_plot_var + " <= %f", yield_cut_u_h);
		} else {
			yield_u_bin_cut.Form( u_plot_var + " > %f && " + u_plot_var +" <= %f",yield_cut_u_l, yield_cut_u_h);
		}

//		cout << yield_u_bin_cut << endl;
		
//		exit(0);



		obj_tree->Draw("phi_pq*180/3.1415926 >> phi_real", yield_u_bin_cut + "&&" + Diamond_Cut + " && " + Cointime_Cut, "goff");
		obj_tree->Draw("phi_pq*180/3.1415926 >> phi_rand", yield_u_bin_cut + "&&" + Diamond_Cut + " && (" + Rand_Cut + ")", "goff");

		obj_tree->Draw("t>>u_check", yield_u_bin_cut + "&&" + Diamond_Cut + " && " + Cointime_Cut, "goff");		
		
		u_check_h->SetOption("E");
		u_check_h->Write("check_u");

		u_phi->cd(u_bin_set_num);



///*--------------------------------------------------*/
/// This section is to compute the normalized yield and it is not needed 
//
// 	 	Double_t good_events_total = event_phi_real->Integral(0, -1); 
// 	 	Double_t rand_events_total = event_phi_rand->Integral(0, -1); 
// 
// 		yield = (float) good_events_total - ((float)rand_events_total)/random_real_ratio;;
// 
// 		if (is_run_dummy) {
// 			real_err = sqrt((float)good_events_total) * dummy_correction;
// 			rand_err = sqrt((float)rand_events_total) * dummy_correction / random_real_ratio;
// 		} else { 
// 			real_err = sqrt((float)good_events_total);
// 			rand_err = sqrt((float)rand_events_total) / random_real_ratio;
// 		}
// 	
// 		yield_err = sqrt( pow(real_err,2) + pow(rand_err, 2) );


	
		TH1F * event_phi_real_sub = (TH1F*) event_phi_real->Clone();
		event_phi_real_sub->Add(event_phi_rand, -1/random_real_ratio);
	

		for (int i = 1;  i <= event_phi_real_sub->GetNbinsX(); i++) {

			real_bin = event_phi_real->GetBinContent(i);
			rand_bin = event_phi_rand->GetBinContent(i);
	
			real_bin_err = sqrt(real_bin);
			rand_bin_err = sqrt(rand_bin) * 1/random_real_ratio;

			sub_err_2 = pow(real_bin_err,2) + pow(rand_bin_err,2);

			norm_yield = event_phi_real_sub->GetBinContent(i);

			if (norm_yield == 0) { 
				norm_yield_err = 0;
			} else {
				norm_yield_err = sub_err_2/pow(norm_yield, 2)  + errrr_temp/pow(acccc_temp, 2);
			}

			norm_yield = norm_yield/acccc_temp;
	        norm_yield = norm_yield / 1000;

			if (is_run_dummy) { 
				norm_yield = norm_yield * dummy_correction;
			}

			norm_yield_err = fabs(norm_yield) * sqrt(norm_yield_err) ;

 			event_phi_real_sub->SetBinContent(i, norm_yield);
 			event_phi_real_sub->SetBinError(i, norm_yield_err);

 		}

		event_phi_real_sub->SetOption("E");
		event_phi_real_sub->DrawClone();
		event_phi_real_sub->Write("phi_real_sub");

		Float_t yield_sum_check = 0.0;

		/*--------------------------------------------------*/
		// Bin in u-phi bins
		
 		for (int iiii = 0;  iiii < phi_bin_num; iiii++) {

			TCanvas * c_phi = new TCanvas();

			///*--------------------------------------------------*/
			// Four phi bin from 0-360, bin center: 45, 135, 225, 315

// 			float phi_cut_u_l = (iiii) * phi_stp + phi_offset;
// 			float phi_cut_u_h = (iiii + 1) * phi_stp + phi_offset;
	

			float phi_cut_u_l = (iiii) * phi_stp + phi_offset;
			float phi_cut_u_h = (iiii + 1) * phi_stp + phi_offset;
	


			///*--------------------------------------------------*/
			// Six phi bin from 0-360, bin center: 0, 60, 120, 180, 240, 300

// 			float phi_cut_u_l = (iiii - 0.5) * 360/float(phi_bin_num);
// 			float phi_cut_u_h = (iiii + 0.5) * 360/float(phi_bin_num);

			TString phi_u_bin_cut;

			/// correcting the negative angle to positive angle
			if(phi_cut_u_l < 0) {

				phi_cut_u_l = 360 + phi_cut_u_l;

				phi_u_bin_cut.Form("(phi_pq*180/3.1415926 > %f || phi_pq*180/3.1415926 <= %f)", phi_cut_u_l, phi_cut_u_h);

			} else if (phi_cut_u_h > 360) {

				phi_cut_u_h =  phi_cut_u_h - 360;

				phi_u_bin_cut.Form("(phi_pq*180/3.1415926 > %f || phi_pq*180/3.1415926 <= %f)", phi_cut_u_l, phi_cut_u_h);

			} else { 

				phi_u_bin_cut.Form("(phi_pq*180/3.1415926 > %f && phi_pq*180/3.1415926 <= %f)", phi_cut_u_l, phi_cut_u_h);

			}








//			cout << "~~~~~~   "<< phi_u_bin_cut << endl;



 			
			obj_tree->Draw("missmass >> h_miss", yield_u_bin_cut + " && " + Diamond_Cut + " && " + phi_u_bin_cut + " && " + Cointime_Cut, "goff");

			obj_tree->Draw("missmass >> h_miss_rand", yield_u_bin_cut + " && " + Diamond_Cut + " && " + phi_u_bin_cut + " && (" + Rand_Cut + ")", "goff");


// 		 	Double_t good_events_total = h_missmass->Integral(0, -1); 
// 		 	Double_t rand_events_total = h_missmass_rand->Integral(0, -1); 
// 
//  			yield = (float) good_events_total - ((float)rand_events_total)/random_real_ratio;
// 
// 			float real_err = 0.0;
// 			float rand_err = 0.0;
// 		
// 			if (is_run_dummy) {
// 				rand_err = sqrt((float)rand_events_total) * dummy_correction / random_real_ratio;
// 			} else { 
// 				rand_err = sqrt((float)rand_events_total) / random_real_ratio;
// 			}
// 		
// 		 	if (is_run_dummy) {
// 				real_err = sqrt((float)good_events_total) * dummy_correction;
// 			} else { 
// 				real_err = sqrt((float)good_events_total);
// 			}
// 		
// 			yield_err = sqrt( pow(real_err,2) + pow(rand_err, 2) );
//
// 			Double_t real_bin;
// 			Double_t rand_bin;	
// 			
// 			Double_t real_bin_err;
// 			Double_t rand_bin_err;
// 		
// 			Double_t sub_err_2;
// 		
// 			Double_t norm_yield;
// 			Double_t norm_yield_err;

			TH1F * h_missmass_sub= (TH1F*) h_missmass->Clone();
 			h_missmass_sub->Add(h_missmass_rand, -1/random_real_ratio);
		
			for (int i = 1;  i <= h_missmass_sub->GetNbinsX(); i++) {
	
				real_bin = h_missmass->GetBinContent(i);
				rand_bin = h_missmass_rand->GetBinContent(i);
		
 				real_bin_err = sqrt(real_bin);
 				rand_bin_err = sqrt(rand_bin) * 1/random_real_ratio;

				sub_err_2 = pow(real_bin_err,2) + pow(rand_bin_err,2);
	
				norm_yield = h_missmass_sub->GetBinContent(i);
				
				if (norm_yield == 0) { 
					norm_yield_err = 0;
				} else {
					norm_yield_err = sub_err_2/pow(norm_yield, 2)  + errrr_temp/pow(acccc_temp, 2);
				}
	
				norm_yield = norm_yield/acccc_temp;
	        	norm_yield = norm_yield / 1000;


				if (is_run_dummy) { 
					norm_yield = norm_yield * dummy_correction;
				}
	
				norm_yield_err = fabs(norm_yield) * sqrt(norm_yield_err) ;

	 			h_missmass_sub->SetBinContent(i, norm_yield);
	 			h_missmass_sub->SetBinError(i, norm_yield_err);
 
	 		}

			h_missmass->SetOption("E");
			h_missmass_sub->SetOption("E");
	
			c_phi->cd();

			h_missmass_sub->Draw();
			
			yield_sum_check = yield_sum_check + h_missmass_sub->Integral(0, -1);
 
 			TString missmass_str;
 			missmass_str.Form("mm_u_%i_phi_%i" , int((yield_cut_u_l+yield_cut_u_h)/2.0*1000), int((phi_cut_u_l + phi_cut_u_h)/2.0 * 10));

  			c_phi->Write(missmass_str); 
 
			Int_t pad_itt = (u_bin_set_num - 1) * phi_bin_num + iiii + 1;

			u_phi_miss->cd(pad_itt);

			c_phi->DrawClonePad();

			u_phi_miss->Update();



			delete c_phi;


 		}

//	    exit(0);

		cout << "phi bin over all check u bin " <<  u_bin_set_num << "    Yield Sum:  " << yield_sum_check << endl;
		
	}

	u_phi->Update();	
	u_phi->Write("u_phi");
	u_phi_miss->Write("u_phi_missmass");


	delete u_phi;
	delete u_phi_miss;

	delete event_phi_real;
	delete event_phi_rand;
	delete u_check_h;
 	delete h_missmass;

//	delete obj_tree;

}














/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Other knematics variables (other than missmass) u-phi bin

void Analysis_omega::Yield_u_Phi_Bin(TTree* obj_tree, TString plot_var) {

	TString plot_name_var;


//	cout << plot_var << endl;

	Double_t real_err = 0.0;
	Double_t rand_err = 0.0;

	Double_t real_bin, real_bin_err, norm_yield    ;
	Double_t rand_bin, rand_bin_err, norm_yield_err;	

	Double_t sub_err_2;

	Float_t lo_limit, hi_limit;
	Int_t bin_num;


// 	 "hsdelta",  80, -10,    10    
// 	 "hsyptar",  80, -0.04,  0.04 
// 	 "hsxptar",  80, -0.1,   0.1  
// 	 "ssdelta",  80, -20,    20   
// 	 "ssyptar",  80, -0.07,  0.07 
// 	 "ssxptar",  80, -0.05,  0.05 




	if (plot_var == "Em") {

		lo_limit = 0.4 ;
		hi_limit = 1.4;	
		bin_num = 50;	

	} else if (plot_var == "Pm"){

		lo_limit = 0.1;
		hi_limit = 1.1;
		bin_num = 50;	

	} else if (plot_var == "hsdelta"){

		lo_limit = -10;
		hi_limit =  10;
		bin_num =   40;	

	} else if (plot_var == "hsyptar"){

		lo_limit = -0.04;
		hi_limit =  0.04;
		bin_num  =  40;	

	} else if (plot_var == "hsxptar"){

		lo_limit = -0.1;
		hi_limit =  0.1;
		bin_num  =  40;	

	} else if (plot_var == "u"){

		lo_limit =  -0.1;
		hi_limit =   0.6;
		bin_num  =   50;	

	} else if (plot_var == "phi_pq"){

		lo_limit =  0.0;
		hi_limit =  6.4;
		bin_num  =  50;	

	} else {

		cout << "Variable: " << plot_var << " is not avaliable " << endl; 
		exit(0);

	}



	if(plot_var == "u") {

		plot_name_var = u_plot_var;

	} else {

		plot_name_var = plot_var;
		
	}






	TH1F* h_missmass      = new TH1F("h_" + plot_var, plot_var + "_real", bin_num,  lo_limit,  hi_limit);
	TH1F* h_missmass_rand = new TH1F("h_" + plot_var + "_rand", plot_var + "_rand", bin_num,  lo_limit,  hi_limit);

	TH1F* u_check_h       = new TH1F("u_check", "u_check", 30, 0, 0.7);

	TCanvas * u_phi_miss  = new TCanvas("u_phi_" + plot_var, "u_phi_" + plot_var, 2000, 1000 );
	u_phi_miss->Divide(phi_bin_num, u_bin_num);

	for(int u_bin_set_num = 1; u_bin_set_num <= u_bin_num; u_bin_set_num++) {

 		float yield_cut_u_l = u_lower_limit[u_bin_set_num-1];
 		float yield_cut_u_h = u_upper_limit[u_bin_set_num-1];

		TString yield_u_bin_cut;
//		yield_u_bin_cut.Form("t > %f && t <= %f",yield_cut_u_l, yield_cut_u_h);

		if(u_bin_set_num == 1) {
			yield_u_bin_cut.Form( u_plot_var + " <= %f", yield_cut_u_h);
		} else {
			yield_u_bin_cut.Form( u_plot_var + " > %f && " + u_plot_var +" <= %f",yield_cut_u_l, yield_cut_u_h);
		}




		
		Float_t yield_sum_check = 0.0;

		/*--------------------------------------------------*/
		// Bin in phi




 		for (int iiii = 0;  iiii < phi_bin_num; iiii++) {

			TCanvas * c_phi = new TCanvas();

			///*--------------------------------------------------*/
			// Four phi bin from 0-360, bin center: 45, 135, 225, 315

// 			float phi_cut_u_l = (iiii) * 360/float(phi_bin_num);
// 			float phi_cut_u_h = (iiii + 1) * 360/float(phi_bin_num);
// 
 
			float phi_cut_u_l = (iiii) * phi_stp + phi_offset;
			float phi_cut_u_h = (iiii + 1) * phi_stp + phi_offset;
	
	



			///*--------------------------------------------------*/
			// Six phi bin from 0-360, bin center: 0, 60, 120, 180, 240, 300

// 			float phi_cut_u_l = (iiii - 0.5) * 360/float(phi_bin_num);
// 			float phi_cut_u_h = (iiii + 0.5) * 360/float(phi_bin_num);

			TString phi_u_bin_cut;

			/// correcting the negative angle to positive angle
			if(phi_cut_u_l < 0) {
				phi_cut_u_l = 360 + phi_cut_u_l;
				phi_u_bin_cut.Form("(phi_pq*180/3.1415926 > %f || phi_pq*180/3.1415926 <= %f)", phi_cut_u_l, phi_cut_u_h);

			} else if (phi_cut_u_h > 360) {

				phi_cut_u_h =  phi_cut_u_h - 360;

				phi_u_bin_cut.Form("(phi_pq*180/3.1415926 > %f || phi_pq*180/3.1415926 <= %f)", phi_cut_u_l, phi_cut_u_h);
			} else { 

				phi_u_bin_cut.Form("(phi_pq*180/3.1415926 > %f && phi_pq*180/3.1415926 <= %f)", phi_cut_u_l, phi_cut_u_h);

			}
 			
			obj_tree->Draw(plot_name_var + " >> h_" + plot_var, yield_u_bin_cut + " && " + Diamond_Cut + " && " + phi_u_bin_cut + " && " + Cointime_Cut, "goff");
			obj_tree->Draw(plot_name_var + " >> h_" + plot_var + "_rand", yield_u_bin_cut + " && " + Diamond_Cut + " && " + phi_u_bin_cut + " && (" + Rand_Cut + ")", "goff");

			TH1F * h_missmass_sub= (TH1F*) h_missmass->Clone();
 			h_missmass_sub->Add(h_missmass_rand, -1/random_real_ratio);

			for (int i = 1;  i <= h_missmass_sub->GetNbinsX(); i++) {

				real_bin = h_missmass->GetBinContent(i);
				rand_bin = h_missmass_rand->GetBinContent(i);
		
 				real_bin_err = sqrt(real_bin);
 				rand_bin_err = sqrt(rand_bin) * 1/random_real_ratio;

				sub_err_2 = pow(real_bin_err,2) + pow(rand_bin_err,2);
	
				norm_yield = h_missmass_sub->GetBinContent(i);
				
				if (norm_yield == 0) { 
					norm_yield_err = 0;
				} else {
					norm_yield_err = sub_err_2/pow(norm_yield, 2)  + errrr_temp/pow(acccc_temp, 2);
				}
	
				norm_yield = norm_yield/acccc_temp;
	        	norm_yield = norm_yield / 1000;

				if (is_run_dummy) { 
					norm_yield = norm_yield * dummy_correction;
				}
	
				norm_yield_err = fabs(norm_yield) * sqrt(norm_yield_err) ;

	 			h_missmass_sub->SetBinContent(i, norm_yield);
	 			h_missmass_sub->SetBinError(i, norm_yield_err);

	 		}

			h_missmass->SetOption("E");
			h_missmass_sub->SetOption("E");
	
			c_phi->cd();

			h_missmass_sub->Draw();
			
			yield_sum_check = yield_sum_check + h_missmass_sub->Integral(0, -1);

 			TString missmass_str;
 			missmass_str.Form(plot_var + "_u_%i_phi_%i" , int((yield_cut_u_l+yield_cut_u_h)/2.0*1000), int((phi_cut_u_l + phi_cut_u_h)/2.0 * 10));

  			c_phi->Write(missmass_str); 
 
			Int_t pad_itt = (u_bin_set_num - 1) * phi_bin_num + iiii + 1;

			u_phi_miss->cd(pad_itt);

			c_phi->DrawClonePad();

			u_phi_miss->Update();

			delete c_phi;

 		}

		cout << "phi bin over all check u bin " <<  u_bin_set_num << "    Yield Sum:  " << yield_sum_check << endl;
		
	}

	u_phi_miss->Write("u_phi_" + plot_var);


	delete u_phi_miss;

	delete u_check_h;
 	delete h_missmass;

}






/*--------------------------------------------------*/
/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Analysis_omega::Sim_Setting_Check_U(TString plot_name , Int_t bin_num, Double_t lower, Double_t upper) {


	TH1F* mm_hist = new TH1F(plot_name, plot_name, bin_num, lower, upper);
//	TH1F* mm_hist = new TH1F();

	mm_hist->SetName(plot_name);

// 	chain->Draw("-Q2 + 2*0.938272*0.938272 + 0.78265*0.78265 + t - W**2 >>" + plot_name, sim_selection);
// 	chain->Draw("Q2 - 2*0.938272*0.938272 - 0.78265*0.78265 - t + W**2 >>" + plot_name, sim_selection);

 	chain->Draw( u_plot_var + ">>" + plot_name, sim_selection);

// 	chain->Draw("-Q2 + 2.87 + 0.61 + t - W**2 >>" + plot_name, sim_selection);

//	cout << "SSS :: " << plot_name << endl;
//	cout << "SSS :: " << plot_name << endl;

	if (plot_name == "thetapq") {

		plot_name = "th_pq";

	} else if (plot_name == "phipq") {

		plot_name = "phi_pq";

	} else if (plot_name == "pmpar") {

		plot_name = "Pmpar";

	} else if (plot_name == "pmper") {

		plot_name = "Pmper";

	} else if (plot_name == "pmoop") {

		plot_name = "Pmoop";

	}

	mm_hist->Write(plot_name);

}




/*--------------------------------------------------*/

void Analysis_omega::Sim_Diamond_Plot() {




	TMultiGraph* mg = new TMultiGraph();

	TCanvas* t1 = new TCanvas();

	t1->Update();
	
//	TString sim_app_selection;
//	sim_app_selection = sim_scale_fac_str + "(" + hms_accept_cut + " && " + sos_accept_cut  + ")";
//	sim_app_selection = sim_scale_fac_str + "(" + full_sim_cut  + ")";

// 	chain->Draw("W:Q2", sim_app_selection, "P");

 	chain->Draw("W:Q2", sim_selection_no_dia, "P");
	

	if(chain->GetSelectedRows() != 0 ) {
	
 		TGraph *diamond_setting_before_tmp = (TGraph*)gPad->GetPrimitive("Graph");	
		TGraph *diamond_setting_before = (TGraph*)diamond_setting_before_tmp->Clone();
		mg->Add(diamond_setting_before);
	
	}



	TCanvas* t2 = new TCanvas();

	// sim_app_selection = sim_scale_fac_str + "(" + hms_accept_cut + " && " + sos_accept_cut + " && " + Diamond_Cut + ")";
	// sim_app_selection = sim_scale_fac_str + "(" + full_sim_cut + " && " + Diamond_Cut + ")";

//	cout << sim_selection << endl;

//	chain->Draw("W:Q2", sim_app_selection, "P");

 	chain->Draw("W:Q2", sim_selection, "P");

	if(chain->GetSelectedRows() != 0 ) {

  		// TGraph *diamond_setting_after = (TGraph*)gPad->GetPrimitive("Graph");	
		TGraph *diamond_setting_after_tmp = (TGraph*)gPad->GetPrimitive("Graph");	
		TGraph *diamond_setting_after = (TGraph*)diamond_setting_after_tmp->Clone();

		diamond_setting_after->SetMarkerColor(2);

 		mg->Add(diamond_setting_after);
	
	}
	
	TCanvas* t3 = new TCanvas();

	t3->cd();

	t3->Clear();	

	t3->cd();

  	mg->Draw("AP");
// 
  	diamond_cut->SetLineColor(2); 
  	diamond_cut->Draw("Lsame");
//  
//   	t1->Update();
   	t3->Write("diamond_after_cut");


//	delete diamond_cut;

	delete mg;
	delete t1;
	delete t2;
	delete t3;


}




/*--------------------------------------------------*/

void Analysis_omega::Sim_2_D_Plot(TString var_1, TString var_2) {

	TString plot_name_str = var_1 + "_" + var_2;

	if ( var_1 == "u") {

		var_1 = u_plot_var;

	} else if (var_1 == "t") { 

		var_1 = t_plot_var;

	} 

	if ( var_2 == "u") {

		var_2 = u_plot_var;

	} else if (var_2 == "t") { 

		var_2 = t_plot_var;

	}


	TString plot_var = var_1 + ":"+ var_2;




	TMultiGraph* mg = new TMultiGraph();

	TCanvas* t1 = new TCanvas();

	t1->Update();
	
 	chain->Draw(plot_var, sim_selection_no_dia, "P");

	if (chain->GetSelectedRows() !=0) { 
	
  		TGraph *diamond_setting_before = (TGraph*)gPad->GetPrimitive("Graph");	
		mg->Add(diamond_setting_before);

	}

	TCanvas* t2 = new TCanvas();

 	chain->Draw(plot_var, sim_selection, "P");


	if (chain->GetSelectedRows() !=0) { 

  		TGraph *diamond_setting_after = (TGraph*)gPad->GetPrimitive("Graph");	
		
		diamond_setting_after->SetMarkerColor(2);

 		mg->Add(diamond_setting_after);

	}

	TCanvas* t3 = new TCanvas();

	t3->cd();

	t3->Clear();	

  	mg->Draw("AP");
// 
  	diamond_cut->SetLineColor(2); 
  	diamond_cut->Draw("Lsame");
//  
//   	t1->Update();
   	t3->Write(plot_name_str);


//	delete diamond_cut;

	delete mg;
	delete t1;
	delete t2;
	delete t3;


}



/*--------------------------------------------------*/

void Analysis_omega::U_Distribution_Check(TTree* obj_tree, TString plot_var, Int_t bin_num, Double_t lower, Double_t upper) {

	TString plot_name_var;

	if (plot_var == "u") {

		plot_name_var = u_plot_var;

	} else if (plot_var == "t") { 

		plot_name_var = t_plot_var;

	} else {

		plot_name_var = plot_var;

	}


	Double_t real_err = 0.0;
	Double_t rand_err = 0.0;

	Double_t real_bin, real_bin_err, norm_yield    ;
	Double_t rand_bin, rand_bin_err, norm_yield_err;	

	Double_t sub_err_2;



	TH1F* event_phi_real  = new TH1F(plot_var + "_real_var", "phi_real", bin_num, lower, upper);
	TH1F* event_phi_rand  = new TH1F(plot_var + "_rand_var", "phi_rand", bin_num, lower, upper);

	TCanvas * u_phi_var  = new TCanvas("u_" + plot_var + "_var" ,"u_phi_var" , 1600, 300);
	u_phi_var->Divide(u_bin_num, 1);

	for(int u_bin_set_num = 1; u_bin_set_num <= u_bin_num; u_bin_set_num++) {

 		float yield_cut_u_l = u_lower_limit[u_bin_set_num-1];
 		float yield_cut_u_h = u_upper_limit[u_bin_set_num-1];

		TString yield_u_bin_cut;
//		yield_u_bin_cut.Form("t > %f && t <= %f",yield_cut_u_l, yield_cut_u_h);

		if(u_bin_set_num == 1) {
			yield_u_bin_cut.Form( u_plot_var + " <= %f", yield_cut_u_h);
		} else {
			yield_u_bin_cut.Form( u_plot_var + " > %f && " + u_plot_var +" <= %f",yield_cut_u_l, yield_cut_u_h);
		}





		obj_tree->Draw(plot_name_var + " >> " + plot_var + "_real_var", yield_u_bin_cut + "&&" + Diamond_Cut + " && " + Cointime_Cut, "goff");
		obj_tree->Draw(plot_name_var + " >> " + plot_var + "_rand_var", yield_u_bin_cut + "&&" + Diamond_Cut + " && (" + Rand_Cut + ")", "goff");
	
	
		TH1F * event_phi_real_sub= (TH1F*) event_phi_real->Clone();
		event_phi_real_sub->Add(event_phi_rand, -1/random_real_ratio);

	
		for (int i = 1;  i <= event_phi_real_sub->GetNbinsX(); i++) {

			real_bin = event_phi_real->GetBinContent(i);
			rand_bin = event_phi_rand->GetBinContent(i);
	
			real_bin_err = sqrt(real_bin);
			rand_bin_err = sqrt(rand_bin) * 1/random_real_ratio;

			sub_err_2 = pow(real_bin_err,2) + pow(rand_bin_err,2);

			norm_yield = event_phi_real_sub->GetBinContent(i);

			if (norm_yield == 0) { 
				norm_yield_err = 0;
			} else {
				norm_yield_err = sub_err_2/pow(norm_yield, 2)  + errrr_temp/pow(acccc_temp, 2);
			}

			norm_yield = norm_yield/acccc_temp;
	        norm_yield = norm_yield / 1000;

			if (is_run_dummy) { 
				norm_yield = norm_yield * dummy_correction;
			}

			norm_yield_err = fabs(norm_yield) * sqrt(norm_yield_err) ;

 			event_phi_real_sub->SetBinContent(i, norm_yield);
 			event_phi_real_sub->SetBinError(i, norm_yield_err);

 		}


		
		cout << " U Bin    " << u_bin_set_num << "    Yield Check:        " << event_phi_real_sub->Integral() << endl;


		u_phi_var->cd(u_bin_set_num);

		event_phi_real_sub->SetOption("E");
		event_phi_real_sub->Draw();

		
	}

	u_phi_var->Write(plot_var + "_real_sub");


	delete event_phi_real;
	delete event_phi_rand;

	delete u_phi_var;

}








/*--------------------------------------------------*/

void Analysis_omega::U_Distribution_Check_Sim(TTree* obj_tree, TString plot_var, Int_t bin_num, Double_t lower, Double_t upper) {

	TString plot_name_var;

	if (plot_var == "u") {

		plot_name_var = u_plot_var;

	} else if (plot_var == "t") { 

		plot_name_var = t_plot_var;

	} else if (plot_var == "phi_pq" ) {

		plot_name_var = "phipq";

	} else {

		plot_name_var = plot_var;

	}


	TH1F* event_phi_real_var     = new TH1F(plot_var + "_real_var",     "phi_real_var", bin_num, lower, upper);
	TH1F* event_phi_real_var_raw = (TH1F*) event_phi_real_var->Clone(plot_var + "_real_var_raw");


	TCanvas * u_phi_var  = new TCanvas("u_" + plot_var + "_var" ,"u_phi_var" , 1600, 300);
	u_phi_var->Divide(u_bin_num, 1);


	for(int u_bin_set_num = 1; u_bin_set_num <= u_bin_num; u_bin_set_num++) {

 		float yield_cut_u_l = u_lower_limit[u_bin_set_num-1];
 		float yield_cut_u_h = u_upper_limit[u_bin_set_num-1];

		TString yield_u_bin_cut;

//		yield_u_bin_cut.Form("t > %f && t <= %f",yield_cut_u_l, yield_cut_u_h);
//		yield_u_bin_cut.Form("(Q2 - 2*0.938272*0.938272 - 0.78265*0.78265 - t + W**2) > %f && (Q2 - 2*0.938272*0.938272 - 0.78265*0.78265 - t + W**2) <= %f",yield_cut_u_l, yield_cut_u_h);


///		yield_u_bin_cut.Form(u_plot_var + " > %f && " + u_plot_var + " <= %f",yield_cut_u_l, yield_cut_u_h);

		//// No Lower limits for the first u bin
		
		if(u_bin_set_num == 1) {
			yield_u_bin_cut.Form( u_plot_var + " <= %f", yield_cut_u_h);
		} else {
			yield_u_bin_cut.Form( u_plot_var + " > %f && " + u_plot_var + " <= %f",yield_cut_u_l, yield_cut_u_h);
		}





//		cout << "asdasdasd   " << u_bin_num << "    " << yield_cut_u_l << "    "  << yield_u_bin_cut + "&&" + Diamond_Cut << endl;

//		exit(0);

		TString sim_app_selection;
//		sim_app_selection.Form("%f*Weight*", normfac/event);
//		sim_app_selection = sim_app_selection + "(" + hms_accept_cut + " && " + sos_accept_cut + " && " + yield_u_bin_cut + " && " + Diamond_Cut + ")";

		sim_app_selection = sim_scale_fac_str + "(" + full_sim_cut + " && " + yield_u_bin_cut + " && " + Diamond_Cut + ")";

		chain->Draw(plot_name_var + " >> " + plot_var + "_real_var",     sim_app_selection, "goff");

		sim_app_selection = full_sim_cut + " && " + yield_u_bin_cut + " && " + Diamond_Cut;

		chain->Draw(plot_name_var + " >> " + plot_var + "_real_var_raw", sim_app_selection, "goff");

		Sim_Set_Histo_Err(event_phi_real_var, event_phi_real_var_raw);

		u_phi_var->cd(u_bin_set_num);






		//event_phi_real->SetOption("E");
		event_phi_real_var->DrawClone();

	}

// 	cout << "aaaaaaaaaaaaaaaaaaa " << endl;
// 	exit(0);

	u_phi_var->Write(plot_var + "_real_sub");

	delete event_phi_real_var;
	delete event_phi_real_var_raw;
	delete u_phi_var;

}






/*--------------------------------------------------*/
/*--------------------------------------------------*/

Analysis_omega::~Analysis_omega() {

	delete data_tree_in_cl;
	
	delete chain;
	
	delete list;

	file_out_ana->Close();
//	delete file_out_ana;

	
	cout << " sadfasdfsdaf " << endl;

 	delete graph_yield_check;
 
 	delete mm_peak_check;
 
 	cout << " sadfasdfsdaf " << endl;
  
 	delete run_tree;
 
  	delete tree_out;
 
 	cout << " sadfasdfsdaf " << endl;
 
}



/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Analysis_omega::Sim_Set_Histo_Err(TH1F* hist, TH1F* hist_raw){

		for (int i = 1;  i <= hist->GetNbinsX(); i++) {
	
			Double_t raw_count  = hist_raw->GetBinContent(i);
			Double_t real_count = hist->GetBinContent(i);

			Double_t bin_err;

			if (real_count != 0.0) {
				Float_t average_weight = real_count / raw_count;
				bin_err = sqrt(raw_count) * average_weight;
			} else { 
				bin_err = 0.0;
			}

			hist->SetBinError(i, bin_err);

		}

}



/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Analysis_omega

Analysis_omega_itt::Analysis_omega_itt() {


} 






/*--------------------------------------------------*/

Analysis_omega_itt::Analysis_omega_itt(TString sim_data_type) {

	sim_ana_name = sim_data_type;
	weight_app = "_mod";

}
 





/*--------------------------------------------------*/




void Analysis_omega::Yield_u_Phi_Bin_Sim(TTree* obj_tree) { 


	TH1F* h_missmass      = new TH1F("h_miss", "missmass_real", 70,  0.4,  1.1);
	TH1F* h_missmass_raw  = (TH1F*) h_missmass->Clone("h_miss_raw");

	TH1F* event_phi_real  = new TH1F("phi_real", "phi_real", phi_bin_num, 0, 360);
	TH1F* u_check_h       = new TH1F("u_check", "u_check", 30, 0, 0.7);


	TCanvas * u_phi  = new TCanvas("u_phi" ,"u_phi" , 1600, 300);
	u_phi->Divide(u_bin_num, 1);

	TCanvas * u_phi_miss  = new TCanvas("u_phi_miss", "u_phi_miss", 2000, 1000 );
	u_phi_miss->Divide(phi_bin_num, u_bin_num);


	TString sim_app_selection;

	for(Int_t u_bin_set_num = 1; u_bin_set_num <= u_bin_num; u_bin_set_num++) {

 		float yield_cut_u_l = u_lower_limit[u_bin_set_num-1];
 		float yield_cut_u_h = u_upper_limit[u_bin_set_num-1];

		TString yield_u_bin_cut;
//		yield_u_bin_cut.Form("t > %f && t <= %f",yield_cut_u_l, yield_cut_u_h);
//		yield_u_bin_cut.Form("(Q2 - 2*0.938272*0.938272 - 0.78265*0.78265 - t + W**2) > %f && (Q2 - 2*0.938272*0.938272 - 0.78265*0.78265 - t + W**2) <= %f",yield_cut_u_l, yield_cut_u_h);

//		yield_u_bin_cut.Form(u_plot_var + " > %f && " + u_plot_var + " <= %f",yield_cut_u_l, yield_cut_u_h);

		//// No Lower limits for the first u bin
		
		if(u_bin_set_num == 1) {
			yield_u_bin_cut.Form( u_plot_var + " <= %f", yield_cut_u_h);
		} else {
			yield_u_bin_cut.Form( u_plot_var + " > %f && " + u_plot_var + " <= %f",yield_cut_u_l, yield_cut_u_h);
		}






		cout << endl <<  endl << " ASd   asdasdasd   " << u_bin_num << "    " << yield_cut_u_l << "    "  << yield_u_bin_cut + "&&" + Diamond_Cut << endl << endl;

//		exit(0);

//		sim_app_selection.Form("%f*Weight*", normfac/event);
//		sim_app_selection = sim_app_selection + "(" + hms_accept_cut + " && " + sos_accept_cut + " && " + yield_u_bin_cut + " && " + Diamond_Cut + ")";

		sim_app_selection = sim_scale_fac_str + "(" + full_sim_cut + " && " + yield_u_bin_cut + " && " + Diamond_Cut + ")";

//		chain->Draw("(Q2 - 2*0.938272*0.938272 - 0.78265*0.78265 - t + W**2)>>u_check", yield_u_bin_cut + "&&" + Diamond_Cut, "goff");		
//		chain->Draw("(Q2 - 2*0.938272*0.938272 - 0.78265*0.78265 - t + W**2)>>u_check", sim_app_selection, "goff");		

//		sim_app_selection = full_sim_cut + " && " + yield_u_bin_cut + " && " + Diamond_Cut;

//		cout << sim_scale_fac_str << "    " << sim_app_selection << endl;

		obj_tree->Draw(u_plot_var + " >> u_check", sim_app_selection, "goff");		
//		chain->Draw( "Weight_mod >> u_check", sim_app_selection, "goff");		
		
		//u_check_h->SetOption("E");
		u_check_h->Write("check_u");

		obj_tree->Draw("phipq*180/3.1415926 >> phi_real", yield_u_bin_cut + "&&" + Diamond_Cut, "goff");
//		chain->Draw("phipq*180/3.1415926 >> phi_real", sim_app_selection, "goff");

		u_phi->cd(u_bin_set_num);
		//event_phi_real->SetOption("E");

		event_phi_real->DrawClone();
 		event_phi_real->Write("u_phi_real");

		/*--------------------------------------------------*/
		// Bin in phi

 		for (int iiii = 0;  iiii < phi_bin_num; iiii++) {

			TCanvas * c_phi = new TCanvas();

//			cout << "Phi    " << iiii << "  " << (iiii) * 360/float(phi_bin_num) << "  " << (iiii + 1) * 360/float(phi_bin_num) << endl;

			///*--------------------------------------------------*/
			// Four phi bin from 0-360, bin center: 45, 135, 225, 315

			//float phi_cut_u_l = (iiii) * 360/float(phi_bin_num);
			//float phi_cut_u_h = (iiii + 1) * 360/float(phi_bin_num);

			float phi_cut_u_l = (iiii) * phi_stp + phi_offset;
			float phi_cut_u_h = (iiii + 1) * phi_stp + phi_offset;
	

		   

	

			///*--------------------------------------------------*/
			// Six phi bin from 0-360, bin center: 0, 60, 120, 180, 240, 300

// 			float phi_cut_u_l = (iiii - 0.5) * 360/float(phi_bin_num);
// 			float phi_cut_u_h = (iiii + 0.5) * 360/float(phi_bin_num);
// 
			TString phi_u_bin_cut;

			/// correcting the negative angle to positive angle
			if(phi_cut_u_l < 0) {

				phi_cut_u_l = 360 + phi_cut_u_l;

//				cout << phi_cut_u_l << endl;				
//				exit(0);

				phi_u_bin_cut.Form("(phipq*180/3.1415926 > %f || phipq*180/3.1415926 <= %f)", phi_cut_u_l, phi_cut_u_h);

			} else if (phi_cut_u_h > 360) {

				phi_cut_u_h =  phi_cut_u_h - 360;

				phi_u_bin_cut.Form("(phipq*180/3.1415926 > %f || phipq*180/3.1415926 <= %f)", phi_cut_u_l, phi_cut_u_h);

			} else { 

				phi_u_bin_cut.Form("(phipq*180/3.1415926 > %f && phipq*180/3.1415926 <= %f)", phi_cut_u_l, phi_cut_u_h);

			}
 

//		 	TString sim_app_selection;

			sim_app_selection = sim_scale_fac_str + "(" + full_sim_cut + " && " + yield_u_bin_cut + " && " + Diamond_Cut + " && " + phi_u_bin_cut + ")";

			obj_tree->Draw("missmass >> h_miss", sim_app_selection, "goff");
			
			sim_app_selection =  full_sim_cut + " && " + yield_u_bin_cut + " && " + Diamond_Cut + " && " + phi_u_bin_cut;
			obj_tree->Draw("missmass >> h_miss_raw", sim_app_selection, "goff");

		 	Sim_Set_Histo_Err(h_missmass, h_missmass_raw);

			c_phi->cd();

			h_missmass->Draw();

 
 			TString missmass_str;
 			missmass_str.Form("mm_u_%i_phi_%i" , int((yield_cut_u_l+yield_cut_u_h)/2.0*1000), int((phi_cut_u_l + phi_cut_u_h)/2.0 * 10));

  			c_phi->Write(missmass_str); 
 
			Int_t pad_itt = (u_bin_set_num - 1) * phi_bin_num + iiii + 1;

			u_phi_miss->cd(pad_itt);

			c_phi->DrawClonePad();

			u_phi_miss->Update();

			delete c_phi;

 		}


	}
	

	u_phi->Write("u_phi");
	u_phi_miss->Write("u_phi_missmass");

	delete u_phi;
	delete u_phi_miss;

	delete event_phi_real;
	delete u_check_h;

 	delete h_missmass;
	delete h_missmass_raw;


}








/*--------------------------------------------------*/

void Analysis_omega::Yield_u_Phi_Bin_Sim(TTree* obj_tree, TString plot_var) { 

	Float_t lo_limit, hi_limit;
	Int_t bin_num;


	TString plot_name_var;

// 
// 	if (plot_name_var == "u") {
// 
// 		plot_var = u_plot_var;
// 
// 	} else if (plot_name_var == "t") { 
// 
// 		plot_var = t_plot_var;
// 
// 	} else {
// 
// 		plot_var = plot_name_var;
// 
// 	}
// 
// 




	if (plot_var == "Em") {

		lo_limit = 0.4;
		hi_limit = 1.4;		
		bin_num  = 50; 		


	} else if (plot_var == "Pm"){

		lo_limit = 0.1;
		hi_limit = 1.1;
		bin_num  = 50; 		

	} else if (plot_var == "hsdelta"){

		lo_limit = -10;
		hi_limit =  10;
		bin_num  =  40; 		

	} else if (plot_var == "hsyptar"){

		lo_limit = -0.04;
		hi_limit =  0.04;
		bin_num  =  40; 		

	} else if (plot_var == "hsxptar"){

		lo_limit = -0.1;
		hi_limit =  0.1;
		bin_num  =  40; 		

	} else if (plot_var == "u"){

		lo_limit =  -0.1;
		hi_limit =   0.6;
		bin_num  =   50;	

	} else if (plot_var == "phi_pq"){

		lo_limit =  0.0;
		hi_limit =  6.4;
		bin_num  =  50;	

	} else {

		cout << "Variable: " << plot_var << " is not avaliable " << endl; 
		exit(0);

	}


	if(plot_var == "u") {

		plot_name_var = u_plot_var;

	} else if (plot_var == "phi_pq" ) {

		plot_name_var = "phipq";

	} else {

		plot_name_var = plot_var;
		
	}


	TH1F* h_missmass      = new TH1F("h_" + plot_var, plot_var + "_real", bin_num, lo_limit, hi_limit);
	TH1F* h_missmass_raw  = (TH1F*) h_missmass->Clone("h_" + plot_var + "_raw");

//	TH1F* event_phi_real  = new TH1F("phi_real", "phi_real", phi_bin_num, 0, 360);
//	TH1F* u_check_h       = new TH1F("u_check", "u_check", 30, 0, 0.7);


//	TCanvas * u_phi  = new TCanvas("u_phi" ,"u_phi" , 1600, 300);
//	u_phi->Divide(u_bin_num, 1);

	TCanvas * u_phi_miss  = new TCanvas("u_phi_" + plot_var, "u_phi_miss", 2000, 1000 );
	u_phi_miss->Divide(phi_bin_num, u_bin_num);


	TString sim_app_selection;

	for(Int_t u_bin_set_num = 1; u_bin_set_num <= u_bin_num; u_bin_set_num++) {

 		float yield_cut_u_l = u_lower_limit[u_bin_set_num-1];
 		float yield_cut_u_h = u_upper_limit[u_bin_set_num-1];

		TString yield_u_bin_cut;
//		yield_u_bin_cut.Form("t > %f && t <= %f",yield_cut_u_l, yield_cut_u_h);
//		yield_u_bin_cut.Form("(Q2 - 2*0.938272*0.938272 - 0.78265*0.78265 - t + W**2) > %f && (Q2 - 2*0.938272*0.938272 - 0.78265*0.78265 - t + W**2) <= %f",yield_cut_u_l, yield_cut_u_h);

//		yield_u_bin_cut.Form(u_plot_var + " > %f && " + u_plot_var + " <= %f",yield_cut_u_l, yield_cut_u_h);


		/*--------------------------------------------------*/
		//// No Lower limits for the first u bin
		
		if(u_bin_set_num == 1) {
			yield_u_bin_cut.Form( u_plot_var + " <= %f", yield_cut_u_h);
		} else {
			yield_u_bin_cut.Form( u_plot_var + " > %f && " + u_plot_var + " <= %f",yield_cut_u_l, yield_cut_u_h);
		}







		cout << endl <<  endl << " ASd   asdasdasd   " << u_bin_num << "    " << yield_cut_u_l << "    "  << yield_u_bin_cut + "&&" + Diamond_Cut << endl << endl;

//		exit(0);

//		sim_app_selection.Form("%f*Weight*", normfac/event);
//		sim_app_selection = sim_app_selection + "(" + hms_accept_cut + " && " + sos_accept_cut + " && " + yield_u_bin_cut + " && " + Diamond_Cut + ")";

		sim_app_selection = sim_scale_fac_str + "(" + full_sim_cut + " && " + yield_u_bin_cut + " && " + Diamond_Cut + ")";

//		chain->Draw("(Q2 - 2*0.938272*0.938272 - 0.78265*0.78265 - t + W**2)>>u_check", yield_u_bin_cut + "&&" + Diamond_Cut, "goff");		
//		chain->Draw("(Q2 - 2*0.938272*0.938272 - 0.78265*0.78265 - t + W**2)>>u_check", sim_app_selection, "goff");		

//		sim_app_selection = full_sim_cut + " && " + yield_u_bin_cut + " && " + Diamond_Cut;

//		cout << sim_scale_fac_str << "    " << sim_app_selection << endl;

//		obj_tree->Draw(u_plot_var + " >> u_check", sim_app_selection, "goff");		
//		chain->Draw( "Weight_mod >> u_check", sim_app_selection, "goff");		
		
		//u_check_h->SetOption("E");
//		u_check_h->Write("check_u");

//		obj_tree->Draw("phipq*180/3.1415926 >> phi_real", yield_u_bin_cut + "&&" + Diamond_Cut, "goff");
//		chain->Draw("phipq*180/3.1415926 >> phi_real", sim_app_selection, "goff");

// 		u_phi->cd(u_bin_set_num);
// 		//event_phi_real->SetOption("E");
// 
// 		event_phi_real->DrawClone();
//  		event_phi_real->Write("u_phi_real");
 
		/*--------------------------------------------------*/
		// Bin in phi

 		for (int iiii = 0;  iiii < phi_bin_num; iiii++) {

			TCanvas * c_phi = new TCanvas();

//			cout << "Phi    " << iiii << "  " << (iiii) * 360/float(phi_bin_num) << "  " << (iiii + 1) * 360/float(phi_bin_num) << endl;

			///*--------------------------------------------------*/
			// Four phi bin from 0-360, bin center: 45, 135, 225, 315

			// float phi_cut_u_l = (iiii) * 360/float(phi_bin_num);
			// float phi_cut_u_h = (iiii + 1) * 360/float(phi_bin_num);
	

			float phi_cut_u_l = (iiii) * phi_stp + phi_offset;
			float phi_cut_u_h = (iiii + 1) * phi_stp + phi_offset;
	


			///*--------------------------------------------------*/
			// Six phi bin from 0-360, bin center: 0, 60, 120, 180, 240, 300

// 			float phi_cut_u_l = (iiii - 0.5) * 360/float(phi_bin_num);
// 			float phi_cut_u_h = (iiii + 0.5) * 360/float(phi_bin_num);
// 
			TString phi_u_bin_cut;

			/// correcting the negative angle to positive angle
			if(phi_cut_u_l < 0) {

				phi_cut_u_l = 360 + phi_cut_u_l;

//				cout << phi_cut_u_l << endl;				
//				exit(0);

				phi_u_bin_cut.Form("(phipq*180/3.1415926 > %f || phipq*180/3.1415926 <= %f)", phi_cut_u_l, phi_cut_u_h);

			} else if (phi_cut_u_h > 360) {

				phi_cut_u_h =  phi_cut_u_h - 360;

				phi_u_bin_cut.Form("(phipq*180/3.1415926 > %f || phipq*180/3.1415926 <= %f)", phi_cut_u_l, phi_cut_u_h);

			} else { 

				phi_u_bin_cut.Form("(phipq*180/3.1415926 > %f && phipq*180/3.1415926 <= %f)", phi_cut_u_l, phi_cut_u_h);

			}
 
//		 	TString sim_app_selection;

			sim_app_selection = sim_scale_fac_str + "(" + full_sim_cut + " && " + yield_u_bin_cut + " && " + Diamond_Cut + " && " + phi_u_bin_cut + ")";

			chain->Draw( plot_name_var + " >> h_" + plot_var, sim_app_selection, "goff");
			
			sim_app_selection =  full_sim_cut + " && " + yield_u_bin_cut + " && " + Diamond_Cut + " && " + phi_u_bin_cut;
			chain->Draw( plot_name_var + " >> h_" + plot_var + "_raw", sim_app_selection, "goff");

		 	Sim_Set_Histo_Err(h_missmass, h_missmass_raw);

			c_phi->cd();

			h_missmass->Draw();

 
 			TString missmass_str;
 			missmass_str.Form("mm_u_%i_phi_%i" , int((yield_cut_u_l+yield_cut_u_h)/2.0*1000), int((phi_cut_u_l + phi_cut_u_h)/2.0 * 10));

  			c_phi->Write(missmass_str); 
 
			Int_t pad_itt = (u_bin_set_num - 1) * phi_bin_num + iiii + 1;

			u_phi_miss->cd(pad_itt);

			c_phi->DrawClonePad();

			u_phi_miss->Update();

			delete c_phi;

 		}

	}
	

//	u_phi->Write("u_phi");
	u_phi_miss->Write("u_phi_" + plot_var);

//	delete u_phi;
	delete u_phi_miss;


//	delete event_phi_real;
//	delete u_check_h; 
 	delete h_missmass;
	delete h_missmass_raw;


}























		



/*--------------------------------------------------*/
Analysis_omega_itt::~Analysis_omega_itt() {

}

