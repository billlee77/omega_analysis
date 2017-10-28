#include <iostream>

#include "analysis_heep.h"
//#include "cut.h"

using namespace std;

/*--------------------------------------------------*/

Analysis_heep::Analysis_heep() {

	chain = new TChain();

}

/*--------------------------------------------------*/

Analysis_heep::Analysis_heep(ReadFile::efficiency_profile eff_struc) {


	
	eff_file = "list.settings.heep";
	off_file = "offset.dat";
// 
// 	cout << off_file << "sdafasfs" << endl;

	eff_ana = eff_struc;
	Init();


	yield_setting     = 0.0;
	yield_setting_err = 0.0;

	chain = new TChain();
	list = new TList();

//	coin_beta_check_offset_setting = new TH2F("coin_beta_check_offset", "coin_beta_check", 150, -15, 15,  150, -1, 1);

}


/*--------------------------------------------------*/
/// Analysis for Heep 
void Analysis_heep::Heep_anna(Int_t run_itt) {

	cout << "Analyzing Heep data !" << endl;

	cout << " sdaf " << run_itt << endl;

 	run_tree = Create_File();
 



 	Para_Run_Def(run_itt);

	Load_Data_Tree();

	/// Setting all offsets and parameters

	Set_Coin_Bin(150);
	Set_Coin_Max(10);
	Set_Coin_Min(-20);

 	// Double_t array_temp[6] = {-0.01293, -0.01293, -0.01293, -0.01293, -0.01293, -0.01293};
 	Double_t array_temp[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

	Set_Expected_MM(0.0);
	Set_MM_Offset(array_temp);
     
 	coin_center = Get_Coin_Center();


	Define_Cuts();

	// Missing_Mass_Plot();


	TH2F* t_phi_check = new TH2F("t_phi_check", "t_phi_check", 65, 0, 6.5,  50, 0, 0.5);
	data_tree_in->Draw("t:phi_pq>>t_phi_check", all_coin_cut, "goff");

	TH2F* coin_beta_check = new TH2F("coin_beta_check", "coin_beta_check", 150, -20, 10,  150, -0.5, 1.8);
//	data_tree_in->Draw("hsbeta:cointime>>coin_beta_check", all_coin_cut, "goff");

	data_tree_in->Draw("hsbeta:cointime>>coin_beta_check", all_coin_cut, "goff");


	coin_center_off = coin_beta_check->GetMean(1);
	beta_center_off = coin_beta_check->GetMean(2);


// 	beta_cut                   = 0.1;
// 	gradi                      = 1./9.;	
// 	global_centered_hsbeta_cut = -0.15;


//	beta_cut                   = 0.1;
//	beta_cut                   = beta_center_off - 0.9;


///*--------------------------------------------------*/
/// Old
// 	gradi                      = 1./7.;	
// 	global_centered_hsbeta_cut = -0.1;
// 

	gradi                      = 1./6.;	
	global_centered_hsbeta_cut = -(beta_center_off - 0.9);


	rand_early_l = 5;
	rand_early_r = 3;
	rand_late_l  = 3;
	rand_late_r  = 9;  

	cointime_cut_l = 1.05;  
	cointime_cut_r = 0.95;  

	top_limit = 0.7;
	bot_limit = 0.6;

	random_real_ratio = (abs(rand_late_l - rand_late_r) + abs(rand_early_l-rand_early_r))/(cointime_cut_l + cointime_cut_r);

	/// Target thickness difference between the dummy target and Loop 1 hydrogen target
	dummy_correction = 0.142;

	out_dir->cd();

//	Calculate_Yield(coin_beta_check, global_centered_hsbeta_cut);


	TH2F* coin_beta_check_offset = new TH2F("coin_beta_check_offset", "coin_beta_check", 150, -15, 15,  150, -0.5, 1.8);

	TString offset_str;
	

	if (coin_center_off < 0){
		offset_str.Form("hsbeta:cointime+%f", abs(coin_center_off));
	} else { 
		offset_str.Form("hsbeta:cointime-%f", coin_center_off);
	} 

	cout << "Offset Check: " << coin_center_off << "     " << beta_center_off << "     " << offset_str << endl;

	data_tree_in->Draw(offset_str + " >> coin_beta_check_offset", all_coin_cut, "goff");

	hsbeta_center = beta_center_off;

	Calculate_Yield_Direct();

//	Calculate_Yield(coin_beta_check_offset, global_centered_hsbeta_cut);


//	cout << coin_bin <<  "       " << coin_max << "    " << coin_min << endl;

// 	cout << "Now analyzing " << target << " Target Heep Run #" << run_num << ": coin center=" 
// 		 << coin_center << endl;

// 	coin_all  = new TH1F("all", "all", coin_bin, coin_min, coin_max);
// 
// 	cout << Get_kset() << endl;
// 
// 	coin_center = Get_Coin_Center();
// 
// 	all_coin_cut = HMS_cut() + " && "+ SOS_cut() + " && " + PID_cut() + " && " + Cointime_all(coin_center) + " && " + Diamond_cut() + " && " + Missingmass_cut(miss_mass_offset[kset]) + " && " + Set_t_limit(t_min, t_max);
// 
// 	cout << "___________________" << endl << endl;
// 
// 	cout << coin_center << "   " << Cointime_all(coin_center) << endl;
// 
// 	all_coin_cut = HMS_cut() + " && "+ SOS_cut() + " && " + PID_cut() + " && "  + Cointime_all(coin_center);
// 
// 	cout << "why this is not working? " << all_coin_cut << endl;
// 
// 	data_tree_in->Draw("cointime>>all",  all_coin_cut, "goff");
// 
// 	coin_all->Write("aaaa");
// 
// 	Para_Run_Clean();
// 



	/*--------------------------------------------------*/
	/// Plot to check: missmass, W, t, EXP/SIM, Q2, th_pq,
	/// 			   phi_pq, PmPar, PmPer, PmOop, hsdelta, hsytar,
	///				   hsyptar, hsxptar, ssdelta, ssytar, ssyptar, ssxptar


 	out_dir->cd();

	Check_Plot_Out("missmass");
	Check_Plot_Out("hsdelta");
	Check_Plot_Out("W");
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

//	Check_Plot_Out("hsbeta", "cointime");



	t_phi_check->Write("t_phi_check");

	coin_beta_check->Write("coin_beta_check");




	// coin_beta_check_offset->Write("coin_beta_check_offset");


	TCanvas* d1 = new TCanvas();

	coin_beta_check_offset->Draw();

	coin_cut->SetLineColor(2); 
 	coin_cut->Draw("Lsame"); 

 	rand_early_cut->SetLineColor(6); 
 	rand_early_cut->Draw("Lsame");
	
 	rand_late_cut->SetLineColor(4); 
 	rand_late_cut->Draw("Lsame");


	d1->Write("coin_check");
	delete d1;


	Correct_Offset_Tree();

	


	list->Add(data_tree_in_cl);


//	Check_Plot_Out("ct_corrected");


	delete t_phi_check;
	delete coin_beta_check;

}


/*--------------------------------------------------*/
/// Dufine cut

void Analysis_heep::Define_Cuts() {

//	all_coin_cut = HMS_cut() + " && "+ SOS_cut()  + " && " + Heep_PID_Cut() + " && " + Cointime_all(coin_center) + "  && "+ Missingmass_heep_cut();


//	cout << all_coin_cut << endl;

//	all_coin_cut = HMS_cut() + " && "+ SOS_cut()  + " && " + Heep_PID_Cut() + " && " + Cointime_all(coin_center);

//	all_coin_cut = HMS_cut() + " && "+ SOS_cut()  + " && " + Heep_PID_Cut() + " && " + Cointime_all(coin_center) + " && " + T_heep_cut();
//	all_coin_cut = HMS_cut() + " && "+ SOS_cut()  + " && " + Heep_PID_Cut();


//	all_coin_cut = HMS_cut() + " && "+ SOS_cut()  + " && " + Heep_PID_Cut() + " && "+ Missingmass_heep_cut();


	hms_accept_cut = Omega_HMS_Cut();
	sos_accept_cut = Omega_SOS_Cut();
//	heep_pid_cut  = Omega_PID_Cut_Garth();
//
//



//  	hms_accept_cut = HMS_cut();
//  	sos_accept_cut = SOS_cut();	


	float Q2_set = Get_Q2();

	heep_pid_cut  = Heep_PID_Cut(Q2_set);

	acceptence_cut = hms_accept_cut + " && " + sos_accept_cut;

	acceptence_missmass_cut = acceptence_cut + " && "+ Missingmass_heep_cut();

	all_coin_cut = acceptence_missmass_cut + " && " + heep_pid_cut;

// 	real_coin_cut = HMS_cut() + " && "+ SOS_cut() + " && " + heep_pid_cut + " && " + Missingmass_heep_cut() + " && "+ Cointime_primary_cut(coin_center);
// 	rand_coin_cut = HMS_cut() + " && "+ SOS_cut() + " && " + heep_pid_cut + " && " + Missingmass_heep_cut() + " && "+ Cointime_random_cut(coin_center);

//  real_coin_cut = acceptence_missmass_cut + " && " + heep_pid_cut + " && "+ Cointime_primary_cut(coin_center);
//  rand_coin_cut = acceptence_missmass_cut + " && " + heep_pid_cut + " && "+ Cointime_random_cut(coin_center);

 	real_coin_cut = all_coin_cut + " && "+ Cointime_primary_cut(coin_center);
 	rand_coin_cut = all_coin_cut + " && "+ Cointime_random_cut(coin_center);

}









/*--------------------------------------------------*/
/*--------------------------------------------------*/
//
Int_t Analysis_heep::Calculate_Yield_Direct() {

 	float x_mean = coin_center;
// 	float y_mean = hist_target->GetMean(2);
 	float y_mean = hsbeta_center;


	cout << coin_center << endl;
	cout << hsbeta_center << endl;

//	exit(0);	

	float y_center = y_mean;

/*--------------------------------------------------*/


	TString coin_str;

	if (coin_center < 0){
		coin_str.Form("cointime+%f", abs(coin_center));
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
		offset_str.Form("hsbeta:cointime+%f", abs(coin_center));
	} else { 
//		offset_str.Form("hsbeta-%f:cointime-%f", y_mean, coin_center);
		offset_str.Form("hsbeta:cointime+%f", abs(coin_center));
	} 

/*--------------------------------------------------*/

	TCanvas* v1 = new TCanvas();
	


	TH1F *frame = v1->DrawFrame(-15, -0.3, 15, 1.7);

// 	v1->Update();
//
//
//	data_tree_in->Draw(offset_str, "(coin_cut || rand_early_cut || rand_late_cut) && " + all_coin_cut);
	data_tree_in->Draw(offset_str, "(coin_cut || rand_early_cut || rand_late_cut) &&" + all_coin_cut, "Psame");
 

 	coin_cut->Draw("Lsame");
 	rand_early_cut->Draw("Lsame");
 	rand_late_cut->Draw("Lsame");

	v1->Update();
	v1->Write("testtest");	
 





// 	data_tree_in->Draw(offset_str, "coin_cut && " + all_coin_cut, "P");
//  	TGraph* coin_gr = (TGraph*)gPad->GetPrimitive("Graph");
// 	TGraph* coin_gr_clone = (TGraph*)coin_gr->Clone();
// 	Int_t good_events_total = coin_gr->GetN(); 
// 
// 	data_tree_in->Draw(offset_str, "rand_early_cut && " + all_coin_cut, "P");
//  	TGraph* rand_early_gr = (TGraph*)gPad->GetPrimitive("Graph");
// 	TGraph* rand_early_clone = (TGraph*)rand_early_gr->Clone();
// 	Int_t random_early_total = rand_early_gr->GetN(); 
// 
// 	data_tree_in->Draw(offset_str, "rand_late_cut && " + all_coin_cut, "P");
//  	TGraph* rand_late_gr = (TGraph*)gPad->GetPrimitive("Graph");
// 	TGraph* rand_late_clone = (TGraph*)rand_late_gr->Clone();
// 	Int_t random_late_total = rand_late_gr->GetN(); 


	data_tree_in->Draw(offset_str, "coin_cut && " + all_coin_cut, "goff");
	Float_t good_events_total = data_tree_in->GetSelectedRows(); 

	data_tree_in->Draw(offset_str, "rand_early_cut && " + all_coin_cut, "goff");
	Float_t random_early_total = data_tree_in->GetSelectedRows(); 

	data_tree_in->Draw(offset_str, "rand_late_cut && " + all_coin_cut, "goff");
	Float_t random_late_total = data_tree_in->GetSelectedRows(); 


	data_tree_in->Draw(offset_str, "coin_cor_cut && " + all_coin_cut, "goff");
	Float_t good_events_total_cor = data_tree_in->GetSelectedRows(); 

	data_tree_in->Draw(offset_str, "rand_early_cor_cut && " + all_coin_cut, "goff");
	Float_t random_early_total_cor = data_tree_in->GetSelectedRows(); 

	data_tree_in->Draw(offset_str, "rand_late_cor_cut && " + all_coin_cut, "goff");
	Float_t random_late_total_cor = data_tree_in->GetSelectedRows(); 


	/// Correct for the hsbeta=0 events

	good_events_total  = good_events_total  + good_events_total_cor;
	random_early_total = random_early_total + random_early_total_cor;
	random_late_total  = random_late_total  + random_late_total_cor;
// 

	/*--------------------------------------------------*/
	/// This part is to account the strange events that has hsbeta=0
	/// Coin window

	TString strange_hsbeta_coin_str;
	strange_hsbeta_coin_str.Form("abs(hsbeta) < %f &&  (" + coin_str + ") < %f && (" + coin_str + ") > %f", 0.1 , cointime_cut_r, -cointime_cut_l);


	cout << coin_str << "   " << cointime_cut_r + 0.1  << "   " << cointime_cut_l - 0.1 << endl;

	data_tree_in->Draw(offset_str,  strange_hsbeta_coin_str + " && " + all_coin_cut);

 	TGraph* strange_hsbeta_coin_gr = (TGraph*)gPad->GetPrimitive("Graph");
	TGraph* strange_hsbeta_coin_gr_clone = (TGraph*)strange_hsbeta_coin_gr->Clone();
	Int_t strange_hsbeta_coin_total = strange_hsbeta_coin_gr->GetN(); 

	cout << "Good event total: " << good_events_total  << endl;
	cout << "Rand early total: " << random_early_total << endl;
	cout << "Rand late total : " << random_late_total  << endl;

	cout << strange_hsbeta_coin_total << endl;

//	exit(0);


//	cout << strange_hsbeta_coin_total << "    " << strange_hsbeta_rand_early_total << "    " 
//         << strange_hsbeta_rand_late_total << endl;


//	exit(0);




///*--------------------------------------------------*/
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
/*--------------------------------------------------*/

	cout << "/*--------------------------------------------------*/" << endl;
	cout << good_events_total << endl; 
	cout << random_early_total << endl; 
	cout << random_late_total << endl; 
	cout << "/*--------------------------------------------------*/" << endl;



//	cout << "Before  " << good_events_total << "   " << dummy_correction * good_events_total << endl;
	
    if (is_run_dummy) {

		cout <<" dummy Correction " << endl;  
		good_events_total  = dummy_correction * good_events_total;
        random_early_total = dummy_correction * random_early_total;
        random_late_total  = dummy_correction * random_late_total;

	} 

//	cout << "After  " << good_events_total << endl;






	yield = (float) good_events_total - ((float)random_early_total + (float)random_late_total)/random_real_ratio;;

	cout <<  good_events_total << "  " << random_early_total << "  " << random_late_total << "  " << random_real_ratio << "  " << dummy_correction << endl;
 

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

	cout << "Yield:   " << yield  
//		 << "   Yield Error:  " << yield_err << endl 
// 		 << random_early_total << "    " << pow(sqrt(random_early_total)/4,2) << endl 
// 		 << random_late_total << "    " << pow(sqrt(random_late_total)/4,2) << endl
		 << "     Efficiency " << Get_Total_Eff()
		 << "     Charge " << Get_Charge() 
		 << "     acccc " << acccc_temp << endl << endl;


//	exit(0);


	run_tree->Fill();

	Print_out_data();

	v1->cd();


//	data_tree_in_cl = data_tree_in->CopyTree(all_coin_cut);

	data_tree_in_cl = data_tree_in->CopyTree("(coin_cut || rand_early_cut || rand_late_cut || coin_cor_cut || rand_early_cor_cut || rand_late_cor_cut) && " + all_coin_cut);

	delete v1;

}


















/*--------------------------------------------------*/
/// Calculate the efficiency

Int_t Analysis_heep::Calculate_Yield(TH2F* hist_target, Float_t hsbeta_limit) {

  	float x_mean = hist_target->GetMean(1);
  	float y_mean = hist_target->GetMean(2);

// 	float x_mean = 0.0;
// 	float y_mean = 0.0;


// 	coin_center_off = x_mean;
// 	beta_center_off = y_mean;
 
 
//	float x_mean = coin_center;
//	float y_mean = hist_target->GetMean(2);

	float y_max  = hist_target->GetYaxis()->GetXmax();
	float y_min  = hist_target->GetYaxis()->GetXmin();

	/*--------------------------------------------------*/
	Int_t good_events_beta_cut = hist_target->Integral(hist_target->GetXaxis()->FindBin(x_mean-1), 
			hist_target->GetXaxis()->FindBin(x_mean + 1), 
			hist_target->GetYaxis()->FindBin(hsbeta_limit), 
			hist_target->GetYaxis()->FindBin(y_max)); 

	Int_t random_late_beta_cut = hist_target->Integral(hist_target->GetXaxis()->FindBin(x_mean + 3), 
			hist_target->GetXaxis()->FindBin(x_mean + 9), 
			hist_target->GetYaxis()->FindBin(hsbeta_limit), 
			hist_target->GetYaxis()->FindBin(y_max));

	Int_t random_early_beta_cut = hist_target->Integral(hist_target->GetXaxis()->FindBin(x_mean -5), 
			hist_target->GetXaxis()->FindBin(x_mean - 3 ), 
			hist_target->GetYaxis()->FindBin(hsbeta_limit), 
			hist_target->GetYaxis()->FindBin(y_max)); 



	/*--------------------------------------------------*/
    /// Correct for hsbeta = 0 events
	Int_t good_events_beta_cor = hist_target->Integral(hist_target->GetXaxis()->FindBin(x_mean-1), 
			hist_target->GetXaxis()->FindBin(x_mean + 1), 
			hist_target->GetYaxis()->FindBin(-0.2), 
			hist_target->GetYaxis()->FindBin(0.2)); 

	Int_t random_late_beta_cor = hist_target->Integral(hist_target->GetXaxis()->FindBin(x_mean + 3), 
			hist_target->GetXaxis()->FindBin(x_mean + 9), 
			hist_target->GetYaxis()->FindBin(-0.2), 
			hist_target->GetYaxis()->FindBin(0.2));

	Int_t random_early_beta_cor = hist_target->Integral(hist_target->GetXaxis()->FindBin(x_mean -5), 
			hist_target->GetXaxis()->FindBin(x_mean - 3 ), 
			hist_target->GetYaxis()->FindBin(-0.2), 
			hist_target->GetYaxis()->FindBin(0.2)); 


	


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

				top_good  = Calculate_miss_x (x_mean+1, global_centered_hsbeta_cut, hist_target->GetYaxis()->GetBinCenter(j));
				down_good = Calculate_miss_x (x_mean-1, global_centered_hsbeta_cut, hist_target->GetYaxis()->GetBinCenter(j));

				top_early  = Calculate_miss_x (x_mean-3, global_centered_hsbeta_cut, hist_target->GetYaxis()->GetBinCenter(j));
				down_early = Calculate_miss_x (x_mean-5, global_centered_hsbeta_cut, hist_target->GetYaxis()->GetBinCenter(j));

				top_late  = Calculate_miss_x (x_mean+9, global_centered_hsbeta_cut, hist_target->GetYaxis()->GetBinCenter(j));
				down_late = Calculate_miss_x (x_mean+3, global_centered_hsbeta_cut, hist_target->GetYaxis()->GetBinCenter(j));


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

 


	
	cout << "Good beta box: " << good_events_beta_cut << "   Good tail: " << tail_count_good << endl;
	cout << "Rand Early: " << random_early_beta_cut << "   Rand Early: " << tail_count_early << endl;
	cout << "Rand Late : " << random_late_beta_cut << "   Rand Late : " << tail_count_late << endl;


//	exit(0);



//	cout << good_events_beta_cut << "   " << good_events_beta_cor << endl;
	 
//	exit(0);

	Int_t good_events_total  = good_events_beta_cut + tail_count_good  ;
	Int_t random_early_total = random_late_beta_cut + tail_count_late  ;
	Int_t random_late_total  = random_early_beta_cut + tail_count_early;

	/*--------------------------------------------------*/
// 	good_events_total  = good_events_total  + good_events_beta_cor;
// 	random_early_total = random_early_total + random_late_beta_cor; 
// 	random_late_total  = random_late_total  + random_early_beta_cor;

	
	Float_t cor_prime_ratio;

	cor_prime_ratio = Float_t(good_events_beta_cor - (random_late_beta_cor+random_early_beta_cor)/4.0) / Float_t((good_events_beta_cut +tail_count_good) - (random_late_beta_cut + tail_count_late + random_early_beta_cut + tail_count_early)/4.0);

	cout << "hsbeta = 0 events correction check: "<<  good_events_beta_cor - (random_late_beta_cor+random_early_beta_cor)/4.0 << "    " << (good_events_beta_cut +tail_count_good) - (random_late_beta_cut + tail_count_late + random_early_beta_cut + tail_count_early)/4.0 << "    " << cor_prime_ratio << endl; 


//	exit(0);

	/*--------------------------------------------------*/

	if (is_run_dummy) {

		cout <<" dummy Correction " << endl;  
		good_events_total  = dummy_correction * good_events_total;
        random_early_total = dummy_correction * random_early_total;
        random_late_total  = dummy_correction * random_late_total;

	} 

	



	yield = (float) good_events_total - ((float)random_early_total + (float)random_late_total)/4.;



	cout << "BB yield           " << yield << endl;

	float real_err = 0.0;
	float rand_err = 0.0;



// 	real_err = sqrt((float)good_events_total);
// 	rand_err = sqrt((float)random_early_total + (float)random_late_total)/4;


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
// 

	yield_setting_err = yield_setting_err + pow(yield_err,2);




	cout << "Yield:   " << yield  
//		 << "   Yield Error:  " << yield_err << endl 
// 		 << random_early_total << "    " << pow(sqrt(random_early_total)/4,2) << endl 
// 		 << random_late_total << "    " << pow(sqrt(random_late_total)/4,2) << endl
		 << "     Efficiency " << Get_Total_Eff()
		 << "     Charge " << Get_Charge() 
		 << "     acccc " << acccc_temp << endl << endl;


	run_tree->Fill();

	Analysis_heep::Print_out_data();

	return 0;

}





/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Analysis_heep::Print_out_data() {



	run_tree->Write("yield");
	run_tree->SetName("tree_out");

	yield = 0.0;

	yield_err = 0.0;


	data_tree_in->SetName("111");

	// cout << Get_Current_Run_File() << endl;

	cout << "asdasdasd " << Get_Current_Run_File() << endl;

 	chain->Add(Get_Current_Run_File()+"/h9500");


}




/*--------------------------------------------------*/
/*--------------------------------------------------*/

TTree* Analysis_heep::Create_File() {

 	TTree* exam = new TTree();

 	exam->Branch("yield",     &yield, "yield/F");
 	exam->Branch("yield_err", &yield_err, "yield_err/F");

	return exam;

}





/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Analysis_heep::Yield_Out() {

	file_out_ana->cd();

	TTree* setting_tree_all = TTree::MergeTrees(list);


	cout << "Check Error: " << sqrt(yield_setting_err)/yield_setting << "     " << sqrt(errrr_temp)/acccc_temp 
		 << "   " << yield_setting_err/pow(yield_setting, 2) << "      "  <<  errrr_temp/pow(acccc_temp, 2) << endl;


	cout << "This is for Sammip " << yield_setting << "    " << yield_setting_err << endl;

	yield_setting_err = yield_setting_err/pow(yield_setting, 2)  + errrr_temp/pow(acccc_temp, 2);

//	yield_setting_err = yield_setting_err/pow(yield_setting, 2)  + errrr_temp/pow(charge_tot, 2);


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


//	exit(0);


	Float_t scale_factor;

	scale_factor = 0.0;

	scale_factor = 1.0/(acccc_temp*1000);

	cout  << " WWWWWWWWWWWW !!!! : " << dummy_correction << endl;

 	if (is_run_dummy)
 		scale_factor =  dummy_correction * scale_factor;


	Set_Scale_Factor(scale_factor);


	SetYield(yield_setting);






// 	tree_out->Fill();
// 
// 	tree_out->Write("yield");
// 	tree_out->SetName("tree_out");
// 
// 	mm_setting->Write("mm_setting"); 
// 	mm_offset_setting->Write("mm_offset_setting"); 
// 
// 	cout << endl << endl;
// 
// 	run_tree->Write("yield");
// 	run_tree->SetName("tree_out");
//
//

// 
// 
// 	chain->Draw("missmass", all_coin_cut);

 	Setting_Check_Plot_Out("missmass", 80, -0.15, 0.15 );
 	Setting_Check_Plot_Out("hsdelta",  80, -10,   5    );
 	Setting_Check_Plot_Out("W",        80,  0.6,  1.8  );
 	Setting_Check_Plot_Out("t",        80, -1,    1    );
 	Setting_Check_Plot_Out("Q2",       80,  3,    8    );
 	Setting_Check_Plot_Out("th_pq",    80, -0.02, 0.1  );
 	Setting_Check_Plot_Out("phi_pq",   80,  0,    6.5  );
 	Setting_Check_Plot_Out("Pmpar",    80, -0.5,  0.1  );
 	Setting_Check_Plot_Out("Pmper",    80, -0.1,  0.3  );
 	Setting_Check_Plot_Out("Pmoop",    80, -0.05, 0.05 );
 	Setting_Check_Plot_Out("hsytar",   80, -2.5,  2.5  );
 	Setting_Check_Plot_Out("hsyptar",  80, -0.04, 0.03 );
 	Setting_Check_Plot_Out("hsxptar",  80, -0.04, 0.03 );
 	Setting_Check_Plot_Out("ssdelta",  80, -20,   10   );
 	Setting_Check_Plot_Out("ssytar",   80, -2.5,  2    );
 	Setting_Check_Plot_Out("ssyptar",  80, -0.07, 0.07 );
 	Setting_Check_Plot_Out("ssxptar",  80, -0.05, 0.05 );
 	Setting_Check_Plot_Out("Em",       600, -0.1,  0.7  );
 	Setting_Check_Plot_Out("Pm",       800, -0.1,  0.7  );
 	Setting_Check_Plot_Out("haero_su", 300, -100.0, 200.0);
 	Setting_Check_Plot_Out("hcer_npe", 150, -5.0,  10.0);
 	Setting_Check_Plot_Out("hsbeta",   170, -0.2,  1.5 );

// 	Setting_Check_Plot_Out("hsbeta:cointime");

	Setting_Beta_Cointime_Check();



//	chain->Draw("Em:missmass", all_coin_cut, "hist");
//	TH1F* htempp = (TH1F*)gPad->GetPrimitive("htemp");
//	htempp ->Write("Em_missmass");
	
	Setting_Beta_Cointime_Check_Setting_Tree(setting_tree_all);

	Proton_Obsorption_Study();
	hsbeta_Omega_Cut_Study();


	Acceptance_check(Get_Acceptence_Vector(), heep_pid_cut, "data");



	Yield_output();

}



/*--------------------------------------------------*/
/*--------------------------------------------------*/

Float_t Analysis_heep::Calculate_miss_x (Float_t x_boundary, Float_t y_boudary, Float_t y_pos) {

	Float_t x_pos;

	x_pos = x_boundary - (y_boudary - ( y_pos )) / gradi;

	return x_pos;

}







void Analysis_heep::Check_Plot_Out(TString plot_name) {
	
 	out_dir->cd();

	data_tree_in->Draw(plot_name, all_coin_cut);

//	data_tree_in->Write(plot_name);

//	data_tree_in->Draw(plot_name, real_coin_cut);
 	TH1F *htemp_real = (TH1F*)gPad->GetPrimitive("htemp");

	if (is_run_dummy) {
		htemp_real->Scale(dummy_correction);
	}
	


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
// 	TH1F *htemp_rand = (TH1F*)gPad->GetPrimitive("htemp");
// 
// //	htemp_real->SetName(plot_name);
// //	htemp_real->Write(plot_name);
//
//
//


// 	Double_t real_bin;
// 	
// 	Double_t real_bin_err;
// 
// 	Double_t sub_err_2;
// 
// 	Double_t norm_yield;
// 	Double_t norm_yield_err;
// 
//  	for (int i = 1;  i <= htemp_real->GetNbinsX(); i++) {
//  		
// // 			real_bin = htemp_real->GetBinContent(i);
// // 	
// // 			Double_t real_bin_err = sqrt(real_bin);
//  	
// 
// // 			/*--------------------------------------------------*/
// // 			/*--------------------------------------------------*/
// // 
// 			norm_yield = htemp_real->GetBinContent(i);
//   			norm_yield_err = sub_err_2;
// // 
// //  			if (norm_yield == 0) {  
// //  				norm_yield_err = 0; 
// //  			} else {
// // //  				norm_yield_err = sub_err_2/pow(norm_yield, 2)  + errrr_temp/pow(acccc_temp, 2);
// //   				norm_yield_err = sub_err_2/pow(norm_yield, 2);
// //   			}
// 
// //			norm_yield = norm_yield;
// 
// //			htemp_real->SetBinError(i, norm_yield_err);
//  
//   	}

	htemp_real->SetOption("E");

	htemp_real->Write(plot_name);

 
}

/*--------------------------------------------------*/
/*--------------------------------------------------*/
// 
// void Analysis_heep::Check_Plot_Out(TString plot_name) {
// 	
//  	out_dir->cd();
// 	
// 	data_tree_in->Draw(plot_name, all_coin_cut);
// 
// //	data_tree_in->Write(plot_name);
// 
// //	data_tree_in->Draw(plot_name, real_coin_cut);
//  	TH1F *htemp_real = (TH1F*)gPad->GetPrimitive("htemp");
// // 	
// // 	htemp_real->SetName(plot_name);
//  	htemp_real->Write(plot_name);
//  
// }

/*--------------------------------------------------*/
/*--------------------------------------------------*/




// void Analysis_heep::Check_Plot_Out(TString plot_var_1, TString plot_var_2) {
// 	
//  	out_dir->cd();
// 	
// 	data_tree_in->Draw(plot_var_1 + ":" + plot_var_2, all_coin_cut, "P");
// 
// //	data_tree_in->Write(plot_name);
// 
// //	data_tree_in->Draw(plot_name, real_coin_cut);
//  	TGraph *gtemp = (TGraph*)gPad->GetPrimitive("Graph");
// // 	
// // 	htemp_real->SetName(plot_name);
// // 	gtemp->Draw("AP");
//  	gtemp->Write(plot_var_1 + "_vs_" + plot_var_2);
// // 
// // 	Double_t min  = htemp_real->GetXaxis()->GetXmin();
// // 	Double_t max  = htemp_real->GetXaxis()->GetXmax();
// // 	Int_t bin_num = htemp_real->GetXaxis()->GetNbins();
// // 
// // 	// cout << min << "   " << max << "   " << bin_num << endl;
// // 
// // 	data_tree_in->Draw(plot_name, rand_coin_cut);
// // 	TH1F *htemp_rand = (TH1F*)gPad->GetPrimitive("htemp");
// // 
// // //	htemp_real->SetName(plot_name);
// // //	htemp_real->Write(plot_name);
// // 
// }















/*--------------------------------------------------*/
/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Analysing Heep Simulation Data


void Analysis_heep::Heep_sim_anna(Int_t run_itt) {




 	TString sim_file;
 	sim_file = "list_sim_normfac.dat"; 

//	cout << setting_str << endl;
 
//	ReadFile* rff_sim = new ReadFile("list_sim_normfac.dat");
	ReadFile* rff_sim = new ReadFile(sim_file);
	sim_ana = rff_sim->sim_pro;

// 	cout << sim_ana.event_num[run_itt] << endl;

	Int_t q2_f, q2_r;
	Int_t ebeam_f, ebeam_r;

	q2_f    = floor(sim_ana.q2_setting[run_itt]);
//	q2_r    = sim_ana.q2_setting[run_itt]*1000. - float(q2_f) *1000.;

	q2_r    = round(sim_ana.q2_setting[run_itt]*1000.) - q2_f *1000.;

	ebeam_f = floor(sim_ana.ebeam[run_itt]);
	ebeam_r = round(sim_ana.ebeam[run_itt]*1000.) - ebeam_f *1000.;
//	ebeam_r = sim_ana.ebeam[run_itt]*1000. - float(ebeam_f) *1000.;

	cout <<"MMMMM " << ebeam_r%10 << endl;

	if (ebeam_r%10==0) {

		ebeam_r = ebeam_r/10;
	}

	TString sim_in_file_name;
	sim_in_file_name.Form("heep_ebeam_%ip%i_q2_%ip%i", ebeam_f, ebeam_r, q2_f, q2_r);
	

	cout << sim_in_file_name << endl;
	
 	Int_t event = sim_ana.event_num[run_itt];
 	Double_t normfac = sim_ana.normfac[run_itt];

//	Int_t event = 5000;
//	Double_t normfac = 5000;

//	TString run_number("heep_ebeam_4p21_q2_2p406");
	TString run_number;

	run_number = sim_in_file_name;

	data_file_dir = "sim_data/";
	data_file = data_file_dir + run_number + ".root";

	if (!rff_sim->File_exist(data_file)) return;

	file_in = TFile::Open(data_file);

	Define_Cuts();

	data_tree_in = (TTree*)file_in->Get("h666");;

	TH1F* mm = new TH1F("mm", run_number + " mm", 100, -0.1, 0.1);
	TH1F* mm_org = new TH1F("mm_org", run_number + " mm", 100, -0.1, 0.1);
	TH1F* mm_try = new TH1F("mm_try", run_number + " mm", 100, -0.1, 0.1);

	sim_selection.Form("%f*Weight*", normfac/event);
	
	sim_selection = sim_selection + "(" + acceptence_missmass_cut + ")";

//	data_tree_in->Draw("missmass>>mm", sim_selection + "(" + acceptence_missmass_cut + ")" , "goff");
	data_tree_in->Draw("missmass>>mm", sim_selection, "goff");

	data_tree_in->Draw("missmass>>mm_org", acceptence_missmass_cut, "goff");
	data_tree_in->Draw("missmass>>mm_try", "Weight * (" + acceptence_missmass_cut +")", "goff");

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


	Sim_Setting_Check_Plot_Out("missmass", 80, -0.15, 0.15 );
 	Sim_Setting_Check_Plot_Out("hsdelta",  80, -10,   5    );
 	Sim_Setting_Check_Plot_Out("W",        80,  0.6,  1.8  );
 	Sim_Setting_Check_Plot_Out("t",        80, -1,    1    );
 	Sim_Setting_Check_Plot_Out("Q2",       80,  3,    8    );
 	Sim_Setting_Check_Plot_Out("thetapq",  80, -0.02, 0.1  );
 	Sim_Setting_Check_Plot_Out("phipq",    80,  0,    6.5  );
	Sim_Setting_Check_Plot_Out("pmpar",    80, -0.5,  0.1  );
	Sim_Setting_Check_Plot_Out("pmper",    80, -0.1,  0.3  );
	Sim_Setting_Check_Plot_Out("pmoop",    80, -0.05, 0.05 );
 	Sim_Setting_Check_Plot_Out("hsytar",   80, -2.5,  2.5  );
 	Sim_Setting_Check_Plot_Out("hsyptar",  80, -0.04, 0.03 );
 	Sim_Setting_Check_Plot_Out("hsxptar",  80, -0.04, 0.03 );
	Sim_Setting_Check_Plot_Out("ssdelta",  80, -20,   10   );
	Sim_Setting_Check_Plot_Out("ssytar",   80, -2.5,  2    );
 	Sim_Setting_Check_Plot_Out("ssyptar",  80, -0.07, 0.07 );
 	Sim_Setting_Check_Plot_Out("ssxptar",  80, -0.05, 0.05 );
 	Sim_Setting_Check_Plot_Out("Em",       600, -0.1,  0.7  );
 	Sim_Setting_Check_Plot_Out("Pm",       800, -0.1,  0.7  );





	chain->Draw("Em:missmass", acceptence_missmass_cut, "hist");

	TH1F* htempp = (TH1F*)gPad->GetPrimitive("htemp");

	htempp ->Write("Em_missmass");





	TString sim_app_selection;

	sim_app_selection.Form("%f*Weight*", normfac/event);
	
	sim_app_selection = sim_app_selection + "(" + hms_accept_cut + " && " + sos_accept_cut + ")";


	cout << sim_selection << endl;
	cout << sim_app_selection << endl;

//	exit(0);


	Acceptance_check(Get_Acceptence_Vector(), sim_app_selection, "sim");

	yield_setting = yield_weighted;
	yield_setting_err = yield_weighted_error;

	Yield_output();


}



void Analysis_heep::Sim_Setting_Check_Plot_Out(TString plot_name , Int_t bin_num, Double_t lower, Double_t upper) {



	TH1F* mm_hist = new TH1F(plot_name, plot_name, bin_num, lower, upper);


 	chain->Draw(plot_name + ">>" + plot_name, sim_selection);

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


// // // 	TString selection_string;
// // // 
// // // 	selection_string.Form("%f*(",1/event);
// // // 	
// // // 	selection_string = selection_string  + all_coin_cut + ")" ;
// // // 	
// // 
// // //	cout << selection_string << endl; 
// // 
// // //	chain->Draw(plot_name);
// // 	chain->Draw(plot_name, sim_selection);
// // 
// // //	data_tree_in->Draw(plot_name, real_coin_cut);
// // 
// // 	TH1F *htemp_real = (TH1F*)gPad->GetPrimitive("htemp");
// // 
// // 	htemp_real->SetName(plot_name);
// // 
// // 	htemp_real->Write(plot_name);


// 	Double_t min  = htemp_real->GetXaxis()->GetXmin();
// 	Double_t max  = htemp_real->GetXaxis()->GetXmax();
// 	Int_t bin_num = htemp_real->GetXaxis()->GetNbins();
// 
// 	// cout << min << "   " << max << "   " << bin_num << endl;
// 
// 	chain->Draw(plot_name, rand_coin_cut);
// 	TH1F *htemp_rand = (TH1F*)gPad->GetPrimitive("htemp");
// 
//	htemp_real->SetName(plot_name);
//	htemp_real->Write(plot_name);


}




/*--------------------------------------------------*/

void Analysis_heep::Setting_Beta_Cointime_Check() {


	TTree *setting_tree = TTree::MergeTrees(list);
	
	TCanvas* beta_coin_can = new TCanvas();
// 	
// 	TH2F* setting_coin_beta_check = new TH2F("setting_coin_beta_check", "coin_beta_check", 150, -20, 10,  150, 0, 2);
// 
// 	chain->Draw("hsbeta:cointime>>setting_coin_beta_check", all_coin_cut, "goff");
// 
// //	setting_coin_beta_check->Write("coin_beta_check");
// 
// 	beta_coin_can->cd();
// 
// //	setting_coin_beta_check->Draw();
// 
// 	float x_cen = setting_coin_beta_check->GetMean(1);
// 	float y_cen = setting_coin_beta_check->GetMean(2) + global_centered_hsbeta_cut;
// 
// 	setting_coin_beta_check->Draw("scat=0.5");
// 
//  	DrawGoodBox(x_cen, y_cen);
//  	DrawRandomBox(x_cen, y_cen);


	float x_cen;
	float y_cen;

	TH2F* setting_coin_beta_check_cor = new TH2F("setting_coin_beta_check_cor", "", 150, -10, 10,  150, 0, 2);

	setting_tree->Draw("hsbeta:ct_corrected>>setting_coin_beta_check_cor", all_coin_cut, "goff");

 	beta_coin_can->Update();

	x_cen = setting_coin_beta_check_cor->GetMean(1);
	y_cen = setting_coin_beta_check_cor->GetMean(2) + global_centered_hsbeta_cut;

	setting_coin_beta_check_cor->SetMarkerColor(2);

	//setting_coin_beta_check_cor->Draw("scat=0.5 same");

	TH1F *frame = beta_coin_can->DrawFrame(-15, -0.3, 15, 1.7);

	setting_coin_beta_check_cor->Draw();

//  	DrawGoodBox(x_cen, y_cen);
//  	DrawRandomBox(x_cen, y_cen);


    setting_coin_beta_check_cor-> GetYaxis()->SetTitle("hsbeta");
    setting_coin_beta_check_cor-> GetYaxis()->CenterTitle();

    setting_coin_beta_check_cor-> GetXaxis()->SetTitle("Coincidence");
    setting_coin_beta_check_cor-> GetXaxis()->CenterTitle();



// 	beta_coin_can->Modified();

 	beta_coin_can->Update();
 	beta_coin_can->Write("coin_beta_check");

	delete beta_coin_can;
//	delete setting_tree;

}


/*--------------------------------------------------*/

void Analysis_heep::Setting_Check_Plot_Out(TString plot_name , Int_t bin_num, Double_t lower, Double_t upper) {


	TString plot_name_var;


	if (plot_name == "Pmpar") {
		plot_name_var = "-Pmpar";
	} else if (plot_name == "Pmper") {
		plot_name_var = "abs(Pmper)";

// 	cout << plot_name_var << endl;
// 
// 	exit(0);



	} else {

		plot_name_var = plot_name;

	}


	TTree *setting_tree = TTree::MergeTrees(list);


 	TH1F* mm_hist_real = new TH1F(plot_name+"_real", plot_name, bin_num, lower, upper);
 
//  	chain->Draw(plot_name + ">>" + plot_name+"_real", real_coin_cut);
  	setting_tree->Draw(plot_name_var + ">>" + plot_name+"_real", real_coin_cut, "goff");
 

 	TH1F* mm_hist_rand = new TH1F(plot_name+"_rand", plot_name, bin_num, lower, upper);
 
	if (is_run_dummy) {
		mm_hist_real->Scale(dummy_correction);
		mm_hist_rand->Scale(dummy_correction);
	}




//  	chain->Draw(plot_name + ">>" + plot_name+"_rand", rand_coin_cut);
  	setting_tree->Draw(plot_name_var + ">>" + plot_name+"_rand", rand_coin_cut, "goff");


	TH1F *mm_hist_real_sub = (TH1F*) mm_hist_real->Clone();
 	mm_hist_real_sub->Add(mm_hist_rand, -0.333333333);




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
			Double_t rand_bin_err = sqrt(rand_bin) * 0.3333333333;
	
			sub_err_2 = pow(real_bin_err,2) + pow(rand_bin_err,2);




			/*--------------------------------------------------*/
			/*--------------------------------------------------*/


			norm_yield = mm_hist_real_sub->GetBinContent(i);

			if (norm_yield == 0) { 


			} else {

				norm_yield_err = sub_err_2/pow(norm_yield, 2)  + errrr_temp/pow(acccc_temp, 2);

			}



			norm_yield = norm_yield/acccc_temp;
				
			norm_yield_err =  norm_yield * sqrt(norm_yield_err);




			norm_yield = norm_yield / 1000;
			
			norm_yield_err = norm_yield_err / 1000;

//			cout  << i << "   "   <<  mm_hist_real_sub->GetBinCenter(i) << "    "<< norm_yield << endl;



 			mm_hist_real_sub->SetBinContent(i, norm_yield);

			if (norm_yield == 0.0) {
 				mm_hist_real_sub->SetBinContent(i, 0.0);
			} else {
  				mm_hist_real_sub->SetBinError(i, norm_yield_err);
			}


// 			if (!isinf(norm_yield_err)) {
//  				mm_hist_real_sub->SetBinError(i, norm_yield_err);
// 			} else {
//  				mm_hist_real_sub->SetBinError(i, 0.0);
// 			}

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


	if (plot_name == "hsbeta" ) {
//	if (plot_name == "hsbeta" || plot_name=="haero_su") {
		mm_hist_real_sub->SetOption("hist");
	}

	mm_hist_real_sub->Write(plot_name);
// 


//    exit(0);

	delete setting_tree;

	delete mm_hist_real;
	delete mm_hist_rand;
	delete mm_hist_real_sub;


}






/*--------------------------------------------------*/
//// Old function

// void Analysis_heep::Setting_Check_Plot_Out(TString plot_name , Int_t bin_num, Double_t lower, Double_t upper) {
// 
// 
// 	TString selection_string;
// 
// 	selection_string.Form("%f*(",1/(acccc_temp*1000));
// 	
// //	selection_string = selection_string  + all_coin_cut + ")" ;
// 	selection_string = selection_string  + real_coin_cut + ")" ;
// 
// 	TH1F* mm_hist_real = new TH1F(plot_name+"_real", plot_name, bin_num, lower, upper);
// 
//  	chain->Draw(plot_name + ">>" + plot_name+"_real", selection_string);
// 
// 	TString selection_string1;
// 
// 	selection_string1.Form("%f*(",1/(acccc_temp*1000));
// 	
// 	selection_string1 = selection_string1  + rand_coin_cut + ")" ;
// 
// 	TH1F* mm_hist_rand = new TH1F(plot_name+"_rand", plot_name, bin_num, lower, upper);
// 
//  	chain->Draw(plot_name + ">>" + plot_name+"_rand", selection_string1);
// 
// 	mm_hist_real->Add(mm_hist_rand, -0.333333333);
// 
// 	mm_hist_real->Write(plot_name);
// 
// }







/*--------------------------------------------------*/

void Analysis_heep::Setting_Check_Plot_Out(TString plot_name) {
	
	
	chain->Draw(plot_name, all_coin_cut);
 	TH1F *htemp = (TH1F*)gPad->GetPrimitive("htemp");

	htemp->Write(plot_name);

}















/*--------------------------------------------------*/
/// Plot boxes


void Analysis_heep::DrawGoodBox(float x_center, float y_center) {

//	cout << x_center - 1 << endl; 

	TLine *line1 = new TLine();
 	line1->SetLineWidth(1.4);
 	line1->SetLineColor(2);
	line1->DrawLine(x_center-1, y_center, x_center-1, y_center+1.15);
	line1->DrawLine(x_center+1, y_center, x_center+1, y_center+1.15);
	line1->DrawLine(x_center-1, y_center, x_center+1, y_center);

	line1->DrawLine(x_center-1, y_center, Calculate_tail_x (x_center-1, y_center ), y_center-0.6);
	line1->DrawLine(x_center+1, y_center, Calculate_tail_x (x_center+1, y_center ), y_center-0.6);


}



void Analysis_heep::DrawRandomBox(float x_center, float y_center) {

	TLine *line1 = new TLine();
 	line1->SetLineStyle(2);
 	line1->SetLineWidth(1.4);
	line1->DrawLine(x_center+3, y_center, x_center+3, y_center+1.15);
	line1->DrawLine(x_center+9, y_center, x_center+9, y_center+1.15);
	line1->DrawLine(x_center+3, y_center, x_center+9, y_center);

	line1->DrawLine(x_center+3, y_center, Calculate_tail_x (x_center+3, y_center ), y_center-0.6);
	line1->DrawLine(x_center+9, y_center, Calculate_tail_x (x_center+9, y_center ), y_center-0.6);


	TLine *line2 = new TLine();
 	line2->SetLineStyle(2);
 	line2->SetLineWidth(1.4);
	line2->DrawLine(x_center-5, y_center, x_center-5, y_center+1.15);
	line2->DrawLine(x_center-3, y_center, x_center-3, y_center+1.15);
	line2->DrawLine(x_center-5, y_center, x_center-3, y_center);

	line2->DrawLine(x_center-3, y_center, Calculate_tail_x (x_center-3, y_center), y_center-0.6);
	line2->DrawLine(x_center-5, y_center, Calculate_tail_x (x_center-5, y_center), y_center-0.6);

}



Float_t Analysis_heep::Calculate_tail_x (Float_t xx, Float_t yy) {
	
	Float_t xx_tail;

	xx_tail = xx - (yy - (yy-0.6)) / gradi;

//	cout << "Calculation Check: "  << xx_tail << "    " << xx << "    " << float (yy - (-0.6)) * gradi << endl;

	return xx_tail;


}















/*--------------------------------------------------*/

void Analysis_heep::Yield_output() {

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


	cout <<  yield_setting << "   " << yield_setting_err << endl;

	yield_file_out << yield_setting  << "   "<< yield_setting_err << endl;


//	exit(0);

	


}





/*--------------------------------------------------*/




void Analysis_heep::Acceptance_check(vector<TString> list_vec, TString cut_str, TString data_type) {


	vector<TString>::iterator it;

	TDirectory* acceptence_dir;

	acceptence_dir = file_out_ana->mkdir("acceptence_check");

	acceptence_dir->cd();

  	for (it = list_vec.begin(); it != list_vec.end(); ++it) {

		TString tmp;
		tmp = *it;

		Apt_Setting_Check(tmp, cut_str, data_type);

//    	cout << "    " << tmp << endl;

	}

//    	cout << ' ' << *it << endl;
// 
// 	exit(0);

}


/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Analysis_heep::Apt_Setting_Check(TString plot_name, TString cut_str, TString data_type) {
	
// 	ssdelta_h     = new TH1F("ssdelta", "ssdelta", 200, -20, 20);
// 	ssxptar_h     = new TH1F("ssxptar", "ssxptar", 200, -0.05, 0.05);
// 	ssyptar_h     = new TH1F("ssyptar", "ssyptar", 200, -0.1, 0.1);
// 
// 	w_h           = new TH1F("w_h",  "w_h",  200, 0.5, 1.5);
// 	q2_h          = new TH1F("q2_h", "q2_h", 200, 4, 6);
// 
//	sszbeam_h     = new TH1F("sszbeam",  "sszbeam" , 200, -5, 10);
//
//

// w_h           = new TH1F("w_h",  "w_h",  200, 0.5, 1.5)



	
	TString plot_op_str;
	

	if (plot_name == "ssxptar") {
	
		plot_op_str = plot_name + ">> htemp(200, -0.05, 0.05)";	

	} else if (plot_name == "ssyptar") {

		plot_op_str = plot_name + ">> htemp(200, -0.1, 0.1)";	
 
	} else if (plot_name == "Q2") {

		plot_op_str = plot_name + ">> htemp(200, 3, 8)";	

	} else if (plot_name == "sszbeam") {

		plot_op_str = plot_name + ">> htemp(200, -5, 10)";	

	} else if (plot_name == "ssytar") {

		plot_op_str = plot_name + ">> htemp(200, -5, 5)";	

	} else if (plot_name == "ssdelta") {

		plot_op_str = plot_name + ">> htemp(200, -40, 40)";	

	} else if (plot_name == "W") {

		plot_op_str = plot_name;

		if( data_type == "data") {

			plot_op_str = "w";

		}

 		plot_op_str = plot_op_str + " >> htemp(200, 0.5, 1.5)";

	}  else if (plot_name == "ssxfp") {

		plot_op_str = plot_name + " >> htemp(200, -20, 15)";

	}  else if (plot_name == "ssxpfp") {

		plot_op_str = plot_name + " >> htemp(200, -0.15, 0.14)";

	} else if (plot_name == "missmass") {

		plot_op_str = plot_name + " >> htemp(200, -0.15, 0.15)";

	} else if (plot_name == "hsdelta") {

		plot_op_str = plot_name + " >> htemp(200, -15, 15)";

	} else if (plot_name == "hsxptar") {

		plot_op_str = plot_name + " >> htemp(200, -0.1, 0.1)";

	} else if (plot_name == "hsyptar") {

		plot_op_str = plot_name + " >> htemp(200, -0.1, 0.1)";

	} else {

		plot_op_str = plot_name;

	}



	// cout << " asdasda    " << cut_str << endl;




 	chain->Draw(plot_op_str, cut_str);
 
   	TH1F* htempp = (TH1F*)gPad->GetPrimitive("htemp");
// 
// 
	if (data_type == "data") {
	
		htempp ->Scale(Get_Scale_Factor());

	} 
// 
// 
//

	cout << plot_name << endl;

 	htempp->Write(plot_name);

	delete htempp;


}

/*--------------------------------------------------*/

vector<TString> Analysis_heep::Get_Acceptence_Vector() {

	static const TString list_arr[] = { "missmass", "Q2","ssdelta", "ssxfp", "ssxpfp", "ssytar", "ssxptar", "ssyptar", "hsdelta", "hsxptar", "hsyptar"};

	vector<TString> vec (list_arr, list_arr + sizeof(list_arr) / sizeof(list_arr[0]) );

	return vec;

}



/*--------------------------------------------------*/

void Analysis_heep::Correct_Offset_Tree() { 

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
/// To study for the proton obsorption correction

void Analysis_heep::Proton_Obsorption_Study() {


	TTree *setting_tree_tmp = TTree::MergeTrees(list);

	Proton_Obsorption_Study(setting_tree_tmp);


	
	cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA " << endl;

	delete setting_tree_tmp;

}


/*--------------------------------------------------*/
/*--------------------------------------------------*/
void Analysis_heep::Proton_Obsorption_Study(TTree* setting_tree) {

	TString plot_name;
	plot_name = "hsbeta";

	Int_t bin_num;
	Double_t lower, upper;

	bin_num = 150;
	lower = -0.2;
	upper =  1.5;




///*--------------------------------------------------*/
///*--------------------------------------------------*/
/// Random subtraction method 

// 
//   	TH1F* mm_hist_real = new TH1F(plot_name+"_real_p", plot_name, bin_num, lower, upper);
//   
//  //  	chain->Draw(plot_name + ">>" + plot_name+"_real", real_coin_cut);
//    	setting_tree->Draw(plot_name + ">>" + plot_name+"_real_p", real_coin_cut);
//   
//  
//   	TH1F* mm_hist_rand = new TH1F(plot_name+"_rand_p", plot_name, bin_num, lower, upper);
//   
//  	if (is_run_dummy) {
//  		mm_hist_real->Scale(dummy_correction);
//  		mm_hist_rand->Scale(dummy_correction);
//  	}
//  
//  //  	chain->Draw(plot_name + ">>" + plot_name+"_rand", rand_coin_cut);
//    	setting_tree->Draw(plot_name + ">>" + plot_name+"_rand_p", rand_coin_cut);
//  
//  
//  	TH1F *mm_hist_real_sub = (TH1F*) mm_hist_real->Clone();
//   	mm_hist_real_sub->Add(mm_hist_rand, -0.333333333);
//  
//  	Double_t real_bin;
//  	Double_t rand_bin;	
//  	
//  	Double_t real_bin_err;
//  	Double_t rand_bin_err;
//  
//  	Double_t sub_err_2;
//  
//  	Double_t norm_yield;
//  	Double_t norm_yield_err;
//  
//  	for (int i = 1;  i <= mm_hist_real_sub->GetNbinsX(); i++) {
//  		
//  			real_bin = mm_hist_real->GetBinContent(i);
//  			rand_bin = mm_hist_rand->GetBinContent(i);
//  	
//  			Double_t real_bin_err = sqrt(real_bin);
//  			Double_t rand_bin_err = sqrt(rand_bin) * 0.3333333333;
//  	
//  			sub_err_2 = pow(real_bin_err,2) + pow(rand_bin_err,2);
//  
//  			/*--------------------------------------------------*/
//  			/*--------------------------------------------------*/
//  
//  
//  			norm_yield = mm_hist_real_sub->GetBinContent(i);
//  
//  			if (norm_yield == 0) { 
//  
//  
//  			} else {
//  
//  				norm_yield_err = sub_err_2/pow(norm_yield, 2)  + errrr_temp/pow(acccc_temp, 2);
//  
//  			}
//  
//  			norm_yield = norm_yield/acccc_temp;
//  				
//  			norm_yield_err =  norm_yield * sqrt(norm_yield_err);
//  
//  			norm_yield = norm_yield / 1000;
//  			
//  			norm_yield_err = norm_yield_err / 1000;
//  
//   			mm_hist_real_sub->SetBinContent(i, norm_yield);
//   			mm_hist_real_sub->SetBinError(i, norm_yield_err);
//  
//   	}











///*--------------------------------------------------*/
///*--------------------------------------------------*/
///  Non random subtraction method 

 	TH1F* mm_hist_real = new TH1F(plot_name+"_real_p", plot_name, bin_num, lower, upper);
  	setting_tree->Draw(plot_name + ">>" + plot_name+"_real_p", all_coin_cut);

 	TH1F *mm_hist_real_sub = (TH1F*) mm_hist_real->Clone();


/*--------------------------------------------------*/


	if (plot_name == "hsbeta") {
		mm_hist_real_sub->SetOption("hist");
	}

	mm_hist_real_sub->GetXaxis()->SetRangeUser(0.85, 1.05);


	TCanvas* v2 = new TCanvas();  

	mm_hist_real_sub->Draw("hist");

	Double_t hsbeta_max = mm_hist_real_sub->GetMean();	

	mm_hist_real_sub->GetXaxis()->SetRangeUser(-0.2, 1.5);

	cout << hsbeta_max << endl;


	Float_t hsbeta_0_events;
	Float_t hsbeta_l_tail;
	Float_t hsbeta_r_tail;
	Float_t hsbeta_peak;
	Float_t hsbeta_tot;

	
// 	hsbeta_0_events = mm_hist_real_sub->Integral(mm_hist_real_sub->FindBin(-0.1), mm_hist_real_sub->FindBin(0.1));
// 	hsbeta_l_tail   = mm_hist_real_sub->Integral(mm_hist_real_sub->FindBin(0.1)+1, mm_hist_real_sub->FindBin(hsbeta_max-0.11)-1);
// 	hsbeta_r_tail   = mm_hist_real_sub->Integral(mm_hist_real_sub->FindBin(hsbeta_max+0.11)+1, mm_hist_real_sub->FindBin(1.40));
// 	hsbeta_peak     = mm_hist_real_sub->Integral(mm_hist_real_sub->FindBin(hsbeta_max-0.11), mm_hist_real_sub->FindBin(hsbeta_max+0.11));
// 

	hsbeta_0_events = mm_hist_real_sub->Integral(mm_hist_real_sub->FindBin(-0.1), mm_hist_real_sub->FindBin(0.1));
	hsbeta_l_tail   = mm_hist_real_sub->Integral(mm_hist_real_sub->FindBin(0.1)+1, mm_hist_real_sub->FindBin(0.9));
	hsbeta_r_tail   = mm_hist_real_sub->Integral(mm_hist_real_sub->FindBin(hsbeta_max+0.11)+1, mm_hist_real_sub->FindBin(1.40));
	hsbeta_peak     = mm_hist_real_sub->Integral(mm_hist_real_sub->FindBin(0.9)+1, mm_hist_real_sub->FindBin(1.5));




	hsbeta_tot =  mm_hist_real_sub->Integral(0, -1);

	cout << mm_hist_real_sub->Integral(0, -1) << endl;   

	cout << hsbeta_0_events << "  " << hsbeta_l_tail << "  " << hsbeta_peak << "  " << hsbeta_r_tail << endl;

	cout << "Ratio " << Get_Q2() << "    " << file_out_ana->GetName() << "   " << hsbeta_0_events/hsbeta_tot *100 << "   " << hsbeta_l_tail/hsbeta_tot *100 << "   "<< hsbeta_peak/hsbeta_tot*100 << endl;

	mm_hist_real_sub->Write(plot_name + "_p_absorption_cor");


//	cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA " << endl;

// 	TLine* line = new TLine();
// 
// 	line->DrawLine(-0.1, 0, -0.1, 1);
// 	line->DrawLine( 0.1, 0,  0.1, 1);
// 	line->DrawLine(hsbeta_max-0.11, 0, hsbeta_max-0.11, 1);
// 	line->DrawLine(hsbeta_max+0.11, 0, hsbeta_max+0.11, 1);
// 	line->DrawLine(1.4, 0, 1.4, 1);
 


	v2->SetLogy();

	v2->Write("p_absorption");


	delete v2;

}







// /*--------------------------------------------------*/
// /*--------------------------------------------------*/
// 
// Double_t Analysis_heep::GetBetaCenter() { 
// 
// 
// 	TString coin_cut_str;
// 
// 	coin_cut_str.Form("cointime > %f && cointime <%f", coin_center - cointime_cut_l, coin_center + cointime_cut_r);
// 
// 
// 	TCanvas* v2 = new TCanvas();  
// 
// 	TH1F* hsbeta_check = new TH1F("hsbeta_check", "hsbeta",  400, 0, 2);
// 
// 	data_tree_in->Draw("hsbeta >> hsbeta_check", all_coin_cut + " && " + coin_cut_str, "goff");
// 
// 	v2->cd();
// 
// 	hsbeta_check->GetXaxis()->SetRangeUser(0.85, 1.05);
// 
// 	hsbeta_check->Draw("hist");
// 
// 	Double_t hsbeta_max = hsbeta_check->GetMaximum();	
// 
// 	v2->Update();
// 
// 	v2->Write("hsbeta_check");
// 
// 	Double_t hsbeta_center_tmp = hsbeta_check->GetMean();
// 
// 	// cout << hsbeta_center_tmp << " asdasdas " << endl;
// 
// 	return hsbeta_center_tmp;
// 
// }
// 



// /*--------------------------------------------------*/
// /*--------------------------------------------------*/
// 
// Double_t Analysis_heep::GetBetaCenterSetting() { 
// 
// 
// 	TString coin_cut_str;
// 
// 	coin_cut_str.Form("cointime > %f && cointime <%f", coin_center - cointime_cut_l, coin_center + cointime_cut_r);
// 
// 
// 	TCanvas* v2 = new TCanvas();  
// 
// 	TH1F* hsbeta_check = new TH1F("hsbeta_check", "hsbeta",  400, 0, 2);
// 
// 	data_tree_in->Draw("hsbeta >> hsbeta_check", all_coin_cut + " && " + coin_cut_str, "goff");
// 
// 	v2->cd();
// 
// 	hsbeta_check->GetXaxis()->SetRangeUser(0.85, 1.05);
// 
// 	hsbeta_check->Draw("hist");
// 
// 	Double_t hsbeta_max = hsbeta_check->GetMaximum();	
// 
// 	v2->Update();
// 
// 	v2->Write("hsbeta_check");
// 
// 	Double_t hsbeta_center_tmp = hsbeta_check->GetMean();
// 
// 	// cout << hsbeta_center_tmp << " asdasdas " << endl;
// 
// 	return hsbeta_center_tmp;
// 
// }










/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Analysis_heep::hsbeta_Omega_Cut_Study() {

	TString plot_name;
	plot_name = "haero_su";

	Int_t bin_num;
	Double_t lower, upper;

	bin_num = 500;
	lower   =  -10;
	upper   =  40;



	TString plot_name_var;
	plot_name_var = plot_name;


	TTree *setting_tree = TTree::MergeTrees(list);


 	TH1F* mm_hist_real = new TH1F(plot_name+"_real", plot_name, bin_num, lower, upper);
 
//  	chain->Draw(plot_name + ">>" + plot_name+"_real", real_coin_cut);
  	setting_tree->Draw(plot_name_var + ">>" + plot_name+"_real", real_coin_cut, "goff");
 

 	TH1F* mm_hist_rand = new TH1F(plot_name+"_rand", plot_name, bin_num, lower, upper);
 
	if (is_run_dummy) {
		mm_hist_real->Scale(dummy_correction);
		mm_hist_rand->Scale(dummy_correction);
	}




//  	chain->Draw(plot_name + ">>" + plot_name+"_rand", rand_coin_cut);
  	setting_tree->Draw(plot_name_var + ">>" + plot_name+"_rand", rand_coin_cut, "goff");


	TH1F *mm_hist_real_sub = (TH1F*) mm_hist_real->Clone();
 	mm_hist_real_sub->Add(mm_hist_rand, -0.333333333);




	Double_t real_bin;
	Double_t rand_bin;	
	
	Double_t real_bin_err;
	Double_t rand_bin_err;

	Double_t sub_err_2;

	Double_t norm_yield;
	Double_t norm_yield_err;

	for (int i = 1;  i <= mm_hist_real_sub->GetNbinsX(); i++) {
		
			real_bin = mm_hist_real->GetBinContent(i);
			rand_bin = mm_hist_rand->GetBinContent(i);
	
			Double_t real_bin_err = sqrt(real_bin);
			Double_t rand_bin_err = sqrt(rand_bin) * 0.3333333333;
	
			sub_err_2 = pow(real_bin_err,2) + pow(rand_bin_err,2);

			/*--------------------------------------------------*/
			/*--------------------------------------------------*/


			norm_yield = mm_hist_real_sub->GetBinContent(i);

			if (norm_yield == 0) { 


			} else {

				norm_yield_err = sub_err_2/pow(norm_yield, 2)  + errrr_temp/pow(acccc_temp, 2);

			}

			norm_yield = norm_yield/acccc_temp;
				
			norm_yield_err =  norm_yield * sqrt(norm_yield_err);

			norm_yield = norm_yield / 1000;
			
			norm_yield_err = norm_yield_err / 1000;

 			mm_hist_real_sub->SetBinContent(i, norm_yield);

			if (norm_yield == 0.0) {
 				mm_hist_real_sub->SetBinContent(i, 0.0);
			} else {
  				mm_hist_real_sub->SetBinError(i, norm_yield_err);
			}

 	}


	mm_hist_real_sub->SetOption("E");





	Float_t haero_cut_left, haero_cut_right, haero_cut, haero_total, haero_eff;
	Float_t haero_cut_left_err, haero_cut_right_err, haero_cut_err, haero_total_err, haero_eff_err;

	haero_cut_left  = mm_hist_real_sub->Integral(0, mm_hist_real_sub->FindBin(-2.5)-1);
	haero_cut       = mm_hist_real_sub->Integral(mm_hist_real_sub->FindBin(-2.5), mm_hist_real_sub->FindBin(2.5));
	haero_cut_right = mm_hist_real_sub->Integral(mm_hist_real_sub->FindBin(2.5)+1, -1);

	
	
	



	haero_cut_err       = 0;
	haero_cut_left_err  = 0;
	haero_cut_right_err = 0;
	haero_total_err     = 0;
	haero_eff_err       = 0;

		
	for (int i = 0;  i < mm_hist_real_sub->FindBin(-2.5)-1; i++) {

  		Float_t cont     = mm_hist_real_sub->GetBinContent(i);
  		Float_t cont_err = mm_hist_real_sub->GetBinError(i);

//		haero_cut_left_err = haero_cut_left_err + pow(cont_err/cont,2);
		haero_cut_left_err = haero_cut_left_err + pow(cont_err,2);

 	}

	haero_cut_left_err = sqrt(haero_cut_left_err);

		
	for (int i = mm_hist_real_sub->FindBin(2.5)+1;  i < mm_hist_real_sub->GetNbinsX(); i++) {

  		Float_t cont     = mm_hist_real_sub->GetBinContent(i);
  		Float_t cont_err = mm_hist_real_sub->GetBinError(i);

//		haero_cut_right_err = haero_cut_right_err + pow(cont_err/cont,2);
		haero_cut_right_err = haero_cut_right_err + pow(cont_err,2);

 	}

	haero_cut_right_err = sqrt(haero_cut_right_err);
		
	for (int i = mm_hist_real_sub->FindBin(-2.5);  i <= mm_hist_real_sub->FindBin(2.5); i++) {

  		Float_t cont     = mm_hist_real_sub->GetBinContent(i);
  		Float_t cont_err = mm_hist_real_sub->GetBinError(i);

//		haero_cut_err = haero_cut_err+ pow(cont_err/cont,2);
		haero_cut_err = haero_cut_err+ pow(cont_err,2);

 	}

	haero_cut_err = sqrt(haero_cut_err);


	haero_total = haero_cut + haero_cut_left + haero_cut_right;
	haero_total_err = sqrt(pow(haero_cut_err,2) + pow(haero_cut_left_err,2) + pow(haero_cut_right_err,2));

	haero_eff = haero_cut/haero_total;
	haero_eff_err = sqrt( pow(haero_cut_err/haero_cut,2) + pow(haero_total_err/haero_total,2) ) * haero_eff;




	TString total_event_str;
	TString total_event_l_str;
	TString total_event_r_str;
	TString total_event_cut_str;
	TString eff_str;

//	total_event.Form("Total: %f #pm %f", haero_cut_left + haero_cut + haero_cut_right, );
	total_event_cut_str.Form("Cut: %f #pm %f", haero_cut, haero_cut_err);
	total_event_l_str.Form("L: %f #pm %f", haero_cut_left, haero_cut_left_err);
	total_event_r_str.Form("R: %f #pm %f", haero_cut_right, haero_cut_right_err);
	total_event_str.Form("Total: %f #pm %f", haero_total, haero_total_err);
	eff_str.Form("Eff: %f #pm %f", haero_eff, haero_eff_err);

//	total_event_cut.Form("Efficiency: %f ", haero_cut/(haero_cut_left + haero_cut + haero_cut_right));




	TCanvas* cc1 = new TCanvas();


	cc1->SetLogy();

	mm_hist_real_sub->Draw("E"); 


	TLatex* latex = new TLatex();
	latex->SetNDC();
	latex->SetTextFont(13);
	latex->SetTextSize(20);

	latex->DrawLatex(.55, .70, total_event_cut_str);
	latex->DrawLatex(.55, .65, total_event_l_str);
	latex->DrawLatex(.55, .60, total_event_r_str);
	latex->DrawLatex(.55, .55, total_event_str);
	latex->DrawLatex(.55, .50, eff_str);

	TLine* line = new TLine();

	line->SetLineColor(2);
	line->SetLineWidth(2);

	line->DrawLine(-2.5, 0, -2.5, mm_hist_real_sub->GetMaximum()/2);
	line->DrawLine( 2.5, 0,  2.5, mm_hist_real_sub->GetMaximum()/2);

//	mm_hist_real_sub->Write("Omega_cut_study");


	mm_hist_real_sub->Write("Omega_cut_study_hist");	

	cc1->Write("Omega_cut_study");

	delete setting_tree;
	delete mm_hist_real;
	delete mm_hist_rand;
	delete mm_hist_real_sub;
	delete latex;
	delete cc1;

}

/*--------------------------------------------------*/

void Analysis_heep::Setting_Beta_Cointime_Check_Setting_Tree(TTree* obj_tree) {
	
	TCanvas* beta_coin_can = new TCanvas();

//	TTree *setting_tree = TTree::MergeTrees(list);
//  	setting_tree->Draw("hsbeta:ct_corrected", Cointime_Cut + " || " + Rand_Cut);
 

	Int_t coin_prime, coin_cor;
	Int_t rand_prime, rand_cor;
	Float_t cor_prime_ratio;

 	obj_tree->Draw("hsbeta:ct_corrected", Cointime_Cut + " && hsbeta > 0.1", "goff");
	coin_prime = obj_tree->GetSelectedRows();

 	obj_tree->Draw("hsbeta:ct_corrected", Cointime_Cut + "&& hsbeta < 0.1", "goff");
	coin_cor = obj_tree->GetSelectedRows();

 	obj_tree->Draw("hsbeta:ct_corrected", Rand_Cut + "&& hsbeta > 0.1", "goff");
	rand_prime = obj_tree->GetSelectedRows();

	obj_tree->Draw("hsbeta:ct_corrected", Rand_Cut + "&& hsbeta < 0.1", "goff");
	rand_cor = obj_tree->GetSelectedRows();

	cor_prime_ratio = Float_t(coin_cor - rand_cor/random_real_ratio) / Float_t(coin_prime - rand_prime/random_real_ratio);

// 	cout << coin_prime << "  " << rand_prime << "  " << coin_cor   << "  " << rand_cor << "  " << cor_prime_ratio << endl;
// 
// 	cout << "Rand subtracted Prime Events: "    << Float_t(coin_prime - rand_prime/random_real_ratio) << endl;
// 	cout << "Rand subtracted hsbeta=0 Events: " << Float_t(coin_cor - rand_cor/random_real_ratio)     << endl;
// 	cout << "(hsbeta=0)/Prime Ratio: " << cor_prime_ratio     << endl;

//	exit(0);


	TH1F *frame = beta_coin_can->DrawFrame(-15, -0.3, 15, 1.7);
	
 	obj_tree->Draw("hsbeta:ct_corrected", "(" + Cointime_Cut + " || "+ Rand_Cut + ")" + " && hsbeta < 1.5 && hsbeta > -1.5", "psame");

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


	beta_coin_can->Clear();


	obj_tree->Draw("ct_corrected >>htemp(300,-15,15)", "hsbeta < 0.2 && hsbeta > -0.1", "");

	beta_coin_can->Write("zero_check");
	

	delete beta_coin_can;

//	delete setting_tree;

}










/*--------------------------------------------------*/

Analysis_heep::~Analysis_heep() {

	file_out_ana->Close();

}





