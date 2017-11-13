TH1F* peak;
TH1F* peak1;
TH1F* peak2;

TH1F* peak_omega;
TH1F* peak_rho;
TH1F* peak_xphsp;
TH1F* peak_eta;
TH1F* peak_etap;

TF1* fun_t1 = new TF1();
TF1* fun_t2 = new TF1();
TF1* fun_t3 = new TF1();

// const Int_t u_bin_num = 3;
// const Int_t phi_bin_num = 8;

// 


// const Int_t u_bin_num = 2;
// const Int_t phi_bin_num = 10;
// const Int_t bin_num = u_bin_num * phi_bin_num;
// 
// 




Int_t u_bin_num;
Int_t phi_bin_num;


// 


//TH1F * peak[6];

// Double_t acc_yield         [bin_num];
// Double_t acc_sim_yield     [bin_num];
// 
// Double_t acc_yield_err     [bin_num];
// Double_t acc_sim_yield_err [bin_num];
 
Double_t* acc_yield;
Double_t* acc_sim_yield;

Double_t* acc_yield_err;
Double_t* acc_sim_yield_err;


ofstream yield_file_out;

bool is_itt;
bool bg_refit;


TString out_dir_str;

out_dir_str = "u_bin_fit/";


TString ave_dir_str;
ave_dir_str = "averages/";

TFile* file_pre;

void bg_background_test () {


//	gROOT->LoadMacro("u_phi_fit_functions.C");

	
	file_pre = new TFile(out_dir_str + "presnetation_plots.root","recreate"); 

	gROOT->ProcessLine(".L u_phi_fit_functions.C");
	gROOT->ProcessLine(".L fitting_functions.C");
	gROOT->ProcessLine(".L recon_kin.C");
	gROOT->ProcessLine(".L get_u_phi_bin.C");

	TString file_name;

	TString q2;
	TString eps;
	TString theta;

//	file_name = "final_160_32_3000_out_put.root";

	is_itt = false;

	bool bg_refit;


	u_bin_num   = Get_u_bin();
	phi_bin_num = Get_phi_bin();
	
	const Int_t u_bin_num_const   = u_bin_num;  
	const Int_t phi_bin_num_const = phi_bin_num;  
	const Int_t bin_num_const     = u_bin_num_const * phi_bin_num_const;

 	acc_yield         =  new Double_t[bin_num_const];
 	acc_sim_yield     =  new Double_t[bin_num_const];
                                           
 	acc_yield_err     = new Double_t[bin_num_const];
 	acc_sim_yield_err = new Double_t[bin_num_const];

	for(Int_t i =0; i < bin_num_const; i++) {

		acc_yield[i]         = 0.0;
        acc_sim_yield[i]     = 0.0;
                            
	    acc_yield_err[i]     = 0.0;
        acc_sim_yield_err[i] = 0.0;

	}


	cout << Get_u_bin() << "    " << Get_phi_bin() << endl;

	TString q2_setting[2]     = {"160", "245"};
	
	TString eps_160[2]        = {"32", "59"};
	TString eps_245[2]        = {"27", "55"};
	
	TString hms_angle_160_l[2]  = {"+0970", "+3000"};
	TString hms_angle_160_h[3]  = {"-2730",  "+0000", "+3000"};
	
	TString hms_angle_245_l[2]  = {"+1350", "+3000"};
	TString hms_angle_245_h[3]  = {"-3000", "+0000", "+3000"};

	for (int i = 0; i < 2; i++) {

		q2 = q2_setting[i];

		if( q2 =="160") {

			/*--------------------------------------------------*/
			/// 2nd order polynomial

			for (int dd = 0; dd < 2; dd++) {

				eps = eps_160[dd];

				if ( eps.Atoi()  < 40) {
					
					yield_file_out.open(ave_dir_str + "aver.pl_" + q2 + "_" + eps + ".dat", std::fstream::out); 

//					yield_file_out << 11111 << endl;

					Init_yield_array();

					for (int iii = 0; iii < 2; iii++) {

//					for (int iii = 0; iii < 1; iii++) {

						theta = hms_angle_160_l[iii]; 

						file_name = q2 + "_" + eps + "_" + theta;

						if (iii == 0) {
//							fit_bg(file_name, 0.6, 0.85, 0.70, 0.9, 0.76, 0.88);
							fit_bg(file_name, 0.63, 0.80, 0.70, 0.83, 0.70, 0.90);
						} else {
//							fit_bg(file_name, 0.6, 0.85, 0.70, 0.9, 0.70, 1.0);
							fit_bg(file_name, 0.63, 0.80, 0.70, 0.83, 0.70, 0.90);
						}

//						fit_bg(file_name, 0.6, 0.80, 0.70, 0.88, 0.75, 0.97);

					}

					Output_yield_array();

					

				} else {

					Init_yield_array();

					yield_file_out.open(ave_dir_str + "aver.pl_"+q2 + "_" + eps + ".dat", std::fstream::out); 
 
					for (int iii = 0; iii < 3; iii++) {
				
						theta = hms_angle_160_h[iii]; 

						file_name = q2 + "_" + eps + "_" + theta;

						if (iii == 0) {
//							fit_bg(file_name, 0.60, 0.85, 0.70, 0.9, 0.70, 1);
							fit_bg(file_name, 0.63, 0.80, 0.70, 0.83, 0.77, 0.95);
						} else if (iii == 1) {
//							fit_bg(file_name, 0.66, 0.85, 0.70, 0.9, 0.76, 0.92);
							fit_bg(file_name, 0.63, 0.80, 0.70, 0.83, 0.77, 0.95);
						} else {
//							fit_bg(file_name, 0.66, 0.85, 0.70, 0.9, 0.76, 0.94);
							fit_bg(file_name, 0.63, 0.80, 0.70, 0.83, 0.77, 0.95);
						}

//						fit_bg(file_name, 0.6, 0.80, 0.74, 0.88, 0.77, 0.95);

					}

					Output_yield_array();
 
 				}

			}

/*--------------------------------------------------*/

		}	else {

			for (int ii = 0; ii < 2; ii++) {

				eps = eps_245[ii];

				if (eps.Atoi() < 40) {

					yield_file_out.open(ave_dir_str + "aver.pl_"+q2 + "_" + eps + ".dat", std::fstream::out); 

					Init_yield_array();

					for (int iii = 0; iii < 2; iii++) {
//					for (int iii = 0; iii < 1; iii++) {

						theta = hms_angle_245_l[iii]; 

						file_name = q2 + "_" + eps + "_" + theta;

						if (iii == 0) {
//							fit_bg(file_name, 0.60, 0.85, 0.73, 0.88, 0.75, 0.94);
							fit_bg(file_name, 0.63, 0.80, 0.70, 0.83, 0.77, 0.95);
						} else {
//							fit_bg(file_name, 0.60, 0.85, 0.71, 0.85, 0.75, 0.94);
							fit_bg(file_name, 0.63, 0.80, 0.70, 0.83, 0.77, 0.95);
						}

//						fit_bg(file_name, 0.6, 0.80, 0.75, 0.88, 0.75, 0.97);

					}


					Output_yield_array();
 
				} else {

					yield_file_out.open(ave_dir_str + "aver.pl_"+q2 + "_" + eps + ".dat", std::fstream::out); 

					Init_yield_array();

					for (int iii = 0; iii < 3; iii++) {

						theta = hms_angle_245_h[iii]; 

						file_name = q2 + "_" + eps + "_" + theta;

						if (iii == 0) {
//							fit_bg(file_name, 0.65, 0.85, 0.68, 0.88, 0.7, 1.00);
							fit_bg(file_name, 0.63, 0.80, 0.70, 0.83, 0.77, 0.95);
						} else if (iii == 1) 
//							fit_bg(file_name, 0.60, 0.85, 0.73, 0.87, 0.70, 0.98);
							fit_bg(file_name, 0.63, 0.80, 0.70, 0.83, 0.77, 0.95);
						else {
//							fit_bg(file_name, 0.65, 0.85, 0.71, 0.85, 0.75, 1.00);
							fit_bg(file_name, 0.63, 0.80, 0.70, 0.83, 0.77, 0.95);
						}
// 
//						fit_bg(file_name, 0.6, 0.80, 0.80, 0.92, 0.75, 0.97);

					}

					

					Output_yield_array();

				}			

			}
			
  		}


// 			for (int ii = 0; ii < 2; ii++) {
// 
// 				eps = eps_245[ii];
// 
// 				if (eps.Atoi() < 40) {
// 
// 					yield_file_out.open(ave_dir_str + "aver.pl_"+q2 + "_" + eps + ".dat", std::fstream::out); 
// 
// 					Init_yield_array();
// 
// 					for (int iii = 0; iii < 2; iii++) {
// //					for (int iii = 0; iii < 1; iii++) {
// 
// 						theta = hms_angle_245_l[iii]; 
// 
// 						file_name = q2 + "_" + eps + "_" + theta;
// 
// 						if (iii == 0) {
// //							fit_bg(file_name, 0.60, 0.85, 0.73, 0.88, 0.75, 0.94);
// 							fit_bg(file_name, 0.65, 0.79, 0.75, 0.85, 0.78, 0.90);
// 						} else {
// //							fit_bg(file_name, 0.60, 0.85, 0.71, 0.85, 0.75, 0.94);
// 							fit_bg(file_name, 0.65, 0.80, 0.75, 0.88, 0.75, 0.97);
// 						}
// 
// //						fit_bg(file_name, 0.6, 0.80, 0.75, 0.88, 0.75, 0.97);
// 
// 					}
// 
// 
// 					Output_yield_array();
//  
// 				} else {
// 
// 					yield_file_out.open(ave_dir_str + "aver.pl_"+q2 + "_" + eps + ".dat", std::fstream::out); 
// 
// 					Init_yield_array();
// 
// 					for (int iii = 0; iii < 3; iii++) {
// 
// 						theta = hms_angle_245_h[iii]; 
// 
// 						file_name = q2 + "_" + eps + "_" + theta;
// 
// 						if (iii == 0) {
// //							fit_bg(file_name, 0.65, 0.85, 0.68, 0.88, 0.7, 1.00);
// 							fit_bg(file_name, 0.6, 0.85, 0.70, 0.90, 0.75, 0.97);
// 						} else if (iii == 1) 
// //							fit_bg(file_name, 0.60, 0.85, 0.73, 0.87, 0.70, 0.98);
// 							fit_bg(file_name, 0.6, 0.80, 0.78, 0.88, 0.55, 0.94);
// 						else {
// //							fit_bg(file_name, 0.65, 0.85, 0.71, 0.85, 0.75, 1.00);
// 							fit_bg(file_name, 0.6, 0.80, 0.75, 0.85, 0.75, 0.97);
// 						}
// // 
// //						fit_bg(file_name, 0.6, 0.80, 0.80, 0.92, 0.75, 0.97);
// 
// 					}
// 
// 					
// 
// 					Output_yield_array();
// 
// 				}			
//
//			}


	}	




// 
// 	if float(setting) == 160:
// 		for eps in eps_160:
// 			if float(eps) == 32:
// 				for angle in hms_angle_160_l:
// 					print setting + "_"+ eps  + "_" + angle
// 					Setting_by_setting(setting + "_"+ eps  + "_" + angle)
// 
// 			else:
// 				for angle in hms_angle_160_h:
// 					print setting + "_"+ eps  + "_" + angle
// 					Setting_by_setting(setting + "_"+ eps  + "_" + angle)
// 		
// 
// 
// 	else:		
// 		for eps in eps_245:
// 			if float(eps) == 27:
// 				for angle in hms_angle_245_l:
// 					print setting + "_"+ eps  + "_" + angle
// 					Setting_by_setting(setting + "_"+ eps  + "_" + angle)
// 			else:
// 				for angle in hms_angle_245_h:
// 					print setting + "_"+ eps  + "_" + angle
// 					Setting_by_setting(setting + "_"+ eps  + "_" + angle)
// 
// 
	




	
	
	// file_name = "final_160_59_0000_out_put.root";

	// fit_bg(file_name);	

	file_pre->Close();

}


/*--------------------------------------------------*/
/*--------------------------------------------------*/
// Version: 1 


void fit_bg (TString file_name, Double_t fit_low_1, Double_t fit_low_2, Double_t fit_mid_1, Double_t fit_mid_2, Double_t fit_high_1, Double_t fit_high_2) {

	TFile *file1 = TFile::Open("final_"+ file_name + "_out_put.root");
	TCanvas* c1 = (TCanvas*)file1->Get("missmass_dummy_sub_sim");

	//	c1->Draw();

   	TH1F* miss_target= (TH1F*)c1->GetPrimitive("missmass_real_dia");
   	peak1= (TH1F*)c1->GetPrimitive("missmass_dia");

	TCanvas* c2 = new TCanvas();

	miss_target->Draw();

	TF1 *ftot = new TF1("ftot",ftotal_all, fit_low_1, fit_high_2 ,5);

//	ftot->SetParameters(3.5, -16.4, 25, -12.2, 20);

	Double_t norm = miss_target->GetMaximum();

	//	ftot->SetParameters(0, 0, 0, 0);
	//	ftot->SetParLimits(0,.1*norm,norm);

	miss_target->Fit("ftot","MUER");   

	//	cout << file_name + ".pdf"<< endl; 

   	c2->Print(out_dir_str + file_name + ".png");
	



// 	/*--------------------------------------------------*/
// 	/// Average kinematics 
// 
// 	Double_t w_ave, w_ave_err;
// 
// 	TCanvas* can_w;
// 	can_w = (TCanvas*)file1->Get("W_u_bin");
//
//	average_kin(can_w, &w_ave, &w_ave_err);

 

//	cout << w_ave << "    " << w_ave_err << endl;

//	exit(0);




	/*--------------------------------------------------*/
	/// t binning 

	TCanvas* c4;
	c4 = (TCanvas*)file1->Get("missmass_u_bin");
//	t_fit_2(c4, file_name, fit_low_1, fit_low_2, fit_mid_1, fit_mid_2, fit_high_1, fit_high_2);
//	u_fit(c4, file_name, fit_low_1, fit_low_2, fit_mid_1, fit_mid_2, fit_high_1, fit_high_2);
//	u_fit_full(c4, file_name, fit_low_1, fit_low_2, fit_mid_1, fit_mid_2, fit_high_1, fit_high_2);
// 
// 
//  
// 	/*--------------------------------------------------*/
// 	/// Individual u_phi binning
// 	// Calculate t-phi bin ratio and average ratio
// 





	TCanvas* c5;
	c5 = (TCanvas*)file1->Get("missmass_u_phi_bin");




//    bg_refit = true;
	bg_refit = false;

	 


	if(bg_refit) {

		average_kin(file1, file_name);
		u_phi_fit_sim(c5, file_name);
	
	} else {

	 	u_phi_int_sim_itt(c5, file_name);
	 	
	 	TCanvas* c6;
	 	c6 = (TCanvas*)file1->Get("Em_u_phi_bin");
	 	recon_kin(c6, "Em", file_name);
	 	
	 	TCanvas* c6;
	 	c6 = (TCanvas*)file1->Get("Pm_u_phi_bin");
	 	recon_kin(c6, "Pm", file_name);
	 	
	 	TCanvas* c6;
	 	c6 = (TCanvas*)file1->Get("hsdelta_u_phi_bin");
	 	recon_kin(c6, "hsbelta", file_name);
	 	
	 	TCanvas* c6;
	 	c6 = (TCanvas*)file1->Get("hsxptar_u_phi_bin");
	 	recon_kin(c6, "hsxptar", file_name);
	 	    
	 	TCanvas* c6;
	 	c6 = (TCanvas*)file1->Get("hsyptar_u_phi_bin");
	 	recon_kin(c6, "hsyptar", file_name);

// 	 	TCanvas* c6;
// 	 	c6 = (TCanvas*)file1->Get("t_u_phi_bin");
// 	 	recon_kin(c6, "t", file_name);
// 
// 	 	TCanvas* c6;
// 	 	c6 = (TCanvas*)file1->Get("u_u_phi_bin");
// 	 	recon_kin(c6, "u", file_name);
// 
// 	 	TCanvas* c6;
// 	 	c6 = (TCanvas*)file1->Get("W_u_phi_bin");
// 	 	recon_kin(c6, "W", file_name);
// 
// 	 	TCanvas* c6;
// 	 	c6 = (TCanvas*)file1->Get("Q2_u_phi_bin");
// 	 	recon_kin(c6, "Q2", file_name);
// 



	}




	/*--------------------------------------------------*/
	/// u_phi bin fitting without common background


//	u_phi_fit_sim(c3, file_name, fit_low_1, fit_low_2, fit_mid_1, fit_mid_2, fit_high_1, fit_high_2);
//	exit(0);

  

   	delete c1;
   	delete c2;


   	file1->Close();

}






/*--------------------------------------------------*/
/*--------------------------------------------------*/
// Version: 2

void fit_bg (TString file_name, Double_t fit_low, Double_t fit_high) {

	TFile *file1 = TFile::Open("final_"+ file_name + "_out_put.root");
	TCanvas* c1 = (TCanvas*)file1->Get("missmass_dummy_sub_sim");

	//	c1->Draw();

   	TH1F* miss_target= (TH1F*)c1->GetPrimitive("missmass_real");
   	peak1= (TH1F*)c1->GetPrimitive("missmass");

	TCanvas* c2 = new TCanvas();

	miss_target->Draw();

	TF1 *ftot = new TF1("ftot",ftotal_all, fit_low, fit_high ,5);

//	ftot->SetParameters(3.5, -16.4, 25, -12.2, 20);

	Double_t norm = miss_target->GetMaximum();

	//	ftot->SetParameters(0, 0, 0, 0);
	//	ftot->SetParLimits(0,.1*norm,norm);

	miss_target->Fit("ftot","MUER");   

	//	cout << file_name + ".pdf"<< endl; 

   	c2->Print(out_dir_str + file_name + ".png");


	
	TCanvas* c3 = new TCanvas();
	c3 = (TCanvas*)file1->Get("testtest4");
	u_phi_fit(c3);
// 
// //
// //
// //
// 
// 
// 	TCanvas* c5 = new TCanvas();
// 	c5 = (TCanvas*)file1->Get("testtest5");
// 	u_fit(c5, file_name, fit_low, fit_high);
// 

   	delete c1;
   	delete c2;

   	file1->Close();

}







// /*--------------------------------------------------*/
// /*--------------------------------------------------*/
// /*--------------------------------------------------*/
// // Just fit missing mass in t
// // Version: 1
// 
// 
// void u_fit(TCanvas* cc3, TString file_name, Double_t fit_low, Double_t fit_mid, Double_t fit_high, Double_t inter_1, Double_t inter_2, Double_t inter_3) {
// 	
// 	cout << "asddddddddddddddddd " << endl;
// 
// // 	cc3->Draw();
// // 
//  	TCanvas* cc5 = new TCanvas("cc5", "cc5", 1600, 800);
//  	cc5->Divide(3,1, 0.003);
// // 
//  	TH1F* h_tgt[3];
//  	TH1F* h_sim[3];
// 
// //	TF1* fit_u_fun[3];
// 
// // 
//  	for (Int_t i = 1; i <=3; i++) {
// // 	
// 		TPad* cur_pad = cc5->GetPad(i);
// 
// 		cur_pad->cd();
// 		cur_pad->SetLeftMargin(0.1);
// 		cur_pad->SetRightMargin(0.002);
// 
// 
// 
// 		cout << i << "   " << i-1 << endl;
// 
// 		const int in =  i -1;
// 
// 		cout << in << endl;
// 
// 		TH1F* h_tgt_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("phi_real_var");
// 		TH1F* h_tgt[in] = (TH1F*) h_tgt_tmp->Clone("h_tgt_clone"); 
// 
// 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("missmass_u_sim");
// 		TH1F* h_sim[in] = (TH1F*) h_sim_tmp->Clone("h_sim_clone");
// 
// 
// 		
// 	 	h_tgt[in]->DrawCopy();
// 
// //	 	h_tgt[in]->DrawClone();
// 
// //		h_sim[in]->Scale(8);
//   
// 		h_sim[in]->SetLineColor(2);
// 
// 	 	h_sim[in]->DrawCopy("same");
// 
// 		peak1 = h_sim[in];
// 		
// // 		TF1 *fit_u_fun;
// // 
// // 		if (i==1) {
// // 			 fit_u_fun = new TF1("fit_u_fun",ftotal_all, fit_low, fit_high ,5);
// // 		} else if (i==2) {
// // 			 fit_u_fun = new TF1("fit_u_fun",ftotal_all, fit_low, fit_high ,5);
// // 		} else {
// // 			 fit_u_fun = new TF1("fit_u_fun",ftotal_all, fit_low, fit_high ,5);
// // 		}
// 
// 
// // 		fit_local_lo = fit_low;
// // 		fit_local_hi = fit_high;
// // 
// 
// 
// 		Double_t fit_local_lo;
// 		Double_t fit_local_hi;
// 
// 		if (i==1) {
// 
// 			fit_local_lo = fit_low;
// 			fit_local_hi = fit_low+inter_1;
// 
// 
// 		} else if (i==2) {
// 
// 			fit_local_lo = fit_mid;
// 			fit_local_hi = fit_mid +inter_2;
// 
// 		} else {
// 
// 			fit_local_lo  = fit_high;
// 			fit_local_hi = fit_high + inter_3;
// 
// 		}
// // 
// 
// 		
// 
// 
// // 		TF1 *fit_u_fun = new TF1("fit_u_fun",ftotal_all, fit_low, fit_high ,5);
// // 
// // 		h_tgt[in]->Fit("fit_u_fun", "MUER");
// 		
// 
// //			TF1 *fit_u_fun = new TF1("fit_u_fun", ftotal, fit_low, fit_high ,5);
// 
// 	
// //		TF1 *fit_u_fun = new TF1("fit_u_fun", ftotal_all, fit_local_lo, fit_local_hi, 5);
// 		TF1 *fit_u_fun = new TF1("fit_u_fun", ftotal_all_2, fit_local_lo, fit_local_hi, 4);
// 
// //		fit_u_fun->SetNpx(400);
// 
// 		fit_u_fun->SetParameters(10.0, -50.0, 75, -36.0, 15);
// 
// //  			fit_u_fun->FixParameter(0, 11.7661);
// //  			fit_u_fun->FixParameter(1, -54.3);
// //  			fit_u_fun->FixParameter(2, 82.7825);
// //  			fit_u_fun->FixParameter(3, -41.5146);
// // 
// //  			fit_u_fun->FixParameter(4, 200);
// 			
// 
// // 			fit_u_fun->SetParLimits(0, 0, 50);
// // 			fit_u_fun->SetParLimits(1, -100, 0);
// // 			fit_u_fun->SetParLimits(2, 0, 200);
// // 			fit_u_fun->SetParLimits(3, -100, 0);
//  			fit_u_fun->SetParLimits(4, 0, 60);
// // 
// // 
// 
// 		h_tgt[in]->Fit("fit_u_fun", "MEURN");
// 
// 
// //			h_tgt[in]->Fit("fit_u_fun", "VR+");
// 
// //			h_sim[in]->Scale(18.868);
// 	 	h_sim[in]->Draw("same");
// 
// //			fit_u_fun->Draw("");
// 
// 	 	TH1F *h_sim_clone = (TH1F*) h_sim[in]->Clone();
// 		
// //		TF1 *fit_u_fun_pol = new TF1("fit_u_fun_pol", ftotal_br, fit_local_lo, fit_local_hi, 4);
// 		TF1 *fit_u_fun_pol = new TF1("fit_u_fun_pol", ftotal_br_2, fit_local_lo, fit_local_hi, 4);
// 
// 		fit_u_fun_pol->SetLineColor(4);
// 
//  		fit_u_fun_pol->FixParameter(0, fit_u_fun->GetParameter(0));
//  		fit_u_fun_pol->FixParameter(1, fit_u_fun->GetParameter(1));
//  		fit_u_fun_pol->FixParameter(2, fit_u_fun->GetParameter(2));
// // 		fit_u_fun_pol->FixParameter(3, fit_u_fun->GetParameter(3));
// 
// 		fit_u_fun_pol->Draw("same");
// 
// 		h_sim_clone->SetLineColor(6);
// //		h_sim_clone->Scale(fit_u_fun->GetParameter(4));		
// 		h_sim_clone->Scale(fit_u_fun->GetParameter(3));		
//  	 	h_sim_clone->DrawCopy("same");
// 
//  		h_sim_clone->Add(fit_u_fun_pol);		
// 		h_sim_clone->SetLineColor(1);
//  	 	h_sim_clone->DrawCopy("same");
// 
// 
// 
// 
//  	}
// 
// 
// 
// 	cc5->Print(out_dir_str + "t_fit_" + file_name+".png");
// 
// }
// 






/*--------------------------------------------------*/
/*--------------------------------------------------*/
/*--------------------------------------------------*/
// Just fit missing mass in t
// Version: 1


void u_fit(TCanvas* cc3, TString file_name, Double_t fit_low_1, Double_t fit_low_2, Double_t fit_mid_1, Double_t fit_mid_2, Double_t fit_high_1, Double_t fit_high_2) {
	
	cout << "asddddddddddddddddd " << endl;

// 	cc3->Draw();
// 
 	TCanvas* cc5 = new TCanvas("cc5", "cc5", 1600, 800);
 	cc5->Divide(u_bin_num,1, 0.003);
// 
 	TH1F* h_tgt[3];
 	TH1F* h_sim[3];
 	TH1F* h_sim_rho[3];
 	TH1F* h_sim_xphsp[3];

//	TF1* fit_u_fun[3];

// 
 	for (Int_t i = 1; i <=3; i++) {
// 	
		TPad* cur_pad = cc5->GetPad(i);

		cur_pad->cd();
		cur_pad->SetLeftMargin(0.1);
		cur_pad->SetRightMargin(0.002);



		cout << i << "   " << i-1 << endl;

		const int in =  i -1;

		cout << in << endl;

		TH1F* h_tgt_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("u_target");
		TH1F* h_tgt[in] = (TH1F*) h_tgt_tmp->Clone("h_tgt_clone"); 

		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("u_sim_omega");
		TH1F* h_sim[in] = (TH1F*) h_sim_tmp->Clone("h_sim_omega_clone");

		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("u_sim_rho");
		TH1F* h_sim_rho[in] = (TH1F*) h_sim_tmp->Clone("h_sim_rho_clone");

		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("u_sim_xphsp");
		TH1F* h_sim_xphsp[in] = (TH1F*) h_sim_tmp->Clone("h_sim_xphsp_clone");




		
	 	h_tgt[in]->DrawCopy();

//	 	h_tgt[in]->DrawClone();

//		h_sim[in]->Scale(8);
  
		h_sim[in]->SetLineColor(2);

	 	h_sim[in]->DrawCopy("same");

		peak_omega = h_sim[in];
		peak_rho   = h_sim_rho[in];
		peak_xphsp = h_sim_xphsp[in];
		
// 		TF1 *fit_u_fun;
// 		if (i==1) {
// 			 fit_u_fun = new TF1("fit_u_fun",ftotal_all, fit_low, fit_high ,5);
// 		} else if (i==2) {
// 			 fit_u_fun = new TF1("fit_u_fun",ftotal_all, fit_low, fit_high ,5);
// 		} else {
// 			 fit_u_fun = new TF1("fit_u_fun",ftotal_all, fit_low, fit_high ,5);
// 		}

// 		fit_local_lo = fit_low;
// 		fit_local_hi = fit_high;

		Double_t fit_local_lo;
		Double_t fit_local_hi;

		if (i==1) {

			fit_local_lo = fit_low_1;
			fit_local_hi = fit_low_2;


		} else if (i==2) {

			fit_local_lo = fit_mid_1;
			fit_local_hi = fit_mid_2;

		} else {

			fit_local_lo = fit_high_1;
			fit_local_hi = fit_high_2;

		}

// 		fit_local_lo = 0.7;
// 		fit_local_lo = 0.85;
	
		TF1 *fit_u_fun = new TF1("fit_u_fun", fun_sim_total, fit_local_lo, fit_local_hi, 3);



 		fit_u_fun->SetParameters(30, 400, 5000);

		fit_u_fun->SetParLimits(0, 0, 100000);
		fit_u_fun->SetParLimits(1, 0, 100000);
		fit_u_fun->SetParLimits(2, 0, 100000);


		h_tgt[in]->Fit("fit_u_fun", "MNR");

//	 	h_sim[in]->Draw("same");

	 	TH1F *h_sim_omega_clone = (TH1F*) h_sim[in]->Clone();
	 	TH1F *h_sim_rho_clone   = (TH1F*) h_sim_rho[in]->Clone();
	 	TH1F *h_sim_xphsp_clone = (TH1F*) h_sim_xphsp[in]->Clone();


//		h_sim_clone->SetLineColor(6);

		h_sim_omega_clone->Scale(fit_u_fun->GetParameter(0));		
		h_sim_rho_clone->Scale(fit_u_fun->GetParameter(1));		
		h_sim_xphsp_clone->Scale(fit_u_fun->GetParameter(2));		

 	 	h_sim_omega_clone->DrawCopy("same");
 	 	h_sim_rho_clone->DrawCopy("same");
 	 	h_sim_xphsp_clone->DrawCopy("same");


 		h_sim_xphsp_clone->Add(h_sim_omega_clone);		
 		h_sim_xphsp_clone->Add(h_sim_rho_clone);		


// 		h_sim_clone->Add(fit_u_fun_pol);		
		h_sim_xphsp_clone->SetLineColor(6);
 	 	h_sim_xphsp_clone->DrawCopy("same");




// 		TF1 *fit_u_fun_eta = new TF1("fit_u_fun_pol", ftotal_eta_prime_peak, 0.9, 1, 3);
// 
//  		fit_u_fun_eta->FixParameter(0, fit_u_fun->GetParameter(5));
//  		fit_u_fun_eta->FixParameter(1, fit_u_fun->GetParameter(6));
//  		fit_u_fun_eta->FixParameter(2, fit_u_fun->GetParameter(7));
// 
// 		fit_u_fun_eta->SetLineColor(2);
//  	 	fit_u_fun_eta->DrawCopy("same");
// 
	

		


		
		if (i==1) {

//			fun_t1 = (TF1*) fit_u_fun->Clone("fun_t1");
//			fun_t1 = (TF1*) fit_u_fun->Clone("fun_t1");

			fit_u_fun->Copy(*fun_t1);


		} else if (i==2) {

//			fun_t2 = (TF1*) fit_u_fun->Clone("fun_t2");

			fit_u_fun->Copy(*fun_t2);

		} else {

//			fun_t3 = (TF1*) fit_u_fun->Clone("fun_t3");
			fit_u_fun->Copy(*fun_t3);

		}
// 



 	}

	cc5->Print(out_dir_str + "u_fit_" + file_name+".png");

	delete cc5;

}
// 



/*--------------------------------------------------*/
/*--------------------------------------------------*/
/*--------------------------------------------------*/
// Just fit missing mass in t
// Version: 3


void u_fit_full(TCanvas* cc3, TString file_name, Double_t fit_low_1, Double_t fit_low_2, Double_t fit_mid_1, Double_t fit_mid_2, Double_t fit_high_1, Double_t fit_high_2) {
	
	cout << "asddddddddddddddddd " << endl;

// 	cc3->Draw();
// 
 	TCanvas* cc5 = new TCanvas("cc5", "cc5", 1600, 800);
 	cc5->Divide(u_bin_num,1, 0.003);
// 
 	TH1F* h_tgt[3];
 	TH1F* h_sim[3];
 	TH1F* h_sim_rho[3];
 	TH1F* h_sim_xphsp[3];
 	TH1F* h_sim_eta[3];
 	TH1F* h_sim_etap[3];

//	TF1* fit_u_fun[3];

// 
 	for (Int_t i = 1; i <=3; i++) {
// 	
		TPad* cur_pad = cc5->GetPad(i);

		cur_pad->cd();
		cur_pad->SetLeftMargin(0.1);
		cur_pad->SetRightMargin(0.002);



		cout << i << "   " << i-1 << endl;

		const int in =  i -1;

		cout << in << endl;

		TH1F* h_tgt_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("u_target");
		TH1F* h_tgt[in] = (TH1F*) h_tgt_tmp->Clone("h_tgt_clone"); 

		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("u_sim_omega");
		TH1F* h_sim[in] = (TH1F*) h_sim_tmp->Clone("h_sim_omega_clone");

		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("u_sim_rho");
		TH1F* h_sim_rho[in] = (TH1F*) h_sim_tmp->Clone("h_sim_rho_clone");

		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("u_sim_xphsp");
		TH1F* h_sim_xphsp[in] = (TH1F*) h_sim_tmp->Clone("h_sim_xphsp_clone");

		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("u_sim_eta");
		TH1F* h_sim_eta[in] = (TH1F*) h_sim_tmp->Clone("h_sim_eta_clone");

		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("u_sim_etap");
		TH1F* h_sim_etap[in] = (TH1F*) h_sim_tmp->Clone("h_sim_etap_clone");




		
	 	h_tgt[in]->DrawCopy();

//	 	h_tgt[in]->DrawClone();

//		h_sim[in]->Scale(8);
  
 		h_sim[in]->SetLineColor(2);
// 
// 	 	h_sim[in]->DrawCopy("same");

		peak_omega = h_sim[in];
		peak_rho   = h_sim_rho[in];
		peak_xphsp = h_sim_xphsp[in];
		peak_eta   = h_sim_eta[in];
		peak_etap  = h_sim_etap[in];
		
// 		TF1 *fit_u_fun;
// 		if (i==1) {
// 			 fit_u_fun = new TF1("fit_u_fun",ftotal_all, fit_low, fit_high ,5);
// 		} else if (i==2) {
// 			 fit_u_fun = new TF1("fit_u_fun",ftotal_all, fit_low, fit_high ,5);
// 		} else {
// 			 fit_u_fun = new TF1("fit_u_fun",ftotal_all, fit_low, fit_high ,5);
// 		}

// 		fit_local_lo = fit_low;
// 		fit_local_hi = fit_high;

		Double_t fit_local_lo;
		Double_t fit_local_hi;

		if (i==1) {

			fit_local_lo = fit_low_1;
			fit_local_hi = fit_low_2;


		} else if (i==2) {

			fit_local_lo = fit_mid_1;
			fit_local_hi = fit_mid_2;

		} else {

			fit_local_lo = fit_high_1;
			fit_local_hi = fit_high_2;

		}

// 		fit_local_lo = 0.7;
// 		fit_local_lo = 0.85;
	
//		TF1 *fit_u_fun = new TF1("fit_u_fun", fun_sim_total, fit_local_lo, fit_local_hi, 3);
		TF1 *fit_u_fun = new TF1("fit_u_fun", fun_sim_total, fit_local_lo, fit_local_hi, 5);



 		fit_u_fun->SetParameters(30, 400, 5000);

		fit_u_fun->SetParLimits(0, 0, 100000);
		fit_u_fun->SetParLimits(1, 0, 100000);
		fit_u_fun->SetParLimits(2, 0, 100000);
		fit_u_fun->SetParLimits(3, 0, 100000);
		fit_u_fun->SetParLimits(4, 0, 100000);


		h_tgt[in]->Fit("fit_u_fun", "MNR");

//	 	h_sim[in]->Draw("same");

	 	TH1F *h_sim_omega_clone = (TH1F*) h_sim[in]->Clone();
	 	TH1F *h_sim_rho_clone   = (TH1F*) h_sim_rho[in]->Clone();
	 	TH1F *h_sim_xphsp_clone = (TH1F*) h_sim_xphsp[in]->Clone();
	 	TH1F *h_sim_eta_clone   = (TH1F*) h_sim_eta[in]->Clone();
	 	TH1F *h_sim_etap_clone  = (TH1F*) h_sim_etap[in]->Clone();




//		h_sim_clone->SetLineColor(6);

		h_sim_omega_clone->Scale(fit_u_fun->GetParameter(0));		
		h_sim_rho_clone->Scale(fit_u_fun->GetParameter(1));		
		h_sim_xphsp_clone->Scale(fit_u_fun->GetParameter(2));		
		h_sim_eta_clone->Scale(fit_u_fun->GetParameter(3));		
		h_sim_etap_clone->Scale(fit_u_fun->GetParameter(4));		

 	 	h_sim_omega_clone->DrawCopy("same");
 	 	h_sim_rho_clone->DrawCopy("same");
 	 	h_sim_xphsp_clone->DrawCopy("same");
		h_sim_eta_clone->DrawCopy("same");		
		h_sim_etap_clone->DrawCopy("same");		


 		h_sim_xphsp_clone->Add(h_sim_omega_clone);		
 		h_sim_xphsp_clone->Add(h_sim_rho_clone);		
 		h_sim_xphsp_clone->Add(h_sim_eta_clone);		
 		h_sim_xphsp_clone->Add(h_sim_etap_clone);		


// 		h_sim_clone->Add(fit_u_fun_pol);		
		h_sim_xphsp_clone->SetLineColor(6);
 	 	h_sim_xphsp_clone->DrawCopy("same");




// 		TF1 *fit_u_fun_eta = new TF1("fit_u_fun_pol", ftotal_eta_prime_peak, 0.9, 1, 3);
// 
//  		fit_u_fun_eta->FixParameter(0, fit_u_fun->GetParameter(5));
//  		fit_u_fun_eta->FixParameter(1, fit_u_fun->GetParameter(6));
//  		fit_u_fun_eta->FixParameter(2, fit_u_fun->GetParameter(7));
// 
// 		fit_u_fun_eta->SetLineColor(2);
//  	 	fit_u_fun_eta->DrawCopy("same");
// 
	

		


		
		if (i==1) {

//			fun_t1 = (TF1*) fit_u_fun->Clone("fun_t1");
//			fun_t1 = (TF1*) fit_u_fun->Clone("fun_t1");

			fit_u_fun->Copy(*fun_t1);


		} else if (i==2) {

//			fun_t2 = (TF1*) fit_u_fun->Clone("fun_t2");

			fit_u_fun->Copy(*fun_t2);

		} else {

//			fun_t3 = (TF1*) fit_u_fun->Clone("fun_t3");
			fit_u_fun->Copy(*fun_t3);

		}
// 



 	}

	file_pre->cd();

	cc5->Print(out_dir_str + "u_fit_" + file_name+".png");
	cc5->Write("u_fit_" + file_name);

	delete cc5;

}
// 




/*--------------------------------------------------*/
/*--------------------------------------------------*/
// Just fit missing mass in t
// Version: 2

void u_fit(TCanvas* cc3, TString file_name, Double_t fit_low, Double_t fit_high) {
	
	cout << "asddddddddddddddddd " << endl;

// 	cc3->Draw();
// 
 	TCanvas* cc5 = new TCanvas("cc5", "cc5", 1600, 800);
 	cc5->Divide(u_bin_num,1, 0.003);
// 
 	TH1F* h_tgt[3];
 	TH1F* h_sim[3];

//	TF1* fit_u_fun[3];

// 
 	for (Int_t i = 1; i <=3; i++) {
// 	
		TPad* cur_pad = cc5->GetPad(i);

		cur_pad->cd();
		cur_pad->SetLeftMargin(0.1);
		cur_pad->SetRightMargin(0.002);



		cout << i << "   " << i-1 << endl;

		const int in =  i -1;

		cout << in << endl;

		TH1F* h_tgt_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("phi_real_var");
		TH1F* h_tgt[in] = (TH1F*) h_tgt_tmp->Clone("h_sim_clone"); 

		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("missmass_u_sim");
		TH1F* h_sim[in] = (TH1F*) h_sim_tmp->Clone("h_sim_clone");


		
	 	h_tgt[in]->DrawCopy();

//	 	h_tgt[in]->DrawClone();

//		h_sim[in]->Scale(8);
  
		h_sim[in]->SetLineColor(2);

	 	h_sim[in]->DrawCopy("same");

		peak1 = h_sim[in];
		
// 		TF1 *fit_u_fun;
// 
// 		if (i==1) {
// 			 fit_u_fun = new TF1("fit_u_fun",ftotal_all, fit_low, fit_high ,5);
// 		} else if (i==2) {
// 			 fit_u_fun = new TF1("fit_u_fun",ftotal_all, fit_low, fit_high ,5);
// 		} else {
// 			 fit_u_fun = new TF1("fit_u_fun",ftotal_all, fit_low, fit_high ,5);
// 		}


// 		fit_local_lo = fit_low;
// 		fit_local_hi = fit_high;
// 


		Double_t fit_local_lo;
		Double_t fit_local_hi;

		if (i==1) {

			fit_local_lo = fit_low;
			fit_local_hi = fit_low+0.27;


		} else if (i==2) {

			fit_local_lo = fit_low + 0.18;
			fit_local_hi = fit_high - 0.12;

		} else {

			fit_local_lo  = fit_high - 0.22;
			fit_local_hi = fit_high-0.01;

		}
// 

		


// 		TF1 *fit_u_fun = new TF1("fit_u_fun",ftotal_all, fit_low, fit_high ,5);
// 
// 		h_tgt[in]->Fit("fit_u_fun", "MUER");
		

//			TF1 *fit_u_fun = new TF1("fit_u_fun", ftotal, fit_low, fit_high ,5);

	
		TF1 *fit_u_fun = new TF1("fit_u_fun", ftotal_all, fit_local_lo, fit_local_hi, 5);

		fit_u_fun->SetNpx(400);

		fit_u_fun->SetParameters(10.0, -50.0, 75, -36.0, 15);

//  			fit_u_fun->FixParameter(0, 11.7661);
//  			fit_u_fun->FixParameter(1, -54.3);
//  			fit_u_fun->FixParameter(2, 82.7825);
//  			fit_u_fun->FixParameter(3, -41.5146);
// 
//  			fit_u_fun->FixParameter(4, 200);
			

// 			fit_u_fun->SetParLimits(0, 0, 50);
// 			fit_u_fun->SetParLimits(1, -100, 0);
// 			fit_u_fun->SetParLimits(2, 0, 200);
// 			fit_u_fun->SetParLimits(3, -100, 0);
 			fit_u_fun->SetParLimits(4, 0, 60);
// 
// 

		h_tgt[in]->Fit("fit_u_fun", "MEURN");


//			h_tgt[in]->Fit("fit_u_fun", "VR+");

//			h_sim[in]->Scale(18.868);
	 	h_sim[in]->Draw("same");

//			fit_u_fun->Draw("");

	 	TH1F *h_sim_clone = (TH1F*) h_sim[in]->Clone();
		
		TF1 *fit_u_fun_pol = new TF1("fit_u_fun_pol", ftotal_br, fit_local_lo, fit_local_hi, 4);

		fit_u_fun_pol->SetLineColor(4);

//  		fit_u_fun_pol->FixParameter(0, fit_u_fun->GetParameter(0));
//  		fit_u_fun_pol->FixParameter(1, fit_u_fun->GetParameter(1));
//  		fit_u_fun_pol->FixParameter(2, fit_u_fun->GetParameter(2));
//  		fit_u_fun_pol->FixParameter(3, fit_u_fun->GetParameter(3));
// 
		fit_u_fun_pol->Draw("same");

		h_sim_clone->SetLineColor(6);
		h_sim_clone->Scale(fit_u_fun->GetParameter(4));		
 	 	h_sim_clone->DrawCopy("same");

 		h_sim_clone->Add(fit_u_fun_pol);		
		h_sim_clone->SetLineColor(1);
 	 	h_sim_clone->DrawCopy("same");

 	}


	cc5->Print(out_dir_str + "u_fit_" + file_name+".png");

}





/*--------------------------------------------------*/
/*--------------------------------------------------*/
/*--------------------------------------------------*/
// Just fit missing mass in t 
// Version: 1
// Description: 2nd order polynomial

void t_fit_2(TCanvas* cc3, TString file_name, Double_t fit_low_1, Double_t fit_low_2, Double_t fit_mid_1, Double_t fit_mid_2, Double_t fit_high_1, Double_t fit_high_2) {
	
	cout << "asddddddddddddddddd " << endl;

// 	cc3->Draw();
// 
 	TCanvas* cc5 = new TCanvas("cc5", "cc5", 1600, 800);
 	cc5->Divide(u_bin_num, 1, 0.003);
// 
 	TH1F* h_tgt[3];
 	TH1F* h_sim[3];

//	TF1* fit_u_fun[3];

// 
 	for (Int_t i = 1; i <=3; i++) {
// 	
		TPad* cur_pad = cc5->GetPad(i);

		cur_pad->cd();
		cur_pad->SetLeftMargin(0.1);
		cur_pad->SetRightMargin(0.002);



		cout << i << "   " << i-1 << endl;

		const int in =  i -1;

		cout << in << endl;

		TH1F* h_tgt_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("phi_real_var");
		TH1F* h_tgt[in] = (TH1F*) h_tgt_tmp->Clone("h_tgt_clone"); 

		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("missmass_u_sim");
		TH1F* h_sim[in] = (TH1F*) h_sim_tmp->Clone("h_sim_clone");


		
	 	h_tgt[in]->DrawCopy();

//	 	h_tgt[in]->DrawClone();

//		h_sim[in]->Scale(8);
  
		h_sim[in]->SetLineColor(2);

	 	h_sim[in]->DrawCopy("same");

		peak1 = h_sim[in];
		
// 		TF1 *fit_u_fun;
// 
// 		if (i==1) {
// 			 fit_u_fun = new TF1("fit_u_fun",ftotal_all, fit_low, fit_high ,5);
// 		} else if (i==2) {
// 			 fit_u_fun = new TF1("fit_u_fun",ftotal_all, fit_low, fit_high ,5);
// 		} else {
// 			 fit_u_fun = new TF1("fit_u_fun",ftotal_all, fit_low, fit_high ,5);
// 		}


// 		fit_local_lo = fit_low;
// 		fit_local_hi = fit_high;
// 


		Double_t fit_local_lo;
		Double_t fit_local_hi;

		if (i==1) {

			fit_local_lo = fit_low_1;
			fit_local_hi = fit_low_2;


		} else if (i==2) {

			fit_local_lo = fit_mid_1;
			fit_local_hi = fit_mid_2;

		} else {

			fit_local_lo = fit_high_1;
			fit_local_hi = fit_high_2;

		}
// 

		


// 		TF1 *fit_u_fun = new TF1("fit_u_fun",ftotal_all, fit_low, fit_high ,5);
// 
// 		h_tgt[in]->Fit("fit_u_fun", "MUER");
		

//			TF1 *fit_u_fun = new TF1("fit_u_fun", ftotal, fit_low, fit_high ,5);

	
//		TF1 *fit_u_fun = new TF1("fit_u_fun", ftotal_all, fit_local_lo, fit_local_hi, 5);
 		TF1 *fit_u_fun = new TF1("fit_u_fun", ftotal_all_2, fit_local_lo, fit_local_hi, 4);
//		fit_u_fun->SetNpx(400);

//		fit_u_fun->SetParameters(10.0, -50.0, 75, 15);

//  			fit_u_fun->FixParameter(0, 11.7661);
//  			fit_u_fun->FixParameter(1, -54.3);
//  			fit_u_fun->FixParameter(2, 82.7825);
//  			fit_u_fun->FixParameter(3, -41.5146);
// 
//  			fit_u_fun->FixParameter(4, 200);
			

// 			fit_u_fun->SetParLimits(0, 0, 50);
// 			fit_u_fun->SetParLimits(1, -100, 0);
// 			fit_u_fun->SetParLimits(2, 0, 200);
// 			fit_u_fun->SetParLimits(3, -100, 0);
 			fit_u_fun->SetParLimits(3, 0, 60);
// 
// 

		h_tgt[in]->Fit("fit_u_fun", "MEURN");


//			h_tgt[in]->Fit("fit_u_fun", "VR+");

//			h_sim[in]->Scale(18.868);
	 	h_sim[in]->Draw("same");

//			fit_u_fun->Draw("");

	 	TH1F *h_sim_clone = (TH1F*) h_sim[in]->Clone();
		
//		TF1 *fit_u_fun_pol = new TF1("fit_u_fun_pol", ftotal_br, fit_local_lo, fit_local_hi, 4);
		TF1 *fit_u_fun_pol = new TF1("fit_u_fun_pol", ftotal_br_2, fit_local_lo, fit_local_hi, 3);

		fit_u_fun_pol->SetLineColor(4);

 		fit_u_fun_pol->FixParameter(0, fit_u_fun->GetParameter(0));
 		fit_u_fun_pol->FixParameter(1, fit_u_fun->GetParameter(1));
 		fit_u_fun_pol->FixParameter(2, fit_u_fun->GetParameter(2));
// 		fit_u_fun_pol->FixParameter(3, fit_u_fun->GetParameter(3));

		fit_u_fun_pol->Draw("same");

		h_sim_clone->SetLineColor(6);
//		h_sim_clone->Scale(fit_u_fun->GetParameter(4));
		h_sim_clone->Scale(fit_u_fun->GetParameter(3));			
 	 	h_sim_clone->DrawCopy("same");

 		h_sim_clone->Add(fit_u_fun_pol);		
		h_sim_clone->SetLineColor(1);
 	 	h_sim_clone->DrawCopy("same");

 	}

	cc5->Print(out_dir_str + "t_fit_" + file_name+".png");


}





/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Reinitialization of the yield array

void Init_yield_array() {

	for(Int_t i = 0; i < sizeof(acc_yield)/sizeof(acc_yield[0]); i++) {

		acc_yield[i]     = 0.0;
		acc_sim_yield[i] = 0.0;

		acc_yield_err[i]     = 0.0;
		acc_sim_yield_err[i] = 0.0;
	}

}



/*--------------------------------------------------*/
/// Output of the yield array

void Output_yield_array() {

// 	for(Int_t i = 0; i < sizeof(acc_yield)/sizeof(acc_yield[0]); i++) {
// 		
// 		cout << acc_yield[i] << "    " << acc_sim_yield[i] << "    " << acc_yield[i]/acc_sim_yield[i] << endl;
// 
// 	}


	if(is_itt == true) {
	
		Int_t in;
	
		for(Int_t i = 0; i < u_bin_num; i++) {
			
			for(Int_t ii = 0; ii < phi_bin_num; ii++) {
	
				in = i * (phi_bin_num) + ii ;
			
				cout << i << "  " << ii << "  " << in << "  " << acc_yield[in] << "  " << acc_sim_yield[in] << "  " << acc_yield_err[in] + acc_sim_yield_err[in] << "  "  << ii+1 << "  " << i+1 << endl;

				Float_t exp_sim_ratio = acc_yield[in]/acc_sim_yield[in];
				Float_t exp_sim_ratio_err = sqrt(acc_yield_err[in]/(acc_yield[in]**2) + acc_sim_yield_err[i]/(acc_sim_yield[in]**2)) * (acc_yield[in]/acc_sim_yield[in]);




				if ( exp_sim_ratio != exp_sim_ratio) {

					exp_sim_ratio     = 0.0;
					exp_sim_ratio_err = 5.0;

					cout << exp_sim_ratio << " !!! is not a number!!! " << endl;

//					exit(0);

				}
				
				yield_file_out << fixed << setw(8) << setprecision(5) << exp_sim_ratio << "    " << exp_sim_ratio_err << "    "  << ii+1 << "    "  << i+1 << "     " << acc_yield[in] << "     " << acc_sim_yield[in] << endl;
	
			}
	
		}

	}





	yield_file_out.close();



 	for(Int_t i = 0; i < u_bin_num * phi_bin_num; i++) {
 
 		acc_yield[i]         = 0.0;
        acc_sim_yield[i]     = 0.0;
                             
 	    acc_yield_err[i]     = 0.0;
        acc_sim_yield_err[i] = 0.0;
 
 	}



//	exit(0);

}




