///*--------------------------------------------------*/
/// 
/// WL 
/// Date: 1/Feb/2016
/// Main function to determine experimental yield 
/// 
/// Omega experimental yield = Data - Background
/// where backgroumd includes Xphasespace and rho production
///
//
//

#include<vector>

Double_t intergration_limit  = 0.04;   // GeV
Double_t fitting_limit       = 0.04;   // percentage / 100


Double_t low_stats_limit     = 0.015;   // Normalized yield
// Double_t low_stats_limit     = 0.005;   // Normalized yield

Double_t int_hi_global = intergration_limit;
Double_t int_lo_global = intergration_limit;

Double_t fit_hi_global = fitting_limit;
Double_t fit_lo_global = fitting_limit;







void u_phi_fit_sim(TCanvas* cc3, TString file_name) {

	const Int_t u_bin_num_const   = u_bin_num;  
	const Int_t phi_bin_num_const = phi_bin_num;  
	const Int_t bin_num_const     = u_bin_num_const * phi_bin_num_const;

	ofstream fit_out_file;
	fit_out_file.open("u_bin_fit/fit_output.pl_" + file_name + ".dat", std::fstream::out); 

	TFile* file_out_root = new TFile(out_dir_str + "u_fit_phi_" + file_name+".root", "RECREATE");

	cout << "asddddddddddddddddd " << endl;

//	cc3->Draw();

//	TPad* cur_pad;

	TCanvas* cc4 = new TCanvas("cc4", "cc4", 1600, 800);
	cc4->Divide(phi_bin_num_const, u_bin_num_const, 0.003, 0.003);

	TCanvas* cc2 = new TCanvas("cc2", "cc2", 1600, 800);
	cc2->Divide(phi_bin_num_const, u_bin_num_const, 0.003, 0.003);

	TCanvas* cc6 = new TCanvas("cc6", "cc6", 1600, 800);
	cc6->Divide(phi_bin_num_const, u_bin_num_const, 0.003, 0.003);

	TCanvas* cc7 = new TCanvas("cc7", "cc7", 1600, 800);
	cc7->Divide(phi_bin_num_const, u_bin_num_const, 0.003, 0.003);


 	TH1F* h_tgt[bin_num_const];
 	TH1F* h_sim[bin_num_const];
 	TH1F* h_sim_rho[bin_num_const];
 	TH1F* h_sim_xphsp[bin_num_const];
 	TH1F* h_sim_eta[bin_num_const];
 	TH1F* h_sim_etap[bin_num_const];
 
	TPad *pp2;
	TPad *pp6;
	TPad *pp7;

	TPad *ppnew;
	TPad *ppnew1;
	TPad *ppnew2;
	TPad *ppnew3;

	TH1F* chi2_omega_acc = new TH1F("chi2_omega", "chi2_omega_check", 40, 0, 10);
	TH1F* chi2_bg_acc    = new TH1F("chi2_bg"   , "chi2_bg_check",    40, 0, 10);

// 	TH1F* data_omega_sum = new TH1F("data_omega_sum", "data_omega_sum", 100, 0.1, 1.1);
// 	TH1F* sim_omega_sum  = new TH1F("sim_omega_sum", "sim_omega_sum", 100, 0.1, 1.1);
// 

	TH1F* raw_data_sum;
	TH1F* data_omega_sum;

    TH1F* sim_omega_sum;
	TH1F* sim_rho_sum;
	TH1F* sim_xphsp_sum;
	TH1F* sim_eta_sum;
	TH1F* sim_etap_sum;

	Int_t iittt = 1;

	bool is_refit = false;

	
	Int_t refit_itt = 0;

	for (Int_t i = 1; i <= u_bin_num*phi_bin_num; i++) {

		// fit_out_file << fixed << setprecision(4) << Int_t((i-1)/phi_bin_num) + 1 << "  " << (i-1)%phi_bin_num + 1 << endl;


		if (is_refit) {

//			i = i - 1;
//			exit(0);
//

			if (refit_itt == 1) {


//				exit(0);
				refit_itt = 0;
				is_refit = false;
				continue;

			}

			refit_itt++;

		}


		TString f_str_temp;
		f_str_temp.Form("_%i_%i", Int_t((i-1)/phi_bin_num) + 1, (i-1)%phi_bin_num + 1);
		f_str_temp = "u_fit_phi_" + file_name + f_str_temp;

		TPad*cur_pad = (TPad*)cc4->GetPad(i);

		cur_pad->cd();
		cur_pad->SetLeftMargin(0.1);
		cur_pad->SetRightMargin(0.002);

		// cout << i << "   " << i-1 << endl;

		const int in =  i -1;

//		cout << in << endl;

		cout << endl << "/*--------------------------------------------------*/" << endl;
		cout << "u_bin: " << Int_t((i-1)/phi_bin_num) + 1 << "   Phi_bin: " << (i-1)%phi_bin_num + 1 
			 << "    Overall counter: " << in<< endl;
	
//		TH1F* h_tgt[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
//		TH1F* h_sim[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");

 		TH1F* h_tgt_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("u_phi_target");
 		TH1F* h_tgt[in] = (TH1F*) h_tgt_tmp->Clone("h_tgt_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("u_phi_sim_omega");
 		TH1F* h_sim[in] = (TH1F*) h_sim_tmp->Clone("h_sim_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("u_phi_sim_rho");
 		TH1F* h_sim_rho[in] = (TH1F*) h_sim_tmp->Clone("h_sim_rho_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("u_phi_sim_xphsp");
 		TH1F* h_sim_xphsp[in] = (TH1F*) h_sim_tmp->Clone("h_sim_xphsp_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("u_phi_sim_eta");
 		TH1F* h_sim_eta[in] = (TH1F*) h_sim_tmp->Clone("h_sim_eta_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("u_phi_sim_etap");
 		TH1F* h_sim_etap[in] = (TH1F*) h_sim_tmp->Clone("h_sim_etap_clone");



		h_tgt[in]->SetTitle(f_str_temp);

	 	h_tgt[in]->DrawCopy();  


//	 	h_sim[in]->DrawCopy("same");

		peak_omega = h_sim[in];
		peak_rho   = h_sim_rho[in];
		peak_xphsp = h_sim_xphsp[in];
		peak_eta   = h_sim_eta[in];
		peak_etap  = h_sim_etap[in];

		Float_t omega_scale_fac = 0.0;
		Float_t rho_scale_fac   = 0.0; 
		Float_t xphsp_scale_fac = 0.0;
		Float_t eta_scale_fac   = 0.0;
		Float_t etap_scale_fac  = 0.0;

		Float_t omega_scale_fac_err = 0.0;
		Float_t rho_scale_fac_err   = 0.0; 
		Float_t xphsp_scale_fac_err = 0.0;
		Float_t eta_scale_fac_err   = 0.0;
		Float_t etap_scale_fac_err  = 0.0;




		Double_t intg_lo        = 0.0;
        Double_t intg_hi        = 0.0;

		Double_t fit_local_lo   = 0.0;
		Double_t fit_local_hi   = 0.0;

	    Double_t omega_yield    = 0.0;   
        Double_t omega_yield_err = 0.0;

		Double_t content_total = 0.0;
		Double_t content_total_error_sq = 0.0;
	
// 		fit_local_lo = 0.5;
// 		fit_local_hi = 0.9;

		TPad* tar_pad = (TPad*) cc3->GetPad(i);

		bool rad_only = false;
		bool low_stat = false;
		bool off_focal = false;
		bool is_excluded = false;

		
		Float_t omega_mass_ratio;

		omega_mass_ratio = h_tgt[in]->Integral(0, h_tgt[in]->FindBin(0.8))/ h_tgt[in]->Integral();


// 		if( (h_tgt[in]->GetMean() - h_tgt[in]->GetRMS()) > 0.783 && omega_mass_ratio < 0.1) {
// 			rad_only = true;
// 		}


		if( (peak_omega->GetMean() - peak_omega->GetRMS()) > 0.783) {
			rad_only = true;
//			rad_only = false;
		}

		if( h_tgt[in]->Integral(0, -1) < low_stats_limit) {

//		if( h_tgt[in]->Integral(0, -1) < 0.04 ) {

//		if( h_tgt[in]->Integral(0, -1) < 0.05 ) {
			low_stat = true;
		}

		if( ( peak_omega->GetRMS() + peak_omega->GetMean()) < 0.783 ) {
			off_focal = true;
		}
// 



		h_sim[in]->DrawCopy("hsame");


//		if( h_tgt[in]->Integral(0,-1) > 0.001)   
//		if( h_tgt[in]->Integral(0,-1) > 0.005 && !rad_only)   
		if( !low_stat && !rad_only && !off_focal) {  


			if (iittt == 1) {

				raw_data_sum   = (TH1F*) h_tgt[in]->Clone("raw_data_sum")  ;
				data_omega_sum = (TH1F*) h_tgt[in]->Clone("data_omega_sum") ;
				sim_omega_sum  = (TH1F*) h_tgt[in]->Clone("sim_omega_sum")  ;
				sim_rho_sum    = (TH1F*) h_tgt[in]->Clone("sim_rho_sum") ;
				sim_xphsp_sum  = (TH1F*) h_tgt[in]->Clone("sim_xphsp_sum")  ;				
				sim_eta_sum    = (TH1F*) h_tgt[in]->Clone("sim_eta_sum") ;
				sim_etap_sum   = (TH1F*) h_tgt[in]->Clone("sim_etap_sum")  ;
	
				raw_data_sum   ->Reset();
				data_omega_sum ->Reset();
				sim_omega_sum  ->Reset();
				sim_rho_sum    ->Reset();
				sim_xphsp_sum  ->Reset();			
				sim_eta_sum    ->Reset();
				sim_etap_sum   ->Reset();

				iittt =0;
			}

//		if( !low_stat && !rad_only)   

			cout << "Target Raw Yield: " << h_tgt[in]->Integral(0,-1) << endl;

//  		fit_local_lo = Find_low_range_xphsp(tar_pad);
//  		fit_local_hi = Find_high_range_xphsp(tar_pad);

			/*--------------------------------------------------*/
// 	 		fit_local_lo = Find_low_range(tar_pad);
//  			fit_local_hi = Find_high_range(tar_pad);


//
//			cout << "asdasdasda "<< fit_local_lo << "    " << fit_local_hi << endl;

//			fit_local_lo = Find_int_lo(tar_pad);
//			fit_local_hi = Find_int_hi(tar_pad);

//		    fit_local_lo = 0.65;
//		    fit_local_hi = 0.75;


// 			fit_local_lo = 0.65;
// 			fit_local_hi = 0.9;



// 			fit_local_lo = Find_low_range_full(tar_pad);
// 			fit_local_hi = Find_high_range_full(tar_pad);
// 
// 			if (is_refit) {
// 				fit_local_lo = Find_low_range_full(tar_pad);
// 				fit_local_hi = Find_high_range_full(tar_pad);
// 			}


//			fit_local_lo = Find_int_lo(tar_pad);
//			fit_local_hi = Find_int_hi(tar_pad);


			
			if (is_refit) {

// 			  fit_local_hi = Find_high_range_full(tar_pad) - 0.02;
//			  fit_local_lo = h_tgt[in]->GetBinCenter(h_tgt[in]->FindFirstBinAbove(0.0)) + 0.02;


	  		   fit_local_lo = Find_low_range_full(tar_pad) + 0.005;
 		   	   fit_local_hi = Find_high_range_full(tar_pad) - 0.005;


//	   	 	  fit_local_lo = Find_int_lo_refit(tar_pad);
//              fit_local_hi = Find_int_hi_refit(tar_pad);



			} else {

//	 			fit_local_lo = 0.65;
// 				fit_local_hi = 0.9;
	
// 				fit_local_lo = Find_int_lo(tar_pad);
// 				fit_local_hi = Find_int_hi(tar_pad);

// 				fit_local_lo = Find_low_range_full(tar_pad);

// 		   	   fit_local_hi = Find_high_range_full(tar_pad) - 0.01;
//	  		   fit_local_lo = h_tgt[in]->GetBinCenter(h_tgt[in]->FindFirstBinAbove(0.0)) + 0.01;



	  		   fit_local_lo = Find_low_range_full(tar_pad);
 		   	   fit_local_hi = Find_high_range_full(tar_pad);





//			  fit_local_hi = h_tgt[in]->GetBinCenter(h_tgt[in]->FindLastBinAbove(0.0));


// 			  if ( h_tgt[in]->GetBinCenter(h_tgt[in]->FindLastBinAbove(0.0)) < 0.9) {
//               
// 				fit_local_hi = h_tgt[in]->GetBinCenter(h_tgt[in]->FindLastBinAbove(0.0));
//               
// 			  } else {
// 
//               	fit_local_hi = 0.9;
// 
// 			  }

			}

	
			TF1 *fit_u_fun = new TF1("fit_u_fun", fun_sim_total, fit_local_lo, fit_local_hi, 5);
	
	 		fit_u_fun->SetParameters(1, 40, 5000, 0.1, 10);
	
			fit_u_fun->SetParLimits(0, 0, 100);
			fit_u_fun->SetParLimits(1, 0, 100000);
			fit_u_fun->SetParLimits(2, 0, 100000);
			fit_u_fun->SetParLimits(3, 0, 100000);
			fit_u_fun->SetParLimits(4, 0, 100000);
			




			TString fit_fun_used_str;
			
			fit_fun_used_str = "Bg:";
 			fit_fun_used_str = fit_fun_used_str + "#rho";
			fit_fun_used_str = fit_fun_used_str + "+xp";




		
// 			fit_u_fun->FixParameter(1, 0);
// 			fit_u_fun->FixParameter(2, 0);
// 			fit_u_fun->FixParameter(3, 0);
// 			fit_u_fun->FixParameter(4, 0);

			if (include_etap(tar_pad)) {
				fit_u_fun->FixParameter(3, 0);
 				fit_fun_used_str = fit_fun_used_str + "+etap";
			} else {
				fit_u_fun->FixParameter(4, 0);
 				fit_fun_used_str = fit_fun_used_str + "+#eta";
			}




// 			if (h_tgt[in]->GetBinCenter(h_tgt[in]->FindFirstBinAbove(0)) < 0.50 ) {
// 
// 				fit_u_fun->ReleaseParameter(3);
// 				fit_fun_used_str = fit_fun_used_str + "+#eta";
// 
// 			} else if (h_tgt[in]->GetMean(1) > 0.85 ) {
// 				
// 				fit_u_fun->ReleaseParameter(4);
// 				fit_fun_used_str = fit_fun_used_str + "+#etap";
// 					
//     		}
//
//
// 
// 			if (include_etap) {
// 
//  				fit_fun_used_str = fit_fun_used_str + "+#etap";
// 			}
 

	 
			h_tgt[in]->Fit("fit_u_fun", "MNR");

			fit_u_fun->SetParameters(fit_u_fun->GetParameter(0), 
									 fit_u_fun->GetParameter(1), 
									 fit_u_fun->GetParameter(2), 
									 fit_u_fun->GetParameter(3), 
									 fit_u_fun->GetParameter(4));
			
			h_tgt[in]->Fit("fit_u_fun", "MNR");
	
		 	TH1F *h_sim_omega_clone = (TH1F*) h_sim[in]->Clone();
		 	TH1F *h_sim_rho_clone   = (TH1F*) h_sim_rho[in]->Clone();
		 	TH1F *h_sim_xphsp_clone = (TH1F*) h_sim_xphsp[in]->Clone();
		 	TH1F *h_sim_eta_clone   = (TH1F*) h_sim_eta[in]->Clone();
		 	TH1F *h_sim_etap_clone  = (TH1F*) h_sim_etap[in]->Clone();
			
			omega_scale_fac = fit_u_fun->GetParameter(0);
			rho_scale_fac   = fit_u_fun->GetParameter(1);
			xphsp_scale_fac = fit_u_fun->GetParameter(2);
			eta_scale_fac   = fit_u_fun->GetParameter(3);
			etap_scale_fac  = fit_u_fun->GetParameter(4);
	
			omega_scale_fac_err = fit_u_fun->GetParError(0);
            rho_scale_fac_err   = fit_u_fun->GetParError(1);
            xphsp_scale_fac_err = fit_u_fun->GetParError(2);
            eta_scale_fac_err   = fit_u_fun->GetParError(3);
            etap_scale_fac_err  = fit_u_fun->GetParError(4);

			h_sim_omega_clone->Scale(omega_scale_fac);		
			h_sim_rho_clone  ->Scale(rho_scale_fac);		
			h_sim_xphsp_clone->Scale(xphsp_scale_fac);		
			h_sim_eta_clone  ->Scale(eta_scale_fac);		
			h_sim_etap_clone ->Scale(etap_scale_fac);		
			
		 	TH1F *h_sim_all = (TH1F*) h_sim_xphsp_clone->Clone();
	
	 		h_sim_all->Add(h_sim_omega_clone);		
	 		h_sim_all->Add(h_sim_rho_clone);		
	 		h_sim_all->Add(h_sim_eta_clone);		
	 		h_sim_all->Add(h_sim_etap_clone);		



			Float_t fit_tgt_diff_per = (h_tgt[in]->Integral(0,-1)-h_sim_all->Integral(0, -1))/h_tgt[in]->Integral(0, -1);

			cout << h_tgt[in]->Integral(0,-1) << "   " << h_sim_all->Integral(0, -1) << "   "<< fit_tgt_diff_per << endl;

//  			if( fabs(fit_tgt_diff_per) > 0.4 ) {
//  
// //				fit_local_lo = h_tgt[in]->GetBinCenter(h_tgt[in]->GetMaximumBin())-0.05;
// //				fit_local_hi = h_tgt[in]->GetBinCenter(h_tgt[in]->GetMaximumBin())+0.05;
// 
// //				fit_local_lo = 0.65;
// //				fit_local_hi = 0.75;
// 
// 
// 				fit_local_lo = fit_local_lo + 0.05;
// 				fit_local_hi = fit_local_hi - 0.05;
// 
// 
// 
// 				fit_u_fun->SetRange(fit_local_lo, fit_local_hi);
//  
//  				h_tgt[in]->Fit("fit_u_fun", "MNR");
//  
// 
// 		 		h_sim_omega_clone = (TH1F*) h_sim[in]->Clone();
// 		 		h_sim_rho_clone   = (TH1F*) h_sim_rho[in]->Clone();
// 		 		h_sim_xphsp_clone = (TH1F*) h_sim_xphsp[in]->Clone();
// 		 		h_sim_eta_clone   = (TH1F*) h_sim_eta[in]->Clone();
// 		 		h_sim_etap_clone  = (TH1F*) h_sim_etap[in]->Clone();
// 				
// 				omega_scale_fac = fit_u_fun->GetParameter(0);
// 				rho_scale_fac   = fit_u_fun->GetParameter(1);
// 				xphsp_scale_fac = fit_u_fun->GetParameter(2);
// 				eta_scale_fac   = fit_u_fun->GetParameter(3);
// 				etap_scale_fac  = fit_u_fun->GetParameter(4);
// 	
// 				h_sim_omega_clone->Scale(fit_u_fun->GetParameter(0));		
// 				h_sim_rho_clone->Scale(fit_u_fun->GetParameter(1));		
// 				h_sim_xphsp_clone->Scale(fit_u_fun->GetParameter(2));		
// 				h_sim_eta_clone->Scale(fit_u_fun->GetParameter(3));		
// 				h_sim_etap_clone->Scale(fit_u_fun->GetParameter(4));		
// 				
// 		 		h_sim_all = (TH1F*) h_sim_xphsp_clone->Clone();
// 	
// 	 			h_sim_all->Add(h_sim_omega_clone);		
// 	 			h_sim_all->Add(h_sim_rho_clone);		
// 	 			h_sim_all->Add(h_sim_eta_clone);		
// 	 			h_sim_all->Add(h_sim_etap_clone);		
//  
// 				fit_fun_used_str = fit_fun_used_str + "+R";
// 
// 
//  			}

			h_sim_omega_clone->SetStats(0);
			h_sim_rho_clone->SetStats(0);
			h_sim_xphsp_clone->SetStats(0);
			h_sim_eta_clone->SetStats(0);
			h_sim_etap_clone->SetStats(0);
	
			h_sim[in]->SetStats(0);
			h_sim[in]->SetFillColor(2);
	
			h_sim[in]->DrawCopy("hsame");
	
	 		h_sim_xphsp_clone->SetLineColor(3);
	 
	 		h_sim_rho_clone->SetLineColor(4); 		
	 
	  	 	h_sim_xphsp_clone->DrawCopy("hsame");
	  	 	h_sim_rho_clone->DrawCopy("hsame");
	  	 	h_sim_omega_clone->DrawCopy("hsame");
	 
	  	 	h_sim_eta_clone->DrawCopy("hsame");
	  	 	h_sim_etap_clone->DrawCopy("hsame");
	
			h_sim_all->SetLineColor(6);
	 	 	h_sim_all->DrawCopy("hsame");

//			exit(0);


	
	   		TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
	   		TH1* h = new TH1F("", "", 1, 0, 1);
	   		leg->AddEntry(h_sim_xphsp_clone, "Sim xphsp", "l");
	   		leg->AddEntry(h_sim_rho_clone,   "Sim rho",   "l");
	   		leg->AddEntry(h_sim_omega_clone, "Sim omega", "l");
	   		leg->AddEntry(h_sim_eta_clone,   "Sim eta",   "l");
	   		leg->AddEntry(h_sim_etap_clone,  "Sim etap",  "l");
	   		leg->AddEntry(h_sim_all,         "Sim All",   "l");
	   		leg->AddEntry(h_tgt[in],         "Data",    "epl");
	
			leg->SetFillColor(0);
			leg->Draw();
	
			// cout << h_tgt_tmp->GetMean(1) << "   "<< fit_local_lo << "   " << fit_local_hi  << endl;
	
			/*--------------------------------------------------*/
			/// Data output
			//
	
			TString range_str;
			range_str.Form("Range: %f - %f ", fit_local_lo, fit_local_hi);
	
			TString mean_str;
			mean_str.Form("Mean: %f RMS: %f Sim: %f ", h_tgt[in]->GetMean(1), h_tgt[in]->GetRMS(1), h_sim[in]->GetBinCenter(h_sim[in]->GetMaximumBin()));
	
			TString ratio_str;
			ratio_str.Form("w_f: %f, r_f: %f, x_f: %f", fit_u_fun->GetParameter(0), fit_u_fun->GetParameter(1), fit_u_fun->GetParameter(2));
	
	 		TText fitting_range;
	
	 		fitting_range.SetNDC();
	 		fitting_range.SetTextSize(0.04);
	
	 		fitting_range.DrawText(0.12, 0.85, range_str);
	 		fitting_range.DrawText(0.12, 0.80, mean_str);
	 		fitting_range.DrawText(0.12, 0.75, ratio_str);
	
			ratio_str.Form("e_f: %f, ep_f: %f", fit_u_fun->GetParameter(3), fit_u_fun->GetParameter(4));
	 		fitting_range.DrawText(0.12, 0.70, ratio_str);
	
			Int_t lo_limit_bin;
			Int_t hi_limit_bin;
	
	// 		///*--------------------------------------------------*/
	// 		//
	// 		Double_t intg_lo = h_sim_omega_clone->GetBinCenter(h_sim_omega_clone->GetMaximumBin()) - 0.035;
	// 		Double_t intg_hi = h_sim_omega_clone->GetBinCenter(h_sim_omega_clone->GetMaximumBin()) + 0.035;
	// 
	// 		if (h_sim_omega_clone->GetMaximumBin() < 10) {
	// 
	// 			cout << "aaaaaaaaaaaa "<< h_sim_omega_clone->GetMaximumBin() << endl; 
	// 		
	// 			intg_lo =  0.783 - 0.05;
	// 			intg_hi =  0.783 + 0.05;
	//  
	// 		}
	
			///*--------------------------------------------------*/
	
	 		// Double_t intg_lo = h_sim_omega_clone->GetBinCenter(h_sim_omega_clone->FindFirstBinAbove(0.0));
	 		// Double_t intg_hi = h_sim_omega_clone->GetBinCenter(h_sim_omega_clone->FindLastBinAbove(0.0));
	
	//  		Double_t intg_lo = fit_local_lo;
	
// 			intg_lo = Find_int_lo(tar_pad);
// 	  		intg_hi = Find_int_hi(tar_pad);

	 		intg_lo =  0.783 - int_lo_global;
	 		intg_hi =  0.783 + int_hi_global;

//	 		intg_hi =   fit_local_hi;
//	 		intg_lo =   fit_local_lo;



			
			cout << "Int_limit "<< intg_lo << "    " << intg_hi << endl;
//			exit(0);



	
			lo_limit_bin = h_sim_omega_clone->FindBin(intg_lo);
	        hi_limit_bin = h_sim_omega_clone->FindBin(intg_hi);
	
			Double_t sim_omega_yield = h_sim[in]->Integral(lo_limit_bin, hi_limit_bin);
	
			omega_yield = h_tgt[in]->Integral(lo_limit_bin, hi_limit_bin) - h_sim_rho_clone->Integral(lo_limit_bin, hi_limit_bin) - h_sim_xphsp_clone->Integral(lo_limit_bin, hi_limit_bin) - h_sim_eta_clone->Integral(lo_limit_bin, hi_limit_bin) - h_sim_etap_clone->Integral(lo_limit_bin, hi_limit_bin);
	




/*--------------------------------------------------*/
/// Error propagation
/*--------------------------------------------------*/


			Float_t data_con         = 0.0;
			Float_t data_con_err_sq  = 0.0;

			Float_t rho_con          = 0.0;
			Float_t rho_con_err_sq   = 0.0;

			Float_t xphsp_con        = 0.0;
			Float_t xphsp_con_err_sq = 0.0;

			Float_t eta_con          = 0.0;
			Float_t eta_con_err_sq   = 0.0;
                                     
			Float_t etap_con         = 0.0;
			Float_t etap_con_err_sq  = 0.0;

			for (Int_t bin_itt = lo_limit_bin; bin_itt <= hi_limit_bin; bin_itt++) {
	
				Double_t bin_con = h_tgt[in]->GetBinContent(bin_itt);
				Double_t bin_err = h_tgt[in]->GetBinError(bin_itt);
	
				data_con         = data_con + bin_con;
				data_con_err_sq  = data_con_err_sq + bin_err**2;

			}
	
			for (Int_t bin_itt = lo_limit_bin; bin_itt <= hi_limit_bin; bin_itt++) {
	
				Double_t bin_con = h_sim_rho[in]->GetBinContent(bin_itt) * rho_scale_fac;
				Double_t bin_err = h_sim_rho[in]->GetBinError(bin_itt) * rho_scale_fac;
	
				rho_con          = rho_con + bin_con;
				rho_con_err_sq   = rho_con_err_sq + bin_err**2;

			}
	
			for (Int_t bin_itt = lo_limit_bin; bin_itt <= hi_limit_bin; bin_itt++) {
	
				Double_t bin_con = h_sim_xphsp[in]->GetBinContent(bin_itt) * xphsp_scale_fac;
				Double_t bin_err = h_sim_xphsp[in]->GetBinError(bin_itt) * xphsp_scale_fac;
	
				xphsp_con        = xphsp_con + bin_con;
				xphsp_con_err_sq = xphsp_con_err_sq + bin_err**2;
	
			}

			for (Int_t bin_itt = lo_limit_bin; bin_itt <= hi_limit_bin; bin_itt++) {
	
				Double_t bin_con = h_sim_eta[in]->GetBinContent(bin_itt) * eta_scale_fac;
				Double_t bin_err = h_sim_eta[in]->GetBinError(bin_itt) * eta_scale_fac;
	
				eta_con        = eta_con + bin_con;
				eta_con_err_sq = eta_con_err_sq + bin_err**2;
	
			}

			for (Int_t bin_itt = lo_limit_bin; bin_itt <= hi_limit_bin; bin_itt++) {
	
				Double_t bin_con = h_sim_etap[in]->GetBinContent(bin_itt) * etap_scale_fac;
				Double_t bin_err = h_sim_etap[in]->GetBinError(bin_itt) * etap_scale_fac;
	
				etap_con        = etap_con + bin_con;
				etap_con_err_sq = etap_con_err_sq + bin_err**2;
	
			}


	

			////*--------------------------------------------------*/
			//// Safe Guard Statement
//			if (omega_yield < 0.00005 || omega_scale_fac < 0.00005 ) 
			if (omega_scale_fac < 0.0001 ) {

//				omega_yield     = 0.0;
//				omega_yield_err = 0.0;
	
			} else { 

	////*--------------------------------------------------*/
	//// Original
	//			omega_yield_err = sqrt(xphsp_con_err_sq/xphsp_con**2)*omega_yield;
	
		
	//			omega_yield_err  = sqrt( data_con_err_sq/pow(data_con,2) + rho_con_err_sq/pow(rho_con,2) + xphsp_con_err_sq/pow(xphsp_con,2) + (xphsp_scale_fac_err/xphsp_scale_fac)**2 +  (rho_scale_fac_err/rho_scale_fac)**2) * omega_yield;
	


	
	//			cout << omega_yield_err << "  "<< sqrt(data_con_err_sq/pow(data_con,2)) << "  " << sqrt(rho_con_err_sq/pow(rho_con,2)) << "  " << sqrt(xphsp_con_err_sq/pow(xphsp_con,2)) << "  "<< rho_scale_fac_err/rho_scale_fac << "  " << xphsp_scale_fac_err/xphsp_scale_fac << endl;
	
	//			cout << omega_yield_err << "  "<< sqrt(data_con_err_sq/pow(data_con,2)) << "  " << sqrt(rho_con_err_sq/pow(rho_con,2)) << "  " << sqrt(xphsp_con_err_sq/pow(xphsp_con,2)) << "  " << omega_scale_fac_err/omega_scale_fac << endl;
	
	//			exit(0);


	//			omega_yield_err  = sqrt( data_con_err_sq/pow(data_con,2) + rho_con_err_sq/pow(rho_con,2) + xphsp_con_err_sq/pow(xphsp_con,2) + (omega_scale_fac_err/omega_scale_fac)**2 ) * omega_yield;

	//			omega_yield_err  = sqrt( data_con_err_sq/pow(data_con,2) + rho_con_err_sq/pow(rho_con,2) + xphsp_con_err_sq/pow(xphsp_con,2) ) * omega_yield;



				cout << endl << "Omega Yield check: " << data_con << "  " << rho_con 
				     << "  " << xphsp_con << "  " << etap_con << "  " << etap_con << endl;



				cout << "Omega error check: " << omega_yield_err << "  " << data_con_err_sq 
					 << "  " << rho_con_err_sq << "  " << xphsp_con_err_sq << "  " << etap_con_err_sq 
					 << "  " << etap_con_err_sq << endl << endl;


				omega_yield_err  = sqrt(data_con_err_sq + rho_con_err_sq + xphsp_con_err_sq + eta_con_err_sq + etap_con_err_sq);


				
		 		TH1F* h_sim_omega_norm_clone = (TH1F*) h_sim[in]->Clone();
		 		TH1F* h_sim_omega_hi_clone   = (TH1F*) h_sim[in]->Clone();
		 		TH1F* h_sim_omega_lo_clone   = (TH1F*) h_sim[in]->Clone();


				h_sim_omega_norm_clone->Scale(omega_scale_fac);
				h_sim_omega_lo_clone  ->Scale(omega_scale_fac + omega_scale_fac_err); 
                h_sim_omega_hi_clone  ->Scale(omega_scale_fac - omega_scale_fac_err); 




				Float_t yield_bg_norm = h_tgt[in]->Integral(lo_limit_bin, hi_limit_bin) -   h_sim_omega_norm_clone->Integral(lo_limit_bin, hi_limit_bin); 
				Float_t yield_bg_lo   = h_tgt[in]->Integral(lo_limit_bin, hi_limit_bin) -   h_sim_omega_lo_clone->Integral(lo_limit_bin, hi_limit_bin);
				Float_t yield_bg_hi   = h_tgt[in]->Integral(lo_limit_bin, hi_limit_bin) -   h_sim_omega_hi_clone->Integral(lo_limit_bin, hi_limit_bin);



				Float_t yield_bg_rho_xphsp = h_sim_rho_clone->Integral(lo_limit_bin, hi_limit_bin) +   h_sim_xphsp_clone->Integral(lo_limit_bin, hi_limit_bin); 


				Float_t bg_error = (yield_bg_hi - yield_bg_lo)/2.0; 


				
//				cout << "Yield Check:    " << h_tgt[in]->Integral(lo_limit_bin, hi_limit_bin) << "  " << h_sim_omega_norm_clone->Integral(lo_limit_bin, hi_limit_bin) << "  " <<  h_sim_omega_lo_clone->Integral(lo_limit_bin, hi_limit_bin) << "  " << h_sim_omega_hi_clone->Integral(lo_limit_bin, hi_limit_bin) << endl;

//				cout << "Error Check:    " <<  omega_scale_fac_err/omega_scale_fac << "    " << bg_error/yield_bg_norm << endl;

//				cout << "BG Check: " << yield_bg_lo << "  " << yield_bg_norm  << "  " << yield_bg_hi << "  "<< yield_bg_rho_xphsp << endl;
 

//				omega_yield_err  = sqrt( data_con_err_sq/pow(data_con,2) + rho_con_err_sq/pow(rho_con,2) + xphsp_con_err_sq/pow(xphsp_con,2) + (omega_scale_fac_err/omega_scale_fac)**2 ) * omega_yield;


//				exit(0);


			}










/*--------------------------------------------------*/
/// End
/*--------------------------------------------------*/


			Double_t x_coord = 0.7;
			Double_t y_coord = 0.7;
	
			TText sim_data;
	 		sim_data.SetNDC();
	 		sim_data.SetTextSize(0.04);
	
			TString sim_data_str;
	
			sim_data_str.Form("Data All: %f", h_tgt[in]->Integral(lo_limit_bin, hi_limit_bin));	
			sim_data.DrawText(x_coord, y_coord, sim_data_str);	
	
			sim_data_str.Form("Sim All: %f", h_sim_xphsp_clone->Integral(lo_limit_bin, hi_limit_bin)+ h_sim_rho_clone->Integral(lo_limit_bin, hi_limit_bin) +  h_sim_eta_clone->Integral(lo_limit_bin, hi_limit_bin) + h_sim_etap_clone->Integral(lo_limit_bin, hi_limit_bin) +  sim_omega_yield );
		
			sim_data.DrawText(x_coord, y_coord-0.05, sim_data_str);	
	
			sim_data_str.Form("xphsp: %f", h_sim_xphsp_clone->Integral(lo_limit_bin, hi_limit_bin));	
			sim_data.DrawText(x_coord, y_coord-0.10, sim_data_str);	
	
			sim_data_str.Form("rho: %f", h_sim_rho_clone->Integral(lo_limit_bin, hi_limit_bin));	
			sim_data.DrawText(x_coord, y_coord-0.15, sim_data_str);	
	
			sim_data_str.Form("eta: %f", h_sim_eta_clone->Integral(lo_limit_bin, hi_limit_bin));	
			sim_data.DrawText(x_coord, y_coord-0.20, sim_data_str);	
	
			sim_data_str.Form("etap: %f", h_sim_etap_clone->Integral(lo_limit_bin, hi_limit_bin));	
			sim_data.DrawText(x_coord, y_coord-0.25, sim_data_str);	
	
			sim_data_str.Form("Omega Sim: %f", sim_omega_yield);
			sim_data.DrawText(x_coord, y_coord-0.30, sim_data_str);
	
			sim_data_str.Form("Omega Exp: %f", omega_yield);	
			sim_data.DrawText(x_coord, y_coord-0.35, sim_data_str);	
	
			///*--------------------------------------------------*/
			/// Outputing the fitting Status on the plot
			/// Color code: OK or Converged fits are in Blue
			/// 			Failed Fits are in Red!
	
			TText fitting_result;
	 		fitting_result.SetNDC();
	 		fitting_result.SetTextSize(0.06);
	
			TString fitting_result_str = "Fit: " + gMinuit->fCstatu;	
	
			if (fitting_result_str.Contains("OK") || fitting_result_str.Contains("CONVERGED")) {
				fitting_result.SetTextColor(4);
			} else {
				fitting_result.SetTextColor(2);
			}

			fitting_result.DrawText(0.05, 0.30, fitting_result_str);	

			TLatex fitting_bg_str;
	 		fitting_bg_str.SetNDC();
	 		fitting_bg_str.SetTextSize(0.06);
	 		fitting_bg_str.SetTextColor(2);
			fitting_bg_str.DrawLatex(0.05, 0.25, fit_fun_used_str);	
	

			TLatex refit_str;
	 		refit_str.SetNDC();
	 		refit_str.SetTextSize(0.06);
			if (is_refit) {
				refit_str.DrawLatex(0.05, 0.20, "Refitted!");
			}	


	//		sim_data_str.Form("Rat: %f", omega_yield/sim_omega_yield);	
	//		sim_data.DrawText(x_coord, y_coord-0.35, sim_data_str);	
	
	// 		sim_data_str.Form("omega: %f", h_sim_omega_clone->Integral(0, -1));	
	// 		sim_data.DrawText(0.9, 0.40, sim_data_str);	
	
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	// 		acc_yield[in]     = acc_yield[in] + omega_yield ;
	// 		acc_sim_yield[in] = acc_sim_yield[in] + sim_omega_yield ;
	// 		acc_yield_err[in] =  acc_yield_err[in] + content_total_error_sq;
	
			/*--------------------------------------------------*/
			/// Draw fitting boundaries
	
			TLine boundary_line_fitting;
			boundary_line_fitting.SetLineWidth(2);
			boundary_line_fitting.SetLineStyle(7);
			boundary_line_fitting.SetLineColor(2);
			boundary_line_fitting.DrawLine(fit_local_lo, 0, fit_local_lo, h_tgt[in]->GetMaximum()*0.8);
			boundary_line_fitting.DrawLine(fit_local_hi, 0, fit_local_hi, h_tgt[in]->GetMaximum()*0.8);
	
			TLine boundary_line_integral;
			boundary_line_integral.SetLineWidth(2);
			boundary_line_integral.SetLineStyle(7);
			boundary_line_integral.SetLineColor(4);
			boundary_line_integral.DrawLine(intg_lo, 0, intg_lo, h_tgt[in]->GetMaximum()*0.8);
			boundary_line_integral.DrawLine(intg_hi, 0, intg_hi, h_tgt[in]->GetMaximum()*0.8);
	
			cout << in << "  " << i << "  " << omega_scale_fac << "  " << rho_scale_fac << "  " 
						 << xphsp_scale_fac << "  " << intg_lo << "  " << intg_hi << "  " 
	                     <<  omega_yield << "  " << fit_local_lo << "   " << fit_local_hi  << endl;
	
	//		fit_out_file << fixed << setprecision(4) << setw(9) << (i-1)%u_bin_num << "  " << (i-1)%phi_bin_num 
	//          << "  " << omega_scale_fac  << "  " << rho_scale_fac << "  " << xphsp_scale_fac << "  " 
	//           << intg_lo << "  " << intg_hi << "  "<<  omega_yield << endl;
	
	//	 	fit_out_file << setiosflags(ios::fixed | ios::showpoint) << setprecision(5);
	
	// 		fit_out_file  << fixed << setprecision(5) << (i-1)%u_bin_num << "  " << (i-1)%phi_bin_num 
	//            << "  " << Form("%4.4f", omega_scale_fac)  << "  " <<  Form("%4.4f", rho_scale_fac) 
	//            << "  " << Form("%4.4f", xphsp_scale_fac) << "  " 
	//            << intg_lo << "  " << intg_hi << "  "<<  omega_yield << endl;



	




			///*--------------------------------------------------*/
			/// Beckground Plotting Check
	
			TH1F* h_tgt_sub = (TH1F*) h_tgt[in]->Clone("h_tgt_sub");
	
	// 		h_tgt_sub->Add(h_sim_rho_clone,-1);
	// 		h_tgt_sub->Add(h_sim_xphsp_clone,-1);
	
	
			h_tgt_sub->Add(h_sim_omega_clone, -1);
	
	// 		h_tgt_sub->Add(h_sim_eta_clone,-1);
	// 		h_tgt_sub->Add(h_sim_etap_clone,-1);
	
	  		pp2 = (TPad*) cc2->GetPad(i);
	 		pp2->cd();

			pp2->SetLeftMargin(0.1);
			pp2->SetRightMargin(0.002);
	
	
			h_tgt_sub->GetXaxis()->SetRangeUser(0.60, 1);
//			h_tgt_sub->GetXaxis()->SetRangeUser(fit_local_lo-0.05, fit_local_hi+0.05);

			TH1F* h_sim_bg = (TH1F*) h_sim_rho_clone->Clone();
			h_sim_bg->Add(h_sim_xphsp_clone, 1);
			h_sim_bg->Add(h_sim_eta_clone, 1);
			h_sim_bg->Add(h_sim_etap_clone, 1);


			TH1F* h_sim_bg_int = (TH1F*) h_sim_bg->Clone("h_tgt_zero_int");
			h_sim_bg_int->Reset();


			Float_t bg_chi2 = 0.0;
			Float_t bg_dof = 0.0;	

			for (Int_t bin_itt = lo_limit_bin; bin_itt <= hi_limit_bin; bin_itt++) {
	
				Double_t bin_con = h_tgt_sub->GetBinContent(bin_itt);
	

				Float_t data_err  = h_tgt[in]->GetBinError(bin_itt);

				if(h_tgt[in]->GetBinContent(bin_itt) == 0.0 && bin_con != 0.0) {
					data_err = 1.0;
				}

				Float_t omega_err = h_sim[in]->GetBinError(bin_itt) * omega_scale_fac;

				Double_t bin_err = sqrt(data_err**2 + omega_err**2);

//				h_tgt_sub->GetBinError(bin_err);

				h_sim_bg_int->SetBinContent(bin_itt, bin_con);
				h_sim_bg_int->SetBinError(bin_itt, bin_err);



				if (h_sim_omega_clone->GetBinContent(bin_itt) != 0.0 && h_tgt[in]->GetBinContent(bin_itt) != 0.0 ) { 
				
					bg_chi2 = bg_chi2 + ((fabs(bin_con) - fabs(h_sim_bg->GetBinContent(bin_itt)))/bin_err)**2;

//					cout << h_sim_omega_clone->GetBinContent(bin_itt) << "   Check chi2:  " << ((bin_con - h_sim_omega_clone->GetBinContent(bin_itt))/bin_err)**2 << endl;

					bg_dof++;

				}


			}

    		h_tgt_sub->Draw();

			h_sim_bg_int->SetLineColor(2);
			h_sim_bg_int->Draw("Esame");
			
			pp2->Update(); 
	

	
	//		h_sim_bg->SetFillColor(6);
			h_sim_bg->SetLineColor(6);
	
			h_sim_rho_clone->SetFillColor(4);
			h_sim_rho_clone->SetFillStyle(3004);
	
	
			h_sim_xphsp_clone->SetFillColor(3);
			h_sim_xphsp_clone->SetFillStyle(3005);
	
		
			h_sim_eta_clone->SetFillColor(1);
			h_sim_eta_clone->SetFillStyle(3006);
	

			h_sim_etap_clone->SetFillColor(1);
			h_sim_etap_clone->SetFillStyle(1001);
	

	
			h_sim_rho_clone->Draw("hsame");
			h_sim_xphsp_clone->Draw("hsame");
			h_sim_eta_clone->Draw("hsame");
			h_sim_etap_clone->Draw("hsame");
	
			h_sim_bg->Draw("hsame");
			
	// 		h_tgt_sub->Draw("same");
	
	
	
	
	// 		h_tgt[in]->Draw();
	
	 		boundary_line_integral.DrawLine(intg_lo, 0, intg_lo, h_tgt_sub->GetMaximum()*0.8);
	 		boundary_line_integral.DrawLine(intg_hi, 0, intg_hi, h_tgt_sub->GetMaximum()*0.8);

			bg_dof = bg_dof - 2;

 			TString chi2_str;
 			chi2_str.Form("#chi^{2}/dof: %.2f", bg_chi2/bg_dof);

			TLatex chi2_lx;
	 		chi2_lx.SetNDC();
	 		chi2_lx.SetTextSize(0.06);
			chi2_lx.DrawLatex(0.6, 0.7, chi2_str);	

			
	 		pp2->Update();
	//		delete h_tgt_sub;



			///*--------------------------------------------------*/
			/// Plot Omega from data and Sim

			TH1F* h_omega_sub = (TH1F*) h_tgt[in]->Clone("h_tgt_sub");
	
			h_omega_sub->Add(h_sim_rho_clone, -1);
			h_omega_sub->Add(h_sim_xphsp_clone, -1);
			h_omega_sub->Add(h_sim_eta_clone, -1);
			h_omega_sub->Add(h_sim_etap_clone, -1);
	
	  		pp7 = (TPad*) cc7->GetPad(i);
	 		pp7->cd();

			pp7->SetLeftMargin(0.1);
			pp7->SetRightMargin(0.002);







	
			h_omega_sub->GetXaxis()->SetRangeUser(0.60, 1);

			TH1F* h_omega_sub_int = (TH1F*) h_omega_sub->Clone("h_tgt_zero_int");
			h_omega_sub_int->Reset();


			Float_t omega_chi2 = 0.0;
			Float_t omega_dof = 0.0;

			for (Int_t bin_itt = lo_limit_bin; bin_itt <= hi_limit_bin; bin_itt++) {
	
				Double_t bin_con = h_omega_sub->GetBinContent(bin_itt);
	

				
				Float_t data_con  = h_tgt[in]->GetBinContent(bin_itt);

				Float_t data_err  = h_tgt[in]->GetBinError(bin_itt);

				if(h_tgt[in]->GetBinContent(bin_itt) == 0.0 && bin_con != 0.0) {
					data_err = 1.0;
				}

				Float_t omega_err = h_sim[in]->GetBinError(bin_itt) * omega_scale_fac;
				Float_t rho_err   = h_sim_rho[in]->GetBinError(bin_itt) * rho_scale_fac;
				Float_t xphsp_err = h_sim_xphsp[in]->GetBinError(bin_itt) * xphsp_scale_fac;
				Float_t eta_err = h_sim_eta[in]->GetBinError(bin_itt) * eta_scale_fac;
				Float_t etap_err = h_sim_etap[in]->GetBinError(bin_itt) * etap_scale_fac;

				Double_t bin_err = sqrt(data_err**2 + rho_err**2 + xphsp_err**2 + eta_err**2 + etap_err**2 +  omega_err**2);
//				Double_t bin_err = sqrt(data_err**2/data_con**2) * bin_err;

//				h_omega_sub->GetBinError(bin_err);

				
				
				h_omega_sub_int->SetBinContent(bin_itt, bin_con);
				h_omega_sub_int->SetBinError(bin_itt, bin_err);


				if (h_sim_omega_clone->GetBinContent(bin_itt) != 0.0 && h_tgt[in]->GetBinContent(bin_itt) != 0.0 ) { 
				
					omega_chi2 = omega_chi2 + ((fabs(bin_con) - fabs(h_sim_omega_clone->GetBinContent(bin_itt)))/bin_err)**2;

					cout << "Check chi2:  " << ((fabs(bin_con) - fabs(h_sim_omega_clone->GetBinContent(bin_itt)))/bin_err)**2 << "  " << fabs(bin_con) << "  " <<  bin_err  << "  " << data_err << "  " << rho_err << "  " << xphsp_err << endl;

					omega_dof++;

				}
			


			}



    		h_omega_sub->Draw();

			h_omega_sub_int->SetLineColor(2);
			h_omega_sub_int->Draw("Esame");

			
			pp7->Update(); 
	
			h_sim_omega_clone->SetFillColor(2);
			h_sim_omega_clone->SetFillStyle(3004);
	
			h_sim_omega_clone->Draw("hsame");
	
	 		boundary_line_integral.DrawLine(intg_lo, 0, intg_lo, h_tgt_sub->GetMaximum()*0.8);
	 		boundary_line_integral.DrawLine(intg_hi, 0, intg_hi, h_tgt_sub->GetMaximum()*0.8);

//			dof = hi_limit_bin - lo_limit_bin - 1;
			omega_dof = omega_dof - 1;

 			TString chi2_str;
 			chi2_str.Form("#chi^{2}/dof: %.2f", omega_chi2/omega_dof);
			
			TLatex chi2_lx;
	 		chi2_lx.SetNDC();
	 		chi2_lx.SetTextSize(0.06);
			chi2_lx.DrawLatex(0.6, 0.7, chi2_str);	


//			cout << "chi2/dof   " <<  omega_chi2/omega_dof << endl;








// 			if (in == 16) {
// 
// 				exit(0);
// 
// 			}

//			chi2_omega_acc->Fill(omega_chi2);

//			cout << "chi2/dof   " <<  omega_chi2/dof << endl;
			
//			exit(0);
	

			if (omega_chi2/omega_dof > 2.5) { 

				is_refit = true;

			} 


			if ( !is_refit || refit_itt!=0) { 

				chi2_bg_acc->Fill(bg_chi2/bg_dof);
				chi2_omega_acc->Fill(omega_chi2/omega_dof);

// 
// 				if (omega_chi2/omega_dof > 3.0)  {
// 					exit(0);
// 				}
// 

			}




	 		pp7->Update();






			///*--------------------------------------------------*/
			/// Plotting only the zero!

			pp6 = (TPad*) cc6->GetPad(i);
			pp6->cd();
			pp6->SetLeftMargin(0.1);
			pp6->SetRightMargin(0.002);


			TH1F* h_tgt_zero = (TH1F*) h_tgt[in]->Clone("h_tgt_sub");

			h_tgt_zero->Add(h_sim_omega_clone, -1);
			h_tgt_zero->Add(h_sim_rho_clone,   -1);
			h_tgt_zero->Add(h_sim_xphsp_clone, -1);
			h_tgt_zero->Add(h_sim_eta_clone,   -1);

//			h_tgt_zero->GetXaxis()->SetRangeUser(fit_local_lo-0.05, fit_local_hi+0.05);
			h_tgt_zero->GetXaxis()->SetRangeUser(0.6, 1.0);

			h_tgt_zero->Draw();

			Float_t bg_all = 0.0;

			TH1F* h_tgt_zero_int = (TH1F*) h_tgt_zero->Clone("h_tgt_zero_int");

			h_tgt_zero_int->Reset();

			Float_t omega_chi2 = 0.0;
			Float_t dof = 0.0;
			Float_t omega_yield_err_2 = 0.0;


			for (Int_t bin_itt = lo_limit_bin; bin_itt <= hi_limit_bin; bin_itt++) {
	
				Double_t bin_con = h_tgt_zero->GetBinContent(bin_itt);
//				Double_t bin_err = h_tgt_zero->GetBinError(bin_itt);
	
			    bg_all        = bg_all + bin_con;
//				xphsp_con_err_sq = xphsp_con_err_sq + bin_err**2;

				Float_t data_err  = h_tgt[in]->GetBinError(bin_itt);

				if(h_tgt[in]->GetBinContent(bin_itt) == 0.0 && bin_con != 0.0) {
//				if(h_tgt[in]->GetBinContent(bin_itt) == 0.0) {
					data_err = 1.0;
				}

//				data_err = 1.0;


				Float_t omega_err = h_sim[in]->GetBinError(bin_itt) * omega_scale_fac;
				Float_t rho_err   = h_sim_rho[in]->GetBinError(bin_itt) * rho_scale_fac;
				Float_t xphsp_err = h_sim_xphsp[in]->GetBinError(bin_itt) * xphsp_scale_fac;
				Float_t eta_err = h_sim_xphsp[in]->GetBinError(bin_itt) * xphsp_scale_fac;

				Double_t bin_err = sqrt(data_err**2 + omega_err**2 + rho_err**2 + xphsp_err**2 + eta_err**2);

//				h_tgt_zero->GetBinError(bin_err);
				
				h_tgt_zero_int->SetBinContent(bin_itt, bin_con);
				h_tgt_zero_int->SetBinError(bin_itt, bin_err);

				cout << "Check: " << h_tgt[in]->GetBinContent(bin_itt) << "  " <<  bin_err << "  " << data_err << "  " << omega_err << endl;  


// 				if (bin_err == 1 ) {
// 					omega_yield_err_2 = omega_yield_err_2 + bin_err**2;
// 				}




// 				if (h_sim_omega_clone->GetBinContent(bin_itt) != 0.0 && h_tgt[in]->GetBinContent(bin_itt) != 0.0 ) { 
// 				
// 					omega_chi2 = omega_chi2 + ((bin_con - h_tgt_zero->GetBinContent(bin_itt))/bin_err)**2;
// 					cout << h_sim_omega_clone->GetBinContent(bin_itt) << "   Check chi2:  " << ((bin_con - h_sim_omega_clone->GetBinContent(bin_itt))/bin_err)**2 << endl;
// 
// 					dof++;
//				}



			}

//			omega_yield_err = sqrt(omega_yield_err_2);


//			cout << omega_yield << "    " << omega_yield_err << endl;

//			exit(0);


			cout << "Total:  " << bg_all << "   " << h_tgt_zero->Integral(lo_limit_bin, hi_limit_bin) << endl;

//    			if (in == 18) {
//    				exit(0);
//    			}

			TF1* check =  new TF1("ch0", "[0]", intg_lo, intg_hi);

			h_tgt_zero_int->SetLineColor(2);
			h_tgt_zero_int->Draw("Esame");
			
			
			h_tgt_zero_int->Fit("ch0", "R");

//			check->Draw("same");

			Float_t check_zero = h_tgt_zero->Integral(lo_limit_bin, hi_limit_bin)/(hi_limit_bin-lo_limit_bin);

			TString zero_bg_str;
			zero_bg_str.Form("Int: %f", check_zero);
//			zero_bg_str.Form("Int: %f", bg_all);


//			cout << check->GetParameter(0) << endl;
//			exit(0);




 			TString fit_bg_str;
 			fit_bg_str.Form("Fit: %f#pm%f", check->GetParameter(0), check->GetParError(0));



			TLatex zero_bg;
	 		zero_bg.SetNDC();
	 		zero_bg.SetTextSize(0.06);
//	 		zero_bg.SetTextColor(6);
			zero_bg.DrawLatex(0.4, 0.25, zero_bg_str);	
			zero_bg.DrawLatex(0.4, 0.20, fit_bg_str);

	 		boundary_line_integral.DrawLine(intg_lo, 0, intg_lo, h_tgt_zero->GetMaximum()*0.8);
	 		boundary_line_integral.DrawLine(intg_hi, 0, intg_hi, h_tgt_zero->GetMaximum()*0.8);

			TLine *line = new TLine(intg_lo, check_zero, intg_hi, check_zero);
			line->SetLineColor(4);
			line->SetLineWidth(2);

			line->Draw("same");

			pp6->Update(); 



	 		TH1F* h_omega_temp = (TH1F*) h_tgt[in]->Clone("h_tgt_temp");
			
			raw_data_sum->Add(h_omega_temp, 1);
			data_omega_sum->Add(h_omega_sub, 1);
			sim_omega_sum->Add(h_sim_omega_clone, 1);

			sim_rho_sum->Add(h_sim_rho_clone, 1);
			sim_xphsp_sum->Add(h_sim_xphsp_clone, 1);
			sim_eta_sum->Add(h_sim_eta_clone, 1);
			sim_etap_sum->Add(h_sim_etap_clone, 1);




			////*--------------------------------------------------*/
			// Safe Guard Check statement

			if (fabs(omega_yield_err) > 100*fabs(omega_yield) ) {

				cout << "The Error is larger than the scale factor:" << endl;


				cout << omega_yield_err << "  "<< sqrt(data_con_err_sq/pow(data_con,2)) << "  " << omega_scale_fac_err/omega_scale_fac << endl;

//				cout << omega_yield_err << "  "<< sqrt(data_con_err_sq/pow(data_con,2)) << "  " << sqrt(rho_con_err_sq/pow(rho_con,2)) << "  " << sqrt(xphsp_con_err_sq/pow(xphsp_con,2)) << "  " << "  " << omega_scale_fac_err/omega_scale_fac << endl;

//				exit(0);

			}










		} else if (low_stat) {

			TText ext_txt;
	 		ext_txt.SetNDC();
	 		ext_txt.SetTextSize(0.1);
			ext_txt.SetTextAngle(45);
			ext_txt.SetTextColor(4);
			ext_txt.DrawText(0.15, 0.15, "Excluded Due to Statistics");	

			is_excluded = true;

		} else if (rad_only) {

			TText ext_txt;
	 		ext_txt.SetNDC();
	 		ext_txt.SetTextSize(0.1);
			ext_txt.SetTextAngle(45);
			ext_txt.SetTextColor(2);
			ext_txt.DrawText(0.15, 0.15, "Excluded Due to Radiative Tail");	

			is_excluded = true;

		} else if (off_focal) {

			TLatex ext_txt_1;
	 		ext_txt_1.SetNDC();
	 		ext_txt_1.SetTextSize(0.1);
			ext_txt_1.SetTextAngle(45);
			ext_txt_1.SetTextColor(6);
			ext_txt_1.DrawLatex(0.15, 0.15, "#omega is off the focal plane");	

			is_excluded = true;

		}






		if (!is_refit || refit_itt!=0) { 


			if (omega_yield < 0.0) {
	
				omega_yield = 0.0;
	
			}
	



			fit_out_file << fixed << setprecision(4) << Int_t((i-1)/phi_bin_num) + 1 << "  " << (i-1)%phi_bin_num + 1
	           << " "  << setw(9) <<  omega_scale_fac << "  " << setw(9) <<  rho_scale_fac 
	           << "  " << setw(9) <<  xphsp_scale_fac << "  " << setw(9) <<  eta_scale_fac  
	           << "  " << setw(9) <<  etap_scale_fac << "  "  << intg_lo << "  " << intg_hi 
			   << "  " << setw(7) <<  omega_yield << " " << setw(7) << omega_yield_err << endl;




//			if (omega_yield_err == 0 && in == 2) {
// 			if (in == 1) {
// 
// 				cout << fit_local_lo << "  " << fit_local_hi << endl; 
// 				cout << "Omega yield error: " << omega_yield_err << endl;
// 
// 				exit(0);
// 				
// 			}




	
			cc4->Update();
	
			TCanvas* cc5 = new TCanvas("cc5", "cc5", 1600, 400);
			cc5->Divide(4, 1);
	
			cc5->cd(1);
	
	 		ppnew = (TPad*) cur_pad->Clone();
	
	//		ppnew->SetPad(0.01, 0.01, 0.01, 0.01);
	//		TPad *pnew = new TPad();
	
			ppnew->SetPad(cc5->GetX1(),cc5->GetY1(),cc5->GetX2(),cc5->GetY2());
	
	//		cout << out_dir_str + "individual_plot/" + f_str_temp + ".png" << endl;
			
			ppnew->Draw();
	
	
			if (!is_excluded) {		
	
	 			cc5->cd(2);
	 			ppnew1  = (TPad*) pp7->Clone();
			    ppnew1->SetPad(cc5->GetX1(),cc5->GetY1(),cc5->GetX2(),cc5->GetY2());
	 			ppnew1->Draw();
	
	 			cc5->cd(3);
	 			ppnew2  = (TPad*) pp2->Clone();
			    ppnew2->SetPad(cc5->GetX1(),cc5->GetY1(),cc5->GetX2(),cc5->GetY2());
	 			ppnew2->Draw();
	
	 			cc5->cd(4);
	 			ppnew3  = (TPad*) pp6->Clone();
			    ppnew3->SetPad(cc5->GetX1(),cc5->GetY1(),cc5->GetX2(),cc5->GetY2());
	 			ppnew3->Draw();
	
			}
	
			cc5->Print(out_dir_str + "individual_plot/" + f_str_temp + ".png");
	
			cc5->Write(f_str_temp);

		}






//		exit(0);
		
//		delete pnew;

// 		TH1F* h_tgt[in] = (TH1F*) h_tgt_tmp->Clone("h_tgt_clone");
//		TH1F* h_tgt_sub = (TH1F*) h_tgt[in]->Clone("h_tgt_sub");

// 		h_tgt_sub->Add(h_sim_rho_clone,-1);
// 		h_tgt_sub->Add(h_sim_xphsp_clone,-1);
// 		h_tgt_sub->Add(h_sim_eta_clone,-1);
// 		h_tgt_sub->Add(h_sim_etap_clone,-1);

//  		pp2 = (TPad*) cc2->GetPad(i);
// 		pp2->cd();
// 		h_tgt_sub->Draw();
// 
// 		boundary_line_fitting.DrawLine(fit_local_lo, 0, fit_local_lo, h_tgt[in]->GetMaximum()*0.8);
// 		boundary_line_fitting.DrawLine(fit_local_hi, 0, fit_local_hi, h_tgt[in]->GetMaximum()*0.8);
// 

//		delete h_tgt_sub;

//		cc4->cd();
 

	

		if (is_refit) {

			i = i - 1;

		}



		delete cc5;

	}

	cc4->Print(out_dir_str + "u_fit_phi_" + file_name+".png");
	cc2->Print(out_dir_str + "u_fit_phi_sub_" + file_name+".png");
	cc6->Print(out_dir_str + "u_fit_phi_zero_" + file_name+".png");
	cc7->Print(out_dir_str + "u_fit_phi_omega_" + file_name+".png");


	cc4->Write(out_dir_str + "u_fit_phi_" + file_name+".root");
	cc2->Write(out_dir_str + "u_fit_phi_sub_" + file_name+".root");
	cc6->Write(out_dir_str + "u_fit_phi_zero_" + file_name+".root");
	cc7->Write(out_dir_str + "u_fit_phi_omega_" + file_name+".png");


	TCanvas* chi2_can = new TCanvas("chi2_can", "Chi2_can", 800, 400);
	chi2_can->Divide(2,1);
	 
	chi2_can->cd(1);
	chi2_omega_acc->Draw();

	chi2_can->cd(2);
	chi2_bg_acc->Draw();
	chi2_can->Update();

	chi2_can->Print(out_dir_str + "chi2_check_" + file_name+".png");

	TFile* chi2_file = new TFile(out_dir_str + "chi2_check_" + file_name+".root", "recreate");

	chi2_omega_acc->Write("chi2_omega");	
	chi2_bg_acc->Write("chi2_bg");	

	chi2_file->Close();
	delete chi2_file;



	TCanvas* omega_check_can = new TCanvas();


	sim_omega_sum->SetLineColor(2);
	
 	data_omega_sum->Draw("E");
 	sim_omega_sum->Draw("hsame");

	omega_check_can->Print(out_dir_str + "Omega_" + file_name+".png");

	

	omega_check_can->Clear();

 	raw_data_sum->Draw("E");

 	sim_omega_sum->Draw("hsame");
  	sim_rho_sum->Draw("hsame");
  	sim_xphsp_sum->Draw("hsame");
  	sim_eta_sum->Draw("hsame");
  	sim_etap_sum->Draw("hsame");

	

	TH1F* sim_sum = (TH1F*)sim_omega_sum->Clone("sim_sum");
 	sim_sum->Add(sim_rho_sum, 1);
 	sim_sum->Add(sim_xphsp_sum, 1);
 	sim_sum->Add(sim_eta_sum, 1);
 	sim_sum->Add(sim_etap_sum, 1);

	sim_sum->Draw("hsame");


	omega_check_can->Print(out_dir_str + "data_check_" + file_name+".png");

//	omega_check_can->Write(file_name + "_presentation");




	delete omega_check_can;

	delete data_omega_sum;
	delete sim_omega_sum;

//	exit(0);


//	cc4->Print(out_dir_str + "u_fit_phi_" + file_name+".root");
//	cc2->Print(out_dir_str + "u_fit_phi_sub_" + file_name+".root");

//	exit(0);

	delete chi2_omega_acc;
	delete chi2_bg_acc;

	delete chi2_can;

	delete pp2;
	delete pp6;
	delete pp7;

	delete cc2;
	delete cc6;
	delete cc7;

	delete cc4;
	delete file_out_root;

//	exit(0);

}











/*--------------------------------------------------*/
/*--------------------------------------------------*/
/*--------------------------------------------------*/
///
/// Date: 26/April/2016
///
/// WL
/// Itteration process to integrate Omega simulation 
///

void u_phi_int_sim_itt(TCanvas* cc3, TString file_name) {

	const Int_t u_bin_num_const   = u_bin_num;  
	const Int_t phi_bin_num_const = phi_bin_num;  
	const Int_t bin_num_const     = u_bin_num_const * phi_bin_num_const;

//	exit(0);

    TString out_int_dir_str;

	out_int_dir_str = "u_int_fit/";

	is_itt = true;

	cout << "asddddddddddddddddd " << endl;

	ifstream file1;

	TString file_str;
	
	file_str = "u_bin_fit/fit_output.pl_" + file_name + ".dat";

	file1.open(file_str);
	
	

	cout << file_str << endl;
	
 
	Int_t   t_bin, phi_bin;
	Float_t  omega_scale_fac, rho_scale_fac, xphsp_scale_fac;
	Float_t  eta_scale_fac, etap_scale_fac;
	Float_t  lo_limit, hi_limit, exp_yield, exp_yield_err;

	vector<Float_t> lo_limit_vec, hi_limit_vec, exp_yield_vec, exp_yield_err_vec;
	vector<Float_t> rho_f_vec, xphsp_f_vec, eta_f_vec, etap_f_vec;

 	string str;

 	if (file1.is_open()) {
// 		while (!file1.eof()) {
 		while (	getline(file1, str)) {

//     		file1 >> t_bin >> phi_bin >> omega_scale_fac >> rho_scale_fac >> xphsp_scale_fac
// 				  >> lo_limit >> hi_limit >> exp_yield; 
//     		cout<< lo_limit << "  " << hi_limit << "  " << exp_yield << endl;
// 			cout << str << endl;

			stringstream ss(str);

     		ss >> t_bin >> phi_bin >> omega_scale_fac >> rho_scale_fac >> xphsp_scale_fac 
               >> eta_scale_fac >> etap_scale_fac >> lo_limit >> hi_limit >> exp_yield >> exp_yield_err;

			lo_limit_vec.push_back(lo_limit);
			hi_limit_vec.push_back(hi_limit);
			exp_yield_vec.push_back(exp_yield);
			exp_yield_err_vec.push_back(exp_yield_err);

			rho_f_vec.push_back(rho_scale_fac);
			xphsp_f_vec.push_back(xphsp_scale_fac);
			eta_f_vec.push_back(eta_scale_fac);
			etap_f_vec.push_back(etap_scale_fac);

 		}

	}

//	cout << exp_yield_vec.size() << endl;
//	exit(0);

	


	TCanvas* cc4 = new TCanvas("cc4", "cc4", 1600, 800);
	cc4->Divide(phi_bin_num_const, u_bin_num_const, 0.003, 0.003);
 
 	TH1F* h_sim       [bin_num_const];
 	TH1F* h_tgt       [bin_num_const];
 	TH1F* h_sim_rho   [bin_num_const];
 	TH1F* h_sim_xphsp [bin_num_const];
	TH1F* h_sim_eta   [bin_num_const];
 	TH1F* h_sim_etap  [bin_num_const];

	for (Int_t i = 1; i <= u_bin_num*phi_bin_num; i++) {

		TString f_str_temp;
		f_str_temp.Form("_%i_%i", Int_t((i-1)/phi_bin_num) + 1, (i-1)%phi_bin_num + 1);
		f_str_temp = "u_fit_phi_" + file_name + f_str_temp;


//		i = 11;
	
		TPad* cur_pad = (TPad*) cc4->GetPad(i);

		cur_pad->cd();
		cur_pad->SetLeftMargin(0.1);
		cur_pad->SetRightMargin(0.002);

//		cur_pad->SetLogy();


//		cout << i << "   " << i-1 << endl;

		const int in =  i -1;

//		cout << in << endl;

  		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("u_phi_sim_omega");
  		TH1F* h_sim[in] = (TH1F*) h_sim_tmp->Clone("h_sim_clone");

		TH1F* h_tgt_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("u_phi_target");
		TH1F* h_tgt[in] = (TH1F*) h_tgt_tmp->Clone("h_tgt_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("u_phi_sim_rho");
 		TH1F* h_sim_rho[in] = (TH1F*) h_sim_tmp->Clone("h_sim_rho_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("u_phi_sim_xphsp");
 		TH1F* h_sim_xphsp[in] = (TH1F*) h_sim_tmp->Clone("h_sim_xphsp_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("u_phi_sim_eta");
 		TH1F* h_sim_eta[in] = (TH1F*) h_sim_tmp->Clone("h_sim_eta_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("u_phi_sim_etap");
 		TH1F* h_sim_etap[in] = (TH1F*) h_sim_tmp->Clone("h_sim_etap_clone");
// 
// 
//
//
//

// 		TH1F* h_omega_sub = (TH1F*) h_tgt[in]->Clone("h_tgt_sub");
// 	
// 		h_omega_sub->Add(h_sim_rho_clone, -1);
// 		h_omega_sub->Add(h_sim_xphsp_clone, -1);
// 		h_omega_sub->Add(h_sim_eta_clone, -1);
// 		h_omega_sub->Add(h_sim_etap_clone, -1);
// 	
//  		pp7->cd();
// 	
// 		h_omega_sub->GetXaxis()->SetRangeUser(0.60, 1);

		
		h_sim_rho[in]->Scale(rho_f_vec[in]);
        h_sim_xphsp[in]->Scale(xphsp_f_vec[in]);
        h_sim_eta[in]->Scale(eta_f_vec[in]);
        h_sim_etap[in]->Scale(etap_f_vec[in]);

		h_tgt[in]->Add(h_sim_rho[in], -1);
		h_tgt[in]->Add(h_sim_xphsp[in], -1);
// 		h_tgt[in]->Add(h_sim_eta[in], -1);
// 		h_tgt[in]->Add(h_sim_etap[in], -1);

 		if (eta_f_vec[in] != 0.0) {
 			h_tgt[in]->Add(h_sim_eta[in], -1);
 		} 
 

//		cout << "!!!  asdasd   " << in << "   " << rho_f_vec[in] << "  " << xphsp_f_vec[in] << "  " << eta_f_vec[in] <<  "  " << etap_f_vec[in] << endl;

 		if (etap_f_vec[in] != Float_t(0)) {
 			h_tgt[in]->Add(h_sim_etap[in], -1);
//			cout << "???  asdasd   " << in <<  "    " << etap_f_vec[in] << endl;
 		}
 








		

		h_tgt[in]->GetXaxis()->SetRangeUser(0.60, 1);


		Double_t fit_local_lo;
		Double_t fit_local_hi;

// 		fit_local_lo = 0.5;
// 		fit_local_hi = 0.9;

		TPad* tar_pad = (TPad*) cc3->GetPad(i);


		TF1 *fit_u_fun = new TF1("fit_u_fun", fun_sim_total, fit_local_lo, fit_local_hi, 3);

 		fit_u_fun->SetParameters(3, 40, 5000);

		fit_u_fun->SetParLimits(0, 0, 100000);
		fit_u_fun->SetParLimits(1, 0, 100000);
		fit_u_fun->SetParLimits(2, 0, 100000);

 
//		h_tgt[in]->Fit("fit_u_fun", "MNR");

		h_sim[in]->SetStats(0);

		h_sim[in]->SetFillColor(2);

 		h_tgt[in]->SetTitle(f_str_temp);
 
 		if(h_tgt[in]->GetMaximum() < h_sim[in]->GetMaximum()) {
 			h_tgt[in]->SetMaximum(h_sim[in]->GetMaximum() * 1.2);
 		}

 		h_tgt[in]->DrawCopy();
		h_sim[in]->DrawCopy("hsame");
		h_tgt[in]->DrawCopy("same");



//		h_tgt[in]->SetMaximum




		
// 		/*--------------------------------------------------*/
// 		/// Data output
// 		//
// 
// 		TString range_str;
// 		range_str.Form("Range: %f - %f ", fit_local_lo, fit_local_hi);
// 
// 		TString mean_str;
// 		mean_str.Form("Mean: %f RMS: %f Sim: %f ", h_tgt[in]->GetMean(1), h_tgt[in]->GetRMS(1), h_sim[in]->GetBinCenter(h_sim[in]->GetMaximumBin()));
// 
// 		TString ratio_str;
// 		ratio_str.Form("w_f: %f, r_f: %f, x_f: %f", fit_u_fun->GetParameter(0), fit_u_fun->GetParameter(1), fit_u_fun->GetParameter(2));
// 

// 
// 
//  		TText fitting_range;
// 
//  		fitting_range.SetNDC();
//  		fitting_range.SetTextSize(0.04);
// 
//  		fitting_range.DrawText(0.12, 0.85, range_str);
//  		fitting_range.DrawText(0.12, 0.80, mean_str);
//  		fitting_range.DrawText(0.12, 0.75, ratio_str);


		Int_t lo_limit_bin;
		Int_t hi_limit_bin;

		Double_t intg_lo = lo_limit_vec[in];
		Double_t intg_hi = hi_limit_vec[in];


		if (intg_lo != 0.0 && intg_hi != 0.0) {
	
			Double_t omega_yield     = exp_yield_vec[in];
			Double_t omega_yield_err = exp_yield_err_vec[in];

//			lo_limit_bin = h_sim[in]->FindBin(intg_lo + 0.014);
//        	hi_limit_bin = h_sim[in]->FindBin(intg_hi - 0.014);


// 			lo_limit_bin = h_sim[in]->FindBin(intg_lo + 0.02);
//         	hi_limit_bin = h_sim[in]->FindBin(intg_hi - 0.02);

 			lo_limit_bin = h_sim[in]->FindBin(intg_lo);
         	hi_limit_bin = h_sim[in]->FindBin(intg_hi);


			Double_t content_total = 0.0;
			Double_t content_total_error_sq = 0.0;


			
	    	bool low_stat = false;


			for (Int_t bin_itt = lo_limit_bin; bin_itt <= hi_limit_bin; bin_itt++) {

				Double_t bin_err = h_sim[in]->GetBinError(bin_itt);
//				cout << "Ased  " << bin_itt  << "  " << bin_con << "  " << bin_con << "  " << bin_err/bin_con << endl;

				content_total_error_sq = content_total_error_sq + bin_err**2;


			}


//			Double_t sim_omega_yield = h_sim_omega_clone->Integral(lo_limit_bin, hi_limit_bin);

			Double_t sim_omega_yield = h_sim[in]->Integral(lo_limit_bin, hi_limit_bin);

//			cout << sim_omega_yield << endl;





// /*--------------------------------------------------*/ 
// //  Low Statistics exclusion creteria is removed 
// //  Date: 20/Sep/2017
//
// 			if (omega_yield < low_stats_limit) { 

 			if (omega_yield < low_stats_limit) { 

//				cout << omega_yield << "    " << low_stats_limit << endl;
//				low_stat = true;

//				exit(0);

			}




			/*--------------------------------------------------*/
			// Oct 26, 2016
			// Garth and Bill is introducing a 10% integration limit systematic error and a 5% fitting limit systematic error
			// This is ony temp, will be revisited.

 			if (!low_stat) { 

				acc_yield[in]     = acc_yield[in] + omega_yield ;
	
				acc_sim_yield[in] = acc_sim_yield[in] + sim_omega_yield ;

// 			acc_yield_err[in] =  acc_yield_err[in] + omega_yield_err**2 + (omega_yield*0.05)**2 + (omega_yield*0.1)**2; 

 				acc_yield_err[in] =  acc_yield_err[in] + omega_yield_err**2; 
	
				acc_sim_yield_err[in] = acc_sim_yield_err[in] + content_total_error_sq;

			}

 

			/*--------------------------------------------------*/
			/// Draw fitting boundaries

			TLine boundary_line_integral;
			boundary_line_integral.SetLineWidth(2);
			boundary_line_integral.SetLineStyle(7);
			boundary_line_integral.SetLineColor(4);
			boundary_line_integral.DrawLine(intg_lo, 0, intg_lo, h_tgt[in]->GetMaximum()*0.8);
			boundary_line_integral.DrawLine(intg_hi, 0, intg_hi, h_tgt[in]->GetMaximum()*0.8);


			TString sim_yield_str;
			sim_yield_str.Form("Sim Yield: %f", sim_omega_yield);
			
			TString exp_yield_str;
			exp_yield_str.Form("Exp Yield: %f", omega_yield);
			


 			TText fitting_range;

 			fitting_range.SetNDC();
 			fitting_range.SetTextSize(0.07);

 			fitting_range.DrawText(0.2, 0.7, sim_yield_str);
 			fitting_range.DrawText(0.2, 0.65, exp_yield_str);





 			if (low_stat) { 

 				TText excluded_txt_1;

 				excluded_txt_1.SetNDC();
 				excluded_txt_1.SetTextSize(0.2);
 				excluded_txt_1.SetTextColor(2);
 				excluded_txt_1.SetTextAngle(-45);

	 			excluded_txt_1.DrawText(0.2, 0.8, "Excluded !");
		
//				exit(0);

			}



			

// 			exit(0);
//


			cc4->Update();


			TCanvas* cc5 = new TCanvas("cc5", "cc5", 800, 600);
			cc5->cd();
//			pnew->SetPad(0.01, 0.01, 0.01, 0.01);

			TPad *pnew = (TPad*) cur_pad->Clone();

			pnew->SetPad(cc5->GetX1(),cc5->GetY1(),cc5->GetX2(),cc5->GetY2());

			cout << out_int_dir_str + "individual_plot/" + f_str_temp + ".png" << endl;
			
			pnew->Draw();





			cc5->Print(out_int_dir_str + "individual_plot/" + f_str_temp + ".png");
//			cc5->Print(out_int_dir_str + "individual_plot/" + f_str_temp + ".C");

//			exit(0);
			
			delete cc5;




		} else {

			h_tgt[in]->Reset();
			h_tgt[in]->Draw();

//			cur_pad->Clear(); 

 			TText excluded_txt;

 			excluded_txt.SetNDC();
 			excluded_txt.SetTextSize(0.2);
 			excluded_txt.SetTextAngle(45);

 			excluded_txt.DrawText(0.2, 0.2, "Excluded !");

		}

		cc4->cd();

	}

	cc4->Print(out_int_dir_str + "u_fit_itt_test_" + file_name+".png");

//	cout << "asdasdasdasdas  asdasd <<<  aEnd Loop" << endl;


	delete cc4;

//	exit(0);
 
}


































/*--------------------------------------------------*/









void u_phi_fit_sim(TCanvas* cc3, TString file_name, Double_t fit_low_1, Double_t fit_low_2, Double_t fit_mid_1, Double_t fit_mid_2, Double_t fit_high_1, Double_t fit_high_2) {
	
	cout << "asddddddddddddddddd " << endl;


//	cc3->Draw();

	const Int_t u_bin_num_const   = u_bin_num;  
	const Int_t phi_bin_num_const = phi_bin_num;  

	TCanvas* cc4 = new TCanvas("cc4", "cc4", 1600, 800);
	cc4->Divide(phi_bin_num_const, u_bin_num_const, 0.003, 0.003);

	TH1F* h_tgt[6];
	TH1F* h_sim[6];
	TH1F* h_sim_rho[6];
	TH1F* h_sim_xphsp[6];

	for (Int_t i = 1; i <= phi_bin_num; i++) {
	
		TPad* cur_pad = cc4->GetPad(i);

		cur_pad->cd();
		cur_pad->SetLeftMargin(0.1);
		cur_pad->SetRightMargin(0.002);



		cout << i << "   " << i-1 << endl;

		const int in =  i -1;

		cout << in << endl;

//		TH1F* h_tgt[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
//		TH1F* h_sim[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");

 		TH1F* h_tgt_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
 		TH1F* h_tgt[in] = (TH1F*) h_tgt_tmp->Clone("h_tgt_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");
 		TH1F* h_sim[in] = (TH1F*) h_sim_tmp->Clone("h_sim_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim_rho");
 		TH1F* h_sim_rho[in] = (TH1F*) h_sim_tmp->Clone("h_sim_rho_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim_xphsp");
 		TH1F* h_sim_xphsp[in] = (TH1F*) h_sim_tmp->Clone("h_sim_xphsp_clone");


	 	h_tgt[in]->DrawCopy();  
	 	h_sim[in]->DrawCopy("same");

		peak_omega = h_sim[in];
		peak_rho   = h_sim_rho[in];
		peak_xphsp = h_sim_xphsp[in];


		Double_t fit_local_lo;
		Double_t fit_local_hi;



  		fit_local_lo = fit_low_1;
  		fit_local_hi = fit_low_2;

		TF1 *fit_u_fun = new TF1("fit_u_fun", fun_sim_total, fit_local_lo, fit_local_hi, 3);

 		fit_u_fun->SetParameters(3, 40, 5000);

		fit_u_fun->SetParLimits(0, 0, 100000);
		fit_u_fun->SetParLimits(1, 0, 100000);
		fit_u_fun->SetParLimits(2, 0, 100000);

 
		h_tgt[in]->Fit("fit_u_fun", "MNR");



	 	TH1F *h_sim_omega_clone = (TH1F*) h_sim[in]->Clone();
	 	TH1F *h_sim_rho_clone   = (TH1F*) h_sim_rho[in]->Clone();
	 	TH1F *h_sim_xphsp_clone = (TH1F*) h_sim_xphsp[in]->Clone();


		h_sim_omega_clone->Scale(fit_u_fun->GetParameter(0));		
		h_sim_rho_clone->Scale(fit_u_fun->GetParameter(1));		
		h_sim_xphsp_clone->Scale(fit_u_fun->GetParameter(2));		

 	 	h_sim_omega_clone->DrawCopy("same");
 	 	h_sim_rho_clone->DrawCopy("same");
 	 	h_sim_xphsp_clone->DrawCopy("same");


 		h_sim_xphsp_clone->Add(h_sim_omega_clone);		
 		h_sim_xphsp_clone->Add(h_sim_rho_clone);		


		h_sim_xphsp_clone->SetLineColor(6);
 	 	h_sim_xphsp_clone->DrawCopy("same");



	}




 	for (Int_t i =  phi_bin_num + 1; i <= 2*phi_bin_num; i++) {
	
		TPad* cur_pad1 = cc4->GetPad(i);

		cur_pad1->cd();
		cur_pad1->SetLeftMargin(0.1);
		cur_pad1->SetRightMargin(0.002);



		cout << i << "   " << i-1 << endl;

		const int in =  i - phi_bin_num;

		cout << in << endl;

		TH1F* h_tgt_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
 		TH1F* h_tgt[in] = (TH1F*) h_tgt_tmp->Clone("h_tgt_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");
 		TH1F* h_sim[in] = (TH1F*) h_sim_tmp->Clone("h_sim_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim_rho");
 		TH1F* h_sim_rho[in] = (TH1F*) h_sim_tmp->Clone("h_sim_rho_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim_xphsp");
 		TH1F* h_sim_xphsp[in] = (TH1F*) h_sim_tmp->Clone("h_sim_xphsp_clone");


	 	h_tgt[in]->DrawCopy();  
	 	h_sim[in]->DrawCopy("same");

		peak_omega = h_sim[in];
		peak_rho   = h_sim_rho[in];
		peak_xphsp = h_sim_xphsp[in];


		Double_t fit_local_lo;
		Double_t fit_local_hi;


		fit_local_lo = fit_mid_1;
		fit_local_hi = fit_mid_2;

//		exit(0);

		TF1 *fit_u_fun = new TF1("fit_u_fun", fun_sim_total, fit_local_lo, fit_local_hi, 3);

 		fit_u_fun->SetParameters(3, 40, 500);

		fit_u_fun->SetParLimits(0, 0, 100000);
		fit_u_fun->SetParLimits(1, 0, 100000);
		fit_u_fun->SetParLimits(2, 0, 100000);

 
		h_tgt[in]->Fit("fit_u_fun", "MNR");


	 	TH1F *h_sim_omega_clone = (TH1F*) h_sim[in]->Clone();
	 	TH1F *h_sim_rho_clone   = (TH1F*) h_sim_rho[in]->Clone();
	 	TH1F *h_sim_xphsp_clone = (TH1F*) h_sim_xphsp[in]->Clone();


		h_sim_omega_clone->Scale(fit_u_fun->GetParameter(0));		
		h_sim_rho_clone->Scale(fit_u_fun->GetParameter(1));		
		h_sim_xphsp_clone->Scale(fit_u_fun->GetParameter(2));		

 	 	h_sim_omega_clone->DrawCopy("same");
 	 	h_sim_rho_clone->DrawCopy("same");
 	 	h_sim_xphsp_clone->DrawCopy("same");


 		h_sim_xphsp_clone->Add(h_sim_omega_clone);		
 		h_sim_xphsp_clone->Add(h_sim_rho_clone);		


		h_sim_xphsp_clone->SetLineColor(6);
 	 	h_sim_xphsp_clone->DrawCopy("same");



	}


 	for (Int_t i =  2*phi_bin_num + 1; i <= 3*phi_bin_num; i++) {

		TPad* cur_pad2 = cc4->GetPad(i);

		cur_pad2->cd();
		cur_pad2->SetLeftMargin(0.1);
		cur_pad2->SetRightMargin(0.002);



		cout << i << "   " << i-1 << endl;

		const int in =  i - (2*phi_bin_num + 1);
	

		cout << in << endl;

		TH1F* h_tgt_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
 		TH1F* h_tgt[in] = (TH1F*) h_tgt_tmp->Clone("h_tgt_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");
 		TH1F* h_sim[in] = (TH1F*) h_sim_tmp->Clone("h_sim_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim_rho");
 		TH1F* h_sim_rho[in] = (TH1F*) h_sim_tmp->Clone("h_sim_rho_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim_xphsp");
 		TH1F* h_sim_xphsp[in] = (TH1F*) h_sim_tmp->Clone("h_sim_xphsp_clone");


	 	h_tgt[in]->DrawCopy();  
	 	h_sim[in]->DrawCopy("same");

		peak_omega = h_sim[in];
		peak_rho   = h_sim_rho[in];
		peak_xphsp = h_sim_xphsp[in];


		Double_t fit_local_lo;
		Double_t fit_local_hi;


		fit_local_lo = fit_high_1;
		fit_local_hi = fit_high_2;

//		exit(0);

		TF1 *fit_u_fun = new TF1("fit_u_fun", fun_sim_total, fit_local_lo, fit_local_hi, 3);

 		fit_u_fun->SetParameters(3, 40, 500);

		fit_u_fun->SetParLimits(0, 0, 100000);
		fit_u_fun->SetParLimits(1, 0, 100000);
		fit_u_fun->SetParLimits(2, 0, 100000);

 
		h_tgt[in]->Fit("fit_u_fun", "MNR");


	 	TH1F *h_sim_omega_clone = (TH1F*) h_sim[in]->Clone();
	 	TH1F *h_sim_rho_clone   = (TH1F*) h_sim_rho[in]->Clone();
	 	TH1F *h_sim_xphsp_clone = (TH1F*) h_sim_xphsp[in]->Clone();


		h_sim_omega_clone->Scale(fit_u_fun->GetParameter(0));		
		h_sim_rho_clone->Scale(fit_u_fun->GetParameter(1));		
		h_sim_xphsp_clone->Scale(fit_u_fun->GetParameter(2));		

 	 	h_sim_omega_clone->DrawCopy("same");
 	 	h_sim_rho_clone->DrawCopy("same");
 	 	h_sim_xphsp_clone->DrawCopy("same");


 		h_sim_xphsp_clone->Add(h_sim_omega_clone);		
 		h_sim_xphsp_clone->Add(h_sim_rho_clone);		


		h_sim_xphsp_clone->SetLineColor(6);
 	 	h_sim_xphsp_clone->DrawCopy("same");



	}


	cc4->Print(out_dir_str + "u_fit_phi_no_const_bg" + file_name+".png");

	delete cc4;

//	exit(0);

}

































/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Invividual fitting with 3rd order poly

void u_phi_fit(TCanvas* cc3, TString file_name, Double_t fit_low_1, Double_t fit_low_2, Double_t fit_mid_1, Double_t fit_mid_2, Double_t fit_high_1, Double_t fit_high_2) {
	
	cout << "asddddddddddddddddd " << endl;

//	cc3->Draw();

	TCanvas* cc4 = new TCanvas("cc4", "cc4", 1600, 800);
	cc4->Divide(phi_bin_num_const, u_bin_num_const, 0.003, 0.003);

	TH1F* h_tgt[6];
	TH1F* h_sim[6];

	for (Int_t i = 1; i <= phi_bin_num; i++) {
	
		TPad* cur_pad = cc4->GetPad(i);

		cur_pad->cd();
		cur_pad->SetLeftMargin(0.1);
		cur_pad->SetRightMargin(0.002);



		cout << i << "   " << i-1 << endl;

		const int in =  i -1;

		cout << in << endl;

//		TH1F* h_tgt[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
//		TH1F* h_sim[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");

 		TH1F* h_tgt_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
 		TH1F* h_tgt[in] = (TH1F*) h_tgt_tmp->Clone("h_tgt_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");
 		TH1F* h_sim[in] = (TH1F*) h_sim_tmp->Clone("h_sim_clone");


	 	h_tgt[in]->DrawCopy();  
	 	h_sim[in]->DrawCopy("same");

		peak1 = h_sim[in];

//		Double_t fit_local_lo;
//		Double_t fit_local_hi;

//  		fit_local_lo = fit_low_1;
//  		fit_local_hi = fit_low_2;



		TPad* tar_pad = (TPad*) cc3->GetPad(i);

  		fit_local_lo = Find_low_range(tar_pad);
  		fit_local_hi = Find_high_range(tar_pad);


		TF1 *fit_u_fun = new TF1("fit_u_fun", ftotal_all, fit_local_lo, fit_local_hi, 5);
		TF1 *fit_u_fun_pol = new TF1("fit_u_fun_pol", ftotal_br, fit_local_lo, fit_local_hi, 4);


// 		TF1 *fit_u_fun = new TF1("fit_u_fun", ftotal_all_skewed_gauss, fit_local_lo, fit_local_hi, 5);
// 		TF1 *fit_u_fun_pol = new TF1("fit_u_fun_pol", skewed_gaussion_fun, fit_local_lo, fit_local_hi, 3);
// 
// 		fit_u_fun->SetParameters(0.001, 0.7, 0.03, 1);
// 
//		TF1 *fit_u_fun = new TF1("fit_u_fun", "gaus", fit_local_lo, fit_local_hi);
//		TF1 *fit_u_fun_pol = new TF1("fit_u_fun_pol", skewed_gaussion_fun, fit_local_lo, fit_local_hi, 3);


//		fit_u_fun->SetParameters(10.0, -50.0, 75, -36.0, 5);
// 	  	fit_u_fun->SetParLimits(4, 0, 60);




// // 		TF1 *fit_u_fun = new TF1("fit_u_fun", "TMath::Dilog", fit_local_lo, fit_local_hi);
//  		TF1 *fit_u_fun = new TF1("fit_u_fun", "", fit_local_lo, fit_local_hi);
// 
		h_tgt[in]->Fit("fit_u_fun", "MEURN");

	 	h_sim[in]->Draw("same");

	 	TH1F *h_sim_clone = (TH1F*) h_sim[in]->Clone();
   	

		fit_u_fun_pol->SetLineColor(4);

 		fit_u_fun_pol->FixParameter(0, fit_u_fun->GetParameter(0));
 		fit_u_fun_pol->FixParameter(1, fit_u_fun->GetParameter(1));
 		fit_u_fun_pol->FixParameter(2, fit_u_fun->GetParameter(2));
 		fit_u_fun_pol->FixParameter(3, fit_u_fun->GetParameter(3));

		fit_u_fun_pol->Draw("same");

		h_sim_clone->SetLineColor(6);
		h_sim_clone->Scale(fit_u_fun->GetParameter(4));		
//		h_sim_clone->Scale(fit_u_fun->GetParameter(3));		
 	 	h_sim_clone->DrawCopy("same");

 		h_sim_clone->Add(fit_u_fun_pol);		
		h_sim_clone->SetLineColor(1);
 	 	h_sim_clone->DrawCopy("same");

	}



 	for (Int_t i =  phi_bin_num + 1; i <= 2*phi_bin_num; i++) {
	
		TPad* cur_pad1 = cc4->GetPad(i);

		cur_pad1->cd();
		cur_pad1->SetLeftMargin(0.1);
		cur_pad1->SetRightMargin(0.002);



		cout << i << "   " << i-1 << endl;

		const int in =  i - phi_bin_num;

		cout << in << endl;

//		TH1F* h_tgt[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
//		TH1F* h_sim[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");

 		TH1F* h_tgt_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
 		TH1F* h_tgt[in] = (TH1F*) h_tgt_tmp->Clone("h_tgt_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");
 		TH1F* h_sim[in] = (TH1F*) h_sim_tmp->Clone("h_sim_clone");


	 	h_tgt[in]->DrawCopy();  
	 	h_sim[in]->DrawCopy("same");

		peak1 = h_sim[in];

		Double_t fit_local_lo;
		Double_t fit_local_hi;

// 		fit_local_lo = fit_mid_1;
// 		fit_local_hi = fit_mid_2;

		TPad* tar_pad = (TPad*) cc3->GetPad(i);

  		fit_local_lo = Find_low_range(tar_pad);
  		fit_local_hi = Find_high_range(tar_pad);

 		TF1 *fit_u_fun = new TF1("fit_u_fun", ftotal_all, fit_local_lo, fit_local_hi, 5);

 		fit_u_fun->SetParameters(10.0, -50.0, 75, -36.0, 5);
		fit_u_fun->SetParLimits(4, 0, 60);

 		h_tgt[in]->Fit("fit_u_fun", "MEURN");
 
 	 	h_sim[in]->Draw("same");
 
 	 	TH1F *h_sim_clone = (TH1F*) h_sim[in]->Clone();
		
		TF1 *fit_u_fun_pol = new TF1("fit_u_fun_pol", ftotal_br, fit_local_lo, fit_local_hi, 4);

		fit_u_fun_pol->SetLineColor(4);

 		fit_u_fun_pol->FixParameter(0, fit_u_fun->GetParameter(0));
 		fit_u_fun_pol->FixParameter(1, fit_u_fun->GetParameter(1));
 		fit_u_fun_pol->FixParameter(2, fit_u_fun->GetParameter(2));
 		fit_u_fun_pol->FixParameter(3, fit_u_fun->GetParameter(3));

		fit_u_fun_pol->Draw("same");

		h_sim_clone->SetLineColor(6);
		h_sim_clone->Scale(fit_u_fun->GetParameter(4));		
 	 	h_sim_clone->DrawCopy("same");

 		h_sim_clone->Add(fit_u_fun_pol);		
		h_sim_clone->SetLineColor(1);
 	 	h_sim_clone->DrawCopy("same");

	}


 	for (Int_t i =  2*phi_bin_num + 1; i <= 3*phi_bin_num; i++) {

		TPad* cur_pad2 = cc4->GetPad(i);

		cur_pad2->cd();
		cur_pad2->SetLeftMargin(0.1);
		cur_pad2->SetRightMargin(0.002);



		cout << i << "   " << i-1 << endl;

		const int in =  i - (2*phi_bin_num + 1);
	

		cout << in << endl;

//		TH1F* h_tgt[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
//		TH1F* h_sim[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");

 		TH1F* h_tgt_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
 		TH1F* h_tgt[in] = (TH1F*) h_tgt_tmp->Clone("h_tgt_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");
 		TH1F* h_sim[in] = (TH1F*) h_sim_tmp->Clone("h_sim_clone");


	 	h_tgt[in]->DrawCopy();  
	 	h_sim[in]->DrawCopy("same");

		peak1 = h_sim[in];

		Double_t fit_local_lo;
		Double_t fit_local_hi;

//  		fit_local_lo = fit_high_1;
//  		fit_local_hi = fit_high_2;
 
		TPad* tar_pad = (TPad*) cc3->GetPad(i);

  		fit_local_lo = Find_low_range(tar_pad);
  		fit_local_hi = Find_high_range(tar_pad);

 		TF1 *fit_u_fun = new TF1("fit_u_fun", ftotal_all, fit_local_lo, fit_local_hi, 5);

 		fit_u_fun->SetParameters(10.0, -50.0, 75, -36.0, 5);
		fit_u_fun->SetParLimits(4, 0, 60);

 		h_tgt[in]->Fit("fit_u_fun", "MEURN");
 
 	 	h_sim[in]->Draw("same");
 
 	 	TH1F *h_sim_clone = (TH1F*) h_sim[in]->Clone();
		
		TF1 *fit_u_fun_pol = new TF1("fit_u_fun_pol", ftotal_br, fit_local_lo, fit_local_hi, 4);

		fit_u_fun_pol->SetLineColor(4);

 		fit_u_fun_pol->FixParameter(0, fit_u_fun->GetParameter(0));
 		fit_u_fun_pol->FixParameter(1, fit_u_fun->GetParameter(1));
 		fit_u_fun_pol->FixParameter(2, fit_u_fun->GetParameter(2));
 		fit_u_fun_pol->FixParameter(3, fit_u_fun->GetParameter(3));

		fit_u_fun_pol->Draw("same");

		h_sim_clone->SetLineColor(6);
		h_sim_clone->Scale(fit_u_fun->GetParameter(4));		
 	 	h_sim_clone->DrawCopy("same");

 		h_sim_clone->Add(fit_u_fun_pol);		
		h_sim_clone->SetLineColor(1);
 	 	h_sim_clone->DrawCopy("same");

	}


	cc4->Print(out_dir_str + "u_fit_phi_no_const_bg" + file_name+".png");

//	exit(0);

}





/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Invividual fitting with 2nd order poly

void u_phi_fit_2(TCanvas* cc3, TString file_name, Double_t fit_low_1, Double_t fit_low_2, Double_t fit_mid_1, Double_t fit_mid_2, Double_t fit_high_1, Double_t fit_high_2) {
	
//void u_phi_fit_2(TCanvas* cc3, TString file_name) {

	cout << "asddddddddddddddddd " << endl;

//	cc3->Draw();

	const Int_t u_bin_num_const   = u_bin_num;  
	const Int_t phi_bin_num_const = phi_bin_num;  

	TCanvas* cc4 = new TCanvas("cc4", "cc4", 1600, 800);
	cc4->Divide(phi_bin_num_const, u_bin_num_const, 0.003, 0.003);

	TH1F* h_tgt[6];
	TH1F* h_sim[6];

	for (Int_t i = 1; i <= phi_bin_num; i++) {
	
		TPad* cur_pad = cc4->GetPad(i);

		cur_pad->cd();
		cur_pad->SetLeftMargin(0.1);
		cur_pad->SetRightMargin(0.002);

		cout << i << "   " << i-1 << endl;

		const int in =  i -1;

		cout << in << endl;




//		TH1F* h_tgt[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
//		TH1F* h_sim[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");

 		TH1F* h_tgt_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
 		TH1F* h_tgt[in] = (TH1F*) h_tgt_tmp->Clone("h_tgt_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");
 		TH1F* h_sim[in] = (TH1F*) h_sim_tmp->Clone("h_sim_clone");
		

//		Find_low_range(tar_pad);


	 	h_tgt[in]->DrawCopy();  
	 	h_sim[in]->DrawCopy("same");

		peak1 = h_sim[in];

		Double_t fit_local_lo;
		Double_t fit_local_hi;

//  		fit_local_lo = fit_low_1;
//  		fit_local_hi = fit_low_2;

		TPad* tar_pad = (TPad*) cc3->GetPad(i);

  		fit_local_lo = Find_low_range(tar_pad);
  		fit_local_hi = Find_high_range(tar_pad);


		cout << "low limit " << fit_local_lo << "    " << fit_local_hi << endl;

//		exit(0);


//		TF1 *fit_u_fun = new TF1("fit_u_fun", ftotal_all, fit_local_lo, fit_local_hi, 5);
 		TF1 *fit_u_fun = new TF1("fit_u_fun", ftotal_all_2, fit_local_lo, fit_local_hi, 4);


// 		fit_u_fun->SetParameters(10.0, -50.0, 75, -36.0, 5);
		fit_u_fun->SetParLimits(3, 0, 60);

 		h_tgt[in]->Fit("fit_u_fun", "MEURN");
 
 	 	h_sim[in]->Draw("same");
 
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



 	for (Int_t i =  phi_bin_num + 1; i <= 2*phi_bin_num; i++) {
	
		TPad* cur_pad1 = cc4->GetPad(i);

		cur_pad1->cd();
		cur_pad1->SetLeftMargin(0.1);
		cur_pad1->SetRightMargin(0.002);



		cout << i << "   " << i-1 << endl;

		const int in =  i - phi_bin_num;

		cout << in << endl;

//		TH1F* h_tgt[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
//		TH1F* h_sim[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");

 		TH1F* h_tgt_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
 		TH1F* h_tgt[in] = (TH1F*) h_tgt_tmp->Clone("h_tgt_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");
 		TH1F* h_sim[in] = (TH1F*) h_sim_tmp->Clone("h_sim_clone");


	 	h_tgt[in]->DrawCopy();  
	 	h_sim[in]->DrawCopy("same");

		peak1 = h_sim[in];

		Double_t fit_local_lo;
		Double_t fit_local_hi;

//  		fit_local_lo = fit_mid_1;
//  		fit_local_hi = fit_mid_2;

		TPad* tar_pad = (TPad*) cc3->GetPad(i);
 
 		fit_local_lo = Find_low_range(tar_pad);
  		fit_local_hi = Find_high_range(tar_pad);

// 		TF1 *fit_u_fun = new TF1("fit_u_fun", ftotal_all, fit_local_lo, fit_local_hi, 5);
 		TF1 *fit_u_fun = new TF1("fit_u_fun", ftotal_all_2, fit_local_lo, fit_local_hi, 4);

 //		fit_u_fun->SetParameters(10.0, -50.0, 75, -36.0, 5);
		fit_u_fun->SetParLimits(4, 0, 60);

 		h_tgt[in]->Fit("fit_u_fun", "MEURN");
 
 	 	h_sim[in]->Draw("same");
 
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


 	for (Int_t i =  2*phi_bin_num + 1; i <= 3*phi_bin_num; i++) {

		TPad* cur_pad2 = cc4->GetPad(i);

		cur_pad2->cd();
		cur_pad2->SetLeftMargin(0.1);
		cur_pad2->SetRightMargin(0.002);



		cout << i << "   " << i-1 << endl;

		const int in =  i - (2*phi_bin_num + 1);
	

		cout << in << endl;

//		TH1F* h_tgt[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
//		TH1F* h_sim[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");

 		TH1F* h_tgt_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
 		TH1F* h_tgt[in] = (TH1F*) h_tgt_tmp->Clone("h_tgt_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");
 		TH1F* h_sim[in] = (TH1F*) h_sim_tmp->Clone("h_sim_clone");


	 	h_tgt[in]->DrawCopy();  
	 	h_sim[in]->DrawCopy("same");

		peak1 = h_sim[in];

		Double_t fit_local_lo;
		Double_t fit_local_hi;

//  		fit_local_lo = fit_high_1;
//  		fit_local_hi = fit_high_2;

		TPad* tar_pad = (TPad*) cc3->GetPad(i);

  		fit_local_lo = Find_low_range(tar_pad);
  		fit_local_hi = Find_high_range(tar_pad);


// 		TF1 *fit_u_fun = new TF1("fit_u_fun", ftotal_all, fit_local_lo, fit_local_hi, 5);
 		TF1 *fit_u_fun = new TF1("fit_u_fun", ftotal_all_2, fit_local_lo, fit_local_hi, 4);

// 		fit_u_fun->SetParameters(10.0, -50.0, 75, -36.0, 5);
		fit_u_fun->SetParLimits(4, 0, 60);

 		h_tgt[in]->Fit("fit_u_fun", "MEURN");
 
 	 	h_sim[in]->Draw("same");
 
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


	cc4->Print(out_dir_str + "u_fit_phi_no_const_bg" + file_name+".png");


}






























/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Individual fit with constant background 3rd order poly version

void u_phi_fit_const_bg(TCanvas* cc3, TString file_name, Double_t fit_low_1, Double_t fit_low_2, Double_t fit_mid_1, Double_t fit_mid_2, Double_t fit_high_1, Double_t fit_high_2) {
	
	cout << "asddddddddddddddddd " << endl;

//	cc3->Draw();

	const Int_t u_bin_num_const   = u_bin_num;  
	const Int_t phi_bin_num_const = phi_bin_num;  

	TCanvas* cc4 = new TCanvas("cc4", "cc4", 1600, 800);
	cc4->Divide(phi_bin_num_const, u_bin_num_const, 0.003, 0.003);

	TH1F* h_tgt[6];
	TH1F* h_sim[6];

	for (Int_t i = 1; i <= phi_bin_num; i++) {
	
		TPad* cur_pad = cc4->GetPad(i);

		cur_pad->cd();
		cur_pad->SetLeftMargin(0.1);
		cur_pad->SetRightMargin(0.002);



		cout << i << "   " << i-1 << endl;

		const int in =  i -1;

		cout << in << endl;

//		TH1F* h_tgt[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
//		TH1F* h_sim[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");

 		TH1F* h_tgt_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
 		TH1F* h_tgt[in] = (TH1F*) h_tgt_tmp->Clone("h_tgt_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");
 		TH1F* h_sim[in] = (TH1F*) h_sim_tmp->Clone("h_sim_clone");


	 	h_tgt[in]->DrawCopy();  
	 	h_sim[in]->DrawCopy("same");

//		fun_t1->Scale(0.01);
//		fun_t1->SetName("fun_t1");
//		fun_t1->Draw();


		peak1 = h_sim[in];

		Double_t fit_local_lo;
		Double_t fit_local_hi;

//  		fit_local_lo = fit_low_1;
//  		fit_local_hi = fit_low_2;

		TPad* tar_pad = (TPad*) cc3->GetPad(i);

  		fit_local_lo = Find_low_range(tar_pad);
  		fit_local_hi = Find_high_range(tar_pad);
		
 		TF1 *fit_u_fun1 = new TF1("fun1111", ftotal_const_bg_fun, fit_local_lo, fit_local_hi, 6);

		cout << fit_local_lo << "    " << fit_local_hi << endl;

//		exit(0);
		
 		fit_u_fun1->FixParameter(0, fun_t1->GetParameter(0));
 		fit_u_fun1->FixParameter(1, fun_t1->GetParameter(1));
 		fit_u_fun1->FixParameter(2, fun_t1->GetParameter(2));
 		fit_u_fun1->FixParameter(3, fun_t1->GetParameter(3));

//		fit_u_fun1->FixParameter(4, 0.02);

// 		fit_u_fun1->FixParameter(5, 0);

// 		fit_u_fun1->SetParLimits(4, 0, 1);
 		fit_u_fun1->SetParLimits(5, 0, fun_t1->GetParameter(4));
 



//		fit_u_fun1->Draw("same");


// 		TF1 *fa1 = new TF1("fa1","sin(x)/x",0,10); 
// 		fa1->Draw();
// 

// 		fit_u_fun->FixParameter(0, 1);
// 		fit_u_fun->FixParameter(1, 1);		
// 		fit_u_fun->FixParameter(2, 1);
// 		fit_u_fun->FixParameter(3, 1);
// 		fit_u_fun->FixParameter(4, 0.01);
//		fit_u_fun->FixParameter(5, 2);


//		fit_u_fun->Draw();


// 		fit_u_fun->SetParameters(10.0, -50.0, 75, -36.0, 5);
//		fit_u_fun->SetParLimits(4, 0, 60);
//		fit_u_fun->SetParLimits(5, 0, 1);
// 
  		h_tgt[in]->Fit("fun1111", "MEURN");

//		fit_u_fun1->Draw("same");

//  
//  	 	h_sim[in]->Draw("same");
//  
  	 	TH1F *h_sim_clone = (TH1F*) h_sim[in]->Clone();
// 		
 		TF1 *fit_u_fun_pol = new TF1("fit_u_fun_pol", ftotal_const_bg_bg_only, fit_local_lo, fit_local_hi, 5);
// 
 		fit_u_fun_pol->SetLineColor(4);
// 
 		fit_u_fun_pol->FixParameter(0, fit_u_fun1->GetParameter(0));
 		fit_u_fun_pol->FixParameter(1, fit_u_fun1->GetParameter(1));
 		fit_u_fun_pol->FixParameter(2, fit_u_fun1->GetParameter(2));
 		fit_u_fun_pol->FixParameter(3, fit_u_fun1->GetParameter(3));
 		fit_u_fun_pol->FixParameter(4, fit_u_fun1->GetParameter(4));
 
 		fit_u_fun_pol->Draw("same");



// 
 		h_sim_clone->SetLineColor(6);
 		h_sim_clone->Scale(fit_u_fun1->GetParameter(5));		
  	 	h_sim_clone->DrawCopy("same");
// 
  		h_sim_clone->Add(fit_u_fun_pol);		
 		h_sim_clone->SetLineColor(1);
  	 	h_sim_clone->DrawCopy("same");
// 
	}



	for (Int_t i = phi_bin_num + 1; i <= 2*phi_bin_num; i++) {

		TPad* cur_pad1 = cc4->GetPad(i);

		cur_pad1->cd();
		cur_pad1->SetLeftMargin(0.1);
		cur_pad1->SetRightMargin(0.002);



		cout << i << "   " << i-1 << endl;

		const int in =  i - (phi_bin_num + 1);

		cout << in << endl;

//		TH1F* h_tgt[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
//		TH1F* h_sim[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");

 		TH1F* h_tgt_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
 		TH1F* h_tgt[in] = (TH1F*) h_tgt_tmp->Clone("h_tgt_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");
 		TH1F* h_sim[in] = (TH1F*) h_sim_tmp->Clone("h_sim_clone");


	 	h_tgt[in]->DrawCopy();  
	 	h_sim[in]->DrawCopy("same");

//		fun_t1->Scale(0.01);
//		fun_t1->SetName("fun_t1");
//		fun_t1->Draw();


		peak1 = h_sim[in];

		Double_t fit_local_lo;
		Double_t fit_local_hi;

// 		fit_local_lo = fit_mid_1;
// 		fit_local_hi = fit_mid_2;

		TPad* tar_pad = (TPad*) cc3->GetPad(i);

  		fit_local_lo = Find_low_range(tar_pad);
  		fit_local_hi = Find_high_range(tar_pad);




 		TF1 *fit_u_fun1 = new TF1("fun1111", ftotal_const_bg_fun, fit_local_lo, fit_local_hi, 6);


 		fit_u_fun1->FixParameter(0, fun_t2->GetParameter(0));
 		fit_u_fun1->FixParameter(1, fun_t2->GetParameter(1));
 		fit_u_fun1->FixParameter(2, fun_t2->GetParameter(2));
 		fit_u_fun1->FixParameter(3, fun_t2->GetParameter(3));
// 		fit_u_fun1->SetParameter(4, 0.01);

// 		fit_u_fun1->SetParLimits(4, 0, 1);
 		fit_u_fun1->SetParLimits(5, 0, fun_t2->GetParameter(4));


//		fit_u_fun1->Draw("same");


// 		TF1 *fa1 = new TF1("fa1","sin(x)/x",0,10); 
// 		fa1->Draw();
// 

// 		fit_u_fun->FixParameter(0, 1);
// 		fit_u_fun->FixParameter(1, 1);		
// 		fit_u_fun->FixParameter(2, 1);
// 		fit_u_fun->FixParameter(3, 1);
// 		fit_u_fun->FixParameter(4, 0.01);
//		fit_u_fun->FixParameter(5, 2);


//		fit_u_fun->Draw();


// 		fit_u_fun->SetParameters(10.0, -50.0, 75, -36.0, 5);
//		fit_u_fun->SetParLimits(4, 0, 60);
//		fit_u_fun->SetParLimits(5, 0, 1);
// 
  		h_tgt[in]->Fit("fun1111", "MEURN");

//		fit_u_fun1->Draw("same");

//  
//  	 	h_sim[in]->Draw("same");
//  
  	 	TH1F *h_sim_clone = (TH1F*) h_sim[in]->Clone();
// 		
 		TF1 *fit_u_fun_pol = new TF1("fit_u_fun_pol", ftotal_const_bg_bg_only, fit_local_lo, fit_local_hi, 5);
// 
 		fit_u_fun_pol->SetLineColor(4);
// 
 		fit_u_fun_pol->FixParameter(0, fit_u_fun1->GetParameter(0));
 		fit_u_fun_pol->FixParameter(1, fit_u_fun1->GetParameter(1));
 		fit_u_fun_pol->FixParameter(2, fit_u_fun1->GetParameter(2));
 		fit_u_fun_pol->FixParameter(3, fit_u_fun1->GetParameter(3));
 		fit_u_fun_pol->FixParameter(4, fit_u_fun1->GetParameter(4));
 
 		fit_u_fun_pol->Draw("same");



// 
 		h_sim_clone->SetLineColor(6);
 		h_sim_clone->Scale(fit_u_fun1->GetParameter(5));		
  	 	h_sim_clone->DrawCopy("same");
// 
  		h_sim_clone->Add(fit_u_fun_pol);		
 		h_sim_clone->SetLineColor(1);
  	 	h_sim_clone->DrawCopy("same");


	}
// 
// 

 	for (Int_t i =  2*phi_bin_num + 1; i <= 3*phi_bin_num; i++) {

		TPad* cur_pad2 = cc4->GetPad(i);

		cur_pad2->cd();
		cur_pad2->SetLeftMargin(0.1);
		cur_pad2->SetRightMargin(0.002);



		cout << i << "   " << i-1 << endl;

		const int in =  i - (2*phi_bin_num + 1);

		cout << in << endl;

//		TH1F* h_tgt[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
//		TH1F* h_sim[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");

 		TH1F* h_tgt_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
 		TH1F* h_tgt[in] = (TH1F*) h_tgt_tmp->Clone("h_tgt_clone");

 		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");
 		TH1F* h_sim[in] = (TH1F*) h_sim_tmp->Clone("h_sim_clone");


	 	h_tgt[in]->DrawCopy();  
	 	h_sim[in]->DrawCopy("same");

//		fun_t1->Scale(0.01);
//		fun_t1->SetName("fun_t1");
//		fun_t1->Draw();


		peak1 = h_sim[in];

		Double_t fit_local_lo;
		Double_t fit_local_hi;

// 		fit_local_lo = fit_high_1;
// 		fit_local_hi = fit_high_2;

		TPad* tar_pad = (TPad*) cc3->GetPad(i);

  		fit_local_lo = Find_low_range(tar_pad);
  		fit_local_hi = Find_high_range(tar_pad);


 		TF1 *fit_u_fun1 = new TF1("fun1111", ftotal_const_bg_fun, fit_local_lo, fit_local_hi, 6);


 		fit_u_fun1->FixParameter(0, fun_t3->GetParameter(0));
 		fit_u_fun1->FixParameter(1, fun_t3->GetParameter(1));
 		fit_u_fun1->FixParameter(2, fun_t3->GetParameter(2));
 		fit_u_fun1->FixParameter(3, fun_t3->GetParameter(3));
 
// 		fit_u_fun1->SetParameter(4, 0.02);
// 		fit_u_fun1->FixParameter(4, 0.025);
// 		fit_u_fun1->SetParameter(5, 0.0);

 		fit_u_fun1->SetParLimits(5, 0.0, fun_t3->GetParameter(4));

//		fit_u_fun1->FixParameter(5, 0);


		
// 		fit_u_fun1->SetParLimits(4, 0, );

// 		fit_u_fun1->SetParLimits(4, 0, 1);

// 		fit_u_fun1->SetParameter(4, 0.01);

//		fit_u_fun1->Draw("same");


// 		TF1 *fa1 = new TF1("fa1","sin(x)/x",0,10); 
// 		fa1->Draw();
// 

// 		fit_u_fun->FixParameter(0, 1);
// 		fit_u_fun->FixParameter(1, 1);		
// 		fit_u_fun->FixParameter(2, 1);
// 		fit_u_fun->FixParameter(3, 1);
// 		fit_u_fun->FixParameter(4, 0.01);
//		fit_u_fun->FixParameter(5, 2);


//		fit_u_fun->Draw();


// 		fit_u_fun->SetParameters(10.0, -50.0, 75, -36.0, 5);
//		fit_u_fun->SetParLimits(4, 0, 60);
//		fit_u_fun->SetParLimits(5, 0, 100);
 
  		h_tgt[in]->Fit("fun1111", "MEURN");

//		fit_u_fun1->Draw("same");

//  
//  	 	h_sim[in]->Draw("same");
//  
  	 	TH1F *h_sim_clone = (TH1F*) h_sim[in]->Clone();
// 		
 		TF1 *fit_u_fun_pol = new TF1("fit_u_fun_pol", ftotal_const_bg_bg_only, fit_local_lo, fit_local_hi, 5);
// 
 		fit_u_fun_pol->SetLineColor(4);
// 
 		fit_u_fun_pol->FixParameter(0, fit_u_fun1->GetParameter(0));
 		fit_u_fun_pol->FixParameter(1, fit_u_fun1->GetParameter(1));
 		fit_u_fun_pol->FixParameter(2, fit_u_fun1->GetParameter(2));
 		fit_u_fun_pol->FixParameter(3, fit_u_fun1->GetParameter(3));
 		fit_u_fun_pol->FixParameter(4, fit_u_fun1->GetParameter(4));
 
 		fit_u_fun_pol->Draw("same");



// 
 		h_sim_clone->SetLineColor(6);
 		h_sim_clone->Scale(fit_u_fun1->GetParameter(5));		
  	 	h_sim_clone->DrawCopy("same");
// 
  		h_sim_clone->Add(fit_u_fun_pol);		
 		h_sim_clone->SetLineColor(1);
  	 	h_sim_clone->DrawCopy("same");

 	}

	
	cc4->Print(out_dir_str + "u_fit_phi_const_bg" + file_name+".png");

}




// 
// /*--------------------------------------------------*/
// /*--------------------------------------------------*/
// /// Individual fit with constant background 2nd order poly version
// 
// 
// void u_phi_fit_const_bg(TCanvas* cc3, TString file_name, Double_t fit_low_1, Double_t fit_low_2, Double_t fit_mid_1, Double_t fit_mid_2, Double_t fit_high_1, Double_t fit_high_2) {
// 	
// 	cout << "asddddddddddddddddd " << endl;
// 
// //	cc3->Draw();
// 
// 	TCanvas* cc4 = new TCanvas("cc4", "cc4", 1600, 800);
// 	cc4->Divide(phi_bin_num, u_bin_num, 0.003, 0.003);
// 
// 	TH1F* h_tgt[6];
// 	TH1F* h_sim[6];
// 
// 	for (Int_t i = 1; i <= phi_bin_num; i++) {
// 	
// 		TPad* cur_pad = cc4->GetPad(i);
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
// //		TH1F* h_tgt[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
// //		TH1F* h_sim[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");
// 
//  		TH1F* h_tgt_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
//  		TH1F* h_tgt[in] = (TH1F*) h_tgt_tmp->Clone("h_tgt_clone");
// 
//  		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");
//  		TH1F* h_sim[in] = (TH1F*) h_sim_tmp->Clone("h_sim_clone");
// 
// 
// 	 	h_tgt[in]->DrawCopy();  
// 	 	h_sim[in]->DrawCopy("same");
// 
// //		fun_t1->Scale(0.01);
// //		fun_t1->SetName("fun_t1");
// //		fun_t1->Draw();
// 
// 
// 		peak1 = h_sim[in];
// 
// 		Double_t fit_local_lo;
// 		Double_t fit_local_hi;
// 
//  		fit_local_lo = fit_low_1;
//  		fit_local_hi = fit_low_2;
// 
//  		TF1 *fit_u_fun1 = new TF1("fun1111", ftotal_const_bg_fun, fit_local_lo, fit_local_hi, 6);
// 
// 
//  		fit_u_fun1->FixParameter(0, fun_t1->GetParameter(0));
//  		fit_u_fun1->FixParameter(1, fun_t1->GetParameter(1));
//  		fit_u_fun1->FixParameter(2, fun_t1->GetParameter(2));
//  		fit_u_fun1->FixParameter(3, fun_t1->GetParameter(3));
// 
// 
// // 		fit_u_fun1->SetParLimits(4, 0, 1);
// // 		fit_u_fun1->SetParLimits(5, 0, fun_t1->GetParameter(5));
//  
// 
// 
// 
// //		fit_u_fun1->Draw("same");
// 
// 
// // 		TF1 *fa1 = new TF1("fa1","sin(x)/x",0,10); 
// // 		fa1->Draw();
// // 
// 
// // 		fit_u_fun->FixParameter(0, 1);
// // 		fit_u_fun->FixParameter(1, 1);		
// // 		fit_u_fun->FixParameter(2, 1);
// // 		fit_u_fun->FixParameter(3, 1);
// // 		fit_u_fun->FixParameter(4, 0.01);
// //		fit_u_fun->FixParameter(5, 2);
// 
// 
// //		fit_u_fun->Draw();
// 
// 
// // 		fit_u_fun->SetParameters(10.0, -50.0, 75, -36.0, 5);
// //		fit_u_fun->SetParLimits(4, 0, 60);
// //		fit_u_fun->SetParLimits(5, 0, 1);
// // 
//   		h_tgt[in]->Fit("fun1111", "MEURN");
// 
// //		fit_u_fun1->Draw("same");
// 
// //  
// //  	 	h_sim[in]->Draw("same");
// //  
//   	 	TH1F *h_sim_clone = (TH1F*) h_sim[in]->Clone();
// // 		
//  		TF1 *fit_u_fun_pol = new TF1("fit_u_fun_pol", ftotal_const_bg_bg_only, fit_local_lo, fit_local_hi, 5);
// // 
//  		fit_u_fun_pol->SetLineColor(4);
// // 
//  		fit_u_fun_pol->FixParameter(0, fit_u_fun1->GetParameter(0));
//  		fit_u_fun_pol->FixParameter(1, fit_u_fun1->GetParameter(1));
//  		fit_u_fun_pol->FixParameter(2, fit_u_fun1->GetParameter(2));
//  		fit_u_fun_pol->FixParameter(3, fit_u_fun1->GetParameter(3));
//  		fit_u_fun_pol->FixParameter(4, fit_u_fun1->GetParameter(4));
//  
//  		fit_u_fun_pol->Draw("same");
// 
// 
// 
// // 
//  		h_sim_clone->SetLineColor(6);
//  		h_sim_clone->Scale(fit_u_fun1->GetParameter(5));		
//   	 	h_sim_clone->DrawCopy("same");
// // 
//   		h_sim_clone->Add(fit_u_fun_pol);		
//  		h_sim_clone->SetLineColor(1);
//   	 	h_sim_clone->DrawCopy("same");
// // 
// 	}
// 
// 
// 
// 	for (Int_t i = phi_bin_num + 1; i <= 2*phi_bin_num; i++) {
// 
// 		TPad* cur_pad1 = cc4->GetPad(i);
// 
// 		cur_pad1->cd();
// 		cur_pad1->SetLeftMargin(0.1);
// 		cur_pad1->SetRightMargin(0.002);
// 
// 
// 
// 		cout << i << "   " << i-1 << endl;
// 
// 		const int in =  i - (phi_bin_num + 1);
// 
// 		cout << in << endl;
// 
// //		TH1F* h_tgt[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
// //		TH1F* h_sim[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");
// 
//  		TH1F* h_tgt_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
//  		TH1F* h_tgt[in] = (TH1F*) h_tgt_tmp->Clone("h_tgt_clone");
// 
//  		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");
//  		TH1F* h_sim[in] = (TH1F*) h_sim_tmp->Clone("h_sim_clone");
// 
// 
// 	 	h_tgt[in]->DrawCopy();  
// 	 	h_sim[in]->DrawCopy("same");
// 
// //		fun_t1->Scale(0.01);
// //		fun_t1->SetName("fun_t1");
// //		fun_t1->Draw();
// 
// 
// 		peak1 = h_sim[in];
// 
// 		Double_t fit_local_lo;
// 		Double_t fit_local_hi;
// 
//  		fit_local_lo = fit_mid_1;
//  		fit_local_hi = fit_mid_2;
// 
//  		TF1 *fit_u_fun1 = new TF1("fun1111", ftotal_const_bg_fun, fit_local_lo, fit_local_hi, 6);
// 
// 
//  		fit_u_fun1->FixParameter(0, fun_t2->GetParameter(0));
//  		fit_u_fun1->FixParameter(1, fun_t2->GetParameter(1));
//  		fit_u_fun1->FixParameter(2, fun_t2->GetParameter(2));
//  		fit_u_fun1->FixParameter(3, fun_t2->GetParameter(3));
// // 		fit_u_fun1->SetParameter(4, 0.01);
// 
// // 		fit_u_fun1->SetParLimits(4, 0, 1);
// // 		fit_u_fun1->SetParLimits(5, 0, fun_t2->GetParameter(5));
// 
// 
// //		fit_u_fun1->Draw("same");
// 
// 
// // 		TF1 *fa1 = new TF1("fa1","sin(x)/x",0,10); 
// // 		fa1->Draw();
// // 
// 
// // 		fit_u_fun->FixParameter(0, 1);
// // 		fit_u_fun->FixParameter(1, 1);		
// // 		fit_u_fun->FixParameter(2, 1);
// // 		fit_u_fun->FixParameter(3, 1);
// // 		fit_u_fun->FixParameter(4, 0.01);
// //		fit_u_fun->FixParameter(5, 2);
// 
// 
// //		fit_u_fun->Draw();
// 
// 
// // 		fit_u_fun->SetParameters(10.0, -50.0, 75, -36.0, 5);
// //		fit_u_fun->SetParLimits(4, 0, 60);
// //		fit_u_fun->SetParLimits(5, 0, 1);
// // 
//   		h_tgt[in]->Fit("fun1111", "MEURN");
// 
// //		fit_u_fun1->Draw("same");
// 
// //  
// //  	 	h_sim[in]->Draw("same");
// //  
//   	 	TH1F *h_sim_clone = (TH1F*) h_sim[in]->Clone();
// // 		
//  		TF1 *fit_u_fun_pol = new TF1("fit_u_fun_pol", ftotal_const_bg_bg_only, fit_local_lo, fit_local_hi, 5);
// // 
//  		fit_u_fun_pol->SetLineColor(4);
// // 
//  		fit_u_fun_pol->FixParameter(0, fit_u_fun1->GetParameter(0));
//  		fit_u_fun_pol->FixParameter(1, fit_u_fun1->GetParameter(1));
//  		fit_u_fun_pol->FixParameter(2, fit_u_fun1->GetParameter(2));
//  		fit_u_fun_pol->FixParameter(3, fit_u_fun1->GetParameter(3));
//  		fit_u_fun_pol->FixParameter(4, fit_u_fun1->GetParameter(4));
//  
//  		fit_u_fun_pol->Draw("same");
// 
// 
// 
// // 
//  		h_sim_clone->SetLineColor(6);
//  		h_sim_clone->Scale(fit_u_fun1->GetParameter(5));		
//   	 	h_sim_clone->DrawCopy("same");
// // 
//   		h_sim_clone->Add(fit_u_fun_pol);		
//  		h_sim_clone->SetLineColor(1);
//   	 	h_sim_clone->DrawCopy("same");
// 
// 
// 	}
// // 
// // 
// 
//  	for (Int_t i =  2*phi_bin_num + 1; i <= 3*phi_bin_num; i++) {
// 
// 		TPad* cur_pad2 = cc4->GetPad(i);
// 
// 		cur_pad2->cd();
// 		cur_pad2->SetLeftMargin(0.1);
// 		cur_pad2->SetRightMargin(0.002);
// 
// 
// 
// 		cout << i << "   " << i-1 << endl;
// 
// 		const int in =  i - (2*phi_bin_num + 1);
// 
// 		cout << in << endl;
// 
// //		TH1F* h_tgt[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
// //		TH1F* h_sim[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");
// 
//  		TH1F* h_tgt_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
//  		TH1F* h_tgt[in] = (TH1F*) h_tgt_tmp->Clone("h_tgt_clone");
// 
//  		TH1F* h_sim_tmp = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");
//  		TH1F* h_sim[in] = (TH1F*) h_sim_tmp->Clone("h_sim_clone");
// 
// 
// 	 	h_tgt[in]->DrawCopy();  
// 	 	h_sim[in]->DrawCopy("same");
// 
// //		fun_t1->Scale(0.01);
// //		fun_t1->SetName("fun_t1");
// //		fun_t1->Draw();
// 
// 
// 		peak1 = h_sim[in];
// 
// 		Double_t fit_local_lo;
// 		Double_t fit_local_hi;
// 
//  		fit_local_lo = fit_high_1;
//  		fit_local_hi = fit_high_2;
// 
//  		TF1 *fit_u_fun1 = new TF1("fun1111", ftotal_const_bg_fun, fit_local_lo+0.3, fit_local_hi, 6);
// 
// 
//  		fit_u_fun1->FixParameter(0, fun_t3->GetParameter(0));
//  		fit_u_fun1->FixParameter(1, fun_t3->GetParameter(1));
//  		fit_u_fun1->FixParameter(2, fun_t3->GetParameter(2));
//  		fit_u_fun1->FixParameter(3, fun_t3->GetParameter(3));
// 
// 
// // 		fit_u_fun1->SetParLimits(4, 0, 1);
// // 		fit_u_fun1->SetParLimits(5, 0, fun_t3->GetParameter(5));
// 
// // 		fit_u_fun1->SetParameter(4, 0.01);
// 
// //		fit_u_fun1->Draw("same");
// 
// 
// // 		TF1 *fa1 = new TF1("fa1","sin(x)/x",0,10); 
// // 		fa1->Draw();
// // 
// 
// // 		fit_u_fun->FixParameter(0, 1);
// // 		fit_u_fun->FixParameter(1, 1);		
// // 		fit_u_fun->FixParameter(2, 1);
// // 		fit_u_fun->FixParameter(3, 1);
// // 		fit_u_fun->FixParameter(4, 0.01);
// //		fit_u_fun->FixParameter(5, 2);
// 
// 
// //		fit_u_fun->Draw();
// 
// 
// // 		fit_u_fun->SetParameters(10.0, -50.0, 75, -36.0, 5);
// //		fit_u_fun->SetParLimits(4, 0, 60);
// //		fit_u_fun->SetParLimits(5, 0, 100);
//  
//   		h_tgt[in]->Fit("fun1111", "MEURN");
// 
// //		fit_u_fun1->Draw("same");
// 
// //  
// //  	 	h_sim[in]->Draw("same");
// //  
//   	 	TH1F *h_sim_clone = (TH1F*) h_sim[in]->Clone();
// // 		
//  		TF1 *fit_u_fun_pol = new TF1("fit_u_fun_pol", ftotal_const_bg_bg_only, fit_local_lo, fit_local_hi, 5);
// // 
//  		fit_u_fun_pol->SetLineColor(4);
// // 
//  		fit_u_fun_pol->FixParameter(0, fit_u_fun1->GetParameter(0));
//  		fit_u_fun_pol->FixParameter(1, fit_u_fun1->GetParameter(1));
//  		fit_u_fun_pol->FixParameter(2, fit_u_fun1->GetParameter(2));
//  		fit_u_fun_pol->FixParameter(3, fit_u_fun1->GetParameter(3));
//  		fit_u_fun_pol->FixParameter(4, fit_u_fun1->GetParameter(4));
//  
//  		fit_u_fun_pol->Draw("same");
// 
// 
// 
// // 
//  		h_sim_clone->SetLineColor(6);
//  		h_sim_clone->Scale(fit_u_fun1->GetParameter(5));		
//   	 	h_sim_clone->DrawCopy("same");
// // 
//   		h_sim_clone->Add(fit_u_fun_pol);		
//  		h_sim_clone->SetLineColor(1);
//   	 	h_sim_clone->DrawCopy("same");
// 
//  	}
// 
// 	
// 	cc4->Print(out_dir_str + "u_fit_phi_const_bg" + file_name+".png");
// 
// }




/*--------------------------------------------------*/
/*--------------------------------------------------*/
///
void u_phi_fit(TCanvas* cc3) {
	
	cout << "asddddddddddddddddd " << endl;

	cc3->Draw();

	TCanvas* cc4 = new TCanvas("cc4", "cc4", 1600, 800);
	cc4->Divide(phi_bin_num, 1, 0.003);

	TH1F* h_tgt[6];
	TH1F* h_sim[6];

	for (Int_t i = 1; i <= phi_bin_num; i++) {
	
		TPad* cur_pad = cc4->GetPad(i);

		cur_pad->cd();
		cur_pad->SetLeftMargin(0.1);
		cur_pad->SetRightMargin(0.002);



		cout << i << "   " << i-1 << endl;

		const int in =  i -1;

		cout << in << endl;

		TH1F* h_tgt[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss");
		TH1F* h_sim[in] = (TH1F*) cc3->GetPad(i)->GetPrimitive("h_miss_sim");

	 	h_tgt[in]->Draw();  
	 	h_sim[in]->Draw("same");

	}




	cc4->cd(1);
	
	TH1F* h_tgt_1 = (TH1F*) cc3->GetPad(3)->GetPrimitive("h_miss");
	TH1F* h_sim_1 = (TH1F*) cc3->GetPad(3)->GetPrimitive("h_miss_sim");

   	peak= (TH1F*) cc3->GetPad(3)->GetPrimitive("h_miss_sim");
	
	
// 	h_tgt_1->Draw();  
// 	h_sim_1->Draw("same");
// 
 	TF1 *uphi = new TF1("uphi", ftotal, 0.56, 0.8 ,5);
 	h_tgt_1->Fit("uphi","MUER");   


// 	cc4->cd(2);
// 	TH1F* h_tgt_2 = (TH1F*) cc3->GetPad(2)->GetPrimitive("h_miss");
// 	TH1F* h_sim = (TH1F*) cc3->GetPad(2)->GetPrimitive("h_miss_sim");
// 
// //   	peak2= (TH1F*) cc3->GetPad(1)->GetPrimitive("h_miss_sim");
// 
// 	h_tgt_2->Draw();  
// 	h_sim->Draw("same");
// 




// // 	TF1 *uphi = new TF1("uphi", ftotal, 0.56, 0.8 ,5);
// // 	h_tgt->Fit("uphi","MUER"); 
// 


// /*--------------------------------------------------*/
// /*--------------------------------------------------*/
// /*--------------------------------------------------*/
// /// Fitting Part
//  	ROOT::Fit::BinData data; 
//  	ROOT::Fit::FillData(data, h_tgt_1); 
//  	ROOT::Fit::FillData(data, h_tgt_2); 
// // 	
// // 	cout << "data size is " << data.Size() << endl;
// // 	
// 	TF1 * f1 = new TF1("f1", ftotal_all_all, 0.56, 0.8 ,6);
// 
// 	f1->SetParameters(0.05, -0.23, 0.43, -0.23, 2, 2.5);
// 
// //	TF1 * f1 = new TF1("f1","gaus",-3,3);
// //	f1->SetParameters(1,0,1);
// // 	
//  	ROOT::Math::WrappedTF1 wf(*f1);
// // 	
//  	ROOT::Fit::Fitter fitter;
//  	fitter.SetFunction(wf);
// // 	
//  	fitter.Fit(data);
//  	ROOT::Fit::FitResult result = fitter.Result();
//  	result.Print(std::cout);
// 
// 
// 
// 

}



/*--------------------------------------------------*/

Double_t Find_low_range_full(TPad* tar_pad) {

  	TH1F* h_tt = (TH1F*) tar_pad->GetPrimitive("u_phi_target");
	
 	Int_t lo_end_bin = h_tt->FindFirstBinAbove(0.0);
 	Int_t hi_end_bin = h_tt->FindLastBinAbove(0.0);

	Int_t low_bin;
	Double_t low_limit;

	Int_t bin_itt=lo_end_bin;	
	
	Double_t yield_all =  h_tt->Integral(lo_end_bin, hi_end_bin);

// 	while( h_tt->Integral(lo_end_bin, bin_itt) < 0.02*yield_all) {
// 		bin_itt++;
// 	}


	///*--------------------------------------------------*/
	///

	while( h_tt->Integral(lo_end_bin, bin_itt) < fit_lo_global * yield_all) {
		bin_itt++;
	}


// 	cout << "Sim Lo Bin: "<< lo_end_bin << endl;
// 	cout << "Sim Hi Bin: "<< hi_end_bin << endl;
// 	cout << "10 percent from Lo: "<< bin_itt << endl;
	
	low_bin = bin_itt;
 
 	low_limit = h_tt->GetBinCenter(low_bin); 
 
	if(low_limit < 0.6){
		
		low_limit = 0.52;

	}





	return low_limit;

}


/*--------------------------------------------------*/

Double_t Find_high_range_full(TPad* tar_pad) {

  	TH1F* h_tt = (TH1F*) tar_pad->GetPrimitive("u_phi_target");
	
 	Int_t lo_end_bin = h_tt->FindFirstBinAbove(0.0);
 	Int_t hi_end_bin = h_tt->FindLastBinAbove(0.0);

	Int_t hi_bin;
	Double_t hi_limit;

	Int_t bin_itt = hi_end_bin;	
	
	Double_t yield_all =  h_tt->Integral(lo_end_bin, hi_end_bin);

// 	while( h_tt->Integral(bin_itt, hi_end_bin) < 0.02*yield_all) {
// 		bin_itt--;
// 	}


	///*--------------------------------------------------*/
	//

// 	while( h_tt->Integral(bin_itt, hi_end_bin) < 0.025*yield_all) {
// 		bin_itt--;
// 	}

	while( h_tt->Integral(bin_itt, hi_end_bin) < fit_hi_global * yield_all) {
		bin_itt--;
	}



	cout << "Sim Lo Bin: "<< lo_end_bin << endl;
	cout << "Sim Hi Bin: "<< hi_end_bin << endl;
	cout << "10 percent from Hi: "<< bin_itt << endl;
	
	hi_bin = bin_itt;
 
 	hi_limit = h_tt->GetBinCenter(hi_bin); 
 

//	exit(0);
 		
	return hi_limit;

}











/*--------------------------------------------------*/

Double_t Find_low_range(TPad* tar_pad) {


  	TH1F* h_tgt_tt = (TH1F*) tar_pad->GetPrimitive("u_phi_target");
// 
  	TH1F* h_sim_tt = (TH1F*) tar_pad->GetPrimitive("u_phi_sim_omega");
// 
// 	
	
	Int_t tgt_lo_end_bin = h_tgt_tt->FindFirstBinAbove(0.0);
 	Int_t tgt_hi_end_bin = h_tgt_tt->FindLastBinAbove(0.0); 

 	Double_t tgt_lo_end = h_tgt_tt->GetBinCenter(tgt_lo_end_bin);
 	Double_t tgt_hi_end = h_tgt_tt->GetBinCenter(tgt_hi_end_bin);

//	Int_t tgt_mid_bin = (h_tgt_tt->FindFirstBinAbove(0.001) + h_tgt_tt->FindLastBinAbove(0.001))/2; 

	Double_t tgt_mid =  h_tgt_tt->GetMean(1); 
	Double_t tgt_RMS =  h_tgt_tt->GetRMS(1); 
	Double_t omega_mass =  0.783; 


	Int_t tgt_max_bin = h_tgt_tt->GetMaximumBin(); 

 	Int_t sim_lo_end_bin = h_sim_tt->FindFirstBinAbove(0.0);
 	Int_t sim_hi_end_bin = h_sim_tt->FindLastBinAbove(0.0);

 	Double_t sim_lo_end = h_sim_tt->GetBinCenter( sim_lo_end_bin );
 	Double_t sim_hi_end = h_sim_tt->GetBinCenter( sim_hi_end_bin );

	Int_t sim_max_bin = h_sim_tt->GetMaximumBin(); 

	
// 	cout << "lllllllllllll     " << tgt_lo_end     << "    " << tgt_hi_end     << "    " << tgt_mid_bin << "    " << tgt_max_bin << endl;
// 	cout << "lllllllllllll     " << sim_lo_end_bin << "    " << sim_hi_end_bin << endl;
// 	cout << "lllllllllllll     " << sim_lo_end     << "    " << sim_hi_end     << "    " << sim_max_bin << endl;
 

	Int_t low_bin;
	Double_t low_limit;

	if (h_sim_tt->GetBinCenter(sim_max_bin) < tgt_mid) {

//		low_bin = sim_lo_end_bin - 1;
		
// 		if ( sim_lo_end - tgt_lo_end <= 0.5) {
// 			low_bin = sim_lo_end_bin + 2;
// 		} else {
// 			low_bin = sim_lo_end_bin + 3 ;
// 		}


//		if ( (tgt_mid - h_sim_tt->GetBinCenter(sim_max_bin)) >= 0.75*tgt_RMS ) {
		if ( (tgt_mid - h_sim_tt->GetBinCenter(sim_max_bin)) >= 1.2*tgt_RMS ) {
			low_bin = sim_max_bin - 2;
		} else {
//			low_bin = sim_lo_end_bin + 10 ;
			low_bin = h_tgt_tt->FindBin(tgt_mid - 1.5*tgt_RMS);
		}


	} else {

		////*--------------------------------------------------*/
		/// If boudary of simulation is close to the data boundary
		
// 		if ( tgt_hi_end_bin - sim_hi_end_bin <= 3) {
// 			low_bin = sim_lo_end_bin + 1 ;
// 		} else {
// 			low_bin = sim_lo_end_bin ;
// 		}

// 		if (h_sim_tt->GetBinCenter(sim_max_bin)-tgt_mid >= tgt_RMS) {
//  			low_bin =  h_tgt_tt->FindBin(tgt_mid-tgt_RMS);
// 		} else {
//  			low_bin = sim_lo_end_bin ;
// 		}
		
		if ((h_sim_tt->GetBinCenter(sim_max_bin)-tgt_mid) >= 0.75*tgt_RMS) {
 			low_bin =  h_tgt_tt->FindBin(tgt_mid - 1.5*tgt_RMS);
		} else {
// 			low_bin = sim_lo_end_bin ;
 			low_bin =  h_tgt_tt->FindBin(tgt_mid - 1.5*tgt_RMS);
		}

	}


	if ( fabs(omega_mass - h_sim_tt->GetBinCenter(sim_max_bin)) > tgt_RMS) {

		low_bin = h_tgt_tt->FindBin(tgt_mid - 1.5*tgt_RMS);


		cout << "!!!!!!!!!!!!1 " << omega_mass << "    " << tgt_mid << "    " << h_sim_tt->GetBinCenter(sim_max_bin) << "     " <<  tgt_RMS << "   ???  " << fabs(omega_mass - h_sim_tt->GetMaximum()) << endl;

//		exit(0);

	} 


// 	if (tgt_mid > 0.91) {
// 		cout << "!!!!!!!!!!!!1 " << omega_mass << "    " << tgt_mid << "    " << h_sim_tt->GetBinCenter(sim_max_bin) << "     " <<  tgt_RMS << "   ???  " << fabs(omega_mass - h_sim_tt->GetMaximum()) << endl;
// 
// 		exit(0);
// 	}
// 

//	cout << "aaaaaaaaaaaaaaaaa " << sim_lo_end_bin << "    " << tgt_lo_end_bin << "   " << sim_lo_end_bin - tgt_lo_end_bin  << endl;

	low_limit = h_sim_tt->GetBinCenter(low_bin); 

//	exit(0);
		
	return low_limit;

}


Double_t Find_high_range(TPad* tar_pad) {


  	TH1F* h_tgt_tt = (TH1F*) tar_pad->GetPrimitive("u_phi_target");
// 
  	TH1F* h_sim_tt = (TH1F*) tar_pad->GetPrimitive("u_phi_sim_omega");
// 
// 	
//  	Double_t tgt_lo_end = h_tgt_tt->GetBinCenter( h_tgt_tt->FindFirstBinAbove(0.0) );
//  	Double_t tgt_hi_end = h_tgt_tt->GetBinCenter( h_tgt_tt->FindLastBinAbove(0.0) );
// 


	Double_t thres = 0.0;	

	Int_t tgt_lo_end_bin = h_tgt_tt->FindFirstBinAbove(thres);
	Int_t tgt_hi_end_bin = h_tgt_tt->FindLastBinAbove(thres); 

	Double_t tgt_lo_end = h_tgt_tt->GetBinCenter(tgt_lo_end_bin);
	Double_t tgt_hi_end = h_tgt_tt->GetBinCenter(tgt_hi_end_bin);

//	Int_t tgt_mid_bin = (h_tgt_tt->FindFirstBinAbove(thres) + h_tgt_tt->FindLastBinAbove(thres))/2; 

	Double_t tgt_mid =  h_tgt_tt->GetMean(1); 
	Double_t tgt_RMS =  h_tgt_tt->GetRMS(1); 
	Double_t omega_mass =  0.783; 

	Double_t tgt_max_bin = h_tgt_tt->GetMaximumBin(); 

 	Int_t sim_lo_end_bin = h_sim_tt->FindFirstBinAbove(thres);
 	Int_t sim_hi_end_bin = h_sim_tt->FindLastBinAbove(thres);

 	Double_t sim_lo_end = h_sim_tt->GetBinCenter( sim_lo_end_bin );
 	Double_t sim_hi_end = h_sim_tt->GetBinCenter( sim_hi_end_bin );

	Int_t sim_max_bin = h_sim_tt->GetMaximumBin(); 

	
// 	cout << "hhhhhhhhhhhhh     " << tgt_lo_end     << "    " << tgt_hi_end << "    " << tgt_mid_bin << "    " << tgt_max_bin << endl;
// 	cout << "hhhhhhhhhhhhh     " << sim_lo_end_bin << "    " << sim_hi_end_bin       << endl;
// 	cout << "hhhhhhhhhhhhh     " << sim_lo_end     << "    " << sim_hi_end << "    " << sim_max_bin << endl;


	Int_t  high_bin;
	Double_t high_limit;

	if (h_sim_tt->GetBinCenter(sim_max_bin) < tgt_mid) {

//		high_bin = tgt_mid_bin + 8;

		if ((tgt_mid - h_sim_tt->GetBinCenter(sim_max_bin)) >= 0.75*tgt_RMS ) {
			high_bin = h_tgt_tt->FindBin(tgt_mid + 1.5*tgt_RMS);
		} else {
//			high_bin = sim_hi_end_bin ;
			high_bin = h_tgt_tt->FindBin(tgt_mid + 1.5*tgt_RMS);
		}



	} else {
		
//		high_bin = sim_hi_end_bin + 1;
		
// 		if ( tgt_hi_end_bin - sim_hi_end_bin <= 3) {
// 			high_bin = sim_hi_end_bin-2;
// 		} else {
// 			high_bin = sim_hi_end_bin-1;
// 		}

//		if ((h_sim_tt->GetBinCenter(sim_max_bin)-tgt_mid) >= 0.75*(tgt_RMS)) {
		if ((h_sim_tt->GetBinCenter(sim_max_bin)-tgt_mid) >= tgt_RMS) {
			high_bin = sim_max_bin+2;
		} else {
//			high_bin = sim_hi_end_bin - 10;
			high_bin = h_tgt_tt->FindBin(tgt_mid + 1.4*tgt_RMS);
		}

	}

	if ( fabs(omega_mass - h_sim_tt->GetBinCenter(sim_max_bin)) > tgt_RMS) {

		high_bin = h_tgt_tt->FindBin(tgt_mid + 1.5*tgt_RMS);

	}

	high_limit = h_sim_tt->GetBinCenter(high_bin); 

	return high_limit;

}



/*--------------------------------------------------*/

Double_t Find_low_range_xphsp(TPad* tar_pad) {

//  	TH1F* h_tgt_tt = (TH1F*) tar_pad->GetPrimitive("h_miss");
//		low_limit = h_tgt_tt->GetMean(1)-0.1; 


//  	TH1F* h_sim_tt = (TH1F*) tar_pad->GetPrimitive("h_miss_sim");

	
//	Int_t lo_end_bin = h_sim_tt->GetMaximumBin();

//	Int_t lo_end_bin = h_sim_tt->GetMaximumBin();


 	TH1F* h_sim_tt = (TH1F*) tar_pad->GetPrimitive("h_miss_sim_xphsp");
	Int_t lo_end_bin = h_sim_tt->FindFirstBinAbove(0.0);

	Double_t low_limit = h_sim_tt->GetBinCenter(lo_end_bin+4); 



	return low_limit;

}




/*--------------------------------------------------*/

Double_t Find_high_range_xphsp(TPad* tar_pad) {

//  	TH1F* h_tgt_tt = (TH1F*) tar_pad->GetPrimitive("h_miss");
//		high_limit = h_tgt_tt->GetMean(1)+0.1;; 

//   	TH1F* h_sim_tt = (TH1F*) tar_pad->GetPrimitive("h_miss_sim");
// 	
// 
// 	Int_t hi_end_bin = h_sim_tt->GetMaximumBin();
// 

  	TH1F* h_sim_tt = (TH1F*) tar_pad->GetPrimitive("h_miss_sim_xphsp");

	Int_t hi_end_bin = h_sim_tt->FindLastBinAbove(0.0);
 
	Double_t high_limit = h_sim_tt->GetBinCenter(hi_end_bin-4); 

	return high_limit;

}





/*--------------------------------------------------*/
/*--------------------------------------------------*/

void average_kin(TFile* file_tmp, TString file_name) {

	ofstream kin_file;

	kin_file.open("kindata/kindata.pl_" + file_name + ".dat", std::fstream::out); 

	kin_file << " 1.000000" << endl;

	Double_t w_ave, w_ave_err, q2_ave，q2_ave_err, u_ave, u_ave_err, t_ave, t_ave_err;

	TCanvas* can_w;
	can_w = (TCanvas*)file_tmp->Get("W_u_bin");

	TCanvas* can_q2;
	can_q2 = (TCanvas*)file_tmp->Get("Q2_u_bin");

	TCanvas* can_u;
	can_u = (TCanvas*)file_tmp->Get("u_u_bin");

	TCanvas* can_t;
	can_t = (TCanvas*)file_tmp->Get("t_u_bin");


	cout << u_bin_num << endl;

//	exit(0);


  	for (Int_t i = 1; i <=u_bin_num; i++) {
 
//  		TH1F* h_w  = (TH1F*) can_w ->GetPad(i)->GetPrimitive("u_target");
//  		TH1F* h_q2 = (TH1F*) can_q2->GetPad(i)->GetPrimitive("u_target");
//  		TH1F* h_u  = (TH1F*) can_u ->GetPad(i)->GetPrimitive("u_target");

 		TH1F* h_w  = (TH1F*) can_w ->GetPad(i)->GetPrimitive("u_sim_omega");
 		TH1F* h_q2 = (TH1F*) can_q2->GetPad(i)->GetPrimitive("u_sim_omega");
 		TH1F* h_u  = (TH1F*) can_u ->GetPad(i)->GetPrimitive("u_sim_omega");
 		TH1F* h_t  = (TH1F*) can_t ->GetPad(i)->GetPrimitive("u_sim_omega");

		w_ave      = h_w->GetMean();
		w_ave_err  = h_w->GetRMS();

		q2_ave     = h_q2->GetMean();
		q2_ave_err = h_q2->GetRMS();

		u_ave      = h_u->GetMean();
		u_ave_err  = h_u->GetRMS();

		t_ave      = h_t->GetMean();
		t_ave_err  = h_t->GetRMS();

		cout << w_ave << " " << w_ave_err << " " << q2_ave << " " << q2_ave_err 
			 << " " << u_ave << " " << u_ave_err << " " << t_ave << " " << t_ave_err << endl;

		
// 		kin_file << fixed << setprecision(5) << " " << w_ave << " " << w_ave_err << " " << q2_ave << " " << q2_ave_err 
// 			     << " " << u_ave << " " << u_ave_err  << " " << Form("%.2f", u_ave) << endl;

		kin_file << fixed << setprecision(5) << " " << w_ave << " " << w_ave_err << " " 
				 << q2_ave << " " << q2_ave_err << " " << t_ave << " " << t_ave_err << " " 
                 << u_ave  << " " << u_ave_err  << " "<< Form("%.2f", u_ave) << endl;

		delete h_w;
		delete h_q2;
  
	}

//	exit(0);


	delete can_w;
	delete can_q2;



	kin_file << " 1.000000" << endl;
	kin_file.close();



}






/*--------------------------------------------------*/
/// Find integration range low limit
Double_t Find_int_lo(TPad* tar_pad) {

  	TH1F* h_sim_tt = (TH1F*) tar_pad->GetPrimitive("u_phi_sim_omega");
	
 	Int_t sim_lo_end_bin = h_sim_tt->FindFirstBinAbove(0.0);
 	Int_t sim_hi_end_bin = h_sim_tt->FindLastBinAbove(0.0);

	Int_t low_bin;
	Double_t low_limit;

	Int_t bin_itt=sim_lo_end_bin;	
	
	Double_t sim_yield_all =  h_sim_tt->Integral(sim_lo_end_bin, sim_hi_end_bin);

	while( h_sim_tt->Integral(sim_lo_end_bin, bin_itt) < 0.02*sim_yield_all) {
		bin_itt++;
	}

	cout << "Sim Lo Bin: "<< sim_lo_end_bin << endl;
	cout << "Sim Hi Bin: "<< sim_hi_end_bin << endl;
	cout << "10 percent from Lo: "<< bin_itt << endl;

	
	low_bin = bin_itt;
 
 	low_limit = h_sim_tt->GetBinCenter(low_bin); 
 
//	exit(0);
 		
	return low_limit;

}


///*--------------------------------------------------*/
/// Find integration range high limit
Double_t Find_int_hi(TPad* tar_pad) {

  	TH1F* h_sim_tt = (TH1F*) tar_pad->GetPrimitive("u_phi_sim_omega");
	
 	Int_t sim_lo_end_bin = h_sim_tt->FindFirstBinAbove(0.0);
 	Int_t sim_hi_end_bin = h_sim_tt->FindLastBinAbove(0.0);

	Int_t hi_bin;
	Double_t hi_limit;

	Int_t bin_itt = sim_hi_end_bin;	
	
	Double_t sim_yield_all =  h_sim_tt->Integral(sim_lo_end_bin, sim_hi_end_bin);

	while( h_sim_tt->Integral(bin_itt, sim_hi_end_bin) < 0.02*sim_yield_all) {
		bin_itt--;
	}

	cout << "Sim Lo Bin: "<< sim_lo_end_bin << endl;
	cout << "Sim Hi Bin: "<< sim_hi_end_bin << endl;
	cout << "10 percent from Hi: "<< bin_itt << endl;
	
	hi_bin = bin_itt;
 
 	hi_limit = h_sim_tt->GetBinCenter(hi_bin); 
 
//	exit(0);
 		
	return hi_limit;

}



/*--------------------------------------------------*/
/// Find integration range low limit
Double_t Find_int_lo_refit(TPad* tar_pad) {

  	TH1F* h_sim_tt = (TH1F*) tar_pad->GetPrimitive("u_phi_sim_omega");
	
 	Int_t sim_lo_end_bin = h_sim_tt->FindFirstBinAbove(0.0);
 	Int_t sim_hi_end_bin = h_sim_tt->FindLastBinAbove(0.0);

	Int_t low_bin;
	Double_t low_limit;

	Int_t bin_itt=sim_lo_end_bin;	
	
	Double_t sim_yield_all =  h_sim_tt->Integral(sim_lo_end_bin, sim_hi_end_bin);

	while( h_sim_tt->Integral(sim_lo_end_bin, bin_itt) < 0.04*sim_yield_all) {
		bin_itt++;
	}

	cout << "Sim Lo Bin: "<< sim_lo_end_bin << endl;
	cout << "Sim Hi Bin: "<< sim_hi_end_bin << endl;
	cout << "10 percent from Lo: "<< bin_itt << endl;

	
	low_bin = bin_itt;
 
 	low_limit = h_sim_tt->GetBinCenter(low_bin); 
 
//	exit(0);
 		
	return low_limit;

}


///*--------------------------------------------------*/
/// Find integration range high limit
Double_t Find_int_hi_refit(TPad* tar_pad) {

  	TH1F* h_sim_tt = (TH1F*) tar_pad->GetPrimitive("u_phi_sim_omega");
	
 	Int_t sim_lo_end_bin = h_sim_tt->FindFirstBinAbove(0.0);
 	Int_t sim_hi_end_bin = h_sim_tt->FindLastBinAbove(0.0);

	Int_t hi_bin;
	Double_t hi_limit;

	Int_t bin_itt = sim_hi_end_bin;	
	
	Double_t sim_yield_all =  h_sim_tt->Integral(sim_lo_end_bin, sim_hi_end_bin);

	while( h_sim_tt->Integral(bin_itt, sim_hi_end_bin) < 0.04*sim_yield_all) {
		bin_itt--;
	}

	cout << "Sim Lo Bin: "<< sim_lo_end_bin << endl;
	cout << "Sim Hi Bin: "<< sim_hi_end_bin << endl;
	cout << "10 percent from Hi: "<< bin_itt << endl;
	
	hi_bin = bin_itt;
 
 	hi_limit = h_sim_tt->GetBinCenter(hi_bin); 
 
//	exit(0);
 		
	return hi_limit;

}




/*--------------------------------------------------*/
/*--------------------------------------------------*/


bool include_etap(TPad* tar_pad) {

	bool inc_etap;
	
  	TH1F* h_tt = (TH1F*) tar_pad->GetPrimitive("u_phi_target");
	
 	Int_t lo_end_bin = h_tt->FindFirstBinAbove(0.0);
 	Int_t hi_end_bin = h_tt->FindLastBinAbove(0.0);

// 	Int_t hi_bin;
// 	Double_t hi_limit;
// 
// 	Int_t bin_itt = hi_end_bin;	
// 	
// 	Double_t yield_all =  h_tt->Integral(lo_end_bin, hi_end_bin);
// 
// 	while( h_tt->Integral(bin_itt, hi_end_bin) < 0.02*yield_all) {
// 		bin_itt--;
// 	}
// 
// 	cout << "Sim Lo Bin: "<< lo_end_bin << endl;
// 	cout << "Sim Hi Bin: "<< hi_end_bin << endl;
// 	cout << "10 percent from Hi: "<< bin_itt << endl;
// 	
// 	hi_bin = bin_itt;
//  
//	hi_limit = h_tt->GetBinCenter(hi_bin); 
 

	

	if (h_tt->Integral(h_tt->FindBin(0.94), hi_end_bin) > 0.1 * h_tt->Integral()) {

//		cout << h_tt->Integral(h_tt->FindBin(0.9), hi_end_bin) << "   " <<  0.1 * h_tt->Integral() << endl;

		inc_etap = true;

//		exit(0);

	} else {

		inc_etap = false;

	}

	return inc_etap;

}







