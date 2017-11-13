
void recon_kin(TCanvas* cc3, TString kin_name,TString file_name) {

    TString out_int_dir_str;

	out_int_dir_str = "recon_kin/";
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


	const Int_t u_bin_num_const   = u_bin_num;  
	const Int_t phi_bin_num_const = phi_bin_num;  
	const Int_t bin_num_const     = u_bin_num_const * phi_bin_num_const;

	TCanvas* cc4 = new TCanvas("cc4", "cc4", 1600, 800);
	cc4->Divide(phi_bin_num, u_bin_num, 0.003, 0.003);
 
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
  		TH1F* h_sim_sum = (TH1F*) h_sim_tmp->Clone("h_sim_clone");

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

		
		h_sim_rho[in]->Scale(rho_f_vec[in]);
        h_sim_xphsp[in]->Scale(xphsp_f_vec[in]);
        h_sim_eta[in]->Scale(eta_f_vec[in]);
        h_sim_etap[in]->Scale(etap_f_vec[in]);

		h_sim_sum->Add(h_sim_rho[in],   1);
		h_sim_sum->Add(h_sim_xphsp[in], 1);
 		h_sim_sum->Add(h_sim_eta[in],   1);
 		h_sim_sum->Add(h_sim_etap[in],  1);


		h_sim_sum->SetLineColor(6); 
		h_sim[in]->SetLineColor(2);
		h_sim_xphsp[in]->SetLineColor(3);
		h_sim_rho[in]->SetLineColor(4);
		h_sim_eta[in]->SetLineColor(1);
		h_sim_etap[in]->SetLineColor(1);


  		h_tgt[in]->DrawCopy();
		h_sim_sum->DrawCopy("hsame");
//		h_tgt[in]->DrawCopy("same");

		h_sim[in]->DrawCopy("hsame");
		h_sim_rho[in]->DrawCopy("hsame"); 
		h_sim_xphsp[in]->DrawCopy("hsame");
		h_sim_eta[in]->DrawCopy("hsame"); 
		h_sim_etap[in]->DrawCopy("hsame"); 
 
 		Double_t intg_lo = lo_limit_vec[in];
 		Double_t intg_hi = hi_limit_vec[in];
 
 
 		if (intg_lo != 0.0 && intg_hi != 0.0) {
 
//  			TCanvas* cc5 = new TCanvas("cc5", "cc5", 800, 600);
//  			cc5->cd();
//  //			pnew->SetPad(0.01, 0.01, 0.01, 0.01);
//  
//  			TPad *pnew = (TPad*) cur_pad->Clone();
//  
//  			pnew->SetPad(cc5->GetX1(),cc5->GetY1(),cc5->GetX2(),cc5->GetY2());
//  
//  			cout << out_int_dir_str + "individual_plot/" + f_str_temp + ".png" << endl;
//  			
//  			pnew->Draw();
//  
//  			cc5->Print(out_int_dir_str + "individual_plot/" + f_str_temp + ".png");
//  
//  //			exit(0);
//  			
//  			delete cc5;
 
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
// 
// 		cc4->cd();
// 
//


 	}
 
 	cc4->Print(out_int_dir_str + "u_fit_itt_" + kin_name + "_"+ file_name+".png");
 	cc4->Print(out_int_dir_str + "u_fit_itt_" + kin_name + "_"+ file_name+".root");
// 
// //	cout << "asdasdasdasdas  asdasd <<<  aEnd Loop" << endl;
// 
// 
// 	delete cc4;
// 
// //	exit(0);
 

//	exit(0);

}

