TGraphErrors* q2_g;
TGraphErrors* q2_ratio_lo;
TGraphErrors* q2_ratio_hi;


TGraphErrors* q2_u_set;

TMultiGraph *mg;

// const Int_t u_bin_num = 2;
// Float_t u_lower_limit[u_bin_num];
// Float_t u_upper_limit[u_bin_num];

Int_t u_bin_num;
Double_t* u_lower_limit;
Double_t* u_upper_limit;

TString outplot_dir;
outplot_dir = "ratio_check/";

TString ave_dir;
ave_dir = "averages/";

Int_t iii = 0;

void plot_R() {

	gROOT->ProcessLine(".L get_u_phi_bin.C");

	u_bin_num = Get_u_bin();

	const Int_t u_bin_num_const = u_bin_num;

	u_lower_limit = new Double_t[u_bin_num_const];
	u_upper_limit = new Double_t[u_bin_num_const];


//	Double file_name_arr[4] = {"160_32", "160_59", "245_27", "245_55"};
//	Double file_name_arr[4] = {"160_32", "160_59", "245_27", "245_55"};

//	file_name = "aver_160_32.dat";

// 	file_name = "aver_160_32.dat";
// 	file_name = "aver_160_32.dat";
// 	file_name = "aver_160_32.dat";
// 	file_name = "aver_160_32.dat";
// 

//	file_name = {"160_32", "160_59", "245_27", "245_55"};


//	TString file_name_str;

 	q2_g   = new TGraphErrors();
 	q2_ratio_lo   = new TGraphErrors();
 	q2_ratio_hi   = new TGraphErrors();

 	TF1* q2_fit = new TF1("q2_fit", "pol1", 150, 260);

	mg = new TMultiGraph();
	
 	q2_u_set  = new TGraphErrors();

// 
// 	for (Int_t i; i < 4; i++) {
// 	
// 
// //		Single_Setting("pl_" + file_name_arr[i]);
// //		file_name_str =  "aver_" + file_name[i] + ".dat";
// 
// 
// 
// 	}

	Single_Setting(160, 32);
	Single_Setting(160, 59);
	Single_Setting(245, 27);
	Single_Setting(245, 55);



	TCanvas* c2 = new TCanvas();


	c2->cd();

 	q2_g->SetTitle("u Fit vs Q{^2} ");
 	
    q2_g->GetXaxis()->SetTitle("Q^{2}");		
    q2_g->GetXaxis()->CenterTitle();		

    q2_g->GetYaxis()->SetTitle("Fitted parameter");		
    q2_g->GetYaxis()->SetTitleOffset(1.2);		
    q2_g->GetYaxis()->CenterTitle();		

	q2_g->Draw("A*");
	
	
	q2_g->Fit("q2_fit");

	c2->Print(outplot_dir+"q2_setting.png");


	c2->Clear();
	c2->cd();


 	q2_u_set->SetTitle("Q^{2} vs -u");
 	
    q2_u_set->GetXaxis()->SetTitle("-#it{u} [GeV^{2}]");		
    q2_u_set->GetXaxis()->CenterTitle();		

    q2_u_set->GetYaxis()->SetTitle("Q^{2}");		
    q2_u_set->GetYaxis()->SetTitleOffset(1.2);		
    q2_u_set->GetYaxis()->CenterTitle();		

	q2_u_set->Draw("A*");

	c2->Print(outplot_dir+"q2_vs_u.png");



	c2->Clear();
	c2->cd();


	q2_ratio_hi->SetMarkerColor(2);
	q2_ratio_hi->Draw("A*");

	q2_ratio_hi->SetTitle("Average Yield Ratio (Yexp/Ysim) vs Q^{2}");
	
    q2_ratio_hi->GetXaxis()->SetTitle("Q^{2}");		
    q2_ratio_hi->GetXaxis()->CenterTitle();		
	
    q2_ratio_hi->GetYaxis()->SetTitle("Average Yield Ratio");		
    q2_ratio_hi->GetYaxis()->SetTitleOffset(1.2);		
    q2_ratio_hi->GetYaxis()->CenterTitle();		


	if(TMath::MaxElement(q2_ratio_lo->GetN(),q2_ratio_lo->GetY()) > TMath::MaxElement(q2_ratio_hi->GetN(),q2_ratio_hi->GetY())) {

		q2_ratio_hi->SetMaximum(TMath::MaxElement(q2_ratio_lo->GetN(),q2_ratio_lo->GetY())*1.1);

	}

	q2_ratio_lo->Draw("*");

	c2->Print(outplot_dir + "q2_ratio.png");
	


	c2->Clear();
	c2->cd();



	TF1* u_fit_total = new TF1("u_total", "pol1", 0.04, 0.6); 



	mg->SetTitle("Average Yield Ratio (Yexp/Ysim) vs -u");
	
	mg->Draw("a");

    mg->GetXaxis()->SetTitle("-#it{u} [GeV^{2}]");		
    mg->GetXaxis()->CenterTitle();		
	
    mg->GetYaxis()->SetTitle("Average Yield Ratio");		
    mg->GetYaxis()->SetTitleOffset(1.2);		
    mg->GetYaxis()->CenterTitle();		

	mg->Fit("u_total", "R");


	TLatex* tt = new TLatex();

	tt->SetNDC(true);
	tt->SetTextSize(0.035);

	tt->SetText(0.15, 0.33, "low Q^{2}, low #epsilon");
	tt->SetTextColor(1);
	tt->DrawClone("same");

	tt->SetText(0.15, 0.28, "low Q^{2}, high #epsilon");
	tt->SetTextColor(2);
	tt->DrawClone("same");

	tt->SetText(0.15, 0.23, "high Q^{2}, low #epsilon");
	tt->SetTextColor(3);
	tt->DrawClone("same");

	tt->SetText(0.15, 0.18, "high Q^{2}, high #epsilon");
	tt->SetTextColor(4);
	tt->DrawClone("same");

	c2->Print( outplot_dir + "total_yield_ratio.png");


}




void Single_Setting(Int_t q2_set, Int_t eps_set) {

	const Int_t u_bin_num_const = u_bin_num;

	ifstream u_bin_file("../u_bin_interval", std::fstream::in); 

	if (!u_bin_file.is_open()) {

		cout << "u file doesn't exits! " << endl; 
		exit(0)
	}

	string dummy_str;

	getline(u_bin_file, dummy_str);

 	if( q2_set == 160) {
 
 	/// t from 0.01 to 0.45
 
 	///	u_lower_limit[0] = 0.01;

	//	u_bin_file >> >> >> ; 

		Float_t temp_num[u_bin_num_const];

		for(Int_t i = 0; i <= u_bin_num; i++) {
			u_bin_file >> temp_num[i];
			cout << temp_num[i] << endl;
		}


		for(Int_t i = 0; i < u_bin_num; i++) {

			u_lower_limit[i] = temp_num[i];
			u_upper_limit[i] = temp_num[i+1];

			cout << 	u_lower_limit[i] << "   " << u_upper_limit[i] << endl;

		}


//		exit(0);

//  		u_lower_limit[0] = 0.0;
//  		u_upper_limit[0] = 0.12;
//  
//  		u_lower_limit[1] = 0.12;
//  		u_upper_limit[1] = 0.2;
//  
//  		u_lower_limit[2] = 0.2;
//  		u_upper_limit[2] = 0.4;
 
 	} else if ( q2_set == 245) {
 
		getline(u_bin_file, dummy_str);
		getline(u_bin_file, dummy_str);

		Float_t temp_num[u_bin_num_const];

		for(Int_t i = 0; i <= u_bin_num; i++) {
			u_bin_file >> temp_num[i];
			cout << temp_num[i] << endl;
		}


		for(Int_t i = 0; i < u_bin_num; i++) {

			u_lower_limit[i] = temp_num[i];
			u_upper_limit[i] = temp_num[i+1];

			cout << 	u_lower_limit[i] << "   " << u_upper_limit[i] << endl;

		}

//  		u_lower_limit[0] = 0.0; 
//  		u_upper_limit[0] = 0.212;
//  
//  		u_lower_limit[1] = 0.212;
//  		u_upper_limit[1] = 0.33;
//  
//  		u_lower_limit[2] = 0.33;
//  		u_upper_limit[2] = 0.60;

//			exit(0);
 
 	}




	Double_t t_set;

//	cout << u_lower_limit[2] << endl;
	
//	exit(0);


	TString file_name_str;
	
	file_name_str.Form("pl_%i_%i", q2_set, eps_set);

	file_name =  ave_dir + "aver." + file_name_str + ".dat";

 	TNtuple* n1 = new TNtuple("n1", "n1", "ratio:r_err:phi_bin:t_bin:y_exp:y_sim");
 
 	n1->ReadFile(file_name);
 
    TGraphErrors* u_g = new TGraphErrors();

	TF1* t_fit = new TF1("t1", "pol1", 0, 4); 

	TCanvas *c1 = new TCanvas("n1", "n1", 1200, 400);
	
	c1->Divide(3,1);

	Float_t q2_sum;




	for (Int_t i = 0; i < u_bin_num; i++) {

		Int_t in = i + 1;	

//		cout << in << endl;

		c1->cd(in);

		TString t_bin_limit;
		t_bin_limit.Form( "t_bin == %i && ratio != 0", in);	

		// n1->Draw("ratio:phi_bin", t_bin_limit, "*");
		n1->Draw("ratio:phi_bin:r_err", t_bin_limit, "off");

    	TGraphErrors* g = new TGraphErrors(n1->GetSelectedRows(), n1->GetV2(), n1->GetV1(), 0, n1->GetV3());


		

		g->SetTitle("Omega Yield Ratio: Yexp/Ysim");
	
        g->GetXaxis()->SetTitle("#phi angle");		
        g->GetXaxis()->CenterTitle();		
	
        g->GetYaxis()->SetTitle("Yield Ratio");		
        g->GetYaxis()->SetTitleOffset(1.8);		
        g->GetYaxis()->CenterTitle();		


		Double_t xx;
		Double_t yy;

		Double_t sum_r = 0;
		Double_t sum_w = 0;

		Double_t weighted_ave_top = 0;
		Double_t weighted_ave_bot = 0;



		for( Int_t ii = 0; ii < g->GetN(); ii++) {

			g->GetPoint(ii, xx, yy);		

			if(yy==0.0) {
				g->RemovePoint(ii);
			}
		}
		


		for( Int_t ii = 0; ii < g->GetN(); ii++) {
			
			g->GetPoint(ii, xx, yy);		
	
			cout << ii  << "    " << xx << "    " << yy << "     " << g->GetErrorY(ii) << endl;

			sum_r =  sum_r + yy;
			sum_w =  sum_w + 1/(g->GetErrorY(ii)**2);

			weighted_ave_top = weighted_ave_top + 1/(g->GetErrorY(ii)**2) * yy;
			weighted_ave_bot = weighted_ave_bot + 1/(g->GetErrorY(ii)**2);

		}

 		g->Draw("AP");


		Int_t graph_it = u_g->GetN(); 

		t_set = (u_lower_limit[i] + u_upper_limit[i])/2;

//		u_g->SetPoint(graph_it, i+1, sum_r/4);

//		u_g->SetPoint(graph_it, t_set, sum_r/graph_it);
//		u_g->SetPointError(graph_it, 0, sqrt(sum_r_err)/graph_it);


		cout << "ppppppppppppppppp   "<< g->GetN() << endl;

//		exit(0);
			

//		u_g->SetPoint(graph_it, t_set, sum_r/g->GetN());
//		u_g->SetPointError(graph_it, 0, sqrt(sum_r_err)/g->GetN());

		u_g->SetPoint(graph_it, t_set, weighted_ave_top/weighted_ave_bot);
		u_g->SetPointError(graph_it, 0, 1/sqrt(sum_w));


		Float_t q2_set_tmp;

//		q2_set_tmp = Get_u_bin_q2(q2_set, i);

		q2_set_tmp = Get_u_bin_q2(q2_set, i);


		q2_u_set->SetPoint(q2_u_set->GetN(), t_set, q2_set_tmp);

		
		cout << t_set << "   " << q2_set_tmp << endl;
//		exit(0);


// 		if (q2_set == 160) {		
// 
// 			if (i == 0) { 
// 	
// 				q2_set_tmp = 150.458;
// 			
// 			} else if (i == 1) { 
// 	
// 				q2_set_tmp = 162.791;
// 	
// 			} else {
// 	
// 				q2_set_tmp = 170.599;
// 	
// 			}
// 
// 		} else{
// 
// 			if (i == 0) { 
// 	
// 				q2_set_tmp = 225.943;
// 			
// 			} else if (i == 1) {
// 	
// 				q2_set_tmp = 243.919;
// 	
// 			} else {
// 	
// 				q2_set_tmp = 261.004;
// 	
// 			}
// 
// 		}





		if (eps_set < 40) { 

			Int_t graph_it = q2_ratio_lo->GetN(); 
//			q2_ratio_lo->SetPoint(graph_it, q2_set_tmp, sum_r/graph_it);
//			q2_ratio_lo->SetPoint(graph_it, q2_set_tmp, sum_r/g->GetN());

			q2_ratio_lo->SetPoint(graph_it, q2_set_tmp, weighted_ave_top/weighted_ave_bot);

		} else { 

			Int_t graph_it = q2_ratio_hi->GetN(); 
//			q2_ratio_hi->SetPoint(graph_it, q2_set_tmp, sum_r/graph_it);
//			q2_ratio_hi->SetPoint(graph_it, q2_set_tmp, sum_r/g->GetN());

			q2_ratio_hi->SetPoint(graph_it, q2_set_tmp, weighted_ave_top/weighted_ave_bot);

		}

		q2_sum =  q2_sum + q2_set_tmp;


	}


//	exit(0);



	c1->Print(outplot_dir + "ratio_check_u_phi_bin" + file_name_str + ".png");
	c1->Print(outplot_dir + "ratio_check_u_phi_bin" + file_name_str + ".root");

	c1->Clear();



	c1->cd();
	

    u_g->SetTitle("u Dependence Plot");
     
    u_g->GetXaxis()->SetTitle("-#it{u} [GeV^{2}]");		
    u_g->GetXaxis()->CenterTitle();		
     
    u_g->GetYaxis()->SetTitle("Yield Ratio");		
    u_g->GetYaxis()->SetTitleOffset(1.8);		
    u_g->GetYaxis()->CenterTitle();		

	u_g->Draw("a*");

	TGraphErrors* u_g_1 = (TGraphErrors*) u_g->Clone();
	
	iii = iii+1;

	u_g_1->SetMarkerColor(iii);
	u_g_1->SetLineColor(iii);

    mg->Add(u_g_1, "*");



	
	gStyle->SetOptFit(); 

	u_g->Fit("t1");



	c1->Print(outplot_dir + "average_ratio" + file_name_str + ".png");


	Int_t graph_it = q2_g->GetN();

//	q2_g->SetPoint(graph_it, q2_set, t_fit->GetParameter(1));
//
	

	cout << "!!!!!!   " <<  q2_g->GetN() << "   " << q2_sum/3 << "   " << t_fit->GetParameter(1) << endl ;

	q2_g->SetPoint(graph_it, q2_sum/3, t_fit->GetParameter(1));

//	q2_g->SetPointError(graph_it, 0, t_fit->GetParError(1));

//	exit(0);

	delete c1;

}



// Double_t t_dependent_fitting_function(Double_t x, Double par) {
// 
// 
// }
// 
// 








/*--------------------------------------------------*/
/*--------------------------------------------------*/

Double_t Get_u_bin_q2(Int_t q2_set, Int_t t_bin_num) {
	
	TString file_str;
	TString dir_str;

	dir_str = "averages/";

	file_str.Form("avek.%i.dat",q2_set);

	cout << file_str << endl;

	TNtuple* avek_ntp = new TNtuple("avek", "avek", "w:dw:q2:dq2:theta_lo:theta_hi:t_bin_in");

	avek_ntp->ReadFile(dir_str + file_str);
	
	Float_t q2_ave;
	
	avek_ntp->SetBranchAddress("q2", &q2_ave);	
	
	avek_ntp->GetEntry(t_bin_num);

	cout << q2_ave << endl;


	delete avek_ntp;

	return q2_ave; 


}




/*--------------------------------------------------*/
/*--------------------------------------------------*/

