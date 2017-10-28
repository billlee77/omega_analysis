#include <iostream>

#include "analysis.h"
//#include "cut.h"

using namespace std;


/*--------------------------------------------------*/
/*--------------------------------------------------*/

Analysis::Analysis() {

// 	eff_file = "list.settings.fpi2";
// 	off_file = "fit_piplus/cointime_pl.dat";
// 
// // 	eff_file = "list.settings.heep";
// // 	off_file = "offset.dat";
// // 
// 	Init();

}


/*--------------------------------------------------*/
/*--------------------------------------------------*/

Analysis::Analysis(ReadFile::efficiency_profile eff_struc) {

	eff_file = "list.settings.fpi2";
	off_file = "fit_piplus/cointime_pl.dat";

// 	eff_file = "list.settings.heep";
// 	off_file = "offset.dat";
// 
	eff_ana = eff_struc;

	Init();

}


/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Analysis::Init() {


//	ReadFile* rff = new ReadFile("list.settings.fpi2", "fit_piplus/cointime_pl.dat");
	ReadFile* rff = new ReadFile(eff_file, off_file);

	kin_ana = rff->kin_pro;
//	eff_ana = rff->eff_pro1;
	cen_ana = rff->cen_pro;

	cen_runs = rff->Get_Cen_Runs();

	list_file = rff->Get_List_File();


	data_file_dir = "data/";

	Para_Init(); 


	tree_out = Ceate_File();

	c1 = new TCanvas();
	c2 = new TCanvas();
	c3 = new TCanvas();
	c4 = new TCanvas();



	c2->Divide(1,2);


	phi_no_sub = new TCanvas();
	phi_no_sub->Divide(2,3);

	phi_sub = new TCanvas();
	phi_sub->Divide(2,3);


// 	real_event[0] = new TH1F("phi_real1", "", phi_bin_num, 0, 360);
// 	real_event[1] = new TH1F("phi_real2", "", phi_bin_num, 0, 360);
// 	real_event[2] = new TH1F("phi_real3", "", phi_bin_num, 0, 360);
// 	real_event[3] = new TH1F("phi_real4", "", phi_bin_num, 0, 360);
// 	real_event[4] = new TH1F("phi_real5", "", phi_bin_num, 0, 360);
// 	real_event[5] = new TH1F("phi_real6", "", phi_bin_num, 0, 360);

	real_event[0] = new TH1F("", "real_event_1", phi_bin_num, 0, 360);
	real_event[1] = new TH1F("", "real_event_2", phi_bin_num, 0, 360);
	real_event[2] = new TH1F("", "real_event_3", phi_bin_num, 0, 360);
	real_event[3] = new TH1F("", "real_event_4", phi_bin_num, 0, 360);
	real_event[4] = new TH1F("", "real_event_5", phi_bin_num, 0, 360);
	real_event[5] = new TH1F("", "real_event_6", phi_bin_num, 0, 360);

	chain = new TChain("chain");

	list = new TList;

	is_run_dummy = false;




	diamond_setting     = new TMultiGraph();
	diamond_setting_cut = new TMultiGraph();

	mm_offset_setting = new TH1F( "", "mm", 50, -0.1, 0.12);
	mm_setting        = new TH1F( "", "mm", 50, -0.1, 0.12);





//	exit(0);


}





/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Analysis::Run_by_Run_Analysis(Int_t run_itt){

	is_run_dummy = false;

	cout << " sadadout_dirout_dir " << run_itt << endl;

	// cout << "Is the run dummy? " << is_run_dummy << endl;
	
	run_tree = Ceate_File();

	Para_Run_Def(run_itt);

	Set_Coin_Min(-6);
	Set_Coin_Max(10);
	Set_Coin_Bin(80);

	Set_Expected_MM(0.939565);

	Double_t array_temp[6] = {0.0, 0.00574, 0.00554, 0.00254, 0.00253, 0.00260};

	Set_MM_Offset(array_temp);
	
	Load_Data_Tree();

	Define_Cuts();
	
	/*--------------------------------------------------*/
	/// From the paw
	//
	//    tit1='missmass' ; nb1=50  ; xlo1=-0.075   ; xhi1=0.125
	//    tit2='W'        ; nb2=30  ; xlo2=2.00     ; xhi2=2.40
	//    tit3='t'        ; nb3=18  ; xlo3=[tmin]   ; xhi3=[tmax]
	//    tit4='Q2'       ; nb4=30  ; xlo4=0.5*[Q2] ; xhi4=1.5*[Q2]
	//    tit5='th_pq'    ; nb5=30  ; xlo5=0.       ; xhi5=10.
	//    tit6='phi_pq'   ; nb6=36  ; xlo6=0.       ; xhi6=360.
	//    tit7='PmPar'    ; nb7=30  ; xlo7=[parmn]  ; xhi7=[parmx]
	//    tit8='PmPer'    ; nb8=30  ; xlo8=-0.32    ; xhi8=0.32
	//    tit9='PmOop'    ; nb9=30  ; xlo9=-0.30    ; xhi9=0.30
	// 
	//    tit10='hsdelta' ; nb10=40 ; xlo10=-10.    ; xhi10=10.
	//    tit11='hsytar'  ; nb11=50 ; xlo11=-2.5    ; xhi11=2.5
	//    tit12='hsyptar' ; nb12=50 ; xlo12=-0.05   ; xhi12=0.05
	//    tit13='hsxptar' ; nb13=50 ; xlo13=-0.10   ; xhi13=0.10
	//    tit14='ssdelta' ; nb14=80 ; xlo14=-20.    ; xhi14=20.
	//    tit15='ssytar'  ; nb15=50 ; xlo15=-2.5    ; xhi15=2.5
	//    tit16='ssyptar' ; nb16=50 ; xlo16=-0.10   ; xhi16=.10
	//    tit17='ssxptar' ; nb17=50 ; xlo17=-0.05   ; xhi17=0.05
 

//	data_tree_in->Draw("missmass >> hist_missmass","missmass > -0.075 && missmass < 0.125", "goff");
// 	Q2_limit.Form("Q2 < %f && Q2 > %f", Q_2*1.5, Q_2*0.5);
// 	data_tree_in->Draw("Q2 >> hist_q2", Q2_limit, "goff");

//	cout << all_coin_cut << endl;


	Missing_Mass_Plot();

	TH1F* phi_real_check[t_bin_set];
	TH1F* phi_rand_check[t_bin_set];

	TH1F* phi_sub_check[t_bin_set];

	// cout << t_bin_set << endl;


// 	for(int t_bin_set_tmp = 1; t_bin_set_tmp <= t_bin_set; t_bin_set_tmp++) {
// 
// 		TString check_name;
// 		check_name.Form("phi_check_%i", t_bin_set_tmp);
// 
// 		cout << check_name << endl;
// 		phi_real_check[t_bin_set_tmp-1] = new TH1F(check_name, check_name, phi_bin_num, 0, 360);
// 
// 
// 
// 
// 	}
// 
// 
// 	phi_real_check[0] ->Fill(10);
// 	phi_real_check[0] ->Fill(20);
// 	phi_real_check[0] ->Fill(30);
// 	phi_real_check[0] ->Fill(60);
// 
// 
// 
// 	phi_real_check[4] ->Fill(10);
// 	phi_real_check[4] ->Fill(20);
// 	phi_real_check[4] ->Fill(30);
// 	phi_real_check[4] ->Fill(60);
// 


// 	for(int t_bin_set_tmp = 1; t_bin_set_tmp <= t_bin_set; t_bin_set_tmp++) {
// 
// 
// 
// 
// 		phi_no_sub->cd(t_bin_set_tmp);
// 
// 		phi_real_check[t_bin_set_tmp-1]->Draw("PE");
// 		
// 
// 		phi_no_sub->Update();
// 
// 
// 	}



// 	phi_no_sub->cd();
// 
// 	phi_real_check[4]->Draw("hist");
// 
// 
// 	phi_no_sub->Update();
// 

//	out_dir->cd();
//	phi_no_sub->Write("testtest");



//exit(0);


	for(int t_bin_set_tmp = 1; t_bin_set_tmp <= t_bin_set; t_bin_set_tmp++) {
		
		// float t_c = t_min + (t_bin_set_tmp-0.5) * t_width;
		float yield_cut_t_l = t_min + (t_max - t_min)/t_bin_set * (t_bin_set_tmp-1);
		float yield_cut_t_h = t_min + (t_max - t_min)/t_bin_set * t_bin_set_tmp;

		yield_t_bin_cut.Form("t > %f && t < %f",yield_cut_t_l, yield_cut_t_h);

		TH1F* event_phi_real = new TH1F("phi_real", "phi_real", phi_bin_num, 0, 360);
		data_tree_in->Draw("phi_pq*180/3.1415926 >> phi_real", yield_t_bin_cut + "&&" + real_coin_cut, "goff");

		TH1F* event_phi_rand = new TH1F("phi_rand", "phi_rand", phi_bin_num, 0, 360);
		data_tree_in->Draw("phi_pq*180/3.1415926 >> phi_rand", yield_t_bin_cut + "&&" + rand_coin_cut, "goff");

		TH1F* event_phi_real_clone = (TH1F*) event_phi_real->Clone();
		TH1F* event_phi_rand_clone = (TH1F*) event_phi_rand->Clone();



//  		event_phi_real->SetOption("PE");
// 		event_phi_real->Scale(0.463);		
// 
// 		event_phi_real->Write();


		/*--------------------------------------------------*/
		/*--------------------------------------------------*/
		/// Note that the scale factor of 0.463 is documented in the report by Dave Meekins,  
		/// the scale factor corrects the difference in the target thickness difference between 
		///	dummy target and real target 

		if (is_run_dummy) {
			event_phi_real->Scale(0.463);		
			event_phi_rand->Scale(0.463);		
		}



		phi_real_check[t_bin_set_tmp-1] = (TH1F*) event_phi_real->Clone(); 

		event_phi_real->Add(event_phi_rand, -0.3333333333);

		event_phi_rand->Scale(0.3333333333);

		phi_sub_check[t_bin_set_tmp-1] = (TH1F*) event_phi_real->Clone(); 
		phi_rand_check[t_bin_set_tmp-1] =  (TH1F*) event_phi_rand->Clone(); 


		
// 		out_dir->cd();
// 		event_phi_real->Draw("PE");
// 		event_phi_rand->Draw("PE");
// 
//		event_phi_real->Write();
	
//		event_phi_rand->Write();

// 		if (is_run_dummy) {
// 			event_phi_real->Add(event_phi_rand, -0.3333333333);
// 		} else {
// 			event_phi_real->Add(event_phi_rand, 1);
// 		}
// 

		real_event[t_bin_set_tmp-1]->Add(event_phi_real);	

		for (int iiii = 0;  iiii < phi_bin_num; iiii++) {

//			int real_err_itt = event_phi_real_clone->GetBinContent(iiii+1); 
//			int rand_err_itt = event_phi_rand_clone->GetBinContent(iiii+1);
//			int real_event_err_itt = real_err_itt ;

			if (is_run_dummy) {
				rand_err = event_phi_rand_clone->GetBinError(iiii+1) * 0.463 * 0.3333333333;
			} else { 
//				rand_err = event_phi_rand_clone->GetBinError(iiii+1) * 0.463;
				rand_err = event_phi_rand_clone->GetBinError(iiii+1) *0.3333333333;
			}


 			if (is_run_dummy) {
//				rand_err = event_phi_rand_clone->GetBinError(iiii+1) * 0.463 * 0.3333333333;
				real_err = event_phi_real_clone->GetBinError(iiii+1) * 0.463;
			} else { 
//				rand_err = event_phi_rand_clone->GetBinError(iiii+1) * 0.463;
				real_err = event_phi_real_clone->GetBinError(iiii+1);
			}

			real_event_error = sqrt( pow(real_err,2) + pow(rand_err,2) );		

//			yield_err_sum[t_bin_set_tmp-1][iiii] = yield_err_sum[t_bin_set_tmp-1][iiii] + pow(real_event_error, 2);
			yield_err_sum[t_bin_set_tmp-1][iiii] = pow(real_event_error, 2);

			event_phi_real->SetBinError(iiii+1, real_err);

			phi_real_check[t_bin_set_tmp-1]->SetBinError(iiii+1, real_err);

			phi_rand_check[t_bin_set_tmp-1]->SetBinError(iiii+1, rand_err);

			phi_sub_check[t_bin_set_tmp-1]->SetBinError(iiii+1, real_event_error);


//			cout << "     " << iiii << "    "  <<  real_event_error << endl;

		}

		TString check_name;
		check_name.Form("%i", t_bin_set_tmp);
		
		phi_real_check[t_bin_set_tmp-1] ->SetName("phi_check_real_"+check_name); 
		phi_rand_check[t_bin_set_tmp-1] ->SetName("phi_check_rand_"+check_name); 
		phi_sub_check[t_bin_set_tmp-1] ->SetName("phi_check_sub_"+check_name); 

		phi_no_sub->cd(t_bin_set_tmp);
		phi_real_check[t_bin_set_tmp-1]->Draw("PE");
		phi_rand_check[t_bin_set_tmp-1]->SetLineColor(2);
		phi_rand_check[t_bin_set_tmp-1]->Draw("same");
		phi_no_sub->Update();

 		phi_sub->cd(t_bin_set_tmp);
 		phi_sub_check[t_bin_set_tmp-1]->Draw("PE");
 		phi_sub->Update();

		delete event_phi_real;
		delete event_phi_rand;

		delete event_phi_real_clone;
		delete event_phi_rand_clone;

	}

	
	out_dir->cd();
 	phi_no_sub->Write("phi_check_no_sub");
 	phi_sub->Write("phi_check_sub");

	
	Print_out_data();
	Para_Run_Clean();

}



/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Analysis::Para_Init() {

	Q_2		= kin_ana.Q2set[0];
	kset	= kin_ana.kset[0];
	t_min	= kin_ana.tmnset[0];
	t_max	= kin_ana.tmxset[0];

	t_bin_set = kin_ana.NBtset[0];
	t_width   = ( t_max - t_min ) / t_bin_set;


	phi_bin_num = 16;

	pstp = 360./phi_bin_num;
	pmn  = pstp/2.;
	pmx  = 360. - pstp/2.;

	acccc_temp = 0.0;
	errrr_temp = 0.0;
	charge_tot = 0.0;


	yield_vec.resize(t_bin_set*phi_bin_num);
	yield_err_vec.resize(t_bin_set*phi_bin_num);
	phi_vec.resize(t_bin_set*phi_bin_num);
	tb_vec.resize(t_bin_set*phi_bin_num);

	// cout << kset << endl;

}




/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Analysis::Para_Run_Def(Int_t num) {

//	cout << "aaaaaaaaaaaaaaaasssssssbbbbbbbbbbbbb" << endl;

	run_num     = eff_ana.run_num[num];

	charge      = eff_ana.charge[num];

	tot_eff     = eff_ana.tot_eff[num];

	charge_err  = eff_ana.charge_err[num];

	tot_eff_err = eff_ana.tot_eff_err[num];

	q2_run      = Float_t(eff_ana.Q2[num]);

	ebeam_run   = Float_t(eff_ana.ebeam[num]);


	int pos = -100000000;

 	for(int ii = 0; ii < cen_runs; ii++) {
 		if (cen_ana.center_run_num[ii] == run_num) {
 			pos = ii;
 		}
 	}


	// cout << is_run_dummy << endl;

	if (is_run_dummy) {
		target = "Dummy"; 
	} else {
		target = "Real";
	}

 	coin_center = cen_ana.center_mean[pos];

	

	cout << " Before "<< acccc_temp << endl; 

	acccc_temp = acccc_temp + (charge * tot_eff);


	charge_tot = charge_tot + charge; 


	cout << "After " << acccc_temp << endl;

	errrr_temp = errrr_temp + pow(charge*tot_eff,2) * (pow(charge_err, 2) + pow(tot_eff_err,2)); 

//	errrr_temp = errrr_temp + pow(charge,2) * (pow(charge_err, 2) + pow(tot_eff_err,2)); 

	cout << "Now analyzing " << target << " Target Run #" << run_num << ": coin center=" 
		 << coin_center << endl;


	cout << "acccc: " << acccc_temp << "   charge: " << charge << "   tot_eff: " 
		 << tot_eff << endl; 

//	cout << "file file file     "<< file_out_ana->GetName() << endl;

	root_out_dir_name.Form("%i", run_num);
	
	out_dir = file_out_ana->mkdir(root_out_dir_name);




	/*--------------------------------------------------*/
	/// Calculate phi and t bins

	yield_err_sum.resize(t_bin_set);

	for (unsigned int i = 0; i < t_bin_set; ++i) {
		yield_err_sum[i].resize(phi_bin_num);
	}



}


/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Analysis::Para_Run_Clean() {

	real_event[0]->Reset();
	real_event[1]->Reset();
	real_event[2]->Reset();
	real_event[3]->Reset();
	real_event[4]->Reset();
	real_event[5]->Reset();

	delete coin_all;
	delete coin_real;
	delete coin_rand;

}



/*--------------------------------------------------*/
/*--------------------------------------------------*/


void Analysis::Load_Data_Tree() {

	///*--------------------------------------------------*/
	// Read in Nturple files

	// cout << " ::::  "  << run_num << endl;

	data_file.Form( data_file_dir + "coin%i.root", run_num);


	// cout << data_file_dir << endl;

	file_in = TFile::Open(data_file);

    data_tree_in = (TTree*) file_in->Get("h9500");
// 	data_tree_in->SetBranchAddress("missmass", &t_missmass);
// 	data_tree_in->SetBranchAddress("W", 	   &t_W       );
// 	data_tree_in->SetBranchAddress("t",        &t_t       );
// 	data_tree_in->SetBranchAddress("Q2",       &t_Q2      );
// 	data_tree_in->SetBranchAddress("th_pq",    &t_th_pq   );
// 	data_tree_in->SetBranchAddress("phi_pq",   &t_phi_pq  );
// 	data_tree_in->SetBranchAddress("Pmpar",    &t_PmPar   );
// 	data_tree_in->SetBranchAddress("Pmper",    &t_PmPer   );
// 	data_tree_in->SetBranchAddress("Pmoop",    &t_PmOop   );
// 
// 	data_tree_in->SetBranchAddress("hsdelta",  &t_hsdelta );
// 	data_tree_in->SetBranchAddress("hsytar",   &t_hsytar  );
// 	data_tree_in->SetBranchAddress("hsyptar",  &t_hsyptar );
// 	data_tree_in->SetBranchAddress("hsxptar",  &t_hsxptar );
// 	data_tree_in->SetBranchAddress("ssdelta",  &t_ssdelta );
// 	data_tree_in->SetBranchAddress("ssytar",   &t_ssytar  );
// 	data_tree_in->SetBranchAddress("ssyptar",  &t_ssyptar );
// 	data_tree_in->SetBranchAddress("ssxptar",  &t_ssxptar ); 

}


/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Analysis::Define_Cuts() {

	all_coin_cut = HMS_cut() + " && "+ SOS_cut() + " && " + PID_cut() + " && " + Cointime_all(coin_center) + " && " + Diamond_cut() + " && " + Missingmass_cut(miss_mass_offset[kset]) + " && " + Set_t_limit(t_min, t_max);

	real_coin_cut = HMS_cut() + " && "+ SOS_cut() + " && " + PID_cut() + " && " + Cointime_primary_cut(coin_center) + " && " + Diamond_cut() + " && " + Missingmass_cut(miss_mass_offset[kset]) + " && " + Set_t_limit(t_min, t_max);

	rand_coin_cut = HMS_cut() + " && "+ SOS_cut() + " && " + PID_cut() + " && " + Cointime_random_cut(coin_center) + " && " + Diamond_cut() + " && " + Missingmass_cut(miss_mass_offset[kset]) + " && " + Set_t_limit(t_min, t_max);



}




/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Analysis::Print_out_data() {


//	cout << endl << endl  << endl << endl << endl << endl;


 	for(unsigned int t_bin_set_tmp = 1; t_bin_set_tmp <= t_bin_set; t_bin_set_tmp++) {
 			
 		float t_c = t_min + (t_bin_set_tmp-0.5) * t_width;
 
 		for (float phi_step_tmp = pmn; phi_step_tmp <= pmx; phi_step_tmp = phi_step_tmp + pstp) {
 

			float yield_run = 0;
			float yield_run_err = 0;


 			int phi_itt = (phi_step_tmp + pstp/2) / 22.5;
 
 			phi = phi_step_tmp;			
  			tb = t_c;
 
 			yield_run = real_event[t_bin_set_tmp-1]->GetBinContent(phi_itt);

// 			if (yield == 0){ 
// 				yield_err = 0;
// 			} else{
// 				yield_err = yield_err_sum[t_bin_set_tmp-1][phi_itt-1]/pow(yield, 2)  + errrr_temp/pow(acccc_temp, 2);
// 			}
// 
//  			yield = yield/acccc_temp;
//  
//  			yield_err =  abs(yield) * sqrt(yield_err);
// 
//  			yield = yield/1000.;
//  			yield_err = yield_err/1000.;



			if (yield_run == 0){ 
				yield_run_err = 0;
			} else{
				yield_run_err = yield_err_sum[t_bin_set_tmp-1][phi_itt-1];
			}

			yield = yield_run;
 			yield_err = yield_run_err;
			
			run_tree->Fill();
			 
			int itt_tmp = (t_bin_set_tmp-1)*phi_bin_num + phi_itt - 1;
			
//			yield_vec[itt_tmp]     = yield;
//			yield_err_vec[itt_tmp] = yield_err;

			yield_vec[itt_tmp]     = yield_vec[itt_tmp] + yield;
			yield_err_vec[itt_tmp] = yield_err_vec[itt_tmp] + yield_err;

			phi_vec[itt_tmp]       = phi; 
			tb_vec[itt_tmp]        = tb;

			cout << "*  "; 
	 		cout << yield << "     " << yield_err << "     " << phi << "     " << tb << endl;

		}
 	}

	out_dir->cd();

	run_tree->Write("yield");
	run_tree->SetName("tree_out");

	chain->Add("tree_out");



	//  run_out->Reset();

	list->Add(run_tree);

	

}



/*--------------------------------------------------*/
/*--------------------------------------------------*/


void Analysis::Missing_Mass_Plot() {

	cointime_setting  = new TH1F( "", "", coin_bin, coin_min, coin_max);


	cout << "bin: " << coin_bin << "    Min: "<<  coin_min << "     Max: " <<  coin_max << endl;

	coin_all  = new TH1F("all", "all",   coin_bin, coin_min, coin_max);
	coin_real = new TH1F("real", "real", coin_bin, coin_min, coin_max);
	coin_rand = new TH1F("rand", "rand", coin_bin, coin_min, coin_max);



	data_tree_in->Draw("cointime>>all",  all_coin_cut, "goff");
	data_tree_in->Draw("cointime>>real", real_coin_cut, "goff");
	data_tree_in->Draw("cointime>>rand", rand_coin_cut, "goff");



// 	cout << kset << "   " << miss_mass_offset[kset] << "   " << miss_mass_offset[1] << endl;
// 	cout << all_coin_cut << endl;
 

	c1->cd();
	c1->Update();
	
	coin_all->Draw("hist");
	coin_all->SetLineColor(1);
	coin_real->Draw("same");
	coin_real->SetLineColor(4);
	coin_rand->Draw("same");
	coin_rand->SetLineColor(2);

	
	TLine *line = new TLine();

  	line->SetLineColor(kBlue);
  	line->SetLineWidth(3);

	line->DrawLine(coin_center+1,0,coin_center+1,25);
	line->DrawLine(coin_center-1,0,coin_center-1,25);
	line->DrawLine(coin_center+7,0,coin_center+7,25);
//	line->DrawLine(1,0,1,25);

// 	test->Update();
// 	test->cd();
// 
// 
// 	TH1D* test_h = new TH1D("test_h", "test_h", 100, -5, 10);
// 
// 	// data_tree_in->Draw("cointime>>test_h", test_cut);
// 	// data_tree_in->Draw("cointime");
// 
// 	data_tree_in->Draw("cointime>>test_h");
// 
// 	test->Update();
// 

	c2->Update();	


	TH1D* mm   = new TH1D("mm", "mm", 50, -0.1, 0.12);
	TH1D* mm_1 = new TH1D("mm_1", "mm1", 50, -0.1, 0.12);
	

	TString expected_mm_str;
	expected_mm_str.Form("%f", expected_mm);
	

	//data_tree_in->Draw("missmass-0.939565 >> mm_1", rand_coin_cut, "goff");

	data_tree_in->Draw("missmass- " + expected_mm_str + " >> mm", real_coin_cut, "goff");
	data_tree_in->Draw("missmass- " + expected_mm_str + " >> mm_1", rand_coin_cut, "goff");

	mm->Add(mm_1, -0.3333333);
 	mm->SetMarkerStyle(3);

	TH1F* mm_fit = (TH1F*) mm->Clone();



	c2->cd(1);

	mm_fit->Draw("PE");


	TGraphErrors* mm_diff_err = new TGraphErrors(mm_fit);

	TF1 *f1 = new TF1("f1", "gaus", -0.025, 0.025);

	f1->SetParameter(0, 0);
	f1->SetParameter(1, 0);
	f1->SetParameter(2, 0);


	mm_fit->Draw("EP");
	mm_fit->Fit("f1", "R");


	c2->Update();	



	/// Print missing mass text on the plot
	TString missing_mass_offset;
	missing_mass_offset.Form("Missing Mass Offset: %f", miss_mass_offset[kset]);

	TText* offset_txt = new TText();

	miss_mass_offset_str.Form("missmass-" + expected_mm_str + "-%f", miss_mass_offset[kset]);

	TH1F* mm_off = new TH1F("mm_off", "mm_off", 50, -0.1, 0.12);
	TH1F* mm_1_off = new TH1F("mm_1_off", "mm1_off", 50, -0.1, 0.12);

	data_tree_in->Draw(miss_mass_offset_str + " >> mm_off", real_coin_cut, "goff");
	data_tree_in->Draw(miss_mass_offset_str + " >> mm_1_off", rand_coin_cut, "goff");

	mm_off->Add(mm_1_off, -0.3333333);
 	mm_off->SetMarkerStyle(3);

	TH1F* mm_off_fit = (TH1F*) mm_off->Clone();

	c2->cd(2);


	TF1 *f2 = new TF1("f2", "gaus", -0.025, 0.025);


	mm_off_fit->SetMarkerColor(2);
	mm_off_fit->SetMarkerStyle(27);
	mm_off_fit->SetLineColor(2);
	mm_off_fit->Draw("PE");
	mm_off_fit->Fit("f2", "R");

	offset_txt->SetTextSize(0.05);
	offset_txt->DrawTextNDC(0.65, 0.7, missing_mass_offset);

	c2->Update();

	

	c4->cd();
	c4->Update();

	mm_off->SetMarkerColor(2);
	mm_off->SetMarkerStyle(27);
	mm_off->SetLineColor(2);

	mm->Draw("ep");
	mm_off->Draw("esame");




	offset_txt->SetTextSize(0.035);
	offset_txt->DrawTextNDC(0.65, 0.7, missing_mass_offset);

	
	/// Just to rescale the mm plot if mm_off points are out of the Y range
	/// Interestingly mm->GetYaxis()->GetXmax() acutually works !
	if ( mm_off->GetMaximum() > mm->GetYaxis()->GetXmax() ) {
		mm -> SetMaximum(mm_off->GetMaximum() + mm_off->GetMaximum() * 0.15);
	}	



	c4->Update();


 	c3->cd();
	c3->Update();	




	///*--------------------------------------------------*/
	//// Methods to access the TGraph from tree
	//
	// TGraph *gr = (TGraph *)gPad->GetPrimitive("Graph")->Clone(); // 2D
	// TGraph *gr = new TGraph(tree->GetSelectedRows(), tree->GetV2(), tree->GetV1());


	without_diamond_cut = HMS_cut() + " && "+ SOS_cut() + " && " + PID_cut() + " && " + Missingmass_cut(miss_mass_offset[kset]) + " && " + Set_t_limit(t_min, t_max);

	data_tree_in->Draw("W:Q2", without_diamond_cut);

	TGraph *gr1 = (TGraph*)gPad->GetPrimitive("Graph")->Clone();
	gr1->SetMarkerStyle(7);	
	gr1->SetMarkerSize(3);	
	gr1->SetMarkerColor(1);	

	data_tree_in->Draw("W:Q2", real_coin_cut);
	TGraph *gr2 = (TGraph* )gPad->GetPrimitive("Graph")->Clone();

	gr2->SetMarkerStyle(7);	
	gr2->SetMarkerSize(1);	
	gr2->SetMarkerColor(2);	

	gr1->Draw("AP");
	gr2->Draw("P");

	
 	TText* number_points = new TText();
 	number_points->SetTextSize(0.035);

 	TString points_str; 
 	points_str.Form("Black Events: %i", gr1->GetN());
 	number_points->DrawTextNDC(0.7, 0.65, points_str);

 	points_str.Form("Red Events: %i", gr2->GetN());
 	number_points->DrawTextNDC(0.7, 0.75, points_str);

 	points_str.Form("Total Events: %i", gr1->GetN()+gr2->GetN());
 	number_points->DrawTextNDC(0.7, 0.85, points_str);

	diamond_setting->Add(gr1);
	diamond_setting_cut->Add(gr2);
 
	file_out_ana->cd();
	out_dir->cd();


	mm_fit->SetOption("ep");
	mm_off_fit->SetOption("ep");

	mm_fit->Write();
	mm_off_fit->Write();




	c1->Write("Cointime");
	c2->Write("mm");
	c3->Write("diamond");
	c4->Write("mm_same");

	
/////	diamond_setting->Add(, 1);   

// 	cointime_setting->Add(coin_all, 1);
 	mm_setting->Add(mm, 1);        
 	mm_offset_setting->Add(mm_off, 1); 
// 

	file_out_ana->cd();


}

		
/*--------------------------------------------------*/
/*--------------------------------------------------*/

TTree* Analysis::Ceate_File() {

 	TTree* exam = new TTree();

	exam->Branch("yield",     &yield, "yield/F");
	exam->Branch("yield_err", &yield_err, "yield_err/F");
	exam->Branch("phi",       &phi, "phi/F");
	exam->Branch("tb",        &tb, "tb/F");

	return exam;

}



/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Analysis::Yield_Out() {



	file_out_ana->cd();

	for (unsigned int itt=0; itt < t_bin_set*phi_bin_num; itt++) {

		// cout << itt << endl;
		yield     = yield_vec[itt];    
		yield_err = yield_err_vec[itt];
		phi       = phi_vec[itt];      
		tb        = tb_vec[itt];       

		if (yield == 0){ 
			yield_err = 0;
 		} else{
			yield_err = yield_err/pow(yield, 2)  + errrr_temp/pow(acccc_temp, 2);
		}
 
		yield  = yield/acccc_temp;
  		yield_err =  abs(yield) * sqrt(yield_err);

		yield = yield/1000;
 		yield_err = yield_err/1000.;

	 	cout << "********  ";
	 	cout << yield << "     " << yield_err << "     " << phi << "     " << tb << endl;

		tree_out->Fill();		

	}


	// TTree* tree_out = new TTree("out_dir", "out_dir");
	// tree_out->Write("yield");
	// delete tree_out;	

	// tree_out->Merge(list);

	//chain->Merge("out_dir");

	tree_out->Write("yield");

//	Overall_Dia_Plot();
		
	cointime_setting->Write("cointime_setting");
	mm_setting->Write("mm_setting"); 
	mm_offset_setting->Write("mm_offset_setting"); 

}



/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Output the overall diamond plot for the setting

void Analysis::Overall_Dia_Plot() {

	TCanvas* dia_canva = new TCanvas();

	dia_canva->Update();
		
	diamond_setting->Draw("AP");
	diamond_setting_cut->Draw("P");

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


	TList * graph_list_cut = diamond_setting_cut->GetListOfGraphs();
 	TIter next_cut(graph_list_cut);
 	TObject* object_cut = 0;
	Int_t event_num_cut = 0;

 	while ((object_cut = next_cut()))
 	{
		event_num_cut = event_num_cut + ((TGraph*)object_cut)->GetN();
 	}






//	cout << "  fsafdsf  " << event_num << endl;

 	TText* number_points = new TText();
 	number_points->SetTextSize(0.035);

	TString points_str; 
 	points_str.Form("Total Events: %i", event_num + event_num_cut);
 	number_points->DrawTextNDC(0.7, 0.85, points_str);

 	points_str.Form("Black Events: %i", event_num);
 	number_points->DrawTextNDC(0.7, 0.75, points_str);

 	points_str.Form("Red Events: %i", event_num_cut);
 	number_points->DrawTextNDC(0.7, 0.65, points_str);




 	dia_canva->Update();
 	dia_canva->Write("diamand");

	delete number_points;
	delete dia_canva;


}



// void Analysis::SetDir(TFile* set_file) {
// 
// 
// 
// 
// }


/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Analysis::Set_MM_Offset(Double_t* array) {
	
	miss_mass_offset = new Double_t[6] ;

	for(int i=0; i<6; ++i) {
		miss_mass_offset[i] = array[i];
	}

}








/*--------------------------------------------------*/
/*--------------------------------------------------*/

Analysis::~Analysis() {



}


