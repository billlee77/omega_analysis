#include <iostream>

#include "analysis_heep_single.h"
//#include "cut.h"

#define PI 3.14159265

using namespace std;


/*--------------------------------------------------*/

Analysis_heep_single::Analysis_heep_single() {


	chain = new TChain();

}

/*--------------------------------------------------*/

Analysis_heep_single::Analysis_heep_single(ReadFile::efficiency_profile eff_struc) {
	
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

	Define_Cuts();

}


/*--------------------------------------------------*/
/// Analysis for Heep 

void Analysis_heep_single::Heep_anna_single(Int_t run_itt) {
	
	Define_Cuts();

	cout << "Analyzing Heep data !" << endl;

	dummy_correction  = 0.142;

 	run_tree = Ceate_File();
 
 	Para_Run_Def(run_itt);

	run_num = Get_Run_Num();

	ssxptar_offset = 3.2;
	ssyptar_offset = 0.0;
	

	Load_Data_Tree();


 	chain->Add(Get_Current_Run_File()+"/h9020");

	Check_Out(data_tree_in);

	Calculate_Yield(w_h);
		
 	out_dir->cd();

	Write_his();

	//Correct_Offset_Tree();

	list->Add(data_tree_in);

}

 
/*--------------------------------------------------*/
/// Dufine cut

void Analysis_heep_single::Define_Cuts() {

	single_cut = SOS_cut();


//	w_cut_up   = 1.05;

	w_cut_up   = 1.05;
	w_cut_down = 0.85;



}


/*--------------------------------------------------*/
 
 Analysis_heep_single::~Analysis_heep_single() {
 
}


/*--------------------------------------------------*/

void Analysis_heep_single::Load_Data_Tree() {

	///*--------------------------------------------------*/
	// Read in Nturple files

	data_file.Form( data_file_dir + "sos%i.root", run_num);

	file_in = TFile::Open(data_file);

    data_tree_in = (TTree*) file_in->Get("h9020");

	cout << data_file << endl;


//	cout << SOS_cut() << endl;
	




// 	data_tree_in->SetBranchAddress("ssdelta",  &ssdelta  );
// 	data_tree_in->SetBranchAddress("ssyptar",  &ssyptar  );
// 	data_tree_in->SetBranchAddress("ssxptar",  &ssxptar  ); 
// 	data_tree_in->SetBranchAddress("w",  	   &w_data   ); 
// 	data_tree_in->SetBranchAddress("Q2",       &q2_data  ); 
// 
// 	data_tree_in->SetBranchAddress("sszbeam",  &sszbeam  ); 
// 	data_tree_in->SetBranchAddress("ssp",      &ssp      ); 
// 	data_tree_in->SetBranchAddress("ssenergy", &ssenergy ); 
// 
// 



}



/*--------------------------------------------------*/

void Analysis_heep_single::Check_Plot_Out(TString plot_name) {
	
 	out_dir->cd();
	
	data_tree_in->Draw(plot_name, single_cut);

 	TH1F *htemp_real = (TH1F*)gPad->GetPrimitive("htemp");
 	
 	htemp_real->Write(plot_name);

}



/*--------------------------------------------------*/

void Analysis_heep_single::Heep_sim_anna_single(Int_t run_itt) {


	Define_Cuts();


 	TString sim_file;
 	sim_file = "list_sim_normfac.dat"; 

//	cout << setting_str << endl;
 
//	ReadFile* rff_sim = new ReadFile("list_sim_normfac.dat");
	ReadFile* rff_sim = new ReadFile(sim_file);
	sim_ana = rff_sim->sim_pro;


	Int_t q2_f, q2_r;
	Int_t ebeam_f, ebeam_r;

	q2_f    = floor(sim_ana.q2_setting[run_itt]);

	q2_r    = round(sim_ana.q2_setting[run_itt]*1000.) - q2_f *1000.;

	ebeam_f = floor(sim_ana.ebeam[run_itt]);
	ebeam_r = round(sim_ana.ebeam[run_itt]*1000.) - ebeam_f *1000.;

	cout <<"MMMMM " << ebeam_r%10 << endl;

	if (ebeam_r%10==0) {

		ebeam_r = ebeam_r/10;
	}

	TString sim_in_file_name;
	sim_in_file_name.Form("heep_ebeam_%ip%i_q2_%ip%i_single", ebeam_f, ebeam_r, q2_f, q2_r);


	

	

	cout << sim_in_file_name << endl;
	
 	Int_t event = sim_ana.event_num[run_itt];
 	Double_t normfac = sim_ana.normfac[run_itt];


	TString run_number;

	run_number = sim_in_file_name;

	data_file_dir = "sim_data/";
	data_file = data_file_dir + run_number + ".root";

	if (!rff_sim->File_exist(data_file)) return;

	file_in = TFile::Open(data_file);

//	acceptence_cut = SOS_cut() + " && W > 0.80 && W < 1.05";

	TString w_cut_tmp; 

	w_cut_tmp.Form("W > %f && W < %f", w_cut_down, w_cut_up);

	cout << w_cut_tmp << endl;

//	exit(0);

	acceptence_cut = Heep_Single_SOS_cut() + " && " + w_cut_tmp ;

	data_tree_in = (TTree*)file_in->Get("h666");

	sim_selection.Form("%f*Weight*", normfac/event);
	
	sim_selection = sim_selection + "(" + acceptence_cut + ")";

	w_h_org       = new TH1F("w_h_org",  "w_h_org",  200, 0.5, 1.5);

	w_h           = new TH1F("w_h",  "w_h",  200, 0.5, 1.5);

	q2_h          = new TH1F("q2_h", "q2_h", 200, 4, 6);

	ssdelta_h     = new TH1F("ssdelta_h", "ssdelta_H", 200, -20, 20);
	ssxptar_h     = new TH1F("ssxptar_h", "ssxptar_h", 200, -0.05, 0.05);
	ssyptar_h     = new TH1F("ssyptar_h", "ssyptar_h", 200, -0.1, 0.1);


	data_tree_in->Draw("W>>w_h",    sim_selection, "goff");
	data_tree_in->Draw("W>>w_h_org", acceptence_cut, "goff");

	data_tree_in->Draw("ssdelta>>ssdelta_h", sim_selection, "goff");
	data_tree_in->Draw("ssxptar>>ssxptar_h", sim_selection, "goff");
	data_tree_in->Draw("ssyptar>>ssyptar_h", sim_selection, "goff");

    
 	file_out_ana->cd();

	w_h->Write("W");

	ssdelta_h->Write("ssdelta");
	ssxptar_h->Write("ssxptar");
	ssyptar_h->Write("ssyptar");





//	TString sim_data;
	
//	sim_data.Form(data_file_dir + "coin%i.root", run_number);

//	cout << sim_data << endl;



	cout << data_file << endl;

 	chain->Add(data_file+"/h666");

	Acceptance_check(Get_Acceptence_Vector(), sim_selection, "sime_data");


//	Calculate_Yield(w_h);




/*--------------------------------------------------*/
/*--------------------------------------------------*/



	cout << "Now Calculating the yield: " << endl;
	
	// cout << "Get Bins: " << mm->GetNbinsX() << endl;

	float yield_raw       = 0.0;		
	float yield_weighted  = 0.0;		

	float yield_err   = 0.0;		
	float yield_err_2 = 0.0;		

	for(int i=0; i <= w_h_org->GetNbinsX();i++) {

//		cout << i << "   " << mm->GetBinContent(i) << endl;
		
		yield_raw      = yield_raw + w_h_org->GetBinContent(i) ;
		yield_weighted = yield_weighted + w_h->GetBinContent(i) ;
		
		
	}

	yield_err = sqrt(yield_err_2);

	cout << "Final raw yield : "   << yield_raw << endl;
	cout << "Weighed raw yield : " << yield_weighted << endl;



	float average_weight = yield_weighted / yield_raw;
	float yield_weighted_error = sqrt(yield_raw) * average_weight;

	cout << "Average Weight: " << average_weight << endl;
	cout << "Garth's Error:  " << yield_weighted_error << endl;

//	cout << "Sim_check : "<< mm_try->Integral(0, -1) << "    " <<  sqrt(mm_try->Integral(0, -1)) * normfac/event << endl;
	cout << "Sim_check : "<< w_h->Integral(0, -1) << "    " << endl;

//	float pp = upppp / event;

//	cout << "New yield error : " <<   sqrt(upppp)/event << endl;

	cout << endl;
	cout << endl;

	yield_setting     = yield_weighted;
	yield_setting_err = yield_weighted_error;



	delete w_h;
	delete ssdelta_h;
	delete ssxptar_h;
	delete ssyptar_h;


	Yield_output();

}

 

/*--------------------------------------------------*/

float Analysis_heep_single::Get_SOS_Momentum(float q2_tmp){

	float sos_p;

	if (q2_tmp == float(1.91)) {
		sos_p = 1.5808;
	} else if (q2_tmp == float(6.53)) {
		sos_p = 1.74;
	} else if (q2_tmp == float(4.42)) {
		sos_p = 1.42;
	} else if (q2_tmp == float(5.42)) {
		sos_p = 1.74;
	} else {
		cout << "asadaedasdas "<< endl;
	}

//	cout << "{{{{{{{{{{{{{{{" << sos_p << endl;
	
	return sos_p;

}


float Analysis_heep_single::Get_SOS_Angle(float q2_tmp){

	float sos_angle;

	if (q2_tmp == float(1.91)) {
		sos_angle = 51.03;
	} else if (q2_tmp == float(6.53)) {
		sos_angle = 50;
	} else if (q2_tmp == float(4.42)) {
		sos_angle = 54;
	} else if (q2_tmp == float(5.42)) {
		sos_angle = 48;
	} else {
		cout << "asadaedasdas "<< endl;
	}

//	cout << "{{{{{{{{{{{{{{{" << sos_p << endl;
	
	return sos_angle;

}



/*--------------------------------------------------*/
/*--------------------------------------------------*/

void Analysis_heep_single::Yield_Out() {



	/*--------------------------------------------------*/
	/*--------------------------------------------------*/

	cout << "Check Error: " << yield_setting_err/pow(yield_setting, 2) << "      " 
		 <<  errrr_temp/pow(acccc_temp, 2) << endl;

	yield_setting_err = yield_setting_err/pow(yield_setting, 2)  + errrr_temp/pow(acccc_temp, 2);

	yield_setting     = yield_setting / acccc_temp ;
  	yield_setting_err =  yield_setting * sqrt(yield_setting_err);

	yield_setting     = yield_setting / 1000.;
	yield_setting_err = yield_setting_err / 1000.;

	cout << "Total yield:       " << yield_setting     << endl;
	cout << "Total yield Error: " << yield_setting_err << endl;
//	cout << "Total charge:      " << charge_setting    << endl;
	cout << "Effective charge:  " << acccc_temp        << endl;

//	cout << "Real yield Unit :  " << yield_setting / charge_setting /1000 << endl;
//	cout << "Real yield Unit :  " << yield_setting / charge_setting /1000 * 2.4 << endl;
//	cout << "Real yield Unit :  " << yield_setting << endl;

	cout << endl;
	cout << endl;


	/*--------------------------------------------------*/
	/*--------------------------------------------------*/


	Check_Out(chain);

	Float_t scale_factor;

	scale_factor = 0.0;

	scale_factor = 1.0/(acccc_temp*1000);

	cout  << " WWWWWWWWWWWW !!!! : " << dummy_correction << endl;

 	if (is_run_dummy)
 		scale_factor =  dummy_correction * scale_factor;
// 		scale_factor =  0.142 * scale_factor;
// //		scale_factor = 0.385 * 0.463 * scale_factor;
// //		scale_factor = 0.15 * scale_factor;


//	cout  << " WWWWWWWWWWWW !!!! : " << scale_factor << endl;
//	exit(0);


	cout << "*** " << is_run_dummy << " :::  " << scale_factor << "  accc " << acccc_temp << endl;

// 	missingmass_h->Scale(scale_factor);
// 
// 	ssdelta_h 	 ->Scale(scale_factor);
// 	ssxptar_h 	 ->Scale(scale_factor);
// 	ssyptar_h 	 ->Scale(scale_factor);
// 
// 	ssenergy_h	 ->Scale(scale_factor);
// 	sszbeam_h 	 ->Scale(scale_factor);
// 
// 	ssp_h     	 ->Scale(scale_factor);
// 
// 	w_h       	 ->Scale(scale_factor);
// 	q2_h      	 ->Scale(scale_factor);
// 
// 	E_m_h     	 ->Scale(scale_factor);
// 
// 	P_mx_h    	 ->Scale(scale_factor);
// 	P_my_h    	 ->Scale(scale_factor);
// 	P_mz_h    	 ->Scale(scale_factor);

	Set_Hist_Err(missingmass_h);
	             
	Set_Hist_Err(ssdelta_h);
	Set_Hist_Err(ssxptar_h); 	 
	Set_Hist_Err(ssyptar_h); 	 
	             
	Set_Hist_Err(ssenergy_h);	 
	Set_Hist_Err(sszbeam_h); 	 
	             
	Set_Hist_Err(ssp_h);     	 
	             
	Set_Hist_Err(w_h);       	 
	Set_Hist_Err(q2_h);      	 
	             	
	Set_Hist_Err(E_m_h);     	 
	             
	Set_Hist_Err(P_mx_h);    	 
	Set_Hist_Err(P_my_h);   	 
	Set_Hist_Err(P_mz_h);    	 


	P_m_h     	 ->Scale(scale_factor);

	w_p       	 ->Scale(scale_factor);

	delta_m   	 ->Scale(scale_factor);



	/*--------------------------------------------------*/
	/*--------------------------------------------------*/


//	data_scale_factor = scale_factor;
//	data_scale_factor = scale_factor;


	Set_Scale_Factor(scale_factor);

	file_out_ana->cd();

	Write_his();

	Acceptance_check(Get_Acceptence_Vector(), Heep_Single_PID_cut() + " && w > 0.8 && w < 1.05", "data");




	cout << "{{{{{{{{{{{{{{{   " << yield_setting << "    " << yield_setting_err << endl;


	Yield_output();


}


/*--------------------------------------------------*/
/*--------------------------------------------------*/



void Analysis_heep_single::Check_Out(TTree* tree_tgt_org) {


	file_out_ana->cd();

//	TTree * tree_tgt = tree_tgt_org->CopyTree(SOS_cut());
//	TTree * tree_tgt = tree_tgt_org->CopyTree(SOS_cut() + " && " + Heep_Single_PID_cut() + " && evtype < 2.5" );

//	TTree * tree_tgt = tree_tgt_org->CopyTree(SOS_cut() + " && " + Heep_Single_PID_cut());
//
//
//
	TString w_cut_tmp; 

	w_cut_tmp.Form("w > %f && w < %f", w_cut_down, w_cut_up);

	TTree * tree_tgt = tree_tgt_org->CopyTree(Heep_Single_SOS_cut() + " && " + Heep_Single_PID_cut() + " && " + w_cut_tmp);

//	TTree * tree_tgt = tree_tgt_org->CopyTree("");

	tree_tgt->SetBranchAddress("ssdelta",  &ssdelta  );
	tree_tgt->SetBranchAddress("ssyptar",  &ssyptar  );
	tree_tgt->SetBranchAddress("ssxptar",  &ssxptar  ); 
	tree_tgt->SetBranchAddress("w",  	   &w_data   ); 
	tree_tgt->SetBranchAddress("Q2",       &q2_data  ); 

	tree_tgt->SetBranchAddress("sszbeam",  &sszbeam  ); 
	tree_tgt->SetBranchAddress("ssp",      &ssp      ); 
	tree_tgt->SetBranchAddress("ssenergy", &ssenergy ); 

	tree_tgt->SetBranchAddress("evtype",   &evtype ); 

// 	tree_tgt->SetBranchAddress("ssxtar",  &ssxtar_data );
// 	tree_tgt->SetBranchAddress("ssytar",  &ssytar_data );
// 

	ssdelta_h     = new TH1F("ssdelta", "ssdelta", 200, -20, 20);
	ssxptar_h     = new TH1F("ssxptar", "ssxptar", 200, -0.05, 0.05);
	ssyptar_h     = new TH1F("ssyptar", "ssyptar", 200, -0.1, 0.1);

//	ssxtar_h     = new TH1F("ssxtar", "ssxtar", 200, -3, 3);
//	ssytar_h     = new TH1F("ssytar", "ssytar", 200, -3, 3);


	w_h           = new TH1F("w_h",  "w_h",  200, 0.5, 1.5);
	q2_h          = new TH1F("q2_h", "q2_h", 200, 4, 6);

	ssenergy_h    = new TH1F("ssenergy", "ssenergy", 200, 0, 5);
	sszbeam_h     = new TH1F("sszbeam",  "sszbeam" , 200, -5, 10);
	
	ssp_h         = new TH1F("ssp", "ssp", 200, 0, 5);

	missingmass_h = new TH1F("missmass", "missmass", 200, 0, 3);

	E_m_h         = new TH1F("E_m_h",  "E_m_h",  200, 3, 5); 
	P_m_h         = new TH1F("P_m_h",  "P_m_h",  200, 0, 5); 

	P_mx_h        = new TH1F("P_mx_h", "P_mx_h", 50, -0.1, 0.1);
	P_my_h        = new TH1F("P_my_h", "P_my_h", 200, 0.5, 2);
	P_mz_h        = new TH1F("P_mz_h", "P_mz_h", 200, 2,   4);

	w_p           = new TH2F("w_p", "w_p", 100, 0, 3, 100, 0, 5);
	delta_m       = new TH2F("delta_m", "delta_m", 100, -40, 40, 100, 0, 5);

	test_h        = new TH1F("test_h", "test_h", 10, 0, 4);

	float E_missing;
	float P_missing;
	
	float E_beam;
	float P_beam;
	float P_x, P_y, P_z;

	float E_e = Get_E_Beam();
	float Q2 = Get_Q2();

	float sos_p_0 = Get_SOS_Momentum(Q2);
	float sos_angle = Get_SOS_Angle(Q2);

	float sos_p;



//	sos_p = sos_p * (1.0 - 0.00218);


	float P;

	float missingmass;

	for (unsigned int i; i < tree_tgt->GetEntries(); i++) {

		tree_tgt->GetEntry(i);

		if (w_data > w_cut_down && w_data < w_cut_up) {

//		if (ssdelta < 40 && ssdelta > -40 && ssxptar > -0.2 && ssxptar < 0.3 && ssyptar > -0.2 && ssyptar < 0.3) {
	
			ssdelta_h->Fill(ssdelta);
			ssxptar_h->Fill(ssxptar);
			ssyptar_h->Fill(ssyptar);

			sszbeam_h ->Fill(sszbeam);
			ssenergy_h->Fill(ssenergy);
			ssp_h->Fill(ssp);

			w_h ->Fill(w_data);
			q2_h->Fill(q2_data);


//			cout << sos_p << endl;


// 			cout  << " 1.0   " <<  Get_Real_SOS_P(1.0)<< endl;
// 			cout  << " 1.7   " <<  Get_Real_SOS_P(1.7)<< endl;
// 

			
			sos_p = Get_Real_SOS_P(sos_p_0);


			P = sos_p * (1. + ssdelta/100.);



//			cout << "1.5808  :   " << Get_Real_SOS_P(1.5808) << endl;
//			cout << "1.4196  :   " << Get_Real_SOS_P(1.4196) << endl;
//			cout << "1.7400  :   " << Get_Real_SOS_P(1.74) << endl;


//			cout << P << endl;

//			exit(0);
			
//			P = ssp;

//			sos_p = ssp;
//			P = sos_p * (1 - ssdelta/100.);

//			E_e = sszbeam;



//			cout << ssxptar << "    " << ssxptar_offset/1000. <<  endl;

//			ssxptar = ssxptar + ssxptar_offset/1000.;
//			ssyptar = ssyptar + ssyptar_offset/1000.;



//			exit(0);


//			float angle_rad_angle_x = asin(ssxptar);
//			float angle_rad_angle_y = asin(ssyptar);

			float angle_rad_angle_x = ssxptar + ssxptar_offset/1000.;
//			float angle_rad_angle_x = ssxptar;
			float angle_rad_angle_y = ssyptar + ssyptar_offset/1000.;



			float p_x, p_y, p_z;

			p_x = P * sin(angle_rad_angle_x) * cos(angle_rad_angle_y);

			p_y = P * cos(angle_rad_angle_x) * sin(angle_rad_angle_y);

			p_z = P * cos(angle_rad_angle_x) * cos(angle_rad_angle_y);

			float p_x_cor = p_x;

			float p_y_cor = p_y * cos(sos_angle* PI / 180.0) + p_z * sin(sos_angle* PI / 180.0);

			float p_z_cor = E_e - (p_y * (-sin(sos_angle* PI / 180.0)) + p_z * cos(sos_angle* PI / 180.0));

//			cout << "  sadf  "  << p_y_cor << "     " << p_z_cor << endl;

			P_missing = sqrt( pow(p_x_cor, 2) + pow(p_y_cor, 2) + pow(p_z_cor, 2) );		

			E_missing = E_e + 0.938272 - P;

			missingmass = sqrt( pow(E_missing, 2) - pow(P_missing, 2) );

// 			cout << angle_rad_angle_x   << "    " << angle_rad_angle_y << endl;
// 
// 			cout << ssdelta   << "    " << ssxptar << "    " << ssyptar << endl;
// 			cout << sos_angle << "    " << E_e     << "    " << sos_p << "    " << P 
// 				 << "    " << Q2 << "    " << q2_data << endl;
// 
// 			cout << "Px:   "           << p_x       << "   " << p_x_cor    << endl;
// 			cout << "Py:   "           << p_y       << "   " << p_y_cor    << endl;
// 			cout << "Pz:   "           << p_z       << "   " << p_z_cor    << endl;
// 			cout << "Missing mass:   " << E_missing << "   " << P_missing  << "   " << missingmass << endl;
// 			cout << "E P Check:   "    << sszbeam   << "   " << ssenergy   << "   " << ssp         << endl;
// 
// 			cout << q2_data/(E_e -P) << endl;
// 			
// 			cout << endl;	
// 
// 
// 		    exit(0);
// // // 
//			missingmass_h->Fill(E_missing-2);

//			cout << w_data << "     " << missingmass << endl;


			missingmass_h->Fill(missingmass);

			E_m_h->Fill(E_missing);

			P_m_h->Fill(P_missing);

			P_mx_h->Fill(p_x_cor);
			P_my_h->Fill(p_y_cor);
			P_mz_h->Fill(p_z_cor);

			P_mz_h->Fill(P_missing);

			w_p->Fill(w_data, P_missing);
			delta_m->Fill(ssdelta, P_missing);
			

			test_h->Fill(evtype);

		}



	}




// 	delete ssdelta_h;
// 	delete ssxptar_h;
// 	delete ssyptar_h;
// 


//	exit(0);


}





void Analysis_heep_single::Write_his() {

	missingmass_h->SetOption("E");

	ssdelta_h ->SetOption("E");
	ssxptar_h ->SetOption("E");
	ssyptar_h ->SetOption("E");

	ssenergy_h->SetOption("E");
	sszbeam_h ->SetOption("E");

	ssp_h     ->SetOption("E");

	w_h       ->SetOption("E");
	q2_h      ->SetOption("E");

	E_m_h     ->SetOption("E");
	P_m_h     ->SetOption("E");
             
	P_mx_h    ->SetOption("E");
	P_my_h    ->SetOption("E");
	P_mz_h    ->SetOption("E");
            





	missingmass_h->Write("missmass");

	ssdelta_h ->Write("ssdelta");
	ssxptar_h ->Write("ssxptar");
	ssyptar_h ->Write("ssyptar");

	ssenergy_h->Write("ssenergy");
	sszbeam_h ->Write("sszbeam");

	ssp_h     ->Write("ssp");

	w_h       ->Write("W");
	q2_h      ->Write("Q2");

	E_m_h     ->Write("E_m_h");
	P_m_h     ->Write("P_m_h");
            
	P_mx_h    ->Write("P_mx_h");
	P_my_h    ->Write("P_my_h");
	P_mz_h    ->Write("P_mz_h");
            

	w_p       ->SetOption("colz");
	w_p       ->Write("w_p");

	delta_m   ->SetOption("colz");
	delta_m   ->Write("delta_m");

	test_h   ->Write("test_h");


	delete missingmass_h;

	delete ssdelta_h ;
	delete ssxptar_h ;
	delete ssyptar_h ;

	delete ssenergy_h;
	delete sszbeam_h ;

	delete ssp_h     ;

	delete w_h       ;
	delete q2_h      ;

	delete E_m_h     ;

	delete P_mx_h    ;
	delete P_my_h    ;
	delete P_mz_h    ;

	delete P_m_h     ;

	delete w_p       ;

	delete delta_m   ;

	delete test_h	 ;

}





float Analysis_heep_single::Get_Real_SOS_P(float p_0) {

	float sos_p_cor;

	float epsilon;

	float a = 0.002;
	float b = -0.17;

	if (p_0  >= 1.5 ){

//		cout << b * pow((p_0 - 1.5),2) << endl;

		epsilon = a + b * pow((p_0 - 1.5),2);

	} else {
		
		epsilon = 0.006 - ( 0.006 - 0.002) / 1.5 * p_0;

	}

//	cout << epsilon << endl;
	
	sos_p_cor = p_0 * (1.0 + epsilon);

	return sos_p_cor;

}


// 
// void Analysis_heep_single::Acceptance_check(vector<TString> list_vec, TString cut_str, TString data_type) {
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
// 
// //    	cout << ' ' << *it << endl;
// // 
// // 	exit(0);
// 
// }
// 

// 
// /*--------------------------------------------------*/
// 
// void Analysis_heep_single::Apt_Setting_Check(TString plot_name, TString cut_str, TString data_type) {
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
// 	}
// 
// 
// 
// 	// cout << " asdasda    " << cut_str << endl;
// 
// 
// 
// 
// 	chain->Draw(plot_op_str, cut_str);
// 
//   	TH1F* htempp = (TH1F*)gPad->GetPrimitive("htemp");
// 
// 
// 	if (data_type == "data") {
// 	
// 		htempp ->Scale(data_scale_factor);
// 
// 	} 
// 
// 
// 	htempp->Write(plot_name);
// 
// }
// 


vector<TString> Analysis_heep_single::Get_Acceptence_Vector() {

	static const TString list_arr[] = {"W", "Q2","ssdelta", "ssxfp", "ssxpfp", "ssytar", "ssxptar", "ssyptar"};

	vector<TString> vec (list_arr, list_arr + sizeof(list_arr) / sizeof(list_arr[0]) );

	return vec;

}







void Analysis_heep_single::Calculate_Yield(TH1F* hist_target) {

	yield = hist_target->Integral(0, -1); 

	yield_err = sqrt((float)yield);



 	if (is_run_dummy) {
		yield     = yield * dummy_correction;
		yield_err = yield_err * dummy_correction;
	} 


	yield_setting     = yield_setting + yield;
	yield_setting_err = yield_setting_err + pow(yield_err,2);


}



void Analysis_heep_single::Set_Hist_Err(TH1F* tgt_hist){


	Double_t real_bin;
	
	Double_t real_bin_err;

	Double_t sub_err_2;

	Double_t norm_yield;
	Double_t norm_yield_err;

	for (int i = 1;  i <= tgt_hist->GetNbinsX(); i++) {
		
		real_bin = tgt_hist->GetBinContent(i);
	
		Double_t real_bin_err = sqrt(real_bin);
	
		sub_err_2 = pow(real_bin_err,2);

		/*--------------------------------------------------*/
		/*--------------------------------------------------*/

		norm_yield = tgt_hist->GetBinContent(i);

		if (norm_yield != 0) { 
			norm_yield_err = sub_err_2/pow(norm_yield, 2)  + errrr_temp/pow(acccc_temp, 2);
		}

		norm_yield = norm_yield/acccc_temp;
			
		norm_yield_err =  norm_yield * sqrt(norm_yield_err);

		norm_yield = norm_yield / 1000;
		
		norm_yield_err = norm_yield_err / 1000;

		cout  << i << "   "   <<  tgt_hist->GetBinCenter(i) << "    "<< norm_yield << endl;
		
 		tgt_hist->SetBinContent(i, norm_yield);
 		tgt_hist->SetBinError(i, norm_yield_err);

 	}

}


