// Program for analysing the omega data 


void root_anna_pl(){


	Float_t x, y, z;

	const int phi_bin_num = 16;


	float pstp = 360./phi_bin_num;
	float pmn  = pstp/2.;
	float pmx  = 360 - pstp/2.;

	// cout << pmn << "   "  << pmx << endl;



	TCanvas *c1 = new TCanvas();
	TCanvas *c2 = new TCanvas();
	TCanvas *c3 = new TCanvas();



	TCanvas *c4 = new TCanvas();

	c2->Divide(1,2);

	Double_t miss_mass_offset[6] = {0.0, 0.00574, 0.00554, 0.00254, 0.00253, 0.00260};


	cout << " this is working" << endl;

//	TNtuple * eff_ntp = new TNtuple("eff","eff", "run_num:charge:charge_err:tot_eff:tot_eff_err:Q_2:epsilon:th0");
	TNtuple * setting_ntp = new TNtuple("setting","setting", "polset:Q2set:epsset:thpqset:tmnset:tmxset:NBtset:kset");

	setting_ntp->ReadFile("list.settings.fpi2");


	float yield;
	float yield_err;
	float phi;
	float tb;

	TTree *tree_out = new TTree("tree_out","Data output");
//	tree_out->Branch("run_num", &run_num_out, "run_num/I");
//	tree_out->Branch("run_random", &run_random_out, "run_random/I");

	tree_out->Branch("yield", &yield, "yield/F");
	tree_out->Branch("yield_err", &yield_err, "yield_err/F");
	tree_out->Branch("phi", &phi, "phi/F");
	tree_out->Branch("tb", &tb, "tb/F");



//	eff_ntp->Draw("run_num", "", "goff");

// 	/*--------------------------------------------------*/
// 	// Put everything into array
// 	// Method 1:
// 	
// 	eff_ntp->SetBranchAddress("run_num", &x);
// 	eff_ntp->SetBranchAddress("v1", &y);
// 	
// 	Int_t nentries = (Int_t)eff_ntp->GetEntries();
//  	eff_ntp->GetEntry(0);
// 
// 	TArrayD* xx = new TArrayD();
// 
// 	xx->Set(nentries);
// 
// 	for (Int_t i=0;i<nentries;i++) {
// 		eff_ntp->GetEntry(i);
// //		cout << x  << "   " << y << endl;
// 		xx->AddAt(x,i);
// 	}
// 
// 
// //	cout << xx->GetAt(0) << endl;
// //	cout << xx->GetAt(1) << endl;
// //	cout << xx->GetAt(2) << endl;
// 	
// 
// 
// 
// 	/*--------------------------------------------------*/
// 	// Put everything into array
// 	// Method 2:
// 
// 
// 	eff_ntp->Draw("run_num:v1", "", "goff");
// 
// 
// 	TH2F *htemp = gROOT->FindObject("htemp");
// 	or
//  TH1F *hnew = (TH1F*)gDirectory->Get("hnew");
//   or
//  TH1F *hnew = (TH1F*)gPad->GetPrimitive("hnew");
// 
// 	Double_t* xxx;
// 
// 	TGraph *gr = new TGraph(eff_ntp->GetSelectedRows(), eff_ntp->GetV1(), eff_ntp->GetV2());
// 	xxx = gr->GetX();
// 
// 	cout << xxx[0] << endl;
// 





	/*--------------------------------------------------*/
	// Put everything into array
	// Method 3:

	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	// Read Setting file
	//
	//
	


	setting_ntp->Draw("polset:Q2set:epsset", "", "goff");

	TArrayD* polset_tmp = new TArrayD(setting_ntp->GetSelectedRows(), setting_ntp->GetV1());
	TArrayD* Q2set_tmp  = new TArrayD(setting_ntp->GetSelectedRows(), setting_ntp->GetV2());
	TArrayD* epsset_tmp = new TArrayD(setting_ntp->GetSelectedRows(), setting_ntp->GetV3());
 
 	Double_t* polset = polset_tmp->GetArray();
 	Double_t* Q2set  = Q2set_tmp->GetArray();
 	Double_t* epsset = epsset_tmp->GetArray();


	setting_ntp->Draw("thpqset:tmnset:tmxset", "", "goff");

	TArrayD* thpqset_tmp  = new TArrayD(setting_ntp->GetSelectedRows(), setting_ntp->GetV1());
	TArrayD* tmnset_tmp   = new TArrayD(setting_ntp->GetSelectedRows(), setting_ntp->GetV2());
	TArrayD* tmxset_tmp   = new TArrayD(setting_ntp->GetSelectedRows(), setting_ntp->GetV3());

 	Double_t* thpqset = thpqset_tmp->GetArray();
 	Double_t* tmnset  = tmnset_tmp->GetArray();
 	Double_t* tmxset  = tmxset_tmp->GetArray();



	setting_ntp->Draw("NBtset:kset", "", "goff");

	TArrayD* NBtset_tmp  = new TArrayD(setting_ntp->GetSelectedRows(), setting_ntp->GetV1());
	TArrayD* kset_tmp    = new TArrayD(setting_ntp->GetSelectedRows(), setting_ntp->GetV2());

 	Double_t* NBtset = NBtset_tmp->GetArray();
 	Double_t* kset   = kset_tmp->GetArray();


	/*--------------------------------------------------*/
	/// Input file name
	TString dir_name ("lists/");  
	TString file_name ("list.dummy_245_27_plus");


	/*--------------------------------------------------*/
	/// Output file name

	TString out_dir_name ("file_out/");  
	TString out_file_name;
	out_file_name =  out_dir_name + file_name + "_out.root"; 


	TFile *file_out = new TFile(out_file_name, "RECREATE");
	


	///*--------------------------------------------------*/
	// Read Efficiency file, and put numbers nicely into array
	//
	/*--------------------------------------------------*/
	// Put everything into array
	// Method 3:


	// TNtuple * eff_ntp = new TNtuple("eff","eff", "run_num:charge:charge_err:tot_eff:tot_eff_err:Q_2:epsilon:th0");
	TNtuple * eff_ntp = new TNtuple("eff","eff", "run_num:charge:charge_err:tot_eff:tot_eff_err:Q_2:epsilon");
	eff_ntp->ReadFile(dir_name + file_name);

 	eff_ntp->Draw("run_num:charge:charge_err", "", "goff");
 //	eff_ntp->Draw("run_num");
 

 	TArrayD* run_num_tmp    = new TArrayD(eff_ntp->GetSelectedRows(), eff_ntp->GetV1());
 	TArrayD* charge_tmp     = new TArrayD(eff_ntp->GetSelectedRows(), eff_ntp->GetV2());
 	TArrayD* charge_err_tmp = new TArrayD(eff_ntp->GetSelectedRows(), eff_ntp->GetV3());
 



	Double_t* run_num_d = run_num_tmp->GetArray();

	



 	// Int_t array_size = sizeof(run_num_d)/sizeof(run_num_d[0]);

 	Int_t array_size = SizeOfArray(run_num_d);



 	Int_t* run_num = new Int_t[array_size];
 	
 	for(int i=0; i < array_size ;i++) {
  		run_num[i] = (int)run_num_d[i];

//		cout << " sadfasdfasdfdsa " << run_num_d[i] << endl;
 	}

	Double_t* charge     = charge_tmp->GetArray();
	Double_t* charge_err = charge_err_tmp->GetArray();


 	eff_ntp->Draw("tot_eff:tot_eff_err", "", "goff");
 
 	TArrayD* tot_eff_tmp = new TArrayD(eff_ntp->GetSelectedRows(), eff_ntp->GetV1());
 	TArrayD* tot_eff_err_tmp = new TArrayD(eff_ntp->GetSelectedRows(), eff_ntp->GetV2());
 
 	Double_t* tot_eff     = tot_eff_tmp->GetArray();
 	Double_t* tot_eff_err = tot_eff_err_tmp->GetArray();
 
 	// eff_ntp->Draw("Q_2:epsilon:th0", "", "goff");
 	eff_ntp->Draw("Q_2:epsilon", "", "goff");
 
 	TArrayD* Q2_tmp = new TArrayD(eff_ntp->GetSelectedRows(), eff_ntp->GetV1());
 	TArrayD* epsilon_tmp = new TArrayD(eff_ntp->GetSelectedRows(), eff_ntp->GetV2());
 	TArrayD* th0_tmp = new TArrayD(eff_ntp->GetSelectedRows(), eff_ntp->GetV3());
 
 
 	Double_t* Q_2 = Q2_tmp->GetArray();
 	Double_t* epsilon = epsilon_tmp->GetArray();
// 	Double_t* th0 = th0_tmp->GetArray();
 
 //	cout << yy->GetAt(0) << endl;
 //	cout << run_num[0] << "     " << tot_eff[2] << endl;
 


	///*--------------------------------------------------*/
	// Read center file
	/*--------------------------------------------------*/

	TNtuple* center_ntp = new TNtuple("center","center", "center_run_num:center_mean:center_sigma");
	center_ntp->ReadFile("fit_piplus/cointime_pl.dat");
 
	center_ntp->Draw("center_run_num:center_mean:center_sigma", "", "goff");

//	center_ntp->Draw("center_run_num:center_mean:center_sigma", "", "");
 
 	TArrayD* center_run_num_tmp = new TArrayD(center_ntp->GetSelectedRows(), center_ntp->GetV1());
 	TArrayD* center_mean_tmp    = new TArrayD(center_ntp->GetSelectedRows(), center_ntp->GetV2());
 	TArrayD* center_sigma_tmp   = new TArrayD(center_ntp->GetSelectedRows(), center_ntp->GetV3());
 


//	cout << center_run_num_tmp->At(5) << "  AAA   " << center_run_num_tmp->GetArray()[9] << endl;

 	Double_t* center_run_num_d = center_run_num_tmp->GetArray();

 	const int array_size_1 = SizeOfArray(center_run_num_d);


 	Int_t* center_run_num = new Int_t[array_size_1];
 	
 	for(int i=0; i < array_size_1 ;i++) {
  		center_run_num[i] = (int)center_run_num_d[i];
//		cout << i << "    " << center_run_num[i] << "  AAA   "  << array_size_1 << endl;
 	}

 	Double_t* center_mean = center_mean_tmp -> GetArray();
 	Double_t* center_sigma = center_sigma_tmp -> GetArray();
 

	Int_t run_itt;



	
	Int_t t_bin_set;
	float t_max, t_min, t_width;


	t_bin_set = NBtset[0];
	t_min = tmnset[0];
	t_max = tmxset[0];

	t_width = ( t_max - t_min ) / t_bin_set;

	cout << "OOOOOOOOOOOAsdassadfas   " << t_max << "   " << t_min << "   " << t_bin_set << "   "<<  t_width << endl;




	TH1F* real_event[6];


	real_event[0] = new TH1F("phi_real1", "", phi_bin_num, 0, 360);
	real_event[1] = new TH1F("phi_real2", "", phi_bin_num, 0, 360);
	real_event[2] = new TH1F("phi_real3", "", phi_bin_num, 0, 360);
	real_event[3] = new TH1F("phi_real4", "", phi_bin_num, 0, 360);
	real_event[4] = new TH1F("phi_real5", "", phi_bin_num, 0, 360);
	real_event[5] = new TH1F("phi_real6", "", phi_bin_num, 0, 360);



	float acccc_temp = 0.0;


	for(int i = 0; i < array_size; i++) {
		
		run_n = run_num[i];	
//		cout << run_n << endl;


	
		int pos = -100000000;
	
		for(int ii = 0; ii < array_size_1; ii++) {
//			if (center_run_num[ii] == run_num[0]) {
			if (center_run_num[ii] == run_n) {
				pos = ii;
	//			cout << " sadfsdaf" << endl;
			}
		}
	
		// cout << "sdafdsafsa " << pos << endl;
	
		
		Double_t coin_center = center_mean[pos];
	
	
		// cout << "sdafdsafsa " << pos << "    " << coin_center << endl;
		
	
	
	
	
	
		///*--------------------------------------------------*/
		// Read in Nturple files
		
		TString data_file_dir("data/");
		TString data_file;
	
		// data_file.Form( "coin%i.root", run_num[0]);
		data_file.Form( "coin%i.root", run_n);
		data_file.Insert(0, data_file_dir);
		
	//	data_file_dir = "sdaf";
		
	//	cout <<  data_file << endl;
	
		TFile* file_in = TFile::Open(data_file);
	
		file_in->ls();
	
	    TTree* t1 = (TTree*) file_in->Get("h9500");
	
	//	t1->Show(100);
	
		Float_t t_missmass, t_W, t_t, t_Q2, t_th_pq, t_phi_pq, t_PmPar, t_PmPer, t_PmOop;
	
		t1->SetBranchAddress("missmass", &t_missmass);
		t1->SetBranchAddress("W", 		 &t_W       );
		t1->SetBranchAddress("t",        &t_t       );
		t1->SetBranchAddress("Q2",       &t_Q2      );
		t1->SetBranchAddress("th_pq",    &t_th_pq   );
		t1->SetBranchAddress("phi_pq",   &t_phi_pq  );
		t1->SetBranchAddress("Pmpar",    &t_PmPar   );
		t1->SetBranchAddress("Pmper",    &t_PmPer   );
		t1->SetBranchAddress("Pmoop",    &t_PmOop   );
	
	
	
	 	Float_t t_hsdelta, t_hsytar, t_hsyptar, t_hsxptar, t_ssdelta, t_ssytar, t_ssyptar, t_ssxptar; 
	
		t1->SetBranchAddress("hsdelta", &t_hsdelta);
		t1->SetBranchAddress("hsytar",  &t_hsytar );
		t1->SetBranchAddress("hsyptar", &t_hsyptar);
		t1->SetBranchAddress("hsxptar", &t_hsxptar);
		t1->SetBranchAddress("ssdelta", &t_ssdelta);
		t1->SetBranchAddress("ssytar",  &t_ssytar );
		t1->SetBranchAddress("ssyptar", &t_ssyptar);
		t1->SetBranchAddress("ssxptar", &t_ssxptar); 
	
	
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/// From the kumac 
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
	 
	 
	
	
		t1->Draw("missmass >> hist_missmass","missmass > -0.075 && missmass < 0.125", "goff");
	
		TString Q2_limit;
		Q2_limit.Form("Q2 < %f && Q2 > %f", Q_2[0]*1.5, Q_2[0]*0.5);
		t1->Draw("Q2 >> hist_q2", Q2_limit, "goff");
	
	
	//	t1->Draw("W >> hist_W"," W > 2 && W < 2.40", "goff");
	
	
	
	
	
	
	//	initial_plot();
	
	
	
	//	TCanvas*c1 = new TCanvas("c1");
	//    TH1F *hnew = (TH1F*)gDirectory->Get("hist_missmass");
	//	hnew->Draw();
	
	
		TH1F * testtest = initial_plot(t1, "missmass", -0.075, 0.125);
		TH1F * testtest1 = initial_plot(t1, "W", 2, 2.40);
	
	//	testtest->Draw();
	
	
	
	
	//	t1->Draw("missmass");
	//	t1->Draw("cointime");
	
	
	
	
	
	
	
	
	
	//	t1->Draw("Q2 >> hist_q2", Q2_limit);
	
	
	
	
	
	//	t1->Draw("W",);
	//	t1->Draw("t",);
	
	//	t1->Draw("Q2",);
	
	//	t1->Draw("th_pq",);
	
	
	//  missmass, W, t, Q2, th_pq, phi_pq, PmPar, PmPer, PmOop
	// 
	// 	t1->Branch("missmass", &missmass);
	// 	t1->Branch("W", &W);
	// 	t1->Branch("t", &t);
	// 	t1->Branch("Q2", &Q2);
	// 	t1->Branch("th_pq", &th_pq);
	// 	t1->Branch("phi_pq", &phi_pq);
	// 	t1->Branch("PmPar", &PmPar);
	// 	t1->Branch("PmPer", &PmPer);
	// 	t1->Branch("PmOop", &PmOop);
	// 
	
	
	
	//   cuts $66 $10.and.$20.and.$30.and.$41.and.$45.and.$50.and.$60
	//   cuts $67 $10.and.$20.and.$30.and.$42.and.$45.and.$50.and.$60
	//   cuts $68 $10.and.$20.and.$30.and.$43.and.$45.and.$50.and.$60
	 
		
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	//	Real and random coincidence cuts
	//
	//	TCut coin_cut;
	//	coin_cut.Form("coin_cut -6/7 ", );
	//	TCut rand_cut 
	
	
	//	TCut cointime_cut_real = "cointime + 6.8 > -1 && cointime + 6.8 < 0";	
	//	TCut cointime_cut_rand = "cointime + 6.8 >  0 && cointime + 6.8 < 7"; 
		
	// 	TCut real_coin_cut = HMS_cut() && SOS_cut() && PID_cut() && Missingmass_cut(miss_mass_offset[int(kset[0])]) && Cointime_primary_cut(coin_center);
	// 
	// 	TCut rand_coin_cut = HMS_cut() && SOS_cut() && PID_cut() && Missingmass_cut(miss_mass_offset[int(kset[0])]) && Cointime_random_cut(coin_center);
	
	
	
	//	cout << HMS_cut() << endl;
	
	
	
	
	
	//	t1->Draw("W", hms_cut());
	//	t1->Draw("missmass");
	
	//	t1->Draw("cointime", real_coin_cut);
	
	
	//	t1->Draw("cointime+6.8", real_coin_cut);
	//	t1->Draw("cointime+6.8", rand_coin_cut);
	
	
	
	// 	Cointime_primary_cut(coin_center);
	// 	Cointime_random_cut(coin_center);
	
	
	//	t1->Draw("cointime", Cointime_primary_cut(coin_center));
	//	t1->Draw("cointime", Cointime_random_cut(coin_center));
	//	t1->Draw("cointime");
	
	
	///*--------------------------------------------------*/
	///*--------------------------------------------------*/
	//  
	//    cuts $66 $10.and.$20.and.$30.and.$41.and.$45.and.$50.and.$60
	//    cuts $67 $10.and.$20.and.$30.and.$42.and.$45.and.$50.and.$60
	//    cuts $68 $10.and.$20.and.$30.and.$43.and.$45.and.$50.and.$60
	
	
	//	TCut primary_cut = HMS_cut() && SOS_cut() && PID_cut() && Missingmass_cut(miss_mass_offset[int(kset[0])]) && Diamond_cut() && cointime_cut_real;
	
	//	TCut real_coin_cut = HMS_cut() && SOS_cut() && PID_cut() && Missingmass_cut(miss_mass_offset[int(kset[0])]) && Cointime_primary_cut(coin_center) && Set_t_limit(tmnset[0], tmxset[0]);
	
	
	//	TCut real_coin_cut = HMS_cut() && SOS_cut() && PID_cut() && Cointime_primary_cut(coin_center) && Set_t_limit(tmnset[0], tmxset[0]) && Diamond_cut() && Missingmass_cut(miss_mass_offset[int(kset[0])]) ;
	
	
	
		TString all_coin_cut;
		all_coin_cut = HMS_cut() + " && "+ SOS_cut() + " && " + PID_cut() + " && " + Cointime_all(coin_center) + " && " + Diamond_cut() + " && " + Missingmass_cut(miss_mass_offset[int(kset[0])]) + " && " + Set_t_limit(tmnset[0], tmxset[0]);
	
	
	
		TString real_coin_cut;
		real_coin_cut = HMS_cut() + " && "+ SOS_cut() + " && " + PID_cut() + " && " + Cointime_primary_cut(coin_center) + " && " + Diamond_cut() + " && " + Missingmass_cut(miss_mass_offset[int(kset[0])]) + " && " + Set_t_limit(tmnset[0], tmxset[0]);
	
	
		TString rand_coin_cut; 
		rand_coin_cut = HMS_cut() + " && "+ SOS_cut() + " && " + PID_cut() + " && " + Cointime_random_cut(coin_center) + " && " + Diamond_cut() + " && " + Missingmass_cut(miss_mass_offset[int(kset[0])]) + " && " + Set_t_limit(tmnset[0], tmxset[0]);
	
	
	//	cout << real_coin_cut << endl;
	
	//	TString teeee = real_coin_cut;
	//	cout << teeee << endl;
	
	
	
		TH1F* coin_all = new TH1F("all", "all", 70, -4, 10);
	
	
		TH1F* coin_real = new TH1F("real", "real", 70, -4, 10);
		TH1F* coin_rand = new TH1F("rand", "rand", 70, -4, 10);
	
	
	
		t1->Draw("cointime>>all", all_coin_cut, "goff");
		t1->Draw("cointime>>real", real_coin_cut, "goff");
		t1->Draw("cointime>>rand", rand_coin_cut, "goff");
	
	


		c1->cd();
		c1->Update();
		
		coin_all->Draw("hist");
		coin_all->SetLineColor(1);
		coin_real->Draw("same");
		coin_real->SetLineColor(4);
		coin_rand->Draw("same");
		coin_rand->SetLineColor(2);
		
	
		
		c2->Update();	
		
		
		c2->cd(1);
	
		TH1F* mm = new TH1F("mm", "mm", 50, -0.1, 0.12);
		TH1F* mm_1 = new TH1F("mm1", "mm1", 50, -0.1, 0.12);
		
		t1->Draw("missmass-0.939565 >> mm", real_coin_cut, "goff");
		t1->Draw("missmass-0.939565 >> mm1", rand_coin_cut, "goff");
	
	
	//  	mm->SetMarkerStyle(5);
	// // 	mm->Draw("PE");
	//  	
	//  	mm1->SetMarkerStyle(4);
	// // 	mm1->Draw("Psame");
	 
	
		TH1F* mm_diff = (TH1F*) mm->Clone();
	
		mm_diff->Add(mm1, -0.3333333);
	 	mm_diff->SetMarkerStyle(3);
		mm_diff->Draw("PE");
	
	
	
	
		TString miss_mass_offset_str;
		miss_mass_offset_str.Form("missmass-0.939565-%f", miss_mass_offset[int(kset[0])]);
	
		TH1F* mm_off = new TH1F("mm_off", "mm_off", 50, -0.1, 0.12);
		TH1F* mm_1_off = new TH1F("mm1_off", "mm1_off", 50, -0.1, 0.12);
	
		t1->Draw(miss_mass_offset_str + " >> mm_off", real_coin_cut, "goff");
		t1->Draw(miss_mass_offset_str + " >> mm1_off", rand_coin_cut, "goff");
	
	
	
		TH1F* mm_offset_diff = (TH1F*) mm_off->Clone();
	
		mm_offset_diff->Add(mm_1_off, -0.3333333);
	 	mm_offset_diff->SetMarkerStyle(3);
	
	
		c2->cd(2);
		mm_offset_diff->Draw("PE");
	
	
	
		c2->Update();
	
	
	 	c3->cd();
		c3->Update();	
	
	
		///*--------------------------------------------------*/
		//// Methods to access the TGraph from tree
		//
		// TGraph *gr = (TGraph *)gPad->GetPrimitive("Graph")->Clone(); // 2D
		// TGraph *gr = new TGraph(tree->GetSelectedRows(), tree->GetV2(), tree->GetV1());
	
	
	
		TString without_diamond_cut;
	
		without_diamond_cut = HMS_cut() + " && "+ SOS_cut() + " && " + PID_cut() + " && " + Missingmass_cut(miss_mass_offset[int(kset[0])]) + " && " + Set_t_limit(tmnset[0], tmxset[0]);
	
	
		t1->Draw("W:Q2", without_diamond_cut);
	
		TGraph *gr1 = (TGraph*)gPad->GetPrimitive("Graph")->Clone();
	
		gr1->SetMarkerStyle(7);	
		gr1->SetMarkerSize(3);	
		gr1->SetMarkerColor(1);	
	
	
	
		t1->Draw("W:Q2", real_coin_cut);
		TGraph *gr2 = (TGraph* )gPad->GetPrimitive("Graph")->Clone();
	
	//	t1->Draw("W:Q2", real_coin_cut, "goff");
		
	
	
		gr2->SetMarkerStyle(7);	
		gr2->SetMarkerSize(1);	
		gr2->SetMarkerColor(2);	
	
	
	
		gr1->Draw("AP");
		gr2->Draw("P");


		
		TString root_out_dir;
		root_out_dir.Form("%i", run_n);

		file_out->cd();
		TDirectory*asd = file_out->mkdir(root_out_dir);
		asd->cd();



		mm_1_off->Write();
		mm_offset_diff->Write();

		c1->Write("Cointime");
		c2->Write("mm");
		c3->Write("diamond");


// 		run_num_out = run_num;
// 		run_random_out = gRandom->Gaus(0,1);
// 	
	
//      tw=([tmax]-[tmin])/[NBt]
//      tc=[tmin]+([j]-0.5)*[tw]
//      ve/cr tb([NBphi]) r [NBphi]*[tc]


		acccc_temp = acccc_temp + (charge[i] * tot_eff[i]);

//		cout << yield_real << "     " << yield_rand  <<  "      " << acccc_temp <<  "     " << yield_yield/acccc_temp << endl;


//		cout << "PPPPP Acc : " << acccc_temp << endl;



		file_out->cd();

		/*--------------------------------------------------*/
		/// Calculate phi and t bins

		for(int t_bin_set_tmp = 1; t_bin_set_tmp <= t_bin_set; t_bin_set_tmp++) {
			
			float t_c = t_min + (t_bin_set_tmp-0.5) * t_width;
			

			/*--------------------------------------------------*/
			/// calculate yield and error


        	//  cut $71 ([tmin]+([tmax]-[tmin])/[NBt]*([j]-1))<t<([tmin]+([tmax]-[tmin])/[NBt]*[j])

			// root_out_dir.Form("%i", run_n);

			TString yield_cut;

			float yield_cut_t_l = t_min + (t_max - t_min)/t_bin_set * (t_bin_set_tmp-1);
			float yield_cut_t_h = t_min + (t_max - t_min)/t_bin_set * t_bin_set_tmp;

			yield_cut.Form("t > %f && t < %f",yield_cut_t_l, yield_cut_t_h);


			TH1F* event_phi_real = new TH1F("phi_real", "phi_real", phi_bin_num, 0, 360);
			t1->Draw("phi_pq*180/3.1415926 >> phi_real", yield_cut + "&&" + real_coin_cut);
//			event_phi_real->Write();
				

			TH1F* event_phi_rand = new TH1F("phi_rand", "phi_rand", phi_bin_num, 0, 360);
			t1->Draw("phi_pq*180/3.1415926 >> phi_rand", yield_cut + "&&" + rand_coin_cut);
//			event_phi_rand->Write();

			
//			TH1F* event_real_diff = (TH1F*) event_phi_real->Clone();
//			event_real_diff->Add(event_phi_rand, -0.3333333333);
//			event_real_diff->Write("asdasd");
			event_phi_real->Add(event_phi_rand, -0.3333333333);

			if (t_bin_set_tmp==1) {			
				real_event[1]->Add(event_phi_real);	
			} else if (t_bin_set_tmp==2) {
				real_event[2]->Add(event_phi_real);	
			}
			

//			cout << yield_cut_t_l << "    " << yield_cut_t_h << endl;	

//			t1->Draw(">>elist", yield_cut + "&&" + real_coin_cut);
//			Int_t yield_real = elist->GetN();

//			yield_real = htemp->GetN();

//		 	TH1F* htemp = (TH1F*) gROOT->FindObject("htemp");
//			htemp->SetTitle(yield_cut);
//			htemp->Write();
//			delete htemp;

//			t1->Draw(">>elist", yield_cut + "&&" + rand_coin_cut);
//			Int_t yield_rand = elist->GetN();

//			yield_rand = htemp->GetN();

//		 	TH1F* htemp = (TH1F*) gROOT->FindObject("htemp");
//			htemp->SetTitle(yield_cut);
//			htemp->Write();
//			delete htemp;


// 			float yield_yield = (float(yield_real) - float(yield_rand) * 0.3333333333) * 0.463;
// 
// 			cout << yield_cut << endl;	
// 			cout << yield_real << "     " << yield_rand << "    " << yield_yield / acccc_temp/1000. << endl;
// 	
 




//   	     	cout << "aaaaaaaaaaaa  " << (1. / charge[i] * tot_eff[i]) / 1000. <<  endl; 


			for (float phi_step_tmp = pmn; phi_step_tmp <= pmx; phi_step_tmp = phi_step_tmp + pstp) {
				// cout << phi_step_tmp << endl;


				int phi_itt = (phi_step_tmp + pstp/2) / 22.5;


//				yield = yield_yield / acccc_temp/1000.;

				phi = phi_step_tmp;			
				// tb = phi_step_tmp * t_c;
// 				tb = t_c;
// 
// 				int yield_real = event_phi_real->GetBinContent(phi_itt);
// 				int yield_rand = event_phi_rand->GetBinContent(phi_itt);
// 

// 				//cout << phi_itt << "   " << yield_real << "  " << yield_rand << "  " << phi << "  " << tb << endl;
// 
// 				float yield_yield = (float(yield_real) - float(yield_rand)*0.3333333333) * 0.463 / acccc_temp / 1000.;
// 
// 				// cout  << yield_real << "   " << yield_rand << endl;
// 
// 				yield = yield_yield;
// 							
// 				cout << yield  << "    " << t_bin_set_tmp << "    " <<  phi << "   " << tb << endl;
//				cout << phi_itt << "    "<< yield  << "    " << t_bin_set_tmp << "    " <<  phi << "   " << tb << endl;

				tree_out->Fill();

			}

			delete event_phi_real;
			delete event_phi_rand;
//			delete event_real_diff; 

		}
	
	} 
	
	
	//	t1->Draw("missmass-0.939565", real_coin_cut);

	

	for (int iii=1; iii < 17; iii++ ) {

		float yield0 = real_event[1]->GetBinContent(iii)/acccc_temp/1000. *0.463;
		float yield1 = real_event[2]->GetBinContent(iii)/acccc_temp/1000. *0.463;
		cout << yield0 << "    " << yield1 << endl;

	}




	file_out->cd();
	tree_out->Write();
	
	
	delete [] run_num;
	

}





TH1F* initial_plot(TTree*t1, TString name, float dn_limit, float up_limit) {


//	cout << "asdf" << endl;



	TString cut_cut;
	cut_cut.Form(" %s > %f && %s < %f", name.Data(), dn_limit, name.Data(), up_limit);

	TString plot_opt;
	plot_opt.Form("%s >> hist_%s", name.Data(), name.Data());
//	plot_opt.Form("%s", name.Data());

//	cout << cut_cut << "   " << name << endl;

	t1->Draw(plot_opt, cut_cut, "goff");

//	t1->Draw("W >> hist_W"," W > 2 && W < 2.40", "goff");


	TString plot_name;
	plot_name.Form("hist_%s", name.Data());

	TH1F *hnew = (TH1F*)gDirectory->Get(plot_name);

	//hnew->Draw();

	return hnew;
} 


// 
// /*--------------------------------------------------*/
// // macro get_yields_exp pol Q2 epsil th_pq tmin tmax NBt kin itar
// 
// 
// void get_yields_exp(Double_t Q_2, Double_t epsilon, Double_t th_pq, Double_t t_min, Double_t t_max) {
// 
// 	cout << "Calculate the yield" << endl;
// 
// 
// 
// 
// 
// }
// 
// 




// /*--------------------------------------------------*/
// Define Cuts
// Define HMS Cuts

/*--------------------------------------------------*/
// Kumac code HMS cut
//    cuts $11 abs(hsdelta)<8.0
// *   cuts $12 abs(hsxptar)<0.060  | cut on Jochen's thesis p.53
//    cuts $12 abs(hsxptar)<0.080   | extra sieve slit fitting should improve things
//    cuts $13 abs(hsyptar)<0.035
// 
//    cuts $10 $11.and.$12.and.$13


TString HMS_cut() {

	TString hsdelta_cut = "abs(hsdelta) < 8.0";

//	TCut* hsdelta_cut = ;

	TString hsxptar_cut = "abs(hsxptar) < 0.080";
	TString hsyptar_cut = "abs(hsyptar) < 0.035";

	TString hms_cut = hsdelta_cut + " && " + hsxptar_cut + " && " + hsyptar_cut;

 	return hms_cut;

}






// /*--------------------------------------------------*/
// Define SOS Cuts
//    cuts $21 abs(ssdelta)<15                       
//    cuts $22 abs(ssxfp)<20	
//    cuts $23 abs(ssxfp+ssxpfp*313.01+5.64)<50.
//    cuts $24 abs(ssyfp+ssypfp*313.01)<30.
//    cuts $25 ssytar<1.5          | cut on Jochen's thesis p.80, for th_SOS>42deg.
//    cuts $26 abs(ssxptar)<0.04
//    cuts $27 abs(ssyptar)<0.065                 



TString SOS_cut() {

	TString ssdelta_cut = "abs(ssdelta) < 15.";                  
	TString ssxfp1_cut  = "abs(ssxfp) < 20.";	
	TString ssxfp2_cut  = "abs(ssxfp+ssxpfp*313.01+5.64) < 50."; 
	TString ssyfp_cut   = "abs(ssyfp+ssypfp*313.01) < 30.";
	TString ssytar_cut  = "ssytar < 1.5";       
	TString ssxptar_cut = "abs(ssxptar) < 0.04";
	TString ssyptar_cut = "abs(ssyptar) < 0.065";               

	TString t1 = "ssyptar > (-0.125+0.00425*ssdelta+0.064*ssytar-0.0017*ssdelta*ssytar)";               
	TString t2 = "ssyptar < (0.125-0.00425*ssdelta+0.064*ssytar-0.0017*ssdelta*ssytar)";               

	TString sos_cut = ssdelta_cut + " && " + ssxfp1_cut + " && " +  ssxfp2_cut + " && " + 
				   ssyfp_cut  + " && " + ssytar_cut + " && " + ssxptar_cut + " && " + 
				   ssyptar_cut + " && " + t1 + " && "+ t2;

   	return sos_cut; 

}

// /*--------------------------------------------------*/
// Define PID Cuts
// 

//    cuts $31 hsbeta>0.95
//    cuts $32 (haero_po+haero_ne)>3.0
//    cuts $33 hcer_npe<2.
//    cuts $34 scer_npe>0.5
//    cuts $35 ssshtrk>0.70           
// 
//    cuts $30 $31.and.$32.and.$33.and.$34.and.$35



TString PID_cut() {

	TString hsbeta_cut    = "hsbeta > 0.95"; 
	TString haero_cut     = "(haero_po + haero_ne) > 3.0";
	TString hcer_npe_cut  = "hcer_npe < 2.";
	TString scer_cut      = "scer_npe > 0.5";
	TString ssshtrk_cut   = "ssshtrk > 0.70";

	TString pid_cut = hsbeta_cut + " && " + haero_cut + " && " + hcer_npe_cut + " && " + scer_cut + " && " + ssshtrk_cut;


	return pid_cut;

}




TString Diamond_cut() {

// ----------------------------------------------------*/
// Kumac diamand cuts
//
//    cuts $51 q2>(-2.916*(w-2.4733)+1.5846)
//    cuts $52 q2>(-4.452*(w-1.94)+3.4199) 
//    cuts $53 q2<(-2.677*(w-2.3308)+2.231)
//    cuts $54 q2<(-4.3748*(w-2.41)+1.78)

	TString diamond_cut_1; 
	diamond_cut_1.Form("Q2 > (-2.916*(W-2.4733)+1.5846)");

	TString diamond_cut_2;
	diamond_cut_2.Form("Q2 > (-4.452*(W-1.94)+3.4199)");

	TString diamond_cut_3;
	diamond_cut_3.Form("Q2 < (-2.677*(W-2.3308)+2.231)");

	TString diamond_cut_4;
	diamond_cut_4.Form("Q2 < (-4.3748*(W-2.41)+1.78)");

	TString combine_cut;

 	TString cut1 = diamond_cut_1 + " && " + diamond_cut_2 + " && " + diamond_cut_3 + " && "+ diamond_cut_4; 

	return cut1;

}





TString Missingmass_cut(Double_t m_m_offset) {


//	cout << "ssssssssss : "<< m_m_offset << endl;


	TString cut_tmp; 

	cut_tmp.Form("(missmass-%f) > 0.875  && (missmass-%f) < 0.98", m_m_offset, m_m_offset);

//	TCut missmass_cut    = "0.875 < missmass"; 
//	$45 0.875<missmass-[mmof]<1.03;

	TString missmass_cut = cut_tmp; 

	return missmass_cut;

}


TString Set_t_limit(Double_t tmin, Double_t tmax) {

	TString cut_tmp; 
	cut_tmp.Form("t > %f && t < %f", tmin, tmax);

	TString t_limit = cut_tmp; 

	return t_limit;
}

// Double_t Get_center(Double* tar_arr, Int_t pos) { 
// 
// 
// 
// 
// 
// 
// 
// 
// }
// 


// Int_t Get_center(Double* tar_arr, Int_t pos) { 
// 
// 
// 
// 
// 
// 
// 
// 
// }
// 

TString Cointime_primary_cut(Double_t center) {
//    cuts $41 -1.<cointime-[center]<1.
	TString cut_tmp;
	cut_tmp.Form("-1 < (cointime - %f) && (cointime - %f) < 1", center, center);
// 	TString cut1 = cut_tmp.Data(); 

//	cout << cut_tmp << endl;
	return cut_tmp;
}




TString Cointime_random_cut(Double_t center) {

//    cuts $42 1.<cointime-[center]<7.

	TString cut_tmp;
	cut_tmp.Form("1 < cointime - %f && cointime - %f < 7", center, center);
 	TString cut1 = cut_tmp.Data(); 

	return cut1;
}

TString Cointime_all(Double_t center) {

//    cuts $42 1.<cointime-[center]<7.

	TString cut_tmp;
	cut_tmp.Form("-6 < cointime - %f && cointime - %f < 10", center, center);
 	TString cut1 = cut_tmp.Data(); 

	return cut1;
}


//	int sizeOfArray = 0;
//	while (center_run_num_d[sizeOfArray] != '\0') { //loop until NULL character is found
//    	sizeOfArray++; //increment size
//	}
//	cout <<  sizeOfArray << "     " << endl;



Double_t SizeOfArray(Double_t* tar_arr) {

	int size_of_array = 0;

// 	while (tar_arr[size_of_array]) { //loop until NULL character is found
// 		size_of_array++; //increment size
// 		cout << "PPPP : "<< tar_arr[size_of_array] << endl;
// 	}

// 	while(true) {
// 
//   		if( (tar_arr[size_of_array]+size_of_array) == '\0') // returns 5...
//     		break;
// 
//  		size_of_array++; //increment size
// 
//  		cout << "PPPP : "<< tar_arr[size_of_array] << endl;
//  	 }
// 
// 	for(int* it = tar_arr.begin(); it!=tar_arr.end(); ++it){
//   		cout<<*it<<endl;
// 	}
// 

// 	for (int& x :tar_arr) {
// 		cout << x << endl;
// 	}


	while( int (*tar_arr) > 0 ) {
//	    cout<<*tar_arr<<endl;
	    tar_arr++;
		size_of_array++;
	}

//	cout << size_of_array << endl;

	return size_of_array;

}





void Clear_Up() {
	delete [] run_num;
}

