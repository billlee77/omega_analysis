
const Float_t m_p     = 0.93827;
const Float_t m_omega = 0.78259;

const Float_t Mpi0    = 134.9766;
const Float_t Mpi02   = Mpi0**2;
const Float_t pi      = 3.141592653589793; 

Int_t neg_count;


TString q2_set_str;

/*--------------------------------------------------*/
/*--------------------------------------------------*/

void wt_test() {

	TString q2;
	TString eps;
	TString theta;

	TString q2_setting[2]     = {"160", "245"};
	
	TString eps_160[2]        = {"32", "59"};
	TString eps_245[2]        = {"27", "55"};
	
// 	TString hms_angle_160_l[2]  = {"+0970", "+3000"};
// 	TString hms_angle_160_h[3]  = {"-2730",  "+0000", "+3000"};
// 	
// 	TString hms_angle_245_l[2]  = {"+1350", "+3000"};
// 	TString hms_angle_245_h[3]  = {"-3000", "+0000", "+3000"};


	TString hms_angle_160_l[2]  = {"+1", "+3"};
	TString hms_angle_160_h[3]  = {"-3",  "00", "+3"};
	
	TString hms_angle_245_l[2]  = {"+1", "+3"};
	TString hms_angle_245_h[3]  = {"-3", "00", "+3"};


	TString dir_str;
	TString file_name;

	dir_str = "../";


	for (int i = 0; i < 2; i++) {

		q2 = q2_setting[i];

		if( q2 =="160") {

			for (int dd = 0; dd < 2; dd++) {

				eps = eps_160[dd];

				if ( eps.Atoi()  < 40) {
					
					for (int iii = 0; iii < 2; iii++) {

						theta = hms_angle_160_l[iii]; 

						file_name = "omega_" + q2 + "_" + eps + "_" + theta + ".root";

						cout << file_name  << endl;

						wt_test_single(dir_str + file_name);

					}

				} else {
 
					for (int iii = 0; iii < 3; iii++) {
				
						theta = hms_angle_160_h[iii]; 

						file_name = "omega_"+ q2 + "_" + eps + "_" + theta + ".root";

						cout << file_name  << endl;

						wt_test_single(dir_str + file_name);

					}

				}

			}

		}	else {

			for (int ii = 0; ii < 2; ii++) {

				eps = eps_245[ii];

				if (eps.Atoi() < 40) {

					for (int iii = 0; iii < 2; iii++) {

						theta = hms_angle_245_l[iii]; 

						file_name = "omega_"+ q2 + "_" + eps + "_" + theta + ".root";

						cout << file_name  << endl;

						wt_test_single(dir_str + file_name);

					}
 
				} else {

					for (int iii = 0; iii < 3; iii++) {

						theta = hms_angle_245_h[iii]; 

						file_name = "omega_"+ q2 + "_" + eps + "_" + theta + ".root";

						cout << file_name  << endl;

						wt_test_single(dir_str + file_name);

					}

				}			

			}
			
  		}

	}	


//	wt_test_single("omega_160_32_+1.root");
//	wt_test_single("omega_test");


}















void wt_test_single(TString root_file_str) {

	
//	TString root_file_str = "omega_160_32_+1.root";

//	TString root_file_str = "omega_test";

// 	TFile* f1     = new TFile(root_file_str + ".root", "READ");
// 	TFile* f1_mod = new TFile(root_file_str + "_mod.root", "RECREATE");

//	TFile* f1     = new TFile(root_file_str + ".root", "UPDATE");

	TFile* f1     = new TFile(root_file_str, "UPDATE");
	
	TTree* t1 = (TTree*)f1->Get("h666");

	cout << root_file_str << "    " << root_file_str.ReplaceAll(".root", "") << endl; 
	
	TFile* f1_mod = new TFile(root_file_str.ReplaceAll(".root", "") + "_mod.root", "RECREATE");

	neg_count = 0;	

//	exit(0);


//	TBranch* b1 =

//	t1->Draw("Weight");
	

	Float_t weight_mod;
	Float_t sig_mod;
	Float_t weight;
	Float_t sigcm;

	Float_t thetacm, thetapq, phicm, t, w, Q2, eps;
	Float_t wi, ti, Q2i, thetacmi, phipqi, epsiloni, nu, q;
	Float_t phicmi, wcmi, tprimei, ti, ui;



//	TBranch* w_mod = t1->Branch("mod", &weight_mod, "mod/F");

	TBranch *b = t1->GetBranch("Weight_mod");
	t1->GetListOfBranches()->Remove(b);
//	b->Delete();



	TBranch* w_mod = t1->Branch("Weight_mod", &weight_mod, "Weight_mod");

	t1->SetBranchAddress("Weight",  &weight  );
	t1->SetBranchAddress("sigcm",   &sigcm   );

	t1->SetBranchAddress("thetacm", &thetacm );

	t1->SetBranchAddress("phipq",   &phicm   );
	t1->SetBranchAddress("t",       &t       );
	t1->SetBranchAddress("W",       &w       );

	t1->SetBranchAddress("Q2",      &Q2      );
	t1->SetBranchAddress("epsilon", &eps     );

	t1->SetBranchAddress("Wi",       &wi       );
	t1->SetBranchAddress("Q2i",      &Q2i      );
	t1->SetBranchAddress("thetacmi", &thetacmi );
	t1->SetBranchAddress("phipqi",   &phipqi   );

	t1->SetBranchAddress("nu",       &nu       );
	t1->SetBranchAddress("q",        &q        );


	t1->SetBranchAddress("epsiloni", &epsiloni );
	t1->SetBranchAddress("phicmi",   &phicmi   );
	t1->SetBranchAddress("wcmi",     &wcmi     );
	t1->SetBranchAddress("tprimei",  &tprimei  );

	t1->SetBranchAddress("ti",       &ti  );
	t1->SetBranchAddress("ui",       &ui  );


//	t1->SetBranchAddress("sigcm", &sigcm);
//	t1->SetBranchAddress("sigcm", &sigcm);
//	t1->SetBranchAddress("sigcm", &sigcm);
//	t1->SetBranchAddress("sigcm", &sigcm);

    Long64_t nentries = t1->GetEntries(); // read the number of entries in the t3





//	t1->SetBranchStatus("mod", 1);

    for (Long64_t i = 0; i < nentries; i++) {
//    for (Long64_t i = 0; i < 1; i++) {

		
//		t1->SetBranchStatus("Weight", 1);

		t1->GetEntry(i);

		Float_t Wsq, qsq, m_psq, invm, tt, uu, theta, phi;
		Float_t mass;

		
		/*--------------------------------------------------*/

// 		invm  = w;
// 		tt    = t;
// //		qsq   = Q2;
// 		theta = thetacm;
// 	 	phi   = phicm;

// 		/*--------------------------------------------------*/
 

		invm  = wcmi;
		tt    = ti;
        qsq   = Q2i;
		theta = thetacmi;
		phi   = phicmi;
		eps   = epsiloni;
		uu    = ui;

// 

//		invm = sqrt((nu+m_p)**2-q**2);
	
		/*--------------------------------------------------*/
		/// Constant		
		Wsq   = invm*invm;
		m_psq = m_p*m_p;
		mass  = m_omega;







		Float_t e_photCM = (Wsq - qsq - m_psq)/invm/2.;
        Float_t e_omCM   = (Wsq + mass**2 - m_psq)/invm/2.;
      
        Float_t t_min    = -qsq + mass**2 -2.*(e_omCM*e_photCM-sqrt((e_omCM**2-mass**2)*(e_photCM**2+qsq)));

//		Float_t tprime   = abs(tt)-abs(t_min);

		Float_t tprime   = tprimei;


//		cout<< "Check tprime:  " << -tt << "     " << tprime << endl;


		/*--------------------------------------------------*/	
		// Cal_New_Mod(thetacm, phicm, t_gev, tprime_gev, q2_gev, w_gev, eps);		/// Calculate the new Cross section
//		 sig_mod = Cal_New_Mod(theta, phi, -tt, tprime, qsq, invm, eps);
//		 sig_mod = sig_gmh_mod(theta, phi, -tt, tprime, qsq, invm, eps);



 		if (root_file_str.Contains("160")) {
 
			q2_set_str = "160";

 		} else if (root_file_str.Contains("245")) {

			q2_set_str = "245";
 
 		} else {
 
// 			sig_mod = 0.0;
 			exit(0);
// 	
 		}

//     /*--------------------------------------------------*/
//     //  Itt_20 Fitted parameterization
// 
// 	if( q2_set_str == "160") {
// 		cout << "Q2 = 1.6 Parameterization is used ! " << endl; 
// 	} else if (q2_set_str == "245") {
// 		cout << "Q2 = 2.45 Parameterization is used !" << endl; 
// 	} else {
// 		cout << "No parameterization found !" << endl;
// 		exit(0);
// 	}

		sig_mod = sig_gmh_mod_u_two_model(theta, phi, -uu, qsq, invm, eps);

		if (sig_mod < 0.0) {

// 			cout << "The sig_mod calculated is below 0 !!!!!   " << sig_mod << endl;
// 			cout << "If you see this message, go back to your fit_in_t and then refit!" << endl;
// 			cout << "Parameterization should not yield sig < 0" << endl;
//			exit(0);

			neg_count++;

//			sig_mod = 0.0;

		}

//		sig_mod = sig_gmh_mod_u(theta, phi, -uu, qsq, invm, eps);

		/*--------------------------------------------------*/	
		/// Calculate the new weight
	

		if( weight*sig_mod/sigcm != weight*sig_mod/sigcm) {

			weight_mod = weight;			
//			cout << "!!! " << weight_mod << endl;	

		} else {

			weight_mod = weight * sig_mod/sigcm;

		}

//		weight_mod = sig_mod / sigcm;

// 		cout << weight << "     " << weight_mod << "    " << sigcm << "    " << sig_mod<< endl;

//		exit(0);

// 		if (Float_t(sigcm/sig_mod) != Float_t(1.0)) {
// 
//  			cout << "  Sig: " << sigcm << "    Sig_mod:" << sig_mod << "  "  << sigcm / sig_mod << endl;
// 			
// 
// 		} 
// 	

//		exit();

//		weight_mod = 1.1;

//		t1->SetBranchStatus("*", 0);
//		t1->SetBranchStatus("mod", 1);

		w_mod->Fill();

//		t1->SetBranchStatus("*", 1);


    }






// 	t1->SetBranchAddress("Weight_mod/F", &weight);
// 
//     for (Long64_t i = 0; i < nentries; i++) {
// 
// 		t1->GetEntry(i);
// 
// 		cout << weight<< "     " << weight_mod << endl;
// 
//     }
//	t1->SetBranchStatus("*",1);
//	t1->Print();

 	cout << "The total number of negtive counts:  " << neg_count << endl;

	/*--------------------------------------------------*/	
	/// Print into a new file 
	
	TTree* t11 = t1->CloneTree();

//	f1->Close();
//	t1->Delete();
//	f1->Clear();


	f1_mod->cd();
//	f1->cd();
//	f1->DeleteAll();
//	gDirectory->Delete("h666;*");
//	f1->ReOpen("RECREATE");

// 	t11->Print();
 	t11->Write();

	f1_mod->Close();
	f1->Close();

}

/*--------------------------------------------------*/
/*--------------------------------------------------*/
/*--------------------------------------------------*/

Float_t sig_gmh_mod_u_two_model(Float_t thetacm, Float_t phicm, Float_t u_gev, Float_t q2_gev, Float_t w_gev, Float_t eps) {

	Float_t sig_gmh_mod;

//	sig_gmh_mod = 100.0;

	Float_t sig, wfactor, Mp, m_p;
	Float_t sigt,sigl,siglt,sigtt;
//	Float_t thetacm,phicm,t_gev,tprime_gev,q2_gev,w_gev,eps,tp;
	Float_t lambda0_sq,lambdapi_sq,alphapi_t,alphapi_0;
	Float_t a,b,c,d,e,f,g,h,fpi,fpi235;
	Float_t m_pi0sq, tp, up;
	Float_t mp = m_p*1000;
	
    Float_t t0, t1, t2, t3 ;
    Float_t l0, l1, l2, l3 ;
    Float_t lt0,lt1,lt2,lt3;
    Float_t tt0,tt1,tt2,tt3;  

	Mp = 938.27231;
	m_p = Mp/1000;


//	cout << endl << "Theta: "<< thetacm << "   phi: " << phicm << "   t_gev: " << t_gev << "   t_prime: " << tprime_gev << "   q2_gev: " << q2_gev  <<  "   w_gev: "<< w_gev  << "  eps: " << eps << endl << endl;   



    m_pi0sq= Mpi02/1.e6;

    up = fabs(u_gev);      // just to make sure it's positive


	///*--------------------------------------------------*/
	/// Checking 
//  	w_gev   =    2.2616900000000002       ;  
//  	q2_gev  =    1.4605500000000000       ; 
//  	up      =   5.00000000000000028E-002  ;
//  	phicm   =   0.39269874999999999       ;	
//  	thetacm =    3.0453998878191055       ; 
//  	eps     =   0.31417277374253511       ; 
	/*--------------------------------------------------*/























/*--------------------------------------------------*/
/// Stright Line
    /*--------------------------------------------------*/
    //  Itt_151 Starting parameterization

  	if( q2_set_str == "160") {
 //		cout << "Q2 = 1.6 Parameterization is used ! " << endl; 
 //		exit(0);
 
 		t0  =        7.73587     ;
 		t1  =        -7.9672     ;
 		t2  =         0.0000     ;
 		t3  =         0.0000     ;
 		l0  =        13.2553     ;
 		l1  =       -47.2633     ;
 		l2  =         0.0000     ;
 		l3  =         0.0000     ;
 		lt0 =        -0.3439        ;
 		lt1 =         5.9217        ;
 		lt2 =         0.0000        ;
 		lt3 =         0.0000        ;
 		tt0 =         8.1221        ;
 		tt1 =      -139.8422        ;
 		tt2 =         0.0000        ;
 		tt3 =         0.0000        ;
 
 	} else if (q2_set_str == "245") {
 //		cout << "Q2 = 2.45 Parameterization is used !" << endl; 
 //		exit(0);

 		t0  =      6.16527     ;
 		t1  =      -4.2124     ;
 		t2  =       0.0000     ;
 		t3  =       0.0000     ;
 		l0  =      12.2546     ;
 		l1  =     -29.8629     ;
 		l2  =       0.0000     ;
 		l3  =       0.0000     ;
 		lt0 =      -0.3620        ;
 		lt1 =       3.1028        ;
 		lt2 =       0.0000        ;
 		lt3 =       0.0000        ;
 		tt0 =      -7.4032        ;
 		tt1 =      63.4705        ;
 		tt2 =       0.0000        ;
 		tt3 =       0.0000        ;
 
 	} else {
 
 		cout << "No parameterization found !" << endl;
 		
 		exit(0);
 
 	}
// 

 	/// Iteration 40 starting parameter
 
//  	t0  =  0.05  ;
//  	t1  =  0.2   ;
//  	t2  =  0.25  ;
//  	t3  =  -0.9  ;
//  	l0  =  0.7   ;
//  	l1  =  -2.7  ;
//  	l2  =  -0.5  ;
//  	l3  =  2.5   ;
//  	lt0 =  0.05  ;
//  	lt1 =  -0.3  ;
//  	lt2 =  -0.15 ;
//  	lt3 =  0.7   ;
//  	tt0 =  0.1   ;
//  	tt1 =  -1.8  ;
//  	tt2 =  -0.3  ;
//  	tt3 =   3    ;

// 
//     /*--------------------------------------------------*/
//     // Sigma T
//       sigt = t0 + t1 * up + t2 * log(q2_gev) + t3 * up * log(q2_gev);

       sigt = t0 / sqrt(q2_gev) + t1 * up / sqrt(q2_gev) + t2 / sqrt(q2_gev) + t3 * up / sqrt(q2_gev);

// 



////     /*--------------------------------------------------*/
////     // Sigma T
//
//  	if( q2_set_str == "160") {
//
// 	  //Double_t f = par[0]*(xx+ par[1])*(xx+ par[1]) + par[2];
//
//      sigt = t0 * (up + t1)*(up + t1) + t2 ;
//
// 	} else if (q2_set_str == "245") {
//
//      // f = par[0]*exp(par[1]*xx) + par[2]
//  
//      sigt = t0 * exp( t1 * up ) + t2;
//
// 	} else {
// 		exit(0);
// 	}





//    sigt = t0;

    /*--------------------------------------------------*/
    // Sigma L
//     sigl = l1* (exp(-l2*(up-l3)) + l2*(up-l3));
//     sigl = l1 + l2 * up;
//     sigl = l1 * exp( l2 * up ) + l3 / up;
//     sigl = l1 * (up-l2)**2 + l3;

//    sigl = l0 + l1 * up + l2 * log(q2_gev) + l3 * up * log(q2_gev);

      sigl = l0/(q2_gev*q2_gev) + l1 * up / (q2_gev*q2_gev) + l2 / q2_gev + l3 * up / q2_gev;


//    sigl = l0 * exp( l1 * up ) + l2;
//      sigl = l0 * (up + l1)*(up + l1) + l2 ;

    /*--------------------------------------------------*/
    // Sigma LT  
//    siglt = (lt0 + lt1 * up + lt2 * log(q2_gev) + lt3 * up * log(q2_gev)) * sin(thetacm);

    siglt = (lt0/q2_gev + lt1 * up / q2_gev + lt2 / q2_gev + lt3 * up / q2_gev) * sin(thetacm);

    //siglt = lt1 * exp( -lt2 * (up-lt3)**2) * sin(thetacm);
     
    /*--------------------------------------------------*/
    // Sigma TT 
//    sigtt = (tt0 + tt1 * up + tt2 * log(q2_gev) + tt3 * up * log(q2_gev)) * sin(thetacm) * sin(thetacm);
//    sigtt = (tt1 * exp( tt2 * up ) + tt3 / up) * sin(thetacm) * sin(thetacm);

    sigtt = (tt0/q2_gev + tt1 * up / q2_gev + tt2 / q2_gev + tt3 * up / q2_gev) * sin(thetacm) * sin(thetacm);


/*--------------------------------------------------*/

//	wfactor=((2.47**2-m_p**2)**2)/((w_gev**2-m_p**2)**2);
	wfactor=1/((w_gev**2-m_p**2)**2);

	sigt = sigt*wfactor;
	sigl = sigl*wfactor;
	siglt= siglt*wfactor;
	sigtt= sigtt*wfactor;
	
	sig = sigt + eps*sigl +eps*cos(2.*phicm)*sigtt +sqrt(2.*eps*(1.+eps))*cos(phicm)*siglt;
	
	sig = sig/2./pi/1.e06;

	sig_gmh_mod = sig;
	
//	cout << "eps=" << eps << " u=" << up << " sigT=" << sigt << " sigL=" << sigl << " sigLT=" << siglt << " sigTT=" << sigtt << " x_mod=" << sig_gmh_mod << endl;
	
//	cout << " wfactor=" << wfactor << " m_p=" << m_p << " w_gev=" << w_gev << endl;


//	exit(0);

	return sig_gmh_mod;



}






















/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// u dependent model original 
///

Float_t sig_gmh_mod_u(Float_t thetacm, Float_t phicm, Float_t u_gev, Float_t q2_gev, Float_t w_gev, Float_t eps) {

	Float_t sig_gmh_mod;

//	sig_gmh_mod = 100.0;

	Float_t sig, wfactor;
	Float_t sigt,sigl,siglt,sigtt;
//	Float_t thetacm,phicm,t_gev,tprime_gev,q2_gev,w_gev,eps,tp;
	Float_t lambda0_sq,lambdapi_sq,alphapi_t,alphapi_0;
	Float_t a,b,c,d,e,f,fpi,fpi235;
	Float_t m_pi0sq, tp, up;
	Float_t mp = m_p*1000;

//	cout << endl << "Theta: "<< thetacm << "   phi: " << phicm << "   t_gev: " << t_gev << "   t_prime: " << tprime_gev << "   q2_gev: " << q2_gev  <<  "   w_gev: "<< w_gev  << "  eps: " << eps << endl << endl;   



    m_pi0sq= Mpi02/1.e6;

    up = fabs(u_gev);      // just to make sure it's positive





	/*--------------------------------------------------*/
    /// Itt 6
	a = 0.0495;
	b = 7;
	c = 0.0;

    sigl = a * exp( -b * up) + c;

	/*--------------------------------------------------*/
    /// Itt 6
	a = 0.088;
	b = 5; 	
	c = 0.00022;

    sigt = a * exp( -b * up) + c;

	/*--------------------------------------------------*/
    /// Itt 6
    a = 0.0015;
    b = 1.04;
    c = 0.000158;

    siglt = a * exp( -b * up) + c;

	/*--------------------------------------------------*/
    /// Itt 6
    a = 0.001;
    b = 1.1789;
    c = 0.00024;

    sigtt = a * exp( -b * up) + c;


	wfactor=((2.47**2-mp**2)**2)/((w_gev**2-mp**2)**2);

	sigl = sigl*wfactor;
	sigt = sigt*wfactor;
	siglt= siglt*wfactor;
	sigtt= sigtt*wfactor;
	
	sig = sigt +eps*sigl +eps*cos(2.*phicm)*sigtt +sqrt(2.*eps*(1.+eps))*cos(phicm)*siglt;
	
	sig = sig/2./pi/1.e06;

	sig_gmh_mod = sig;






	
	return sig_gmh_mod;

}



/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// u dependent model for Q2=1.60
///

Float_t sig_gmh_mod_u_160(Float_t thetacm, Float_t phicm, Float_t u_gev, Float_t q2_gev, Float_t w_gev, Float_t eps) {

	Float_t sig_gmh_mod;

//	sig_gmh_mod = 100.0;

	Float_t sig, wfactor;
	Float_t sigt,sigl,siglt,sigtt;
//	Float_t thetacm,phicm,t_gev,tprime_gev,q2_gev,w_gev,eps,tp;
	Float_t lambda0_sq,lambdapi_sq,alphapi_t,alphapi_0;
	Float_t a,b,c,d,e,f,g,h,fpi,fpi235;
	Float_t m_pi0sq, tp, up;
	Float_t mp = m_p*1000;

    Float_t l1,l2,l3;
    Float_t t1,t2,t3;
    Float_t lt1,lt2,lt3;
    Float_t tt1,tt2,tt3;

//	cout << endl << "Theta: "<< thetacm << "   phi: " << phicm << "   t_gev: " << t_gev << "   t_prime: " << tprime_gev << "   q2_gev: " << q2_gev  <<  "   w_gev: "<< w_gev  << "  eps: " << eps << endl << endl;   

    m_pi0sq= Mpi02/1.e6;

    up = fabs(u_gev);      // just to make sure it's positive

//       /*--------------------------------------------------*/
//       /// Itt_18
//       t1  =   0.1       ;
//       t2  =  -5.0       ;
//       t3  =   0.0002    ; 
//       l1  =   0.05      ;
//       l2  =   0.6       ;
//       l3  =   0.02      ;
//       lt1 =   0.000000  ;
//       lt2 =   0.000000  ; 
//       lt3 =   0.000000  ; 
//       tt1 =   0.000000  ;
//       tt2 =   0.000000  ;
//       tt3 =   0.000000  ;

//       /*--------------------------------------------------*/
      /// Itt_19
      t1  =   0.23417  ;
      t2  =  -12.0748  ;
      t3  =   -0.0108  ; 
      l1  =    0.0000  ;
      l2  =    0.0000  ;
      l3  =    0.0616  ;
      lt1 =   -0.0237  ;
      lt2 =    1.2987  ; 
      lt3 =    0.0014  ; 
      tt1 =    0.0166  ;
      tt2 =    0.0000  ;
      tt3 =    0.0000  ;




//     /*--------------------------------------------------*/
//     // Sigma T
//     sigt = t1 * exp( -t2 * up) + t3;
//       
//     /*--------------------------------------------------*/
//     // Sigma L
//     sigl = l1 * exp(l2 * up) + l3;
//       
//     /*--------------------------------------------------*/
//     // Sigma TT  
//     sigtt = tt1 * exp( -tt2 * up) * sin(thetacm) * sin(thetacm);
//       
//     /*--------------------------------------------------*/
//     // Sigma LT  
//     siglt = lt1 * exp( -lt2 * up) * sin(thetacm);
//     //siglt = lt1 * exp( -lt2 * (up-lt3)**2) * sin(thetacm);
// 

     /*--------------------------------------------------*/
     // Sigma T
     sigt = t1 * exp(t2 * up) + t3;
       

     /*--------------------------------------------------*/
     // Sigma L
//     sigl = l1* (exp(-l2*(up-l3)) + l2*(up-l3));
//     sigl = l1 + l2 * up;
//     sigl = l1 * exp( l2 * up ) + l3 / up;
     sigl = l1 * (up-l2)**2 + l3;
 


     /*--------------------------------------------------*/
     // Sigma LT  
     siglt = (lt1 * exp( lt2 * up ) +lt3 / up) * sin(thetacm);
     //siglt = lt1 * exp( -lt2 * (up-lt3)**2) * sin(thetacm);


// 	if (siglt < 0.0) {
// 		
// 		cout << "Negative Sigma: sig check: " << sigl << " " << sigt << " " << siglt << " " << sigtt << " " << sigtt << " " << thetacm << " "<< up << endl;
// //		cout << "Negative Sigma: sig check: " << siglt << " " << lt1 << " " << lt2 << " " << lt3 << endl;
// //		cout << "Negative Sigma: thetacm check: " << up << lt1 * exp( lt2 * up ) +lt3 / up << " " << thetacm << " " << sin(thetacm) << endl;
// 
// 		exit(0);
// 
// 	}
 

     /*--------------------------------------------------*/
     // Sigma TT  
     sigtt = tt1 * sin(thetacm) * sin(thetacm);

	wfactor=((2.47**2-mp**2)**2)/((w_gev**2-mp**2)**2);

	sigl = sigl*wfactor;
	sigt = sigt*wfactor;
	siglt= siglt*wfactor;
	sigtt= sigtt*wfactor;
	
	sig = sigt +eps*sigl +eps*cos(2.*phicm)*sigtt +sqrt(2.*eps*(1.+eps))*cos(phicm)*siglt;
	
	sig = sig/2./pi/1.e06;

	sig_gmh_mod = sig;
	

	if (sig_gmh_mod < 0.0) {
		
//		cout << "Negative Sigma: sig check: " << sigl << " " << sigt << " " << siglt << " " << sigtt << " " << sigtt << " " << thetacm << " "<< up << endl;
		cout << "Negative Sigma: sig check: " << siglt << " " << lt1 << " " << lt2 << " " << lt3 << endl;
//		cout << "Negative Sigma: thetacm check: " << up << lt1 * exp( lt2 * up ) +lt3 / up << " " << thetacm << " " << sin(thetacm) << endl;

//		exit(0);

	}

	return sig_gmh_mod;

}


/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// u dependent model for Q2=2.45
///

Float_t sig_gmh_mod_u_245(Float_t thetacm, Float_t phicm, Float_t u_gev, Float_t q2_gev, Float_t w_gev, Float_t eps) {

	Float_t sig_gmh_mod;

//	sig_gmh_mod = 100.0;

	Float_t sig, wfactor;
	Float_t sigt,sigl,siglt,sigtt;
//	Float_t thetacm,phicm,t_gev,tprime_gev,q2_gev,w_gev,eps,tp;
	Float_t lambda0_sq,lambdapi_sq,alphapi_t,alphapi_0;
	Float_t a,b,c,d,e,f,g,h,fpi,fpi235;
	Float_t m_pi0sq, tp, up;
	Float_t mp = m_p*1000;
	
    Float_t l1,l2,l3;
    Float_t t1,t2,t3;
    Float_t lt1,lt2,lt3;
    Float_t tt1,tt2,tt3;  


//	cout << endl << "Theta: "<< thetacm << "   phi: " << phicm << "   t_gev: " << t_gev << "   t_prime: " << tprime_gev << "   q2_gev: " << q2_gev  <<  "   w_gev: "<< w_gev  << "  eps: " << eps << endl << endl;   



    m_pi0sq= Mpi02/1.e6;

    up = fabs(u_gev);      // just to make sure it's positive

//     /*--------------------------------------------------*/
//     /// Itt_18
//     t1  =     0.1           ; 
//     t2  =    -5.0           ;
//     t3  =     0.0002        ; 
//     l1  =     0.05          ;
//     l2  =     0.6           ;
//     l3  =     0.02          ;
//     lt1 =     0.000000      ;
//     lt2 =     0.000000      ;
//     lt3 =     0.000000      ;
//     tt1 =     0.000000      ;
//     tt2 =     0.000000      ;
//     tt3 =     0.000000      ;
 

    /*--------------------------------------------------*/
    /// Itt_19
    t1  =       0.23417        ; 
    t2  =      -12.0748        ;
    t3  =       -0.0108        ; 
    l1  =        0.0000        ;
    l2  =        0.0000        ;
    l3  =        0.0616        ;
    lt1 =       -0.0237        ;
    lt2 =        1.2987        ;
    lt3 =        0.0014        ;
    tt1 =        0.0166        ;
    tt2 =        0.0000        ;
    tt3 =        0.0000        ;

             
//     /*--------------------------------------------------*/
//     // Sigma T
//     sigt = t1 * exp( -t2 * up) + t3;
//       
//     /*--------------------------------------------------*/
//     // Sigma L
//     sigl = l1 * exp(l2 * up) + l3;
//       
//     /*--------------------------------------------------*/
//     // Sigma TT  
//     sigtt = tt1 * exp( -tt2 * up) * sin(thetacm) * sin(thetacm);
//       
//     /*--------------------------------------------------*/
//     // Sigma LT  
//     siglt = lt1 * exp( -lt2 * up) * sin(thetacm);
//     //siglt = lt1 * exp( -lt2 * (up-lt3)**2) * sin(thetacm);

     /*--------------------------------------------------*/
     // Sigma T
     sigt = t1 * exp(t2 * up) + t3;

     /*--------------------------------------------------*/
     // Sigma L
//     sigl = l1* (exp(-l2*(up-l3)) + l2*(up-l3));
//     sigl = l1 + l2 * up;
//     sigl = l1 * exp( l2 * up ) + l3 / up;
     sigl = l1 * (up-l2)**2 + l3;




     /*--------------------------------------------------*/
     // Sigma LT  
     siglt = (lt1 * exp( lt2 * up ) + lt3 / up) * sin(thetacm);
     //siglt = lt1 * exp( -lt2 * (up-lt3)**2) * sin(thetacm);
      
     /*--------------------------------------------------*/
     // Sigma TT  
     sigtt = tt1 * sin(thetacm) * sin(thetacm);





	wfactor=((2.47**2-mp**2)**2)/((w_gev**2-mp**2)**2);

	sigl = sigl*wfactor;
	sigt = sigt*wfactor;
	siglt= siglt*wfactor;
	sigtt= sigtt*wfactor;
	
	sig = sigt +eps*sigl +eps*cos(2.*phicm)*sigtt +sqrt(2.*eps*(1.+eps))*cos(phicm)*siglt;
	
	sig = sig/2./pi/1.e06;

	sig_gmh_mod = sig;
	
	return sig_gmh_mod;

}


















/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// WL Newly modified model
//

Float_t Cal_New_Mod(Float_t thetacm, Float_t phicm, Float_t t_gev, Float_t tprime_gev, Float_t q2_gev, Float_t w_gev, Float_t eps) {

	Float_t sig_gmh_mod;

//	sig_gmh_mod = 100.0;

	Float_t sig, wfactor;
	Float_t sigt,sigl,siglt,sigtt;
//	Float_t thetacm,phicm,t_gev,tprime_gev,q2_gev,w_gev,eps,tp;
	Float_t lambda0_sq,lambdapi_sq,alphapi_t,alphapi_0;
	Float_t a,b,c,d,e,f,g,h,fpi,fpi235;
	Float_t m_pi0sq, tp;
	Float_t mp = m_p*1000;

//	cout << endl << "Theta: "<< thetacm << "   phi: " << phicm << "   t_gev: " << t_gev << "   t_prime: " << tprime_gev << "   q2_gev: " << q2_gev  <<  "   w_gev: "<< w_gev  << "  eps: " << eps << endl << endl;   



    m_pi0sq= Mpi02/1.e6;



	///*--------------------------------------------------*/
	/// Important !! 
	/// Must be fabs instead of abs!!!!

	tp = fabs(tprime_gev);

	lambda0_sq= 0.462;



//	cout << m_pi0sq << endl;

	if (t_gev>m_pi0sq) { 
	   alphapi_t=0.7*(t_gev-m_pi0sq);
	} else {
	   alphapi_t=tanh(0.7*(t_gev-m_pi0sq));
	}
	
	alphapi_0=0.7*(0.-m_pi0sq);
	lambdapi_sq=lambda0_sq*((1.+alphapi_t)/(1.+alphapi_0))**2;
	fpi    = 1./(1.+q2_gev/lambdapi_sq);
	fpi235 = 1./(1.+2.35/lambdapi_sq);
	
//	cout << fpi <<  "  " << fpi235 << endl;

	a = 0.16675;
//	b = 0.89524;

	b = 0.89524;
	c = 3.5991;
	d = 6.4562;
	e = -9.6199;
	f = 5.8319;

	
	
	/*--------------------------------------------------*/
	/// Iteration: 1  
	// sigl = a*exp(-b*tp)+c*exp(-d*(tp**0.5))+e*exp(-f*(tp**0.25));
	// sigl = sigl *(fpi/fpi235)**2;

	


	/*--------------------------------------------------*/
	/// Iteration: 2  
	// sigl = a*exp(-b*tp) + c*exp(-d*(tp**0.5));
	// sigl = sigl * (1 + 3 * q2_gev)/5;



//	tp =4;
//	tp =5;

	/*--------------------------------------------------*/
	/// Iteration: 3  
	sigl = a*exp(0.08*tp) + c*exp(-d*(tp**0.5));
	sigl = sigl * (1 + 6 * q2_gev);


	



	a = -0.12286;
//	b = 0.56383;
	b = 1;
	c = 1.4673;
	d = 2.1988;
	e = 0.65170;
	f = 18.501;

	/*--------------------------------------------------*/
	/// Iteration: 1
	// sigt = a*exp(-b*tp)+c*exp(-d*(tp**0.5))+e*exp(-f*(tp**2));
	// sigt = sigt *(fpi/fpi235)**2;


	/*--------------------------------------------------*/
	/// Iteration: 2
	// sigt = a*exp(-b*tp) + c*exp(-d*(tp**0.5));
	// sigt = sigt * (1 + 4 * q2_gev)/5;

//	sigt = a*exp(-b*tp) + c*exp(-d*(tp**0.5));
	sigt = 0.005*exp(0.1*tp);

	sigt = sigt * (1 + 8 * q2_gev);


// 
// 	cout << tp << "Sigma L: " << sigl << endl;;
// 	cout << "Sigma T: "<< sigt << endl;;
// 
//	exit(0);



// 	/*--------------------------------------------------*/
// 	/// 
// 	a = 0.46397;
// 	b = 4.9179;
// 	c = 3.3023;
// 	d = 3.1741;
// 	siglt = a*exp(-b*tp)+c*exp(-d*(tp**0.5));
// 	siglt = siglt*(fpi/fpi235)**2*sin(thetacm);
// 	siglt = -siglt/sqrt(2.);
// 	
// //	cout << "Sigma LT: "<< siglt << endl;;
// 
// 
// 	/*--------------------------------------------------*/
// 	/// 
// 
// 	a = -0.26497;
// 	b = 2.7655;
// 	c = 2.8034;
// 	d = 3.8586;
// 	sigtt = a*exp(-b*tp)+c*exp(-d*(tp**0.5));
// 	sigtt= sigtt*(fpi/fpi235)**2*(sin(thetacm))**2;
// 	
// //	cout << "Sigma TT: " << sigtt << endl;;


	siglt = 0;
	sigtt = 0;
	






	/*--------------------------------------------------*/
	/*--------------------------------------------------*/


	wfactor=((2.47**2-mp**2)**2)/((w_gev**2-mp**2)**2);
	
	sigl = sigl*wfactor;
	sigt = sigt*wfactor;
	siglt= siglt*wfactor;
	sigtt= sigtt*wfactor;
	
	sig = sigt +eps*sigl +eps*cos(2.*phicm)*sigtt +sqrt(2.*eps*(1.+eps))*cos(phicm)*siglt;
	
	sig = sig/2./pi/1.e06;

	sig_gmh_mod = sig;
	
	
	
	return sig_gmh_mod;

}


















































/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// GMH model original 
//
/// Do Not Change!

Float_t sig_gmh_mod(Float_t thetacm, Float_t phicm, Float_t t_gev, Float_t tprime_gev, Float_t q2_gev, Float_t w_gev, Float_t eps) {

	Float_t sig_gmh_mod;

//	sig_gmh_mod = 100.0;

	Float_t sig, wfactor;
	Float_t sigt,sigl,siglt,sigtt;
//	Float_t thetacm,phicm,t_gev,tprime_gev,q2_gev,w_gev,eps,tp;
	Float_t lambda0_sq,lambdapi_sq,alphapi_t,alphapi_0;
	Float_t a,b,c,d,e,f,fpi,fpi235;
	Float_t m_pi0sq, tp;
	Float_t mp = m_p*1000;

//	cout << endl << "Theta: "<< thetacm << "   phi: " << phicm << "   t_gev: " << t_gev << "   t_prime: " << tprime_gev << "   q2_gev: " << q2_gev  <<  "   w_gev: "<< w_gev  << "  eps: " << eps << endl << endl;   



    m_pi0sq= Mpi02/1.e6;

	tp = fabs(tprime_gev);

	lambda0_sq= 0.462;



//	cout << m_pi0sq << endl;

	if (t_gev>m_pi0sq) { 
	   alphapi_t=0.7*(t_gev-m_pi0sq);
	} else {
	   alphapi_t=tanh(0.7*(t_gev-m_pi0sq));
	}
	
	alphapi_0=0.7*(0.-m_pi0sq);
	lambdapi_sq=lambda0_sq*((1.+alphapi_t)/(1.+alphapi_0))**2;
	fpi    = 1./(1.+q2_gev/lambdapi_sq);
	fpi235 = 1./(1.+2.35/lambdapi_sq);
	
	
	a = 0.16675;
	b = 0.89524;
	c = 3.5991;
	d = 6.4562;
	e = -9.6199;
	f = 5.8319;
	sigl = a*exp(-b*tp)+c*exp(-d*(tp**0.5))+e*exp(-f*(tp**0.25));
	sigl = sigl *(fpi/fpi235)**2;
	

//	cout << "Sigma L: "<< fpi <<  "  " << fpi235 << "   "<< tp << "  " << sigl << endl;;


	a = -0.12286;
	b = 0.56383;
	c = 1.4673;
	d = 2.1988;
	e = 0.65170;
	f = 18.501;
	sigt = a*exp(-b*tp)+c*exp(-d*(tp**0.5))+e*exp(-f*(tp**2));
	sigt = sigt *(fpi/fpi235)**2;

//	cout << "Sigma T: "<< sigt << endl;;

	
	a = 0.46397;
	b = 4.9179;
	c = 3.3023;
	d = 3.1741;
	siglt = a*exp(-b*tp)+c*exp(-d*(tp**0.5));
	siglt = siglt*(fpi/fpi235)**2*sin(thetacm);
	siglt = -siglt/sqrt(2.);
	
//	cout << "Sigma LT: "<< siglt << endl;;



	a = -0.26497;
	b = 2.7655;
	c = 2.8034;
	d = 3.8586;
	sigtt = a*exp(-b*tp)+c*exp(-d*(tp**0.5));
	sigtt= sigtt*(fpi/fpi235)**2*(sin(thetacm))**2;
	
//	cout << "Sigma TT: " << sigtt << endl;;

	wfactor=((2.47**2-mp**2)**2)/((w_gev**2-mp**2)**2);
	
//	cout << "wfactor: " << mp << "   " << wfactor << endl;;


	sigl = sigl*wfactor;
	sigt = sigt*wfactor;
	siglt= siglt*wfactor;
	sigtt= sigtt*wfactor;
	
	sig = sigt +eps*sigl +eps*cos(2.*phicm)*sigtt +sqrt(2.*eps*(1.+eps))*cos(phicm)*siglt;
	
	sig = sig/2./pi/1.e06;

	sig_gmh_mod = sig;
	
	
	
	return sig_gmh_mod;

}

