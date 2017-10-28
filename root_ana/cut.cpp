#include <iostream>
#include <stdlib.h> 

#include "TString.h"


// /*--------------------------------------------------*/
// Define Cuts
// Define HMS Cuts
//
/*--------------------------------------------------*/
// Kumac code HMS cut from Paw script
//    cuts $11 abs(hsdelta)<8.0
// *   cuts $12 abs(hsxptar)<0.060  | cut on Jochen's thesis p.53
//    cuts $12 abs(hsxptar)<0.080   | extra sieve slit fitting should improve things
//    cuts $13 abs(hsyptar)<0.035
// 
//    cuts $10 $11.and.$12.and.$13

using namespace std;


TString HMS_cut() {

	TString hsdelta_cut = "abs(hsdelta) < 8.0";

    /*--------------------------------------------------*/
    /// *   cuts $12 abs(hsxptar)<0.060  | cut on Jochen's thesis p.53
    /// cuts $12 abs(hsxptar)<0.080   | extra sieve slit fitting should improve things
    /*--------------------------------------------------*/

	TString hsxptar_cut = "abs(hsxptar) < 0.080";
	TString hsyptar_cut = "abs(hsyptar) < 0.035";

	TString hms_cut = hsdelta_cut + " && " + hsxptar_cut + " && " + hsyptar_cut;

 	return hms_cut;

}

/*--------------------------------------------------*/
// /*--------------------------------------------------*/
// Define SOS Cuts from PAW script
//    cuts $21 abs(ssdelta)<15                       
//    cuts $22 abs(ssxfp)<20	
//    cuts $23 abs(ssxfp+ssxpfp*313.01+5.64)<50.
//    cuts $24 abs(ssyfp+ssypfp*313.01)<30.
//    cuts $25 ssytar<1.5          | cut on Jochen's thesis p.80, for th_SOS>42deg.
//    cuts $26 abs(ssxptar)<0.04
//    cuts $27 abs(ssyptar)<0.065                 

//    Jochen acceptance cut (see p.17 of Blok paper)
//    t1=-0.125+0.00425*ssdelta+0.064*ssytar-0.0017*ssdelta*ssytar
//    t2= 0.125-0.00425*ssdelta+0.064*ssytar-0.0017*ssdelta*ssytar 
//    cuts $28 [t1]<ssyptar<[t2]
//
//    cuts $20 $21.and.$22.and.$23.and.$24.and.$25.and.$26.and.$27.and.$28

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



/*--------------------------------------------------*/
// /*--------------------------------------------------*/
// Define PID Cuts from PAW Script
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

//	TString hsbeta_cut    = "hsbeta > 0.1 && hsbeta < 1.5"; 

	TString hcer_npe_cut  = "hcer_npe < 2.";

//	TString hcer_npe_cut  = "hcer_npe < 1.5";

	TString haero_cut     = "(haero_po + haero_ne) > 3.0";
	TString scer_cut      = "scer_npe > 0.5";
	TString ssshtrk_cut   = "ssshtrk > 0.70";

	TString pid_cut = hsbeta_cut + " && " + haero_cut + " && " + hcer_npe_cut + " && " + scer_cut + " && " + ssshtrk_cut;

	return pid_cut;

}





/*--------------------------------------------------*/
// /*--------------------------------------------------*/
// Kumac diamand cuts from PAW Script
//
//    cuts $51 q2>(-2.916*(w-2.4733)+1.5846)
//    cuts $52 q2>(-4.452*(w-1.94)+3.4199) 
//    cuts $53 q2<(-2.677*(w-2.3308)+2.231)
//    cuts $54 q2<(-4.3748*(w-2.41)+1.78)


TString Diamond_cut() {

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

	TString cut_tmp; 

	cut_tmp.Form("(missmass-%f) > 0.875  && (missmass-%f) < 0.98", m_m_offset, m_m_offset);

	TString missmass_cut = cut_tmp; 

	return missmass_cut;

}


TString Set_t_limit(Double_t tmin, Double_t tmax) {

	TString cut_tmp; 
	cut_tmp.Form("t > %f && t < %f", tmin, tmax);

	TString t_limit = cut_tmp; 

	return t_limit;
}




TString Cointime_primary_cut(Double_t center) {

	TString cut_tmp;
	// cut_tmp.Form("-1 < (cointime - %f) && (cointime - %f) < 1", center, center);

	if(center > 0) {

		cut_tmp.Form("-1 < (cointime - %f) && (cointime - %f) < 1", center, center);

	} else {

		cut_tmp.Form("-1 < (cointime + %f) && (cointime + %f) < 1", -center, -center);

	}




	return cut_tmp;
}




TString Cointime_random_cut(Double_t center) {


	TString cut_tmp;



	if(center > 0) {

		cut_tmp.Form("1 < (cointime - %f) && (cointime - %f) < 7", center, center);

	} else {

		cut_tmp.Form("1 < (cointime + %f) && (cointime + %f) < 7", -center, -center);

	}

 	TString cut1 = cut_tmp.Data(); 

	return cut1;
}


TString Cointime_all(Double_t center) {


	TString cut_tmp;

	if(center > 0) {

		cut_tmp.Form("-20 < (cointime - %f) && (cointime - %f) < 10", center, center);

	} else {

		cut_tmp.Form("-20 < (cointime + %f) && (cointime + %f) < 10", -center, -center);

	}

 	TString cut1 = cut_tmp.Data(); 

	return cut1;
}







TString Heep_PID_Cut(float Q2_set) {

//    TString haero_cut 	  = "abs(haero_su) < 50";

     TString haero_cut;
     TString haero_cut_1;



///*--------------------------------------------------*/

 	if(Q2_set == float(1.91)) {
  
 		haero_cut.Form("haero_su > %f", -50.0);
 		haero_cut_1.Form("haero_su < %f", 90.0);
 
 		haero_cut = haero_cut + " && "+ haero_cut_1;
 
 	} else {
 
 		haero_cut.Form("haero_su < %f", 29.0);
 
 	}

// 
// 	cout << Q2_set << endl;
// 
///*--------------------------------------------------*/
// 	if(Q2_set == float(1.91)) {
// 
// 		haero_cut.Form("abs(haero_su) < %f", 25.0);
// 
// 	} else if (Q2_set == float(6.53)) {
// 
// 		haero_cut.Form("haero_su > %f", 0.5);
// 		haero_cut_1.Form("haero_su < %f", 10.0);
// 
// 		haero_cut = haero_cut + " && "+ haero_cut_1;
// 
// 	} else {
// 
// 		haero_cut.Form("abs(haero_su) < %f", 2.5);
// 
// 	}







//	cout << haero_cut << endl;

// 	TString hsbeta_cut    = "hsbeta > 0.1 && hsbeta < 1.5"; 
// 	TString hcer_npe_cut  = "hcer_npe < 1.5";
// 	TString haero_cut     = "(haero_po + haero_ne) > 3.0";
// 	TString scer_cut      = "scer_npe > 0.5";
// 	TString ssshtrk_cut   = "ssshtrk > 0.70";
//
//

//	TString pid_cut = hsbeta_cut + " && " + haero_cut + " && " + hcer_npe_cut + " && " + scer_cut + " && " + ssshtrk_cut;

//	return pid_cut;

	// Beta Cut
// 	TString hsbeta_cut    = "hsbeta > 0.1 && hsbeta < 1.5"; 


///*--------------------------------------------------*/

// 	/*--------------------------------------------------*/
// 	///
// 	// Beta Cut
//  	TString hsbeta_cut    = "hsbeta > -0.1 && hsbeta < 1.5"; 
//  	TString hcer_npe_cut  = "hcer_npe < 0.5";
// // 	TString haero_cut  	  = "abs(haero_su) < 30";


	/*--------------------------------------------------*/
	/// Omega PID cuts
// 	TString hsbeta_cut    = "hsbeta > -0.1 && hsbeta < 1.5"; 
 	TString hcer_npe_cut  = "hcer_npe < 2.0";
// 	TString haero_cut  	  = "abs(haero_su) < 200";



	// SOS Showercounter cut
	TString ssshtrk_cut = "ssshtrk > 0.70";           
	TString scer_npe_cut = "scer_npe > 0.50";           

	// TString pid_cut = hsbeta_cut + " && " + hcer_npe_cut + " && " + haero_cut + " && "+ ssshtrk_cut + " && " + scer_npe_cut;

	// TString pid_cut = hsbeta_cut + " && " + hcer_npe_cut + " && " + haero_cut + " && "+ ssshtrk_cut + " && " + scer_npe_cut;
	//TString pid_cut = hsbeta_cut + " && "+ hcer_npe_cut + " && " + haero_cut + " && "+ ssshtrk_cut + " && " + scer_npe_cut;


	TString pid_cut = hcer_npe_cut + " && " + haero_cut + " && "+ ssshtrk_cut + " && " + scer_npe_cut;
//	TString pid_cut = hcer_npe_cut + " && "+ ssshtrk_cut + " && " + scer_npe_cut;

	return pid_cut;

}


TString T_heep_cut() {

	TString cut_tmp; 

	// cut_tmp.Form("(missmass-%f) > 0.875  && (missmass-%f) < 0.98", m_m_offset, m_m_offset);
	cut_tmp = "t < 0.15";

	TString missmass_cut = cut_tmp; 

	return missmass_cut;

}





TString Missingmass_heep_cut() {

//	TString cut_tmp; 

	// cut_tmp.Form("(missmass-%f) > 0.875  && (missmass-%f) < 0.98", m_m_offset, m_m_offset);
//	cut_tmp = "missmass > -0.11  && missmass < 0.13";
//	cut_tmp = "missmass > -0.1 && missmass < 0.6 && Em > -0.02 && Em <0.1";
//	cut_tmp = " Em > -0.02 && Em <0.1";
//	cut_tmp = " Em <0.1 && Em > -0.03";


///*--------------------------------------------------*/
/// Garth's propose cut A:
//	cut_tmp = " missmass > -0.14 && Em < 0.1";


///*--------------------------------------------------*/
/// Garth's propose cut B:
///	cut_tmp = " missmass > -0.14 && missmass < 0.08";


///*--------------------------------------------------*/
/// Garth's propose cut C (A+B) :
//	cut_tmp = " missmass > -0.14 && Em < 0.1 && missmass < 0.08";


//	cut_tmp = " missmass > -0.14";

//	cut_tmp = " missmass > -0.14 && Em < 0.35 && missmass < 0.06";
//	cut_tmp = " missmass > -0.14 && Em < 0.3 && missmass < 0.08";
//	cut_tmp = " missmass > -0.14 && Em < 0.4 && missmass < 0.08";




//	cut_tmp = " missmass > -0.14 && Em < 0.3 && missmass < 0.04";
//	cut_tmp = " missmass > -0.14 && Em < 0.3 && missmass < 0.05";
//	cut_tmp = " missmass > -0.14 && Em < 0.3 && missmass < 0.06";



//	cut_tmp = " missmass > -0.14 && missmass < 0.6";


//	cut_tmp = " missmass > -0.14 && Em < 0.35 && missmass < 0.08";


/*--------------------------------------------------*/
	// Good cut 
//	cut_tmp = "Em < 0.27";

//  Old acceptance cut	
//	cut_tmp = " missmass > -0.085 && Em < 0.27 && missmass < 0.04";


//  Garth acceptance cut	
//	cut_tmp = " missmass > -0.06 && Em < 0.27 && missmass < 0.04";
//	cut_tmp = " missmass > -0.15 && Em < 0.27 && missmass < 0.04";
//
//

	TString E_cut;
	TString missmass_cut;
	TString total_cut;


//	E_cut = "Em < 0.7";


//	missmass_cut = "missmass > -0.07 && missmass < 0.05";
// 	missmass_cut = "missmass > -0.20 && missmass < 0.8";
// 	missmass_cut = "missmass > -0.05 && missmass < 0.40";


//	missmass_cut = "missmass > -0.06 && missmass < 0.05";
//	missmass_cut = "missmass > -0.04 && missmass < 0.05";



//	E_cut = "Em < 0.27";
//	missmass_cut = "missmass > -0.5 && missmass < 0.3";



//	E_cut = "Em < 0.27";
//	missmass_cut = "missmass > -0.05 && missmass < 0.03";



// 	/*--------------------------------------------------*/
 	E_cut = "Em < 0.10";
// 	missmass_cut = "missmass > -0.03 && missmass < 0.02";


// 	missmass_cut = "missmass > -0.035 && missmass < 0.025";
//
//

// 	missmass_cut = "missmass > -0.04 && missmass < 0.03";
// 	missmass_cut = "missmass > -0.045 && missmass < 0.035";
// 	missmass_cut = "missmass > -0.050 && missmass < 0.040";

// 	missmass_cut = "missmass > -0.055 && missmass < 0.045";

//	E_cut = "Em < 0.27";
//	missmass_cut = "missmass > -0.1 && missmass < 0.1";


// 	missmass_cut = "missmass > -0.027 && missmass < 0.013";
 	missmass_cut = "missmass > -0.032 && missmass < 0.018";
// 	missmass_cut = "missmass > -0.037 && missmass < 0.023";

//
//
//
// 	missmass_cut = "missmass > -0.03 && missmass < 0.025";

// 	missmass_cut = "missmass > -0.0365 && missmass < 0.0235";
// 	missmass_cut = "missmass > -0.0415 && missmass < 0.0285";
// 	missmass_cut = "missmass > -0.0465 && missmass < 0.0335";


	total_cut = missmass_cut + " && " +  E_cut;
//	total_cut =   E_cut;


//	total_cut = "missmass > -10";

	return total_cut;

}




/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Singles check cut on W

TString Heep_single_W_cut(float w_down, float w_up) {

	TString cut_tmp; 

	// cut_tmp.Form("(missmass-%f) > 0.875  && (missmass-%f) < 0.98", m_m_offset, m_m_offset);
	// cut_tmp = "missmass > -0.11  && missmass < 0.13";


	cut_tmp.Form("W > %f && W < %f", w_down, w_up);

	TString missmass_cut = cut_tmp; 

	return missmass_cut;

}

/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Singles PID cut

TString Heep_Single_PID_cut() {

	TString scer_cut      = "scer_npe > 0.5";
	TString ssshtrk_cut   = "ssshtrk > 0.70";

	TString pid_cut = scer_cut + " && " + ssshtrk_cut;

	return pid_cut;

} 


/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Singles SOS cut
TString Heep_Single_SOS_cut() {

	TString ssdelta_cut = "abs(ssdelta) < 15.";                  
	TString ssxfp1_cut  = "abs(ssxfp) < 20.";	

	TString ssxfp2_cut  = "abs(ssxfp+ssxpfp*313.01+5.64) < 50."; 
	TString ssyfp_cut   = "abs(ssyfp+ssypfp*313.01) < 30.";
	TString ssytar_cut  = "ssytar < 1.5";       

	TString ssxptar_cut = "abs(ssxptar) < 0.04";
	TString ssyptar_cut = "abs(ssyptar) < 0.065";               

	TString t1 = "ssyptar > (-0.125+0.00425*ssdelta+0.064*ssytar-0.0017*ssdelta*ssytar)";               
	TString t2 = "ssyptar < (0.125-0.00425*ssdelta+0.064*ssytar-0.0017*ssdelta*ssytar)";               

// 	TString sos_cut = ssdelta_cut + " && " + ssxfp1_cut + " && " +  ssxfp2_cut + " && " + 
// 				   ssyfp_cut  + " && " + ssytar_cut + " && " + ssxptar_cut + " && " + 
// 				   ssyptar_cut + " && " + t1 + " && "+ t2;
//
//
//

	TString sos_cut = ssdelta_cut + " && " + ssxfp1_cut + " && " + ssxptar_cut + " && " + 
				   ssyptar_cut + " && " + t1 + " && "+ t2;

   	return sos_cut; 



}


/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Garth's Omega cuts form teleconf report dated 05/April/2005

TString Omega_HMS_Cut_Garth() {

	TString hsdelta_cut = "abs(hsdelta) <= 8.0";

	TString hsxptar_cut = "abs(hsxptar) <= 0.090";
	TString hsyptar_cut = "abs(hsyptar) <= 0.055";

	TString hms_cut = hsdelta_cut + " && " + hsxptar_cut + " && " + hsyptar_cut;

 	return hms_cut;

}

TString Omega_SOS_Cut_Garth() {

	TString ssdelta_cut = "abs(ssdelta) <= 15.";                  

	TString ssxptar_cut = "abs(ssxptar) <= 0.04";
	TString ssyptar_cut = "abs(ssyptar) <= 0.08";               

	TString ssxfp1_cut  = "ssxfp >= -20.";	
	TString ssxfp2_cut  = "ssxfp <=  22.";	

	TString sos_cut = ssdelta_cut + " && " + ssxfp1_cut + " && " +  ssxfp2_cut + " && " + 
				   ssxptar_cut + " && " + ssyptar_cut;

   	return sos_cut; 

}

TString Omega_PID_Cut_Garth() {

 	TString haero_cut  	 = "abs(haero_su) <= 2.5";
 	TString hcer_npe_cut = "hcer_npe <= 2.0";

	// SOS Showercounter cut
	TString ssshtrk_cut  = "ssshtrk >= 0.70";           
	TString scer_npe_cut = "scer_npe >= 0.50";           

	TString pid_cut = hcer_npe_cut + " && " + haero_cut + " && "+ ssshtrk_cut + " && " + scer_npe_cut;

//	TString pid_cut = hcer_npe_cut + " && "+ ssshtrk_cut + " && " + scer_npe_cut;

	return pid_cut;

}


/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Modified cut based on Garth's Omega cuts and stardard pion analysis cut
//  01/March/2016

TString Omega_HMS_Cut() {

	
	TString hsytar_cut  = "abs(hsytar) <= 1.75";

	TString hsdelta_cut = "abs(hsdelta) <= 8.0";

    /*--------------------------------------------------*/
    /// *   cuts $12 abs(hsxptar)<0.060  | cut on Jochen's thesis p.53
    /// cuts $12 abs(hsxptar)<0.080   | extra sieve slit fitting should improve things
    /*--------------------------------------------------*/

	TString hsxptar_cut = "abs(hsxptar) <= 0.080";
	TString hsyptar_cut = "abs(hsyptar) <= 0.035";

//	TString hms_cut = hsdelta_cut + " && " + hsxptar_cut + " && " + hsyptar_cut;

	TString hms_cut = hsytar_cut + " && "+ hsdelta_cut + " && " + hsxptar_cut + " && " + hsyptar_cut;

 	return hms_cut;

}


TString Omega_SOS_Cut() {



	TString ssytar_cut  = "ssytar < 1.5";       

	TString ssdelta_cut = "abs(ssdelta) <= 15.";                  
	TString ssxfp1_cut  = "abs(ssxfp) <= 20.";	

	TString ssxptar_cut = "abs(ssxptar) <= 0.04";
	TString ssyptar_cut = "abs(ssyptar) <= 0.065";               

//	TString ssxfp2_cut  = "abs(ssxfp+ssxpfp*313.01+5.64) < 50."; 
//	TString ssyfp_cut   = "abs(ssyfp+ssypfp*313.01) < 30.";

//	TString t1 = "ssyptar > (-0.125+0.00425*ssdelta+0.064*ssytar-0.0017*ssdelta*ssytar)";               
//	TString t2 = "ssyptar < (0.125-0.00425*ssdelta+0.064*ssytar-0.0017*ssdelta*ssytar)";               



	TString sos_cut = ssdelta_cut + " && " + ssxfp1_cut + " && " + ssytar_cut + " && " 
					+ ssxptar_cut + " && " + ssyptar_cut;

   	return sos_cut; 

}


TString Omega_PID_Cut() {

	// A loose beta Cut
 	TString hsbeta_cut    = "hsbeta > 0.2 && hsbeta < 1.5";
 
 	TString haero_cut  	  = "abs(haero_su) < 4";
 	TString hcer_npe_cut  = "hcer_npe < 0.5";

	// SOS Showercounter cut
	TString ssshtrk_cut = "ssshtrk >= 0.70";           
	TString scer_npe_cut = "scer_npe >= 0.50";           

	TString pid_cut = hsbeta_cut + " && " + hcer_npe_cut + " && " + haero_cut + " && "+ ssshtrk_cut + " && " + scer_npe_cut;

	return pid_cut;

}


TString Omega_U_Cut() {

	// A loose beta Cut 
 	TString u_cut = "(Pm**2-(0.938272-Em)**2) > -100";
 	//TString u_cut = "(Pm**2-(0.938272-Em)**2) > 0";
 
	return u_cut;

}



TString Omega_U_Sim_Cut() {

	// A loose beta Cut 
 	TString u_cut = "minus_u > -100";
 	//TString u_cut = "minus_u > 0";
 
	return u_cut;

}



/*--------------------------------------------------*/
// Missmass cut

TString Omega_Missmass_Cut() {

 	TString mm_cut = "missmass > 0.65";
 
	return mm_cut;

}

