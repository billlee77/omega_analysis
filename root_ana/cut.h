#ifndef __CUT_H_INCLUDED__
#define __CUT_H_INCLUDED__

#include "TString.h"

/*--------------------------------------------------*/
/// Cut definition
//
extern TString HMS_cut();

extern TString SOS_cut();

extern TString PID_cut();
 
extern TString Diamond_cut();

extern TString Missingmass_cut(Double_t m_m_offset);

extern TString Set_t_limit(Double_t tmin, Double_t tmax);
 
extern TString Cointime_primary_cut(Double_t center); 

extern TString Cointime_random_cut(Double_t center);

extern TString Cointime_all(Double_t center);

extern TString Q2_limit;

extern TString Heep_PID_Cut(float);

extern TString Missingmass_heep_cut();

extern TString T_heep_cut();

extern TString Heep_single_W_cut(float , float);

extern TString Heep_Single_PID_cut();

extern TString Heep_Single_SOS_cut();


/*--------------------------------------------------*/
extern TString Omega_HMS_Cut_Garth();
extern TString Omega_SOS_Cut_Garth();
extern TString Omega_PID_Cut_Garth();


/*--------------------------------------------------*/
extern TString Omega_HMS_Cut();
extern TString Omega_SOS_Cut();
extern TString Omega_PID_Cut();

extern TString Omega_U_Cut();
extern TString Omega_U_Sim_Cut();
extern TString Omega_Missmass_Cut();

#endif


