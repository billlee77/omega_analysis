/*--------------------------------------------------*/
/*--------------------------------------------------*/
/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Fitting function defined

Double_t ftotal(Double_t *x, Double_t *par) {

   Double_t xx = x[0];
   Int_t bin = peak->GetXaxis()->FindBin(xx);

//   Double_t br = par[3]*background->GetBinContent(bin);
//   Double_t arg = (xx-par[1])/par[2];
//   Double_t sr = par[0]*TMath::Exp(-0.5*arg*arg);

   Double_t sr = par[4] * peak->GetBinContent(bin);
   Double_t br = par[0]  + par[1] * xx  + par[2]* xx * xx + par[3]* xx * xx *  xx;

   return sr + br;

}


Double_t ftotal_all(Double_t *x, Double_t *par) {
	return ftotal_br(x, par) + ftotal_sg1(x, &par[4]);
}



Double_t ftotal_all_all(Double_t *x, Double_t *par) {
	return ftotal_br(x, par) + ftotal_sg1(x, &par[4]) + ftotal_sg2(x, &par[5]);
}


/*--------------------------------------------------*/

Double_t ftotal_br(Double_t *x, Double_t *par) {

   Double_t xx = x[0];

//   Double_t br = par[3]*background->GetBinContent(bin);
//   Double_t arg = (xx-par[1])/par[2];
//   Double_t sr = par[0]*TMath::Exp(-0.5*arg*arg);

   Double_t br = par[0]  + par[1] * xx  + par[2]* xx * xx + par[3]* xx * xx *  xx;

   return br;

}

/*--------------------------------------------------*/

Double_t ftotal_sg1(Double_t *x, Double_t *par) {

   Double_t xx = x[0];
   Int_t bin = peak1->GetXaxis()->FindBin(xx);

   Double_t sr = par[0] * peak1->GetBinContent(bin);

   return sr;

}

/*--------------------------------------------------*/

Double_t ftotal_sg2(Double_t *x, Double_t *par) {

   Double_t xx = x[0];
   Int_t bin = peak2->GetXaxis()->FindBin(xx);

   Double_t sr = par[0] * peak2->GetBinContent(bin);

   return sr;

}




// /*--------------------------------------------------*/
// /*--------------------------------------------------*/
// 
// Double_t ftotal_1(Double_t *x, Double_t *par) {
// 
//    Double_t xx = x[0];
//    Int_t bin = peak1->GetXaxis()->FindBin(xx);
// 
// //   Double_t br = par[3]*background->GetBinContent(bin);
// //   Double_t arg = (xx-par[1])/par[2];
// //   Double_t sr = par[0]*TMath::Exp(-0.5*arg*arg);
// 
//    Double_t sr = par[4] * peak->GetBinContent(bin);
//    Double_t br = par[0]  + par[1] * xx  + par[2]* xx * xx + par[3]* xx * xx *  xx;
// 
//    return sr + br;
// 
// }
// 
// /*--------------------------------------------------*/
// /*--------------------------------------------------*/
// 
// Double_t ftotal_2(Double_t *x, Double_t *par) {
// 
//    Double_t xx = x[0];
//    Int_t bin = peak1->GetXaxis()->FindBin(xx);
// 
// //   Double_t br = par[3]*background->GetBinContent(bin);
// //   Double_t arg = (xx-par[1])/par[2];
// //   Double_t sr = par[0]*TMath::Exp(-0.5*arg*arg);
// 
//    Double_t sr = par[4] * peak->GetBinContent(bin);
//    Double_t br = par[0]  + par[1] * xx  + par[2]* xx * xx + par[3]* xx * xx *  xx;
// 
//    return sr + br;
// 
// }
//
//
//
//
//
//


/*--------------------------------------------------*/
/*--------------------------------------------------*/

Double_t ftotal_all_2(Double_t *x, Double_t *par) {
	return ftotal_br_2(x, par) + ftotal_sg1(x, &par[3]);
}



/*--------------------------------------------------*/
/*--------------------------------------------------*/

Double_t ftotal_br_2(Double_t *x, Double_t *par) {

   Double_t xx = x[0];

//   Double_t br = par[3]*background->GetBinContent(bin);
//   Double_t arg = (xx-par[1])/par[2];
//   Double_t sr = par[0]*TMath::Exp(-0.5*arg*arg);

   Double_t br = par[0]  + par[1] * xx  + par[2]* xx * xx;

   return br;

}



/*--------------------------------------------------*/
/*--------------------------------------------------*/


Double_t ftotal_const_bg_fun(Double_t *x, Double_t *par) {

//   return  ftotal_br(x, par) + ftotal_sg1(x, &par[4]);

//  return  ftotal_br(x, par) + ftotal_sg1(x, &par[4]) + ftotal_bg_scale(x, &par[5]);



//   return  ftotal_br(x, par) * 0.01;
   return  ftotal_br(x, par) * ftotal_bg_scale(x, &par[4]) + ftotal_sg1(x, &par[5]);

}



Double_t ftotal_bg_scale(Double_t *x, Double_t *par) {

   Double_t xx = x[0];

//   Double_t br = par[3]*background->GetBinContent(bin);
//   Double_t arg = (xx-par[1])/par[2];
//   Double_t sr = par[0]*TMath::Exp(-0.5*arg*arg);

   Double_t br = par[0];

   return br;

}


/*--------------------------------------------------*/

Double_t ftotal_const_bg_bg_only(Double_t *x, Double_t *par) {

   return ftotal_br(x, par) * ftotal_bg_scale(x, &par[4]);

}


/*--------------------------------------------------*/

Double_t skewed_gaussion_fun(Double_t *x, Double_t *par){

		Double_t xx = x[0];

//		Double_t sgf = par[0] * exp(-0.5 * pow((x-par[1])/(par[2]+(x-par[1])*par[3]*(x-par[1])),2) );
		// p0*exp(-0.5*((x-p1)/p2)^2)


/*--------------------------------------------------*/
//// Suggested skewed gaussian function
//		Double_t gausss = par[0] * exp(-0.5 * pow((xx-par[1])/(par[2]+(xx-par[1])*par[3]*(xx-par[1])), 2) );

 

		Double_t gauss_pt = par[0] * exp(-0.5 * pow((xx-par[1])/par[2], 2) );

		Double_t erf_pt = TMath::Erf(xx/par[3]-par[4]);

		return gauss_pt * erf_pt;		

}



/*--------------------------------------------------*/

Double_t ftotal_all_skewed_gauss(Double_t *x, Double_t *par) {
//	return skewed_gaussion_fun(x, par) + + ftotal_sg1(x, &par[4]);
	return skewed_gaussion_fun(x, par) + ftotal_sg1(x, &par[4]);
//	return skewed_gaussion_fun(x, par);

}



/*--------------------------------------------------*/
/// Eta prime peak function
Double_t ftotal_eta_prime_peak(Double_t *x, Double_t *par) {

  return (0.5*par[0]*par[1]/TMath::Pi()) / TMath::Max( 1.e-10,(x[0]-par[2])*(x[0]-par[2]) + .25*par[1]*par[1]);

}



Double_t ftotal_bg_sg_eta(Double_t *x, Double_t *par) {
	return ftotal_br(x, par) + ftotal_sg1(x, &par[4]) + ftotal_eta_prime_peak(x, &par[5]);
}







/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Simulation background



/*--------------------------------------------------*/

Double_t fun_omega(Double_t *x, Double_t *par) {

   Double_t xx = x[0];
   Int_t bin = peak_omega->GetXaxis()->FindBin(xx);

   Double_t sr = par[0] * peak_omega->GetBinContent(bin);

   return sr;

}


/*--------------------------------------------------*/

Double_t fun_rho(Double_t *x, Double_t *par) {

   Double_t xx = x[0];
   Int_t bin = peak_rho->GetXaxis()->FindBin(xx);

   Double_t sr = par[0] * peak_rho->GetBinContent(bin);

   return sr;

}


/*--------------------------------------------------*/

Double_t fun_xphsp(Double_t *x, Double_t *par) {

   Double_t xx = x[0];
   Int_t bin = peak_xphsp->GetXaxis()->FindBin(xx);

   Double_t sr = par[0] * peak_xphsp->GetBinContent(bin);

   return sr;

}

/*--------------------------------------------------*/

Double_t fun_eta(Double_t *x, Double_t *par) {

   Double_t xx = x[0];
   Int_t bin = peak_eta->GetXaxis()->FindBin(xx);

   Double_t sr = par[0] * peak_eta->GetBinContent(bin);

   return sr;

}

/*--------------------------------------------------*/

Double_t fun_etap(Double_t *x, Double_t *par) {

   Double_t xx = x[0];
   Int_t bin = peak_etap->GetXaxis()->FindBin(xx);

   Double_t sr = par[0] * peak_etap->GetBinContent(bin);

   return sr;

}

/*--------------------------------------------------*/

Double_t fun_sim_total(Double_t *x, Double_t *par) {
//	return fun_omega(x, par) + fun_rho(x, &par[1]) + fun_xphsp(x, &par[2]);
	return fun_omega(x, par) + fun_rho(x, &par[1]) + fun_xphsp(x, &par[2]) + fun_eta(x, &par[3]) + fun_etap(x, &par[4]);
}


