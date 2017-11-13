{


	TFile *f1 = TFile::Open("yields.omega_160_32_+0970.target.root");
	TFile *f2 = TFile::Open("yields.omega_160_32_+3000.target.root");

  	TH1F* Q2_160_l_l_angle_tmp = (TH1F*) f1->Get("t_dia");
  	TH1F* Q2_160_l_l_angle = (TH1F*) Q2_160_l_l_angle_tmp->Clone();
// 
  	TH1F* Q2_160_l_h_angle_tmp = (TH1F*) f2->Get("t_dia");
  	TH1F* Q2_160_l_h_angle = (TH1F*) Q2_160_l_h_angle_tmp->Clone();

  	TCanvas * d1 = new TCanvas();

 	Q2_160_l_l_angle->SetLineColor(4);
 	Q2_160_l_h_angle->SetLineColor(2);

//  	TH1F* Q2_160_l_l_angle_sum = (TH1F*)Q2_160_l_l_angle->Clone();
// 
// 	Q2_160_l_l_angle_sum->Add(Q2_160_l_h_angle, 1);
// 
// 	Q2_160_l_l_angle_sum->Draw();
 
 	Q2_160_l_l_angle->Draw();
  	Q2_160_l_h_angle->Draw("same");


  	TH1F* Q2_160_l_l_angle_ph_tmp = (TH1F*) f1->Get("phi_pq_dia");
  	TH1F* Q2_160_l_l_angle_ph = (TH1F*) Q2_160_l_l_angle_ph_tmp->Clone();
// 
  	TH1F* Q2_160_l_h_angle_ph_tmp = (TH1F*) f2->Get("phi_pq_dia");
  	TH1F* Q2_160_l_h_angle_ph = (TH1F*) Q2_160_l_h_angle_ph_tmp->Clone();

  	TCanvas * t1 = new TCanvas("t1", "t1");

 	Q2_160_l_l_angle_ph->SetLineColor(4);
 	Q2_160_l_h_angle_ph->SetLineColor(2);

  	Q2_160_l_h_angle_ph->Draw();
 	Q2_160_l_l_angle_ph->Draw("same");


 	TCanvas* c1 = (TCanvas*)f1->Get("t_phi_check");
  	TGraphPolar* Q2_160_l_l_gp = (TGraphPolar*) c1->GetListOfPrimitives()->FindObject("Graph");

	TCanvas* c2 = (TCanvas*)f2->Get("t_phi_check");
  	TGraphPolar* Q2_160_l_h_gp = (TGraphPolar*) c2->GetListOfPrimitives()->FindObject("Graph");


  	TCanvas * p1 = new TCanvas("p1", "p1", 600, 600);

	Q2_160_l_l_gp->SetMarkerColor(4);
	Q2_160_l_h_gp->SetMarkerColor(2);

	Q2_160_l_l_gp->Draw("AOP");

	p1->Update();

	Q2_160_l_l_gp->GetPolargram()->SetNdivPolar(4);
 	Q2_160_l_l_gp->GetPolargram()->SetNdivRadial(4);
 	Q2_160_l_l_gp->GetPolargram()->SetRangeRadial(0, 0.4); 
 
	Q2_160_l_h_gp->GetPolargram()->SetRangeRadial(0, 0.4); 

	Q2_160_l_h_gp->Draw("OPsame");
 
	p1->Update();

	


// 	TMultiGraph * mg1_tmp = (TMultiGraph*) c1->GetListOfPrimitives()->FindObject("diamond_gr");
// 	TMultiGraph * mg1 = (TMultiGraph*) mg1_tmp->Clone();
 
//  	TGraphPolar*  = (TGraphPolar*) f2->Get("t_dia");

 
// 	/*--------------------------------------------------*/
// 	/*--------------------------------------------------*/
// 
// 
	TFile *f1 = TFile::Open("yields.omega_160_59_-2730.target.root");
	TFile *f2 = TFile::Open("yields.omega_160_59_+0000.target.root");
	TFile *f3 = TFile::Open("yields.omega_160_59_+3000.target.root");
 

  	TH2F* Q2_160_h_l_angle_tmp = (TH2F*) f1->Get("t_dia");
  	TH2F* Q2_160_h_l_angle = (TH2F*) Q2_160_h_l_angle_tmp->Clone();
 
  	TH2F* Q2_160_h_0_angle_tmp = (TH2F*) f2->Get("t_dia");
  	TH2F* Q2_160_h_0_angle = (TH2F*) Q2_160_h_0_angle_tmp->Clone();


  	TH2F* Q2_160_h_h_angle_tmp = (TH2F*) f3->Get("t_dia");
  	TH2F* Q2_160_h_h_angle = (TH2F*) Q2_160_h_h_angle_tmp->Clone();


  	TCanvas * d2 = new TCanvas();

 	Q2_160_h_0_angle->SetLineColor(4);
  	Q2_160_h_h_angle->SetLineColor(2);

  	Q2_160_h_0_angle->Draw();

  	Q2_160_h_l_angle->Draw("same");
   	Q2_160_h_h_angle->Draw("same");

  	TH1F* Q2_160_h_l_angle_ph_tmp = (TH1F*) f1->Get("phi_pq_dia");
  	TH1F* Q2_160_h_l_angle_ph = (TH1F*) Q2_160_h_l_angle_ph_tmp->Clone();

  	TH1F* Q2_160_h_0_angle_ph_tmp = (TH1F*) f2->Get("phi_pq_dia");
  	TH1F* Q2_160_h_0_angle_ph = (TH1F*) Q2_160_h_0_angle_ph_tmp->Clone();
 
  	TH1F* Q2_160_h_h_angle_ph_tmp = (TH1F*) f3->Get("phi_pq_dia");
  	TH1F* Q2_160_h_h_angle_ph = (TH1F*) Q2_160_h_h_angle_ph_tmp->Clone();

  	TCanvas * t2 = new TCanvas("t2", "t2");
 
 	Q2_160_h_0_angle_ph->SetLineColor(4);
	Q2_160_h_h_angle_ph->SetLineColor(2);

 	Q2_160_h_h_angle_ph->Draw();
  	Q2_160_h_l_angle_ph->Draw("same");
 	Q2_160_h_0_angle_ph->Draw("same");

 	TCanvas* c1 = (TCanvas*)f1->Get("t_phi_check");
  	TGraphPolar* Q2_160_h_l_gp = (TGraphPolar*) c1->GetListOfPrimitives()->FindObject("Graph");

	TCanvas* c2 = (TCanvas*)f2->Get("t_phi_check");
  	TGraphPolar* Q2_160_h_0_gp = (TGraphPolar*) c2->GetListOfPrimitives()->FindObject("Graph");

	TCanvas* c3 = (TCanvas*)f3->Get("t_phi_check");
  	TGraphPolar* Q2_160_h_h_gp = (TGraphPolar*) c3->GetListOfPrimitives()->FindObject("Graph");

  	TCanvas * p2 = new TCanvas("p2", "p2", 600, 600);

	Q2_160_h_0_gp->SetMarkerColor(4);
	Q2_160_h_h_gp->SetMarkerColor(2);

	Q2_160_h_0_gp->Draw("AOP");

	p2->Update();

	Q2_160_h_0_gp->GetPolargram()->SetNdivPolar(4);
 	Q2_160_h_0_gp->GetPolargram()->SetNdivRadial(4);
 	Q2_160_h_0_gp->GetPolargram()->SetRangeRadial(0, 0.4); 

	Q2_160_h_h_gp->GetPolargram()->SetRangeRadial(0, 0.4); 
	Q2_160_h_l_gp->GetPolargram()->SetRangeRadial(0, 0.4); 

	Q2_160_h_h_gp->Draw("OPsame");
	Q2_160_h_l_gp->Draw("OPsame");

	p2->Update();


/*--------------------------------------------------*/
/*--------------------------------------------------*/
/*--------------------------------------------------*/

// 
// 	TCanvas * d2 = new TCanvas();
// 
// 
// //	TFile* file_out = new TFile("", );
// 
	TFile *f1 = TFile::Open("yields.omega_245_27_+1350.target.root");
	TFile *f2 = TFile::Open("yields.omega_245_27_+3000.target.root");


  	TH2F* Q2_245_l_l_angle_tmp = (TH2F*) f1->Get("t_dia");
  	TH2F* Q2_245_l_l_angle = (TH2F*) Q2_245_l_l_angle_tmp->Clone();
// 
  	TH2F* Q2_245_l_h_angle_tmp = (TH2F*) f2->Get("t_dia");
  	TH2F* Q2_245_l_h_angle = (TH2F*) Q2_245_l_h_angle_tmp->Clone();

  	TCanvas * d3 = new TCanvas();

 	Q2_245_l_l_angle->SetLineColor(4);
 	Q2_245_l_h_angle->SetLineColor(2);

 	Q2_245_l_l_angle->Draw();
  	Q2_245_l_h_angle->Draw("same");

  	TH1F* Q2_245_l_l_angle_ph_tmp = (TH1F*) f1->Get("phi_pq_dia");
  	TH1F* Q2_245_l_l_angle_ph = (TH1F*) Q2_245_l_l_angle_ph_tmp->Clone();
// 
  	TH1F* Q2_245_l_h_angle_ph_tmp = (TH1F*) f2->Get("phi_pq_dia");
  	TH1F* Q2_245_l_h_angle_ph = (TH1F*) Q2_245_l_h_angle_ph_tmp->Clone();

  	TCanvas * t3 = new TCanvas("t3", "t3");

 	Q2_245_l_l_angle_ph->SetLineColor(4);
 	Q2_245_l_h_angle_ph->SetLineColor(2);
 
 	Q2_245_l_h_angle_ph->Draw();
 	Q2_245_l_l_angle_ph->Draw("same");


 	TCanvas* c1 = (TCanvas*)f1->Get("t_phi_check");
  	TGraphPolar* Q2_245_l_l_gp = (TGraphPolar*) c1->GetListOfPrimitives()->FindObject("Graph");

	TCanvas* c2 = (TCanvas*)f2->Get("t_phi_check");
  	TGraphPolar* Q2_245_l_h_gp = (TGraphPolar*) c2->GetListOfPrimitives()->FindObject("Graph");


  	TCanvas * p3 = new TCanvas("p3", "p3", 600, 600);

	Q2_245_l_l_gp->SetMarkerColor(4);
	Q2_245_l_h_gp->SetMarkerColor(2);

	Q2_245_l_l_gp->Draw("AOP");

	p3->Update();

	Q2_245_l_l_gp->GetPolargram()->SetNdivPolar(4);
 	Q2_245_l_l_gp->GetPolargram()->SetNdivRadial(5);
 	Q2_245_l_l_gp->GetPolargram()->SetRangeRadial(0, 0.5); 
 	Q2_245_l_h_gp->GetPolargram()->SetRangeRadial(0, 0.5); 


	Q2_245_l_h_gp->Draw("OPsame");
	p3->Update();




/*--------------------------------------------------*/
/*--------------------------------------------------*/
/*--------------------------------------------------*/


	TFile *f1 = TFile::Open("yields.omega_245_55_-3000.target.root");
	TFile *f2 = TFile::Open("yields.omega_245_55_+0000.target.root");
	TFile *f3 = TFile::Open("yields.omega_245_55_+3000.target.root");
 

  	TH2F* Q2_245_h_l_angle_tmp = (TH2F*) f1->Get("t_dia");
  	TH2F* Q2_245_h_l_angle = (TH2F*) Q2_245_h_l_angle_tmp->Clone();
 
  	TH2F* Q2_245_h_0_angle_tmp = (TH2F*) f2->Get("t_dia");
  	TH2F* Q2_245_h_0_angle = (TH2F*) Q2_245_h_0_angle_tmp->Clone();


  	TH2F* Q2_245_h_h_angle_tmp = (TH2F*) f3->Get("t_dia");
  	TH2F* Q2_245_h_h_angle = (TH2F*) Q2_245_h_h_angle_tmp->Clone();


  	TCanvas * d4 = new TCanvas();

  	Q2_245_h_0_angle->SetLineColor(4);
 	Q2_245_h_h_angle->SetLineColor(2);
 
 	Q2_245_h_0_angle->Draw();
  	Q2_245_h_l_angle->Draw("same");
   	Q2_245_h_h_angle->Draw("same");


  	TH1F* Q2_245_h_l_angle_ph_tmp = (TH1F*) f1->Get("phi_pq_dia");
  	TH1F* Q2_245_h_l_angle_ph = (TH1F*) Q2_245_h_l_angle_ph_tmp->Clone();

  	TH1F* Q2_245_h_0_angle_ph_tmp = (TH1F*) f2->Get("phi_pq_dia");
  	TH1F* Q2_245_h_0_angle_ph = (TH1F*) Q2_245_h_0_angle_ph_tmp->Clone();
 
  	TH1F* Q2_245_h_h_angle_ph_tmp = (TH1F*) f3->Get("phi_pq_dia");
  	TH1F* Q2_245_h_h_angle_ph = (TH1F*) Q2_245_h_h_angle_ph_tmp->Clone();

  	TCanvas * t4 = new TCanvas("t4", "t4");

 
 	Q2_245_h_0_angle_ph->SetLineColor(4);
	Q2_245_h_h_angle_ph->SetLineColor(2);
 
 	Q2_245_h_h_angle_ph->Draw();
  	Q2_245_h_l_angle_ph->Draw("same");
 	Q2_245_h_0_angle_ph->Draw("same");



 	TCanvas* c1 = (TCanvas*)f1->Get("t_phi_check");
  	TGraphPolar* Q2_245_h_l_gp = (TGraphPolar*) c1->GetListOfPrimitives()->FindObject("Graph");

	TCanvas* c2 = (TCanvas*)f2->Get("t_phi_check");
  	TGraphPolar* Q2_245_h_0_gp = (TGraphPolar*) c2->GetListOfPrimitives()->FindObject("Graph");

	TCanvas* c3 = (TCanvas*)f3->Get("t_phi_check");
  	TGraphPolar* Q2_245_h_h_gp = (TGraphPolar*) c3->GetListOfPrimitives()->FindObject("Graph");


  	TCanvas * p4 = new TCanvas("p4", "p4", 600, 600);

	Q2_245_h_0_gp->SetMarkerColor(4);
	Q2_245_h_h_gp->SetMarkerColor(2);

	Q2_245_h_0_gp->Draw("AOP");

	p4->Update();

	Q2_245_h_0_gp->GetPolargram()->SetNdivPolar(4);
 	Q2_245_h_0_gp->GetPolargram()->SetNdivRadial(5);
 	Q2_245_h_0_gp->GetPolargram()->SetRangeRadial(0, 0.5); 

 	Q2_245_h_h_gp->GetPolargram()->SetRangeRadial(0, 0.5); 
 	Q2_245_h_l_gp->GetPolargram()->SetRangeRadial(0, 0.5); 

	Q2_245_h_h_gp->Draw("OPsame");
	Q2_245_h_l_gp->Draw("OPsame");

	p4->Update();





	







/*--------------------------------------------------*/
/*--------------------------------------------------*/

  	TCanvas * t_can = new TCanvas("t", "t", 800, 600);

	t_can->Divide(2,2);

	t_can->cd(1);
	d1->DrawClonePad();

	t_can->cd(2);
	d2->DrawClonePad();

	t_can->cd(3);
	d3->DrawClonePad();

	t_can->cd(4);
	d4->DrawClonePad();

	t_can->Update();



/*--------------------------------------------------*/
/*--------------------------------------------------*/

  	TCanvas * phi_can = new TCanvas("phi", "phi", 800, 600);

	phi_can->Divide(2,2);

 	phi_can->cd(1);
 	t1->DrawClonePad();
 
 	phi_can->cd(2);
 	t2->DrawClonePad();
 
 	phi_can->cd(3);
 	t3->DrawClonePad();
 
 	phi_can->cd(4);
 	t4->DrawClonePad();
 
 	phi_can->Update();



/*--------------------------------------------------*/
/*--------------------------------------------------*/

  	TCanvas * polar = new TCanvas("polar", "polar", 800, 800);

	polar->Divide(2,2);

	polar->cd(1);
	p1->DrawClonePad();

	polar->cd(2);
	p2->DrawClonePad();

	polar->cd(3);
	p3->DrawClonePad();

	polar->cd(4);
	p4->DrawClonePad();

	polar->Update();



}




