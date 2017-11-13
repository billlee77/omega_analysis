{


//	TFile* file_out = new TFile("", );

	TFile *f1 = TFile::Open("yields.omega_160_32_+0970.target.root");
	TFile *f2 = TFile::Open("yields.omega_160_32_+3000.target.root");

	TCanvas* c1 = (TCanvas*)f1->Get("diamand");

	TCanvas* c2 = (TCanvas*)f2->Get("diamand");

	c1->Draw();
	c2->Draw();

	TMultiGraph * mg1_tmp = (TMultiGraph*) c1->GetListOfPrimitives()->FindObject("diamond_gr");
	TMultiGraph * mg1 = (TMultiGraph*) mg1_tmp->Clone();

	

	TMultiGraph * mg2_tmp = (TMultiGraph*) c2->GetListOfPrimitives()->FindObject("diamond_gr");
	TMultiGraph * mg2 = (TMultiGraph*) mg2_tmp->Clone();



 	delete c1;
 	delete c2;



	TCanvas * d1 = new TCanvas();
	
// 	mg1->Draw("AP");
// 
// 	//mg2->SetMarkerColor(2);
// 
// 	mg2->Draw("P");
// 	
// 
	cout << mg1->GetListOfGraphs()->GetSize() << endl;


	TList* list1 = (TList*) mg1->GetListOfGraphs();
	TList* list2 = (TList*) mg2->GetListOfGraphs();

	d1->Update();

	TMultiGraph* m_160_low = new TMultiGraph();

	TGraph* gg1 = new TGraph();
	gg1->Merge(list1);
//	gg1->Draw("P");

	TGraph* gg2 = new TGraph();
	gg2->Merge(list2);
//	gg2->Draw("P");

	gg1->SetMarkerColor(2);
	gg2->SetMarkerColor(2);
	
	m_160_low->Add(gg1);
	m_160_low->Add(gg2);



	
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/


	TFile *f1 = TFile::Open("yields.omega_160_59_-2730.target.root");
	TFile *f2 = TFile::Open("yields.omega_160_59_+0000.target.root");
	TFile *f3 = TFile::Open("yields.omega_160_59_+3000.target.root");

	TCanvas* c1 = (TCanvas*)f1->Get("diamand");

	TCanvas* c2 = (TCanvas*)f2->Get("diamand");

	TCanvas* c3 = (TCanvas*)f3->Get("diamand");

	c1->Draw();
	c2->Draw();

	TMultiGraph * mg1_tmp = (TMultiGraph*) c1->GetListOfPrimitives()->FindObject("diamond_gr");
	TMultiGraph * mg1 = (TMultiGraph*) mg1_tmp->Clone();

	TMultiGraph * mg2_tmp = (TMultiGraph*) c2->GetListOfPrimitives()->FindObject("diamond_gr");
	TMultiGraph * mg2 = (TMultiGraph*) mg2_tmp->Clone();

	TMultiGraph * mg3_tmp = (TMultiGraph*) c3->GetListOfPrimitives()->FindObject("diamond_gr");
	TMultiGraph * mg3 = (TMultiGraph*) mg2_tmp->Clone();

 
	delete c1;
	delete c2;
 	delete c3;


	
	TMultiGraph* m_160_high = new TMultiGraph();




	TList* list1 = (TList*) mg1->GetListOfGraphs();
	TList* list2 = (TList*) mg2->GetListOfGraphs();
	TList* list3 = (TList*) mg3->GetListOfGraphs();

	d1->Update();

	TGraph* ggg1 = new TGraph();
	ggg1->Merge(list1);
//	ggg1->Draw("AP");

	TGraph* ggg2 = new TGraph();
	ggg2->Merge(list2);
//	ggg2->Draw("P");

	TGraph* ggg3 = new TGraph();
	ggg3->Merge(list2);
//	ggg3->Draw("P");


	m_160_high->Add(ggg1);
	m_160_high->Add(ggg2);
	m_160_high->Add(ggg3);


	d1->Update();

	m_160_high->Draw("AP");
	m_160_low->Draw("P");
	
	d1->Update();







/*--------------------------------------------------*/
/*--------------------------------------------------*/


	TCanvas * d2 = new TCanvas();


//	TFile* file_out = new TFile("", );

	TFile *f1 = TFile::Open("yields.omega_245_27_+1350.target.root");
	TFile *f2 = TFile::Open("yields.omega_245_27_+3000.target.root");

	TCanvas* c1 = (TCanvas*)f1->Get("diamand");
	TCanvas* c2 = (TCanvas*)f2->Get("diamand");

	c1->Draw();
	c2->Draw();

	TMultiGraph * mg1_tmp = (TMultiGraph*) c1->GetListOfPrimitives()->FindObject("diamond_gr");
	TMultiGraph * mg1 = (TMultiGraph*) mg1_tmp->Clone();

	

	TMultiGraph * mg2_tmp = (TMultiGraph*) c2->GetListOfPrimitives()->FindObject("diamond_gr");
	TMultiGraph * mg2 = (TMultiGraph*) mg2_tmp->Clone();



 	delete c1;
 	delete c2;



	
// 	mg1->Draw("AP");
// 
// 	//mg2->SetMarkerColor(2);
// 
// 	mg2->Draw("P");
// 	
// 
	cout << mg1->GetListOfGraphs()->GetSize() << endl;


	TList* list1 = (TList*) mg1->GetListOfGraphs();
	TList* list2 = (TList*) mg2->GetListOfGraphs();


	TMultiGraph* m_245_low = new TMultiGraph();

	TGraph* gg11 = new TGraph();
	gg11->Merge(list1);
//	gg1->Draw("P");

	TGraph* gg22 = new TGraph();
	gg22->Merge(list2);
//	gg2->Draw("P");

	gg11->SetMarkerColor(2);
	gg22->SetMarkerColor(2);
	
	m_245_low->Add(gg11);
	m_245_low->Add(gg22);




/*--------------------------------------------------*/
/*--------------------------------------------------*/



	TFile *f1 = TFile::Open("yields.omega_245_55_-3000.target.root");
	TFile *f2 = TFile::Open("yields.omega_245_55_+0000.target.root");
	TFile *f3 = TFile::Open("yields.omega_245_55_+3000.target.root");

	TCanvas* c1 = (TCanvas*)f1->Get("diamand");

	TCanvas* c2 = (TCanvas*)f2->Get("diamand");

	TCanvas* c3 = (TCanvas*)f3->Get("diamand");

	c1->Draw();
	c2->Draw();

	TMultiGraph * mg1_tmp = (TMultiGraph*) c1->GetListOfPrimitives()->FindObject("diamond_gr");
	TMultiGraph * mg1 = (TMultiGraph*) mg1_tmp->Clone();

	TMultiGraph * mg2_tmp = (TMultiGraph*) c2->GetListOfPrimitives()->FindObject("diamond_gr");
	TMultiGraph * mg2 = (TMultiGraph*) mg2_tmp->Clone();

	TMultiGraph * mg3_tmp = (TMultiGraph*) c3->GetListOfPrimitives()->FindObject("diamond_gr");
	TMultiGraph * mg3 = (TMultiGraph*) mg2_tmp->Clone();

 
	delete c1;
	delete c2;
 	delete c3;




	TMultiGraph* m_245_high = new TMultiGraph();

	TList* list1 = (TList*) mg1->GetListOfGraphs();
	TList* list2 = (TList*) mg2->GetListOfGraphs();
	TList* list3 = (TList*) mg3->GetListOfGraphs();


	TGraph* gggg1 = new TGraph();
	gggg1->Merge(list1);
//	gggg1->Draw("P");

	TGraph* gggg2 = new TGraph();
	gggg2->Merge(list2);
//	gggg2->Draw("P");

	TGraph* gggg3 = new TGraph();
	gggg3->Merge(list2);
//	gggg3->Draw("P");


	m_245_high->Add(gggg1);
	m_245_high->Add(gggg2);
	m_245_high->Add(gggg3);


	d2->cd();
	d2->Update();




	m_245_high->Draw("AP");
	m_245_low->Draw("P");
	
	d2->Update();




	TMultiGraph* m_245_high_cl = (TMultiGraph*)m_245_high->Clone();


	TCanvas* d3 = new TCanvas();



	m_245_high_cl->Draw("AP");

	m_245_high_cl->GetXaxis()->SetTitle("Q^{2}");
	m_245_high_cl->GetXaxis()->CenterTitle();

	m_245_high_cl->GetYaxis()->SetTitle("W");
	m_245_high_cl->GetYaxis()->CenterTitle();

	
	m_245_high_cl->GetXaxis()->SetLimits(0.8,3.5);

	m_245_low->Draw("P");
	m_160_high->Draw("P");
	m_160_low->Draw("P");



	d3->Update();





}
