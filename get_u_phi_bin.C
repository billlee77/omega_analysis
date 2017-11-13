

TString bin_file = "../u_bin_interval" ;

Int_t Get_u_bin(){

	Float_t Q2_set;
	Int_t u_bin, phi_bin;


	ifstream u_phi_file;
	u_phi_file.open(bin_file);


	if (u_phi_file.is_open()) {
		u_phi_file >> Q2_set >> u_bin >> phi_bin;
	} else {
		cerr << "File u_bin_interval doesn't exist! " << endl;
		exit(0);
	}

//	cout << Q2_set << "  " << u_bin << "  " << phi_bin << endl; 
	u_phi_file.close();
	return u_bin;

}

Int_t Get_phi_bin(){

	Float_t Q2_set;
	Int_t u_bin, phi_bin;

	ifstream u_phi_file;
	u_phi_file.open(bin_file);

	if (u_phi_file.is_open()) {
		u_phi_file >> Q2_set >> u_bin >> phi_bin;
	} else {
		cerr << "File u_bin_interval doesn't exist! " << endl;
		exit(0);
	}

//	cout << Q2_set << "  " << u_bin << "  " << phi_bin << endl; 
	u_phi_file.close();
	return phi_bin;

}
