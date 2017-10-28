#include <iostream>
#include "print_data.h"



using namespace std;




File_Output::File_Output() {



	static int tsttst = 0;

	if (tsttst == 0) {
		cout << "/*--------------------------------------------------*/" << endl;
		cout << "                  Output File Created                 " << endl << endl;
//		cout << "/*--------------------------------------------------*/" << endl;
		tsttst++;
	}

	out_dir = "file_out/";

}

File_Output::~File_Output() {}


TFile* File_Output::Create_Out_File(TString fname){

	/*--------------------------------------------------*/
	/// Output file name

	TFile* file_out;


	out_file_name =  out_dir + fname + ".root"; 
	file_out = new TFile(out_file_name, "RECREATE");

	return file_out;

}





