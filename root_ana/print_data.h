#ifndef PRINT_DATA_H
#define PRINT_DATA_H

#include "TFile.h"



class File_Output {



	public:

		File_Output();
		~File_Output();

		TFile* Create_Out_File(TString fname);

		TString Get_Dir() {return out_dir;}
		// TFile* file_out;


	private:
	
		/*--------------------------------------------------*/
		// File out and firectory
		TString out_file_name;
		TString out_dir; 
		
		

};

#endif
