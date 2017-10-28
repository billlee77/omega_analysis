#include <iostream>

#include "general_utility.h"





Int_t SizeOfArray(Double_t* tar_arr) {

	int size_of_array = 0;

	while( int (*tar_arr) > 0 ) {
	    tar_arr++;
		size_of_array++;
	}

	return size_of_array;

}



Int_t* Array_D_to_I (Double_t* tar_arr, Int_t array_size) {

	Int_t* int_array_tmp = new Int_t[array_size];
 	
 	for(int i=0; i < array_size ;i++) {
  
		int_array_tmp[i] = (int)tar_arr[i];

 	}

	return int_array_tmp;

}





