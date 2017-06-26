{

	gROOT->ProcessLine(".x $CUORE_INSTALL/lib/diana_root.C");

	QChain * qtree = new QChain("qtree");
	
	//ulite
	TString file_path = "file path in the cluster";
	TString process_type = "OfficialProcessed_v02.30";
	file_path += process_type;
	file_path += "/";

	TString ds_name = "2079";
	
	Int_t ds_index = 0;
	Bool_t b_reduced = 1;

	if(ds_index == 11) TString list_name = "All_reduced_calib.list";
	else {
	        file_path += "ds2079/";
		if(b_reduced ==1){
			TString list_name = "calib_BlindedReduced_ds2079.list";
		}
		else TString list_name = "calib_Blinded_ds2079.list";
	}

	TString file_name =  file_path+list_name; 

	qtree->Add(file_name.Data());
	cout<<"Adding Root Files ... "<<file_name.Data()<<endl;
}
