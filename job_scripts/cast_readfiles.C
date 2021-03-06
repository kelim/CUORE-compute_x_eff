{

	gROOT->ProcessLine(".x $CUORE_INSTALL/lib/diana_root.C");

	QChain * qtree = new QChain("qtree");
	
	//ulite
	TString file_path = "file path in the cluster";
	TString process_type = "OfficialProcessed_v02.30";
	file_path += process_type;
	file_path += "/";

	TString ds_name = "_ds_name_";
	
	Int_t ds_index = _ds_index_;
	Bool_t b_reduced = _b_reduced_;

	if(ds_index == _ds_index_all_) TString list_name = "All_reduced__data_type_.list";
	else {
	        file_path += "ds_ds_name_/";
		if(b_reduced ==1){
			TString list_name = "_data_type__BlindedReduced_ds_ds_name_.list";
		}
		else TString list_name = "_data_type__Blinded_ds_ds_name_.list";
	}

	TString file_name =  file_path+list_name; 

	qtree->Add(file_name.Data());
	cout<<"Adding Root Files ... "<<file_name.Data()<<endl;
}
