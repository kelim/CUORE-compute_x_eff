{
	TString current_path = "path in the cluster";
	TString s_readfile = current_path+"job_scripts/"+"readfiles_1.C";

	gROOT->Macro(s_readfile.Data());
        gROOT->SetBatch(1);

	////// Set Variables to draw histogram
	Double_t energy_max, energy_min;
	Double_t h_energy_binsize;
	Int_t ch_number;
	TString draw_branch;
	TString peak_name;
	Int_t n_multiplicity;
	Int_t n_sigma_for_peak_count, n_sigma_for_bg_count;
	Int_t peak_index, draw_branch_index;
	
	peak_name = "K-40";
	energy_max = 1480;	energy_min = 1440;
	ch_number = 2;
	h_energy_binsize = 0.500000;
	draw_branch = "Energy";
	n_multiplicity = 1;
	n_sigma_for_peak_count = 5;
	n_sigma_for_bg_count = 15;
	peak_index = 2;
	draw_branch_index = 0;

	//Default analysis cuts
	TString s_setcuts = current_path+"setcuts.C";
	gROOT->Macro(s_setcuts.Data());

	////////// Set the cuts to draw the histogram

        TCut x_energy = "Energy > 1440 && Energy< 1480";
	TCut x_tot_energy = "TotalEnergy_RawTime > 500.000000";

	//Good event selection with OT_ChiSquare selection
	TCut x_ot_chisq = "Multiplicity_RawTime > 0";

	//Mutliplicity cut
	TCut x_multiplicity = "Multiplicity_RawTime == 1";

	//Timing cut for each time bin
	TCut x_time ="";

	// Quality cuts to choose the clean sample for the evaluation of the cut efficiency
	//	TCut x_base = x_signal+x_bi+x_energy+x_channel+x_run+x_ot_chisq;
	TCut x_base = x_signal+x_bi+x_energy+x_run+x_ot_chisq;
	TCut x_test = x_base+x_multiplicity;

	cout<<"x_base: "<<x_base.Print();cout<<endl;
        cout<<"x_test: "<<x_test.Print(); cout<<endl;

	//Output file name to write
        TString output_file_name = current_path + "results/";
        output_file_name += "x_eff_";
        output_file_name += process_type; //from readfiles.C
	output_file_name +="_";
        output_file_name += list_name; //from readfiles.C
        cout<<output_file_name.Data() <<endl;

        //open file to write the outputs
	FILE *fout = fopen(output_file_name.Data(), "a+");

	//	TCut x_draw = x_channel+x_energy+x_signal+x_bi+x_multiplicity;
	TCut x_draw = x_base;
	if(n_multiplicity>1) TCut x_draw = x_base+x_multiplicity+x_tot_energy+x_ot_chisq;
	cout<<"x_draw: "<<x_draw.Print(); cout<<endl;


	// fit energy spectrum for each channel
        Double_t h_energy_xbin_size = h_energy_binsize;//keV
        Double_t h_energy_xmin = energy_min;    Double_t h_energy_xmax = energy_max;
        Int_t n_h_energy_xbin = (Int_t) ((h_energy_xmax - h_energy_xmin)/h_energy_xbin_size);
        cout<<"h_energy_xmin: "<< h_energy_xmin <<endl;
        cout<<"h_energy_xbin_size: "<<h_energy_xbin_size<<" keV"<<endl;

	/////////// Set Histogram w/ base cuts
        TString h_energy_name = "h_energy_";
        h_energy_name += peak_name;
        h_energy_name += "; "; h_energy_name += draw_branch; h_energy_name += " [keV]";
        h_energy_name += "; Rate [cts/";
        h_energy_name += h_energy_binsize;
        h_energy_name += "keV]";
        TString h_energy_draw_branch = "";
        h_energy_draw_branch += draw_branch;
        h_energy_draw_branch += ">>h_energy";

        TH1D * h_energy = new TH1D("h_energy", h_energy_name.Data(), n_h_energy_xbin, 
				       h_energy_xmin, h_energy_xmax);
        qtree->Draw(h_energy_draw_branch.Data(),x_base);

        //Set Histogram w/ test cuts
        TString h_energy_cut_name = "h_energy_cut_";
        h_energy_cut_name += peak_name;
        h_energy_cut_name += ";"; h_energy_cut_name += draw_branch; h_energy_cut_name += " [keV]";
        h_energy_cut_name += "; Rate [cts/";
        h_energy_cut_name += h_energy_binsize;
        h_energy_cut_name += "keV]";
        TString h_energy_cut_draw_branch = "";
        h_energy_cut_draw_branch += draw_branch;
        h_energy_cut_draw_branch += ">>h_energy_cut";

        TH1D * h_energy_cut = new TH1D("h_energy_cut",h_energy_cut_name.Data(),n_h_energy_xbin, h_energy_xmin,h_energy_xmax);
        qtree->Draw(h_energy_cut_draw_branch.Data(),x_test);

	//Get Efficiency histogram
        TEfficiency * pEff = new TEfficiency(*h_energy_cut,*h_energy);

        Double_t h_energy_zoom_xmin = energy_min;
        Double_t h_energy_zoom_xmax = energy_max;
        Int_t n_h_energy_zoom_xbin = (Int_t) ((h_energy_zoom_xmax - h_energy_zoom_xmin)/h_energy_xbin_size);
        cout<<"h_energy_zoom_xmin: "<<h_energy_zoom_xmin<<endl;
        cout<<"h_energy_zoom_xmax: "<<h_energy_zoom_xmax<<endl;

        //Set Histogram for zoom
        // w/ basic cuts
        TString h_energy_zoom_name = "h_energy_zoom_";
        h_energy_zoom_name += peak_name;
        h_energy_zoom_name += ";"; h_energy_zoom_name += draw_branch; h_energy_zoom_name += " [keV]";
        h_energy_zoom_name += "; Rate [cts/";
        h_energy_zoom_name += h_energy_binsize;
        h_energy_zoom_name += "keV]";
        TString h_energy_zoom_draw_branch = "";
        h_energy_zoom_draw_branch += draw_branch;
        h_energy_zoom_draw_branch +=">>h_energy_zoom";

        TH1D * h_energy_zoom = new TH1D("h_energy_zoom",h_energy_zoom_name.Data(),n_h_energy_zoom_xbin, h_energy_zoom_xmin,h_energy_zoom_xmax);
        qtree->Draw(h_energy_zoom_draw_branch.Data(),x_base);
        h_energy_zoom->Sumw2();

	// w/ test cuts
        TString h_energy_zoom_cut_name = "h_energy_zoom_cut_";
        h_energy_zoom_cut_name += peak_name;
        h_energy_zoom_cut_name += ";"; h_energy_zoom_cut_name += draw_branch; h_energy_zoom_cut_name += " [keV]";
        h_energy_zoom_cut_name += "; Rate [cts/";
        h_energy_zoom_cut_name += h_energy_binsize;
        h_energy_zoom_cut_name += "keV]";
        TString h_energy_zoom_cut_draw_branch = "";
        h_energy_zoom_cut_draw_branch += draw_branch;
        h_energy_zoom_cut_draw_branch +=">>h_energy_zoom_cut";

        TH1D * h_energy_zoom_cut = new TH1D("h_energy_zoom_cut",h_energy_zoom_cut_name.Data(),n_h_energy_zoom_xbin, h_energy_zoom_xmin,h_energy_zoom_xmax);
        qtree->Draw(h_energy_zoom_cut_draw_branch.Data(),x_test);
        h_energy_zoom_cut->Sumw2();


	//Set Fit Ranges
        Double_t n_sigma_fit = 3.;
        Double_t h_energy_fit_min = h_energy_xmin + (h_energy->GetMaximumBin()*h_energy_xbin_size) - 0.03*(h_energy->GetMaximumBin()*h_energy_xbin_size)*n_sigma_fit;
        Double_t h_energy_fit_max = h_energy_xmin + (h_energy->GetMaximumBin()*h_energy_xbin_size) + 0.03*(h_energy->GetMaximumBin()*h_energy_xbin_size)*n_sigma_fit;
													   
        cout<<"h_energy_fit_min: "<<h_energy_fit_min<<endl;
        cout<<"h_energy_fit_max: "<<h_energy_fit_max<<endl;

	TF1 *f_gaus = new TF1("f_gaus","gaus",h_energy_fit_min,h_energy_fit_max);
        f_gaus->SetLineColor(kMagenta);

        //fitting histgoram w/ basic cuts
        h_energy->Fit(f_gaus,"R");
        Double_t f_gaus_getchisq = f_gaus->GetChisquare();
        Double_t f_gaus_getndf = f_gaus->GetNDF();
        Double_t f_gaus_reduced_chisq = f_gaus_getchisq/f_gaus_getndf;
        Double_t f_gaus_prob = f_gaus->GetProb();
        Double_t f_gaus_fit_par0 = f_gaus->GetParameter(0);
        Double_t f_gaus_fit_par0_err = f_gaus->GetParError(0);
        Double_t f_gaus_fit_par1 = f_gaus->GetParameter(1);
        Double_t f_gaus_fit_par1_err = f_gaus->GetParError(1);
        Double_t f_gaus_fit_par2 = f_gaus->GetParameter(2);
        Double_t f_gaus_fit_par2_err = f_gaus->GetParError(2);

        //fitting histgoram w/ test cuts
        TF1 *f_gaus_cut = new TF1("f_gaus_cut","gaus",h_energy_fit_min,h_energy_fit_max);
        f_gaus_cut->SetLineColor(kOrange+8);

        h_energy_zoom_cut->Fit(f_gaus_cut,"R");
        Double_t f_gaus_cut_getchisq = f_gaus_cut->GetChisquare();
        Double_t f_gaus_cut_getndf = f_gaus_cut->GetNDF();
        Double_t f_gaus_cut_reduced_chisq = f_gaus_cut_getchisq/f_gaus_cut_getndf;
        Double_t f_gaus_cut_pvalue = f_gaus_cut->GetProb();
        Double_t f_gaus_cut_fit_par0 = f_gaus_cut->GetParameter(0);
        Double_t f_gaus_cut_fit_par0_err = f_gaus_cut->GetParError(0);
        Double_t f_gaus_cut_fit_par1 = f_gaus_cut->GetParameter(1);
        Double_t f_gaus_cut_fit_par1_err = f_gaus_cut->GetParError(1);
        Double_t f_gaus_cut_fit_par2 = f_gaus_cut->GetParameter(2);
        Double_t f_gaus_cut_fit_par2_err = f_gaus_cut->GetParError(2);


        //Compute residuals of the fit
	//	Int_t np = (Int_t) ((h_energy_fit_max - h_energy_fit_min)/h_energy_xbin_size);
	Int_t np = (Int_t) ((h_energy_fit_max - h_energy_fit_min)/h_energy_xbin_size)+1;
	cout<<"np: "<<np<<endl;
        TGraph * g_residual = new TGraph();     TGraph * g_residual_cut = new TGraph();
        g_residual->SetTitle("Residual");       g_residual_cut->SetTitle("Residual_cut");
        Double_t get_f_gaus, get_f_gaus_cut;
        Double_t get_h_energy_zoom, get_h_energy_zoom_err;      Double_t get_h_energy_zoom_cut, get_h_energy_zoom_cut_err;
        Double_t g_residual_x , g_residual_y;   Double_t g_residual_cut_x , g_residual_cut_y;

        int ij = 0;     int ik = 0;
        for (int ii =0; ii<np;ii++){
                g_residual_x = h_energy_zoom->GetBinLowEdge(h_energy_zoom->FindBin(h_energy_fit_min))+(ii+0.5)*h_energy_xbin_size;
                get_f_gaus = f_gaus->Eval(g_residual_x);
                get_h_energy_zoom = h_energy_zoom->GetBinContent(h_energy_zoom->FindBin(h_energy_fit_min)+ii);
                get_h_energy_zoom_err = h_energy_zoom->GetBinError(h_energy_zoom->FindBin(h_energy_fit_min)+ii);
                if(get_h_energy_zoom_err>0){
                        g_residual_y = (get_h_energy_zoom-get_f_gaus)/get_h_energy_zoom_err;
                        g_residual->SetPoint(ij, g_residual_x,g_residual_y);
                        //                      cout<<"ii :" <<ii;
                        //                      cout<<"ij :" <<ij;
                        //                      cout<<"g_residual_x: "<<g_residual_x<<endl;
                        //                      cout<<"get_f_gaus: "<<get_f_gaus <<endl;
                        //                      cout<<"get_h_energy_zoom: "<<get_h_energy_zoom<<endl;
                        //                      cout<<"g_residual_y: "<<g_residual_y<<endl;
                        //                      cin.get();
                        ij++;
                }//if

                g_residual_cut_x = h_energy_zoom_cut->GetBinLowEdge(h_energy_zoom_cut->FindBin(h_energy_fit_min))+(ii+0.5)*h_energy_xbin_size;
                get_f_gaus_cut = f_gaus_cut->Eval(g_residual_cut_x);
                get_h_energy_zoom_cut = h_energy_zoom->GetBinContent(h_energy_zoom_cut->FindBin(h_energy_fit_min)+ii);
                get_h_energy_zoom_cut_err = h_energy_zoom->GetBinError(h_energy_zoom_cut->FindBin(h_energy_fit_min)+ii);
                if(get_h_energy_zoom_cut_err>0){
                        g_residual_cut_y = (get_h_energy_zoom_cut-get_f_gaus_cut)/get_h_energy_zoom_cut_err;
                        g_residual_cut->SetPoint(ik, g_residual_cut_x,g_residual_cut_y);
                        //                      cout<<"ii :" <<ii;
                        //                      cout<<"ij :" <<ij;
                        //                      cout<<"g_residual_x: "<<g_residual_x<<endl;
                        //                      cout<<"get_f_gaus: "<<get_f_gaus <<endl;
                        //                      cout<<"get_h_energy_zoom: "<<get_h_energy_zoom<<endl;
                        //                      cout<<"g_residual_y: "<<g_residual_y<<endl;
                        ik++;
                }//if

        }//for graph

        g_residual->GetXaxis()->SetRangeUser(h_energy_fit_min,h_energy_fit_max);
        g_residual_cut->GetXaxis()->SetRangeUser(h_energy_fit_min,h_energy_fit_max);

        //To compute the number of signal and noise events in the peak, noise is assumed to be flat/slwoly varying
        //-----------------  Method 1, count number of events and subtract width-normalized flat background
        Int_t h_energy_peak_min_bin = h_energy_zoom->GetMaximumBin() - (Int_t) (n_sigma_for_peak_count*f_gaus_fit_par2/h_energy_xbin_size);
        Int_t h_energy_peak_max_bin = h_energy_zoom->GetMaximumBin() + (Int_t) (n_sigma_for_peak_count*f_gaus_fit_par2/h_energy_xbin_size);
        cout<<"h_energy_peak_min_bin: "<<h_energy_peak_min_bin<<endl;
        cout<<"h_energy_peak_max_bin: "<<h_energy_peak_max_bin<<endl;

        Double_t h_energy_peak_min = h_energy_zoom->GetBinLowEdge(h_energy_peak_min_bin);
        Double_t h_energy_peak_max = h_energy_zoom->GetBinLowEdge(h_energy_peak_max_bin+1);
        cout<<"h_energy_peak_min: "<<h_energy_peak_min<<endl;
        cout<<"h_energy_peak_max: "<<h_energy_peak_max<<endl;

        Double_t n_h_energy_peak = h_energy_zoom->Integral(h_energy_peak_min_bin, h_energy_peak_max_bin);
        Double_t n_h_energy_cut_peak = h_energy_zoom_cut->Integral(h_energy_peak_min_bin, h_energy_peak_max_bin);
        cout<<"n_h_energy_peak: "<<n_h_energy_peak<<endl;
        cout<<"n_h_energy_cut_peak: "<<n_h_energy_cut_peak<<endl;

        //width of peak before/after
        Double_t n_sigma = n_sigma_for_bg_count;
        cout<<"h_energy_peak_max+ n_sigma*f_gaus_fit_par2: "<<h_energy_peak_max+ n_sigma*f_gaus_fit_par2<<endl;
        cout<<"h_energy_zoom_max: "<<energy_max<<endl;

        if(h_energy_peak_max+ n_sigma*f_gaus_fit_par2 < energy_max){
                Int_t nbin_width_peak_before = (Int_t) (n_sigma*f_gaus_fit_par2/h_energy_xbin_size);
                Int_t nbin_width_peak_after = nbin_width_peak_before;
        }
        else{
                Int_t nbin_width_peak_before = (Int_t) ((energy_max - h_energy_peak_max)/h_energy_xbin_size)-1;
                Int_t nbin_width_peak_after = nbin_width_peak_before;
        }

        Double_t h_energy_peak_before_min = h_energy_zoom->GetBinLowEdge(h_energy_peak_min_bin-nbin_width_peak_before);
        Double_t h_energy_peak_before_max = h_energy_zoom->GetBinLowEdge(h_energy_peak_min_bin);
        Double_t h_energy_peak_after_min = h_energy_zoom->GetBinLowEdge(h_energy_peak_max_bin+1);
        Double_t h_energy_peak_after_max = h_energy_zoom->GetBinLowEdge(h_energy_peak_max_bin+nbin_width_peak_after+1);

        Int_t n_h_energy_peak_before = h_energy_zoom->Integral(h_energy_peak_min_bin - nbin_width_peak_before,h_energy_peak_min_bin-1);
        Int_t n_h_energy_peak_after = h_energy_zoom->Integral(h_energy_peak_max_bin+1, h_energy_peak_max_bin+nbin_width_peak_after);
        Int_t n_h_energy_cut_peak_before = h_energy_zoom_cut->Integral(h_energy_peak_min_bin - nbin_width_peak_before,h_energy_peak_min_bin-1);
        Int_t n_h_energy_cut_peak_after = h_energy_zoom_cut->Integral(h_energy_peak_max_bin+1, h_energy_peak_max_bin+nbin_width_peak_after);

        cout<<"n_h_energy_peak_before: "<<n_h_energy_peak_before<<endl;
        cout<<"n_h_energy_peak_after: "<<n_h_energy_peak_after<<endl;
        cout<<"n_h_energy_cut_peak_before: "<<n_h_energy_cut_peak_before<<endl;
        cout<<"n_h_energy_cut_peak_after: "<<n_h_energy_cut_peak_after<<endl;

        const Double_t low_e_count = h_energy_peak_before_min;
        const Double_t high_e_count = h_energy_peak_after_max;
        TString tmp = draw_branch; tmp += " > "; tmp += low_e_count; tmp += " && "; tmp += draw_branch; tmp += " < "; tmp +=high_e_count;
        TCut x_e_peak_count = TCut(tmp.Data());
        cout<<"x_e_peak_count: "<<x_e_peak_count.Print();cout<<endl;
        Double_t n_get_entries = qtree->GetEntries(x_base+x_e_peak_count);
        Double_t n_sum = n_h_energy_peak_before+n_h_energy_peak+n_h_energy_peak_after;
        if(n_get_entries!=n_sum){
                cout<<"WARNING: n_get_entries and n_sum are NOT consistent!!!"<<endl;
                cout<<"n_get_entries: "<<n_get_entries <<", n_sum: "<<n_sum<<endl;
        }

        //normalized to the width
        Double_t n_h_energy_peak_bg=
                (n_h_energy_peak_before+n_h_energy_peak_after)*(h_energy_peak_max-h_energy_peak_min)/((Double_t) (nbin_width_peak_before+nbin_width_peak_after)*h_energy_xbin_size);
        Double_t n_h_energy_peak_signal = n_h_energy_peak- n_h_energy_peak_bg;
        Double_t h_energy_signal_to_total = n_h_energy_peak_signal/n_h_energy_peak;
        cout<<"n_h_energy_peak_bg: "<<n_h_energy_peak_bg<<endl;
        cout<<"n_h_energy_peak_signal: "<<n_h_energy_peak_signal<<endl;
        cout<<"h_energy_signal_to_total: "<<h_energy_signal_to_total<<endl;

        Double_t n_h_energy_cut_peak_bg=
                (n_h_energy_cut_peak_before+n_h_energy_cut_peak_after)*(h_energy_peak_max-h_energy_peak_min)/((Double_t) (nbin_width_peak_before+nbin_width_peak_after)*h_energy_xbin_size);
        Double_t n_h_energy_cut_peak_signal = n_h_energy_cut_peak- n_h_energy_cut_peak_bg;
        Double_t h_energy_cut_signal_to_total = n_h_energy_cut_peak_signal/n_h_energy_cut_peak;
        cout<<"n_h_energy_cut_peak_bg: "<<n_h_energy_cut_peak_bg<<endl;
        cout<<"n_h_energy_cut_peak_signal: "<<n_h_energy_cut_peak_signal<<endl;
        cout<<"h_energy_cut_signal_to_total: "<<h_energy_cut_signal_to_total<<endl;
        //      cin.get();

        Double_t n_h_energy_peak_signal_x_eff = n_h_energy_cut_peak_signal/n_h_energy_peak_signal;
        Double_t n_h_energy_peak_bg_x_eff = n_h_energy_cut_peak_bg/n_h_energy_peak_bg;
        Double_t n_h_energy_peak_bg_rej_eff =(n_h_energy_peak_bg-n_h_energy_cut_peak_bg)/n_h_energy_peak_bg;

        //compute uncertainty on the efficiency
        TEfficiency * x_eff = new TEfficiency();
        Double_t cl = 0.683;
        Double_t n_h_energy_peak_signal_x_eff_l = x_eff->ClopperPearson(n_h_energy_peak_signal,n_h_energy_cut_peak_signal,cl,0);
        Double_t n_h_energy_peak_signal_x_eff_h = x_eff->ClopperPearson(n_h_energy_peak_signal,n_h_energy_cut_peak_signal,cl,1);
        //      Double_t n_h_energy_peak_bg_x_eff_l = x_eff->ClopperPearson(n_h_energy_peak_bg,n_h_energy_cut_peak_bg,cl,0);
        //      Double_t n_h_energy_peak_bg_x_eff_h = x_eff->ClopperPearson(n_h_energy_peak_bg,n_h_energy_cut_peak_bg,cl,1);
        // The real measurement on the background have more stat
        Double_t n_h_energy_peak_bg_x_eff_l = x_eff->ClopperPearson(n_h_energy_peak_before+n_h_energy_peak_after,n_h_energy_cut_peak_before+n_h_energy_cut_peak_after,cl,0);
        Double_t n_h_energy_peak_bg_x_eff_h = x_eff->ClopperPearson(n_h_energy_peak_before+n_h_energy_peak_after,n_h_energy_cut_peak_before+n_h_energy_cut_peak_after,cl,1);

        Double_t n_h_energy_peak_bg_rej_eff_l = x_eff->ClopperPearson(n_h_energy_peak_bg,
                                                                      n_h_energy_peak_bg-n_h_energy_cut_peak_bg,cl,0);
        Double_t n_h_energy_peak_bg_rej_eff_h = x_eff->ClopperPearson(n_h_energy_peak_bg,
                                                                      n_h_energy_peak_bg-n_h_energy_cut_peak_bg,cl,1);
        cout<<"n_h_energy_peak_bg_x_eff: "<<n_h_energy_peak_bg_x_eff<< "[ "
            << n_h_energy_peak_bg_x_eff_h <<", "
            << n_h_energy_peak_bg_x_eff_l <<" ]"
            <<endl;

        cout<<"n_h_energy_peak_signal_x_eff: "<<n_h_energy_peak_signal_x_eff<< " [ "
            << n_h_energy_peak_signal_x_eff_h  <<", "
            << n_h_energy_peak_signal_x_eff_l <<" ]"
            <<endl;

        //--------------- Method 2, use gaus + flat background fit -----------------------//
        TF1 * f_gaus_pol1 = new TF1("f_gaus_pol1","gaus(0)+pol1(3)",h_energy_peak_before_min,h_energy_peak_after_max);
        Double_t f_gaus_pol1_par[5] = {f_gaus_fit_par0*0.5,f_gaus_fit_par1, f_gaus_fit_par2,n_h_energy_peak_bg/(h_energy_peak_max-h_energy_peak_min), 0.};
        f_gaus_pol1->SetParameters(f_gaus_pol1_par);
        f_gaus_pol1->SetLineColor(kMagenta);
        f_gaus_pol1->SetLineStyle(2);

        TH1D *h_energy_zoom_2 = (TH1D*)h_energy_zoom->Clone("h_energy_zoom_2");
        h_energy_zoom_2->Scale(1./h_energy_xbin_size);
        h_energy_zoom_2->SetXTitle("Energy [keV]");
        h_energy_zoom_2->SetYTitle("Rate [cts/keV]");
        h_energy_zoom_2->Fit(f_gaus_pol1,"R");
        Double_t f_gaus_pol1_getchisq = f_gaus_pol1->GetChisquare();
        Double_t f_gaus_pol1_getndf = f_gaus_pol1->GetNDF();
        Double_t f_gaus_pol1_reduced_chisq = f_gaus_pol1_getchisq/f_gaus_pol1_getndf;
        Double_t f_gaus_pol1_pvalue = f_gaus_pol1->GetProb();
        Double_t f_gaus_pol1_fit_par0 = f_gaus_pol1->GetParameter(0);
        Double_t f_gaus_pol1_fit_par0_err = f_gaus_pol1->GetParError(0);
        Double_t f_gaus_pol1_fit_par1 = f_gaus_pol1->GetParameter(1);
        Double_t f_gaus_pol1_fit_par1_err = f_gaus_pol1->GetParError(1);
        Double_t f_gaus_pol1_fit_par2 = f_gaus_pol1->GetParameter(2);
        Double_t f_gaus_pol1_fit_par2_err = f_gaus_pol1->GetParError(2);
        Double_t f_gaus_pol1_fit_par3 = f_gaus_pol1->GetParameter(3);
        Double_t f_gaus_pol1_fit_par3_err = f_gaus_pol1->GetParError(3);
        Double_t f_gaus_pol1_fit_par4 = f_gaus_pol1->GetParameter(4);
        Double_t f_gaus_pol1_fit_par4_err = f_gaus_pol1->GetParError(4);

        Double_t n_f_gaus_pol1_signal = TMath::Abs(f_gaus_pol1_fit_par0*f_gaus_pol1_fit_par2*TMath::Power(2.*TMath::Pi(),0.5));
        cout<<"n_f_gaus_pol1_signal: "<<n_f_gaus_pol1_signal<<endl;

	TF1* f_pol1 = new TF1("f_pol1","pol1(0)",h_energy_peak_min, h_energy_peak_max);
        f_pol1->SetParameters(f_gaus_pol1_fit_par3, f_gaus_pol1_fit_par4);
        Double_t n_f_gaus_pol1_bg = f_pol1->Integral(h_energy_peak_min,h_energy_peak_max);
        cout<<"n_f_gaus_pol1_bg: "<<n_f_gaus_pol1_bg<<endl;

        //fit function for the cut histogram
        TF1 * f_gaus_pol1_cut = new TF1("f_gaus_pol1_cut","gaus(0)+pol1(3)",h_energy_peak_before_min,h_energy_peak_after_max);
        Double_t f_gaus_pol1_cut_par[5] = {f_gaus_fit_par0*0.5,f_gaus_fit_par1, f_gaus_fit_par2,n_h_energy_peak_bg/(h_energy_peak_max-h_energy_peak_min), 0.};
        f_gaus_pol1_cut->SetParameters(f_gaus_pol1_cut_par);
        f_gaus_pol1_cut->SetLineColor(kOrange+8);
        f_gaus_pol1_cut->SetLineStyle(2);

        //histogram w/ cuts
        TH1D *h_energy_zoom_2_cut = (TH1D*)h_energy_zoom_cut->Clone("h_energy_zoom_2_cut");
        h_energy_zoom_2_cut->Scale(1./h_energy_xbin_size);
        h_energy_zoom_2_cut->SetXTitle("Energy [keV]");
        h_energy_zoom_2_cut->SetYTitle("Rate [cts/keV]");

        h_energy_zoom_2_cut->Fit(f_gaus_pol1_cut,"R");
        Double_t f_gaus_pol1_cut_getchisq = f_gaus_pol1_cut->GetChisquare();
        Double_t f_gaus_pol1_cut_getndf = f_gaus_pol1_cut->GetNDF();
        Double_t f_gaus_pol1_cut_reduced_chisq = f_gaus_pol1_cut_getchisq/f_gaus_pol1_cut_getndf;
        Double_t f_gaus_pol1_cut_pvalue = f_gaus_pol1_cut->GetProb();
        Double_t f_gaus_pol1_cut_fit_par0 = f_gaus_pol1_cut->GetParameter(0);
        Double_t f_gaus_pol1_cut_fit_par0_err = f_gaus_pol1_cut->GetParError(0);
        Double_t f_gaus_pol1_cut_fit_par1 = f_gaus_pol1_cut->GetParameter(1);
        Double_t f_gaus_pol1_cut_fit_par1_err = f_gaus_pol1_cut->GetParError(1);
        Double_t f_gaus_pol1_cut_fit_par2 = f_gaus_pol1_cut->GetParameter(2);
        Double_t f_gaus_pol1_cut_fit_par2_err = f_gaus_pol1_cut->GetParError(2);
        Double_t f_gaus_pol1_cut_fit_par3 = f_gaus_pol1_cut->GetParameter(3);
        Double_t f_gaus_pol1_cut_fit_par3_err = f_gaus_pol1_cut->GetParError(3);
        Double_t f_gaus_pol1_cut_fit_par4 = f_gaus_pol1_cut->GetParameter(4);
        Double_t f_gaus_pol1_cut_fit_par4_err = f_gaus_pol1_cut->GetParError(4);

        Double_t n_f_gaus_pol1_cut_signal = TMath::Abs(f_gaus_pol1_cut_fit_par0*f_gaus_pol1_cut_fit_par2*TMath::Power(2.*TMath::Pi(),0.5));
        cout<<"n_f_gaus_pol1_cut_signal: "<<n_f_gaus_pol1_cut_signal<<endl;

        TF1* f_pol1_cut = new TF1("f_pol1_cut","pol1(0)",h_energy_peak_min, h_energy_peak_max);
        f_pol1_cut->SetParameters(f_gaus_pol1_cut_fit_par3, f_gaus_pol1_cut_fit_par4);
        Double_t n_f_gaus_pol1_cut_bg = f_pol1_cut->Integral(h_energy_peak_min,h_energy_peak_max);
        cout<<"n_f_gaus_pol1_cut_bg: "<<n_f_gaus_pol1_cut_bg<<endl;

        Double_t f_f_gaus_pol1_signal_x_eff = n_f_gaus_pol1_cut_signal/n_f_gaus_pol1_signal;
        Double_t f_f_gaus_pol1_bg_x_eff = n_f_gaus_pol1_cut_bg/n_f_gaus_pol1_bg;
        Double_t f_f_gaus_pol1_bg_rej_eff = (n_f_gaus_pol1_bg-n_f_gaus_pol1_cut_bg)/n_f_gaus_pol1_bg;

        Double_t f_f_gaus_pol1_signal_x_eff_l = x_eff->ClopperPearson(n_f_gaus_pol1_signal,n_f_gaus_pol1_cut_signal,cl,0);
        Double_t f_f_gaus_pol1_signal_x_eff_h = x_eff->ClopperPearson(n_f_gaus_pol1_signal,n_f_gaus_pol1_cut_signal,cl,1);
        Double_t f_f_gaus_pol1_bg_x_eff_l = x_eff->ClopperPearson(n_f_gaus_pol1_bg,n_f_gaus_pol1_cut_bg,cl,0);
        Double_t f_f_gaus_pol1_bg_x_eff_h = x_eff->ClopperPearson(n_f_gaus_pol1_bg,n_f_gaus_pol1_cut_bg,cl,1);
        Double_t f_f_gaus_pol1_bg_rej_eff_l = x_eff->ClopperPearson(n_f_gaus_pol1_bg,
                                                                    n_f_gaus_pol1_bg-n_f_gaus_pol1_cut_bg,cl,0);
        Double_t f_f_gaus_pol1_bg_rej_eff_h = x_eff->ClopperPearson(n_f_gaus_pol1_bg,
                                                                    n_f_gaus_pol1_bg-n_f_gaus_pol1_cut_bg,cl,1);

        cout<<"f_f_gaus_pol1_bg_x_eff: "<<f_f_gaus_pol1_bg_x_eff<< "[ "
            << f_f_gaus_pol1_bg_x_eff_h <<", "
            << f_f_gaus_pol1_bg_x_eff_l <<" ]"
            <<endl;

        cout<<"f_f_gaus_pol1_signal_x_eff: "<<f_f_gaus_pol1_signal_x_eff<< " [ "
            << f_f_gaus_pol1_signal_x_eff_h  <<", "
            << f_f_gaus_pol1_signal_x_eff_l <<" ]"
            <<endl;

        //compute residuals
        Int_t npp = nbin_width_peak_before+ (h_energy_peak_max_bin - h_energy_peak_min_bin + 1) + nbin_width_peak_after;
        cout<<"npp: "<<npp<<endl;
        TGraph * g_residual_gaus_pol1 = new TGraph();   TGraph * g_residual_gaus_pol1_cut = new TGraph();
        g_residual_gaus_pol1->SetTitle("Residual_gaus_pol1");
        g_residual_gaus_pol1_cut->SetTitle("Residual_gaus_pol1_cut");
        Double_t get_f_gaus_pol1, get_f_gaus_pol1_cut;
        Double_t get_h_energy_zoom_2, get_h_energy_zoom_2_err;  Double_t get_h_energy_zoom_2_cut, get_h_energy_zoom_2_cut_err;
        Double_t g_residual_gaus_pol1_x , g_residual_gaus_pol1_y;       Double_t g_residual_gaus_pol1_cut_x , g_residual_gaus_pol1_cut_y;

        int iij = 0; int iik = 0;
        for (int iii =0; iii<npp;iii++){
                //w/o cuts
                g_residual_gaus_pol1_x = h_energy_peak_before_min+(iii+0.5)*h_energy_xbin_size;
                get_f_gaus_pol1 = f_gaus_pol1->Eval(g_residual_gaus_pol1_x);
                get_h_energy_zoom_2 = h_energy_zoom_2->GetBinContent(h_energy_peak_min_bin-nbin_width_peak_before+iii);
                get_h_energy_zoom_2_err = h_energy_zoom_2->GetBinError(h_energy_peak_min_bin-nbin_width_peak_before+iii);


                if(get_h_energy_zoom_2_err>0){
                        g_residual_gaus_pol1_y = (get_h_energy_zoom_2-get_f_gaus_pol1)/get_h_energy_zoom_2_err;
                        g_residual_gaus_pol1->SetPoint(iij, g_residual_gaus_pol1_x,g_residual_gaus_pol1_y);
                        //                      cout<<"iii :" <<iii;
                        //                      cout<<"iij :" <<iij;
                        //                      cout<<"g_residual_gaus_pol1_x: "<<g_residual_gaus_pol1_x<<endl;
                        //                      cout<<"get_f_gaus_pol1: "<<get_f_gaus_pol1 <<endl;
                        //                      cout<<"get_h_energy_zoom_2: "<<get_h_energy_zoom_2<<endl;
                        //                      cout<<"g_residual_gaus_pol1_y: "<<g_residual_gaus_pol1_y<<endl;
                        //                      cin.get();
                        iij++;
                }//if

                //w cuts
                g_residual_gaus_pol1_cut_x = h_energy_peak_before_min+(iii+0.5)*h_energy_xbin_size;
                get_f_gaus_pol1_cut = f_gaus_pol1_cut->Eval(g_residual_gaus_pol1_cut_x);
                get_h_energy_zoom_2_cut = h_energy_zoom_2_cut->GetBinContent(h_energy_peak_min_bin-nbin_width_peak_before+iii);
                get_h_energy_zoom_2_cut_err = h_energy_zoom_2_cut->GetBinError(h_energy_peak_min_bin-nbin_width_peak_before+iii);

                if(get_h_energy_zoom_2_cut_err>0){
                        g_residual_gaus_pol1_cut_y = (get_h_energy_zoom_2_cut-get_f_gaus_pol1_cut)/get_h_energy_zoom_2_cut_err;
                        g_residual_gaus_pol1_cut->SetPoint(iik, g_residual_gaus_pol1_cut_x,g_residual_gaus_pol1_cut_y);
                        //                      cout<<"iii :" <<iii;
                        //                      cout<<"iij :" <<iij;
                        //                      cout<<"g_residual_gaus_pol1_cut_x: "<<g_residual_gaus_pol1_cut_x<<endl;
                        //                      cout<<"get_f_gaus_pol1_cut: "<<get_f_gaus_pol1_cut <<endl;
                        //                      cout<<"get_h_energy_zoom_2: "<<get_h_energy_zoom_2<<endl;
                        //                      cout<<"g_residual_gaus_pol1_cut_y: "<<g_residual_gaus_pol1_cut_y<<endl;
                        iik++;
                }//if

        }//for graph
        cout<<"h_energy_peak_before_min: " << h_energy_peak_before_min <<" , h_energy_peak_after_max: " <<h_energy_peak_after_max << endl;
        g_residual_gaus_pol1->GetXaxis()->SetRangeUser(h_energy_peak_before_min,h_energy_peak_after_max);
        g_residual_gaus_pol1_cut->GetXaxis()->SetRangeUser(h_energy_peak_before_min,h_energy_peak_after_max);

        //Compute uncertainty
        //x_eff is obtained as a sum of 1) statistical uncertainty obtinaed by using ClopperPearson 1sigma C.I
        //and 2) difference of the value obtained from the fit as systematic uncertainty
        Double_t x_eff_uncertainty_stat = (n_h_energy_peak_signal_x_eff_h - n_h_energy_peak_signal_x_eff_l)*0.5;
        Double_t x_eff_uncertainty_sys = TMath::Abs((n_h_energy_peak_signal_x_eff-f_f_gaus_pol1_signal_x_eff)*0.5);

        //Comparison with the binomial distribution
        Double_t x_eff_uncertainty_stat_bi = TMath::Power(n_h_energy_cut_peak_signal*(1.-n_h_energy_peak_signal_x_eff),0.5)/n_h_energy_peak_signal;

        Double_t x_eff_uncertainty =
                TMath::Power(TMath::Power(x_eff_uncertainty_stat,2)+TMath::Power(x_eff_uncertainty_sys,2),0.5);

        //x_eff_bg uncertainty
        Double_t x_eff_bg_uncertainty_stat = (n_h_energy_peak_bg_x_eff_h - n_h_energy_peak_bg_x_eff_l)*0.5;
        Double_t x_eff_bg_uncertainty_sys = TMath::Abs((n_h_energy_peak_bg_x_eff-f_f_gaus_pol1_bg_x_eff)*0.5);

        //Comparison with the binomial distribution
        Double_t x_eff_bg_uncertainty_stat_bi = TMath::Power((n_h_energy_cut_peak_before+n_h_energy_cut_peak_after)*(1.-n_h_energy_peak_bg_x_eff),0.5)/(n_h_energy_peak_before+n_h_energy_peak_after);

        Double_t x_eff_bg_uncertainty =
                TMath::Power(TMath::Power(x_eff_bg_uncertainty_stat,2)+TMath::Power(x_eff_bg_uncertainty_sys,2),0.5);


	//open file to write the fit results
	//        FILE *fout = fopen(output_file_name.Data(), "a+");
	//write output
        fprintf(fout,"%s %s %d %s %d  %1.1f   %1.3f %1.3f %d %d %1.3f %1.3f %1.3f %1.3f  %d %d  %1.3f %1.3f %1.3f  %1.3f %1.3f %1.3f   %1.3f %1.3f %1.3f %1.3f   %1.3f %1.3f %1.3f   %1.3f %1.3f %1.3f    %1.3f %1.3f %1.3f %1.3f\n",
		ds_name.Data(), 
		draw_branch.Data(), draw_branch_index, peak_name.Data(), peak_index,
		//		ch_number,
		h_energy_binsize, 
		h_energy_fit_min, h_energy_fit_max,
		n_sigma_for_peak_count, n_sigma_for_bg_count,
		h_energy_peak_min, h_energy_peak_max, h_energy_peak_before_min, h_energy_peak_after_max,
		n_h_energy_cut_peak_signal, n_h_energy_cut_peak_bg,
		n_h_energy_peak_signal_x_eff, n_h_energy_peak_signal_x_eff_l, n_h_energy_peak_signal_x_eff_h,
		f_f_gaus_pol1_signal_x_eff, f_f_gaus_pol1_signal_x_eff_l, f_f_gaus_pol1_signal_x_eff_h,
		//		f_gaus_pol1_pvalue, f_gaus_pol1_cut_pvalue,
		x_eff_uncertainty_stat, x_eff_uncertainty_stat_bi, x_eff_uncertainty_sys, x_eff_uncertainty,
		n_h_energy_peak_bg_x_eff, n_h_energy_peak_bg_x_eff_l, n_h_energy_peak_bg_x_eff_h,
		f_f_gaus_pol1_bg_x_eff, f_f_gaus_pol1_bg_x_eff_l, f_f_gaus_pol1_bg_x_eff_h,
		x_eff_bg_uncertainty_stat, x_eff_bg_uncertainty_stat_bi, x_eff_bg_uncertainty_sys, x_eff_bg_uncertainty);
		//		f_gaus_fit_par0, f_gaus_fit_par0_err, 
		//		f_gaus_fit_par1, f_gaus_fit_par1_err, 
		//		f_gaus_fit_par2, f_gaus_fit_par2_err,
		//		f_gaus_reduced_chisq, f_gaus_prob);

	//close output file
        fclose(fout);


	/// Draw Canvas
        TCanvas * can = new TCanvas("can","",0,0,1200,800);
        can->Divide(3,2);
        can->SetGridy(0);
        gStyle->SetOptFit(111111111);

	can->cd(1);
        h_energy->SetStats(0);
        h_energy->Draw();
        h_energy_cut->SetLineColor(kBlue);
        h_energy_cut->SetMarkerColor(kBlue);
        h_energy_cut->Draw("same");

        can->cd(2);
        h_energy_zoom->SetStats(0);
        h_energy_zoom->SetMinimum(0);
        h_energy_zoom->Draw();

        Double_t f_gaus_ext_min = f_gaus_fit_par1-f_gaus_fit_par2*6.;
        Double_t f_gaus_ext_max = f_gaus_fit_par1+f_gaus_fit_par2*6.;
        cout<<"f_gaus_ext_min: "<<f_gaus_ext_min<< " f_gaus_ext_max: "<<f_gaus_ext_max<<endl;

        TF1 *f_gaus_ext = new TF1("f_gaus_ext", "gaus",f_gaus_ext_min,f_gaus_ext_max);
        Double_t f_gaus_ext_par[3] = {f_gaus_fit_par0, f_gaus_fit_par1, f_gaus_fit_par2};
        f_gaus_ext->SetParameters(f_gaus_ext_par);
        f_gaus_ext->SetLineStyle(2);
        f_gaus_ext->SetLineColor(kMagenta);
        f_gaus_ext->Draw("same");

        Double_t l_ymax = h_energy_zoom->GetMaximum()+h_energy_zoom->GetMinimum();
        TLine *l = new TLine();
        l->SetLineColor(kGreen);
        l->SetLineWidth(2);
        l->DrawLine(h_energy_peak_before_max,0,h_energy_peak_before_max,l_ymax);
        l->DrawLine(h_energy_peak_after_min,0,h_energy_peak_after_min,l_ymax);

        TLine *l_before = new TLine();
        l_before->SetLineColor(kGreen);
        l_before->SetLineWidth(2);
        l_before->SetLineStyle(2);
        l_before->DrawLine(h_energy_peak_before_min,0,h_energy_peak_before_min,l_ymax);

        TLine *l_after = new TLine();
        l_after->SetLineColor(kGreen);
        l_after->SetLineWidth(2);
        l_after->SetLineStyle(2);
        l_after->DrawLine(h_energy_peak_after_max,0,h_energy_peak_after_max,l_ymax);

        Double_t tex_x = 0.25;
        Double_t tex_y = 0.81;
        Double_t tex_size = 0.065;

        TLatex * tex_n_signal_x_eff = new TLatex(tex_x, tex_y, Form("%1.3f [%1.3f, %1.3f]", n_h_energy_peak_signal_x_eff, n_h_energy_peak_signal_x_eff_l, n_h_energy_peak_signal_x_eff_h));
        tex_n_signal_x_eff->SetNDC();
        tex_n_signal_x_eff->SetTextFont(132);
        tex_n_signal_x_eff->SetTextSize(tex_size);
        tex_n_signal_x_eff->Draw("same");

        tex_y = 0.72;
        TLatex * tex_n_bg_x_eff = new TLatex(tex_x, tex_y, Form("(BG rej: %1.3f [%1.3f, %1.3f])", n_h_energy_peak_bg_rej_eff, n_h_energy_peak_bg_rej_eff_l, n_h_energy_peak_bg_rej_eff_h));
        tex_n_bg_x_eff->SetNDC();
        tex_n_bg_x_eff->SetTextFont(132);
        tex_n_bg_x_eff->SetTextSize(tex_size);
        tex_n_bg_x_eff->Draw("same");

        can->cd(3);
        g_residual->SetMarkerStyle(20);
        g_residual->SetMarkerSize(0.9);
        g_residual->GetYaxis()->SetTitle("Residual");
        g_residual->GetXaxis()->SetTitle("Energy [keV]");
        g_residual->Draw("ap");

        can->cd(4);
        pEff->Draw("ap");

        can->cd(5);
        h_energy_zoom_cut->SetLineColor(kBlue);
        h_energy_zoom_cut->SetMarkerColor(kBlue);
        h_energy_zoom_cut->SetStats(0);
        h_energy_zoom_cut->SetMinimum(0);
        h_energy_zoom_cut->Draw();
        l->DrawLine(h_energy_peak_min,0,h_energy_peak_min,l_ymax);
        l->DrawLine(h_energy_peak_max,0,h_energy_peak_max,l_ymax);
        l_before->DrawLine(h_energy_peak_before_min,0,h_energy_peak_before_min,l_ymax);
        l_after->DrawLine(h_energy_peak_after_max,0,h_energy_peak_after_max,l_ymax);

        Double_t f_gaus_cut_ext_min = f_gaus_cut_fit_par1-f_gaus_cut_fit_par2*6.;
        Double_t f_gaus_cut_ext_max = f_gaus_cut_fit_par1+f_gaus_cut_fit_par2*6.;
        cout<<"f_gaus_cut_ext_min: "<<f_gaus_cut_ext_min<< " f_gaus_cut_ext_max: "<<f_gaus_cut_ext_max<<endl;

        TF1 *f_gaus_cut_ext = new TF1("f_gaus_cut_ext", "gaus",f_gaus_cut_ext_min,f_gaus_cut_ext_max);
        Double_t f_gaus_cut_ext_par[3] = {f_gaus_cut_fit_par0, f_gaus_cut_fit_par1, f_gaus_cut_fit_par2};
        f_gaus_cut_ext->SetParameters(f_gaus_cut_ext_par);
        f_gaus_cut_ext->SetLineStyle(2);
        f_gaus_cut_ext->SetLineColor(kOrange+8);
        f_gaus_cut_ext->Draw("same");

        can->cd(6);
        g_residual_cut->SetMarkerStyle(20);
        g_residual_cut->SetMarkerSize(0.9);
        g_residual_cut->SetMarkerColor(kBlue);
        g_residual_cut->GetYaxis()->SetTitle("Residual");
        g_residual_cut->GetXaxis()->SetTitle("Energy [keV]");
        g_residual_cut->Draw("ap");

        //Save
        TString can_fig_name = current_path +"fig/";  can_fig_name += process_type; can_fig_name += "/"; can_fig_name +=draw_branch; can_fig_name +="/"; 
        TString can_C_name = current_path +"C/"; can_C_name += process_type; can_C_name +="/"; can_C_name +=draw_branch; can_C_name +="/"; 

	can_fig_name += peak_name; 	can_fig_name += "_count_";
        can_fig_name += h_energy_xbin_size;     can_fig_name += "keVbin_";
	can_fig_name += 5;       can_fig_name +="sigma_peak_range_";
        can_fig_name += list_name;//list_name comes from read_files.C
	//        can_fig_name +="_ch";  can_fig_name += 2;
        can_fig_name +="_m";  can_fig_name += 1;
        can_fig_name +=".png";

	can_C_name += peak_name; 	can_C_name += "_count_";
        can_C_name += h_energy_xbin_size;       can_C_name += "keVbin_";
	can_C_name += 5; can_C_name +="sigma_peak_range_";
        can_C_name += list_name;//list_name comes from read_files.C
	//        can_C_name +="_ch";  can_C_name += 2;
        can_C_name +="_m";  can_C_name += 1;
        can_C_name +=".C";

        can->SaveAs(can_fig_name.Data());
        can->SaveAs(can_C_name.Data());


        TCanvas *can2 = new TCanvas("can2","",1200,600);
        can2->Divide(2,2);
        can2->SetGridy(0);
        gStyle->SetOptFit(11111111);

        can2->cd(1);
        h_energy_zoom_2->SetStats(0);
        h_energy_zoom_2->SetMinimum(0);
        h_energy_zoom_2->Draw();

        Double_t l_ymax = h_energy_zoom_2->GetMaximum()+h_energy_zoom_2->GetMinimum();
        l->DrawLine(h_energy_peak_before_max,0,h_energy_peak_before_max,l_ymax);
        l->DrawLine(h_energy_peak_after_min,0,h_energy_peak_after_min,l_ymax);
        l_before->DrawLine(h_energy_peak_before_min,0,h_energy_peak_before_min,l_ymax);
        l_after->DrawLine(h_energy_peak_after_max,0,h_energy_peak_after_max,l_ymax);

        tex_x = 0.52;
        tex_y = 0.81;
        TLatex * tex_f_signal_x_eff = new TLatex(tex_x,tex_y, Form("%1.3f [%1.3f, %1.3f]",f_f_gaus_pol1_signal_x_eff,f_f_gaus_pol1_signal_x_eff_l, f_f_gaus_pol1_signal_x_eff_h));
        tex_f_signal_x_eff->SetNDC();
        tex_f_signal_x_eff->SetTextFont(132);
        tex_f_signal_x_eff->SetTextSize(tex_size);
        tex_f_signal_x_eff->Draw("same");

        tex_y = 0.72;
        TLatex * tex_f_bg_x_eff = new TLatex(tex_x,tex_y, Form("(BG rej: %1.3f [%1.3f, %1.3f])",f_f_gaus_pol1_bg_rej_eff,f_f_gaus_pol1_bg_rej_eff_l, f_f_gaus_pol1_bg_rej_eff_h));
        tex_f_bg_x_eff->SetNDC();
        tex_f_bg_x_eff->SetTextFont(132);
        tex_f_bg_x_eff->SetTextSize(tex_size);
        tex_f_bg_x_eff->Draw("same");

        can2->cd(2);
        g_residual_gaus_pol1->SetMarkerStyle(20);
        g_residual_gaus_pol1->SetMarkerSize(0.9);
        g_residual_gaus_pol1->GetXaxis()->SetTitle("Energy [keV]");
        g_residual_gaus_pol1->GetYaxis()->SetTitle("Residual");
        g_residual_gaus_pol1->Draw("ap");

	can2->cd(3);
        h_energy_zoom_2_cut->SetLineColor(kBlue);
        h_energy_zoom_2_cut->SetMarkerColor(kBlue);
        h_energy_zoom_2_cut->SetStats(0);
        h_energy_zoom_2_cut->SetMinimum(0);
        h_energy_zoom_2_cut->Draw("");

        l->DrawLine(h_energy_peak_before_max,0,h_energy_peak_before_max,l_ymax);
        l->DrawLine(h_energy_peak_after_min,0,h_energy_peak_after_min,l_ymax);
        l_before->DrawLine(h_energy_peak_before_min,0,h_energy_peak_before_min,l_ymax);
        l_after->DrawLine(h_energy_peak_after_max,0,h_energy_peak_after_max,l_ymax);

        can2->cd(4);
        g_residual_gaus_pol1_cut->SetMarkerColor(kBlue);
        g_residual_gaus_pol1_cut->SetMarkerStyle(20);
        g_residual_gaus_pol1_cut->SetMarkerSize(0.9);
        g_residual_gaus_pol1_cut->GetXaxis()->SetTitle("Energy [keV]");
        g_residual_gaus_pol1_cut->GetYaxis()->SetTitle("Residual");
        g_residual_gaus_pol1_cut->Draw("ap");

        //Save
        TString can2_fig_name = current_path +"fig/"; can2_fig_name += process_type; can2_fig_name += "/"; can2_fig_name +=draw_branch; can2_fig_name +="/"; 
        TString can2_C_name = current_path +"C/"; can2_C_name += process_type; can2_C_name += "/"; can2_C_name +=draw_branch; can2_C_name +="/"; 

        can2_fig_name += peak_name;
        can2_fig_name +="_fit_";
        can2_fig_name += h_energy_xbin_size;    can2_fig_name += "keVbin_";
        can2_fig_name += 5;      can2_fig_name +="sigma_peak_range_";
        can2_fig_name += list_name;
        can2_fig_name +="_m";  can2_fig_name += 1;
        can2_fig_name +=".png";

        can2_C_name += peak_name;
        can2_C_name +="_fit_";
        can2_C_name += h_energy_xbin_size;      can2_C_name += "keVbin_";
        can2_C_name += 5;        can2_C_name +="sigma_peak_range_";
        can2_C_name += list_name;
        can2_C_name +="_m";  can2_C_name += 1;
        can2_C_name +=".C";

        can2->SaveAs(can2_fig_name.Data());
        can2->SaveAs(can2_C_name.Data());

}
