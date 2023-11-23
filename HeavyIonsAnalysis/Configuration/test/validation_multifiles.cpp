// This is called after RhoConeAnaCondor for Rho and random cone variables
// RhoConeAnaCondor requires the Rho and random cone (0.2 and 0.4) trees

void centrality_plot(std::string &pdf_file, std::vector<TFile *> in_files, std::vector<std::string> processes, std::vector<Int_t> colours){
	TCanvas *c = new TCanvas();
	c->cd();
	TPad *pad1 = new TPad(("tpad" + std::to_string(1)).c_str(), ("tpad" + std::to_string(1)).c_str(), 0, 0.29, 1, 1);
	pad1->SetBottomMargin(0); // Upper and lower plot are joined
	pad1->SetLogy();
	// dtpad->SetLogy();
   	pad1->Draw();        // Go back to the main canvas before defining pad2
   	TPad *pad2 = new TPad(("bpad" + std::to_string(1)).c_str(), ("bpad" + std::to_string(1)).c_str(), 0, 0, 1, 0.3);
   	pad2->SetTopMargin(0);
   	pad2->SetBottomMargin(0.3);
   	pad2->SetGridy();
   	pad2->Draw();
	// c->SetLogy();
	TLegend *legend;
	legend = new TLegend( 0.7, 0.7, 0.9, 0.9 );
	for(Int_t i{0}; i < in_files.size(); ++i){
		TH1 *h = (TH1*)in_files[i]->Get( "centrality" );
		TH1F *a = (TH1F*) h->Clone();
	   	a->SetDirectory(0);
		a->SetMarkerColor(colours.at(i));
		legend->AddEntry( a, processes.at(i).c_str(), "p" );
		a->SetMarkerStyle(20);
		a->GetXaxis()->SetRange(0,180);
		a->Scale(1./a->Integral(0,180));
		a->SetMaximum(1.);
		a->SetMinimum(0.00001);
		if(i == 1){
			
		}
		pad1->cd();
		a->Draw("SAME,e0e2x0p");
	}
	pad1->cd();
	legend->Draw();
	// c->Update();
	pad2->cd();
	for(Int_t i{1};i<in_files.size();++i){
		TH1F *rat = (TH1F*) in_files[0]->Get("centrality")->Clone();
		rat->SetDirectory(0);
		rat->GetXaxis()->SetRange(0,180);
		rat->Scale(1./rat->Integral(0,180));
		TH1F *other_hist = (TH1F*) in_files[i]->Get("centrality")->Clone();
		other_hist->SetDirectory(0);
		other_hist->GetXaxis()->SetRange(0,180);
		other_hist->Scale(1./other_hist->Integral(0,180));
		rat ->Divide(other_hist,rat);
		rat->SetTitle("");
		rat->SetMarkerStyle(20);
		rat->SetMarkerColor(colours[i]);
		rat->SetLineColor(colours[i]);
		rat->GetYaxis()->SetTitle( ( " X / " + processes.at(0) ).c_str());
	   	rat->GetYaxis()->CenterTitle();
	   	rat->GetYaxis()->SetTitleSize(13);
	   	rat->GetYaxis()->SetTitleFont(43);
	   	rat->GetYaxis()->SetTitleOffset(1.5);
		rat->Draw("SAME,e0x0p");
		rat->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
	   	rat->GetYaxis()->SetLabelSize(12);
	   	rat->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
	   	rat->GetXaxis()->SetLabelSize(12);
	   	rat->GetXaxis()->SetTitleSize(13);
	   	rat->GetXaxis()->SetTitleFont(43);
	   	rat->GetXaxis()->SetTitleOffset(2.5);
	   	rat->SetMaximum(1.5);
	   	rat->SetMinimum(0.5);
	   	rat->GetYaxis()->SetNdivisions(5, kTRUE);
   }
   c->Print( ( pdf_file ).c_str(),"pdf");
}

void separate_fits(std::string &pdf_file, std::string fits_file, std::vector<TFile *> in_files, std::string radius, std::vector<std::string> processes, std::vector<Int_t> colours ){
	std::vector < Int_t > band_lowbin = {0, 20, 40, 60, 80, 100, 120, 140};
    std::vector < Int_t > band_highbin = {20, 40, 60, 80, 100, 120, 140, 160};
    TCanvas *c = new TCanvas();
	c->SetLogy();
	std::cout << "before canvases" << std::endl;
	std::vector<TCanvas*> vc = {};
	std::vector<TPad*> tp = {};
	std::vector<TPad*> bp = {};
	for( int i=0; i<4; ++i){
		TCanvas *dummyC = new TCanvas();
		dummyC->cd();
		TPad *dtpad = new TPad(("tpad" + std::to_string(i)).c_str(), ("tpad" + std::to_string(i)).c_str(), 0, 0.29, 1, 1);
   		dtpad->SetBottomMargin(0); // Upper and lower plot are joined
   		if( i == 1 || i == 2 )
	   		dtpad->SetLogy();	
   		tp.push_back(dtpad);
	   	dtpad->Draw();             // Draw the upper pad: pad1
		dummyC->cd();          // Go back to the main canvas before defining pad2
	   	TPad *dbpad = new TPad(("bpad" + std::to_string(i)).c_str(), ("bpad" + std::to_string(i)).c_str(), 0, 0, 1, 0.3);
	   	dbpad->SetTopMargin(0);
	   	dbpad->SetBottomMargin(0.3);
	   	dbpad->SetGridy();
	   	dbpad->Draw();
	   	// dbpad->SetLogy();
	   	bp.push_back(dbpad);
		vc.push_back(dummyC);
	}
	// TCanvas *c1 = new TCanvas();
	// c1->SetLogy();
	// TCanvas *c2 = new TCanvas();
	// TCanvas *c3 = new TCanvas();
	// c3->SetLogy();
	// TCanvas *c4 = new TCanvas();
	TLegend * legend;
	legend = new TLegend( 0.7, 0.7, 0.9, 0.9 );

    //loop over files
    std::vector <TH1F*> to_save_shift;
	std::vector <TH1F*> to_save_sigma;
	std::vector <TH1F*> to_save_dist_shift;
	std::vector <TH1F*> to_save_dist_sigma;
	//book histograms
	Int_t fit_iter = 7;	
    for(Int_t i{0}; i < in_files.size(); ++i ){
    	//histograms to contain the fit parameters as fns of centrality
    	TH1F *a = new TH1F( (radius + std::to_string(i) + processes.at(i) + "shift").c_str(),( radius + "  #mu_{LHS};HiBin;#mu_{LHS} [GeV]").c_str(), 8, 0, 160);
		TH1F *b = new TH1F( (radius + std::to_string(i) + processes.at(i) + "sigma").c_str(),( radius + "  #sigma_{LHS};HiBin;#sigma_{LHS} [GeV]").c_str(), 8, 0, 160);
		TH1F *a_dist = new TH1F( (radius + std::to_string(i) + processes.at(i) + "distribution shift").c_str(),( radius + "  #mu_{total};HiBin;#mu_{total} [GeV]").c_str(), 8, 0, 160);
		TH1F *b_dist = new TH1F( (radius + std::to_string(i) + processes.at(i) + "distribution sigma").c_str(),( radius + "  #sigma_{total};HiBin;#sigma_{total} [GeV]").c_str(), 8, 0, 160);
		a->SetMaximum(5);
		a->SetMinimum(-5);
		a_dist->SetMaximum(5);
		a_dist->SetMinimum(-5);
		b->SetMaximum(50);
		b->SetMinimum(0.1);
		b_dist->SetMaximum(50);
		b_dist->SetMinimum(0.1);
		to_save_shift.push_back(a);
		to_save_sigma.push_back(b);
		to_save_dist_shift.push_back(a_dist);
		to_save_dist_sigma.push_back(b_dist);
    	//loop over histograms
    	Double_t a1_sigma = 10;
		Double_t a1_mean = 0;
    	for(size_t j{0}; j<band_lowbin.size();++j){
    		std::string histogram = radius + " axis corr. pt hiBin [" + std::to_string(band_lowbin.at(j)) + "-" + std::to_string(band_highbin.at(j)) + "]";
    		TH1 *af = (TH1*)in_files.at(i)->Get( histogram.c_str() );
    		std::string h_title = af->GetTitle();
    		af->SetTitle((processes.at(i) + " " + h_title + ";#Delta p_{T};").c_str());
			Double_t final_mean = 0;
			Double_t final_sigma = 0;
			Double_t final_mean_err = 0;
			Double_t final_sigma_err = 0;
			TF1 *g1 = new TF1("g1","gaus",-100,100);
			g1->SetParameters(1,a1_mean,a1_sigma);
			c->cd();
			af->Draw();
			std::cout <<"Fit of " << af ->GetTitle() << std::endl;
			for(Int_t i{1}; i<=fit_iter;++i){
				af->Fit("g1","Q","",g1->GetParameter(1)-2*g1->GetParameter(2), g1->GetParameter(1)+g1->GetParameter(2));
				TF1* fit = af->GetFunction("g1");
				TF1* fit_copy = (TF1*)fit->IsA()->New();
				fit->Copy(*fit_copy);
				fit_copy->SetRange( g1->GetParameter(1)-2*g1->GetParameter(2), g1->GetParameter(1)+g1->GetParameter(2) );
				fit_copy->SetLineColor(kRed);
				fit_copy->SetFillColorAlpha(kRed, 0.15);
				fit_copy->SetLineWidth(1.5);
				fit_copy->SetLineStyle(4);
				// fit_copy->SetFillStyle(3805+s*20);
				fit_copy->SetFillStyle(1001);
				// fit->SetLineColor(3);
				if(i == fit_iter){
					fit_copy->Draw("SAME FC");
					final_mean = g1->GetParameter(1);
					final_mean_err = g1->GetParError(1);
					final_sigma = g1->GetParameter(2);
					final_sigma_err = g1->GetParError(2);
					TF1* fit_copy2 = (TF1*)fit->IsA()->New();
					fit->Copy(*fit_copy2);
					fit_copy2->SetRange(-100,100);
					fit_copy2->SetLineColor(kBlue);
					fit_copy2->SetLineWidth(1);
					fit_copy2->Draw("SAME FC");
					a1_mean = g1->GetParameter(1);
					a1_sigma = g1->GetParameter(2);
					TLatex l;
			   		l.SetTextSize(0.03);
			      	l.DrawLatexNDC(0.55,0.85,Form("#mu_{LHS}=%.3g, #sigma_{LHS}=%.3g, #sigma_{hist}=%.3g",final_mean,final_sigma, a->GetStdDev()));
				}
			}
			c->Update();
			c->Print( ( fits_file ).c_str(),"pdf");
			c->Clear();
			std::cout << "after fits" << std::endl;
			// std::cout << projY->GetMean() << " at bin " << j << " [" << band_lowbin.at(j) << "-" << band_highbin.at(j) << "]" << std::endl;
			to_save_shift.at(i)->SetBinContent(j+1, final_mean);
			to_save_shift.at(i)->SetBinError(j+1, final_mean_err);
			to_save_sigma.at(i)->SetBinContent(j+1, final_sigma);
			to_save_sigma.at(i)->SetBinError(j+1, final_sigma_err);
			to_save_dist_shift.at(i)->SetBinContent(j+1, af->GetMean());
			to_save_dist_sigma.at(i)->SetBinContent(j+1, af->GetStdDev());
    	}
    	to_save_shift.at(i)->SetMarkerColor(colours.at(i));
		to_save_shift.at(i)->SetFillColor(colours.at(i));
		to_save_sigma.at(i)->SetMarkerColor(colours.at(i));
		to_save_sigma.at(i)->SetFillColor(colours.at(i));
		to_save_dist_shift.at(i)->SetMarkerColor(colours.at(i));
		to_save_dist_shift.at(i)->SetFillColor(colours.at(i));
		to_save_dist_sigma.at(i)->SetMarkerColor(colours.at(i));
		to_save_dist_sigma.at(i)->SetFillColor(colours.at(i));
		legend->AddEntry( to_save_shift.at(i), processes.at(i).c_str(), "p" );
			to_save_shift.at(i)->SetMarkerStyle(20);
			to_save_sigma.at(i)->SetMarkerStyle(20);
			to_save_dist_shift.at(i)->SetMarkerStyle(20);
			to_save_dist_sigma.at(i)->SetMarkerStyle(20);
		
		// else{
		// 	to_save_shift.at(i)->SetMarkerStyle(20);
		// 	to_save_sigma.at(i)->SetMarkerStyle(20);
		// 	to_save_dist_shift.at(i)->SetMarkerStyle(20);
		// 	to_save_dist_sigma.at(i)->SetMarkerStyle(20);
		// }
		// vc.at(1)->cd();
		tp.at(0)->cd();
		tp.at(0)->Update();
		// c2->cd();
		to_save_shift.at(i)->Draw("SAME,e0e2x0p");
		// tp.at(0)->Update();
		// c1->cd();
		vc.at(1)->cd();
		tp.at(1)->cd();
		to_save_sigma.at(i)->Draw("SAME,e0e2x0p");
		// c3->cd();
		vc.at(2)->cd();
		tp.at(2)->cd();
		to_save_dist_sigma.at(i)->Draw("SAME,e0e2x0p");
		// c4->cd();
		vc.at(3)->cd();
		tp.at(3)->cd();
		to_save_dist_shift.at(i)->Draw("SAME,e0e2x0p");
		std::cout << "end" << std::endl;
    }
    std::vector<std::vector <TH1F*> > four_histos;
    four_histos.push_back(to_save_shift);
    four_histos.push_back(to_save_sigma);
    four_histos.push_back(to_save_dist_sigma);
    four_histos.push_back(to_save_dist_shift);
    for( int i{0};i<4;++i){
	    vc.at(i)->cd();
	    tp.at(i)->Update();
	    tp.at(i)->cd();
	    tp.at(i)->Draw();
		legend->Draw();
		vc.at(i)->Update();
		bp.at(i)->cd();
		for(Size_t f{1};f<in_files.size();++f){
			TH1F *rat = (TH1F*) four_histos.at(i).at(0)->Clone();
			rat->SetDirectory(0);
			rat ->Divide(four_histos.at(i).at(f),rat);
			rat->SetTitle("");
			rat->SetMarkerStyle(20);
			rat->SetMarkerColor(colours.at(f));
			rat->SetLineColor(colours.at(f));
			rat->GetYaxis()->SetTitle( ( " X / " + processes.at(0) ).c_str());
		   	rat->GetYaxis()->CenterTitle();
		   	rat->GetYaxis()->SetTitleSize(13);
		   	rat->GetYaxis()->SetTitleFont(43);
		   	rat->GetYaxis()->SetTitleOffset(1.5);
			rat->Draw("SAME,e0x0p");
			rat->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
		   	rat->GetYaxis()->SetLabelSize(12);
		   	rat->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
		   	rat->GetXaxis()->SetLabelSize(12);
		   	rat->GetXaxis()->SetTitleSize(13);
		   	rat->GetXaxis()->SetTitleFont(43);
		   	rat->GetXaxis()->SetTitleOffset(2.5);
		   	rat->SetMaximum(1.5);
		   	rat->SetMinimum(0.5);
		   	if( radius.find("R=0.2")!=std::string::npos && i==0 ){
		   		rat->SetMaximum(5);
		   		rat->SetMinimum(0.3);
		   	}
		   	rat->GetYaxis()->SetNdivisions(5, kTRUE);
		}
		vc.at(i)->Update();
		vc.at(i)->Print( ( pdf_file ).c_str(),"pdf");
	   	vc.at(i)->Clear();
	}
}

void scatter_proj( std::string &pdf_file, std::vector<TFile *> in_files, std::string histogram, std::vector<std::string> processes, std::vector<Int_t> colours, Bool_t fitGaussian ){
	std::vector < Int_t > band_lowbin = {1, 21, 41, 61, 81, 101, 121, 141};
	std::vector < Int_t > band_highbin = {20, 40, 60, 80, 100, 120, 140, 160};
	TH2 *test = (TH2*)in_files.at(0)->Get(histogram.c_str());
	Double_t histXMin = test->GetXaxis()->GetXmin();
	Double_t histYMax = test->GetXaxis()->GetXmax();
	// std::cout << histXMin << " " << histYMax << " ranges" << std::endl;
	TCanvas *c = new TCanvas();
	TPad *pad1 = new TPad(("tpad" + std::to_string(1)).c_str(), ("tpad" + std::to_string(1)).c_str(), 0, 0.29, 1, 1);
	pad1->SetBottomMargin(0); // Upper and lower plot are joined
	// dtpad->SetLogy();
   	pad1->Draw();        // Go back to the main canvas before defining pad2
   	TPad *pad2 = new TPad(("bpad" + std::to_string(1)).c_str(), ("bpad" + std::to_string(1)).c_str(), 0, 0, 1, 0.3);
   	pad2->SetTopMargin(0);
   	pad2->SetBottomMargin(0.3);
   	pad2->SetGridy();
   	pad2->Draw();
	c->SetLogy();
	TLegend * legend;
	legend = new TLegend( 0.7, 0.7, 0.9, 0.9 );
	std::vector <TH1F*> to_save_shift;
	TH2 *testing = (TH2*)in_files[0]->Get( histogram.c_str() );
	std::string axis_name = "";
	std::string hist_title(testing->GetTitle());
	if((hist_title.find( "R=" )) != std::string::npos)
		axis_name = "<p_{T}> [GeV]";
	if((hist_title.find( "#rho" )) != std::string::npos)
		axis_name = "<#rho> [GeV]";
	for(Int_t i{0}; i < in_files.size(); ++i){
		TH1F *a = new TH1F( ("tosave " + hist_title + processes.at(i)).c_str(), (hist_title + ";hiBin;" + axis_name).c_str(), 8, 0, 160);
		a->SetMaximum(250);
		to_save_shift.push_back(a);
	}
	// TH1F *to_save_shift = new TH1F(("mean " + histogram).c_str(), ("mean " + histogram).c_str());
	for(size_t i{0}; i < in_files.size(); ++i){
		TH2 *a = (TH2*)in_files[i]->Get( histogram.c_str() );
		// TH1F *to_save_shift = new TH1F(("mean " + histogram).c_str(), ("mean " + histogram).c_str(), 10, histXMin, histYMax);
		for(size_t j{0}; j < band_lowbin.size(); ++j){
			TCanvas *c1 = new TCanvas();
			TH1D* projY;
			projY = a->ProjectionY("",band_lowbin.at(j),band_highbin.at(j));
			projY->SetTitle((processes.at(i) + histogram + " [" + std::to_string(band_lowbin.at(j)) + "-" + std::to_string(band_highbin.at(j)) + "]").c_str());
			c1->cd();
			c1->SetLogy();
			// projY->Draw();
			// c1->Print( ( pdf_file ).c_str(),"pdf");
			c1->Clear();
			// std::cout << projY->GetMean() << " at bin " << j << " [" << band_lowbin.at(j) << "-" << band_highbin.at(j) << "]" << std::endl;
			to_save_shift.at(i)->SetBinContent(j+1,projY->GetMean());
			to_save_shift.at(i)->SetBinError(j+1,projY->GetStdDev());
		}
		c->cd();
		pad1->cd();
		pad1->SetLogy();
		to_save_shift.at(i)->SetMarkerColor(colours.at(i));
		to_save_shift.at(i)->SetLineColor(colours.at(i));
		legend->AddEntry( to_save_shift.at(i), processes.at(i).c_str(), "p" );
		to_save_shift.at(i)->SetMarkerStyle(20);
		to_save_shift.at(i)->Draw("SAME,p");
		c->Update();
	}
	pad1->cd();
	legend->Draw();
	pad2->cd();
	//ratio plots on pad2
	for(Size_t i{1}; i< in_files.size();++i){
		TH1F *rat = (TH1F*) to_save_shift.at(0)->Clone();
		rat->SetDirectory(0);
		rat ->Divide(to_save_shift.at(i),rat);
		rat->SetTitle("");
		rat->SetMarkerStyle(20);
		rat->SetMarkerColor(colours[i]);
		rat->SetLineColor(colours[i]);
		rat->SetFillColor(colours[i]);
		rat->GetYaxis()->SetTitle( ( "X / " + processes.at(0) ).c_str());
	   	rat->GetYaxis()->CenterTitle();
	   	rat->GetYaxis()->SetTitleSize(13);
	   	rat->GetYaxis()->SetTitleFont(43);
	   	rat->GetYaxis()->SetTitleOffset(1.5);
		rat->Draw("SAME,e0x0p");
		rat->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
	   	rat->GetYaxis()->SetLabelSize(12);
	   	rat->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
	   	rat->GetXaxis()->SetLabelSize(12);
	   	rat->GetXaxis()->SetTitleSize(13);
	   	rat->GetXaxis()->SetTitleFont(43);
	   	rat->GetXaxis()->SetTitleOffset(2.5);
	   	rat->SetMaximum(1.5);
	   	rat->SetMinimum(0.5);
	   	rat->GetYaxis()->SetNdivisions(5, kTRUE);
		
	}
	c->Print( ( pdf_file ).c_str(),"pdf");
	c->Clear();
}

void validation_multifiles(){
	TH1::SetDefaultSumw2( );
	//name and destinations of output files
	std::string nameOfFile = "HIEHIERP_comparison.pdf";
	std::string pdf_file = "/afs/cern.ch/user/v/vavladim/public/pdfFolder/" + nameOfFile;
	std::string fits_file = "/afs/cern.ch/user/v/vavladim/public/pdfFolder/fits_" + nameOfFile;
	Int_t ci = 2555;
	// colours and names of the histograms
	// new TColor(++ci, 54/256., 161/256., 24/256.); // dark green
	new TColor(++ci, 54/256., 161/256., 24/256.); // dark green
	new TColor(++ci, 252/256., 161/256., 3/256.);	// dark orange
	new TColor(++ci, 173/256., 16/256., 204/256.); // purple
	new TColor(++ci, 242/256., 235/256., 36/256.); // yellow
	new TColor(++ci, 36/256., 116/256., 156/256.); // deep blue
	new TColor(++ci, 132/256., 240/256., 209/256.); // cyan
	new TColor(++ci, 252/256., 3/256., 190/256.); //pink
	new TColor(++ci, 121/256., 133/256., 78/256.); // olive

	ci = 2555;
	std::vector < Int_t > colours = { ci+1, ci+2, ci+3, ci+4, ci+5, ci+6, ci+7, ci+8 };
	// std::vector <std::string> processes = {"Run3 data", "HYDIJET", "PbPb data", "Pythia embedded", "pp data", "pp Pythia", "PbPb CaloTrig" };
	std::vector<TFile *> inf = {};
	//location of directory containing the combined files for each run
	std::string dir = "/afs/cern.ch/user/v/vavladim/public/rootFolder/";
	std::vector<std::string> filenames = {	
											// "Mbias_combined.root",
											// "trig_filt_new_data.root", 
											// "half_second_run_processed.root", 
											// "374345.root",
											// "374354.root",
											"HIExpress_374322.root",
											"HIExpressRawPrime_374322.root",
											"HIExpress_374345.root",
											"HIExpressRawPrime_374345.root",
											"HIExpress_374354.root",
											"HIExpressRawPrime_374354.root",
											"HIExpress_374289.root",
											"HIExpressRawPrime_374289.root",

										};
	for(Size_t s{0}; s<filenames.size();++s){
		TFile *f = new TFile((dir + filenames.at(s)).c_str(),"READ");
		inf.push_back(f);
	}
	//The names corresponding to the files (could use the filenames and cut the names)
	std::vector<std::string> processes = { "374322 HIE", "374322 HIERP", "374345 HIE", "374345 HIERP", "374354 HIE", "374354 HIERP","374289 HIE", "374289 HIERP"  };
	// std::vector<std::string> processes ={"Run 2", "HIExpressRawPrime 374354"};
	gROOT->SetBatch(kTRUE);
	TIter gonext( inf.at(0)->GetListOfKeys() );
	TKey *key;
	TCanvas *c = new TCanvas();
	c->cd();
	c->Print( ( pdf_file + "(").c_str(),"pdf");
	c->Print( ( fits_file + "(").c_str(),"pdf");
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0, 1.0, 1.0);
	pad1->cd();
	gStyle->SetOptStat(0);
	scatter_proj(pdf_file, inf, "rho vs centrality acceptance", processes, colours, false );
	separate_fits(pdf_file, fits_file, inf, "R=0.4", processes, colours);
	separate_fits(pdf_file, fits_file, inf, "R=0.2", processes, colours);
	scatter_proj(pdf_file, inf, "raw R=0.4 cone pt vs centrality", processes, colours, false );
	scatter_proj(pdf_file, inf, "raw R=0.2 cone pt vs centrality", processes, colours, false );
	centrality_plot(pdf_file, inf, processes, colours);
	c->Clear();
	c->Print( ( pdf_file+ ")" ).c_str(),"pdf");
	c->Print( ( fits_file+ ")" ).c_str(),"pdf");
}