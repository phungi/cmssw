#include "TCanvas.h"
#include <iterator>
#include "TMath.h"
#include <cmath>
#include <cstdio>
#include <ctime>
#include "THStack.h"
#include "TFile.h"
#include "TLegend.h"
#include <typeinfo>
#include "TAxis.h"
// #include <stdlib.h>

Int_t histogram_number(std::vector <Int_t> &band_lowbin, std::vector <Int_t> &band_highbin, Int_t hiBin){
    for(Int_t i{0}; i<band_lowbin.size(); ++i){
        if(hiBin >= band_lowbin.at(i) && hiBin <= band_highbin.at(i)){
            return i;
        }
    }
    std::cout << "can't find where to put this " << hiBin << std::endl;
    return 0;
}

void Reduced_rho(std::string input_file, std::string output_file){
    Bool_t isMC = false;
    Bool_t ispp = false;
    if( ( input_file.find( "embedded" ) ) != std::string::npos || ( input_file.find( "pythia" ) ) != std::string::npos || ( input_file.find( "HYDJET" ) ) != std::string::npos ){
        isMC = true;
    }
    if( ( input_file.find( "pp" ) ) != std::string::npos ){
        ispp = true;
    }

    std::cout << "Starting the processing..." << std::endl;
	const Int_t MAXJETS=100;
	TH1::SetDefaultSumw2( );
	std::string filename;
    std::cout << "Opening file..." << std::endl;
	TFile *file = new TFile(input_file.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cout << "ERROR opening file" << std::endl;
        return;
    }
	TTree *t; //tree
	TTree *w_t; //weight tree
    TTree *hlt_t;
    TTree *add_sel_t;
    TTree *randomCone4_t;
    TTree *randomCone2_t;
    TTree *rho_t;
	std::string path = "akCs4PFJetAnalyzer/t";
    if(ispp)
        path = "ak4PFJetAnalyzer/t";
	std::string weight_path = "hiEvtAnalyzer/HiTree";
    std::string trigger_path = "hltanalysis/HltTree";
    std::string add_sel_tree_path = "skimanalysis/HltTree";
    std::string rho_path = "rhoAnalysis/t";
    std::string cone2_path = "randomConeAnalysisR2/t";
    std::string cone4_path = "randomConeAnalysisR4/t";
    Bool_t hasRhoAnaR4 = false;
    Bool_t hasRhoAnaR2 = false;
    std::cout << "Getting the trees..." << std::endl;
	file->GetObject(path.c_str(),t);
	file->GetObject(weight_path.c_str(),w_t);
    file->GetObject(rho_path.c_str(), rho_t);
    if( ( input_file.find( "Run3_data" ) ) != std::string::npos || ( input_file.find( "Run2_data" ) ) != std::string::npos || ( input_file.find( "PbPb" ) ) != std::string::npos || ( input_file.find( "embedded" ) ) != std::string::npos || ( input_file.find( "HYDIJET" ) ) != std::string::npos || ( input_file.find( "rho" ) ) != std::string::npos )
        std::cout << "Opening the tree" << std::endl;
        file->GetObject( add_sel_tree_path.c_str(), add_sel_t );
	if( !isMC )
    	file->GetObject( trigger_path.c_str(), hlt_t );
    TDirectory *dir = file->GetDirectory("randomConeAnalysisR4");
    if (dir){
        file->GetObject( cone4_path.c_str(), randomCone4_t );
        hasRhoAnaR4 = true;
    }
    TDirectory *dir1 = file->GetDirectory("randomConeAnalysisR2");
    if (dir1){
        file->GetObject( cone2_path.c_str(), randomCone2_t );
        hasRhoAnaR2 = true;
    }

    std::cout << "Starting to book histograms..." << std::endl;

    TH1F *rhoAverage = new TH1F("average rho all", "<#rho> in |#eta|<5.191;<#rho>_{tot};", 50, 0, 250);
    TH1F *rhoAverageEta = new TH1F("average rho acceptance", "<#rho> in |#eta|<2.043;<#rho>_{acc};", 50, 0, 250);
    
    TH1F *randomConeRawPt4 = new TH1F("random cone R=0.4 raw pt", "random cone R=0.4 raw p_{T};p^{cone}_{T,raw};", 100, 0, 200);
    TH1F *randomConePtAxis4 = new TH1F("random cone R=0.4 pt axis corr", "random cone R=0.4 axis corrected p_{T};p^{cone}_{T,axis};", 100, -100, 100);
    TH1F *randomConePtBands4 = new TH1F("random cone R=0.4 pt band corr", "random cone R=0.4 band corrected p_{T};p^{cone}_{T,band};", 100, -100, 100);

    TH1F *randomConeRawPt2 = new TH1F("random cone R=0.2 raw pt", "random cone R=0.2 raw p_{T};p^{cone}_{T,raw};", 100, 0, 200);
    TH1F *randomConePtAxis2 = new TH1F("random cone R=0.2 pt axis corr", "random cone R=0.2 axis corrected p_{T};p^{cone}_{T,axis};", 100, -100, 100);
    TH1F *randomConePtBands2 = new TH1F("random cone R=0.2 pt band corr", "random cone R=0.2 band corrected p_{T};p^{cone}_{T,band};", 100, -100, 100);

    TH1F *hist_vz = new TH1F("vz", "v_z;v_z;", 100, -20, 20);

    
    TH2F *rho_vs_centrTot = new TH2F("rho vs centrality all", "#rho (|#eta|<5) vs centrality;HiBin;#rho", 201, 0, 201, 100, 0, 250);
    TH2F *rho_vs_centrAcc = new TH2F("rho vs centrality acceptance", "#rho (|#eta|<2) vs centrality;HiBin;#rho", 201, 0, 201, 100, 0, 250);
    TH2F *raw_cone_pt_vs_centr4 = new TH2F("raw R=0.4 cone pt vs centrality", "raw R=0.4 cone p_{T} vs centrality;HiBin;p^{cone}_{T,raw}", 201, 0, 201, 100, 0, 200);
    TH2F *axis_corr_pt_vs_centr4 = new TH2F("axis corr. R=0.4 cone pt vs centrality", "axis corr. R=0.4 cone p_{T} vs centrality;HiBin;p^{cone}_{T,axis}", 201, 0, 201, 100, -100, 100);
    TH2F *raw_cone_pt_vs_centr2 = new TH2F("raw R=0.2 cone pt vs centrality", "raw R=0.2 cone p_{T} vs centrality;HiBin;p^{cone}_{T,raw}", 201, 0, 201, 100, 0, 200);
    TH2F *axis_corr_pt_vs_centr2 = new TH2F("axis corr. R=0.2 cone pt vs centrality", "axis corr. R=0.2 cone p_{T} vs centrality;HiBin;p^{cone}_{T,axis}", 201, 0, 201, 100, -100, 100);
    
    TH1F *centrality = new TH1F("centrality", "centrality;centrality;", 200, 0, 200);

    std::vector < Int_t > band_lowbin = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180};
    std::vector < Int_t > band_highbin = {20, 40, 60, 80, 100, 120, 140, 160, 180, 200};

    std::vector<TH1F*> R2_cone_axis_corr_pt;
    std::vector<TH1F*> R4_cone_axis_corr_pt;

    std::cout << "Before vector of histograms..." << std::endl;

    std::vector<TH1F*> rhoAllBins;
    for(int i{0}; i < 82; ++i){
        TH1F *f = new TH1F(("rho_in_bin_" + std::to_string(i)).c_str(), ("rho_in_bin_" + std::to_string(i)).c_str(), 100, 0, 300);
        rhoAllBins.push_back(f);
    }
    for(int i{0}; i<band_lowbin.size(); ++i){
        TH1F *a;
        TH1F *b;
        if(i > 3 && i <= 5)
            a = new TH1F(("R=0.4 axis corr. pt hiBin [" + std::to_string(band_lowbin.at(i)) + "-" + std::to_string(band_highbin.at(i)) + "]").c_str(), ("R=0.4 axis corr. pt hiBin [" + std::to_string(band_lowbin.at(i)) + "-" + std::to_string(band_highbin.at(i)) + "];p^{cone}_{T,axis};").c_str(), 100, -40, 40);
        else if( i > 5)
            a = new TH1F(("R=0.4 axis corr. pt hiBin [" + std::to_string(band_lowbin.at(i)) + "-" + std::to_string(band_highbin.at(i)) + "]").c_str(), ("R=0.4 axis corr. pt hiBin [" + std::to_string(band_lowbin.at(i)) + "-" + std::to_string(band_highbin.at(i)) + "];p^{cone}_{T,axis};").c_str(), 100, -20, 20);
        else
            a = new TH1F(("R=0.4 axis corr. pt hiBin [" + std::to_string(band_lowbin.at(i)) + "-" + std::to_string(band_highbin.at(i)) + "]").c_str(), ("R=0.4 axis corr. pt hiBin [" + std::to_string(band_lowbin.at(i)) + "-" + std::to_string(band_highbin.at(i)) + "];p^{cone}_{T,axis};").c_str(), 100, -100, 100);
        R4_cone_axis_corr_pt.push_back(a);
        if(i > 5)
            b = new TH1F(("R=0.2 axis corr. pt hiBin [" + std::to_string(band_lowbin.at(i)) + "-" + std::to_string(band_highbin.at(i)) + "]").c_str(), ("R=0.2 axis corr. pt hiBin [" + std::to_string(band_lowbin.at(i)) + "-" + std::to_string(band_highbin.at(i)) + "];p^{cone}_{T,axis};").c_str(), 100, -5, 5);
        else
            b = new TH1F(("R=0.2 axis corr. pt hiBin [" + std::to_string(band_lowbin.at(i)) + "-" + std::to_string(band_highbin.at(i)) + "]").c_str(), ("R=0.2 axis corr. pt hiBin [" + std::to_string(band_lowbin.at(i)) + "-" + std::to_string(band_highbin.at(i)) + "];p^{cone}_{T,axis};").c_str(), 100, -30, 30);
        R2_cone_axis_corr_pt.push_back(b);
    }


    std::cout << "All histograms booked, loading branches..." << std::endl;
    Int_t   trigger_MinBias = 1;
    Int_t   trigger_ZeroBias = 1;
    Int_t   trigger=1;
    Int_t   hiBin=0;
    float   vz=0;
    Int_t   nref= 0;
    float   weight=1;
    float   jtdyn_kt[MAXJETS];
    float   jtdyn_theta[MAXJETS];
    float   jtdyn_var[MAXJETS];
    Int_t   jtdyn_split[MAXJETS];
    float   jtpt[MAXJETS];
    float   jteta[MAXJETS];
    float   jtphi[MAXJETS];

    std::vector<double> *etaMax = 0;
    std::vector<double> *etaMin = 0;
    std::vector<double> *rho = 0;
    Double_t rcone_eta;
    Double_t rcone_phi;
    Double_t rcone_pt_raw;
    Double_t rcone_pt_axisCorr;
    Double_t rcone_pt_bandsCorr;

    Double_t rcone_eta2;
    Double_t rcone_phi2;
    Double_t rcone_pt_raw2;
    Double_t rcone_pt_axisCorr2;
    Double_t rcone_pt_bandsCorr2;

    Int_t pclusterCompatibilityFilter = -1;
    Int_t pphfCoincFilter2Th4 = -1;
    Int_t pprimaryVertexFilter = -1;

    t->SetBranchAddress("jtpt", &jtpt);

	w_t->SetBranchAddress("hiBin", &hiBin);
    w_t->SetBranchAddress("vz", &vz);
	//if it's MC, load weights and truth level varables
    rho_t->SetBranchAddress("etaMax", &etaMax);
    rho_t->SetBranchAddress("etaMin", &etaMin);
    rho_t->SetBranchAddress("rho", &rho);
    if(hasRhoAnaR4){
        randomCone4_t->SetBranchAddress("rcone_eta", &rcone_eta);
        randomCone4_t->SetBranchAddress("rcone_phi", &rcone_phi);
        randomCone4_t->SetBranchAddress("rcone_pt_raw", &rcone_pt_raw);
        randomCone4_t->SetBranchAddress("rcone_pt_axisCorr", &rcone_pt_axisCorr);
        randomCone4_t->SetBranchAddress("rcone_pt_bandsCorr", &rcone_pt_bandsCorr);
    }
    if(hasRhoAnaR2){
        randomCone2_t->SetBranchAddress("rcone_eta", &rcone_eta2);
        randomCone2_t->SetBranchAddress("rcone_phi", &rcone_phi2);
        randomCone2_t->SetBranchAddress("rcone_pt_raw", &rcone_pt_raw2);
        randomCone2_t->SetBranchAddress("rcone_pt_axisCorr", &rcone_pt_axisCorr2);
        randomCone2_t->SetBranchAddress("rcone_pt_bandsCorr", &rcone_pt_bandsCorr2);
    }

    hlt_t->SetBranchAddress("HLT_HIMinimumBiasHF1AND_v1",&trigger_MinBias);
    if( (input_file.find( "PbPb" )) != std::string::npos || ( input_file.find("embedded")) != std::string::npos || (input_file.find("HYDIJET"))!= std::string::npos || ( input_file.find("data")) != std::string::npos ){
        std::cout << "Loading filters" << std::endl;
        add_sel_t->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter);
        add_sel_t->SetBranchAddress("pphfCoincFilter2Th4", &pphfCoincFilter2Th4);
        add_sel_t->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter);
    }

    for (Int_t j{0}; j < t->GetEntries(); j++ ){
        if(j%1000 == 0){
            std::cout << "Processing entry " << j << " / " << t->GetEntries() << std::endl;
        }
        t->GetEntry(j);
        w_t->GetEntry(j);
        rho_t->GetEntry(j);
        if(hasRhoAnaR4)
            randomCone4_t->GetEntry(j);
        if(hasRhoAnaR2)
            randomCone2_t->GetEntry(j);
        if( !isMC ){
        	hlt_t->GetEntry(j);
        }
        if( trigger == -1 ){
            std::cout << "Selected trigger is not filled/value is -1! Exiting" << std::endl;
            return;
        }
        if( trigger != 1 )
            continue;
        if( trigger_MinBias != 1 )
            continue;
        if( ( input_file.find( "data" ) ) != std::string::npos || ( input_file.find( "PbPb" ) ) != std::string::npos || ( input_file.find( "embedded" ) ) != std::string::npos || ( input_file.find( "HYDIJET" ) ) != std::string::npos)
            add_sel_t->GetEntry(j);
        if( pprimaryVertexFilter != 1 || pclusterCompatibilityFilter != 1 || pphfCoincFilter2Th4 != 1 ){
            std::cout << "Throwing out event due to filter" << std::endl;
            std::cout << pprimaryVertexFilter << " " << pclusterCompatibilityFilter << " " << pphfCoincFilter2Th4 << std::endl;
            continue;
        }
        //average rho over whole range and acceptance range eta in ~(-2;2)
        Double_t rho_accu = 0;
        Double_t rho_accu_acc = 0;
        Int_t accepted_bins = 1;
        for(size_t r{0};r < rho->size(); ++r){
            rho_accu += rho->at(r);
            if(etaMin->at(r)>= -2.043 && etaMax->at(r) <= 2.043){
                rho_accu_acc += rho->at(r);
                accepted_bins++;
            }
        }
        rho_accu = rho_accu/rho->size();
        rho_accu_acc = rho_accu_acc/accepted_bins;
        //filling rho histogram for every eta band
        for(size_t r{0};r < rho->size(); ++r){
            rhoAllBins.at(r)->Fill(rho->at(r), weight);
            // rhoAllBins.at(r)->SetName(("rho_in_bin_" + std::to_string(r)).c_str());
            std::stringstream s1;
            s1 << std::fixed << std::setprecision(3) << etaMin->at(r);
            std::stringstream s2;
            s2 << std::fixed << std::setprecision(3) << etaMax->at(r);
            rhoAllBins.at(r)->SetTitle(("#rho in bin " + std::to_string(r) + " #eta #in (" + s1.str() + "," + s2.str() + ")").c_str());
            rhoAllBins.at(r)->GetXaxis()->SetTitle("#rho");
        }
        centrality->Fill(hiBin, weight);
        rhoAverage->Fill(rho_accu, weight);
        rhoAverageEta->Fill(rho_accu_acc, weight);
        randomConeRawPt4->Fill(rcone_pt_raw, weight);
        randomConePtAxis4->Fill(rcone_pt_axisCorr, weight);
        randomConePtBands4->Fill(rcone_pt_bandsCorr, weight);
        randomConeRawPt2->Fill(rcone_pt_raw2, weight);
        randomConePtAxis2->Fill(rcone_pt_axisCorr2, weight);
        randomConePtBands2->Fill(rcone_pt_bandsCorr2, weight);
        rho_vs_centrTot->Fill(hiBin, rho_accu, weight);
        raw_cone_pt_vs_centr2 ->Fill(hiBin, rcone_pt_raw2, weight);
        axis_corr_pt_vs_centr2->Fill(hiBin, rcone_pt_axisCorr2, weight);
        raw_cone_pt_vs_centr4 ->Fill(hiBin, rcone_pt_raw, weight);
        axis_corr_pt_vs_centr4 ->Fill(hiBin, rcone_pt_axisCorr, weight);

        Int_t hist_number = histogram_number(band_lowbin, band_highbin, hiBin);
        if(rcone_pt_raw2 != 0){
            R2_cone_axis_corr_pt.at(hist_number)->Fill(rcone_pt_axisCorr2, weight);
        }
        if(rcone_pt_raw != 0){
            R4_cone_axis_corr_pt.at(hist_number)->Fill(rcone_pt_axisCorr, weight);
        }
        hist_vz->Fill(vz, weight);
	}

	TFile *of = new TFile(output_file.c_str(),"RECREATE");
	of->cd();
	for(size_t i{0}; i<rhoAllBins.size();++i){
        rhoAllBins.at(i)->Write();
    }
    for(size_t i{0}; i<R2_cone_axis_corr_pt.size();++i){
        R2_cone_axis_corr_pt.at(i)->Write();
    }
    for(size_t i{0}; i<R4_cone_axis_corr_pt.size();++i){
        R4_cone_axis_corr_pt.at(i)->Write();
    }
    rhoAverage->Write();
    rhoAverageEta->Write();
    randomConeRawPt4->Write();
    randomConePtAxis4->Write();
    randomConePtBands4->Write();
    randomConeRawPt2->Write();
    randomConePtAxis2->Write();
    randomConePtBands2->Write();

    hist_vz->Write();
    centrality->Write();
    rho_vs_centrTot->Write();
    raw_cone_pt_vs_centr2 ->Write();
    axis_corr_pt_vs_centr2->Write();
    raw_cone_pt_vs_centr4 ->Write();
    axis_corr_pt_vs_centr4 ->Write();

	of->Close();
	delete t;

	std::cout << "All done!" << std::endl;
	gSystem->Exit(0);
}
