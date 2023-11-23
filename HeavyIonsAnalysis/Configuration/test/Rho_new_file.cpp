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

float findNcoll(int hiBin) {
   const int nbins = 200;
   const vector<float> Ncoll = {1976.95, 1944.02, 1927.29, 1891.9, 1845.3, 1807.2, 1760.45, 1729.18, 1674.8, 1630.3, 1590.52, 1561.72, 1516.1, 1486.5, 1444.68, 1410.88, 1376.4, 1347.32, 1309.71, 1279.98, 1255.31, 1219.89, 1195.13, 1165.96, 1138.92, 1113.37, 1082.26, 1062.42, 1030.6, 1009.96, 980.229, 955.443, 936.501, 915.97, 892.063, 871.289, 847.364, 825.127, 806.584, 789.163, 765.42, 751.187, 733.001, 708.31, 690.972, 677.711, 660.682, 640.431, 623.839, 607.456, 593.307, 576.364, 560.967, 548.909, 530.475, 519.575, 505.105, 490.027, 478.133, 462.372, 451.115, 442.642, 425.76, 416.364, 405.154, 392.688, 380.565, 371.167, 360.28, 348.239, 340.587, 328.746, 320.268, 311.752, 300.742, 292.172, 281.361, 274.249, 267.025, 258.625, 249.931, 240.497, 235.423, 228.63, 219.854, 214.004, 205.425, 199.114, 193.618, 185.644, 180.923, 174.289, 169.641, 161.016, 157.398, 152.151, 147.425, 140.933, 135.924, 132.365, 127.017, 122.127, 117.817, 113.076, 109.055, 105.16, 101.323, 98.098, 95.0548, 90.729, 87.6495, 84.0899, 80.2237, 77.2201, 74.8848, 71.3554, 68.7745, 65.9911, 63.4136, 61.3859, 58.1903, 56.4155, 53.8486, 52.0196, 49.2921, 47.0735, 45.4345, 43.8434, 41.7181, 39.8988, 38.2262, 36.4435, 34.8984, 33.4664, 31.8056, 30.351, 29.2074, 27.6924, 26.7754, 25.4965, 24.2802, 22.9651, 22.0059, 21.0915, 19.9129, 19.1041, 18.1487, 17.3218, 16.5957, 15.5323, 14.8035, 14.2514, 13.3782, 12.8667, 12.2891, 11.61, 11.0026, 10.3747, 9.90294, 9.42648, 8.85324, 8.50121, 7.89834, 7.65197, 7.22768, 6.7755, 6.34855, 5.98336, 5.76555, 5.38056, 5.11024, 4.7748, 4.59117, 4.23247, 4.00814, 3.79607, 3.68702, 3.3767, 3.16309, 2.98282, 2.8095, 2.65875, 2.50561, 2.32516, 2.16357, 2.03235, 1.84061, 1.72628, 1.62305, 1.48916, 1.38784, 1.28366, 1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};
   return Ncoll.at(hiBin);
}

Double_t findRhoWeight(Double_t rho){
    const vector<Double_t> RatioBinLowEdge = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245};
    const vector<Double_t> RatioBinContent = {0, 0, 0, 0, 0, 0.00145813, 0, 0, 0.00242926, 0.0849748, 0.338255, 0.737106, 0.942886, 1.02698, 1.04522, 0.986484, 1.05535, 1.14385, 1.18673, 1.31314, 1.30534, 1.48038, 1.41024, 1.48927, 1.59019, 1.72233, 1.84837, 1.86471, 1.92942, 2.01884, 2.06943, 1.973, 2.22433, 2.28652, 2.35816, 2.23266, 2.25113, 1.90066, 1.5782, 1.01255, 0.804272, 0.438775, 0.286191, 0.160913, 0.170254, 0.111196, 0.124253, 0.300739, 0.226383, 0.675163};
    // no trigger, 300 GeV jet req. 
    // const vector<Double_t> RatioBinContent = {0, 0, 0, 0, 0, 0.000209806, 0, 0, 0.00406503, 0.0696658, 0.422128, 0.904297, 1.24741, 1.40316, 1.51074, 1.64154, 1.75954, 1.96003, 2.22357, 2.43028, 2.62443, 2.75196, 2.92921, 3.24549, 3.31155, 3.46401, 3.71232, 3.86192, 3.93444, 4.16549, 4.3434, 4.19178, 4.52838, 4.66675, 4.81219, 4.64101, 4.52978, 3.73813, 3.14476, 2.24314, 1.57588, 1.05741, 0.588536, 0.332049, 0.268398, 0.125907, 0.173472, 0.271159, 0.571597, 0.96026};

    // no trigger ratios
    // const vector<Double_t> RatioBinContent = {0, 0, 0, 0, 1.8713e-05, 3.98462e-05, 7.36656e-05, 0.000356741, 0.0106945, 0.139252, 0.580562, 1.02557, 1.22079, 1.25531, 1.23058, 1.19959, 1.16536, 1.21154, 1.27643, 1.27865, 1.26867, 1.27884, 1.31107, 1.30472, 1.27864, 1.31418, 1.32454, 1.31515, 1.33055, 1.3335, 1.32517, 1.3505, 1.37704, 1.40912, 1.40224, 1.38693, 1.25505, 1.0807, 0.869455, 0.648113, 0.446237, 0.287623, 0.177387, 0.105519, 0.0676355, 0.0509648, 0.056936, 0.0805727, 0.148572, 0.288392};
    // old ratio weights with 100GeV trigger
    // const vector<Double_t> RatioBinContent = {0, 0, 0, 0, 3.09562e-05, 6.93856e-05, 0.000115449, 0.00045351, 0.0117345, 0.168105, 0.747199, 1.35108, 1.60675, 1.63157, 1.57346, 1.50544, 1.43521, 1.46527, 1.51567, 1.48493, 1.45047, 1.43291, 1.44665, 1.40814, 1.35938, 1.36593, 1.34149, 1.29774, 1.27605, 1.24172, 1.19238, 1.16877, 1.14627, 1.12507, 1.07244, 1.01417, 0.881134, 0.735451, 0.568816, 0.411977, 0.275729, 0.174529, 0.103914, 0.0584471, 0.0344026, 0.0206139, 0.016094, 0.0171954, 0.0312582, 0.0441074};
    //assume bins are of equal width
    const Double_t BinWidth = RatioBinLowEdge.at(1) - RatioBinLowEdge.at(0);
    for(size_t i{0}; i<RatioBinContent.size(); ++i){
        if(rho > RatioBinLowEdge.at(i) && rho < RatioBinLowEdge.at(i)+BinWidth){
            return RatioBinContent.at(i);
        }
    }
    return 0;
}

void Rho_new_file(std::string input_file, std::string output_file){
    Bool_t isMC = false;
    Bool_t ispp = false;
    if( ( input_file.find( "embedded" ) ) != std::string::npos || ( input_file.find( "pythia" ) ) != std::string::npos || ( input_file.find( "HYDJET" ) ) != std::string::npos ){
        isMC = true;
    }
    if( ( input_file.find( "pp" ) ) != std::string::npos ){
        ispp = true;
    }
	// if( R!="2" && R!="3" && R!="4" ){
	// 	std::cout<<"Expect argument 2,3 or 4, exiting" << std::endl;
	// 	return;
	// }
    std::cout << "Starting the processing..." << std::endl;
	const Int_t MAXJETS=100;
	// Bool_t isMC = false;
	// options my_opt;
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

    TH1F *rho300400AverageEtaAcc = new TH1F("AE average rho acceptance 300-400", "Jet acc. <#rho> in |#eta|<2.043 300-400;<#rho>_{acc};", 50, 0, 250);
    TH1F *rhoAverageAcc = new TH1F("AE average rho all", "Jet acc. <#rho> in |#eta|<5.191;<#rho>_{tot};", 50, 0, 250);
    TH1F *rhoAverageEtaAcc = new TH1F("AE average rho acceptance", "Jet acc. <#rho> in |#eta|<2.043;<#rho>_{acc};", 50, 0, 250);
    TH1F *randomConeRawPtAcc4 = new TH1F("AE random cone R=0.4 raw pt", "Jet acc. random cone R=0.4 raw p_{T};p^{cone}_{T,raw};", 100, 0, 200);
    TH1F *randomConePtAxisAcc4 = new TH1F("AE random cone R=0.4 pt axis corr", "Jet acc. random cone R=0.4 axis corrected p_{T};p^{cone}_{T,axis};", 100, -100, 100);
    TH1F *randomConePtBandsAcc4 = new TH1F("AE random cone R=0.4 pt band corr", "Jet acc. random cone R=0.4 band corrected p_{T};p^{cone}_{T,band};", 100, -100, 100);
    
    TH1F *randomConeRawPtAcc2 = new TH1F("AE random cone R=0.2 raw pt", "Jet acc. random cone R=0.2 raw p_{T};p^{cone}_{T,raw};", 100, 0, 200);
    TH1F *randomConePtAxisAcc2 = new TH1F("AE random cone R=0.2 pt axis corr", "Jet acc. random cone R=0.2 axis corrected p_{T};p^{cone}_{T,axis};", 100, -100, 100);
    TH1F *randomConePtBandsAcc2 = new TH1F("AE random cone R=0.2 pt band corr", "Jet acc. random cone R=0.2 band corrected p_{T};p^{cone}_{T,band};", 100, -100, 100);
    

    TH1F *hist_vz = new TH1F("vz", "v_z;v_z;", 100, -20, 20);
    TH1F *hist_vzAcc = new TH1F("AE vz", "Jet acc. v_z;v_z;", 100, -20, 20);

    
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

    std::vector<TH1F*> rho_in_pt;
    for(Int_t i{1}; i < 4; ++i){
        TH1F *rpt = new TH1F(("Rho for pT>" + std::to_string(i*100)).c_str(), ("#rho for p_{T}>" + std::to_string(i*100)).c_str(), 50, 0, 250);
        rho_in_pt.push_back(rpt);
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

    float   refdyn_kt[MAXJETS];
    float   refdyn_theta[MAXJETS];
    float   refdyn_var[MAXJETS];
    Int_t   refdyn_split[MAXJETS];
    float   refpt[MAXJETS];
    float   refeta[MAXJETS];
    float   refphi[MAXJETS];

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
    // t->SetBranchAddress("jteta", &jteta);
    // t->SetBranchAddress("jtphi", &jtphi);
    // t->SetBranchAddress("nref", &nref);
    // t->SetBranchAddress("jtdyn_kt", &jtdyn_kt);
    // t->SetBranchAddress("jtdyn_theta", &jtdyn_theta);
    // t->SetBranchAddress("jtdyn_var", &jtdyn_var);
	// t->SetBranchAddress("jtdyn_split", &jtdyn_split);

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
    if(isMC){
    	w_t->SetBranchAddress("weight",&weight);
	    t->SetBranchAddress("refpt", &refpt);
	    t->SetBranchAddress("refeta", &refeta);
	    t->SetBranchAddress("refphi", &refphi);
	    t->SetBranchAddress("refdyn_kt", &refdyn_kt);
	    t->SetBranchAddress("refdyn_theta", &refdyn_theta);
	    t->SetBranchAddress("refdyn_var", &refdyn_var);
	    t->SetBranchAddress("refdyn_split", &refdyn_split);
	}
	//check whether the filepath contains a hint of what the input dataset is
	// if( ( input_file.find( "PbPb" ) ) != std::string::npos ){
	// 	// hlt_t->SetBranchAddress("HLT_HIPuAK4CaloJet100Eta5p1_v1",&trigger);
        
    // }
    // if( ( input_file.find( "Run3_data" ) ) != std::string::npos ){
    //     hlt_t->SetBranchAddress("HLT_HIMinimumBias_v2",&trigger_MinBias);  
    // }
    hlt_t->SetBranchAddress("HLT_HIMinimumBiasHF1AND_v1",&trigger_MinBias);  
    if( (input_file.find( "PbPb" )) != std::string::npos || ( input_file.find("embedded")) != std::string::npos || (input_file.find("HYDIJET"))!= std::string::npos || ( input_file.find("data")) != std::string::npos ){
        std::cout << "Loading filters" << std::endl;
        add_sel_t->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter);
        add_sel_t->SetBranchAddress("pphfCoincFilter2Th4", &pphfCoincFilter2Th4);
        add_sel_t->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter);
    }

	else if( ( input_file.find( "pp" ) ) != std::string::npos )
		hlt_t->SetBranchAddress("HLT_HIAK4CaloJet100_v1", &trigger);
	// else if(( input_file.find( "MC" )) != std::string::npos || ( input_file.find( "pythia" )) != std::string::npos || ( input_file.find( "embedded" )) != std::string::npos || || ( input_file.find( "HYDIJET" )) != std::string::npos){
	// 	//do nothing
	// }
	// else{
	// 	std::cout << "Not pp, PbPb, embedded (or MC), exiting."<<std::endl;
	// 	return;
	// }

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
        // std::cout << "Event " << j << " weight " << weight << std::endl;
        if( trigger == -1 ){
            std::cout << "Selected trigger is not filled/value is -1! Exiting" << std::endl;
            return;
        }
        // if( weight > 2 ){
        //     std::cout << "event weight is " << weight << std::endl;
        // }
        if( trigger != 1 )
            continue;
        if( trigger_MinBias != 1 )
            continue;
        if( (input_file.find("HIHardProbes")) != std::string::npos && hiBin > 60 )
            continue;
        if(isMC){
           weight = weight*findNcoll(hiBin);
        }
        if ( (input_file.find("Run3_data")) != std::string::npos ){
        }
        if( ( input_file.find( "data" ) ) != std::string::npos || ( input_file.find( "PbPb" ) ) != std::string::npos || ( input_file.find( "embedded" ) ) != std::string::npos || ( input_file.find( "HYDIJET" ) ) != std::string::npos)
            add_sel_t->GetEntry(j);
        if( pprimaryVertexFilter != 1 || pclusterCompatibilityFilter != 1 || pphfCoincFilter2Th4 != 1 ){
            std::cout << "Throwing out event due to filter" << std::endl;
            std::cout << pprimaryVertexFilter << " " << pclusterCompatibilityFilter << " " << pphfCoincFilter2Th4 << std::endl;
            continue;
        }
        
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

        if(isMC){
            weight = weight*findRhoWeight(rho_accu_acc);
        }
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
        rho_vs_centrAcc->Fill(hiBin,rho_accu_acc, weight);
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

        Bool_t jt_above_100 = false;
        Bool_t jt_above_200 = false;
        Bool_t jt_above_300 = false;
        hist_vz->Fill(vz, weight);
        Bool_t eventAccepted = false;
		for(Int_t k = 0; k < nref; k++){
			if ( TMath::Abs(jteta[k]) > 1.6 ) continue;
			if ( TMath::Abs(refeta[k]) > 2 ) continue;
            if ( weight > 0.001 && isMC ) continue;
            if ( jtpt[k] > 100 )
                jt_above_100 = true;
            if( jtpt[k] > 200 )
                jt_above_200 = true;
            if( jtpt[k] > 300 )
                jt_above_300 = true;
            if( jtpt[k]>300 && jtpt[k]<400)
                rho300400AverageEtaAcc->Fill(rho_accu_acc, weight);
			if ( jtpt[k] < 300 ) continue;
            eventAccepted = true;
            
		}

        if ( jt_above_100 )
            rho_in_pt.at(0)->Fill(rho_accu_acc, weight);
        if ( jt_above_200 )
            rho_in_pt.at(1)->Fill(rho_accu_acc, weight);
        if ( jt_above_300 )
            rho_in_pt.at(2)->Fill(rho_accu_acc, weight);
        
        if(eventAccepted){
            rhoAverageAcc->Fill(rho_accu, weight);
            rhoAverageEtaAcc->Fill(rho_accu_acc, weight);
            randomConeRawPtAcc4->Fill(rcone_pt_raw, weight);
            randomConePtAxisAcc4->Fill(rcone_pt_axisCorr, weight);
            randomConePtBandsAcc4->Fill(rcone_pt_bandsCorr, weight);
            randomConeRawPtAcc2->Fill(rcone_pt_raw2, weight);
            randomConePtAxisAcc2->Fill(rcone_pt_axisCorr2, weight);
            randomConePtBandsAcc2->Fill(rcone_pt_bandsCorr2, weight);
            hist_vzAcc->Fill(vz, weight);
        }
	}
	// std::string outfile = "";
	// outfile = outfile + "_kappa_switch.root";
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
    for(size_t i{0}; i<rho_in_pt.size();++i){
        rho_in_pt.at(i)->Write();
    }
    rhoAverage->Write();
    rhoAverageEta->Write();
    rhoAverageAcc->Write();
    rhoAverageEtaAcc->Write();
    rho300400AverageEtaAcc->Write();
    randomConeRawPt4->Write();
    randomConePtAxis4->Write();
    randomConePtBands4->Write();
    randomConeRawPtAcc4->Write();
    randomConePtAxisAcc4->Write();
    randomConePtBandsAcc4->Write();
    randomConeRawPt2->Write();
    randomConePtAxis2->Write();
    randomConePtBands2->Write();
    randomConeRawPtAcc2->Write();
    randomConePtAxisAcc2->Write();
    randomConePtBandsAcc2->Write();
    hist_vz->Write();
    hist_vzAcc->Write();
    centrality->Write();
    rho_vs_centrAcc->Write();
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
