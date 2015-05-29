#include <cstdlib>
#include <fstream>
#include <vector>
#include <set>
#include <map>

#include "TString.h"
#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

#include "../interface/Retina.h"
#include "../interface/Track.h"


struct GenPart_t {
  int    pdg;
  double c;
  double pt;
  double d0;
  double phi;
  double eta;
  double z0;
  double x_P_XY;
  double x_M_XY;
  double x_P_RZ;
  double x_M_RZ;
};

struct GenPart_compare {
  bool operator() (const GenPart_t &a, const GenPart_t &b) const {
    return a.pt < b.pt;
  }
};


struct Pattern_t {
  float pt;
  float pt_rms;
  float phi;
  float phi_rms;
  float eta;
  float eta_rms;
};

std::map<const char*,Pattern_t> road_par;


// --- Function paradigms:
void rotateHits(vector<Hit_t*> hits, double angle);
void confTrans(vector<Hit_t*> hits);
void retina_initialize(std::map <const char*,double> *config);


// --- Retina configuration parameters:
std::map <const char*,double> config[3];

// For -pi < phi < pi
const double rot_angle[8] = {  0.39269908169872414,   //  pi/8
			      -0.39269908169872414,   // -pi/8
			      -1.17809724509617242,   // -3/8 pi
			      -1.96349540849362070,   // -5/8 pi
			       3.53429173528851726,   //  9/8 pi
			       2.74889357189106898,   //  7/8 pi
			       1.96349540849362070,   //  5/8 pi
			       1.17809724509617242 }; //  3/8 pi
// For 0 < phi < 2 pi
//const double rot_angle[8] = {  0.39269908169872414,   //  pi/8
//			      -0.39269908169872414,   // -pi/8
//			      -1.17809724509617242,   // -3/8 pi
//			      -1.96349540849362070,   // -5/8 pi
//			      -2.74889357189106898,   // -7/8 pi
//			      -3.53429173528851726,   // -9/8 pi
//			      -4.31968989868596509,   // -11/8 pi
//			      -5.10508806208341426 }; // -13/8 pi


const int active_layers[16] = { 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 18, 19, 20, 21, 22 };

const double mMagneticField = 3.8;

const int n_scan_points = 1;

int Nroads_trigT[48];
int Ncomb_trigT[48];


// --- Histogram booking:

TH1F * h_Nstubs_tot    = new TH1F("h_Nstubs_tot","Total stubs per event",1000,0,20000);
TH1F * h_Nroads_tot    = new TH1F("h_Nroads_tot","Total roads per event",1000,0,60000);
TH1F * h_Nstubs_tot_PR = new TH1F("h_Nstubs_tot_PR","Total stubs per event after PR",1000,0,20000);
TH1F * h_Nstubs_occ_PR = new TH1F("h_Nstubs_occ_PR","Number of stubs vs trigger tower",48,0,48);
TH1F * h_Nstubs_road   = new TH1F("h_Nstubs_road","Total stubs per road",400,0,400);
TH1F * h_trigT_occ     = new TH1F("h_trigT_occ","Number of roads vs trigger tower",48,0,48);
TH1F * h_Ncomb_road_56    = new TH1F("h_Ncomb_road_56","Number of stub combinations per road",3000,0,3000);
TH1F * h_Ncomb_road_56_XY = new TH1F("h_Ncomb_road_56_XY","Number of stub combinations per road (after XY retina)",3000,0,3000);
TH1F * h_Ncomb_road_56_RZ = new TH1F("h_Ncomb_road_56_RZ","Number of stub combinations per road (after #rhoZ retina)",3000,0,3000);
TH1F * h_Ncomb_road_56_PE = new TH1F("h_Ncomb_road_56_PE","Number of stub combinations per road (after #phi#eta retina)",3000,0,3000);
TH1F * h_Ncomb_road_66    = new TH1F("h_Ncomb_road_66","Number of stub combinations per road",3000,0,3000);
TH1F * h_Ncomb_road_66_XY = new TH1F("h_Ncomb_road_66_XY","Number of stub combinations per road (after XY retina)",3000,0,3000);
TH1F * h_Ncomb_road_66_RZ = new TH1F("h_Ncomb_road_66_RZ","Number of stub combinations per road (after #rhoZ retina)",3000,0,3000);
TH1F * h_Ncomb_road_66_PE = new TH1F("h_Ncomb_road_66_PE","Number of stub combinations per road (after #phi#eta retina)",3000,0,3000);
TH1F * h_Nmax_road_56_XY = new TH1F("h_Nmax_road_56_XY","Number of maxima (XY retina)",100,0,100);
TH1F * h_Nmax_road_56_RZ = new TH1F("h_Nmax_road_56_RZ","Number of maxima (#rhoZ retina)",100,0,100);
TH1F * h_Nmax_road_56_PE = new TH1F("h_Nmax_road_56_PE","Number of maxima (#phi#eta retina)",100,0,100);
TH1F * h_Nmax_road_66_XY = new TH1F("h_Nmax_road_66_XY","Number of maxima (XY retina)",100,0,100);
TH1F * h_Nmax_road_66_RZ = new TH1F("h_Nmax_road_66_RZ","Number of maxima (#rhoZ retina)",100,0,100);
TH1F * h_Nmax_road_66_PE = new TH1F("h_Nmax_road_66_PE","Number of maxima (#phi#eta retina)",100,0,100);


TH1F * h_Ncomb_trigT   = new TH1F("h_Ncomb_trigT","Number of combinations per trigger tower",100,0,3000);

TH1F * h_Nstubs_road_layer  = new TH1F("h_Nstubs_road_layer","Different layers per road",30,0,30);
TH1F * h_Nstubs_trigT_layer = new TH1F("h_Nstubs_trigT_layer","Different layers per trigger tower",30,0,30);

TH1F * h_Nstubs_trigT[48];
TH1F * h_Nroads_trigT[48];

TH1F * h_NstubsLayer_road[16];
TH1F * h_NstubsLayer_trigT[16];


TH1F * h_gen_pdg   = new TH1F("h_gen_pdg","gen pdg ID;pdg ID", 1000, -500, 500.);
TH1F * h_gen_pt    = new TH1F("h_gen_pt","gen p_{T};p_{T}^{gen} [GeV/c]", 100, 0.,100.);
TH1F * h_gen_d0    = new TH1F("h_gen_d0","gen d_{0};d_{0}^{gen} [cm]", 100,-0.05,0.05);
TH1F * h_gen_phi   = new TH1F("h_gen_phi","gen #phi;#phi^{gen} [rad]", 100,-TMath::Pi(),TMath::Pi());
TH1F * h_gen_eta   = new TH1F("h_gen_eta","gen #eta;#eta^{gen}", 100,-4.,4.);
TH1F * h_gen_theta = new TH1F("h_gen_theta","gen #theta;#theta^{gen}", 100,0.,TMath::Pi());
TH1F * h_gen_z0    = new TH1F("h_gen_z0","gen z_{0};z_{0}^{gen} [cm]", 100, -20.,20.);

TH2F * h_max_XY_gen[3];  
TH2F * h_max_RZ_gen[3];  
TH2F * h_max_PE_gen[3];  

TH2F * h_max_XY1[n_scan_points][3];  
TH2F * h_max_XY2[n_scan_points][3];  
TH2F * h_max_RZ1[n_scan_points][3];  
TH2F * h_max_RZ2[n_scan_points][3];  
TH2F * h_max_PE1[n_scan_points][3];  
TH2F * h_max_PE2[n_scan_points][3];  

TH1F * h_Ntrks[n_scan_points];
TH1F * h_trk_c[n_scan_points];
TH1F * h_trk_pt[n_scan_points];
TH1F * h_trk_phi[n_scan_points];
TH1F * h_trk_eta[n_scan_points];
TH1F * h_trk_z0[n_scan_points];

TH1F * h_wxy_1[n_scan_points];
TH1F * h_wxy_2[n_scan_points];
TH1F * h_wrz_1[n_scan_points];
TH1F * h_wrz_2[n_scan_points];

TH1F * h_res_pt_rel[n_scan_points];
TH1F * h_res_phi[n_scan_points];
TH1F * h_res_eta[n_scan_points];
TH1F * h_res_z0[n_scan_points];


int main(int argc, char* argv[]) {

  if ( argc<2 ){
    cout << "\n*** No arguments provided! Run \"./Retina -h\" for a list of available options.\n" << endl;
    return EXIT_SUCCESS;
  }

  long long int n_events = 1;
  bool fit_roads = true;
  string input_file = "";
  int verboseLevel = 1;

  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i],"--fit-roads")) {
      fit_roads = (bool) atoi(argv[i+1]);
    }
    if (!strcmp(argv[i],"-n")) {
      n_events = atoi(argv[i+1]);
    }
    if (!strcmp(argv[i],"-f")) {
      input_file = argv[i+1];
    }
    if (!strcmp(argv[i],"-v")) {
      verboseLevel = atoi(argv[i+1]);
    }
    if (!strcmp(argv[i],"-h")) {
      cout << "-h \t print this help and exit" << endl 
	   << "-n # \t set number of events (default = 1)" << endl 
	   << "-f <filename> \t a root file or a text file containing a list of root files (default = \"\")" 
	   << endl
	   << "--fit-roads # \t fit per road (default = 1)" << endl
	   << "-v # \t verbosity level (default = 1):" << endl
	   << "         0: no printouts" << endl
	   << "         1: print the fitted tracks and the matched gen particles" << endl
	   << "         2: print all the MC gen particles" << endl
	   << "         3: print the stub coordinates and the Retina maxima" << endl
	   << "         4: dump the Retina into a file" << endl
	   << "         5: switch on the fit diagnostics" << endl;
      return EXIT_SUCCESS;
    }
  }


  // --- Book the histogram arrays:
  for (int iscan=0; iscan<n_scan_points; ++iscan){

    TString hname  = Form("h_Ntrks_%d",iscan);
    TString htitle = Form("Number of fitted track (%d)",iscan);
    h_Ntrks[iscan] = new TH1F(hname.Data(),htitle.Data(),500,0,500); 

    hname  = Form("h_trk_c_%d",iscan);
    htitle = Form("track c (%d);c^{trk} [cm^{-1}]",iscan);
    h_trk_c[iscan] = new TH1F(hname.Data(),htitle.Data(), 150, -0.0076, 0.0076);

    hname  = Form("h_trk_pt_%d",iscan);
    htitle = Form("track p_{T} (%d);p_{T}^{trk} [GeV/c]",iscan);
    h_trk_pt[iscan] = new TH1F(hname.Data(),htitle.Data(), 150, 0.,150.);

    hname  = Form("h_trk_phi_%d",iscan);
    htitle = Form("track #phi (%d);#phi^{trk} [rad]",iscan);
    h_trk_phi[iscan] = new TH1F(hname.Data(),htitle.Data(), 100, -TMath::Pi(), TMath::Pi());

    hname  = Form("h_trk_eta_%d",iscan);
    htitle = Form("track #eta (%d);#eta^{trk}",iscan);
    h_trk_eta[iscan] = new TH1F(hname.Data(),htitle.Data(), 100,-3.,3.);

    hname  = Form("h_trk_z0_%d",iscan);
    htitle = Form("track z_{0} (%d);z_{0}^{trk} [cm]",iscan);
    h_trk_z0[iscan] = new TH1F(hname.Data(),htitle.Data(), 100,-25.,25.);

    hname  = Form("h_wxy_1_%d",iscan);
    htitle = Form("XY-1 fit weight (%d)",iscan);
    h_wxy_1[iscan] = new TH1F(hname.Data(),htitle.Data(), 100, 0.,20.);

    hname  = Form("h_wxy_2_%d",iscan);
    htitle = Form("XY-2 fit weight (%d)",iscan);
    h_wxy_2[iscan] = new TH1F(hname.Data(),htitle.Data(), 100, 0.,20.);

    hname  = Form("h_wrz_1_%d",iscan);
    htitle = Form("RZ-1 fit weight (%d)",iscan);
    h_wrz_1[iscan] = new TH1F(hname.Data(),htitle.Data(), 100, 0.,20.);

    hname  = Form("h_wrz_2_%d",iscan);
    htitle = Form("RZ-2 fit weight (%d)",iscan);
    h_wrz_2[iscan] = new TH1F(hname.Data(),htitle.Data(), 100, 0.,20.);

    hname  = Form("h_res_pt_rel_%d",iscan);
    htitle = Form("p_{T} resolution (%d);(p_{T}^{fit}-p_{T}^{gen})/p_{T}^{gen}",iscan);
    h_res_pt_rel[iscan] = new TH1F(hname.Data(),htitle.Data(), 100,-1.,1.);

    hname  = Form("h_res_phi_%d",iscan);
    htitle = Form("#phi resolution (%d);#phi^{fit}-#phi^{gen} [rad]",iscan);
    h_res_phi[iscan] = new TH1F(hname.Data(),htitle.Data(), 100,-0.05,0.05);

    hname  = Form("h_res_eta_%d",iscan);
    htitle = Form("#eta resolution (%d);#eta^{fit}-#eta^{gen}",iscan);
    h_res_eta[iscan] = new TH1F(hname.Data(),htitle.Data(), 100,-0.05,0.05);

    hname  = Form("h_res_z0_%d",iscan);
    htitle = Form("z_{0} resolution (%d);z_{0}^{fit}-z_{0}^{gen} [cm]",iscan);
    h_res_z0[iscan] = new TH1F(hname.Data(),htitle.Data(), 100,-2.,2.);

    for (int ihist=0; ihist<3; ++ihist){
      hname  = Form("h_max_XY1_%d_%d",iscan,ihist);
      htitle = Form("XY-1 maxima (%d);x_{+};x_{-}", iscan);
      h_max_XY1[iscan][ihist] = new TH2F(hname.Data(),htitle.Data(),100,-0.1,0.1,100,-0.1,0.1);
      hname  = Form("h_max_XY2_%d_%d",iscan,ihist);
      htitle = Form("XY-2 maxima (%d);x_{+};x_{-}", iscan);
      h_max_XY2[iscan][ihist] = new TH2F(hname.Data(),htitle.Data(),100,-0.1,0.1,100,-0.1,0.1);
      hname  = Form("h_max_RZ1_%d_%d",iscan,ihist);
      htitle = Form("RZ-1 maxima (%d);x_{+};x_{-}", iscan);
      h_max_RZ1[iscan][ihist] = new TH2F(hname.Data(),htitle.Data(),100,-20.,420.,100,-10.,260);
      hname  = Form("h_max_RZ2_%d_%d",iscan,ihist);
      htitle = Form("RZ-2 maxima (%d);x_{+};x_{-}", iscan);
      h_max_RZ2[iscan][ihist] = new TH2F(hname.Data(),htitle.Data(),100,-20.,420.,100,-10.,260);
    }

  } // iscan loop

  h_max_XY_gen[0] = new TH2F("h_max_XY_gen_0","XY GenPart distribution;x_{+};x_{-}",1000,-2.5,2.5,1000,-1.7,1.7);
  h_max_XY_gen[1] = new TH2F("h_max_XY_gen_1","XY GenPart distribution;x_{+};x_{-}",1000,-2.5,2.5,1000,-1.7,1.7);
  h_max_XY_gen[2] = new TH2F("h_max_XY_gen_2","XY GenPart distribution;x_{+};x_{-}",1000,-2.5,2.5,1000,-1.7,1.7);
    
  h_max_RZ_gen[0] = new TH2F("h_max_RZ_gen_0","RZ GenPart distribution;x_{+};x_{-}",100,-20.,420.,100,-10.,260);
  h_max_RZ_gen[1] = new TH2F("h_max_RZ_gen_1","RZ GenPart distribution;x_{+};x_{-}",100,-20.,420.,100,-10.,260);
  h_max_RZ_gen[2] = new TH2F("h_max_RZ_gen_2","RZ GenPart distribution;x_{+};x_{-}",100,-20.,420.,100,-10.,260);
    
  for (int itow=0; itow<48; ++itow){
    
    TString hname  = Form("h_Nstubs_trigT_%d",itow);
    TString htitle = Form("Total stubs per trigger tower (tower %d)",itow);
    h_Nstubs_trigT[itow] = new TH1F(hname.Data(),htitle.Data(),1000,0,10000);

    hname  = Form("h_Nroads_trigT_%d",itow);
    htitle = Form("Total roads per trigger tower (tower %d)",itow);
    h_Nroads_trigT[itow] = new TH1F(hname.Data(),htitle.Data(),1000,0,10000);

  }


  for (int ilayer=0; ilayer<16; ++ilayer){
    
    TString hname  = Form("h_NstubsLayer_road_%d",ilayer);
    TString htitle = Form("Stubs per layer per road (layer %d)",active_layers[ilayer]);
    h_NstubsLayer_road[ilayer]  = new TH1F(hname.Data(),htitle.Data(),30,0,30);

    hname  = Form("h_NstubsLayer_trigT_%d",ilayer);
    htitle = Form("Stubs per layer per trigger tower (layer %d)",active_layers[ilayer]);
    h_NstubsLayer_trigT[ilayer]  = new TH1F(hname.Data(),htitle.Data(),1000,0,1000);

  }


  // =============================================================================================
  //  Load the pattern parameters
  // =============================================================================================

  cout << "Loading pattern parameters ... ";
  std::ifstream infile("pattern_par_fountain_sec25.txt");
  UInt_t  counter = 0;
  if ( infile.is_open() ){

    UInt_t  index;
    TString pattern_id;
    Float_t pt_mean;
    Float_t pt_rms;
    Float_t phi_mean;
    Float_t phi_rms;
    Float_t eta_mean;
    Float_t eta_rms;
    UInt_t  n;


    while (infile >> index >> pattern_id 
	          >> pt_mean >> pt_rms 
                  >> phi_mean >> phi_rms
	          >> eta_mean >> eta_rms >> n) {

      road_par[pattern_id].pt      =  pt_mean;
      road_par[pattern_id].pt_rms  =  pt_rms;
      road_par[pattern_id].phi     =  phi_mean;
      road_par[pattern_id].phi_rms =  phi_rms;
      road_par[pattern_id].eta     =  eta_mean;
      road_par[pattern_id].eta_rms =  eta_rms;

      //cout << index   << " " << pattern_id  << " " 
      //	   << pt_mean << " " << pt_rms << " " 
      //	   << phi_mean << " " << phi_rms << " " 
      //	   << eta_mean << " " << eta_rms << " " 
      //	   << n
      //	   << endl;

      counter++;

    }

  }

  cout << counter << " lines read." << endl;

 
 
  // =============================================================================================
  //  Get the trees from the root files
  // =============================================================================================

  //bool is_extractedTree = true;

  TChain *m_L1TT = new TChain("TkStubs");   // Tree containing the stub info 
  TChain *m_PATT = new TChain("L1tracks");  // Tree containing the L1 track reco info
  TChain *m_MC   = new TChain("MC");        // Tree containing the true MC info 

  //TChain *m_Patterns = new TChain("Patterns");
  //f->GetListOfKeys()->Contains("graph")


  std::size_t found = input_file.find(".root");
  // Case 1: it's a root file
  if ( found != std::string::npos ) {
    std::cout << "Opening file " << input_file.c_str() << std::endl;
    m_L1TT->Add(input_file.c_str());
    m_PATT->Add(input_file.c_str());
    m_MC->Add(input_file.c_str());
    //m_Patterns->Add(input_file.c_str());
    //if ( m_Patterns->GetEntries()>0 )
    //  is_extractedTree = false;
  }
  // Case 2: it's a list provided into a text file
  else { 
    std::string STRING;
    std::ifstream in(input_file.c_str());
    if (!in) {
      std::cout << "Please provide a valid data filename list" << std::endl; 
      return EXIT_FAILURE;
    }    
    while (!in.eof()) {
      getline(in,STRING);
      found = STRING.find(".root");
      if ( found != std::string::npos ) {
	std::cout << "Opening file " << STRING.c_str() << std::endl;
	m_L1TT->Add(STRING.c_str());   
	m_PATT->Add(STRING.c_str());   
	m_MC->Add(STRING.c_str());   
	//m_Patterns->Add(input_file.c_str());
	//if ( m_Patterns->GetEntries()>0 )
	//  is_extractedTree = false;
      }
    }
    in.close();
  }


  int n_gen = 0; 
  std::vector<int>     *gen_pdg = NULL;
  std::vector<float>   *gen_px  = NULL;
  std::vector<float>   *gen_py  = NULL;
  std::vector<float>   *gen_pz  = NULL;
  std::vector<float>   *gen_vx  = NULL;
  std::vector<float>   *gen_vy  = NULL;
  std::vector<float>   *gen_vz  = NULL;
  int n_pu = 0; 
  std::vector<int>     *pu_pdg  = NULL;
  std::vector<float>   *pu_px   = NULL;
  std::vector<float>   *pu_py   = NULL;
  std::vector<float>   *pu_pz   = NULL;
  std::vector<float>   *pu_eta  = NULL;
  std::vector<float>   *pu_phi  = NULL;
  std::vector<float>   *pu_vx   = NULL;
  std::vector<float>   *pu_vy   = NULL;
  std::vector<float>   *pu_vz   = NULL;

  m_MC->SetBranchAddress("gen_n",   &n_gen);
  m_MC->SetBranchAddress("gen_pdg", &gen_pdg);
  m_MC->SetBranchAddress("gen_px",  &gen_px);
  m_MC->SetBranchAddress("gen_py",  &gen_py);
  m_MC->SetBranchAddress("gen_pz",  &gen_pz);
  m_MC->SetBranchAddress("gen_x",   &gen_vx);
  m_MC->SetBranchAddress("gen_y",   &gen_vy);
  m_MC->SetBranchAddress("gen_z",   &gen_vz);
  m_MC->SetBranchAddress("subpart_n",     &n_pu);
  m_MC->SetBranchAddress("subpart_pdgId", &pu_pdg);
  m_MC->SetBranchAddress("subpart_px",    &pu_px);
  m_MC->SetBranchAddress("subpart_py",    &pu_py);
  m_MC->SetBranchAddress("subpart_pz",    &pu_pz);
  m_MC->SetBranchAddress("subpart_eta",   &pu_eta);
  m_MC->SetBranchAddress("subpart_phi",   &pu_phi);
  m_MC->SetBranchAddress("subpart_x",     &pu_vx);
  m_MC->SetBranchAddress("subpart_y",     &pu_vy);
  m_MC->SetBranchAddress("subpart_z",     &pu_vz);


  int eventId = -1;
  int n_stubs =  0;
  std::vector<int>   *stub_clu1 = NULL;
  std::vector<int>   *stub_clu2 = NULL;
  std::vector<int>   *stub_layer = NULL;
  std::vector<int>   *stub_ladder  = NULL;
  std::vector<int>   *stub_module  = NULL;
  std::vector<int>   *stub_segment = NULL;
  std::vector<float> *stub_strip   = NULL;
  std::vector<float> *stub_x = NULL;
  std::vector<float> *stub_y = NULL;
  std::vector<float> *stub_z = NULL;
  std::vector<float> *stub_ds = NULL;
  std::vector<float> *stub_off = NULL;
  std::vector<float> *stub_phi0 = NULL;   ////////////////////////////////////////////
  std::vector<float> *stub_etaGEN = NULL; ////////////////////////////////////////////
  std::vector<float> *stub_pxGEN = NULL;  ////////////////////////////////////////////
  std::vector<float> *stub_pyGEN = NULL;  ////////////////////////////////////////////
  std::vector<float> *clu_x = NULL;      
  std::vector<float> *clu_y = NULL;
  std::vector<float> *clu_z = NULL;

  
  m_L1TT->SetBranchAddress("L1Tkevt",         &eventId); 
  m_L1TT->SetBranchAddress("L1TkSTUB_n",      &n_stubs);
  m_L1TT->SetBranchAddress("L1TkSTUB_clust1", &stub_clu1);
  m_L1TT->SetBranchAddress("L1TkSTUB_clust2", &stub_clu2);
  m_L1TT->SetBranchAddress("L1TkSTUB_layer",  &stub_layer);
  m_L1TT->SetBranchAddress("L1TkSTUB_ladder", &stub_ladder);
  m_L1TT->SetBranchAddress("L1TkSTUB_module", &stub_module);
  m_L1TT->SetBranchAddress("L1TkSTUB_seg",    &stub_segment);
  m_L1TT->SetBranchAddress("L1TkSTUB_strip",  &stub_strip);
  m_L1TT->SetBranchAddress("L1TkSTUB_x",      &stub_x);
  m_L1TT->SetBranchAddress("L1TkSTUB_y",      &stub_y);
  m_L1TT->SetBranchAddress("L1TkSTUB_z",      &stub_z);
  m_L1TT->SetBranchAddress("L1TkSTUB_deltas", &stub_ds);
  m_L1TT->SetBranchAddress("L1TkSTUB_cor",    &stub_off);
  m_L1TT->SetBranchAddress("L1TkSTUB_PHI0",   &stub_phi0);
  m_L1TT->SetBranchAddress("L1TkSTUB_etaGEN", &stub_etaGEN);
  m_L1TT->SetBranchAddress("L1TkSTUB_pxGEN",  &stub_pxGEN);
  m_L1TT->SetBranchAddress("L1TkSTUB_pyGEN",  &stub_pyGEN);
  m_L1TT->SetBranchAddress("L1TkCLUS_x",      &clu_x);
  m_L1TT->SetBranchAddress("L1TkCLUS_y",      &clu_y);
  m_L1TT->SetBranchAddress("L1TkCLUS_z",      &clu_z);
 
  int n_patterns = 0;                                 // Number of patterns in event eventtId
  std::vector< std::vector<int> > *patt_links = NULL; // Links to the stub ids of the pattern 
  std::vector<int>                *patt_secId = NULL; // Link to the sector ids of the pattern

  m_PATT->SetBranchAddress("L1PATT_n",     &n_patterns);
  m_PATT->SetBranchAddress("L1PATT_links", &patt_links);
  m_PATT->SetBranchAddress("L1PATT_secid", &patt_secId);
 

  // --- Initialize the Retina configuration parameters:
  retina_initialize(config);


  // =============================================================================================
  //  Loop over the tree entries
  // =============================================================================================

  set<int> PR_stubs;                      // stub indices of the event after PR
  map <int, set<int> > road_stubs;        // stub indices per road
  map <int, set<int> > road_stubs_layer;  // occupied layers per road
  map <int, set<int> > trigT_stubs;       // stub indices per trigger tower
  map <int, set<int> > trigT_stubs_layer; // occupied layers per trigger tower
  map <int, int> road_ncomb;

  map <int, vector<int> > v_road_stubs_layer;  // stub layers per road
  map <int, vector<int> > v_trigT_stubs_layer; // stub layers per trigger tower

  map <int, unsigned int> layer_mask;

  n_events = std::min(m_L1TT->GetEntries(),n_events);
  for (long long int ientry=0; ientry<n_events; ientry++) {

    m_L1TT->GetEntry(ientry);
    m_PATT->GetEntry(ientry);
    m_MC->GetEntry(ientry);

    int event_counter = ientry;

    // --- Skip events with no roads:
    if ( n_patterns == 0 ) continue;

    if ( ientry % 1000 == 0 )
      cout << "*** Processing event " << ientry << "/" << n_events << endl;

    if ( verboseLevel > 0 )
      std::cout << "Event = " <<  ientry 
		<< " ----------------------------------------------------------------------" 
		<< std::endl;

    for (int itow=0; itow<48; ++itow){
      Nroads_trigT[itow] = 0;
      Ncomb_trigT[itow]  = 0;
    }


    // ---  x_+ - x_- coverage for gen particles
    
    int n_part = 0;
    std::vector <GenPart_t> genP;
    for (int ipart=0; ipart<n_gen; ++ipart) {
	
      double y0 =  0.0217391304347826081; // 0.5/23.
      double y1 =  0.0046296296296296294; // 0.5/108.

      // --- Keep only charged stable particles:
      if ( fabs(gen_pdg->at(ipart)) !=   11 &&
	   fabs(gen_pdg->at(ipart)) !=   13 &&
	   fabs(gen_pdg->at(ipart)) !=  211 &&
	   fabs(gen_pdg->at(ipart)) !=  321 &&
	   fabs(gen_pdg->at(ipart)) != 2212 ) continue;
    
      double gen_charge = TMath::Sign(1., (Double_t) gen_pdg->at(ipart));
      if ( fabs(gen_pdg->at(ipart)) < 20. )
	gen_charge *= -1.; 

      double gen_pt  = sqrt(gen_px->at(ipart)*gen_px->at(ipart)+gen_py->at(ipart)*gen_py->at(ipart));
      double gen_phi = atan2(gen_py->at(ipart),gen_px->at(ipart));

      GenPart_t tmpPart;
      tmpPart.pdg = gen_pdg->at(ipart);
      tmpPart.pt  = gen_pt;
      tmpPart.phi = gen_phi;

      int phi_sector =  (int) (gen_phi/TMath::PiOver4());
      if ( gen_phi < 0. )
	phi_sector += 7;
      gen_phi += rot_angle[phi_sector];

      double gen_px_rot = gen_px->at(ipart)*cos(rot_angle[phi_sector])- 
	                  gen_py->at(ipart)*sin(rot_angle[phi_sector]);
      double gen_py_rot = gen_px->at(ipart)*sin(rot_angle[phi_sector])+
	                  gen_py->at(ipart)*cos(rot_angle[phi_sector]);
      double gen_vx_rot = gen_vx->at(ipart)*cos(rot_angle[phi_sector])-
                          gen_vy->at(ipart)*sin(rot_angle[phi_sector]);
      double gen_vy_rot = gen_vx->at(ipart)*sin(rot_angle[phi_sector])+
                          gen_vy->at(ipart)*cos(rot_angle[phi_sector]);

      double gen_theta = atan2(gen_pt,gen_pz->at(ipart));
      tmpPart.eta = -log(tan(0.5*gen_theta));
       
      if ( gen_pz->at(ipart) < 0. )
	gen_theta = atan2(gen_pt, -gen_pz->at(ipart));

      double gen_eta = -log(tan(0.5*gen_theta));
    
      // curvature and helix radius:
      double c = gen_charge*0.003*mMagneticField/gen_pt;
      double R = gen_pt/(0.003*mMagneticField);

      tmpPart.c = c;

      // helix center:
      double x_0 = gen_vx_rot - gen_charge*R*gen_py_rot/gen_pt;
      double y_0 = gen_vy_rot + gen_charge*R*gen_px_rot/gen_pt;
    
      // transverse and longitudinal impact parameters:
      double gen_d0 = gen_charge*(sqrt(x_0*x_0+y_0*y_0)-R);
      tmpPart.d0  = gen_d0;
      double diff = gen_vx_rot*gen_vx_rot+gen_vy_rot*gen_vy_rot-gen_d0*gen_d0;
      if ( diff < 0. ) diff = 0.;
      double gen_z0 = gen_vz->at(ipart) - 2./c*gen_pz->at(ipart)/gen_pt*asin(0.5*c*sqrt(diff));
      tmpPart.z0  = gen_z0;
      if (gen_pz->at(ipart)<0.)
	gen_z0 = -gen_vz->at(ipart) - 2./c*(-gen_pz->at(ipart))/gen_pt*asin(0.5*c*sqrt(diff));

      //if ( fabs(gen_pt)  < 1.5 ) continue; 
      if ( fabs(gen_eta) > 2. ) continue; 


      //h_gen_pdg   ->Fill(gen_pdg->at(ipart));
      //h_gen_pt    ->Fill(gen_pt);
      //h_gen_d0    ->Fill(gen_d0);
      //h_gen_phi   ->Fill(gen_phi);
      //h_gen_eta   ->Fill(gen_eta);
      //h_gen_theta ->Fill(gen_theta);
      //h_gen_z0    ->Fill(gen_z0);
      
      cout << " " << n_part << "  -  " << gen_pdg->at(ipart) << " "
	   << gen_pt << " " << gen_phi-rot_angle[phi_sector] << " " << gen_eta << " " << gen_z0 << endl;

      // -- XY view

      double p = -x_0/y_0;
      double q =  0.5/y_0;

      double x_plus  = 0.5*(y0+y1-2.*q)/p;
      double x_minus = 0.5*(y1-y0)/p; 

      tmpPart.x_P_XY = x_plus;
      tmpPart.x_M_XY = x_minus;

      if ( fabs(gen_eta) < 0.95 ) 
	h_max_XY_gen[0]->Fill(x_plus,x_minus);
      if ( fabs(gen_eta)>0.7 &&  fabs(gen_eta)<1.6 ) 
	h_max_XY_gen[1]->Fill(x_plus,x_minus);
      if ( fabs(gen_eta)>1.35 &&  fabs(gen_eta)<2.5 ) 
	h_max_XY_gen[0]->Fill(x_plus,x_minus);

      
      // -- RZ view

      y0 = 0.5/y0;
      y1 = 0.5/y1;

      p = tan(gen_theta);
      q = -p*gen_z0;

      x_plus  = 0.5*(y0+y1-2.*q)/p;
      x_minus = 0.5*(y1-y0)/p; 
      
      tmpPart.x_P_RZ = x_plus;
      tmpPart.x_M_RZ = x_minus;

      genP.push_back(tmpPart);

      if ( fabs(gen_eta) < 0.95 ) 
	h_max_RZ_gen[0]->Fill(x_plus,x_minus);
      if ( fabs(gen_eta)>0.7 &&  fabs(gen_eta)<1.6 ) 
	h_max_RZ_gen[1]->Fill(x_plus,x_minus);
      if ( fabs(gen_eta)>1.35 &&  fabs(gen_eta)<2.5 ) 
	h_max_RZ_gen[2]->Fill(x_plus,x_minus);

      n_part++;
    
    } // ipart loop


    // --- Get the stubs per road and trigger tower:

    // Clear the stub containers:
    PR_stubs.clear();
    road_stubs.clear();
    road_stubs_layer.clear();
    trigT_stubs.clear();
    trigT_stubs_layer.clear();
    road_ncomb.clear();

    v_road_stubs_layer.clear();
    v_trigT_stubs_layer.clear();

    // Loop through the roads to fill the stub containers:
    for (int iroad=0; iroad<n_patterns; ++iroad) {

      const int sector_id = patt_secId->at(iroad);

      Nroads_trigT[sector_id-1]++;
      
      // Loop through the stubs of the road:
      layer_mask[iroad] = 0x0;
      for (unsigned int istub=0; istub<patt_links->at(iroad).size(); ++istub ){
	
      	int stubRef = patt_links->at(iroad).at(istub);

	layer_mask[iroad] |= ( 0x1 << (stub_layer->at(stubRef)-5) );

	PR_stubs.insert(stubRef);
     	road_stubs[iroad].insert(stubRef);
     	road_stubs_layer[iroad].insert(stub_layer->at(stubRef));

      	trigT_stubs[sector_id].insert(stubRef);
     	trigT_stubs_layer[sector_id].insert(stub_layer->at(stubRef));

     	v_road_stubs_layer[iroad].push_back(stub_layer->at(stubRef));
     	v_trigT_stubs_layer[sector_id].push_back(stub_layer->at(stubRef));

      } // istub loop

    } // iroad loop


    h_Nstubs_tot->Fill(n_stubs);
    h_Nroads_tot->Fill(n_patterns);
    h_Nstubs_tot_PR->Fill(PR_stubs.size());

    for (map<int,set<int> >::iterator itow = trigT_stubs.begin(); itow != trigT_stubs.end(); ++itow){
      h_Nstubs_trigT[itow->first]->Fill((itow->second).size());

      for (set<int>::iterator istub=(itow->second).begin(); istub!=(itow->second).end(); ++istub)
	h_Nstubs_occ_PR->Fill(itow->first-1);

    }
    for (map<int,set<int> >::iterator istub = trigT_stubs_layer.begin(); istub != trigT_stubs_layer.end(); ++istub)
      h_Nstubs_trigT_layer->Fill((istub->second).size());
    for (map<int,set<int> >::iterator istub = road_stubs.begin(); istub != road_stubs.end(); ++istub)
      h_Nstubs_road->Fill((istub->second).size());
    for (map<int,set<int> >::iterator istub = road_stubs_layer.begin(); istub != road_stubs_layer.end(); ++istub)
      h_Nstubs_road_layer->Fill((istub->second).size());


    for (map<int,vector<int> >::iterator istub = v_road_stubs_layer.begin(); istub != v_road_stubs_layer.end(); ++istub){
      
      int n_comb_road = 1;

      for ( int ilayer=0; ilayer<16; ++ilayer ){
	int stub_count = std::count((istub->second).begin(), (istub->second).end(), active_layers[ilayer]);
	if (stub_count > 0){
	  h_NstubsLayer_road[ilayer]->Fill(stub_count);
	  n_comb_road *= stub_count;
	}      
      }
      

      road_ncomb[istub->first] =  n_comb_road;

      Ncomb_trigT[patt_secId->at(istub->first)-1] += n_comb_road;

    }
    

    for (map<int,vector<int> >::iterator istub = v_trigT_stubs_layer.begin(); istub != v_trigT_stubs_layer.end(); ++istub){
      
      for ( int ilayer=0; ilayer<16; ++ilayer ){
	int stub_count = std::count((istub->second).begin(), (istub->second).end(), active_layers[ilayer]);
	if (stub_count > 0)
	  h_NstubsLayer_trigT[ilayer]->Fill(stub_count);
      }

    }

    for (int itow=0; itow<48; ++itow) {
      //for (int itow=15; itow<31; ++itow) { ///////////////////////////////////////////////////////////////////

      //if ( Nroads_trigT[itow] > 0 ){
	h_Nroads_trigT[itow]->Fill(Nroads_trigT[itow]);
	h_Ncomb_trigT->Fill(Ncomb_trigT[itow]);
	//}
    }

    for (int iscan=0; iscan<n_scan_points; ++iscan) {

      if ( n_scan_points > 1 && verboseLevel > 0 ) {
	
	cout << endl;
	cout << " ================================================================================ " << endl;
	cout << "   Scan point #" << iscan << endl; 
	cout << " ================================================================================ " << endl;

      }

      // --- Apply the retina filters:
      map<int,set<int> >::iterator istub;
      map<int,set<int> >::iterator istub_begin =  trigT_stubs.begin();
      map<int,set<int> >::iterator istub_end   =  trigT_stubs.end();
      if ( fit_roads ) {
	istub_begin =  road_stubs.begin();
	istub_end   =  road_stubs.end();
      }

      for (istub = istub_begin; istub != istub_end; ++istub){

	int road_id = istub->first; // NB: When fitting per trigger tower, this is the tower id.

	// Constants used in X+-X- transformation:
	double y0 =  0.0217391304347826081; // 0.5/23.
	double y1 =  0.0046296296296296294; // 0.5/108.


	// --- Determine in which phi sector and eta range the trigger tower is:
	const int sector_id = ( fit_roads ? patt_secId->at(istub->first) : istub->first );
	const int phi_sector = sector_id % 8;
	const int eta_range  = sector_id / 8;

	int trigTow_type = 0;
	if ( eta_range==1 || eta_range==4 )
	  trigTow_type = 1;
	else if ( eta_range==0 || eta_range==5 )
	  trigTow_type = 2;

	
	bool is_6outof6 = ( layer_mask[road_id] == 0x3F );
	bool is_5outof6 = ( layer_mask[istub->first] == 0x1F ||
			    layer_mask[istub->first] == 0x2F ||
			    layer_mask[istub->first] == 0x37 ||
			    layer_mask[istub->first] == 0x3B ||
			    layer_mask[istub->first] == 0x3D ||
			    layer_mask[istub->first] == 0x3E );


	// Consider only central towers
	if ( sector_id<16 || sector_id>31 ) continue;

	cout << " sec = " <<  sector_id << "  road = " << road_id << " ----------------------------------------" << endl;

	// --- Fill the stubs vector:
	std::vector <Hit_t*> hits;
	for ( set<int>::iterator ihit=(istub->second).begin(); ihit!=(istub->second).end(); ++ihit ){

	  int stubRef = *ihit;

	  Hit_t* hit = new Hit_t();
	  hit->x     = stub_x->at(stubRef);
	  hit->y     = stub_y->at(stubRef);
	  // NB: We flip the z sign for the trigger towers with negative eta in
	  //     order to use the same retina setup as for eta > 0 
	  hit->z     = ( eta_range > 2 ? stub_z->at(stubRef) : -stub_z->at(stubRef) );
	  hit->rho   = sqrt( stub_x->at(stubRef)*stub_x->at(stubRef) +
			     stub_y->at(stubRef)*stub_y->at(stubRef) );
	  hit->layer = (short) stub_layer->at(stubRef);
	  hit->id    = stubRef;

	  double pitch = 0.01;
	  if ( hit->rho > 60. )
	    pitch = 0.009;

	  double dR = sqrt( (clu_x->at(stub_clu2->at(stubRef))-clu_x->at(stub_clu1->at(stubRef))) *
			    (clu_x->at(stub_clu2->at(stubRef))-clu_x->at(stub_clu1->at(stubRef))) +
			    (clu_y->at(stub_clu2->at(stubRef))-clu_y->at(stub_clu1->at(stubRef))) *
			    (clu_y->at(stub_clu2->at(stubRef))-clu_y->at(stub_clu1->at(stubRef))) );
	  hit->dphi = asin((stub_ds->at(stubRef)+stub_off->at(stubRef))*pitch/dR);
	  //hit->dphi = 0.;
     
	  hits.push_back(hit);
     	
	  if (verboseLevel==3 )
	    cout << std::distance((istub->second).begin(),ihit) << "  - " 
		 << " x = " << stub_x->at(stubRef)
		 << " y = " << stub_y->at(stubRef)
		 << " z = " << stub_z->at(stubRef)
		 << " --> R = " << hit->rho
		 << "  -  id = " <<  hit->id
		 << endl;

	} // ihit loop


	// --- Total combinations before the retina filtering:
	if ( is_6outof6 )
	  h_Ncomb_road_66->Fill(road_ncomb[road_id]);
	if ( is_5outof6 )
	  h_Ncomb_road_56->Fill(road_ncomb[road_id]);

 
	// ===========================================================================
	//  X-Y view
	// ===========================================================================

	// --- Phi sector rotation:
	rotateHits(hits, rot_angle[phi_sector]);

	// --- Conformal transformation: 
	confTrans(hits);

	
	// --------------------------------------------------------------------
	//  First find the retina maximum (to center the PR retina around it)
	// --------------------------------------------------------------------

	// --- Setup the retina:
	double pbins = config[trigTow_type]["xy_pbins_step1"];
	double qbins = config[trigTow_type]["xy_qbins_step1"];
	double pmin  = config[trigTow_type]["xy_pmin_step1"];
	double pmax  = config[trigTow_type]["xy_pmax_step1"];
	double qmin  = config[trigTow_type]["xy_qmin_step1"];
	double qmax  = config[trigTow_type]["xy_qmax_step1"];
  	
	double minWeight = config[trigTow_type]["xy_threshold_step1"];

	double pstep = (pmax-pmin)/pbins;
	double qstep = (qmax-qmin)/qbins;

	vector <double> sigma(8,sqrt(pstep*pstep+qstep*qstep));
	
	//  sigma scan
	//for (int i=0; i<sigma.size(); ++i){
	//  sigma[i] += 0.1*sqrt(pstep*pstep+qstep*qstep) * (iscan-5+1);
	//  //cout << i << " " << sigma[i] << endl;
	//}

	if ( config[trigTow_type]["xy_sigma1_step1"] != 0. ) 
	  sigma[0] = config[trigTow_type]["xy_sigma1_step1"];
	if ( config[trigTow_type]["xy_sigma2_step1"] != 0. ) 
	  sigma[1] = config[trigTow_type]["xy_sigma2_step1"];
	if ( config[trigTow_type]["xy_sigma3_step1"] != 0. ) 
	  sigma[2] = config[trigTow_type]["xy_sigma3_step1"];
	if ( config[trigTow_type]["xy_sigma4_step1"] != 0. ) 
	  sigma[3] = config[trigTow_type]["xy_sigma4_step1"];
	if ( config[trigTow_type]["xy_sigma5_step1"] != 0. ) 
	  sigma[4] = config[trigTow_type]["xy_sigma5_step1"];
	if ( config[trigTow_type]["xy_sigma6_step1"] != 0. ) 
	  sigma[5] = config[trigTow_type]["xy_sigma6_step1"];
	if ( config[trigTow_type]["xy_sigma7_step1"] != 0. ) 
	  sigma[6] = config[trigTow_type]["xy_sigma7_step1"];
	if ( config[trigTow_type]["xy_sigma8_step1"] != 0. ) 
	  sigma[7] = config[trigTow_type]["xy_sigma8_step1"];


	Retina retinaXY(hits, pbins+2, qbins+2, pmin-pstep, pmax+pstep, 
			qmin-qstep, qmax+qstep, sigma, minWeight, 1, XY);


	// --- Fill the retina and find maxima:
	retinaXY.fillGrid();
	retinaXY.findMaxima();
	if ( verboseLevel==6 )
	  retinaXY.dumpGrid(event_counter,1,road_id,iscan);
	if ( verboseLevel==5 )
	  retinaXY.printMaxima();

	// --- Get first step maxima:
	//vector <pqPoint> maximaXY = retinaXY.getMaxima();
	pqPoint bestMaxumum_XY = retinaXY.getBestPQ();


	// --------------------------------------------------------------------
	//  Apply now a new retina to reduce stub combinations
	// --------------------------------------------------------------------

	pbins = 8.;
	qbins = 8.;
	pmin  = bestMaxumum_XY.p*(1.-0.05);
	pmax  = bestMaxumum_XY.p*(1.+0.05);
	qmin  = bestMaxumum_XY.q*(1.+0.05);
	qmax  = bestMaxumum_XY.q*(1.-0.05);
  
	minWeight = 4.5;

	// Safety check:
	if ( pmax < pmin ){
	  double tmp = pmin;
	  pmin = pmax;
	  pmax = tmp;
	}
	if ( qmax < qmin ){
	  double tmp = qmin;
	  qmin = qmax;
	  qmax = tmp;
	}

	pstep = (pmax-pmin)/pbins;
	qstep = (qmax-qmin)/qbins;

	for (unsigned int i=0; i<sigma.size(); ++i)
	  sigma[i] = sqrt(pstep*pstep+qstep*qstep);


	Retina retinaXY_PR(hits, pbins+2, qbins+2, pmin-pstep, pmax+pstep, 
			   qmin-qstep, qmax+qstep, sigma, minWeight, 1, XY);

	// --- Fill the retina and find maxima:
	retinaXY_PR.fillGrid();
	retinaXY_PR.findMaxima();
	if ( verboseLevel==4 )
	  retinaXY_PR.dumpGrid(event_counter,1,road_id,iscan);
	if ( verboseLevel==3 )
	  retinaXY_PR.printMaxima();

	// --- Get first step maxima:
	vector <pqPoint> maximaXY_PR = retinaXY_PR.getMaxima();
	//pqPoint bestMaxumum_XY_PR = retinaXY_PR.getBestPQ();

	if ( is_6outof6 )
	  h_Nmax_road_66_XY->Fill(maximaXY_PR.size());
	if ( is_5outof6 )
	  h_Nmax_road_56_XY->Fill(maximaXY_PR.size());

	// --- Loop through the retina maxima:
	int n_comb_XY = 0;
	for (unsigned int imax_XY=0; imax_XY<maximaXY_PR.size(); ++imax_XY){

	  h_max_XY1[iscan][trigTow_type]->Fill( maximaXY_PR[imax_XY].p, maximaXY_PR[imax_XY].q);
	  h_wxy_1[iscan]->Fill(maximaXY_PR[imax_XY].w);

	  // --- Invert the X+-X- transformation:
	  double p = 0.5*(y1 - y0)/maximaXY_PR[imax_XY].q;
	  double q = y0 - p*(maximaXY_PR[imax_XY].p-maximaXY_PR[imax_XY].q);

	  // --- Associate stubs to this maximum:
	  vector <Hit_t*> hits_max;
	  vector <int> hits_max_layer;
	  for (unsigned int ihit=0; ihit<hits.size(); ++ihit){
	    
	    double dist   = (p*hits[ihit]->x-hits[ihit]->y+q)/p;
	    //double dist   = fabs(hits[ihit]->y-p*hits[ihit]->x-q)/sqrt(1.+p*p);
	    double weight = exp(-0.5*dist*dist/(sigma[0]*sigma[0]));
	    
	    if ( weight > 0.5 ){
	      hits_max.push_back(hits[ihit]);
	      hits_max_layer.push_back(hits[ihit]->layer);
	    }

	  } // ihit loop


	  int n_comb_road = 1;

	  for ( int ilayer=0; ilayer<16; ++ilayer ){
	    int stub_count = std::count(hits_max_layer.begin(), hits_max_layer.end(), active_layers[ilayer]);
	    if (stub_count > 0)
	      n_comb_road *= stub_count;
	  }      

	  //if (hits_max.size()==0)
	  //  n_comb_road = road_ncomb[road_id];

	  n_comb_XY += n_comb_road;

	} // imax_XY loop
	
	if (n_comb_XY==0)
	  n_comb_XY = road_ncomb[road_id];

	cout << "     XY:     Nmax = " << maximaXY_PR.size() << " ---->  " << n_comb_XY
	     << " (" << road_ncomb[road_id] << ")" << endl;

	if ( is_6outof6 )
	  h_Ncomb_road_66_XY->Fill(n_comb_XY);
	if ( is_5outof6 )
	  h_Ncomb_road_56_XY->Fill(n_comb_XY);


	
	// =========================================================================
	//  R-Z view
	// =========================================================================

	y0 = 0.5/y0;
	y1 = 0.5/y1;

	// --------------------------------------------------------------------
	//  First find the retina maximum (to center the PR retina around it)
	// --------------------------------------------------------------------

	// --- Setup the retina:
	pbins = config[trigTow_type]["rz_pbins_step1"];
	qbins = config[trigTow_type]["rz_qbins_step1"];
	pmin  = config[trigTow_type]["rz_pmin_step1"];
	pmax  = config[trigTow_type]["rz_pmax_step1"];
	qmin  = config[trigTow_type]["rz_qmin_step1"];
	qmax  = config[trigTow_type]["rz_qmax_step1"];
	
	minWeight = config[trigTow_type]["rz_threshold_step1"];

	pstep = (pmax-pmin)/pbins;
	qstep = (qmax-qmin)/qbins;

	for (unsigned int ilayer=0; ilayer<8; ++ilayer)
	  sigma[ilayer] = sqrt(pstep*pstep+qstep*qstep);

	if ( config[trigTow_type]["rz_sigma1_step1"] != 0. ) 
	  sigma[0] = config[trigTow_type]["rz_sigma1_step1"];
	if ( config[trigTow_type]["rz_sigma2_step1"] != 0. ) 
	  sigma[1] = config[trigTow_type]["rz_sigma2_step1"];
	if ( config[trigTow_type]["rz_sigma3_step1"] != 0. ) 
	  sigma[2] = config[trigTow_type]["rz_sigma3_step1"];
	if ( config[trigTow_type]["rz_sigma4_step1"] != 0. ) 
	  sigma[3] = config[trigTow_type]["rz_sigma4_step1"];
	if ( config[trigTow_type]["rz_sigma5_step1"] != 0. ) 
	  sigma[4] = config[trigTow_type]["rz_sigma5_step1"];
	if ( config[trigTow_type]["rz_sigma6_step1"] != 0. ) 
	  sigma[5] = config[trigTow_type]["rz_sigma6_step1"];
	if ( config[trigTow_type]["rz_sigma7_step1"] != 0. ) 
	  sigma[6] = config[trigTow_type]["rz_sigma7_step1"];
	if ( config[trigTow_type]["rz_sigma8_step1"] != 0. ) 
	  sigma[7] = config[trigTow_type]["rz_sigma8_step1"];

	
	Retina retinaRZ(hits, pbins+2, qbins+2, pmin-pstep, pmax+pstep, 
			qmin-qstep, qmax+qstep,	sigma, minWeight, 1, RZ);


	// --- Fill the retina and find maxima:
	retinaRZ.fillGrid();
	retinaRZ.findMaxima();
	if ( verboseLevel==6 )
	  retinaRZ.dumpGrid(event_counter,1,road_id,iscan);
	if ( verboseLevel==5 )
	  retinaRZ.printMaxima();

	// --- Get first step maximum:
	//vector <pqPoint> maximaRZ = retinaRZ.getMaxima();
	pqPoint bestMaxumum_RZ = retinaRZ.getBestPQ();


	// --------------------------------------------------------------------
	//  Apply now a new retina to reduce stub combinations
	// --------------------------------------------------------------------

	pbins = 8.;
	qbins = 8.;
	pmin  = bestMaxumum_RZ.p*(1.-0.1);
	pmax  = bestMaxumum_RZ.p*(1.+0.1);
	qmin  = bestMaxumum_RZ.q*(1.-0.1);
	qmax  = bestMaxumum_RZ.q*(1.+0.1);

	minWeight = 4.;

	// Safety check:
	if ( pmax < pmin ){
	  double tmp = pmin;
	  pmin = pmax;
	  pmax = tmp;
	}
	if ( qmax < qmin ){
	  double tmp = qmin;
	  qmin = qmax;
	  qmax = tmp;
	}

	pstep = (pmax-pmin)/pbins;
	qstep = (qmax-qmin)/qbins;

	for (unsigned int i=0; i<sigma.size(); ++i)
	  sigma[i] = sqrt(pstep*pstep+qstep*qstep);

	
	Retina retinaRZ_PR(hits, pbins+2, qbins+2, pmin-pstep, pmax+pstep, 
			   qmin-qstep, qmax+qstep, sigma, minWeight, 1, RZ);

	// --- Fill the retina and find maxima:
	retinaRZ_PR.fillGrid();
	retinaRZ_PR.findMaxima();
	if ( verboseLevel==4 )
	  retinaRZ_PR.dumpGrid(event_counter,1,road_id,iscan);
	if ( verboseLevel==3 )
	  retinaRZ_PR.printMaxima();


	// --- Get first step maxima:
	vector <pqPoint> maximaRZ_PR = retinaRZ_PR.getMaxima();

	if ( is_6outof6 )
	  h_Nmax_road_66_RZ->Fill(maximaRZ_PR.size());
	if ( is_5outof6 )
	  h_Nmax_road_56_RZ->Fill(maximaRZ_PR.size());


	// --- Loop through the retina maxima:
	int n_comb_RZ = 0;
	for (unsigned int imax_RZ=0; imax_RZ<maximaRZ_PR.size(); ++imax_RZ){
	    
	  h_max_RZ1[iscan][trigTow_type]->Fill( maximaRZ_PR[imax_RZ].p, maximaRZ_PR[imax_RZ].q);
	  h_wrz_1[iscan]->Fill(maximaRZ_PR[imax_RZ].w);

	  // --- Invert the X+-X- transformation:
	  double p = 0.5*(y1 - y0)/maximaRZ_PR[imax_RZ].q;
	  double q = y0 - p*(maximaRZ_PR[imax_RZ].p-maximaRZ_PR[imax_RZ].q);

	  // --- Associate stubs to this maximum:
	  vector <Hit_t*> hits_max;
	  vector <int> hits_max_layer;
	  for (unsigned int ihit=0; ihit<hits.size(); ++ihit){
	    
	    double dist   = (p*hits[ihit]->x-hits[ihit]->y+q)/p;
	    //double dist   = fabs(hits[ihit]->y-p*hits[ihit]->x-q)/sqrt(1.+p*p);
	    double weight = exp(-0.5*dist*dist/(sigma[0]*sigma[0]));
	    
	    if ( weight > 0.5 ){
	      hits_max.push_back(hits[ihit]);
	      hits_max_layer.push_back(hits[ihit]->layer);
	    }

	  } // ihit loop


	  int n_comb_road = 1;

	  for ( int ilayer=0; ilayer<16; ++ilayer ){
	    int stub_count = std::count(hits_max_layer.begin(), hits_max_layer.end(), active_layers[ilayer]);
	    if (stub_count > 0)
	      n_comb_road *= stub_count;
	  }      

	  //if (hits_max.size()==0)
	  //  n_comb_road = road_ncomb[road_id];

	  //cout << "     RZ:  road " << road_id << " ---->  " << n_comb_road 
	  //     << " (" << road_ncomb[road_id] << ")" << endl;

	  n_comb_RZ += n_comb_road;

	} // imax_RZ loop

	if (n_comb_RZ==0)
	  n_comb_RZ = road_ncomb[road_id];

	cout << "     RZ:     Nmax = " << maximaRZ_PR.size() << " ---->  " << n_comb_RZ
	     << " (" << road_ncomb[road_id] << ")" << endl;

	if ( is_6outof6 )
	  h_Ncomb_road_66_RZ->Fill(n_comb_RZ);
	if ( is_5outof6 )
	  h_Ncomb_road_56_RZ->Fill(n_comb_RZ);

	  
	// =========================================================================
	//  Phi-Eta view
	// =========================================================================

	// --- Invert the conformal transformation: 
	//confTrans(hits);

	// --------------------------------------------------------------------
	//  First find the retina maximum (to center the PR retina around it)
	// --------------------------------------------------------------------

	// --- Setup the retina:
	double phi_bins = 50;
	double eta_bins = 50;
	double phi_min  = 0.15;
	double phi_max  = 1.25;
	double eta_min  = 0.;
	double eta_max  = 1.;

	minWeight = 2.;


	double phi_step = (phi_max-phi_min)/phi_bins;
	double eta_step = (eta_max-eta_min)/eta_bins;

	for (unsigned int ilayer=0; ilayer<8; ++ilayer){
	  sigma[ilayer] = 2.*sqrt(phi_step*phi_step+eta_step*eta_step);
	  //if (ilayer>2)
	  //  sigma[ilayer] = sqrt(phi_step*phi_step+10.*eta_step*eta_step);
	}

	Retina retinaPhiEta(hits, phi_bins+2, eta_bins+2, 
			    phi_min-phi_step, phi_max+phi_step, 
			    eta_min-eta_step, eta_max+eta_step, 
			    sigma, minWeight, 2, PhiEta);

	retinaPhiEta.fillGrid();
	retinaPhiEta.findMaxima();
	if ( verboseLevel==6 )
	  retinaPhiEta.dumpGrid(event_counter,1,road_id,iscan);
	if ( verboseLevel==5 )
	  retinaPhiEta.printMaxima();

	
	// --- Get first step maximum:
	//vector <pqPoint> maximaPhiEta = retinaPhiEta.getMaxima();
	pqPoint bestMaxumum_PhiEta = retinaPhiEta.getBestPQ();


	// --------------------------------------------------------------------
	//  Apply now a new retina to reduce stub combinations
	// --------------------------------------------------------------------


	phi_bins = 8.;
	eta_bins = 8.;
	phi_min  = bestMaxumum_PhiEta.p-3.*0.015;
	phi_max  = bestMaxumum_PhiEta.p+3.*0.015;
	eta_min  = bestMaxumum_PhiEta.q-3.*0.04;
	eta_max  = bestMaxumum_PhiEta.q+3.*0.04;

	minWeight = 3.;

	// Safety check:
	if ( phi_max < phi_min ){
	  double tmp = phi_min;
	  phi_min = phi_max;
	  phi_max = tmp;
	}
	if ( eta_max < eta_min ){
	  double tmp = eta_min;
	  eta_min = eta_max;
	  eta_max = tmp;
	}


	phi_step = (phi_max-phi_min)/phi_bins;
	eta_step = (eta_max-eta_min)/eta_bins;
	
	for (unsigned int ilayer=0; ilayer<8; ++ilayer){
	  sigma[ilayer] = 1.5*sqrt(phi_step*phi_step+eta_step*eta_step);
	  //if (ilayer>2)
	  //  sigma[ilayer] = sqrt(phi_step*phi_step+10.*eta_step*eta_step);
	}

	Retina retinaPhiEta_PR(hits, phi_bins+2, eta_bins+2, 
			       phi_min-phi_step, phi_max+phi_step, 
			       eta_min-eta_step, eta_max+eta_step, 
			       sigma, minWeight, 2, PhiEta);

	retinaPhiEta_PR.fillGrid();
	retinaPhiEta_PR.findMaxima();
	if ( verboseLevel==4 )
	  retinaPhiEta_PR.dumpGrid(event_counter,1,road_id,iscan);
	if ( verboseLevel==3 )
	  retinaPhiEta_PR.printMaxima();


	vector <pqPoint> maximaPhiEta_PR = retinaPhiEta_PR.getMaxima();

	if ( is_6outof6 )
	  h_Nmax_road_66_PE->Fill(maximaPhiEta_PR.size());
	if ( is_5outof6 )
	  h_Nmax_road_56_PE->Fill(maximaPhiEta_PR.size());


	// --- Loop through the retina maxima:
	int n_comb_PhiEta = 0;
	for (unsigned int imax_PhiEta=0; imax_PhiEta<maximaPhiEta_PR.size(); ++imax_PhiEta){
	    
	  //h_max_PE1[iscan][trigTow_type]->Fill( maximaPhiEta[imax_PhiEta].p, maximaPhiEta[imax_PhiEta].q);
	  //h_wpe_1[iscan]->Fill(maximaPhiEta[imax_PhiEta].w);

	  // --- Associate stubs to this maximum:
	  vector <Hit_t*> hits_max;
	  vector <int> hits_max_layer;
	  for (unsigned int ihit=0; ihit<hits.size(); ++ihit){
	    
	    double phi_stub   = atan2(hits[ihit]->y,hits[ihit]->x) + hits[ihit]->dphi;
	    double theta_stub = atan2(hits[ihit]->rho,hits[ihit]->z); 
	    double eta_stub   = -log(tan(0.5*theta_stub));

	    double delta_phi = phi_stub - maximaPhiEta_PR[imax_PhiEta].p;
	    double delta_eta = eta_stub - maximaPhiEta_PR[imax_PhiEta].q;

	    double dist2 = delta_phi*delta_phi + delta_eta*delta_eta;

	    double weight = exp(-0.5*dist2/(sigma[0]*sigma[0]));
	    
	    if ( weight > 0.5 ){
	      hits_max.push_back(hits[ihit]);
	      hits_max_layer.push_back(hits[ihit]->layer);
	    }

	  } // ihit loop


	  int n_comb_road = 1;

	  for ( int ilayer=0; ilayer<16; ++ilayer ){
	    int stub_count = std::count(hits_max_layer.begin(), hits_max_layer.end(), active_layers[ilayer]);
	    if (stub_count > 0)
	      n_comb_road *= stub_count;
	  }      

	  //if (hits_max.size()==0)
	  //  n_comb_road = road_ncomb[road_id];

	  //cout << "     PhiEta:  road " << road_id << " ---->  " << n_comb_road 
	  //     << " (" << road_ncomb[road_id] << ")" << endl;

	  n_comb_PhiEta += n_comb_road;

	} // imax_PhiEta loop

	if (n_comb_PhiEta==0)
	  n_comb_PhiEta = road_ncomb[road_id];

	cout << "     PhiEta: Nmax = " << maximaPhiEta_PR.size() << " ---->  " << n_comb_PhiEta
	     << " (" << road_ncomb[road_id] << ")" << endl;

	if ( is_6outof6 )
	  h_Ncomb_road_66_PE->Fill(n_comb_PhiEta);
	if ( is_5outof6 )
	  h_Ncomb_road_56_PE->Fill(n_comb_PhiEta);

      } // istub loop

    } // iscan loop

  } // ientry loop


  // --- Save histograms:
  TFile f("RetinaPRStandalone_histos.root","recreate");

  h_Nstubs_tot->Write();  
  h_Nroads_tot->Write();  
  h_Nstubs_tot_PR->Write();  
  h_Nstubs_occ_PR->Write();  
  h_Nstubs_road->Write();  
  h_Nstubs_road_layer->Write();  
  h_Nstubs_trigT_layer->Write();
  h_Ncomb_road_66->Write();
  h_Ncomb_road_66_XY->Write();
  h_Ncomb_road_66_RZ->Write();
  h_Ncomb_road_66_PE->Write();
  h_Ncomb_road_56->Write();
  h_Ncomb_road_56_XY->Write();
  h_Ncomb_road_56_RZ->Write();
  h_Ncomb_road_56_PE->Write();
  h_Nmax_road_66_XY->Write();
  h_Nmax_road_66_RZ->Write();
  h_Nmax_road_66_PE->Write();
  h_Nmax_road_56_XY->Write();
  h_Nmax_road_56_RZ->Write();
  h_Nmax_road_56_PE->Write();
  h_Ncomb_trigT->Write();
  h_trigT_occ->Write();

  for (int itow=0; itow<48; ++itow){
    h_Nstubs_trigT[itow]->Write();
    h_Nroads_trigT[itow]->Write();
  }
  
  for (int ilayer=0; ilayer<16; ++ilayer){
    
    h_NstubsLayer_road[ilayer]->Write();
    h_NstubsLayer_trigT[ilayer]->Write();
  
  }
  
  
  h_gen_pdg->Write();
  h_gen_pt->Write();
  h_gen_d0->Write();
  h_gen_phi->Write();
  h_gen_eta->Write();
  h_gen_theta->Write();
  h_gen_z0->Write();
  
  for (int ihist=0; ihist<3; ++ihist){
    h_max_XY_gen[ihist]->Write();
    h_max_RZ_gen[ihist]->Write();
  }
  
  for (int iscan=0; iscan<n_scan_points; ++iscan){
  
    for (int ihist=0; ihist<3; ++ihist){
      h_max_XY1[iscan][ihist]->Write();
      h_max_XY2[iscan][ihist]->Write();
      h_max_RZ1[iscan][ihist]->Write();
      h_max_RZ2[iscan][ihist]->Write();
    }
  
    h_Ntrks[iscan]->Write();
    h_trk_c[iscan]->Write();
    h_trk_pt[iscan]->Write();
    h_trk_phi[iscan]->Write();
    h_trk_eta[iscan]->Write();
    h_trk_z0[iscan]->Write();
  
    h_wxy_1[iscan]->Write();
    h_wxy_2[iscan]->Write();
    h_wrz_1[iscan]->Write();
    h_wrz_2[iscan]->Write();
  
    h_res_pt_rel[iscan]->Write();
    h_res_phi[iscan]->Write();
    h_res_eta[iscan]->Write();
    h_res_z0[iscan]->Write();
    
  }

  f.Close();


  // --- Clean-up the heap:
  delete m_L1TT;
  delete m_PATT;

  return EXIT_SUCCESS;

}

// ===============================================================================================
void rotateHits(vector<Hit_t*> hits, double angle){
  
  for (unsigned int ihit=0; ihit<hits.size(); ihit++) {
    double x = hits[ihit]->x*cos(angle) - hits[ihit]->y*sin(angle);
    double y = hits[ihit]->x*sin(angle) + hits[ihit]->y*cos(angle);
    hits[ihit]->x = x;
    hits[ihit]->y = y;
  }
  
}

// ===============================================================================================
void confTrans(vector<Hit_t*> hits){
  
  for (unsigned int ihit=0; ihit<hits.size(); ihit++) {
    double R2 = hits[ihit]->x*hits[ihit]->x + hits[ihit]->y*hits[ihit]->y;
    hits[ihit]->x /= R2;
    hits[ihit]->y /= R2;
  }
  
}

// ===============================================================================================
void retina_initialize(std::map <const char*,double> *config){

  // Enter all the retina parameters
  // (we refer to the detector geometry in
  //  http://sviret.web.cern.ch/sviret/Images/CMS/Upgrade/Eta6_Phi8.jpg)

  // --- Central trigger tower:
  config[0]["xy_pbins_step1"]     = 100.;
  config[0]["xy_qbins_step1"]     = 100.;
  config[0]["xy_pmin_step1"]      =  0.;
  config[0]["xy_pmax_step1"]      =  0.05;
  config[0]["xy_qmin_step1"]      = -0.05;
  config[0]["xy_qmax_step1"]      =  0.;
  config[0]["xy_threshold_step1"] =  4.5;
  config[0]["xy_sigma1_step1"]    =  0.;
  config[0]["xy_sigma2_step1"]    =  0.;
  config[0]["xy_sigma3_step1"]    =  0.;
  config[0]["xy_sigma4_step1"]    =  0.;
  config[0]["xy_sigma5_step1"]    =  0.;
  config[0]["xy_sigma6_step1"]    =  0.;
  config[0]["xy_sigma7_step1"]    =  0.;
  config[0]["xy_sigma8_step1"]    =  0.;
  config[0]["xy_pbins_step2"]     = 100.;
  config[0]["xy_qbins_step2"]     = 100.;
  config[0]["xy_zoom_step2"]      = 1.;
  config[0]["xy_threshold_step2"] =  4.5;
  config[0]["xy_sigma1_step2"]    =  0.;
  config[0]["xy_sigma2_step2"]    =  0.;
  config[0]["xy_sigma3_step2"]    =  0.;
  config[0]["xy_sigma4_step2"]    =  0.;
  config[0]["xy_sigma5_step2"]    =  0.;
  config[0]["xy_sigma6_step2"]    =  0.;
  config[0]["xy_sigma7_step2"]    =  0.;
  config[0]["xy_sigma8_step2"]    =  0.;

  config[0]["rz_pbins_step1"]     =  100.;
  config[0]["rz_qbins_step1"]     =  100.;
  config[0]["rz_pmin_step1"]      = -15.;
  config[0]["rz_pmax_step1"]      =  85.;
  config[0]["rz_qmin_step1"]      =  -2.;
  config[0]["rz_qmax_step1"]      =  46.;
  config[0]["rz_threshold_step1"] =  4.5;
  config[0]["rz_sigma1_step1"]    =  0.;
  config[0]["rz_sigma2_step1"]    =  0.;
  config[0]["rz_sigma3_step1"]    =  0.;
  config[0]["rz_sigma4_step1"]    =  0.;
  config[0]["rz_sigma5_step1"]    =  0.;
  config[0]["rz_sigma6_step1"]    =  0.;
  config[0]["rz_sigma7_step1"]    =  0.;
  config[0]["rz_sigma8_step1"]    =  0.;
  config[0]["rz_pbins_step2"]     = 80.;
  config[0]["rz_qbins_step2"]     = 80.;
  config[0]["rz_zoom_step2"]      = 1.5;
  config[0]["rz_threshold_step2"] =  4.;
  config[0]["rz_sigma1_step2"]    =  0.;
  config[0]["rz_sigma2_step2"]    =  0.;
  config[0]["rz_sigma3_step2"]    =  0.;
  config[0]["rz_sigma4_step2"]    =  0.;
  config[0]["rz_sigma5_step2"]    =  0.;
  config[0]["rz_sigma6_step2"]    =  0.;
  config[0]["rz_sigma7_step2"]    =  0.;
  config[0]["rz_sigma8_step2"]    =  0.;

  // --- Hybrid trigger tower:
  config[1]["xy_pbins_step1"]     = config[0]["xy_pbins_step1"];
  config[1]["xy_qbins_step1"]     = config[0]["xy_qbins_step1"];
  config[1]["xy_pmin_step1"]      = config[0]["xy_pmin_step1"];
  config[1]["xy_pmax_step1"]      = config[0]["xy_pmax_step1"];
  config[1]["xy_qmin_step1"]      = config[0]["xy_qmin_step1"];
  config[1]["xy_qmax_step1"]      = config[0]["xy_qmax_step1"];
  config[1]["xy_threshold_step1"] = config[0]["xy_threshold_step1"];
  config[1]["xy_sigma1_step1"]    = config[0]["xy_sigma1_step1"];
  config[1]["xy_sigma2_step1"]    = config[0]["xy_sigma2_step1"];
  config[1]["xy_sigma3_step1"]    = config[0]["xy_sigma3_step1"];
  config[1]["xy_sigma4_step1"]    = config[0]["xy_sigma4_step1"];
  config[1]["xy_sigma5_step1"]    = config[0]["xy_sigma5_step1"];
  config[1]["xy_sigma6_step1"]    = config[0]["xy_sigma6_step1"];
  config[1]["xy_sigma7_step1"]    = config[0]["xy_sigma7_step1"];
  config[1]["xy_sigma8_step1"]    = config[0]["xy_sigma8_step1"];
  config[1]["xy_pbins_step2"]     = config[0]["xy_pbins_step2"];
  config[1]["xy_qbins_step2"]     = config[0]["xy_qbins_step2"];
  config[1]["xy_zoom_step2"]      = config[0]["xy_zoom_step2"];
  config[1]["xy_threshold_step2"] = config[0]["xy_threshold_step2"];
  config[1]["xy_sigma1_step2"]    = config[0]["xy_sigma1_step2"];
  config[1]["xy_sigma2_step2"]    = config[0]["xy_sigma2_step2"];
  config[1]["xy_sigma3_step2"]    = config[0]["xy_sigma3_step2"];
  config[1]["xy_sigma4_step2"]    = config[0]["xy_sigma4_step2"];
  config[1]["xy_sigma5_step2"]    = config[0]["xy_sigma5_step2"];
  config[1]["xy_sigma6_step2"]    = config[0]["xy_sigma6_step2"];
  config[1]["xy_sigma7_step2"]    = config[0]["xy_sigma7_step2"];
  config[1]["xy_sigma8_step2"]    = config[0]["xy_sigma8_step2"];

  config[1]["rz_pbins_step1"]     =  100.;
  config[1]["rz_qbins_step1"]     =  100.;
  config[1]["rz_pmin_step1"]      =  35.;
  config[1]["rz_pmax_step1"]      = 170.;
  config[1]["rz_qmin_step1"]      =  30.;
  config[1]["rz_qmax_step1"]      = 105.;
  config[1]["rz_threshold_step1"] =  4.;
  config[1]["rz_sigma1_step1"]    =  0.;
  config[1]["rz_sigma2_step1"]    =  0.;
  config[1]["rz_sigma3_step1"]    =  0.;
  config[1]["rz_sigma4_step1"]    =  0.;
  config[1]["rz_sigma5_step1"]    =  0.;
  config[1]["rz_sigma6_step1"]    =  0.;
  config[1]["rz_sigma7_step1"]    =  0.;
  config[1]["rz_sigma8_step1"]    =  0.;
  config[1]["rz_pbins_step2"]     = 80.;
  config[1]["rz_qbins_step2"]     = 80.;
  config[1]["rz_zoom_step2"]      = 1.5;
  config[1]["rz_threshold_step2"] =  3.;
  config[1]["rz_sigma1_step2"]    =  0.;
  config[1]["rz_sigma2_step2"]    =  0.;
  config[1]["rz_sigma3_step2"]    =  0.;
  config[1]["rz_sigma4_step2"]    =  0.;
  config[1]["rz_sigma5_step2"]    =  0.;
  config[1]["rz_sigma6_step2"]    =  0.;
  config[1]["rz_sigma7_step2"]    =  0.;
  config[1]["rz_sigma8_step2"]    =  0.;

  // --- Forward trigger tower:
  config[2]["xy_pbins_step1"]     = config[0]["xy_pbins_step1"];
  config[2]["xy_qbins_step1"]     = config[0]["xy_qbins_step1"];
  config[2]["xy_pmin_step1"]      = config[0]["xy_pmin_step1"];
  config[2]["xy_pmax_step1"]      = config[0]["xy_pmax_step1"];
  config[2]["xy_qmin_step1"]      = config[0]["xy_qmin_step1"];
  config[2]["xy_qmax_step1"]      = config[0]["xy_qmax_step1"];
  config[2]["xy_threshold_step1"] = config[0]["xy_threshold_step1"];
  config[2]["xy_sigma1_step1"]    = config[0]["xy_sigma1_step1"];
  config[2]["xy_sigma2_step1"]    = config[0]["xy_sigma2_step1"];
  config[2]["xy_sigma3_step1"]    = config[0]["xy_sigma3_step1"];
  config[2]["xy_sigma4_step1"]    = config[0]["xy_sigma4_step1"];
  config[2]["xy_sigma5_step1"]    = config[0]["xy_sigma5_step1"];
  config[2]["xy_sigma6_step1"]    = config[0]["xy_sigma6_step1"];
  config[2]["xy_sigma7_step1"]    = config[0]["xy_sigma7_step1"];
  config[2]["xy_sigma8_step1"]    = config[0]["xy_sigma8_step1"];
  config[2]["xy_pbins_step2"]     = config[0]["xy_pbins_step2"];
  config[2]["xy_qbins_step2"]     = config[0]["xy_qbins_step2"];
  config[2]["xy_zoom_step2"]      = config[0]["xy_zoom_step2"];
  config[2]["xy_threshold_step2"] = config[0]["xy_threshold_step2"];
  config[2]["xy_sigma1_step2"]    = config[0]["xy_sigma1_step2"];
  config[2]["xy_sigma2_step2"]    = config[0]["xy_sigma2_step2"];
  config[2]["xy_sigma3_step2"]    = config[0]["xy_sigma3_step2"];
  config[2]["xy_sigma4_step2"]    = config[0]["xy_sigma4_step2"];
  config[2]["xy_sigma5_step2"]    = config[0]["xy_sigma5_step2"];
  config[2]["xy_sigma6_step2"]    = config[0]["xy_sigma6_step2"];
  config[2]["xy_sigma7_step2"]    = config[0]["xy_sigma7_step2"];
  config[2]["xy_sigma8_step2"]    = config[0]["xy_sigma8_step2"];

  config[2]["rz_pbins_step1"]     = 100.;
  config[2]["rz_qbins_step1"]     = 100.;
  config[2]["rz_pmin_step1"]      = 100.;
  config[2]["rz_pmax_step1"]      = 410.;
  config[2]["rz_qmin_step1"]      =  75.;
  config[2]["rz_qmax_step1"]      = 260.;
  config[2]["rz_threshold_step1"] =  4.;
  config[2]["rz_sigma1_step1"]    =  1.5;
  config[2]["rz_sigma2_step1"]    =  1.5;
  config[2]["rz_sigma3_step1"]    =  1.5;
  config[2]["rz_sigma4_step1"]    =  1.5;
  config[2]["rz_sigma5_step1"]    =  1.5;
  config[2]["rz_sigma6_step1"]    =  1.5;
  config[2]["rz_sigma7_step1"]    =  1.5;
  config[2]["rz_sigma8_step1"]    =  1.5;
  config[2]["rz_pbins_step2"]     = 80.;
  config[2]["rz_qbins_step2"]     = 80.;
  config[2]["rz_zoom_step2"]      = 1.5;
  config[2]["rz_threshold_step2"] = 3.;
  config[2]["rz_sigma1_step2"]    = 0.;
  config[2]["rz_sigma2_step2"]    = 0.;
  config[2]["rz_sigma3_step2"]    = 0.;
  config[2]["rz_sigma4_step2"]    = 0.;
  config[2]["rz_sigma5_step2"]    = 0.;
  config[2]["rz_sigma6_step2"]    = 0.;
  config[2]["rz_sigma7_step2"]    = 0.;
  config[2]["rz_sigma8_step2"]    = 0.;


}
