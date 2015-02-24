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


// --- Function paradigms:
void rotateHits(vector<Hit_t*> hits, double angle);
void confTrans(vector<Hit_t*> hits);
void retina_initialize(std::map <const char*,double> *config);


// --- Global variables:
const double mMagneticField = 3.8;


// --- Retina configuration parameters:
std::map <const char*,double> config[3];

const double rot_angle[8] = {  0.39269908169872414,   //  pi/8
			      -0.39269908169872414,   // -pi/8
			      -1.17809724509617242,   // -3/8 pi
			      -1.96349540849362070,   // -5/8 pi
			      -2.74889357189106898,   // -7/8 pi
			      -3.53429173528851726,   // -9/8 pi
			      -4.31968989868596509,   // -11/8 pi
			      -5.10508806208341426 }; // -13/8 pi


// --- Histograms:
TH1F * h_Nstubs_tot;  
TH1F * h_Nstubs_tot_PR;  
TH1F * h_Nstubs_road;  
TH1F * h_Nstubs_layer;  
TH1F * h_Nstubs_trigT;  
TH2F * h_max_XY1[3];  
TH2F * h_max_XY2[3];  
TH2F * h_max_RZ1[3];  
TH2F * h_max_RZ2[3];  


int main(int argc, char* argv[]) {

  if ( argc<2 ){
    cout << "\n*** No arguments provided! Run \"./Retina -h\" for a list of available options.\n" << endl;
    return EXIT_SUCCESS;
  }

  long long int n_events = 1;
  bool fit_roads = false;
  string input_file = "";
  int verboseLevel = 0;

  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i],"--fit-roads")) {
      fit_roads = true;
      cout << "Fitting single roads." << endl;
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
	   << "-v \t verbosity level \t\t default = 0" << endl
	   << "-n # \t set number of events \t\t\t default = 1" << endl 
	   << "-f <filename> \t a root file or a text file containing a list of root files: \t default = \"\"" 
	   << endl
	   << "--fit-roads \t fit per road: \t default = false" << endl;
      return EXIT_SUCCESS;
    }
  }


  // =============================================================================================
  //  Book the histograms:
  // =============================================================================================

  h_Nstubs_tot    = new TH1F("h_Nstubs_tot","Total stubs per event",1000,0,20000);  
  h_Nstubs_tot_PR = new TH1F("h_Nstubs_tot_PR","Total stubs per event after PR",1000,0,20000);  
  h_Nstubs_trigT  = new TH1F("h_Nstubs_trigT","Total stubs per trigger tower",400,0,400);  
  h_Nstubs_road   = new TH1F("h_Nstubs_road","Total stubs per road",30,0,30);  
  h_Nstubs_layer  = new TH1F("h_Nstubs_layer","Different layers per road",30,0,30);  

  for (int ihist=0; ihist<3; ++ihist){
    TString hname = Form("h_max_XY1_%d",ihist);
    h_max_XY1[ihist] = new TH2F(hname.Data(),"XY-1 maxima",100,-0.05,0.05,100,-0.05,0.05);
    hname = Form("h_max_XY2_%d",ihist);
    h_max_XY2[ihist] = new TH2F(hname.Data(),"XY-2 maxima",100,-0.05,0.05,100,-0.05,0.05);
    hname = Form("h_max_RZ1_%d",ihist);
    h_max_RZ1[ihist] = new TH2F(hname.Data(),"RZ-1 maxima",100,-20.,240.,100,-60.,180);
    hname = Form("h_max_RZ2_%d",ihist);
    h_max_RZ2[ihist] = new TH2F(hname.Data(),"RZ-2 maxima",100,-20.,240.,100,-60.,180);
  }

  // =============================================================================================
  //  Get the trees from the root files
  // =============================================================================================

  TChain *m_L1TT = new TChain("TkStubs");   // Tree containing the stub info 
  TChain *m_PATT = new TChain("L1tracks");  // Tree containing the L1 track reco info
  TChain *m_MC   = new TChain("MC");        // Tree containing the true MC info 

  std::size_t found = input_file.find(".root");
  // Case 1: it's a root file
  if ( found != std::string::npos ) {
    m_L1TT->Add(input_file.c_str());
    m_PATT->Add(input_file.c_str());
    m_MC->Add(input_file.c_str());
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
	m_L1TT->Add(STRING.c_str());   
	m_PATT->Add(STRING.c_str());   
	m_MC->Add(STRING.c_str());   
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

  m_MC->SetBranchAddress("gen_n",   &n_gen);
  m_MC->SetBranchAddress("gen_pdg", &gen_pdg);
  m_MC->SetBranchAddress("gen_px",  &gen_px);
  m_MC->SetBranchAddress("gen_py",  &gen_py);
  m_MC->SetBranchAddress("gen_pz",  &gen_pz);
  m_MC->SetBranchAddress("gen_x",   &gen_vx);
  m_MC->SetBranchAddress("gen_y",   &gen_vy);
  m_MC->SetBranchAddress("gen_z",   &gen_vz);


  int eventId = -1;
  int n_stubs =  0;
  std::vector<int>   *stub_layer   = NULL;
  //std::vector<int>   *stub_ladder  = NULL;
  //std::vector<int>   *stub_module  = NULL;
  //std::vector<int>   *stub_segment = NULL;
  //std::vector<float> *stub_strip   = NULL;
  std::vector<float> *stub_x = NULL;
  std::vector<float> *stub_y = NULL;
  std::vector<float> *stub_z = NULL;

  m_L1TT->SetBranchAddress("L1Tkevt",            &eventId); 
  m_L1TT->SetBranchAddress("L1TkSTUB_n",         &n_stubs);
  m_L1TT->SetBranchAddress("L1TkSTUB_layer",     &stub_layer);
  //m_L1TT->SetBranchAddress("L1TkSTUB_ladder",    &stub_ladder);
  //m_L1TT->SetBranchAddress("L1TkSTUB_module",    &stub_module);
  //m_L1TT->SetBranchAddress("L1TkSTUB_seg",       &stub_segment);
  //m_L1TT->SetBranchAddress("L1TkSTUB_strip",     &stub_strip);
  m_L1TT->SetBranchAddress("L1TkSTUB_x",         &stub_x);
  m_L1TT->SetBranchAddress("L1TkSTUB_y",         &stub_y);
  m_L1TT->SetBranchAddress("L1TkSTUB_z",         &stub_z);

 
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

  n_events = std::min(m_L1TT->GetEntries(),n_events);
  for (long long int ientry=0; ientry<n_events; ientry++) {
    
    int event_counter = ientry;

    m_L1TT->GetEntry(ientry);
    m_PATT->GetEntry(ientry);
    m_MC->GetEntry(ientry);

   
    // --- Get the stubs per road and trigger tower:
    map <int, set<int> > road_stubs;
    map <int, set<int> > road_stubs_layer;
    map <int, set<int> > trigT_stubs;
    
    for (int iroad=0; iroad<n_patterns; ++iroad) {
  
      const int sector_id = patt_secId->at(iroad);

      for (unsigned int istub=0; istub<patt_links->at(iroad).size(); ++istub ){
	
      	int stubRef =  patt_links->at(iroad).at(istub);
     	road_stubs[iroad].insert(stubRef);
     	road_stubs_layer[iroad].insert(stub_layer->at(stubRef));
      	trigT_stubs[sector_id].insert(stubRef);

      } // istub loop

    } // iroad loop
      

    int n_stubs_PR = 0;
    for (map<int,set<int> >::iterator istub = trigT_stubs.begin(); istub != trigT_stubs.end(); ++istub){
      h_Nstubs_trigT->Fill((istub->second).size());
      n_stubs_PR += (istub->second).size();
    }
    for (map<int,set<int> >::iterator istub = road_stubs.begin(); istub != road_stubs.end(); ++istub)
      h_Nstubs_road->Fill((istub->second).size());
    for (map<int,set<int> >::iterator istub = road_stubs_layer.begin(); istub != road_stubs_layer.end(); ++istub)
      h_Nstubs_layer->Fill((istub->second).size());
    h_Nstubs_tot->Fill(n_stubs);
    h_Nstubs_tot_PR->Fill(n_stubs_PR);


    // --- Do the fit:
    map<int,set<int> >::iterator istub;
    map<int,set<int> >::iterator istub_begin =  trigT_stubs.begin();
    map<int,set<int> >::iterator istub_end   =  trigT_stubs.end();
    if ( fit_roads ) {
      istub_begin =  road_stubs.begin();
      istub_end   =  road_stubs.end();
    }

    vector <Track*> tracks;

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


      // --- Fill the stubs vector:
      vector <Hit_t*> hits;
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
     
     	hits.push_back(hit);
     	
	if (verboseLevel==2 )
	  cout << std::distance((istub->second).begin(),ihit) << "  -  " 
	       << " x = " << hit->x << "   "
	       << " y = " << hit->y << "   "
	       << " z = " << hit->z << " ---> " 
	       << " R = " << hit->rho
	       << "  -  id = " <<  hit->id
	       << endl;
	
      } // ihit loop


      // ===========================================================================
      //  XY fit
      // ===========================================================================

      // --- Phi sector rotation:
      rotateHits(hits, rot_angle[phi_sector]);

      // --- Conformal transformation: 
      confTrans(hits);

      //
      // --- First step ------------------------------------------------------------
      //

      // --- Setup the retina:
      double pbins_step1 = config[trigTow_type]["xy_pbins_step1"];
      double qbins_step1 = config[trigTow_type]["xy_qbins_step1"];
      double pmin_step1  = config[trigTow_type]["xy_pmin_step1"];
      double pmax_step1  = config[trigTow_type]["xy_pmax_step1"];
      double qmin_step1  = config[trigTow_type]["xy_qmin_step1"];
      double qmax_step1  = config[trigTow_type]["xy_qmax_step1"];
  
      double minWeight_step1 = config[trigTow_type]["xy_threshold_step1"];

      double pstep_step1 = (pmax_step1-pmin_step1)/pbins_step1;
      double qstep_step1 = (qmax_step1-qmin_step1)/qbins_step1;

      vector <double> sigma_step1(8,sqrt(pstep_step1*pstep_step1+qstep_step1*qstep_step1));
      if ( config[trigTow_type]["xy_sigma1_step1"] != 0. ) 
	sigma_step1[0] = config[trigTow_type]["xy_sigma1_step1"];
      if ( config[trigTow_type]["xy_sigma2_step1"] != 0. ) 
	sigma_step1[1] = config[trigTow_type]["xy_sigma2_step1"];
      if ( config[trigTow_type]["xy_sigma3_step1"] != 0. ) 
	sigma_step1[2] = config[trigTow_type]["xy_sigma3_step1"];
      if ( config[trigTow_type]["xy_sigma4_step1"] != 0. ) 
	sigma_step1[3] = config[trigTow_type]["xy_sigma4_step1"];
      if ( config[trigTow_type]["xy_sigma5_step1"] != 0. ) 
	sigma_step1[4] = config[trigTow_type]["xy_sigma5_step1"];
      if ( config[trigTow_type]["xy_sigma6_step1"] != 0. ) 
	sigma_step1[5] = config[trigTow_type]["xy_sigma6_step1"];
      if ( config[trigTow_type]["xy_sigma7_step1"] != 0. ) 
	sigma_step1[6] = config[trigTow_type]["xy_sigma7_step1"];
      if ( config[trigTow_type]["xy_sigma8_step1"] != 0. ) 
	sigma_step1[7] = config[trigTow_type]["xy_sigma8_step1"];


      Retina retinaXY_step1(hits, pbins_step1+2, qbins_step1+2, 
			    pmin_step1-pstep_step1, pmax_step1+pstep_step1, 
			    qmin_step1-qstep_step1, qmax_step1+qstep_step1, 
			    sigma_step1, minWeight_step1, 1, XY);


      // --- Fill the retina and find maxima:
      retinaXY_step1.fillGrid();
      retinaXY_step1.findMaxima();
      if ( verboseLevel==1 )
	retinaXY_step1.dumpGrid(event_counter,1,road_id);
      if ( verboseLevel==2 )
	retinaXY_step1.printMaxima();


      // --- Get first step maxima:
      vector <pqPoint> maximaXY_step1 = retinaXY_step1.getMaxima();

      if ( maximaXY_step1.size()==0 && verboseLevel>0 ){
	cout << "*** WARNING in RetinaTrackFitter::fit() at event/road = " << event_counter << "/" 
	     << road_id << ": no maximum found in XY fit-step 1." << endl;
	retinaXY_step1.dumpGrid(event_counter,1,0);
      }

      if ( maximaXY_step1.size() > 10 ){
	cout << "*** ERROR in RetinaTrackFitter::fit() at event/road = " << event_counter 
	     << "/" << road_id << ": " << maximaXY_step1.size() 
	     << " maxima found in XY fit-step 1, fit aborted!" << endl;
	if ( verboseLevel>0 ) 
	  retinaXY_step1.dumpGrid(event_counter,1,0);
	continue;
      }


      //
      // --- Second step -----------------------------------------------------------
      //

      double pbins_step2 = config[trigTow_type]["xy_pbins_step2"];
      double qbins_step2 = config[trigTow_type]["xy_qbins_step2"];

      // --- Zoom around first step maxima:
      for (unsigned int imax=0; imax<maximaXY_step1.size(); ++imax){

	h_max_XY1[trigTow_type]->Fill( maximaXY_step1[imax].p, maximaXY_step1[imax].q);

	// --- Retina setup:
	double pmin_step2 = maximaXY_step1[imax].p - config[trigTow_type]["xy_zoom_step2"]*pstep_step1;
	double pmax_step2 = maximaXY_step1[imax].p + config[trigTow_type]["xy_zoom_step2"]*pstep_step1;
	double qmin_step2 = maximaXY_step1[imax].q - config[trigTow_type]["xy_zoom_step2"]*qstep_step1;
	double qmax_step2 = maximaXY_step1[imax].q + config[trigTow_type]["xy_zoom_step2"]*qstep_step1;
   
	double pstep_step2 = (pmax_step2-pmin_step2)/pbins_step2;
	double qstep_step2 = (qmax_step2-qmin_step2)/qbins_step2;
    
	double minWeight_step2 = config[trigTow_type]["xy_threshold_step2"];

	vector <double> sigma_step2(8,sqrt(pstep_step2*pstep_step2+qstep_step2*qstep_step2));
	if ( config[trigTow_type]["xy_sigma1_step2"] != 0. ) 
	  sigma_step2[0] = config[trigTow_type]["xy_sigma1_step2"];
	if ( config[trigTow_type]["xy_sigma2_step2"] != 0. ) 
	  sigma_step2[1] = config[trigTow_type]["xy_sigma2_step2"];
	if ( config[trigTow_type]["xy_sigma3_step2"] != 0. ) 
	  sigma_step2[2] = config[trigTow_type]["xy_sigma3_step2"];
	if ( config[trigTow_type]["xy_sigma4_step2"] != 0. ) 
	  sigma_step2[3] = config[trigTow_type]["xy_sigma4_step2"];
	if ( config[trigTow_type]["xy_sigma5_step2"] != 0. ) 
	  sigma_step2[4] = config[trigTow_type]["xy_sigma5_step2"];
	if ( config[trigTow_type]["xy_sigma6_step2"] != 0. ) 
	  sigma_step2[5] = config[trigTow_type]["xy_sigma6_step2"];
	if ( config[trigTow_type]["xy_sigma7_step2"] != 0. ) 
	  sigma_step2[6] = config[trigTow_type]["xy_sigma7_step2"];
	if ( config[trigTow_type]["xy_sigma8_step2"] != 0. ) 
	  sigma_step2[7] = config[trigTow_type]["xy_sigma8_step2"];


	Retina retinaXY_step2(hits, pbins_step2+2, qbins_step2+2, 
			      pmin_step2-pstep_step2, pmax_step2+pstep_step2, 
			      qmin_step2-qstep_step2, qmax_step2+qstep_step2, 
			      sigma_step2, minWeight_step2, 1, XY);


	// --- Fill the retina and find maxima:
	retinaXY_step2.fillGrid();
	retinaXY_step2.findMaxima();
	if ( verboseLevel==1 )
	  retinaXY_step2.dumpGrid(event_counter,2,road_id*100+imax);
	if ( verboseLevel==2 )
	  retinaXY_step2.printMaxima();


	// --- Get second step maxima:
	vector <pqPoint> maximaXY_step2 = retinaXY_step2.getMaxima();

	if ( maximaXY_step2.size()==0 && verboseLevel>0 ){
	  cout << "*** WARNING in RetinaTrackFitter::fit() at event/road = " << event_counter << "/" 
	       << road_id << ": no maximum found in XY fit-step 2 (step-1 XY maximum #" << imax << ")." 
	       << endl;
	  retinaXY_step2.dumpGrid(event_counter,2,imax);
	}

	if ( maximaXY_step2.size() > 10 ){
	  cout << "*** ERROR in RetinaTrackFitter::fit() at event/road = " << event_counter 
	       << "/" << road_id << ": " << maximaXY_step2.size() 
	       << " maxima found in XY fit-step 2 (step-1 XY maximum #" << imax << "), fit aborted!" 
	       << endl;
	  if ( verboseLevel>0 ) 
	    retinaXY_step2.dumpGrid(event_counter,2,imax);
	  continue;
	}


	// --- Loop over XY second step maxima:
	for (unsigned int itrk=0; itrk<maximaXY_step2.size(); ++itrk){

	  h_max_XY2[trigTow_type]->Fill( maximaXY_step2[itrk].p, maximaXY_step2[itrk].q);

	  // --- Invert the X+-X- transformation:
	  double p = 0.5*(y1 - y0)/maximaXY_step2[itrk].q;
	  double q = y0 - p*(maximaXY_step2[itrk].p-maximaXY_step2[itrk].q);


	  // --- Associate stubs to this maximum:
	  vector <Hit_t*> hits_RZ;
	  unsigned int n_stubsPS = 0;
	  for (unsigned int ihit=0; ihit<hits.size(); ++ihit){
      
	    double dist   = (p*hits[ihit]->x-hits[ihit]->y+q)/p;
	    //double dist   = fabs(hits[ihit]->y-p*hits[ihit]->x-q)/sqrt(1.+p*p);
	    double weight = exp(-0.5*dist*dist/(sigma_step2[0]*sigma_step2[0]));

	    if ( weight > 0.5 ){
	      hits_RZ.push_back(hits[ihit]);
	      if ( hits[ihit]->rho > 60. )
		n_stubsPS++;
	    }
	    else
	      if ( verboseLevel>0 )
		cout << "*** WARNING in RetinaTrackFitter::fit() at event/road = " << event_counter 
		     << "/" << road_id << ": stub " << hits[ihit]->id << " with weight = " << weight 
		     << " has not been associated to the XY step-2 maximum #" << imax << "." << endl;
	  } // ihit loop
	  if ( hits_RZ.size() < 3 ){
	    if ( verboseLevel>0 )
	      cout << "*** ERROR in RetinaTrackFitter::fit() at event/road = " << event_counter 
		   << "/" << road_id << ": only " <<  hits_RZ.size() 
		   << " stubs associated to the XY step-2 maximum #" << imax << ", fit aborted!" << endl;
	    continue;
	  }

	  // --- Rotate back the original phi sector:
	  //q = q/(cos(rot_angle[phi_sector])+p*sin(rot_angle[phi_sector]));
	  //p = (p*cos(rot_angle[phi_sector])-sin(rot_angle[phi_sector]))/
	  //  (cos(rot_angle[phi_sector])+p*sin(rot_angle[phi_sector]));


	  // --- Invert the conformal transformation and get the track parameters:
	  double a = -0.5*p/q;
	  double b =  0.5/q;
      
	  //  TTTrack might not be ready for a signed curvature!!!!!!!!!!!!!!!!!
	  double charge = TMath::Sign(1.,b);
	  double c   = charge/sqrt(a*a+b*b);
	  //double c   = 1./sqrt(a*a+b*b);

	  double phi = atan(p)-rot_angle[phi_sector];
	  if ( phi>TMath::Pi() )
	    phi -= TMath::TwoPi();
      
	  // =========================================================================
	  //  RZ fit
	  // =========================================================================

	  y0 = 0.5/y0;
	  y1 = 0.5/y1;

	  double eta = -9999.;
	  double z0  = -9999.;
    

	  //
	  // --- First step ----------------------------------------------------------
	  //

	  pbins_step1 = config[trigTow_type]["rz_pbins_step1"];
	  qbins_step1 = config[trigTow_type]["rz_qbins_step1"];
	  pmin_step1  = config[trigTow_type]["rz_pmin_step1"];
	  pmax_step1  = config[trigTow_type]["rz_pmax_step1"];
	  qmin_step1  = config[trigTow_type]["rz_qmin_step1"];
	  qmax_step1  = config[trigTow_type]["rz_qmax_step1"];

	  minWeight_step1 = config[trigTow_type]["rz_threshold_step1"];

	  pstep_step1 = (pmax_step1-pmin_step1)/pbins_step1;
	  qstep_step1 = (qmax_step1-qmin_step1)/qbins_step1;

	  for (unsigned int ilayer=0; ilayer<8; ++ilayer)
	    sigma_step1[ilayer] = sqrt(pstep_step1*pstep_step1+qstep_step1*qstep_step1);

	  if ( config[trigTow_type]["rz_sigma1_step1"] != 0. ) 
	    sigma_step1[0] = config[trigTow_type]["rz_sigma1_step1"];
	  if ( config[trigTow_type]["rz_sigma2_step1"] != 0. ) 
	    sigma_step1[1] = config[trigTow_type]["rz_sigma2_step1"];
	  if ( config[trigTow_type]["rz_sigma3_step1"] != 0. ) 
	    sigma_step1[2] = config[trigTow_type]["rz_sigma3_step1"];
	  if ( config[trigTow_type]["rz_sigma4_step1"] != 0. ) 
	    sigma_step1[3] = config[trigTow_type]["rz_sigma4_step1"];
	  if ( config[trigTow_type]["rz_sigma5_step1"] != 0. ) 
	    sigma_step1[4] = config[trigTow_type]["rz_sigma5_step1"];
	  if ( config[trigTow_type]["rz_sigma6_step1"] != 0. ) 
	    sigma_step1[5] = config[trigTow_type]["rz_sigma6_step1"];
	  if ( config[trigTow_type]["rz_sigma7_step1"] != 0. ) 
	    sigma_step1[6] = config[trigTow_type]["rz_sigma7_step1"];
	  if ( config[trigTow_type]["rz_sigma8_step1"] != 0. ) 
	    sigma_step1[7] = config[trigTow_type]["rz_sigma8_step1"];

      
	  Retina retinaRZ_step1(hits_RZ, pbins_step1+2, qbins_step1+2, 
				pmin_step1-pstep_step1, pmax_step1+pstep_step1, 
				qmin_step1-qstep_step1, qmax_step1+qstep_step1, 
				sigma_step1, minWeight_step1, 1, RZ);


	  retinaRZ_step1.fillGrid();
	  retinaRZ_step1.findMaxima();
	  if ( verboseLevel==1 )
	    retinaRZ_step1.dumpGrid(event_counter,1,imax);
	  if ( verboseLevel==2 )
	    retinaRZ_step1.printMaxima();


	  // --- Get first step maximum:
	  vector <pqPoint> maximaRZ_step1 = retinaRZ_step1.getMaxima();

	  if ( maximaRZ_step1.size()==0 && verboseLevel>0 ){
	    cout << "*** WARNING in RetinaTrackFitter::fit() at event/road = " << event_counter 
		 << "/" << road_id << ": no maximum found in RZ fit-step 1 (step-1 maximum #" 
		 << imax << ")." << endl;
	    retinaXY_step2.dumpGrid(event_counter,2,imax);
	  }

	  if ( maximaRZ_step1.size() > 10 ){
	    cout << "*** ERROR in RetinaTrackFitter::fit() at event/road = " << event_counter 
		 << "/" << road_id << ": " << maximaRZ_step1.size() 
		 << " maxima found in RZ fit-step 1 (step-1 maximum #" << imax << "), fit aborted!" 
		 << endl;
	    if ( verboseLevel>0 ) 
	      retinaRZ_step1.dumpGrid(event_counter,1,imax);
	    continue;
	  }

      
	  //
	  // --- Second step ---------------------------------------------------------
	  //

	  pbins_step2 = config[trigTow_type]["rz_pbins_step2"];
	  qbins_step2 = config[trigTow_type]["rz_qbins_step2"];

	  // Zoom around first step maxima
	  for (unsigned int imax_RZ=0; imax_RZ<maximaRZ_step1.size(); ++imax_RZ){
	    
	    h_max_RZ1[trigTow_type]->Fill( maximaRZ_step1[imax_RZ].p, maximaRZ_step1[imax_RZ].q);

	    double pmin_step2 = maximaRZ_step1[imax_RZ].p - config[trigTow_type]["rz_zoom_step2"]*pstep_step1;
	    double pmax_step2 = maximaRZ_step1[imax_RZ].p + config[trigTow_type]["rz_zoom_step2"]*pstep_step1;
	    double qmin_step2 = maximaRZ_step1[imax_RZ].q - config[trigTow_type]["rz_zoom_step2"]*qstep_step1;
	    double qmax_step2 = maximaRZ_step1[imax_RZ].q + config[trigTow_type]["rz_zoom_step2"]*qstep_step1;
   
	    double pstep_step2 = (pmax_step2-pmin_step2)/pbins_step2;
	    double qstep_step2 = (qmax_step2-qmin_step2)/qbins_step2;
    
	    double minWeight_step2 = config[trigTow_type]["rz_threshold_step2"];
	    // If less than 3 PS stubs are available, adjust the maximum-finder threshold:
	    if (  n_stubsPS < 3 )
	      minWeight_step2 *= 0.66;

	    vector <double> sigma_step2(8,sqrt(pstep_step2*pstep_step2+qstep_step2*qstep_step2));
	    for (unsigned int ilayer=3; ilayer<6; ++ilayer)
	      sigma_step2[ilayer] = 8.*sqrt(pstep_step2*pstep_step2+qstep_step2*qstep_step2);

	    if ( config[trigTow_type]["rz_sigma1_step2"] != 0. ) 
	      sigma_step2[0] = config[trigTow_type]["rz_sigma1_step2"];
	    if ( config[trigTow_type]["rz_sigma2_step2"] != 0. ) 
	      sigma_step2[1] = config[trigTow_type]["rz_sigma2_step2"];
	    if ( config[trigTow_type]["rz_sigma3_step2"] != 0. ) 
	      sigma_step2[2] = config[trigTow_type]["rz_sigma3_step2"];
	    if ( config[trigTow_type]["rz_sigma4_step2"] != 0. ) 
	      sigma_step2[3] = config[trigTow_type]["rz_sigma4_step2"];
	    if ( config[trigTow_type]["rz_sigma5_step2"] != 0. ) 
	      sigma_step2[4] = config[trigTow_type]["rz_sigma5_step2"];
	    if ( config[trigTow_type]["rz_sigma6_step2"] != 0. ) 
	      sigma_step2[5] = config[trigTow_type]["rz_sigma6_step2"];
	    if ( config[trigTow_type]["rz_sigma7_step2"] != 0. ) 
	      sigma_step2[6] = config[trigTow_type]["rz_sigma7_step2"];
	    if ( config[trigTow_type]["rz_sigma8_step2"] != 0. ) 
	      sigma_step2[7] = config[trigTow_type]["rz_sigma8_step2"];


	    Retina retinaRZ_step2(hits_RZ, pbins_step2+2, qbins_step2+2, 
				  pmin_step2-pstep_step2, pmax_step2+pstep_step2, 
				  qmin_step2-qstep_step2, qmax_step2+qstep_step2, 
				  sigma_step2, minWeight_step2, 1, RZ);


	    retinaRZ_step2.fillGrid();
	    retinaRZ_step2.findMaxima();
	    if ( verboseLevel==1 )
	      retinaRZ_step2.dumpGrid(event_counter,2,imax*100+imax_RZ);
	    if ( verboseLevel==2 )
	      retinaRZ_step2.printMaxima();

	    pqPoint bestpqRZ_step2 = retinaRZ_step2.getBestPQ();

	    h_max_RZ2[trigTow_type]->Fill(bestpqRZ_step2.p,bestpqRZ_step2.q);


	    // --- If no RZ maximum is found, skip the road:
	    if ( bestpqRZ_step2.w == -1. ) {
	      if ( verboseLevel>0 )
		cout << "*** ERROR in RetinaTrackFitter::fit() at event/road = " << event_counter 
		     << "/" << road_id << ": no maximum found in RZ fit-step 2 (step-1 XY maximum #" 
		     << imax << ", step-1 RZ maximum #" << imax_RZ << "), fit aborted!" << endl;
	      continue;
	    }

	    if ( retinaRZ_step2.getMaxima().size() > 10 ){
	      cout << "*** ERROR in RetinaTrackFitter::fit() at event/road = " << event_counter 
		   << "/" << road_id << ": " << retinaRZ_step2.getMaxima().size()
		   << " maxima found in RZ fit-step 2 (step-1 maximum #" << imax 
		   << ", step-1 RZ maximum #" << imax_RZ << "), fit aborted!" << endl;
	      if ( verboseLevel>0 ) 
		retinaRZ_step2.dumpGrid(event_counter,2,imax*100+imax_RZ);
	      continue;
	    }


	    // --- Invert the X+-X- transformation:
	    double p = 0.5*(y1 - y0)/bestpqRZ_step2.q;
	    double q = y0 - p*(bestpqRZ_step2.p-bestpqRZ_step2.q);


	    // --- Get the track parameters:
	    double theta = atan(p);
	    if ( theta < 0. ){
	      if ( verboseLevel>0 ) 
		cout << "*** WARNING in RetinaTrackFitter::fit() at event/road = " << event_counter 
		     << "/" << road_id << ": theta corrected from " << theta << " to " << theta+TMath::Pi() 
		     << endl;
	      theta += TMath::Pi();
	    }
	    eta = -log(tan(0.5*theta));
	    z0  = -q/p;


	    // --- Invert eta and z0 signs if we are fitting a negative-eta tower:
	    if ( eta_range < 3 ){
	      eta = -eta;
	      z0  = -z0;
	    }


	    // --- Save the track:
	    Track* trk = new Track(c, 0., phi, eta, z0, maximaXY_step2[itrk].w, bestpqRZ_step2.w);
	    for(unsigned int ihit=0; ihit<hits_RZ.size(); ihit++)
	      trk->addStubIndex(hits_RZ[ihit]->id);

	    tracks.push_back(trk);

	  } // imax_RZ loop


	  hits_RZ.clear();

	} // itrk loop
 

      } // imax loop


      // --- Clean-up pointers:
      for(vector<Hit_t*>::iterator it=hits.begin(); it!=hits.end(); ++it)
	delete *it;
      hits.clear();

    } // istub loop



    if ( verboseLevel == 3 ){

      std::cout << "Event = " <<  ientry 
		<< " ----------------------------------------------------------------------" << std::endl;

      //  --- Printout the generated particles:

      std::cout << " Generated particles:" << std::endl;
      for (int ipart=0; ipart<n_gen; ++ipart) {

	// --- Keep only charged stable particles:
	if ( fabs(gen_pdg->at(ipart)) == 15   ||
	     fabs(gen_pdg->at(ipart)) == 12   ||
	     fabs(gen_pdg->at(ipart)) == 14   ||
	     fabs(gen_pdg->at(ipart)) == 16   ||
	     fabs(gen_pdg->at(ipart)) == 22   ||
	     fabs(gen_pdg->at(ipart)) == 111  ||
	     fabs(gen_pdg->at(ipart)) == 113  ||
	     fabs(gen_pdg->at(ipart)) == 130  ||
	     fabs(gen_pdg->at(ipart)) == 221  ||
	     fabs(gen_pdg->at(ipart)) == 310  ||
	     fabs(gen_pdg->at(ipart)) == 331  ||
	     fabs(gen_pdg->at(ipart)) == 2112 ||
	     fabs(gen_pdg->at(ipart)) > 3000  ) continue;

	double gen_charge = TMath::Sign(1., (Double_t) gen_pdg->at(ipart));

	double gen_pt  = sqrt(gen_px->at(ipart)*gen_px->at(ipart)+gen_py->at(ipart)*gen_py->at(ipart));
	double gen_phi = atan2(gen_py->at(ipart),gen_px->at(ipart));
	//if (gen_phi<0.)
	//gen_phi += TMath::TwoPi();
	double gen_theta = atan2(gen_pt,gen_pz->at(ipart));
	double gen_eta = -log(tan(0.5*gen_theta));

	// curvature and helix radius:
	double c = gen_charge*0.003*mMagneticField/gen_pt;
	double R = gen_pt/(0.003*mMagneticField);
	  
	// helix center:
	double x0 = gen_vx->at(ipart) - gen_charge*R*gen_py->at(ipart)/gen_pt;
	double y0 = gen_vy->at(ipart) + gen_charge*R*gen_px->at(ipart)/gen_pt;

	// transverse and longitudinal impact parameters:
	double gen_d0 = gen_charge*(sqrt(x0*x0+y0*y0)-R);
	double gen_z0 = gen_vz->at(ipart) - 2./c*gen_pz->at(ipart)/gen_pt*
	  asin(0.5*c*sqrt((gen_vx->at(ipart)*gen_vx->at(ipart)+
			   gen_vy->at(ipart)*gen_vy->at(ipart)-gen_d0*gen_d0)/(1+c*gen_d0)));


	// --- Print out the generated particles: 
	std::cout << "  " << ipart << "  -  ID = " << gen_pdg->at(ipart)
		  << "  pt = "  << gen_pt
		  << "  d0 = "  << gen_d0
		  << "  phi = " << gen_phi
		  << "  eta = " << gen_eta
		  << "  z0 = "  << gen_z0 
		  << std::endl; 


      } // ipart loop



      //  --- Printout the fitted tracks:
      cout << " Fitted tracks (no duplicate removal):" << endl;
      for ( std::vector<Track*>::iterator itrk=tracks.begin(); itrk!=tracks.end(); ++itrk ){

	double pt = 0.003*mMagneticField/(*itrk)->getCurve();

	cout << "  "
	     << std::distance(tracks.begin(),itrk)
	     << "  -  c = " << (*itrk)->getCurve()
	     << "  pt = "  << pt
	     << "  phi = " << (*itrk)->getPhi0()
	     << "  eta = " << (*itrk)->getEta0()
	     << "  z0 = "  << (*itrk)->getZ0()
	     << endl;

      } // itrk loop

    } // if ( verboseLevel == 3 )
    
    // --- Clean-up the heap:
    for ( std::vector<Track*>::iterator itrk=tracks.begin(); itrk!=tracks.end(); ++itrk )
      delete *itrk;
    tracks.clear();


  } // ientry loop



  // --- Save histograms:
  TFile f("histos.root","recreate");

  h_Nstubs_tot->Write();  
  h_Nstubs_tot_PR->Write();  
  h_Nstubs_road->Write();  
  h_Nstubs_layer->Write();  
  h_Nstubs_trigT->Write();  
  for (int ihist=0; ihist<3; ++ihist){
    h_max_XY1[ihist]->Write();
    h_max_XY2[ihist]->Write();
    h_max_RZ1[ihist]->Write();
    h_max_RZ2[ihist]->Write();
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
  config[0]["xy_pbins_step1"]     = 40.;
  config[0]["xy_qbins_step1"]     = 40.;
  config[0]["xy_pmin_step1"]      = -0.05;
  config[0]["xy_pmax_step1"]      =  0.05;
  config[0]["xy_qmin_step1"]      = -0.05;
  config[0]["xy_qmax_step1"]      =  0.05;
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

  config[0]["rz_pbins_step1"]     = 20.;
  config[0]["rz_qbins_step1"]     = 20.;
  config[0]["rz_pmin_step1"]      = -20.;
  config[0]["rz_pmax_step1"]      =  60.;
  config[0]["rz_qmin_step1"]      = -60.;
  config[0]["rz_qmax_step1"]      =  60.;
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
  config[1]["xy_pbins_step1"]     = 40.;
  config[1]["xy_qbins_step1"]     = 40.;
  config[1]["xy_pmin_step1"]      = -0.05;
  config[1]["xy_pmax_step1"]      =  0.05;
  config[1]["xy_qmin_step1"]      = -0.05;
  config[1]["xy_qmax_step1"]      =  0.05;
  config[1]["xy_threshold_step1"] =  4.;
  config[1]["xy_sigma1_step1"]    =  0.;
  config[1]["xy_sigma2_step1"]    =  0.;
  config[1]["xy_sigma3_step1"]    =  0.;
  config[1]["xy_sigma4_step1"]    =  0.;
  config[1]["xy_sigma5_step1"]    =  0.;
  config[1]["xy_sigma6_step1"]    =  0.;
  config[1]["xy_sigma7_step1"]    =  0.;
  config[1]["xy_sigma8_step1"]    =  0.;
  config[1]["xy_pbins_step2"]     = 100.;
  config[1]["xy_qbins_step2"]     = 100.;
  config[1]["xy_zoom_step2"]      = 1.;
  config[1]["xy_threshold_step2"] =  4.;
  config[1]["xy_sigma1_step2"]    =  0.;
  config[1]["xy_sigma2_step2"]    =  0.;
  config[1]["xy_sigma3_step2"]    =  0.;
  config[1]["xy_sigma4_step2"]    =  0.;
  config[1]["xy_sigma5_step2"]    =  0.;
  config[1]["xy_sigma6_step2"]    =  0.;
  config[1]["xy_sigma7_step2"]    =  0.;
  config[1]["xy_sigma8_step2"]    =  0.;

  config[1]["rz_pbins_step1"]     = 20.;
  config[1]["rz_qbins_step1"]     = 20.;
  config[1]["rz_pmin_step1"]      = 40.;
  config[1]["rz_pmax_step1"]      = 140.;
  config[1]["rz_qmin_step1"]      = 0.;
  config[1]["rz_qmax_step1"]      = 120.;
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
  config[2]["xy_pbins_step1"]     = 40.;
  config[2]["xy_qbins_step1"]     = 40.;
  config[2]["xy_pmin_step1"]      = -0.05;
  config[2]["xy_pmax_step1"]      =  0.05;
  config[2]["xy_qmin_step1"]      = -0.05;
  config[2]["xy_qmax_step1"]      =  0.05;
  config[2]["xy_threshold_step1"] =  4.;
  config[2]["xy_sigma1_step1"]    =  0.;
  config[2]["xy_sigma2_step1"]    =  0.;
  config[2]["xy_sigma3_step1"]    =  0.;
  config[2]["xy_sigma4_step1"]    =  0.;
  config[2]["xy_sigma5_step1"]    =  0.;
  config[2]["xy_sigma6_step1"]    =  0.;
  config[2]["xy_sigma7_step1"]    =  0.;
  config[2]["xy_sigma8_step1"]    =  0.;
  config[2]["xy_pbins_step2"]     = 100.;
  config[2]["xy_qbins_step2"]     = 100.;
  config[2]["xy_zoom_step2"]      = 1.;
  config[2]["xy_threshold_step2"] =  4.;
  config[2]["xy_sigma1_step2"]    =  0.;
  config[2]["xy_sigma2_step2"]    =  0.;
  config[2]["xy_sigma3_step2"]    =  0.;
  config[2]["xy_sigma4_step2"]    =  0.;
  config[2]["xy_sigma5_step2"]    =  0.;
  config[2]["xy_sigma6_step2"]    =  0.;
  config[2]["xy_sigma7_step2"]    =  0.;
  config[2]["xy_sigma8_step2"]    =  0.;

  config[2]["rz_pbins_step1"]     = 20.;
  config[2]["rz_qbins_step1"]     = 20.;
  config[2]["rz_pmin_step1"]      = 140.;
  config[2]["rz_pmax_step1"]      = 240.;
  config[2]["rz_qmin_step1"]      = 80.;
  config[2]["rz_qmax_step1"]      = 180.;
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
