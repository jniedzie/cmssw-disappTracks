#include <vector>

#include <TTree.h>

using namespace std;

class FriendTreeProcessor
{
public:
  FriendTreeProcessor();
  ~FriendTreeProcessor();
  
  void clearVectors();
  
  void fill();
  
  TTree *outputTree;
  
  uint lumi;
  uint run;
  unsigned long long event;
  
  // GEN-SIM level pion information
  vector<double> pion_vx;
  vector<double> pion_vy;
  vector<double> pion_vz;
  vector<double> pion_px;
  vector<double> pion_py;
  vector<double> pion_pz;
  vector<double> pion_simHits_x;
  vector<double> pion_simHits_y;
  vector<double> pion_simHits_z;
  vector<double> pion_simHits_t;
  vector<int> pion_simHits_subDet;
  vector<int> pion_charge;
  
  // GEN-SIM level chargino information
  vector<double>  chargino_eta;
  vector<double>  chargino_phi;
  vector<double>  chargino_pt;
  vector<double>  chargino_px;
  vector<double>  chargino_py;
  vector<double>  chargino_pz;
  vector<int>     chargino_charge;
  vector<int>     chargino_nTrackerLayers;
  vector<double>  chargino_simHits_x;
  vector<double>  chargino_simHits_y;
  vector<double>  chargino_simHits_z;
  vector<int>     chargino_simHits_subDet;
  
  // Tracker clusters not assigned to any track
  vector<double> pixelCluster_x;
  vector<double> pixelCluster_y;
  vector<double> pixelCluster_z;
  vector<double> pixelCluster_charge;
  vector<double> pixelCluster_subDet;
  
  vector<double> stripCluster_x;
  vector<double> stripCluster_y;
  vector<double> stripCluster_z;
  vector<double> stripCluster_ex;
  vector<double> stripCluster_ey;
  vector<double> stripCluster_ez;
  vector<double> stripCluster_charge;
  vector<double> stripCluster_subDet;
  
  vector<double> pionCluster_x;
  vector<double> pionCluster_y;
  vector<double> pionCluster_z;
  vector<double> pionCluster_ex;
  vector<double> pionCluster_ey;
  vector<double> pionCluster_ez;
  vector<double> pionCluster_charge;
  vector<double> pionCluster_subDet;
  
  vector<double> generalTrack_px;
  vector<double> generalTrack_py;
  vector<double> generalTrack_pz;
  vector<int>    generalTrack_nLoops;
  vector<bool>   generalTrack_isLooper;
  vector<int>    generalTrack_nPionHits;
  vector<double> generalTrack_d0;
  vector<double> generalTrack_charge;
  vector<double> generalTrack_chi2;
  vector<double> generalTrack_eta;
  vector<double> generalTrack_phi;
  vector<double> generalTrack_nHits;
  vector<double> generalTrack_nMissingHits;
};
