
#include "FriendTreeProcessor.hpp"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

using namespace std;
using namespace edm;

FriendTreeProcessor::FriendTreeProcessor()
{
  Service<TFileService> fs;
  
  outputTree = fs->make<TTree>("tree","tree");
  outputTree->Branch("pion_vx", &pion_vx);
  outputTree->Branch("pion_vy", &pion_vy);
  outputTree->Branch("pion_vz", &pion_vz);
  outputTree->Branch("pion_px", &pion_px);
  outputTree->Branch("pion_py", &pion_py);
  outputTree->Branch("pion_pz", &pion_pz);
  outputTree->Branch("pion_charge", &pion_charge);
  outputTree->Branch("pion_simHits_x", &pion_simHits_x);
  outputTree->Branch("pion_simHits_y", &pion_simHits_y);
  outputTree->Branch("pion_simHits_z", &pion_simHits_z);
  outputTree->Branch("pion_simHits_t", &pion_simHits_t);
  outputTree->Branch("pion_simHits_subDet", &pion_simHits_subDet);
  
  outputTree->Branch("chargino_eta", &chargino_eta);
  outputTree->Branch("chargino_phi", &chargino_phi);
  outputTree->Branch("chargino_pt", &chargino_pt);
  outputTree->Branch("chargino_px", &chargino_px);
  outputTree->Branch("chargino_py", &chargino_py);
  outputTree->Branch("chargino_pz", &chargino_pz);
  outputTree->Branch("chargino_charge", &chargino_charge);
  outputTree->Branch("chargino_nTrackerLayers", &chargino_nTrackerLayers);
  outputTree->Branch("chargino_simHits_x", &chargino_simHits_x);
  outputTree->Branch("chargino_simHits_y", &chargino_simHits_y);
  outputTree->Branch("chargino_simHits_z", &chargino_simHits_z);
  outputTree->Branch("chargino_simHits_subDet", &chargino_simHits_subDet);
  
  outputTree->Branch("pixelCluster_x", &pixelCluster_x);
  outputTree->Branch("pixelCluster_y", &pixelCluster_y);
  outputTree->Branch("pixelCluster_z", &pixelCluster_z);
  outputTree->Branch("pixelCluster_charge", &pixelCluster_charge);
  outputTree->Branch("pixelCluster_subDet", &pixelCluster_subDet);
  
  outputTree->Branch("stripCluster_x", &stripCluster_x);
  outputTree->Branch("stripCluster_y", &stripCluster_y);
  outputTree->Branch("stripCluster_z", &stripCluster_z);
  outputTree->Branch("stripCluster_ex", &stripCluster_ex);
  outputTree->Branch("stripCluster_ey", &stripCluster_ey);
  outputTree->Branch("stripCluster_ez", &stripCluster_ez);
  outputTree->Branch("stripCluster_charge", &stripCluster_charge);
  outputTree->Branch("stripCluster_subDet", &stripCluster_subDet);
  
  outputTree->Branch("pionCluster_x", &pionCluster_x);
  outputTree->Branch("pionCluster_y", &pionCluster_y);
  outputTree->Branch("pionCluster_z", &pionCluster_z);
  outputTree->Branch("pionCluster_ex", &pionCluster_ex);
  outputTree->Branch("pionCluster_ey", &pionCluster_ey);
  outputTree->Branch("pionCluster_ez", &pionCluster_ez);
  outputTree->Branch("pionCluster_charge", &pionCluster_charge);
  outputTree->Branch("pionCluster_subDet", &pionCluster_subDet);
  
  outputTree->Branch("generalTrack_px" , &generalTrack_px);
  outputTree->Branch("generalTrack_py" , &generalTrack_py);
  outputTree->Branch("generalTrack_pz" , &generalTrack_pz);
  outputTree->Branch("generalTrack_nLoops" , &generalTrack_nLoops);
  outputTree->Branch("generalTrack_isLooper" , &generalTrack_isLooper);
  outputTree->Branch("generalTrack_nPionHits" , &generalTrack_nPionHits);
  outputTree->Branch("generalTrack_d0" , &generalTrack_d0);
  outputTree->Branch("generalTrack_charge" , &generalTrack_charge);
  outputTree->Branch("generalTrack_chi2" , &generalTrack_chi2);
  outputTree->Branch("generalTrack_eta" , &generalTrack_eta);
  outputTree->Branch("generalTrack_phi" , &generalTrack_phi);
  outputTree->Branch("generalTrack_nHits" , &generalTrack_nHits);
  outputTree->Branch("generalTrack_nMissingHits" , &generalTrack_nMissingHits);
  
  
  outputTree->Branch("runNumber", &run);
  outputTree->Branch("lumiBlock", &lumi);
  outputTree->Branch("eventNumber", &event);
  
}

FriendTreeProcessor::~FriendTreeProcessor()
{
  
}

void FriendTreeProcessor::clearVectors()
{
  pion_vx.clear();
  pion_vy.clear();
  pion_vz.clear();
  pion_px.clear();
  pion_py.clear();
  pion_pz.clear();
  pion_charge.clear();
  
  pion_simHits_x.clear();
  pion_simHits_y.clear();
  pion_simHits_z.clear();
  pion_simHits_t.clear();
  pion_simHits_subDet.clear();
  
  chargino_eta.clear();
  chargino_phi.clear();
  chargino_pt.clear();
  chargino_px.clear();
  chargino_py.clear();
  chargino_pz.clear();
  chargino_charge.clear();
  chargino_nTrackerLayers.clear();
  chargino_simHits_x.clear();
  chargino_simHits_y.clear();
  chargino_simHits_z.clear();
  chargino_simHits_subDet.clear();
  
  pixelCluster_x.clear();
  pixelCluster_y.clear();
  pixelCluster_z.clear();
  pixelCluster_charge.clear();
  pixelCluster_subDet.clear();
  
  stripCluster_x.clear();
  stripCluster_y.clear();
  stripCluster_z.clear();
  stripCluster_ex.clear();
  stripCluster_ey.clear();
  stripCluster_ez.clear();
  stripCluster_charge.clear();
  stripCluster_subDet.clear();
  
  pionCluster_x.clear();
  pionCluster_y.clear();
  pionCluster_z.clear();
  pionCluster_ex.clear();
  pionCluster_ey.clear();
  pionCluster_ez.clear();
  pionCluster_charge.clear();
  pionCluster_subDet.clear();
  
  generalTrack_px.clear();
  generalTrack_py.clear();
  generalTrack_pz.clear();
  generalTrack_nLoops.clear();
  generalTrack_isLooper.clear();
  generalTrack_nPionHits.clear();
  generalTrack_d0.clear();
  generalTrack_charge.clear();
  generalTrack_chi2.clear();
  generalTrack_eta.clear();
  generalTrack_phi.clear();
  generalTrack_nHits.clear();
  generalTrack_nMissingHits.clear();
}

void FriendTreeProcessor::fill()
{
  outputTree->Fill();
}
