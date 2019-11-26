
// -*- C++ -*-
//
// Package:    CharginoAnalysis/CharginoAnalyzer
// Class:      CharginoAnalyzer
//
/**\class CharginoAnalyzer CharginoAnalyzer.cc CharginoAnalysis/CharginoAnalyzer/plugins/CharginoAnalyzer.cc
 
 Description: [one line class summary]
 
 Implementation:
 [Notes on implementation]
 */
//
// Original Author:  Jeremi Niedziela
//         Created:  Tue, 19 Feb 2019 12:48:21 GMT
//
//

#include <memory>
#include <string>
#include <iomanip>
#include <fstream>

#include <TLorentzVector.h>
#include <TH1D.h>
#include <TTree.h>

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/CoreSimTrack.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/GeomDetEnumerators.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"

#include "RecoLocalTracker/SiPixelRecHits/interface/PixelCPEBase.h"
#include "RecoLocalTracker/Records/interface/TkPixelCPERecord.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"

using namespace reco;
using namespace std;
using namespace edm;

class CharginoAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit CharginoAnalyzer(const edm::ParameterSet&);
  ~CharginoAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  void print(const TrackingParticle *particle);
  void FillHitsForTrackID(const PSimHit &hit, uint trackID, bool isPixel, bool isPion, bool isChargino);
  void FillSimHitsForTrack(uint trackID, bool isPion, bool isChargino);
  void FillTrackerClusters(bool filterTrackClusters);
  void PrintTrackingParticles();
  
  // Tokens to get necessary objects
  EDGetTokenT<vector<GenParticle>> genParticlesToken;
  EDGetTokenT<vector<SimTrack>> simTracksToken;
  EDGetTokenT<vector<PSimHit>> simHitsPixelHighToken;
  EDGetTokenT<vector<PSimHit>> simHitsPixelLowToken;
  EDGetTokenT<vector<PSimHit>> simHitsPixelEndcapLowToken;
  EDGetTokenT<vector<PSimHit>> simHitsPixelEndcapHighToken;
  EDGetTokenT<vector<PSimHit>> simHitsTEClowToken;
  EDGetTokenT<vector<PSimHit>> simHitsTEChighToken;
  EDGetTokenT<vector<PSimHit>> simHitsTIBlowToken;
  EDGetTokenT<vector<PSimHit>> simHitsTIBhighToken;
  EDGetTokenT<vector<PSimHit>> simHitsTIDlowToken;
  EDGetTokenT<vector<PSimHit>> simHitsTIDhighToken;
  EDGetTokenT<vector<PSimHit>> simHitsTOBhighToken;
  EDGetTokenT<vector<PSimHit>> simHitsTOBlowToken;
  EDGetTokenT<OwnVector<TrackingRecHit,ClonePolicy<TrackingRecHit>>> trackingRecHitsToken;
  EDGetTokenT<vector<Track>> recTracksToken;
  EDGetTokenT<vector<TrackingParticle>> trackingParticlesToken;
  EDGetTokenT<edmNew::DetSetVector<SiPixelCluster>> pixelClusterToken;
  EDGetTokenT<edmNew::DetSetVector<SiStripCluster>> stripClusterToken;
  EDGetTokenT<std::vector<reco::Track>> generalTracksToken;
  
  // Handles to get objects
  Handle<vector<GenParticle>> genParticles;
  Handle<vector<SimTrack>> simTracks;
  Handle<vector<PSimHit>> simHitsPixelHigh;
  Handle<vector<PSimHit>> simHitsPixelLow;
  Handle<vector<PSimHit>> simHitsPixelEndcapLow;
  Handle<vector<PSimHit>> simHitsPixelEndcapHigh;
  Handle<vector<PSimHit>> simHitsTEClow;
  Handle<vector<PSimHit>> simHitsTEChigh;
  Handle<vector<PSimHit>> simHitsTIBlow;
  Handle<vector<PSimHit>> simHitsTIBhigh;
  Handle<vector<PSimHit>> simHitsTIDlow;
  Handle<vector<PSimHit>> simHitsTIDhigh;
  Handle<vector<PSimHit>> simHitsTOBhigh;
  Handle<vector<PSimHit>> simHitsTOBlow;
  Handle<OwnVector<TrackingRecHit,ClonePolicy<TrackingRecHit>>> trackingRecHitsRef;
  Handle<vector<Track>> recTracks;
  Handle<vector<TrackingParticle>> trackingParticles;
  Handle<edmNew::DetSetVector<SiPixelCluster>> pixelClustersHandle;
  Handle<edmNew::DetSetVector<SiStripCluster>> stripClustersHandle;
  Handle<vector<Track>> generalTracksHandle;
  
  ESHandle<TrackerGeometry> trackerGeometry;
  
  TTree *outputTree;
  
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
  
  uint lumi;
  uint run;
  unsigned long long event;
  
  bool verbose;
  int nCharginoHits;
  
  TrackerHitAssociator::Config trackerHitAssociatorConfig;
  TrackerHitAssociator *trackerHitAssociator;
  vector<uint> pionTrackIDs;
  
  const ParameterSet config;
  const PixelCPEBase* cpe;
  
  void TestClusters(edm::Handle<edmNew::DetSetVector<SiPixelCluster>> inputhandle, edm::ESHandle<TrackerGeometry> & geom);

  vector<tuple<uint, uint, unsigned long long>> pickedEvents;
};

CharginoAnalyzer::CharginoAnalyzer(const ParameterSet& iConfig) :

genParticlesToken(            consumes<vector<GenParticle>>(InputTag("genParticles"))),
simTracksToken(               consumes<vector<SimTrack>>(   InputTag("g4SimHits"))),
simHitsPixelHighToken(        consumes<vector<PSimHit>>(    InputTag("g4SimHits", "TrackerHitsPixelBarrelHighTof"))),
simHitsPixelLowToken(         consumes<vector<PSimHit>>(    InputTag("g4SimHits", "TrackerHitsPixelBarrelLowTof"))),
simHitsPixelEndcapLowToken(   consumes<vector<PSimHit>>(    InputTag("g4SimHits", "TrackerHitsPixelEndcapLowTof"))),
simHitsPixelEndcapHighToken(  consumes<vector<PSimHit>>(	  InputTag("g4SimHits", "TrackerHitsPixelEndcapHighTof"))),
simHitsTEClowToken(           consumes<vector<PSimHit>>(    InputTag("g4SimHits", "TrackerHitsTECLowTof"))),
simHitsTEChighToken(          consumes<vector<PSimHit>>(    InputTag("g4SimHits", "TrackerHitsTECHighTof"))),
simHitsTIBlowToken(           consumes<vector<PSimHit>>(    InputTag("g4SimHits", "TrackerHitsTIBLowTof"))),
simHitsTIBhighToken(          consumes<vector<PSimHit>>(    InputTag("g4SimHits", "TrackerHitsTIBHighTof"))),
simHitsTIDlowToken(           consumes<vector<PSimHit>>(    InputTag("g4SimHits", "TrackerHitsTIDLowTof"))),
simHitsTIDhighToken(          consumes<vector<PSimHit>>(    InputTag("g4SimHits", "TrackerHitsTIDHighTof"))),
simHitsTOBhighToken(          consumes<vector<PSimHit>>(    InputTag("g4SimHits", "TrackerHitsTOBHighTof"))),
simHitsTOBlowToken(           consumes<vector<PSimHit>>(    InputTag("g4SimHits", "TrackerHitsTOBLowTof"))),

trackingRecHitsToken(         consumes<OwnVector<TrackingRecHit,ClonePolicy<TrackingRecHit>>>(InputTag("generalTracks"))),
recTracksToken(               consumes<vector<Track>>(InputTag("generalTracks"))),
trackingParticlesToken(       consumes<vector<TrackingParticle>>(InputTag("mix", "MergedTrackTruth"))),
pixelClusterToken(            consumes<edmNew::DetSetVector<SiPixelCluster>>(InputTag("siPixelClusters"))),
stripClusterToken(            consumes<edmNew::DetSetVector<SiStripCluster>>(InputTag("siStripClusters"))),
generalTracksToken(           consumes<std::vector<reco::Track>>(InputTag("generalTracks"))),

trackerHitAssociatorConfig( iConfig.getParameter<ParameterSet>("ClusterRefiner"), consumesCollector()),
config(iConfig)
{
  verbose = false;
  
  edm::Service<TFileService> fs;
  
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
  
  outputTree->Branch("runNumber", &run);
  outputTree->Branch("lumiBlock", &lumi);
  outputTree->Branch("eventNumber", &event);
  
  
  ifstream infile("/afs/cern.ch/work/j/jniedzie/private/disapp_tracks/friendTreeProducer/survivingSignalEventsAfterL1_all_500_10.txt");

  uint l, r;
  unsigned long long e;
  char junk;
  
  cout<<"Loading picked event numbers from file"<<endl;
  
//  cout<<infile.rdbuf()<<endl;
  
  while(infile >> r >> junk >> l >> junk >> e){
//    cout<<"run: "<<run<<"\tlumi: "<<lumi<<"\tevent: "<<event<<endl;
    pickedEvents.push_back({r, l, e});
  }
}

CharginoAnalyzer::~CharginoAnalyzer()
{
  
}

void CharginoAnalyzer::FillHitsForTrackID(const PSimHit &hit, uint trackID, bool isPixel, bool isPion, bool isChargino)
{
  if(hit.trackId() != trackID) return;
  if(hit.detUnitId() == 0) return; // but not sure why sometimes it's zero...
  
  GlobalPoint globalPoint;
  GeomDetEnumerators::SubDetector subDet;
  
  if(isPixel){
    auto detUnit = dynamic_cast<const PixelGeomDetUnit*>(trackerGeometry->idToDetUnit(hit.detUnitId()));
    subDet = detUnit->specificType().subDetector();
  }
  else{
    auto detUnit = dynamic_cast<const StripGeomDetUnit*>(trackerGeometry->idToDetUnit(hit.detUnitId()));
    subDet = detUnit->specificType().subDetector();
  }
  
  globalPoint = trackerGeometry->idToDet(hit.detUnitId())->surface().toGlobal(hit.localPosition());
  
  if(verbose) cout<<"\t\t\t{"<<globalPoint.x()<<", "<<globalPoint.y()<<", "<<globalPoint.z()<<"},"<<endl;
  
  if(isPion){
    pion_simHits_x.push_back(globalPoint.x());
    pion_simHits_y.push_back(globalPoint.y());
    pion_simHits_z.push_back(globalPoint.z());
    pion_simHits_t.push_back(hit.timeOfFlight());
    pion_simHits_subDet.push_back((int)subDet);
  }
  if(isChargino){
    chargino_simHits_x.push_back(globalPoint.x());
    chargino_simHits_y.push_back(globalPoint.y());
    chargino_simHits_z.push_back(globalPoint.z());
    chargino_simHits_subDet.push_back((int)subDet);
    nCharginoHits++;
  }
}

void CharginoAnalyzer::FillSimHitsForTrack(uint trackID, bool isPion, bool isChargino)
{
  // For each sim track, print it's sim hits
  if(verbose) cout<<"\t\tPXB low hits:"<<endl;
  for(auto hit : *simHitsPixelLow){FillHitsForTrackID(hit, trackID, true, isPion, isChargino);}
  
  if(verbose) cout<<"\t\tPXB high hits:"<<endl;
  for(auto hit : *simHitsPixelHigh){FillHitsForTrackID(hit, trackID, true, isPion, isChargino);}
  
  if(verbose) cout<<"\t\tPXE low hits:"<<endl;
  for(auto hit : *simHitsPixelEndcapLow){FillHitsForTrackID(hit, trackID, true, isPion, isChargino);}
  if(verbose) cout<<"\t\tPXE high hits:"<<endl;
  for(auto hit : *simHitsPixelEndcapHigh){FillHitsForTrackID(hit, trackID, true, isPion, isChargino);}
  
  if(verbose) cout<<"\t\tTIB low hits:"<<endl;
  for(auto hit : *simHitsTIBlow){FillHitsForTrackID(hit, trackID, false, isPion, isChargino);}
  if(verbose) cout<<"\t\tTIB high hits:"<<endl;
  for(auto hit : *simHitsTIBhigh){FillHitsForTrackID(hit, trackID, false, isPion, isChargino);}
  
  if(verbose) cout<<"\t\tTOB low hits:"<<endl;
  for(auto hit : *simHitsTOBlow){FillHitsForTrackID(hit, trackID, false, isPion, isChargino);}
  if(verbose) cout<<"\t\tTOB high hits:"<<endl;
  for(auto hit : *simHitsTOBhigh){FillHitsForTrackID(hit, trackID, false, isPion, isChargino);}
  
  if(verbose) cout<<"\t\tTEC low hits:"<<endl;
  for(auto hit : *simHitsTEClow){FillHitsForTrackID(hit, trackID, false, isPion, isChargino);}
  if(verbose) cout<<"\t\tTEC high hits:"<<endl;
  for(auto hit : *simHitsTEChigh){FillHitsForTrackID(hit, trackID, false, isPion, isChargino);}
  
  if(verbose) cout<<"\t\tTID low hits:"<<endl;
  for(auto hit : *simHitsTIDlow){FillHitsForTrackID(hit, trackID, false, isPion, isChargino);}
  if(verbose) cout<<"\t\tTID high hits:"<<endl;
  for(auto hit : *simHitsTIDhigh){FillHitsForTrackID(hit, trackID, false, isPion, isChargino);}
}

void CharginoAnalyzer::print(const TrackingParticle *particle)
{
  // Print tracking particle properties
  cout<<"------------------------------------------------------------------"<<endl;
  cout<<"\tPID:"<<particle->pdgId()<<"\tstatus:"<<particle->status()<<endl;
  cout<<"\tCharge:"<<particle->charge()<<endl;
  cout<<"\tEnergy:"<<particle->energy()<<"\tEt:"<<particle->et()<<"\tmass:"<<particle->mass()<<endl;
  cout<<"\tPt:"<<particle->pt()<<"\tpz:"<<particle->pz()<<endl;
  cout<<"\tEta:"<<particle->eta()<<"\tphi:"<<particle->phi()<<endl;
  cout<<"\tLong lived:"<<particle->longLived()<<endl;
  cout<<"\tNumber of hits:"<<particle->numberOfHits()<<"\ttracker hits:"<<particle->numberOfTrackerHits()<<"\tlayers:"<<particle->numberOfTrackerLayers()<<endl;
  cout<<"\tParent vertex:"<<particle->vx()<<", "<<particle->vy()<<", "<<particle->vz()<<endl;
  
  // Print vertices
  const TrackingVertexRefVector &verticesRef = particle->decayVertices();
  
  cout<<"\tVertices:"<<endl;
  for(uint iVertex=0; iVertex<verticesRef.size(); iVertex++){
    const TrackingVertex *vertex = &(*verticesRef[iVertex]);
    cout<<"\t\t"<<vertex->position().x()<<", "<<vertex->position().y()<<", "<<vertex->position().z()<<endl;
  }
  
  // Print geant tracks
  const vector<SimTrack> &geantTracks = particle->g4Tracks();
  
  cout<<"\tSim tracks:"<<endl;
  for(uint iTrack=0; iTrack<geantTracks.size(); iTrack++){
    const SimTrack track = geantTracks[iTrack];
    uint trackID = track.trackId();
    
    cout<<"\t\tEta:"<<track.momentum().eta()<<"\tphi:"<<track.momentum().phi()<<"\tpt:"<<track.momentum().pt()<<endl;
    cout<<"\t\ttrackId:"<<trackID<<endl;
  }
  
  // Print gen particles
  const GenParticleRefVector &genParticles = particle->genParticles();
  
  cout<<"\tGen particles:"<<endl;
  for(uint iPart=0; iPart<genParticles.size(); iPart++){
    const GenParticle *part = &(*genParticles[iPart]);
    cout<<"\t\tEta:"<<part->momentum().eta()<<"\tphi:"<<part->momentum().phi()<<"\tpt:"<<sqrt(pow(part->momentum().x(),2)+pow(part->momentum().y(),2))<<endl;
  }
  cout<<"------------------------------------------------------------------"<<endl;
}

void CharginoAnalyzer::PrintTrackingParticles()
{
  vector<vector<double>> charginoDecayVertices;
  
  //------------------------------------------------------
  // Print charginos and neutralinos
  //------------------------------------------------------
  for(size_t i=0; i<trackingParticles->size(); i++){
    const TrackingParticle *particle = &(*trackingParticles)[i];
    
    if(abs(particle->pdgId()) == 1000022 || abs(particle->pdgId()) == 1000024){
      if(verbose){
        cout<<"Tracking particle "<<i<<endl;
        print(particle);
      }
      
      // Save decay vertices of charginos
      const TrackingVertexRefVector &verticesRef = particle->decayVertices();
      
      for(uint iVertex=0; iVertex<verticesRef.size(); iVertex++){
        const TrackingVertex *vertex = &(*verticesRef[iVertex]);
        
        vector<double> decayVertex = { vertex->position().x(), vertex->position().y(), vertex->position().z() };
        charginoDecayVertices.push_back(decayVertex);
      }
      
      // Get sim track of chargino and store sim hits assosiated with it
      const vector<SimTrack> &geantTracks = particle->g4Tracks();
      
      for(uint iTrack=0; iTrack<geantTracks.size(); iTrack++){
        const SimTrack track = geantTracks[iTrack];
        uint trackID = track.trackId();
        bool isPion = false;
        bool isChargino = true;
        
        nCharginoHits = 0;
        FillSimHitsForTrack(trackID, isPion, isChargino);
        
        double eta = track.momentum().eta();
        double phi = track.momentum().phi();
        double pt = track.momentum().Pt();
        double px = track.momentum().Px();
        double py = track.momentum().Py();
        double pz = track.momentum().Pz();
        int charge = track.charge();
        
        chargino_eta.push_back(eta);
        chargino_phi.push_back(phi);
        chargino_pt.push_back(pt);
        chargino_px.push_back(px);
        chargino_py.push_back(py);
        chargino_pz.push_back(pz);
        chargino_nTrackerLayers.push_back(nCharginoHits);
        chargino_charge.push_back(charge);
      }
    }
  }
  
  if(verbose){
    cout<<"\n\n chargino decay vertices:"<<endl;
    for(auto v : charginoDecayVertices){
      cout<<v[0]<<"\t"<<v[1]<<"\t"<<v[2]<<endl;
    }
  }
  
  //------------------------------------------------------
  // Print pions from decays
  //------------------------------------------------------
  for(size_t i=0; i<trackingParticles->size(); i++){
    const TrackingParticle *particle = &(*trackingParticles)[i];
    
    if(abs(particle->pdgId()) == 211){
      for(auto v : charginoDecayVertices){
        if(particle->vx() == v[0] &&
           particle->vy() == v[1] &&
           particle->vz() == v[2] ){
          
          if(verbose){
            cout<<"Tracking particle "<<i<<endl;
            print(particle);
          }
          
          pion_vx.push_back(particle->vx());
          pion_vy.push_back(particle->vy());
          pion_vz.push_back(particle->vz());
          pion_px.push_back(particle->px());
          pion_py.push_back(particle->py());
          pion_pz.push_back(particle->pz());
          pion_charge.push_back(particle->charge());
          
          // Get sim track of pion and store sim hits assosiated with it
          const vector<SimTrack> &geantTracks = particle->g4Tracks();
          
          for(uint iTrack=0; iTrack<geantTracks.size(); iTrack++){
            const SimTrack track = geantTracks[iTrack];
            uint trackID = track.trackId();
            bool isPion = true;
            bool isChargino = false;
            
            FillSimHitsForTrack(trackID, isPion, isChargino);
            pionTrackIDs.push_back(trackID);
          }
          
          break;
        }
      }
    }
  }
}

void CharginoAnalyzer::FillTrackerClusters(bool filterTrackClusters)
{
  //------------------------------------------------------
  // Find IDs of clusters that belong to some tracks
  //------------------------------------------------------
  const std::vector<reco::Track> *tracks = generalTracksHandle.product();
  
  vector<uint> trackClustersIDs;
  
  for(uint iTrack=0;iTrack<tracks->size();iTrack++){
    Track track = tracks->at(iTrack);
    
    for(uint iHit=0;iHit<track.recHitsSize();iHit++){
      const TrackingRecHit *hit = track.recHit(iHit).get();
      trackClustersIDs.push_back(hit->rawId());
    }
  }
  
  //------------------------------------------------------
  // Fill in the strip clusters
  //------------------------------------------------------
  const edmNew::DetSetVector<SiStripCluster> *stripClusters = stripClustersHandle.product();
  
  for (auto clusterSet = stripClusters->begin(); clusterSet!=stripClusters->end(); ++clusterSet) {
    
    DetId detId(clusterSet->detId());
    
    auto detUnit = dynamic_cast<const StripGeomDetUnit*>(trackerGeometry->idToDet(detId));
    auto stripTopology = dynamic_cast<const StripTopology*>(&(detUnit->specificTopology()));
    
    for(auto cluster = clusterSet->begin(); cluster != clusterSet->end(); cluster++){
      
      LocalPoint localPoint   = stripTopology->localPosition(cluster->barycenter());
      GlobalPoint globalPoint = detUnit->surface().toGlobal(localPoint);
      GeomDetEnumerators::SubDetector subDet = detUnit->specificType().subDetector();
      
      double stripLength = stripTopology->localStripLength(localPoint);
      LocalError localError = LocalError(pow(stripTopology->localPitch(localPoint), 2)/12,
                                         0,
                                         pow(stripTopology->localPitch(localPoint), 2)/12);
      
      stripCluster_x.push_back(globalPoint.x());
      stripCluster_y.push_back(globalPoint.y());
      stripCluster_z.push_back(globalPoint.z());
      stripCluster_ex.push_back(sqrt(localError.xx()));
      stripCluster_ey.push_back(sqrt(localError.yy()));
      stripCluster_ez.push_back(stripLength/2.);
      stripCluster_charge.push_back(cluster->charge());
      stripCluster_subDet.push_back((int)subDet);
      
      vector<SimHitIdpr> simtrackid;
      vector<PSimHit> simhit;
      
      trackerHitAssociator->associateCluster(cluster, detId, simtrackid, simhit);
      
      for(auto stripSimTrackId : simtrackid){
        bool found = false;
        for(auto pionTrackID : pionTrackIDs){
          if(stripSimTrackId.first == pionTrackID){
            found = true;
            break;
          }
        }
        if(found){
          pionCluster_x.push_back(globalPoint.x());
          pionCluster_y.push_back(globalPoint.y());
          pionCluster_z.push_back(globalPoint.z());
          pionCluster_ex.push_back(sqrt(localError.xx()));
          pionCluster_ey.push_back(sqrt(localError.yy()));
          pionCluster_ez.push_back(stripLength/2.);
          pionCluster_charge.push_back(cluster->charge());
          pionCluster_subDet.push_back((int)subDet);
          
          break;
        }
      }
    }
  }
  
  //------------------------------------------------------
  // Fill in the pixel clusters
  //------------------------------------------------------
  const edmNew::DetSetVector<SiPixelCluster> *pixelClusters = pixelClustersHandle.product();
  
  for (auto clusterSet = pixelClusters->begin(); clusterSet!=pixelClusters->end(); ++clusterSet) {
    
    DetId detId(clusterSet->detId());
    
    auto detUnit = dynamic_cast<const PixelGeomDetUnit*>(trackerGeometry->idToDet(detId));
    auto pixelTopology = dynamic_cast<const PixelTopology*>(&(detUnit->specificTopology()));
    //
    const GeomDetUnit *genericDet = trackerGeometry->idToDetUnit(detId);
    auto pixDet = dynamic_cast<const PixelGeomDetUnit*>(genericDet);
    //
    
    for(auto cluster = clusterSet->begin(); cluster != clusterSet->end(); cluster++){
      
      LocalPoint localPoint = pixelTopology->localPosition(MeasurementPoint(cluster->x(),cluster->y()));
      GlobalPoint globalPoint = detUnit->surface().toGlobal(localPoint);
      GeomDetEnumerators::SubDetector subDet = detUnit->specificType().subDetector();
      
      pixelCluster_x.push_back(globalPoint.x());
      pixelCluster_y.push_back(globalPoint.y());
      pixelCluster_z.push_back(globalPoint.z());
      pixelCluster_charge.push_back(cluster->charge());
      pixelCluster_subDet.push_back((int)subDet);
    
      /*
      //
      tuple<LocalPoint, LocalError, SiPixelRecHitQuality::QualWordType> tuple = cpe->getParameters(*cluster, *genericDet);
      LocalPoint lp( std::get<0>(tuple) );
      LocalError le( std::get<1>(tuple) );
      SiPixelRecHitQuality::QualWordType rqw( std::get<2>(tuple) );
      
      Ref< edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster > clusterRef = edmNew::makeRefTo(pixelClustersHandle, cluster);
      
      const SiPixelRecHit *hit = new SiPixelRecHit( lp, le, rqw, *genericDet, clusterRef);
      if(!hit){
        cout<<"No pixel rec hit could be created"<<endl;
        continue;
      }
      //
      
      vector<SimHitIdpr> simtrackid;
      trackerHitAssociator->associatePixelRecHit(hit, simtrackid);
      */
    }
  }
}

void CharginoAnalyzer::TestClusters(edm::Handle<edmNew::DetSetVector<SiPixelCluster> >  inputhandle,
                           edm::ESHandle<TrackerGeometry> & geom) {
  if (!cpe){
    edm::LogError("SiPixelRecHitConverter") << " at least one CPE is not ready -- can't run!";
    // TO DO: throw an exception here?  The user may want to know...
    assert(0);
    return;   // clusterizer is invalid, bail out
  }
  
//  SiPixelRecHitCollectionNew output;
  
  int numberOfDetUnits = 0;
  int numberOfClusters = 0;
  
  const edmNew::DetSetVector<SiPixelCluster>& input = *inputhandle;
  
  edmNew::DetSetVector<SiPixelCluster>::const_iterator DSViter=input.begin();
  
  for ( ; DSViter != input.end() ; DSViter++) {
    numberOfDetUnits++;
    unsigned int detid = DSViter->detId();
    DetId detIdObject( detid );
    const GeomDetUnit * genericDet = geom->idToDetUnit( detIdObject );
    const PixelGeomDetUnit * pixDet = dynamic_cast<const PixelGeomDetUnit*>(genericDet);
    assert(pixDet);
//    SiPixelRecHitCollectionNew::FastFiller recHitsOnDetUnit(output,detid);
    
    edmNew::DetSet<SiPixelCluster>::const_iterator clustIt = DSViter->begin(), clustEnd = DSViter->end();
    
    for ( ; clustIt != clustEnd; clustIt++) {
      numberOfClusters++;
      std::tuple<LocalPoint, LocalError,SiPixelRecHitQuality::QualWordType> tuple = cpe->getParameters( *clustIt, *genericDet );
      LocalPoint lp( std::get<0>(tuple) );
      LocalError le( std::get<1>(tuple) );
      SiPixelRecHitQuality::QualWordType rqw( std::get<2>(tuple) );
      // Create a persistent edm::Ref to the cluster
      edm::Ref< edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster > cluster = edmNew::makeRefTo( inputhandle, clustIt);
      // Make a RecHit and add it to the DetSet
      // old : recHitsOnDetUnit.push_back( new SiPixelRecHit( lp, le, detIdObject, &*clustIt) );
      SiPixelRecHit *hit = new SiPixelRecHit( lp, le, rqw, *genericDet, cluster);
      //
      // Now save it =================
//      recHitsOnDetUnit.push_back(hit);
      // =============================

      vector<SimHitIdpr> simtrackid;
      vector<TrackerHitAssociator::simhitAddr> simhitCFPos;
      
      trackerHitAssociator->associatePixelRecHit(hit, simtrackid, &simhitCFPos);
      
      std::cout << "SiPixelRecHitConverterVI " << numberOfClusters << ' '<< lp << " " << le << std::endl;
    } //  <-- End loop on Clusters
    
    
    LogDebug("SiPixelRecHitConverter");
    std::cout << "SiPixelRecHitConverterVI " << " Found RecHits on " << detid << std::endl;
    
    
  } //    <-- End loop on DetUnits
  
  LogDebug ("SiPixelRecHitConverter");
  std::cout << "SiPixelRecHitConverterVI converted " << numberOfClusters << "SiPixelClusters into SiPixelRecHits, in " << numberOfDetUnits << " DetUnits." << std::endl;
}


// ------------ method called for each event  ------------
void CharginoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  bool filterTrackClusters = true;
 
  ESHandle<PixelClusterParameterEstimator> hCPE;
  iSetup.get<TkPixelCPERecord>().get(config.getParameter<string>("CPE"), hCPE);
  cpe = dynamic_cast<const PixelCPEBase*>(&(*hCPE));
  
  run   = iEvent.id().run();
  lumi  = iEvent.id().luminosityBlock();
  event = iEvent.id().event();
  
  tuple<uint, uint, unsigned long long> eventTuple = {run, lumi, event};
  
  if(find(pickedEvents.begin(), pickedEvents.end(), eventTuple) == pickedEvents.end()){
    cout<<"Event not in the picked events list. Skipping..."<<endl;
    return;
  }
  
  //  if(event!=2331) return;
  
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
  
  
  if(verbose){
    cout<<"\n\n================================================================="<<endl;
    cout<<"Event:"<<event<<"\trun:"<<run<<"\tlumi:"<<lumi<<endl;
    cout<<"\n"<<endl;
  }
  
  iEvent.getByToken(genParticlesToken,          genParticles);
  iEvent.getByToken(simTracksToken,             simTracks);
  iEvent.getByToken(simHitsPixelHighToken,      simHitsPixelHigh);
  iEvent.getByToken(simHitsPixelLowToken,       simHitsPixelLow);
  iEvent.getByToken(simHitsPixelEndcapLowToken, simHitsPixelEndcapLow);
  iEvent.getByToken(simHitsPixelEndcapHighToken,simHitsPixelEndcapHigh);
  iEvent.getByToken(simHitsTEClowToken,         simHitsTEClow);
  iEvent.getByToken(simHitsTEChighToken,        simHitsTEChigh);
  iEvent.getByToken(simHitsTIBlowToken,         simHitsTIBlow);
  iEvent.getByToken(simHitsTIBhighToken,        simHitsTIBhigh);
  iEvent.getByToken(simHitsTIDlowToken,         simHitsTIDlow);
  iEvent.getByToken(simHitsTIDhighToken, 	      simHitsTIDhigh);
  iEvent.getByToken(simHitsTOBhighToken,        simHitsTOBhigh);
  iEvent.getByToken(simHitsTOBlowToken,         simHitsTOBlow);
  
  iEvent.getByToken(trackingRecHitsToken,       trackingRecHitsRef);
  iEvent.getByToken(recTracksToken,             recTracks);
  iEvent.getByToken(trackingParticlesToken,     trackingParticles);
  iEvent.getByToken(generalTracksToken,         generalTracksHandle);
  iEvent.getByToken(pixelClusterToken,          pixelClustersHandle);
  iEvent.getByToken(stripClusterToken,          stripClustersHandle);
  
  iSetup.get<TrackerDigiGeometryRecord>().get(trackerGeometry);
  
  trackerHitAssociator = new TrackerHitAssociator(iEvent, trackerHitAssociatorConfig);
  
  // pion id:       211
  // chargino id:   1000024
  // neutralino id: 1000022
  
  // Print gen particles (coming from primary chargino/neutralino)
  /*
   for(size_t i=0; i<genParticles->size(); i++){
   const GenParticle *p = &(*genParticles)[i];
   const Candidate *mom = p->mother();
   if(!mom) continue;
   
   int id = p->pdgId();
   int momId = mom->pdgId();
   
   // check if this is the first chargino/neutralino in the chain
   if(abs(id)!=1000022 && abs(id)!=1000024) continue;
   if(abs(momId) == 1000022 || abs(momId) == 1000024) continue;
   
   if(verbose){
   cout<<"\n\nParticle chain:\n\n"<<endl;
   printDeeply(p);
   }
   }
   */
  PrintTrackingParticles();
  FillTrackerClusters(filterTrackClusters);
  
//  TestClusters(pixelClustersHandle, trackerGeometry);
  
  outputTree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void
CharginoAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
CharginoAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CharginoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CharginoAnalyzer);
