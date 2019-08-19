/**
 * \class TrkEval
 * \brief General purposed evaluation module
 * \author Haiwang Yu, yuhw@nmsu.edu
 *
 * Created: 08-27-2018
 *
 *
 */


#include "TrkEval.h"

#include <interface_main/SQHit.h>
#include <interface_main/SQHit_v1.h>
#include <interface_main/SQMCHit_v1.h>
#include <interface_main/SQHitMap_v1.h>
#include <interface_main/SQHitVector_v1.h>
#include <interface_main/SQEvent_v1.h>
#include <interface_main/SQRun_v1.h>
#include <interface_main/SQSpill_v1.h>
#include <interface_main/SQSpillMap_v1.h>

#include <ktracker/SRecEvent.h>

#include <geom_svc/GeomSvc.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4HitDefs.h>
#include <g4main/PHG4VtxPoint.h>

#include <TFile.h>
#include <TTree.h>

#include <cstring>
#include <cmath>
#include <cfloat>
#include <stdexcept>
#include <limits>
#include <tuple>
#include <typeinfo>

#include <boost/lexical_cast.hpp>

#define LogDebug(exp)		std::cout<<"DEBUG: "  <<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl
#define LogError(exp)		std::cout<<"ERROR: "  <<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl
#define LogWarning(exp)	    std::cout<<"WARNING: "<<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl

using namespace std;

TrkEval::TrkEval(const std::string& name) :
SubsysReco(name),
_hit_container_type("Vector"),
_event(0),
_run_header(nullptr),
_spill_map(nullptr),
_event_header(nullptr),
_hit_map(nullptr),
_hit_vector(nullptr),
_out_name("eval.root")
{
}

int TrkEval::Init(PHCompositeNode* topNode) {
	return Fun4AllReturnCodes::EVENT_OK;
}

int TrkEval::InitRun(PHCompositeNode* topNode) {

	ResetEvalVars();
	InitEvalTree();

	p_geomSvc = GeomSvc::instance();

	int ret = GetNodes(topNode);
	if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

	return Fun4AllReturnCodes::EVENT_OK;
}

int TrkEval::process_event(PHCompositeNode* topNode) {
	int ret = Fun4AllReturnCodes::ABORTRUN;

	ret = TruthEval(topNode);
	if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;
	/*
	if(_recEvent) {
		ret = RecoEval(topNode);
		if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;
	}
	*/
	++_event;

	return ret;
}

int TrkEval::TruthEval(PHCompositeNode* topNode)
{
        //Incase it is not clear: these eval functions run per event. All the limits only limit
        //how many hits/tracks are in 1 single event. 

	if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT)
		std::cout << "Entering TrkEval::TruthEval: " << _event << std::endl;

	ResetEvalVars();

	if(_spill_map) {
	  auto spill_info = _spill_map->get(spill_id);
	  if(spill_info) {
	    target_pos = spill_info->get_target_pos();
	  } else {
	    LogWarning("");
	  }
	}
	
	if(_event_header) {
	  event_id    = _event_header->get_event_id();    
	  emu_trigger = _event_header->get_trigger();
	  spill_id    = _event_header->get_spill_id();
	  run_id      = _event_header->get_run_id();
	}

	std::map<int, int> parID_nhits_dc;
	std::map<int, int> parID_nhits_hodo;
	std::map<int, int> parID_nhits_prop;

	//add info for the number of dark photon detector hits
	std::map<int, int> parID_nhits_dp;
	std::map<int, int> parID_nhits_H4Y;
	std::map<int, std::map<int, int> > parID_detid_elmid;

	typedef std::tuple<int, int> ParDetPair;
	std::map<ParDetPair, int> parID_detID_ihit;
	std::map<int, int> hitID_ihit;

	if(_hit_vector) {
	  n_hits = 0;
	  for(int ihit=0; ihit<_hit_vector->size(); ++ihit) {
	    SQHit *hit = _hit_vector->at(ihit);
	    
	    if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT) {
	      LogInfo(hit->get_detector_id());
	      hit->identify();
	    }
	    
	    int hitID = hit->get_hit_id();
	    hit_id[n_hits]         = hitID;
	    hitID_ihit[hitID]      = ihit;
	    drift_distance[n_hits] = hit->get_drift_distance();
	    pos[n_hits]            = hit->get_pos();
	    detector_z[n_hits]     = p_geomSvc->getPlanePosition(hit->get_detector_id());
	    detector_id[n_hits]    = hit->get_detector_id();
	    element_id[n_hits]     = hit->get_element_id();
	    hodo_mask[n_hits]      = hit->is_hodo_mask();

	    if(_truth) {
	      int track_id = hit->get_track_id();
	      int det_id = hit->get_detector_id();
	      
	      //NOTE HERE:
	      //parID is not the particle ID, but the trackID.
	      //The parID indexes starting from 1, if you are looking at info later on,
	      //anything indexed by n_tracks starts at 0. Make sure you pay attention
	      // what is indexed where.
	      
	      parID_detID_ihit[std::make_tuple(track_id, det_id)] = ihit;

	      auto detid_elmid_iter = parID_detid_elmid.find(track_id);
	      if(detid_elmid_iter != parID_detid_elmid.end()) {
      		detid_elmid_iter->second.insert(std::pair<int, int>(det_id, hit->get_element_id()));
	      } else {
      		std::map<int, int> detid_elmid;
      		detid_elmid.insert(std::pair<int, int>(det_id, hit->get_element_id()));
      		parID_detid_elmid[track_id] = detid_elmid;
	      }
	      
	      // dc hits
	      if(hit->get_detector_id() >= 1 and hit->get_detector_id() <=30) {
        	if(parID_nhits_dc.find(track_id)!=parID_nhits_dc.end())
		  parID_nhits_dc[track_id] = parID_nhits_dc[track_id]+1;
        	else
		  parID_nhits_dc[track_id] = 1;
	      }
	      // hodo hits
	      if(hit->get_detector_id() >= 31 and hit->get_detector_id() <=46) {
		if(parID_nhits_hodo.find(track_id)!=parID_nhits_hodo.end())
		  parID_nhits_hodo[track_id] = parID_nhits_hodo[track_id]+1;
        	else
		  parID_nhits_hodo[track_id] = 1;
	      }
	      // prop tube hits
	      if(hit->get_detector_id() >= 47 and hit->get_detector_id() <=54) {
        	if(parID_nhits_prop.find(track_id)!=parID_nhits_prop.end())
		  parID_nhits_prop[track_id] = parID_nhits_prop[track_id]+1;
        	else
		  parID_nhits_prop[track_id] = 1;
	      }

	      // dp hits
	      if((hit->get_detector_id() >= 55 and hit->get_detector_id() <=62)) {
		if(parID_nhits_dp.find(track_id)!=parID_nhits_dp.end())
		  parID_nhits_dp[track_id] = parID_nhits_dp[track_id]+1;
		else
		  parID_nhits_dp[track_id] = 1;
	      }


	      if((hit->get_detector_id() == 43 or hit->get_detector_id() ==44)) {
		//it is possible that a hit will hit both det 43 and 44 since they over lap,
		//but we are only interested in if it hits 1, so as long as it hits, take that 
		// to be 1 hit.
		parID_nhits_H4Y[track_id] = 1;
	      }
	      
	      
	      //NOTE:
	      //this is not the PHG4 truth information ("true truth") information.
	      //SQHit truth information =/= the information that you give the simulation.
	      //For example, if the information you give says hit 5 should be at 1,1,1 
	      //the "truth" information from SQHit might not give 1,1,1 below.
	      //
	      //PHG4Hit is what has the true truth information.
	      truth_x[n_hits] = hit->get_truth_x();
	      truth_y[n_hits] = hit->get_truth_y();
	      truth_z[n_hits] = hit->get_truth_z();
	      
	      double uVec[3] = {
		p_geomSvc->getPlane(hit->get_detector_id()).uVec[0],
		p_geomSvc->getPlane(hit->get_detector_id()).uVec[1],
		p_geomSvc->getPlane(hit->get_detector_id()).uVec[2]
	      };
	      
	      truth_pos[n_hits] =
		//      			(truth_x[n_hits] - p_geomSvc->getPlane(hit->get_detector_id()).xc)*uVec[0] +
		//      			(truth_y[n_hits] - p_geomSvc->getPlane(hit->get_detector_id()).yc)*uVec[1] +
		//      			(truth_z[n_hits] - p_geomSvc->getPlane(hit->get_detector_id()).zc)*uVec[2];
		(truth_x[n_hits])*uVec[0] +
		(truth_y[n_hits])*uVec[1] +
		(truth_z[n_hits]-p_geomSvc->getPlane(hit->get_detector_id()).zc)*uVec[2];
	      
	      if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT) {
        	LogInfo("");
        	std::cout << truth_pos[n_hits] << " => { "
			  << truth_x[n_hits] << ", "
			  << truth_y[n_hits] << ", "
			  << truth_z[n_hits] << "} {"
			  << uVec[0] << ", "
			  << uVec[1] << ", "
			  << uVec[2] << "}"
			  << std::endl;
	      }
	      
	      //LogDebug("detector_id: " << detector_id[n_hits]);
	    }
	    ++n_hits;
	    if(n_hits>=10000) break;
	  }
	}

	if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT) LogInfo("ghit eval finished");
	
	if(_truth) {
	  for(auto iter=_truth->GetPrimaryParticleRange().first; iter!=_truth->GetPrimaryParticleRange().second; ++iter) {
	
	    std::cout << "looping over phg4particle " << std::endl;
	    //PHG4Particle is the particle information from geant. This has all the inforation
	    //that you put into the simulation and all the information the geant created when
	    //propogating the particle.
	    PHG4Particle * par = iter->second;
	    
	    pid[n_tracks] = par->get_pid();
	    int vtx_id =  par->get_vtx_id();
	    PHG4VtxPoint* vtx = _truth->GetVtx(vtx_id);
	    gvx[n_tracks] = vtx->get_x();
	    gvy[n_tracks] = vtx->get_y();
	    gvz[n_tracks] = vtx->get_z();
	    
	    TVector3 mom(par->get_px(), par->get_py(), par->get_pz());
	    gpx[n_tracks] = par->get_px();
	    gpy[n_tracks] = par->get_py();
	    gpz[n_tracks] = par->get_pz();
	    gpt[n_tracks] = mom.Pt();
	    geta[n_tracks] = mom.Eta();
	    gphi[n_tracks] = mom.Phi();
	    
	    //Not particle ID, trackID
	    int parID = par->get_track_id();
	    par_id[n_tracks] = parID;
	    
	    //The next few blocks are grabbing information at certain stations. For this module
	    //we get the information at stations 1,2 and 3, dark photon detector 1 and 2, and hodoscope
	    //H4Y2. The detector ID and names are listed in e1039-core/packages/geom_svc/GeomSvc.cxx.
	    
	    PHG4HitContainer *D1X_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D0X");
	    if (!D1X_hits)
	      D1X_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D1X");
	    
	    if (!D1X_hits)
	      {
		if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT)
		  cout << Name() << " Could not locate g4 hit node " << "G4HIT_D0X or G4HIT_D1X" << endl;
	      }
	    
	    // trackID + detID -> SQHit -> PHG4Hit -> momentum
	    //detID 1-6 deal with D0, 7-12 D1
	    for(int det_id=1; det_id<=12; ++det_id) {
	      auto iter = parID_detID_ihit.find(std::make_tuple(parID, det_id));
	      if(iter != parID_detID_ihit.end()) {
		SQHit *hit = _hit_vector->at(iter->second);
		if(hit and D1X_hits) {
		  PHG4Hit* g4hit =  D1X_hits->findHit(hit->get_g4hit_id());
		  if (g4hit) {
		    gx_st1[n_tracks]  = g4hit->get_x(0);
		    gy_st1[n_tracks]  = g4hit->get_y(0);
		    gz_st1[n_tracks]  = g4hit->get_z(0);
		    
		    gpx_st1[n_tracks] = g4hit->get_px(0)/1000.;
		    gpy_st1[n_tracks] = g4hit->get_py(0)/1000.;
		    gpz_st1[n_tracks] = g4hit->get_pz(0)/1000.;
		    
		    if(gpz_st1[n_tracks] <0){
		      std::cout << "WARNING:: Negative z-momentum at Station 1!" << std::endl;
		    }
		    break;
		  }
		}
	      }
	    }// end st1 det id loop  
	    if(verbosity>=2) std::cout << "station 1 truth info done." << std::endl;

	    PHG4HitContainer *D2X_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D2X");
	    //if (!D2X_hits)
	    //D1X_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D1X");
	    if(!D2X_hits){
	      D2X_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D2Xp");
	    }
	    if (!D2X_hits)
	      {
		if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT)
		  cout << Name() << " Could not locate g4 hit node " <<  "G4HIT_D2X" << endl;
	      }
	    
	    //detID 13-18 are D2
	    // trackID + detID -> SQHit -> PHG4Hit -> momentum
	    for(int det_id=13; det_id<=18; ++det_id) {
	      auto iter = parID_detID_ihit.find(std::make_tuple(parID, det_id));
	      if(iter != parID_detID_ihit.end()) {
		SQHit *hit = _hit_vector->at(iter->second);
		if(hit and D2X_hits) {
		  PHG4Hit* g4hit =  D2X_hits->findHit(hit->get_g4hit_id());
		  if (g4hit) {
		    gx_st2[n_tracks]  = g4hit->get_x(0);
		    gy_st2[n_tracks]  = g4hit->get_y(0);
		    gz_st2[n_tracks]  = g4hit->get_z(0);
		    
		    gpx_st2[n_tracks] = g4hit->get_px(0)/1000.;
		    gpy_st2[n_tracks] = g4hit->get_py(0)/1000.;
		    gpz_st2[n_tracks] = g4hit->get_pz(0)/1000.;
		    if(gpz_st2[n_tracks] <0){
		      std::cout << "WARNING:: Negative z-momentum at Station 2!" << std::endl;
		    }
		    break;
		  }
		}
	      }
	    }// end st2 det id loop  
	    if(verbosity>=2) std::cout << "station 2 truth info done." << std::endl;
	    
	    PHG4HitContainer *D3X_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D3pX");
	    
	    if (!D3X_hits)
	      {
		D3X_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D3pXp");
	      }
	    
	    
	    //detID 19-24 are D3p, 25-30 are D3m
	    for(int det_id=19; det_id<=24; ++det_id) {
	      auto iter = parID_detID_ihit.find(std::make_tuple(parID, det_id));
	      if(iter != parID_detID_ihit.end()) {
		SQHit *hit = _hit_vector->at(iter->second);
		if(hit and D3X_hits) {
		  PHG4Hit* g4hit =  D3X_hits->findHit(hit->get_g4hit_id());
		  if (g4hit) {
		    gx_st3[n_tracks]  = g4hit->get_x(0);
		    gy_st3[n_tracks]  = g4hit->get_y(0);
		    gz_st3[n_tracks]  = g4hit->get_z(0);
		    
		    gpx_st3[n_tracks] = g4hit->get_px(0)/1000.;
		    gpy_st3[n_tracks] = g4hit->get_py(0)/1000.;
		    gpz_st3[n_tracks] = g4hit->get_pz(0)/1000.;
		    if(gpz_st3[n_tracks] <0){
		      std::cout << "WARNING:: Negative z-momentum at Station 3!" << std::endl;
		    }
		    
		    break;
		  }
		}
	      }
	    } // end st3 det id loop

	    D3X_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D3mX"); 
	    if (!D3X_hits)
	      {
		D3X_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_D3mXp");
	      }
	    
	    if(!D3X_hits){
	      std::cout << "Could not locate D3X container." << std::endl;
	    }
	    
	    //detID  25-30 are D3m
	    for(int det_id=25; det_id<=30; ++det_id) {
	      auto iter = parID_detID_ihit.find(std::make_tuple(parID, det_id));
	      if(iter != parID_detID_ihit.end()) {
		SQHit *hit = _hit_vector->at(iter->second);
		
		if(hit and D3X_hits) {
		  PHG4Hit* g4hit =  D3X_hits->findHit(hit->get_g4hit_id());
		  if (g4hit) {
		    gx_st3[n_tracks]  = g4hit->get_x(0);
		    gy_st3[n_tracks]  = g4hit->get_y(0);
		    gz_st3[n_tracks]  = g4hit->get_z(0);
		    
		    gpx_st3[n_tracks] = g4hit->get_px(0)/1000.;
		    gpy_st3[n_tracks] = g4hit->get_py(0)/1000.;
		    gpz_st3[n_tracks] = g4hit->get_pz(0)/1000.;
		    if(gpz_st3[n_tracks] <0){
		      std::cout << "WARNING:: Negative z-momentum at Station 3!" << std::endl;
		    }
		    
		    break;
		  }
		}
	      }
	    }
	    if(verbosity>=2) std::cout << "station 3 truth info done." << std::endl;
	    
	    //detID 33 and 34 are for H1L/R. This is to check possible roads, if comment
	    //out if you do not want this info.
	    PHG4HitContainer *H1_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H1LL");
	    if(!H1_hits){
	      H1_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H1R");
	    }
	    
	    for(int det_id=33; det_id<=34; ++det_id) {
	      auto iter = parID_detID_ihit.find(std::make_tuple(parID, det_id));
	      if(iter != parID_detID_ihit.end()) {
		SQHit *hit = _hit_vector->at(iter->second);
		if(hit and H1_hits) {
		  PHG4Hit* g4hit =  H1_hits->findHit(hit->get_g4hit_id());
		  if (g4hit) {
		    int h1barID = hit->get_element_id();
		    if(h1barID<11) H1BarID[n_tracks] = hit->get_element_id();
		    if(h1barID>10) H1BarID[n_tracks] = hit->get_element_id() - 8;
		    break;
		  }
		}
	      }
	    }
	    
	    //detID 35 and 36 are for H2L/R. This is to check possible roads, if comment
	    //out if you do not want this info.
	    PHG4HitContainer *H2_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H2LL");
	    if(!H2_hits){
	      H2_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H2R");
	    }
	    
	    for(int det_id=35; det_id<=36; ++det_id) {
	      auto iter = parID_detID_ihit.find(std::make_tuple(parID, det_id));
	      if(iter != parID_detID_ihit.end()) {
		SQHit *hit = _hit_vector->at(iter->second);
		if(hit and H2_hits) {
		  PHG4Hit* g4hit =  H2_hits->findHit(hit->get_g4hit_id());
		  if (g4hit) {
		    int h2barID = hit->get_element_id();
		    if(h2barID<11) H2BarID[n_tracks] = hit->get_element_id();
		    if(h2barID>9) H2BarID[n_tracks] = hit->get_element_id() - 8;
		    break;
		  }
		}
	      }
	    }

	    //detID 43 and 44 are for H4Y2, which is the one used in DP tracking.
		
	    PHG4HitContainer *H4Y_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H4Y2L");
            
	    for(int det_id=43; det_id<44; ++det_id) {
	      auto iter = parID_detID_ihit.find(std::make_tuple(parID, det_id));
	      if(iter != parID_detID_ihit.end()) {
		SQHit *hit = _hit_vector->at(iter->second);
		if(hit and H4Y_hits) {
		  if(verbosity >= Fun4AllBase::VERBOSITY_A_LOT) {
		    LogDebug("h4y2lhit: " << iter->second);
		  }
		  PHG4Hit* g4hit =  H4Y_hits->findHit(hit->get_g4hit_id());
		  if (g4hit) {
		    //if(verbosity >= 2) {
		    //  std::cout << "g4hit h4y2l " << std::endl;
		    //  g4hit->identify();
		    //}
		    int h4ybarID = hit->get_element_id();
		    
		    if(h4ybarID<9) H4YBarID[n_tracks] = 9 - hit->get_element_id();
		    if(h4ybarID>8) H4YBarID[n_tracks] = hit->get_element_id() - 8;
		    //quadH4Y[n_tracks] = (-1)*hit->get_element_id();
		    if (h4ybarID > 8) {
		      quadH4Y[n_tracks] = 8;
		    } else {
		      quadH4Y[n_tracks] = 10;
		    }
		    //std::cout << "h4y2l " << H4YBarID[n_tracks]  << " at " << n_tracks << std::endl;
		    break;
		  }
		}
	      }
	    }
	    
	    H4Y_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_H4Y2R");
	    
	    for(int det_id=44; det_id<45; ++det_id) {
	      auto iter = parID_detID_ihit.find(std::make_tuple(parID, det_id));
	      if(iter != parID_detID_ihit.end()) {
		SQHit *hit = _hit_vector->at(iter->second);
		if(hit and H4Y_hits) {
		  PHG4Hit* g4hit =  H4Y_hits->findHit(hit->get_g4hit_id());
		  if (g4hit) {
		    //if(verbosity >= 2) {
		    //  std::cout << "g4hit h4y2r " << std::endl;
		    //  g4hit->identify();
		    //}
		    int h4ybarID = hit->get_element_id();
		    
		    if(h4ybarID<9) H4YBarID[n_tracks] = 9 - hit->get_element_id();
		    if(h4ybarID>8) H4YBarID[n_tracks] = hit->get_element_id() - 8;
		    if (h4ybarID > 8) {
                      quadH4Y[n_tracks] = 9;
                    } else {
                      quadH4Y[n_tracks] = 11;
                    }
		    //std::cout << "h4y2r " << H4YBarID[n_tracks] << " at " << n_tracks << std::endl;
		    break;
		  }
		}
	      }
	    }
	    if(verbosity>=2) std::cout << "H4Y Truth info done." << std::endl;

	    //Find the hit position and barID for the first Dark Photon Detector
	    PHG4HitContainer* DP1Container;
	    for(int det_id = 55; det_id<=58; det_id++){
	      switch(det_id){
	      case 55:
		DP1Container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_DP1TL");
		break;
		
	      case 56:
		DP1Container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_DP1TR");
		break;
		
	      case 57:
		DP1Container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_DP1BL");
		break;
		
	      case 58:
		DP1Container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_DP1BR");
		break;
		
	      }
	      
	      auto iter = parID_detID_ihit.find(std::make_tuple(parID, det_id));
	      if(iter != parID_detID_ihit.end()) {
		SQHit *hit = _hit_vector->at(iter->second);
		if(hit and DP1Container) {
		  PHG4Hit* g4hit =  DP1Container->findHit(hit->get_g4hit_id());
		  if (g4hit) {
		    //if(verbosity >= 2) {
		    // LogDebug("dp1hit: " << iter->second);
		    //  g4hit->identify();
		    //}
		    GDP1_X[n_tracks] = g4hit->get_x(0);
		    GDP1_Y[n_tracks] = g4hit->get_y(0);
		    GDP1_Z[n_tracks] = g4hit->get_z(0);
		    DP1BarID[n_tracks] = hit->get_element_id();
		    if(hit->get_detector_id()==det_id){
		      switch(det_id){
		      case 55:
			quadDP1[n_tracks] = (80+hit->get_element_id())*(-1);
			break;
		      case 56:
			quadDP1[n_tracks] = (80+hit->get_element_id());
			break;
		      case 57:
			quadDP1[n_tracks] = (81-hit->get_element_id())*(-1);
			break;
		      case 58:
			quadDP1[n_tracks] = 81-hit->get_element_id();
			break;
		      }
		    }
		    break;
		  }
		}
	      }
	    }
	    if(verbosity>=2) std::cout << "DP1 Truth info done." << std::endl;
	    
	    //Find the hit position and barID for the second Dark Photon Detector
	    PHG4HitContainer* DP2Container;
	    
	    for(int det_id = 59; det_id<=62; det_id++){
	      switch(det_id){
	      case 59:
		DP2Container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_DP2TL");
		break;
		
	      case 60:
		DP2Container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_DP2TR");
		break;
		
	      case 61:
		DP2Container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_DP2BL");
		break;
		
	      case 62:
		DP2Container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_DP2BR");
		break;
		
	      }
	      
	      auto iter = parID_detID_ihit.find(std::make_tuple(parID, det_id));
	      if(iter != parID_detID_ihit.end()) {
		SQHit *hit = _hit_vector->at(iter->second);
		//if(verbosity >= Fun4AllBase::VERBOSITY_A_LOT) {
                //  LogDebug("dp2hit: " << iter->second);
		//  hit->identify();
		//}
		if(hit and DP2Container) {
		  PHG4Hit* g4hit =  DP2Container->findHit(hit->get_g4hit_id());
		  if (g4hit) {
		    GDP2_X[n_tracks] = g4hit->get_x(0);
		    GDP2_Y[n_tracks] = g4hit->get_y(0);
		    GDP2_Z[n_tracks] = g4hit->get_z(0);
		    DP2BarID[n_tracks] = hit->get_element_id();
		    
		    if(hit->get_detector_id()==det_id){
		      switch(det_id){
		      case 59:
			quadDP2[n_tracks] = (50+hit->get_element_id())*(-1);
			break;
		      case 60:
			quadDP2[n_tracks] = (50+hit->get_element_id());
			break;
		      case 61:
			quadDP2[n_tracks] = (51-hit->get_element_id())*(-1);
			break;
		      case 62:
			quadDP2[n_tracks] = 51-hit->get_element_id();
			break;
		      }
		    }
		    
		    break;
		  }
		}
	      }
	    }
	    if(verbosity>=2) std::cout << "DP2 Truth info done." << std::endl;
	    
	    gnhits[n_tracks] =
	      parID_nhits_dc[parID] +
	      parID_nhits_hodo[parID] +
	      parID_nhits_prop[parID];
	    
	    gndc[n_tracks] = parID_nhits_dc[parID];
	    gnhodo[n_tracks] = parID_nhits_hodo[parID];
	    gnprop[n_tracks] = parID_nhits_prop[parID];
	    gnDP[n_tracks] = parID_nhits_dp[parID];
	    gnH4Y[n_tracks] = parID_nhits_H4Y[parID];
	    
	    for(auto detid_elmid : parID_detid_elmid[parID]) {
	      int detid = detid_elmid.first;
	      int elmid = detid_elmid.second;
	      if(detid>62) {
		LogWarning("detid>62");
		continue;
	      }
	      gelmid[n_tracks][detid] = elmid;
	    }
	    if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT) LogInfo("particle eval finished");
            ++n_tracks;
            if(n_tracks>=1000) break;
	  }

	  //All the info in the particle eval section above is looking at the truth info from geant4. 
	  //This info is not affected by reconstruction of tracks.
	
          for(auto iter=_truth->GetPrimaryParticleRange().first; iter!=_truth->GetPrimaryParticleRange().second; ++iter) {
	    PHG4Particle * par = iter->second;

	    const double mu_mass = 0.106;
	    if (abs(par->get_pid()) == 13) {
	      auto iter2 = iter;
	      iter2++;
	      for (;iter2 != _truth->GetPrimaryParticleRange().second; ++iter2) {
		PHG4Particle* par2 = iter2->second;

		//std::cout << "getting muons " << std::endl;
		// Un-like charged
		if(par->get_pid()+par2->get_pid()!=0) continue;
		
		// same vtx
		if(par->get_vtx_id() != par2->get_vtx_id()) continue;
		
		TLorentzVector par1_mom;
		par1_mom.SetXYZM(
				 par->get_px(),
				 par->get_py(),
				 par->get_pz(),
				 mu_mass
				 );
		
		TLorentzVector par2_mom;
		par2_mom.SetXYZM(
				 par2->get_px(),
				 par2->get_py(),
				 par2->get_pz(),
				 mu_mass
				 );
		
		TLorentzVector vphoton = par1_mom + par2_mom;
		dimu_gpx[gndimu] = vphoton.Px();
		dimu_gpy[gndimu] = vphoton.Py();
		dimu_gpz[gndimu] = vphoton.Pz();
		dimu_gpt[gndimu] = vphoton.Pt();
		dimu_gmass[gndimu] = vphoton.M();
		dimu_geta[gndimu] = vphoton.Eta();
		dimu_gphi[gndimu] = vphoton.Phi();
		
		int trkID1 = par->get_track_id();//std::get<0>(parID_bestRecID[par->get_track_id()]);
		int trkID2 = par2->get_track_id();//std::get<0>(parID_bestRecID[par2->get_track_id()]);

		par1_H1[gndimu] = H1BarID[trkID1-1];
		par1_H2[gndimu] = H2BarID[trkID1-1];
		par1_H4[gndimu] = H4YBarID[trkID1-1];
		par1_DP1[gndimu] = DP1BarID[trkID1-1];
		par1_DP2[gndimu] = DP2BarID[trkID1-1];
		par1_quadDP1[gndimu] = quadDP1[trkID1-1];
		par1_quadDP2[gndimu] = quadDP2[trkID1-1];
		par1_quadH4Y[gndimu] = quadH4Y[trkID1-1];
		
		par2_H1[gndimu] = H1BarID[trkID2-1];
		par2_H2[gndimu] = H2BarID[trkID2-1];
		par2_H4[gndimu] = H4YBarID[trkID2-1];
		par2_DP1[gndimu] = DP1BarID[trkID2-1];
		par2_DP2[gndimu] = DP2BarID[trkID2-1];
		par2_quadDP1[gndimu] = quadDP1[trkID2-1];
		par2_quadDP2[gndimu] = quadDP2[trkID2-1];
		par2_quadH4Y[gndimu] = quadH4Y[trkID2-1];

		trkID1++;
                trkID2++;
				
		dimu_gnDP[gndimu] = parID_nhits_dp[trkID1]
		  +parID_nhits_dp[trkID2]
		  +parID_nhits_H4Y[trkID1]
		  +parID_nhits_H4Y[trkID2];
		
		par1_DP[gndimu] = parID_nhits_dp[trkID1];
		par1_H4Y[gndimu] = parID_nhits_H4Y[trkID1];
		par1_hodo[gndimu] = parID_nhits_hodo[trkID1];
		
		par2_DP[gndimu] = parID_nhits_dp[trkID2];
		par2_hodo[gndimu] = parID_nhits_hodo[trkID2];
		par2_H4Y[gndimu] = parID_nhits_H4Y[trkID2];
		  
		++gndimu;
		if(gndimu>=1000) break;
	      }
	    }
	    if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT) LogInfo("dimu eval finished");
	  }
	}
	
	_tout_truth->Fill();
	
	if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT)
	  std::cout << "Leaving TrkEval::TruthEval: " << _event << std::endl;
	
	return Fun4AllReturnCodes::EVENT_OK;
}

int TrkEval::End(PHCompositeNode* topNode) {
  if(Verbosity() >= Fun4AllBase::VERBOSITY_A_LOT)
    std::cout << "TrkEval::End" << std::endl;
  
  PHTFileServer::get().cd(_out_name.c_str());
  _tout_truth->Write();
  //_tout_reco->Write();
  return Fun4AllReturnCodes::EVENT_OK;
}

int TrkEval::InitEvalTree() {
  PHTFileServer::get().open(_out_name.c_str(), "RECREATE");

  _tout_truth = new TTree("Truth", "Truth Eval");
  _tout_truth->Branch("runID",         &run_id,          "runID/I");
  _tout_truth->Branch("spillID",       &spill_id,        "spillID/I");
  _tout_truth->Branch("liveProton",    &target_pos,      "liveProton/F");
  _tout_truth->Branch("eventID",       &event_id,        "eventID/I");
  _tout_truth->Branch("emu_trigger",   &emu_trigger,      "emu_trigger/s");
  _tout_truth->Branch("krecstat",      &krecstat,        "krecstat/I");

  _tout_truth->Branch("nHits",         &n_hits,          "nHits/I");
  _tout_truth->Branch("hitID",         hit_id,           "hitID[nHits]/I");
  _tout_truth->Branch("detectorID",    detector_id,      "detectorID[nHits]/I");
  _tout_truth->Branch("elementID",     element_id,       "elementID[nHits]/I");
  _tout_truth->Branch("hodo_mask",     hodo_mask,        "hodo_mask[nHits]/I");
  _tout_truth->Branch("detectorZ",     detector_z,       "detectorZ[nHits]/F");
  _tout_truth->Branch("truth_x",       truth_x,          "truth_x[nHits]/F");
  _tout_truth->Branch("truth_y",       truth_y,          "truth_y[nHits]/F");
  _tout_truth->Branch("truth_z",       truth_z,          "truth_z[nHits]/F");
  _tout_truth->Branch("truth_pos",     truth_pos,        "truth_pos[nHits]/F");
  _tout_truth->Branch("driftDistance", drift_distance,   "driftDistance[nHits]/F");
  _tout_truth->Branch("pos",           pos,              "pos[nHits]/F");

  _tout_truth->Branch("n_tracks",      &n_tracks,           "n_tracks/I");
  _tout_truth->Branch("par_id",        par_id,              "par_id[n_tracks]/I");
  _tout_truth->Branch("rec_id",        rec_id,              "rec_id[n_tracks]/I");
  _tout_truth->Branch("pid",           pid,                 "pid[n_tracks]/I");
  _tout_truth->Branch("gvx",           gvx,                 "gvx[n_tracks]/F");
  _tout_truth->Branch("gvy",           gvy,                 "gvy[n_tracks]/F");
  _tout_truth->Branch("gvz",           gvz,                 "gvz[n_tracks]/F");
  _tout_truth->Branch("gpx",           gpx,                 "gpx[n_tracks]/F");
  _tout_truth->Branch("gpy",           gpy,                 "gpy[n_tracks]/F");
  _tout_truth->Branch("gpz",           gpz,                 "gpz[n_tracks]/F");
  _tout_truth->Branch("gpt",           gpt,                 "gpt[n_tracks]/F");
  _tout_truth->Branch("geta",          geta,                "geta[n_tracks]/F");
  _tout_truth->Branch("gphi",          gphi,                "gphi[n_tracks]/F");
  _tout_truth->Branch("gx_st1",        gx_st1,              "gx_st1[n_tracks]/F");
  _tout_truth->Branch("gy_st1",        gy_st1,              "gy_st1[n_tracks]/F");
  _tout_truth->Branch("gz_st1",        gz_st1,              "gz_st1[n_tracks]/F");

  _tout_truth->Branch("gx_st2",        gx_st2,              "gx_st2[n_tracks]/F");
  _tout_truth->Branch("gy_st2",        gy_st2,              "gy_st2[n_tracks]/F");
  _tout_truth->Branch("gz_st2",        gz_st2,              "gz_st2[n_tracks]/F");


  _tout_truth->Branch("gx_st3",        gx_st3,              "gx_st3[n_tracks]/F");
  _tout_truth->Branch("gy_st3",        gy_st3,              "gy_st3[n_tracks]/F");
  _tout_truth->Branch("gz_st3",        gz_st3,              "gz_st3[n_tracks]/F");

  _tout_truth->Branch("gpx_st1",       gpx_st1,             "gpx_st1[n_tracks]/F");
  _tout_truth->Branch("gpy_st1",       gpy_st1,             "gpy_st1[n_tracks]/F");
  _tout_truth->Branch("gpz_st1",       gpz_st1,             "gpz_st1[n_tracks]/F");

  _tout_truth->Branch("gpx_st2",       gpx_st2,             "gpx_st2[n_tracks]/F");
  _tout_truth->Branch("gpy_st2",       gpy_st2,             "gpy_st2[n_tracks]/F");
  _tout_truth->Branch("gpz_st2",       gpz_st2,             "gpz_st2[n_tracks]/F");

  _tout_truth->Branch("gpx_st3",       gpx_st3,             "gpx_st3[n_tracks]/F");
  _tout_truth->Branch("gpy_st3",       gpy_st3,             "gpy_st3[n_tracks]/F");
  _tout_truth->Branch("gpz_st3",       gpz_st3,             "gpz_st3[n_tracks]/F");

  _tout_truth->Branch("gnhits",        gnhits,              "gnhits[n_tracks]/I");
  _tout_truth->Branch("gndc",          gndc,                "gndc[n_tracks]/I");
  _tout_truth->Branch("gnhodo",        gnhodo,              "gnhodo[n_tracks]/I");
  _tout_truth->Branch("gnprop",        gnprop,              "gnprop[n_tracks]/I");
  _tout_truth->Branch("gnH4Y",         gnH4Y,               "gnH4Y[n_tracks]/F");
  _tout_truth->Branch("gelmid",        gelmid,              "gelmid[n_tracks][54]/I");

  _tout_truth->Branch("ntruhits",      ntruhits,            "ntruhits[n_tracks]/I");
  _tout_truth->Branch("nhits",         nhits,               "nhits[n_tracks]/I");
  _tout_truth->Branch("charge",        charge,              "charge[n_tracks]/I");
  _tout_truth->Branch("vx",            vx,                  "vx[n_tracks]/F");
  _tout_truth->Branch("vy",            vy,                  "vy[n_tracks]/F");
  _tout_truth->Branch("vz",            vz,                  "vz[n_tracks]/F");
  _tout_truth->Branch("px",            px,                  "px[n_tracks]/F");
  _tout_truth->Branch("py",            py,                  "py[n_tracks]/F");
  _tout_truth->Branch("pz",            pz,                  "pz[n_tracks]/F");
  _tout_truth->Branch("pt",            pt,                  "pt[n_tracks]/F");
  _tout_truth->Branch("eta",           eta,                 "eta[n_tracks]/F");
  _tout_truth->Branch("phi",           phi,                 "phi[n_tracks]/F");

  _tout_truth->Branch("x_st1",         x_st1,               "x_st1[n_tracks]/F");
  _tout_truth->Branch("y_st1",         y_st1,               "y_st1[n_tracks]/F");
  _tout_truth->Branch("x_st3",         x_st3,               "x_st3[n_tracks]/F");
  _tout_truth->Branch("y_st3",         y_st3,               "y_st3[n_tracks]/F");
  _tout_truth->Branch("z_st3",         z_st3,               "z_st3[n_tracks]/F");

  _tout_truth->Branch("px_st1",        px_st1,              "px_st1[n_tracks]/F");
  _tout_truth->Branch("py_st1",        py_st1,              "py_st1[n_tracks]/F");
  _tout_truth->Branch("pz_st1",        pz_st1,              "pz_st1[n_tracks]/F");
  _tout_truth->Branch("px_st3",        px_st3,              "px_st3[n_tracks]/F");
  _tout_truth->Branch("py_st3",        py_st3,              "py_st3[n_tracks]/F");
  _tout_truth->Branch("pz_st3",        pz_st3,              "pz_st3[n_tracks]/F");

  _tout_truth->Branch("gndimu",        &gndimu,              "gndimu/I");
  _tout_truth->Branch("dimu_gpx",      dimu_gpx,             "dimu_gpx[gndimu]/F");
  _tout_truth->Branch("dimu_gpy",      dimu_gpy,             "dimu_gpy[gndimu]/F");
  _tout_truth->Branch("dimu_gpz",      dimu_gpz,             "dimu_gpz[gndimu]/F");
  _tout_truth->Branch("dimu_gpt",      dimu_gpt,             "dimu_gpt[gndimu]/F");
  _tout_truth->Branch("dimu_gmass",    dimu_gmass,           "dimu_gmass[gndimu]/F");
  _tout_truth->Branch("dimu_geta",     dimu_geta,            "dimu_geta[gndimu]/F");
  _tout_truth->Branch("dimu_gphi",     dimu_gphi,            "dimu_gphi[gndimu]/F");

  _tout_truth->Branch("dimu_nrec",     dimu_nrec,            "dimu_nrec[gndimu]/I");
  _tout_truth->Branch("dimu_px",       dimu_px,              "dimu_px[gndimu]/F");
  _tout_truth->Branch("dimu_py",       dimu_py,              "dimu_py[gndimu]/F");
  _tout_truth->Branch("dimu_pz",       dimu_pz,              "dimu_pz[gndimu]/F");
  _tout_truth->Branch("dimu_pt",       dimu_pt,              "dimu_pt[gndimu]/F");
  _tout_truth->Branch("dimu_mass",     dimu_mass,            "dimu_mass[gndimu]/F");
  _tout_truth->Branch("dimu_eta",      dimu_eta,             "dimu_eta[gndimu]/F");
  _tout_truth->Branch("dimu_phi",      dimu_phi,             "dimu_phi[gndimu]/F");
  _tout_truth->Branch("dimu_gnDP",     dimu_gnDP,            "dimu_gnDP[gndimu]/F");

  _tout_truth->Branch("par1_H1",       par1_H1,             "par1_H1[gndimu]/F");
  _tout_truth->Branch("par1_H2",       par1_H2,             "par1_H2[gndimu]/F");
  _tout_truth->Branch("par1_H4",       par1_H4,             "par1_H4[gndimu]/F");
  _tout_truth->Branch("par1_DP1",      par1_DP1,            "par1_DP1[gndimu]/F");
  _tout_truth->Branch("par1_DP2",      par1_DP2,            "par1_DP2[gndimu]/F");
  _tout_truth->Branch("par1_quadDP1",  par1_quadDP1,        "par1_quadDP1[gndimu]/F");
  _tout_truth->Branch("par1_quadDP2",  par1_quadDP2,        "par1_quadDP2[gndimu]/F");
  _tout_truth->Branch("par1_quadH4Y",  par1_quadH4Y,        "par1_quadH4Y[gndimu]/F");

  _tout_truth->Branch("par2_H1",       par2_H1,             "par2_H1[gndimu]/F");
  _tout_truth->Branch("par2_H2",       par2_H2,             "par2_H2[gndimu]/F");
  _tout_truth->Branch("par2_H4",       par2_H4,             "par2_H4[gndimu]/F");
  _tout_truth->Branch("par2_DP1",      par2_DP1,            "par2_DP1[gndimu]/F");
  _tout_truth->Branch("par2_DP2",      par2_DP2,            "par2_DP2[gndimu]/F");
  _tout_truth->Branch("par2_quadDP1",  par2_quadDP1,        "par2_quadDP1[gndimu]/F");
  _tout_truth->Branch("par2_quadDP2",  par2_quadDP2,        "par2_quadDP2[gndimu]/F");
  _tout_truth->Branch("par2_quadH4Y",  par2_quadH4Y,        "par2_quadH4Y[gndimu]/F");

  _tout_truth->Branch("GDP1_X",        GDP1_X,               "GDP1_X[n_tracks]/F");
  _tout_truth->Branch("GDP1_Y",        GDP1_Y,               "GDP1_Y[n_tracks]/F");
  _tout_truth->Branch("GDP1_Z",        GDP1_Z,               "GDP1_Z[n_tracks]/F");
  _tout_truth->Branch("DP1BarID",      DP1BarID,             "DP1BarID[n_tracks]/F");

  _tout_truth->Branch("GDP2_X",        GDP2_X,               "GDP2_X[n_tracks]/F");
  _tout_truth->Branch("GDP2_Y",        GDP2_Y,               "GDP2_Y[n_tracks]/F");
  _tout_truth->Branch("GDP2_Z",        GDP2_Z,               "GDP2_Z[n_tracks]/F");
  _tout_truth->Branch("DP2BarID",      DP2BarID,             "DP2BarID[n_tracks]/F");

  _tout_truth->Branch("H4YBarID",      H4YBarID,             "H4YBarID[n_tracks]/F");
  _tout_truth->Branch("gnDP",          gnDP,                 "gnDP[n_tracks]/F");
  _tout_truth->Branch("par1_DP",       par1_DP,              "par1_DP[n_tracks]/F");
  _tout_truth->Branch("par2_DP",       par2_DP,              "par2_DP[n_tracks]/F");
  _tout_truth->Branch("par1_hodo",     par1_hodo,            "par1_hodo[n_tracks]/F");
  _tout_truth->Branch("par2_hodo",     par2_hodo,            "par2_hodo[n_tracks]/F");
  _tout_truth->Branch("par1_H4Y",      par1_H4Y,             "par1_H4Y[n_tracks]/F");
  _tout_truth->Branch("par2_H4Y",      par2_H4Y,             "par2_H4Y[n_tracks]/F");
  _tout_truth->Branch("H1BarID",       H1BarID,              "H1BarID[n_tracks]/F");
  _tout_truth->Branch("H2BarID",       H2BarID,              "H2BarID[n_tracks]/F");

  _tout_truth->Branch("quadDP1",       quadDP1,              "quadDP1[n_tracks]/F");
  _tout_truth->Branch("quadDP2",       quadDP2,              "quadDP2[n_tracks]/F");
  _tout_truth->Branch("quadH4Y",       quadH4Y,              "quadH4Y[n_tracks]/F");

  /*
  _tout_reco = new TTree("Reco", "Reco Eval");
  _tout_reco->Branch("krecstat",      &krecstat,           "krecstat/I");
  _tout_reco->Branch("n_tracks",      &n_tracks,           "n_tracks/I");
  _tout_reco->Branch("par_id",        par_id,              "par_id[n_tracks]/I");
  _tout_reco->Branch("rec_id",        rec_id,              "rec_id[n_tracks]/I");
  _tout_reco->Branch("pid",           pid,                 "pid[n_tracks]/I");
  _tout_reco->Branch("gvx",           gvx,                 "gvx[n_tracks]/F");
  _tout_reco->Branch("gvy",           gvy,                 "gvy[n_tracks]/F");
  _tout_reco->Branch("gvz",           gvz,                 "gvz[n_tracks]/F");
  _tout_reco->Branch("gpx",           gpx,                 "gpx[n_tracks]/F");
  _tout_reco->Branch("gpy",           gpy,                 "gpy[n_tracks]/F");
  _tout_reco->Branch("gpz",           gpz,                 "gpz[n_tracks]/F");
  _tout_reco->Branch("gpt",           gpt,                 "gpt[n_tracks]/F");
  _tout_reco->Branch("geta",          geta,                "geta[n_tracks]/F");
  _tout_reco->Branch("gphi",          gphi,                "gphi[n_tracks]/F");
  _tout_reco->Branch("gx_st1",        gx_st1,              "gx_st1[n_tracks]/F");
  _tout_reco->Branch("gy_st1",        gy_st1,              "gy_st1[n_tracks]/F");
  _tout_reco->Branch("gz_st1",        gz_st1,              "gz_st1[n_tracks]/F");
  _tout_reco->Branch("gx_st3",        gx_st3,              "gx_st3[n_tracks]/F");
  _tout_reco->Branch("gy_st3",        gy_st3,              "gy_st3[n_tracks]/F");
  _tout_reco->Branch("gz_st3",        gz_st3,              "gz_st3[n_tracks]/F");

  _tout_reco->Branch("gpx_st1",       gpx_st1,             "gpx_st1[n_tracks]/F");
  _tout_reco->Branch("gpy_st1",       gpy_st1,             "gpy_st1[n_tracks]/F");
  _tout_reco->Branch("gpz_st1",       gpz_st1,             "gpz_st1[n_tracks]/F");

  _tout_reco->Branch("gpx_st3",       gpx_st3,             "gpx_st3[n_tracks]/F");
  _tout_reco->Branch("gpy_st3",       gpy_st3,             "gpy_st3[n_tracks]/F");
  _tout_reco->Branch("gpz_st3",       gpz_st3,             "gpz_st3[n_tracks]/F");
 
  _tout_reco->Branch("gnhits",        gnhits,              "gnhits[n_tracks]/I");
  _tout_reco->Branch("gndc",          gndc,                "gndc[n_tracks]/I");
  _tout_reco->Branch("gnhodo",        gnhodo,              "gnhodo[n_tracks]/I");
  _tout_reco->Branch("gnprop",        gnprop,              "gnprop[n_tracks]/I");
 
  //_tout_reco->Branch("gelmid",        gelmid,              "gelmid[n_tracks][54]/I");

  _tout_reco->Branch("ntruhits",      ntruhits,            "ntruhits[n_tracks]/I");
  _tout_reco->Branch("nhits",         nhits,               "nhits[n_tracks]/I");
  _tout_reco->Branch("charge",        charge,              "charge[n_tracks]/I");
  _tout_reco->Branch("vx",            vx,                  "vx[n_tracks]/F");
  _tout_reco->Branch("vy",            vy,                  "vy[n_tracks]/F");
  _tout_reco->Branch("vz",            vz,                  "vz[n_tracks]/F");
  _tout_reco->Branch("px",            px,                  "px[n_tracks]/F");
  _tout_reco->Branch("py",            py,                  "py[n_tracks]/F");
  _tout_reco->Branch("pz",            pz,                  "pz[n_tracks]/F");
  _tout_reco->Branch("pt",            pt,                  "pt[n_tracks]/F");
  _tout_reco->Branch("eta",           eta,                 "eta[n_tracks]/F");
  _tout_reco->Branch("phi",           phi,                 "phi[n_tracks]/F");

  _tout_reco->Branch("x_st1",         x_st1,               "x_st1[n_tracks]/F");
  _tout_reco->Branch("y_st1",         y_st1,               "y_st1[n_tracks]/F");


  _tout_reco->Branch("x_st3",         x_st3,               "x_st3[n_tracks]/F");
  _tout_reco->Branch("y_st3",         y_st3,               "y_st3[n_tracks]/F");
  _tout_reco->Branch("z_st3",         z_st3,               "z_st3[n_tracks]/F");

  _tout_reco->Branch("px_st1",        px_st1,              "px_st1[n_tracks]/F");
  _tout_reco->Branch("py_st1",        py_st1,              "py_st1[n_tracks]/F");
  _tout_reco->Branch("pz_st1",        pz_st1,              "pz_st1[n_tracks]/F");

  _tout_reco->Branch("px_st3",        px_st3,              "px_st3[n_tracks]/F");
  _tout_reco->Branch("py_st3",        py_st3,              "py_st3[n_tracks]/F");
  _tout_reco->Branch("pz_st3",        pz_st3,              "pz_st3[n_tracks]/F");
  */
  return 0;
}

int TrkEval::ResetEvalVars() {
  run_id = std::numeric_limits<int>::max();
  spill_id = std::numeric_limits<int>::max();
  target_pos = std::numeric_limits<float>::max();
  event_id = std::numeric_limits<int>::max();
  emu_trigger = 0;
  krecstat = std::numeric_limits<int>::max();
 
  n_hits = 0;
  for(int i=0; i<10000; ++i) {
    detector_id[i]    = std::numeric_limits<short>::max();
    element_id[i]     = std::numeric_limits<short>::max();
    hodo_mask[i]      = std::numeric_limits<short>::max();
    drift_distance[i] = std::numeric_limits<float>::max();
    pos[i]            = std::numeric_limits<float>::max();
    detector_z[i]     = std::numeric_limits<float>::max();

    truth_x[i]       = std::numeric_limits<float>::max();
    truth_y[i]       = std::numeric_limits<float>::max();
    truth_z[i]       = std::numeric_limits<float>::max();
    truth_pos[i]     = std::numeric_limits<float>::max();
  }

  n_tracks = 0;
  for(int i=0; i<1000; ++i) {
    rec_id[i]     = std::numeric_limits<int>::max();
    par_id[i]     = std::numeric_limits<int>::max();
    pid[i]        = std::numeric_limits<int>::max();
    gvx[i]        = std::numeric_limits<float>::max();
    gvy[i]        = std::numeric_limits<float>::max();
    gvz[i]        = std::numeric_limits<float>::max();
    gpx[i]        = std::numeric_limits<float>::max();
    gpy[i]        = std::numeric_limits<float>::max();
    gpz[i]        = std::numeric_limits<float>::max();
    gpt[i]        = std::numeric_limits<float>::max();
    geta[i]       = std::numeric_limits<float>::max();
    gphi[i]       = std::numeric_limits<float>::max();
    gnhits[i]     = std::numeric_limits<int>::max();
    gx_st1[i]     = std::numeric_limits<float>::max();
    gy_st1[i]     = std::numeric_limits<float>::max();
    gz_st1[i]     = std::numeric_limits<float>::max();
    gpx_st1[i]    = std::numeric_limits<float>::max();
    gpy_st1[i]    = std::numeric_limits<float>::max();
    gpz_st1[i]    = std::numeric_limits<float>::max();

    gx_st3[i]     = std::numeric_limits<float>::max();
    gy_st3[i]     = std::numeric_limits<float>::max();
    gz_st3[i]     = std::numeric_limits<float>::max();
    gpx_st3[i]    = std::numeric_limits<float>::max();
    gpy_st3[i]    = std::numeric_limits<float>::max();
    gpz_st3[i]    = std::numeric_limits<float>::max();

    GDP1_X[i] = std::numeric_limits<float>::max();
    GDP1_Y[i] = std::numeric_limits<float>::max();
    GDP1_Z[i] = std::numeric_limits<float>::max();

    GDP2_X[i] = std::numeric_limits<float>::max();
    GDP2_Y[i] = std::numeric_limits<float>::max();
    GDP2_Z[i] = std::numeric_limits<float>::max();

    DP1BarID[i] = std::numeric_limits<float>::max();
    DP2BarID[i] = std::numeric_limits<float>::max();
    H4YBarID[i] = std::numeric_limits<float>::max();
    H1BarID[i] = std::numeric_limits<float>::max();
    H2BarID[i] = std::numeric_limits<float>::max();

    quadDP1[i] = std::numeric_limits<float>::max();
    quadDP2[i] = std::numeric_limits<float>::max();
    quadH4Y[i] = std::numeric_limits<float>::max();

    gnDP[i] = std::numeric_limits<float>::max();
    gnH4Y[i] = std::numeric_limits<float>::max();

    gndc[i]       = std::numeric_limits<int>::max();
    gnhodo[i]     = std::numeric_limits<int>::max();
    gnprop[i]     = std::numeric_limits<int>::max();

    for(int j=0; j<55; ++j) {
    	gelmid[i][j] = std::numeric_limits<int>::max();
    }

    ntruhits[i]   = std::numeric_limits<int>::max();
    nhits[i]      = std::numeric_limits<int>::max();
    charge[i]     = std::numeric_limits<int>::max();
    vx[i]         = std::numeric_limits<float>::max();
    vy[i]         = std::numeric_limits<float>::max();
    vz[i]         = std::numeric_limits<float>::max();
    px[i]         = std::numeric_limits<float>::max();
    py[i]         = std::numeric_limits<float>::max();
    pz[i]         = std::numeric_limits<float>::max();
    pt[i]         = std::numeric_limits<float>::max();
    eta[i]        = std::numeric_limits<float>::max();
    phi[i]        = std::numeric_limits<float>::max();
    x_st1[i]     = std::numeric_limits<float>::max();
    y_st1[i]     = std::numeric_limits<float>::max();
    px_st1[i]     = std::numeric_limits<float>::max();
    py_st1[i]     = std::numeric_limits<float>::max();
    pz_st1[i]     = std::numeric_limits<float>::max();

    x_st3[i]     = std::numeric_limits<float>::max();
    y_st3[i]     = std::numeric_limits<float>::max();
    z_st3[i]     = std::numeric_limits<float>::max();
    px_st3[i]     = std::numeric_limits<float>::max();
    py_st3[i]     = std::numeric_limits<float>::max();
    pz_st3[i]     = std::numeric_limits<float>::max();
  }

  gndimu = 0;
  for(int i=0; i<1000; ++i) {
    dimu_gpx[i]        = std::numeric_limits<float>::max();
    dimu_gpy[i]        = std::numeric_limits<float>::max();
    dimu_gpz[i]        = std::numeric_limits<float>::max();
    dimu_gpt[i]        = std::numeric_limits<float>::max();
    dimu_gmass[i]      = std::numeric_limits<float>::max();
    dimu_geta[i]       = std::numeric_limits<float>::max();
    dimu_gphi[i]       = std::numeric_limits<float>::max();
    
    dimu_nrec[i]       = 0;
    dimu_px[i]         = std::numeric_limits<float>::max();
    dimu_py[i]         = std::numeric_limits<float>::max();
    dimu_pz[i]         = std::numeric_limits<float>::max();
    dimu_pt[i]         = std::numeric_limits<float>::max();
    dimu_mass[i]       = std::numeric_limits<float>::max();
    dimu_eta[i]        = std::numeric_limits<float>::max();
    dimu_phi[i]        = std::numeric_limits<float>::max();
    dimu_gnDP[i]       = std::numeric_limits<float>::max();
    par1_DP[i]         = std::numeric_limits<float>::max();
    par1_hodo[i]       = std::numeric_limits<float>::max();
    
    par2_DP[i]         = std::numeric_limits<float>::max(); 
    par2_hodo[i]       = std::numeric_limits<float>::max();
    
    par1_H1[i]         = std::numeric_limits<float>::max();
    par1_H2[i]         = std::numeric_limits<float>::max();
    par1_H4[i]         = std::numeric_limits<float>::max();
    par1_DP1[i]        = std::numeric_limits<float>::max();
    par1_DP2[i]        = std::numeric_limits<float>::max();
    par1_quadDP1[i]    = std::numeric_limits<float>::max();
    par1_quadDP2[i]    = std::numeric_limits<float>::max();
    par1_quadH4Y[i]    = std::numeric_limits<float>::max();
    
    par2_H1[i]         = std::numeric_limits<float>::max();
    par2_H2[i]         = std::numeric_limits<float>::max();
    par2_H4[i]         = std::numeric_limits<float>::max();
    par2_DP1[i]        = std::numeric_limits<float>::max();
    par2_DP2[i]        = std::numeric_limits<float>::max();
    par2_quadDP1[i]    = std::numeric_limits<float>::max();
    par2_quadDP2[i]    = std::numeric_limits<float>::max();
    par2_quadH4Y[i]    = std::numeric_limits<float>::max();
  }

  return 0;
}

int TrkEval::GetNodes(PHCompositeNode* topNode) {

  _run_header = findNode::getClass<SQRun>(topNode, "SQRun");
  if (!_run_header) {
    LogError("!_run_header");
    //return Fun4AllReturnCodes::ABORTEVENT;
  }

  _spill_map = findNode::getClass<SQSpillMap>(topNode, "SQSpillMap");
  if (!_spill_map) {
    LogError("!_spill_map");
    //return Fun4AllReturnCodes::ABORTEVENT;
  }

  _event_header = findNode::getClass<SQEvent>(topNode, "SQEvent");
  if (!_event_header) {
    LogError("!_event_header");
    //return Fun4AllReturnCodes::ABORTEVENT;
  }

  if(_hit_container_type.find("Map") != std::string::npos) {
    _hit_map = findNode::getClass<SQHitMap>(topNode, "SQHitMap");
    if (!_hit_map) {
      LogError("!_hit_map");
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if(_hit_container_type.find("Vector") != std::string::npos) {
    _hit_vector = findNode::getClass<SQHitVector>(topNode, "SQHitVector");
    if (!_hit_vector) {
      LogError("!_hit_vector");
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  _truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_truth) {
    LogError("!_truth");
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _recEvent = findNode::getClass<SRecEvent>(topNode, "SRecEvent");
  if (!_recEvent) {
    LogError("!_recEvent");
    //return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}







