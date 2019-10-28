#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>
#include <set>
#include <unordered_set>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TClonesArray.h>

#include "GeomSvc.h"
#include "SRawEvent.h"
#include "SRecEvent.h"
#include "FastTracklet.h"

using namespace std;

// constants
float wDP1 = 1.03; float wDP2 = 1.90; float wH4Y = 23;
float xDP1 = 798; float xDP2 = 1471; float xH4Y = 2200;

// hit triplet
struct Triplet
{ 
  int x, y, z; 
}; 

double getDelay(int quadID, int barID) {
    switch(quadID) {
        case 0:
            return 1728.0;
        case 1:
            return 1720.0;
        case 2:
            return 1730.0;
        case 3:
            return (barID>15)?1716.0:1721.0;
        case 4:
            return (barID>32)?1695.0:1703.0;
        case 5:
            return 1697.0;
        case 6:
            return 1704.0;
        case 7:
            return (barID>32)?1695.0:1689.0;
        case 8:
            return 1092;
        case 9:
            return 1091;
        case 10:
            return 1093;
        case 11:
            return 1092;
        default:
            return 1700;
    }
}
// given position
int findh4ybar(float y1){
  int barid = int(y1/wH4Y-.5);
  return barid;
}
int roadHash(int st1, int st2, int h4) {
    return st1*1000 + st2*10 + h4;
}
int main(int argc, char *argv[])
{
    // roads
    unordered_set<int> roadset;
    //ifstream roadsfile("code_400_real_roads.txt"); // sho's
    //ifstream roadsfile("allRoads_n1993.txt"); // dylan's
    ifstream roadsfile("allRoads_n260.txt"); // dylan's  new version  brem mass 1.05 eps -6
    //ifstream roadsfile("roads_3hits_mu12.txt"); //mine with 3 hits
    //ifstream roadsfile("roads_6hits_mu12.txt");
    //ifstream roadsfile("roads.txt")//mine
    if (!roadsfile.is_open()) return 1;
    char line[100];
    while (roadsfile.getline(line,100)) {
        //printf("line=%s\n",line);
        char * tok;
        int st1, st2, h4;
        st1 = atoi(strtok(line," "));
        st2 = atoi(strtok(NULL," "));
        h4 = atoi(strtok(NULL," "));
        int roadhash = roadHash(st1+1, st2+1, h4+1); //VHDL uses zero-indexed vectors, so we need to add 1 to match the elementIDs
        //printf("road=%d\n",roadhash);
        roadset.insert(roadhash);
    }

    //Initialization of geometry and tracked data
    SRawEvent* rawEvent = new SRawEvent();

    // Reading data
    TFile* dataFile = new TFile(argv[1], "READ");
    TTree* dataTree = (TTree*)dataFile->Get("save");
    dataTree->SetBranchAddress("rawEvent", &rawEvent);

    // hitsets is just an array of 12 sets (for each quadrant) with the barIDs where you have hits -> per event
    set<int> hitsets[12];
    set<int> hitsetsNoTimecut[12];
    set<int> hitsetsOffByOne[12];

    vector<float> hitsy1[12];
    vector<float> hits[12];
    vector<float> hitsbars[12];

    // Output
    TFile* saveFile = new TFile(argv[2], "recreate");
    TTree* saveTree = new TTree("save", "save");

    const int NROADS = 10;
    int dp_mask(-1),nim3_mask(-1),nRoads(-1);
    int quadrant[NROADS];
    int baridDP1[NROADS],baridDP2[NROADS],baridH4Y[NROADS],baridH4Y_t[NROADS];
    float yDP1[NROADS],yDP2[NROADS],yH4Y[NROADS],yh[NROADS];
    saveTree->Branch("nRoads", &nRoads, "nRoads/I");
    saveTree->Branch("dp_mask", &dp_mask, "dp_mask/I");
    saveTree->Branch("nim3_mask", &nim3_mask, "nim3_mask/I");
    saveTree->Branch("quadrant", &quadrant, TString::Format("quadrant[%i]/I", NROADS));
    saveTree->Branch("baridDP1", &baridDP1, TString::Format("baridDP1[%i]/I", NROADS));
    saveTree->Branch("baridDP2", &baridDP2, TString::Format("baridDP2[%i]/I", NROADS));
    saveTree->Branch("baridH4Y", &baridH4Y, TString::Format("baridH4Y[%i]/I", NROADS));
    saveTree->Branch("yDP1", &yDP1, TString::Format("yDP1[%i]/F", NROADS));
    saveTree->Branch("yDP2", &yDP2, TString::Format("yDP2[%i]/F", NROADS));
    saveTree->Branch("yH4Y", &yH4Y, TString::Format("yH4Y[%i]/F", NROADS));
    saveTree->Branch("yh", &yh, TString::Format("yh[%i]/F", NROADS));
    saveFile->cd();

    // actual road check
    TH2D* roadsVsTrig = new TH2D("roadsVsTrig","roadsVsTrig",200,-0.5,199.5,16,-0.5,15.5);
    TH2D* roadsVsTrigNoTimecut = new TH2D("roadsVsTrigNoTimecut","roadsVsTrigNoTimecut",200,-0.5,199.5,16,-0.5,15.5);
    TH2D* nimlikeVsTrig = new TH2D("nimlikeVsTrig","nimlikeVsTrig",200,-0.5,199.5,16,-0.5,15.5);

    // hists that save dp1 and dp2 data hits with a hit in h4y quad X
    // all in the same quad
    TH2D* DP1DP2Hists[8][4];
    TH2D* DP1DP2Hists_all[8];
    char name[500];
    for (int bar=0; bar<8; bar++) {
      sprintf(name, "DP1DP2Hist_h4ybar%i_allquad", bar+1);
      DP1DP2Hists_all[bar] = new TH2D(name,name,80,0,80,50,0,50);
      for (int quad=0; quad<4; quad++) {
	sprintf(name, "DP1DP2Hist_h4ybar%i_quad%i", bar+1, quad);
	DP1DP2Hists[bar][quad] = new TH2D(name,name,80,0,80,50,0,50);
      }
    }

    // timing info
    double rfclock = 18.8;
    double timeCut = rfclock/2;
    int dpTriggerMask = 64;
    int nim1TriggerMask = 32;
    int nim3TriggerMask = 128;

    // to compute rates
    int nBitsTrigDP(0),nBitsTrigNIM3(0);
    int nTrigDP(0),nTrigNIM3(0);
    int nMatchedDP(0),nMatchedNIM3(0),nMatchedDPlineH4Y(0),nMatchedNIM3lineH4Y(0); // just matched roads
    int nSameQuadDP(0),nSameQuadNIM3(0),nSameQuadDPlineH4Y(0),nSameQuadNIM3lineH4Y(0); // at least 2 roads in same quad
    int nSame2QuadDP(0),nSame2QuadNIM3(0),nSame2QuadDPlineH4Y(0),nSame2QuadNIM3lineH4Y(0); // exactly 2 roads in same quad
    int nDiffQuadDP(0),nDiffQuadNIM3(0),nDiffQuadDPlineH4Y(0),nDiffQuadNIM3lineH4Y(0); // at least 2 roads with diff quad
    int nTBQuadDP(0),nTBQuadNIM3(0),nTBQuadDPlineH4Y(0),nTBQuadNIM3lineH4Y(0); // at least 2 roads with diff quad (top and bottom)
    int nDiff2hitsDP(0),nDiff2hitsNIM3(0),nDiff2hitsDPlineH4Y(0),nDiff2hitsNIM3lineH4Y(0);
    int nDiff2quadDP(0),nDiff2quadNIM3(0),nDiff2quadDPlineH4Y(0),nDiff2quadNIM3lineH4Y(0);

    int nBitsTrigDPCheck(0),nBitsTrigNIM3Check(0);
    int nBitsTrigDPCheck2(0),nBitsTrigNIM3Check2(0);
    int nBitsTrigDPCheck3(0),nBitsTrigNIM3Check3(0);
    int nBitsTrigDPCheck4(0),nBitsTrigNIM3Check4(0);
    int nBitsTrigDPCheck5(0),nBitsTrigNIM3Check5(0);
    int nBitsTrigDPCheck6(0),nBitsTrigNIM3Check6(0);

    // loop over events
    vector<Triplet> vquad_data_samequad,vquad_data_h4y;
    for(Int_t i = 0; i < dataTree->GetEntries(); ++i)
    //for(Int_t i = 0; i < 1000; ++i)
    {
        dataTree->GetEntry(i);
        //if(i % 1000 == 0) cout << i << endl;

        if (rawEvent->getTriggerBits()>0) {
	  dp_mask = 0; nim3_mask = 0;
	  if (rawEvent->getTriggerBits() == dpTriggerMask) {
	    nTrigDP += 1;
	    dp_mask = 1;
	  }
	  if (rawEvent->getTriggerBits() == nim3TriggerMask) {
	    nTrigNIM3 += 1;
	    nim3_mask = 1;
	  }
	  // loop over event hits
	  for(Int_t k = 0; k < rawEvent->getNHitsAll(); ++k) {
	    Hit h = rawEvent->getHit(k);

	    // get hits detector ID and quadrant
	    int dpQuadID = -1;
	    int barID = -1;
	    float y1 = -1;
	    //if (h.detectorID >= 41 && h.detectorID < 45) dpQuadID = h.detectorID-33; //H4Y1L = 8, H4Y1R = 9, etc.
	    if (h.detectorID >= 55 && h.detectorID <= 62) {
	      dpQuadID = h.detectorID-55;
	      barID = h.elementID;
	      if (h.detectorID >= 55 && h.detectorID <= 58){
		y1 = 7.5 + (barID + .5)*wDP1;
	      }
	      if (h.detectorID >= 59 && h.detectorID <= 62){
                y1 = 7.5 + (barID + .5)*wDP2;
              }
	    }
	    else if (h.detectorID == 43) { // H4Y2L
	      if (h.elementID > 8) {
		dpQuadID = 8;
		barID = h.elementID-8;
	      } else {
		dpQuadID = 10;
		barID = 9-h.elementID;
	      }
	      //barID [0,8] should be the ones to use when computing y and is saved in hitsets
	      //elmid or hits[] ranges from [0,16] and should be used 
	      // hits have elmID and hitsets have barIDs
	      y1 = (barID+.5)*wH4Y;
	    }
	    else if (h.detectorID == 44) { // H4Y2R
	      if (h.elementID > 8) {
		dpQuadID = 9;
		barID = h.elementID-8;
	      } else {
		dpQuadID = 11;
		barID = 9-h.elementID;
	      }
              y1 = (barID+.5)*wH4Y;
	    }

	    // select only hits with h4y2 or dp hits     
	    if(dpQuadID<0 || dpQuadID>11) continue; 

	    // check time
	    double deltaT = h.tdcTime - getDelay(dpQuadID,barID);
	    //tweak deltaT so the NIM triggers are in time (the delays are tuned for the DP trigger)
	    if (rawEvent->getTriggerBits() == nim1TriggerMask) {
	      deltaT -= 5.0;
	    } else if (rawEvent->getTriggerBits() == nim3TriggerMask) {
	      deltaT -= 15.0;
	    }

	    //deltat[k] = deltaT;
	    //if (h.isInTime()) hit_InTime[k] = 1;
	    //else hit_InTime[k] = 0;

	    // making hitsets: there are hitsets 0 - 11 depending on dpquadid filled with the hit barID
	    // 0-7 dp and 8-11 h4y
	    hitsetsNoTimecut[dpQuadID].insert(barID);
	    if (TMath::Abs(deltaT)<timeCut){
	      hitsets[dpQuadID].insert(barID);
	      hitsy1[dpQuadID].push_back(y1);
	      hits[dpQuadID].push_back(h.elementID); // this is the actual barID
	      hitsbars[dpQuadID].push_back(barID); // this is the transformed barID ( only diff for H4Y)
	    }
	    if (TMath::Abs(deltaT+rfclock)<timeCut) hitsetsOffByOne[dpQuadID].insert(barID);

	  } // end hit loop

	  // hitsets is just the hits in each quadrant and the barIDs where you have hits -> per event
	  // now we have collected hits in each quadrant
	  // I want to see if the pattern matches the one in the roads file
	  int roadBits = 0; 
	  int roadBitsNoTimecut = 0;
	  int nimlikeBits = 0;

	  int counthit(0);
	  int testbit(-1), testbit2(-1), testbit3(-1), testbit4(-1);
	  int test(-1);
	  
	  vector<int> quads;
          vector<int> vquad,vquad_h4y,vbaridDP1,vbaridDP2,vbaridH4Y,vbaridH4Y_t;
	  vector<float> vyDP1,vyDP2,vyH4Y,vyh;

	  for (int quad=0;quad<4;quad++) {
	    // so here you are always taking same quad hits
	    int fID = (quad)%4; //forward dp
	    int bID = quad+4; // back dp
	    int hID = quad+8; // h4ys

	    //one beam-left and one beam-right quadrant 
	    // fID: 0 1 2 3
	    // fID: 55 56 57 58 
	    // bID: 4 5 6 7
	    // bID: 59 60 61 62 e.g. 62-55=7
	    // bID: TL TR BL BR
	    // bID: 59
	    // hID: 8 9 10 11
	    // hID:H4 TL TR BL BR

	    // one beam left: 0 4 8 or 2 6 10
	    // one beam right: 1 5 9 or 3 7 11

	    // roadbits
	    // 0101 0110 0111 1001 1010 1011 1101 1110 1111
	    // at least one left and one right - no this is not right
	    // it should be: at least 2 different quadrants with the 3 hits in the DP1,DP2 and H4Y in the same quad and matching to one of dylan's roads
	    // 5: 0 and 2: TL - BL
	    // 6: 1 and 2: TR - BL
	    // 7: 0 1 and 2: TL - TR - BL
	    // 9: 0 and 3: TL - BR
	    // 10: 1 and 3: TR - BR
	    // 11: 0 1 and 3: TL - TR - BR
	    // 13: 0 2 and 3: TL - BL - BR
	    // 14: 1 2 and 3: TR - BL - BR
	    // 15: 0 1 2 and 3: TL - TR - BL - BR

	    // testbit: at least 2 different quadrants with the 3 hits in the DP1,DP2 and H4Y in the same quad and matching to one of dylan's roads but one T and one B
	    // testbit2:  at least 2 different quadrants with the 3 hits in the DP1,DP2 and H4Y in the same quad and matching to one of dylan's roads

	    // quad 0 TL 0 4 8  
            // quad 1 TR 1 5 9 
            // quad 2 BL 2 6 10 
            // quad 3 BR 3 7 11 
	    //cout << " position size " << hitsetsy1.size() << " for fId " << fID << " is " hitsetsy1[fID].size() << endl;
	    // now we are looping over the barIDs
	    //cout << " size set " << hitsets[fID].size() << " vec " << hits[fID].size() << endl; same size!
	    for (set<int>::iterator fIt = hitsets[fID].begin(); fIt!=hitsets[fID].end();++fIt) {
	      for (set<int>::iterator bIt = hitsets[bID].begin(); bIt!=hitsets[bID].end();++bIt) {
		for (set<int>::iterator hIt = hitsets[hID].begin(); hIt!=hitsets[hID].end();++hIt) {
		  int roadhash = roadHash(*fIt,*bIt,*hIt); // for the hits withe st1*1000 + st2*10 + h4
		  // e..g 2 17 2  => 2*1000 + 17*10 + 2 = 2172 barid 2 17 2
		  if (roadset.count(roadhash)) {
		    if(counthit>0){
		      if(testbit==test){ // until now all of them have been in the same quadrant
			if( (quad==2 || quad==3) && (test==0 || test==1)){
			  testbit = quad;
			} 
			if( (quad==0 || quad==1) && (test==2 || test==3)){
			  testbit = quad;
			}
		      }
		      if(testbit2==test){
			testbit2 = quad;
		      }
		      else{
			testbit2 = -1;
		      }
		    }
		    if(counthit==0){
		      test = quad;
		      testbit= quad;
		      testbit2=quad;
		    }
		    counthit+=1;
		    roadBits |= (1 << quad); // << shift bits left by quad 0001 => 0010 or 0100 or 1
		    quads.push_back(quad);
		    //Assign x | y to roadBits
		    // 0 or 0001 => 0001
		    // 0001 or 0010 => 0011
		  }
		}
	      }
	    }

	    // now loop over the hits
	    for(int if0 = 0; if0 < int(hits[fID].size()); if0++) { 
	      for(int ib0 = 0; ib0 < int(hits[bID].size()); ib0++) {
		for(int ih0 = 0; ih0 < int(hits[hID].size()); ih0++) {
		  int fIt = hitsbars[fID][if0]; float fIty = hitsy1[fID][if0];
                  int bIt = hitsbars[bID][ib0]; float bIty = hitsy1[bID][ib0];
                  int hIt = hitsbars[hID][ih0]; float hIty = hitsy1[hID][ih0];
                  int roadhash = roadHash(fIt,bIt,hIt);

		  if(fIt>30){
		    std::cout << " hit bar " << hIt-1 << " quad " << quad  << " fit " << fIt << " bIt " << bIt << std::endl;
		  }
		  if (rawEvent->getTriggerBits() == nim3TriggerMask) {
		    DP1DP2Hists_all[hIt-1]->Fill(fIt,bIt);
                    DP1DP2Hists[hIt-1][quad]->Fill(fIt,bIt);
		  }

                  if (roadset.count(roadhash)) {
		    vquad_data_samequad.push_back({fIt,bIt,hIt});
		    //std::cout << "y hit positions " << fIty << " " << fIt << " " << if0 << " " << bIty << " " << bIt << " " << hIty << " " << hIt << std::endl;
		    
		    // need to work with all ys and x in the first quad  (all positives)
		    float m = (bIty-fIty)/(xDP2-xDP1);
		    float b = bIty-m*xDP2;
		    float yh = m*xH4Y + b;
		    //cout << "yhits: " << fIty << " " << bIty << " " << hIty << " if in line " <<  yh << endl;
		    //cout << "ybars: " << fIt << " " << bIt << " " << hIt << " hbar " << hitsbars[hID][ih0] << " predict " << findh4ybar(hIty) << " if in line " << findh4ybar(yh) << endl;
		    //if(findh4ybar(yh) >= findh4ybar(hIty)-1 && findh4ybar(yh) <= findh4ybar(hIty)+1){
		    if(findh4ybar(yh) == findh4ybar(hIty)){ 
		      vquad_data_h4y.push_back({fIt,bIt,hIt});
		    //if(findh4ybar(yh) == findh4ybar(hIty)){
		      testbit4 = 1;
		      //cout << "yhits: " << fIty << " " <<bIty <<" " << hIty << " if in line " <<  yh <<endl;
		      //cout << "ybars: " << fIt <<" " << bIt << " " << hIt << " hbar " << hitsbars[hID][ih0] << " predict " << findh4ybar(hIty) << " if in line " << findh4ybar(yh) << endl;
		      //cout << " in line " << findh4ybar(hIty)-2 << " "<<findh4ybar(yh) <<" "<<  findh4ybar(hIty)+2 << endl;
		      //cout << " slope!! " << m << " intercept " << b << " yH4Y " << yh << " yhith4y " << hIty << endl;
		      //cout << " true-barh4 " << hIt << " fromy-barh4 " << findh4ybar(yh) << " fromtrue-barh4 " << findh4ybar(hIty) << endl;
		      vquad_h4y.push_back(quad);
		    }
		    //cout << " bars " << hitsbars[hID][if0] << " id " << hIt << endl;
		    vquad.push_back(quad);
		    vbaridDP1.push_back(fIt);
                    vbaridDP2.push_back(bIt);
                    vbaridH4Y.push_back(hIt);
		    vyDP1.push_back(fIty);
                    vyDP2.push_back(bIty);
                    vyH4Y.push_back(hIty);
		    vyh.push_back(yh);
		  }
		}
	      }
	    } // end loop over hits

	    for (set<int>::iterator fIt = hitsetsNoTimecut[fID].begin(); fIt!=hitsetsNoTimecut[fID].end();++fIt) {
	      for (set<int>::iterator bIt = hitsetsNoTimecut[bID].begin(); bIt!=hitsetsNoTimecut[bID].end();++bIt) {
		for (set<int>::iterator hIt = hitsetsNoTimecut[hID].begin(); hIt!=hitsetsNoTimecut[hID].end();++hIt) {
		  int roadhash = roadHash(*fIt,*bIt,*hIt);
		  if (roadset.count(roadhash)) {
		    roadBitsNoTimecut |= (1 << quad);
		  }
		}
	      }
	    }
	    if (!hitsets[bID].empty()) nimlikeBits |= (1 << quad);
	  } // end loop over quads

	  // save info in a tree
	  for(int i0 = 0; i0 < 10; i0++) {
	    quadrant[i0] = -1;
	    baridDP1[i0] = -1;
	    baridDP2[i0] = -1;
	    baridH4Y[i0] = -1;
	    yDP1[i0] = -1;
	    yDP2[i0] = -1;
	    yH4Y[i0] = -1;
	    yh[i0] = -1;
	  }
	  nRoads = int(vquad.size());
	  
	  // save tree and check if its top and bottom
	  int firsttb(-1),tmptb(-1);
	  int istb(-1),istb_h4y(-1),issamequad(-1),issamequad_h4y(-1),is2diffquad(-1),is2diffquad_h4y(-1),is2diffhits(-1),is2diffhits_h4y(-1);
	  int is2samequad(-1),is2samequad_h4y(-1);
	  if(vquad.size()>0) {
	    //std::cout << "nquads " << int(vquad.size()) << " ";
	    for(int i0 = 0; i0 < int(vquad.size()); i0++) {      
	      //std::cout << vquad[i0];
	      if(i0<10){
		quadrant[i0] = vquad[i0];
		baridDP1[i0] = vbaridDP1[i0];
		baridDP2[i0] = vbaridDP2[i0];
		baridH4Y[i0] = vbaridH4Y[i0];
		yDP1[i0] = vyDP1[i0];
		yDP2[i0] = vyDP2[i0];
		yH4Y[i0] = vyH4Y[i0];
		yh[i0] = vyh[i0];
	      }
	      if(i0==0){// check first quad
		if(vquad[i0]==0 || vquad[i0]==1) firsttb = 0; //top
		else if(vquad[i0]==2 || vquad[i0]==3) firsttb = 1; //bottom
	      }
	      else{
		if(vquad[i0]==0 || vquad[i0]==1) tmptb = 0; //top
		else if(vquad[i0]==2 || vquad[i0]==3) tmptb = 1; //bottom
		if((firsttb==0 && tmptb==1)||(firsttb==1 && tmptb==0)) istb+=1;
	      }
	    }

	    // check for istb
	    /*
	    if(istb>=0){
	      cout << "istb " << endl;
	      for(int i0 = 0; i0 < int(vquad.size()); i0++) {
	        cout << "quad " << vquad[i0] << endl;  
		cout << "yhits: " << vyDP1[i0] << " " <<vyDP2[i0] <<" " << vyH4Y[i0] << " if in line " <<  vyh[i0] <<endl;
		cout << "ybars: " << vbaridDP1[i0] <<" " << vbaridDP2[i0] << " " << vbaridH4Y[i0]  << " predict " << findh4ybar(vyH4Y[i0]) << " if in line " << findh4ybar(vyh[i0]) << endl;
	      }
	    }
	    */

	    //std::cout << " \n " << std::endl;
	    //std::cout << "and nhits dp1 2 h4 " << std::endl;
	    //for(int i0 = 0; i0 < 4; i0++){
	    //std::cout << hits[i0%4].size() << " " << hits[i0+4].size() << " " << hits[i0+8].size() << std::endl;
	    //}
	    //std::cout << " \n " << std::endl;

	    // count at least two same quads
	    std::vector<int> myotherquads;
	    for(int i0 = 0; i0 < 4; i0++){
	      if(i0==0) { myotherquads = {1,2,3};}
	      else if(i0==1) { myotherquads = {0,2,3};}
	      else if(i0==2) { myotherquads = {1,0,3};}
	      else if(i0==3) { myotherquads = {1,2,0};}
	      
	      int q0 = std::count (vquad.begin(), vquad.end(), i0);
	      // roads with same quad
	      if(q0>1) {
		issamequad += 1;
		//std::cout << "same quad " << i0 << std::endl;
	      }
	      
	      // exactly 2 roads with same quad
              if(q0==2) {
                is2samequad += 1;
	      }

	      // 2 diff quads but at least 2 hits in each detector 
	      if(q0>0 && hits[i0%4].size()>1 && hits[i0+4].size()>1 && hits[i0+8].size()>1){
		for(int i1 = 0; i1 < myotherquads.size(); i1++){
		  if(std::count (vquad.begin(), vquad.end(), myotherquads[i1])>0 && hits[i1%4].size()>1 && hits[i1+4].size()>1 && hits[i1+8].size()>1){
		    is2diffquad = 1;
		    //std::cout << "diff quad " << i0 << " and " << i1 << " but 2 hits in each" << std::endl;
		  }
		}
	      }

	      // at least 2 roads in same quad but at least 2 hits in each detector
	      if(q0>1 && hits[i0%4].size()>1 && hits[i0+4].size()>1 && hits[i0+8].size()>1){
		is2diffhits += 1;
		//std::cout << "same quad " << i0 << " but 2 hits in each" << std::endl;
	      }
	      
	    }
	    //if((issamequad>=1 && is2diffhits<1) || (issamequad<1 &&is2diffhits>=1)){
	    //  std::cout << " SAMES " << issamequad << " diffhits " << is2diffhits << std::endl;
	    //}
	    // if(is2samequad==1){
	    // std::cout << " exactly 2 roads with same quad " << vquad.size() << std::endl;
	    //}

	    // h4y in line
	    firsttb=-1;tmptb=-1;
	    //cout << "size " << vquad.size() << " in line " << vquad_h4y.size() << endl;
	    for(int i0 = 0; i0 < int(vquad_h4y.size()); i0++) {
	      if(i0==0){
		if(vquad_h4y[i0]==0 || vquad_h4y[i0]==1) firsttb = 0; 
		else if(vquad_h4y[i0]==2 || vquad_h4y[i0]==3) firsttb = 1;
	      }
	      else{
		if(vquad_h4y[i0]==0 || vquad_h4y[i0]==1) tmptb = 0; 
		else if(vquad_h4y[i0]==2 || vquad_h4y[i0]==3) tmptb = 1; 
		if((firsttb==0 && tmptb==1)||(firsttb==1 && tmptb==0)) istb_h4y+=1;
	      }
	    }
	    
	    // count at least two same quads
	    if(vquad_h4y.size()>0){
	      for(int i0 = 0; i0 < 4; i0++){
		if(i0==0) { myotherquads = {1,2,3};}
		else if(i0==1) { myotherquads = {0,2,3};}
		else if(i0==2) { myotherquads = {1,0,3};}
		else if(i0==3) { myotherquads = {1,2,0};}
	      
		int q0 = std::count (vquad_h4y.begin(), vquad_h4y.end(), i0);
		//std::cout<< " counts in quad " << i0 << " is " << q0 << std::endl;

		// at least 2 roads with same quad
		if(q0>1) issamequad_h4y += 1;
		
		if(q0==2) {
		  is2samequad_h4y += 1;
		}
		// 2 diff quads but at least 2 hits in each detector
		if(q0>0 && hits[i0%4].size()>1 && hits[i0+4].size()>1 && hits[i0+8].size()>1){
		  for(int i1 = 0; i1 < myotherquads.size(); i1++){
		    if(std::count (vquad_h4y.begin(), vquad_h4y.end(), myotherquads[i1])>0 && hits[i1%4].size()>1 && hits[i1+4].size()>1 && hits[i1+8].size()>1)
		      is2diffquad_h4y = 1;
		  }
		}
		
		// at least 2 roads in same quad but at least 2 hits in each detector                                                                                                                            
		if(q0>1 && hits[i0%4].size()>1 && hits[i0+4].size()>1 && hits[i0+8].size()>1){
		  is2diffhits_h4y += 1;
		}
	      }
	      /*
	      std::cout << " size of h4y matched " << vquad_h4y.size() << std::endl;
	      std::cout << " quadh4y ";
	      for(int i1 = 0; i1 < int(vquad_h4y.size()); i1++) {
		std::cout << vquad_h4y[i1] << " ";
	      }
	      std::cout<< " \n ";
	      if(issamequad_h4y>=0){
		std::cout<< " at least 2 roads with same quad " << issamequad_h4y << std::endl;
	      }
	      if(is2samequad_h4y>=0){
		std::cout<< " exactly 2 roads with same quad " << is2samequad_h4y << std::endl;
	      }
	      */
	    }
	  }	

	  // *** new check ***
	  if(rawEvent->getTriggerBits() == dpTriggerMask){
	    if(roadBits>0 && vquad.size()>=2){
	      nMatchedDP+=1;

	      if(issamequad>=0){ // at least 2 roads with same quad
		nSameQuadDP +=1;
		if(is2samequad==0) nSame2QuadDP+=1; // exactly 2 roads with same quad => just 1 quad that matches 2 same roads
	      }
	      else nDiffQuadDP +=1; // no same quad roads

	      if(vquad_h4y.size()>=2){
		nMatchedDPlineH4Y+=1;
		if(issamequad_h4y>=0){
		  nSameQuadDPlineH4Y+=1;
		  if(is2samequad_h4y==0) nSame2QuadDPlineH4Y+=1;
		}
		else nDiffQuadDPlineH4Y +=1;
	      }

	      if(istb>=0) {
		nTBQuadDP+=1;
		if(istb_h4y>=0) nTBQuadDPlineH4Y+=1;
	      }

	      if(is2diffhits>=1) {
		nDiff2hitsDP+=1;
		if(is2diffhits_h4y>=1) nDiff2hitsDPlineH4Y+1;
	      }
	      if(is2diffquad==1) {
		nDiff2quadDP+=1;
                if(is2diffquad_h4y==1) nDiff2quadDPlineH4Y+1;
	      }
	    }
          }
	  if(rawEvent->getTriggerBits() == nim3TriggerMask){
            if(roadBits>0 && vquad.size()>=2){
              nMatchedNIM3+=1;

              if(issamequad>=0){
                nSameQuadNIM3 +=1;
                if(is2samequad==0) nSame2QuadNIM3+=1;
              }
              else nDiffQuadNIM3 +=1;

	      if(vquad_h4y.size()>=2){
		nMatchedNIM3lineH4Y+=1;
		if(issamequad_h4y>=0){
		  nSameQuadNIM3lineH4Y+=1;
		  if(is2samequad_h4y==0) nSame2QuadNIM3lineH4Y+=1;
		}
		else nDiffQuadNIM3lineH4Y +=1;
              }

              if(istb>=0) {
                nTBQuadNIM3+=1;
                if(istb_h4y>=0) nTBQuadNIM3lineH4Y+=1;
	      } 
	      if(is2diffhits>=1) {
                nDiff2hitsNIM3+=1;
                if(is2diffhits_h4y>=1) nDiff2hitsNIM3lineH4Y+1;
              }
	      if(is2diffquad==1) {
		nDiff2quadNIM3+=1;
		if(is2diffquad_h4y==1) nDiff2quadNIM3lineH4Y+1;
	      }
            }
          }

	  // select hits
	  // at least one top and one bottom quad
	  if(roadBits>0 && testbit>-1 && testbit!=test && quads.size()>1){
	    if(rawEvent->getTriggerBits() == dpTriggerMask) nBitsTrigDPCheck +=1;
	    if(rawEvent->getTriggerBits() == nim3TriggerMask) nBitsTrigNIM3Check +=1;
	  }
	  // at least two different quad
	  if(roadBits>0 && testbit2>-1 && testbit2!=test && quads.size()>1){
            if(rawEvent->getTriggerBits() == dpTriggerMask) nBitsTrigDPCheck2 +=1;
            if(rawEvent->getTriggerBits() == nim3TriggerMask) nBitsTrigNIM3Check2 +=1;
          }
	  // all of them in the same quad   
	  if(roadBits>0 && testbit2>-1 && testbit2==test && quads.size()>1){
            if(rawEvent->getTriggerBits() == dpTriggerMask) {
	      //cout << " testbit2 " << testbit2 << "  test " << test << std::endl;
	      //for(int i0 = 0; i0 < int(quads.size()); i0++) {
	      //	cout << " quads in same quad " << quads[i0] << endl;
	      //}
	      nBitsTrigDPCheck3 +=1;
	    }
            if(rawEvent->getTriggerBits() == nim3TriggerMask) nBitsTrigNIM3Check3 +=1;
          }
	  // all of them in the same quad but only 2 quads ( not 4..)
	  int quad0(0),quad1(0),quad2(0),quad3(0);
	  if(roadBits>0 && quads.size()>1){
	    quad0 = std::count (quads.begin(), quads.end(), 0);
            quad1 = std::count (quads.begin(), quads.end(), 1);
            quad2 = std::count (quads.begin(), quads.end(), 2);
            quad3 = std::count (quads.begin(), quads.end(), 3);
	    //cout << "quad 0 appears " << std::count (quads.begin(), quads.end(), 0)  << " times.\n";
            //cout << "quad 1 appears " << std::count (quads.begin(), quads.end(), 1)  << " times.\n";
            //cout << "quad 2 appears " << std::count (quads.begin(), quads.end(), 2)  << " times.\n";
            //cout << "quad 3 appears " << std::count (quads.begin(), quads.end(), 3)  << " times.\n";
	    if((quad0==2 && (quad1==0&&quad2==0&&quad3==0)) || (quad1==2 && (quad2==0&&quad2==0&&quad3==0)) || (quad2==2 && (quad1==0&&quad0==0&&quad3==0)) || (quad3==2 && (quad1==0&&quad2==0&&quad0==0))){
	      if(rawEvent->getTriggerBits() == dpTriggerMask) nBitsTrigDPCheck4 +=1;
	      if(rawEvent->getTriggerBits() == nim3TriggerMask) nBitsTrigNIM3Check4 +=1;
	    }
          }
	  
	  // h4y is in line +/- 2 bars
	  if(roadBits>0 && quads.size()>1 && testbit4>-1){
	    if(rawEvent->getTriggerBits() == dpTriggerMask) nBitsTrigDPCheck5 +=1;
	    if(rawEvent->getTriggerBits() == nim3TriggerMask) nBitsTrigNIM3Check5 +=1;
	  }

	  if(roadBits>0 && quads.size()>1 && testbit4>-1 && testbit>-1 && testbit!=test){
            if(rawEvent->getTriggerBits() == dpTriggerMask) nBitsTrigDPCheck6 +=1;
            if(rawEvent->getTriggerBits() == nim3TriggerMask) {
	      nBitsTrigNIM3Check6 +=1;
	      /*
	      std::cout << "testbit4 " << testbit4 <<  " first quad " << test << " and testbit " << testbit << std::endl;
	      std::cout << " MATCHED TO 9 " << std::endl;
	      std::cout << " size of h4y matched " << vquad_h4y.size() << std::endl;
	      std::cout << " quadh4y ";
              for(int i1 = 0; i1 < int(vquad_h4y.size()); i1++) {
		std::cout << vquad_h4y[i1] << " ";
              }
	      for(int i1 = 0; i1 < 4; i1++){
		std::cout << " counts in quad " << i1  << " IS " << std::count (vquad_h4y.begin(), vquad_h4y.end(), i1) << std::endl;
	      }
	      std::cout << " size of matched " << vquad.size() << std::endl;
	      std::cout << " quad ";
	      for(int i1 = 0; i1 < int(vquad.size()); i1++) {
		std::cout << vquad[i1] << " ";
		cout << "yhits: " << vyDP1[i1] << " " <<vyDP2[i1] <<" " << vyH4Y[i1] << " if in line " <<  vyh[i1] <<endl;
		cout << "ybars: " << vbaridDP1[i1] <<" " << vbaridDP2[i1] << " " << vbaridH4Y[i1]  << " predict " << findh4ybar(vyH4Y[i1]) << " if in line " << findh4ybar(vyh[i1]) << endl; 
              }
	      */
	    }
          }

	  // old
	  roadsVsTrig->Fill(rawEvent->getTriggerBits(),roadBits);
          roadsVsTrigNoTimecut->Fill(rawEvent->getTriggerBits(),roadBitsNoTimecut);
          nimlikeVsTrig->Fill(rawEvent->getTriggerBits(),nimlikeBits);

	  if(roadBits == 5 || roadBits == 6 || roadBits == 7 || roadBits == 9 || roadBits == 10 || roadBits == 11 || roadBits == 13 || roadBits == 14 || roadBits == 15){
	    if (rawEvent->getTriggerBits() == dpTriggerMask){
		nBitsTrigDP+=1;
	    }
	    if (rawEvent->getTriggerBits() == nim3TriggerMask){
              nBitsTrigNIM3+=1;
            }
	  }

	  int onlyquadF = -1;
	  int onlyquadB = -1;
	  int onlyquadH = -1;
	  int nhitsF, nhitsB, nhitsH;
	  int nquadsF = 0;
	  int nquadsB = 0;
	  int nquadsH = 0;
	  for (int quad=0;quad<4;quad++) {
	    int fID = (quad)%4;
	    int bID = quad+4;
	    int hID = quad+8;
	    int fHits = hitsets[fID].size();
	    int bHits = hitsets[bID].size();
	    int hHits = hitsets[hID].size();
	    if (fHits>0) {
	      nquadsF++;
	      if (onlyquadF==-1) {
		onlyquadF = quad;
		nhitsF = fHits;
	      } else
		onlyquadF = -2;
	    }
	    if (bHits>0) {
	      nquadsB++;
	      if (onlyquadB==-1) {
		onlyquadB = quad;
		nhitsB = bHits;
	      } else
		onlyquadB = -2;
	    }
	    if (hHits>0) {
	      nquadsH++;
	      if (onlyquadH==-1) {
		onlyquadH = quad;
		nhitsH = hHits;
	      } else
		onlyquadH = -2;
	    }
	  }
	  
	  saveTree->Fill();
	
	  // clearing hitsets
	  for (int j=0;j<12;j++) {
	    hitsets[j].clear();
	    hitsetsNoTimecut[j].clear();
	    hitsetsOffByOne[j].clear();
	    hits[j].clear();
	    hitsbars[j].clear();
	    hitsy1[j].clear();
	  }
	}
    
	rawEvent->clear();
    }

    // write new roads to file
    ofstream roads_samequad;
    roads_samequad.open("roads_samequad.txt");
    for (int i0=0;i0<int(vquad_data_samequad.size());i0++) {
      roads_samequad << vquad_data_samequad[i0].x << " " << vquad_data_samequad[i0].y << "  " << vquad_data_samequad[i0].z << "\n";
    }
    roads_samequad.close();

    ofstream roads_h4y;
    roads_h4y.open("roads_matchh4y.txt");
    for (int i0=0;i0<int(vquad_data_h4y.size());i0++) {
      roads_h4y <<vquad_data_h4y[i0].x << " " << vquad_data_h4y[i0].y << "  " << vquad_data_h4y[i0].z << "\n";
    }
    roads_h4y.close();

    cout << "rates " << endl;
    cout << "L1 DP: " << nTrigDP << " NIM3 " << nTrigNIM3 << endl;
    cout << "DP fires: " << nBitsTrigDP << " NIM3 fires: " << nBitsTrigNIM3 << endl;
    cout << "\n" << endl;

    cout << "all matched roads will contain 3 hits (DP1,DP2, H4Y) in the same quadrant" << endl;
    cout << " DP NIM3 " << endl;
    cout << "at least 2 roads matched per evt =>" << endl;
    cout << " " << nMatchedDP << " "  << nMatchedNIM3 << endl;
    cout << "at least 2 roads matched per evt with same quad" << endl;
    cout << " " << nSameQuadDP << " "  << nSameQuadNIM3 << endl;
    cout << "only 2 roads matched per evt with same quad "<< endl;
    cout << " " << nSame2QuadDP << " "  << nSame2QuadNIM3 << endl; 
    cout << "at least 2 roads matched per evt with different quad " << endl;
    cout << " " << nDiffQuadDP << " "  << nDiffQuadNIM3 << endl;
    cout << "at least 2 roads matched per evt with one top quad and one bottom quad " << endl;
    cout << " " << nTBQuadDP << " "  << nTBQuadNIM3 << endl;
    cout << "at least 2 roads matched per evt in same quad but at least 2 hits in each detector " << endl;
    cout << " " << nDiff2hitsDP << " " << nDiff2hitsNIM3 << endl;
    cout << "at least 2 roads matched in diff quads but at least 2 hits in each detector " << endl;
    cout << " " << nDiff2quadDP << " " << nDiff2quadNIM3 << endl;
    cout << "\n" << endl;

    cout << "all matched roads will contain 3 hits (DP1,DP2, H4Y) in the same quadrant and the H4Y hit will be in line"<< endl;
    cout << " DP NIM3 " << endl;
    cout << "at least 2 roads matched per evt" << endl;
    cout << " " << nMatchedDPlineH4Y << "  "  << nMatchedNIM3lineH4Y << endl;
    cout << "at least 2 roads matched per evt with same quad" << endl;
    cout << "  " << nSameQuadDPlineH4Y << "  "  << nSameQuadNIM3lineH4Y << endl;
    cout << "only 2 roads matched per evt with same quad "<< endl;
    cout << "  " << nSame2QuadDPlineH4Y << "  "  << nSame2QuadNIM3lineH4Y << endl;
    cout << "at least 2 roads matched per evt with different quad " << endl;
    cout << "  " << nDiffQuadDPlineH4Y << "  "  << nDiffQuadNIM3lineH4Y << endl;
    cout << "at least 2 roads matched per evt with one top quad and one bottom quad " << endl;
    cout << "  " << nTBQuadDPlineH4Y << "  "  << nTBQuadNIM3lineH4Y << endl;
    cout << "at least 2 roads matched per evt in same quad but at least 2 hits in each detector " << endl;
    cout << " " << nDiff2hitsDPlineH4Y << " " << nDiff2hitsNIM3lineH4Y << endl;
    cout << "at least 2 roads matched in diff quads but at least 2 hits in each detector " << endl;
    cout << " " << nDiff2quadDPlineH4Y << " " << nDiff2quadNIM3lineH4Y << endl;
    cout << "\n" << endl;

    cout << "OLD " << endl;
    cout << "at least 2 roads with diff quad T and B => DP fires: " << nBitsTrigDPCheck << " NIM3 fires " << nBitsTrigNIM3Check << endl;
    cout << "at least 2 roads with diff quad => DP fires: " << nBitsTrigDPCheck2 << " NIM3 fires " << nBitsTrigNIM3Check2 << endl; // why is this less??
    cout << "at least 2 roads in same quad (6 hits) => DP fires: " << nBitsTrigDPCheck3 << " NIM3 fires " << nBitsTrigNIM3Check3 << endl;
    cout << "only 2 roads in same quad (6 hits) => DP fires: " << nBitsTrigDPCheck4 << " NIM3 fires " << nBitsTrigNIM3Check4 << endl;
    cout << "only hits with h4y in line => DP fires: " << nBitsTrigDPCheck5 << " NIM3 fires " << nBitsTrigNIM3Check5 << endl;
    cout << "only hits with h4y in line and at least 2 diff quad T and B  => DP fires: " << nBitsTrigDPCheck6 << " NIM3 fires " << nBitsTrigNIM3Check6 << endl;

    saveTree->Write();
    saveFile->Write();
    saveFile->Close();
    
    return EXIT_SUCCESS;
}
