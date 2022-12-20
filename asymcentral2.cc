#include "Pythia8/Pythia.h"

#include "Pythia8/HeavyIons.h"
// ROOT, for saving Pythia events as trees in a file.
#include <stdio.h>
#include <stdlib.h>
#include "TH1.h"
#include "TH2.h"
#include "TVirtualPad.h"
#include "TApplication.h"
#include <vector> 
#include "TTree.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGraph.h"
#include <bits/stdc++.h>


using namespace Pythia8;


int main() {

  //const Int_t n = 100000;

  //Int_t ntrack;

  Pythia pythia;
 
  pythia.readString("Beams:idA = 1000822080");
  pythia.readString("Beams:idB = 1000822080"); // The lead ion.
  pythia.readString("Beams:eCM = 2760.0");
  pythia.readString("Beams:frameType = 1"); 
  
 
  pythia.readString("Random:seed=12343669");
  // 
   
  

 //pythia.readString("PhaseSpace:pTHatMin = 20.");
 
 
 pythia.settings.listAll(); 
 pythia.init();
 
 
 TFile* outFile = new TFile("lastcentrality.root", "RECREATE");
 TH1F *mcen = new TH1F("M1","Total number of final state charged particles", 50,0,7000);
 
 TH1F *pt1b = new TH1F("pt1b","pt of region 1, 0-20%",100,0,50);
 TH1F *pt2b = new TH1F("pt2b","pt of region 2, 0-20%",100,0,50);
 TH1F *pt3b = new TH1F("pt3b","pt of region 3, 0-20%",100,0,50);
 TH1F *eta1b = new TH1F("eta1b","eta of region 1, 0-20%",20,-4,4);
 TH1F *eta2b = new TH1F("eta2b","eta of region 2, 0-20%",20,-4,4);
 TH1F *eta3b = new TH1F("eta3b","eta of region 3, 0-20%",20,-4,4);
 TH1F *phi1b = new TH1F("phi1b","phi of region 1, 0-20%",20,-M_PI,M_PI);
 TH1F *phi2b = new TH1F("phi2b","phi of region 2, 0-200%",20,-M_PI,M_PI);
 TH1F *phi3b = new TH1F("phi3b","phi of region 3, 0-20%",20,-M_PI,M_PI);
 
 
  int nEvent=1000;
 
   Event& event = pythia.event;

  for ( int iEvent = 0; iEvent < nEvent; ++iEvent ) {
    if ( !pythia.next() ) continue;
    
    int abstarget=pythia.info.hiInfo->nAbsTarg();
    int difftarget=pythia.info.hiInfo->nDiffTarg();
    double A = (double)(abstarget+difftarget);
    int absproj=pythia.info.hiInfo->nAbsProj();
    int diffproj=pythia.info.hiInfo->nDiffProj();
    double B = (double)(absproj+diffproj);
    double alpha = (double)((A-B)/(A+B));
    double yi=(1.0/2.0)*TMath::Log(A/B);
    double nCharged = 0.0;
    for (int i = 0; i < pythia.event.size(); ++i){
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged() && abs(pythia.event[i].eta()) < 2.5)
      	{
	        nCharged=nCharged+1.0;
	      }
    }
    
    if(nCharged>3979.7 && nCharged<7000.0){
    if(alpha<-0.1)
    {
        for(int i=0;i<pythia.event.size();i++)
        {
          if (pythia.event[i].isFinal() && pythia.event[i].isCharged() && abs(pythia.event[i].eta()) < 2.5){
          pt1b->Fill(pythia.event[i].pT());
          eta1b->Fill(pythia.event[i].eta());
          phi1b->Fill(pythia.event[i].phi());
          }
        }
    }
    if(alpha>0.1)
    {
        for(int i=0;i<pythia.event.size();i++)
        {
          if (pythia.event[i].isFinal() && pythia.event[i].isCharged() && abs(pythia.event[i].eta()) < 2.5){
          pt2b->Fill(pythia.event[i].pT());
          eta2b->Fill(pythia.event[i].eta());
          phi2b->Fill(pythia.event[i].phi());
          }
        }
    }
    if(alpha>-0.1 && alpha<0.1)
    {
        for(int i=0;i<pythia.event.size();i++)
        {
          if (pythia.event[i].isFinal() && pythia.event[i].isCharged() && abs(pythia.event[i].eta()) < 2.5){
          pt3b->Fill(pythia.event[i].pT());
          eta3b->Fill(pythia.event[i].eta());
          phi3b->Fill(pythia.event[i].phi());
          }
        }

    }
    }
  }

 double ret1b, ret2b, rp1b, rp2b, rph1b, rph2b;

  std::ofstream myfile;
  myfile.open ("pseudorapidityData.csv");
    for(int i=1;i<=20;i++){
      double et1b=eta1b->GetBinContent(i);
      double et2b=eta2b->GetBinContent(i);
      double et3b=eta3b->GetBinContent(i);
      if(et1b==0 || et2b==0 || et3b==0) continue;
      ret1b = et1b/et3b;
      ret2b = et2b/et3b;
      myfile << setprecision(9);
    myfile << ret1b << "," << ret2b << "\n";
    }
    myfile.close();

  myfile.open ("pTdata.csv");
    for(int i=1;i<=100;i++){
      double p1b= pt1b->GetBinContent(i);
      double p2b= pt2b->GetBinContent(i);
      double p3b= pt3b->GetBinContent(i);
      if(p1b==0 || p2b==0 || p3b==0) continue;
      rp1b= p1b/p3b;
      rp2b = p2b/p3b;
      myfile << setprecision(9);
    myfile << rp1b << "," << rp2b << "\n";
    }
    myfile.close();
   
  myfile.open ("phidata.csv");
    for(int i=1;i<=20;i++){
      double ph1b= phi1b->GetBinContent(i);
      double ph2b= phi2b->GetBinContent(i);
      double ph3b= phi3b->GetBinContent(i);
      if(ph1b==0 || ph2b==0 || ph3b==0) continue;
      rph1b= ph1b/ph3b;
      rph2b = ph2b/ph3b;
      myfile << setprecision(9);
    myfile << rph1b << "," << rph2b << "\n";
    }
    myfile.close();

    pythia.stat();
    mcen->Write();
    pt1b->Write();
    pt2b->Write();
    pt3b->Write();
    eta1b->Write();
    eta2b->Write();
    eta3b->Write();
    phi1b->Write();
    phi2b->Write();
    phi3b->Write();


    delete outFile;
	}
