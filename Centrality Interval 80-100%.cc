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
 
 
 TFile* outFile = new TFile("firstcentrality.root", "RECREATE");
 TH1F *mcen = new TH1F("M1","Total number of final state charged particles", 50,0,7000);
 TH1F *pt1a = new TH1F("pt1a","pt of region 1, 80-100%",100,0,50);
 TH1F *pt2a = new TH1F("pt2a","pt of region 2, 80-100%",100,0,50);
 TH1F *pt3a = new TH1F("pt3a","pt of region 3, 80-100%",100,0,50);
 TH1F *eta1a = new TH1F("eta1a","eta of region 1, 80-100%",20,-4,4);
 TH1F *eta2a = new TH1F("eta2a","eta of region 2, 80-100%",20,-4,4);
 TH1F *eta3a = new TH1F("eta3a","eta of region 3, 80-100%",20,-4,4);
 TH1F *phi1a = new TH1F("phi1a","phi of region 1, 80-100%",20,-M_PI,M_PI);
 TH1F *phi2a = new TH1F("phi2a","phi of region 2, 80-100%",20,-M_PI,M_PI);
 TH1F *phi3a = new TH1F("phi3a","phi of region 3, 80-100%",20,-M_PI,M_PI);
 
 
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
    
    if(nCharged>0 && nCharged<173.2){
    if(alpha<-0.1)
    {
        for(int i=0;i<pythia.event.size();i++)
        {
          if (pythia.event[i].isFinal() && pythia.event[i].isCharged() && abs(pythia.event[i].eta()) < 2.5){
          pt1a->Fill(pythia.event[i].pT());
          eta1a->Fill(pythia.event[i].eta());
          phi1a->Fill(pythia.event[i].phi());
          }
        }
    }
    if(alpha>0.1)
    {
        for(int i=0;i<pythia.event.size();i++)
        {
          if (pythia.event[i].isFinal() && pythia.event[i].isCharged() && abs(pythia.event[i].eta()) < 2.5){
          pt2a->Fill(pythia.event[i].pT());
          eta2a->Fill(pythia.event[i].eta());
          phi2a->Fill(pythia.event[i].phi());
          }
        }
    }
    if(alpha>-0.1 && alpha<0.1)
    {
        for(int i=0;i<pythia.event.size();i++)
        {
          if (pythia.event[i].isFinal() && pythia.event[i].isCharged() && abs(pythia.event[i].eta()) < 2.5){
          pt3a->Fill(pythia.event[i].pT());
          eta3a->Fill(pythia.event[i].eta());
          phi3a->Fill(pythia.event[i].phi());
          }
        }

    }
    }

  }
  vector<double> reta1a, reta2a, rpt1a, rpt2a, rphi1a, rphi2a;
  
  double ret1a, ret2a, rp1a, rp2a, rph1a, rph2a;

  std::ofstream myfile;
  myfile.open ("pseudorapidityData.csv");
  //cout << "eta bins, 80-100%" << endl;
    for(int i=1;i<=20;i++){
      double et1a=eta1a->GetBinContent(i);
      double et2a=eta2a->GetBinContent(i);
      double et3a=eta3a->GetBinContent(i);
      if(et1a==0 || et2a==0 || et3a==0) continue;
      ret1a = et1a/et3a;
      ret2a = et2a/et3a;
    myfile << setprecision(9);
    myfile << ret1a << "," << ret2a << "\n";
    }
    myfile.close();
  
  myfile.open ("pTData.csv");
  //cout << "pt bins, 80-100%" << endl;
    for(int i=1;i<=100;i++){
      double p1a= pt1a->GetBinContent(i);
      double p2a= pt2a->GetBinContent(i);
      double p3a= pt3a->GetBinContent(i);
      if(p1a==0 || p2a==0 || p3a==0) continue;
      rp1a= p1a/p3a;
      rp2a = p2a/p3a;
    //   rpt1a.push_back(rp1a);
    //   rpt2a.push_back(rp2a);
    myfile << setprecision(9);
    myfile << rp1a << "," << rp2a << "\n";
      }
    myfile.close();
    
    myfile.open ("phiData.csv");
    //cout << "phi bins, 80-100%" << endl;
    for(int i=1;i<=20;i++){
      double ph1a= phi1a->GetBinContent(i);
      double ph2a= phi2a->GetBinContent(i);
      double ph3a= phi3a->GetBinContent(i);
      if(ph1a==0 || ph2a==0 || ph3a==0) continue;
      rph1a= ph1a/ph3a;
      rph2a = ph2a/ph3a;
    //   rphi1a.push_back(rph1a);
    //   rphi2a.push_back(rph2a);
    myfile << setprecision(9);
    myfile << rph1a << "," << rph2a << "\n";
    }
    myfile.close();

    // std::ofstream myfile;
    //   myfile.open ("LongAssymData.csv");
    //   for (int i=0; i< 100; i++)
    //   {
    //     myfile << reta1a[i] << "," << reta2a[i] << "," << rpt1a[i] << "," << rpt2a[i] << "," << rphi1a[i] << "," << rphi2a[i] << "\n";
    //   }
      
      
    //   myfile.close();
    
    pythia.stat();
    mcen->Write();
    pt1a->Write();
    pt2a->Write();
    pt3a->Write();
    eta1a->Write();
    eta2a->Write();
    eta3a->Write();
    phi1a->Write();
    phi2a->Write();
    phi3a->Write();

    delete outFile;
	}
