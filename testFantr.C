#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1F.h>
#include <TSystem.h>
#include <TMath.h>
#include <TParticle.h>
#include <TRandom.h>
#include <TROOT.h>
#include "AliESDVertex.h"
#include "AliExternalTrackParam.h"
#include "FT2.h"
#include <TTree.h>
#include <TFile.h>
#include <TParticlePDG.h>

FT2* det=0;
#endif

TH1F* hdca=0,*hdcaN=0;
TH1F* hdcaZ=0,*hdcaZN=0;
TH1F* hFdca=0,*hFdcaN=0;
TH1F* hFdcaZ=0,*hFdcaZN=0;
TH1F* hPatternITS=0;
TH1F* hFPatternITS=0;


typedef struct {
  Int_t pdg;
  Int_t qMC;
  Float_t ptMC;
  Float_t etaMC;
  Float_t phiMC;
  Float_t drMC;
  Float_t dzMC;
  //
  Int_t nrec;
  //
  Int_t q;
  Float_t zvSPD;
  Float_t zvTrc;
  Float_t pt;
  Float_t eta;
  Float_t phi;
  Float_t dca[2];
  Float_t dcaErr[2];
  Float_t sigpti;
  Float_t tpcChi2;
  Float_t trdChi2;
  //
  UChar_t nclITS;
  UChar_t nclITSF;
  UChar_t nclTPC;
  UChar_t nclTRD;
  UChar_t clITS;
  UChar_t clITSF;
  UChar_t clTRD;
  //
  Int_t mlt05;
  //
  Int_t lbl;
  Int_t lblTPC;
  Int_t lblITS;
  //
  ULong64_t   fFlags;
  ULong64_t   fMyBits;
  //
} trSumm;

trSumm tSum;

void testFantr(const char* inpTreeName="trsf.root",int ntrials=-1, double dndy=2100., Bool_t useKalmanOut=kTRUE);
void testFantr(TTree* inpTree ,int ntrials=-1, double dndy=2100.*3, Bool_t useKalmanOut=kTRUE);

void SetInpTree(TTree* inpTree);
void BookTree(const char* treeFile, const char* treeName, const char* treeTitle="SummaryTree");
void CloseTree();
TFile* treeOutFile = 0;
TTree* treeOut = 0;

void testFantr(const char* inpTreeName ,int ntrials, double dndy, Bool_t useKalmanOut)
{
  TFile* ff = TFile::Open(inpTreeName);
  TTree *tin = (TTree*)ff->Get("sum");
  testFantr(tin,ntrials,dndy,useKalmanOut);
}

void testFantr(TTree* inpTree ,int ntrials, double dndy, Bool_t useKalmanOut)
{
#if defined(__CINT__) && !defined(__MAKECINT__)
  gSystem->Load("libITSUpgradeBase.so");
  gSystem->Load("libITSUpgradeSim.so");
  gSystem->Load("libITSUpgradeRec.so");
  gROOT->ProcessLine(".L FT2.cxx+");
  FT2* det=0;
#endif
  //
  det = new FT2();
  det->InitEnvLocal();
  det->InitDetector();
  det->SetSimMat(kTRUE);
  det->SetMaxStepTGeo(1.);
  det->SetdNdY(dndy);
  det->SetUseKalmanOut(useKalmanOut);
  //
  SetInpTree(inpTree);
  BookTree("ft2anres.root","tsumft2");
  int nent = inpTree->GetEntries();
  //
  AliESDVertex *vtx = new AliESDVertex();
  double vcov[6] = {1e-4, 0, 1e-4, 0, 0, 1e-4};
  vtx->SetCovarianceMatrix(vcov);
  vtx->Print();
  //
  TParticle prt;
  int ntrAcc = 0;
  for (int ntr=0;ntr<nent;ntr++) {
    inpTree->GetEntry(ntr);
    if (tSum.qMC==0 || tSum.nrec==0) continue;
    //
    vtx->SetXv(gRandom->Gaus(-0.1, 50e-4));
    vtx->SetYv(gRandom->Gaus(0.2,  50e-4));
    vtx->SetZv(tSum.zvTrc);
    //
    double phi = tSum.phiMC;
    double eta = tSum.etaMC;
    double theta = 2*TMath::ATan(TMath::Exp(-eta));
    double pt = tSum.ptMC;
    double p = pt/TMath::Sin(theta);
    double pz = p*TMath::Cos(theta);
    //    double pt = p*TMath::Sin(theta);
    double pxyz[3]={pt*TMath::Cos(phi),pt*TMath::Sin(phi),pz};
    prt.SetPdgCode(tSum.pdg);
    TParticlePDG* ppdg = prt.GetPDG();
    if (!ppdg) continue;
    double m = ppdg->Mass();
    double en = TMath::Sqrt(p*p+m*m);
    //
    prt.SetMomentum(pxyz[0],pxyz[1],pxyz[2],en);
    prt.SetProductionVertex(vtx->GetX(),vtx->GetY(),vtx->GetZ(),0);
    //
    if (det->ProcessTrack(&prt,vtx)) {
      //const AliExternalTrackParam& prob = det->GetProbe();
      //      printf("%d %d\n",det->GetNClITS(),det->GetNClTPC());
      //    prob.Print();
      //      /*
      const double* dca = det->GetDCA();
      const double* cov = det->GetDCACov();
      tSum.dca[0] = dca[0];
      tSum.dca[1] = dca[1]; 
      tSum.dcaErr[0] = TMath::Sqrt(cov[0]);
      tSum.dcaErr[1] = TMath::Sqrt(cov[2]);
      //
      tSum.nclITSF = det->GetNClITSFakes();
      tSum.nclITS  = det->GetNClITS();
      //
      tSum.clITS = (UChar_t) det->GetITSPattern();
      tSum.clITSF = det->GetITSPatternFakes();
      //
      const FTProbe& probe = det->GetProbe();
      tSum.pt = probe.Pt();
      tSum.eta = probe.Eta();
      tSum.phi = probe.Phi();
      tSum.sigpti = TMath::Sqrt(probe.GetSigma1Pt2());
      tSum.tpcChi2 = det->GetChi2TPC();
      tSum.trdChi2 = det->GetChi2ITS();
      tSum.nclTPC =  det->GetNClTPC();
      ntrAcc++;
      treeOut->Fill();
      if ( (ntrAcc%5000)==0 ) printf("%d tracks accepted from %d read, total %d\n",ntrAcc,ntr,nent);
    }
    else {
      printf("Failed on track %d\n",ntr);
    }
    if (ntrials>0 && ntrAcc>ntrials) break;
  }
  //
  CloseTree();
}


void SetInpTree(TTree* inpTree)
{
  inpTree->SetBranchAddress("pdg", &tSum.pdg);
  inpTree->SetBranchAddress("qMC", &tSum.qMC);
  inpTree->SetBranchAddress("ptMC", &tSum.ptMC);
  inpTree->SetBranchAddress("etaMC", &tSum.etaMC);
  inpTree->SetBranchAddress("phiMC", &tSum.phiMC);
  inpTree->SetBranchAddress("drMC", &tSum.drMC);
  inpTree->SetBranchAddress("dzMC", &tSum.dzMC);
  //
  inpTree->SetBranchAddress("nrec", &tSum.nrec);
  //
  inpTree->SetBranchAddress("q", &tSum.q);
  inpTree->SetBranchAddress("zvSPD", &tSum.zvSPD);
  inpTree->SetBranchAddress("zvTrc", &tSum.zvTrc);
  inpTree->SetBranchAddress("pt", &tSum.pt);
  inpTree->SetBranchAddress("eta", &tSum.eta);
  inpTree->SetBranchAddress("phi", &tSum.phi);
  inpTree->SetBranchAddress("dca", &tSum.dca);
  inpTree->SetBranchAddress("dcaErr", &tSum.dcaErr);
  inpTree->SetBranchAddress("sigpti", &tSum.sigpti);
  inpTree->SetBranchAddress("tpcChi2", &tSum.tpcChi2);
  inpTree->SetBranchAddress("trdChi2", &tSum.trdChi2);
  //
  inpTree->SetBranchAddress("nclITS", &tSum.nclITS);
  inpTree->SetBranchAddress("nclITSF", &tSum.nclITSF);
  inpTree->SetBranchAddress("nclTPC", &tSum.nclTPC);
  inpTree->SetBranchAddress("nclTRD", &tSum.nclTRD);
  //
  inpTree->SetBranchAddress("clITS", &tSum.clITS);
  inpTree->SetBranchAddress("clITSF", &tSum.clITSF);
  inpTree->SetBranchAddress("clTRD", &tSum.clTRD);
  //
  inpTree->SetBranchAddress("mlt05", &tSum.mlt05);
  //
  inpTree->SetBranchAddress("lbl", &tSum.lbl);
  inpTree->SetBranchAddress("lblTPC", &tSum.lblTPC);
  inpTree->SetBranchAddress("lblITS", &tSum.lblITS);
  //
  inpTree->SetBranchAddress("fFlags", &tSum.fFlags);
  inpTree->SetBranchAddress("fMyBits",&tSum.fMyBits);  
}


//_______________________________________________
void BookTree(const char* treeFile, const char* treeName, const char* treeTitle)
{
  treeOutFile = TFile::Open(treeFile,"recreate");
  treeOut = new TTree(treeName,treeTitle);
  //
  treeOut->Branch("pdg", &tSum.pdg ,"pdg/I");
  treeOut->Branch("qMC", &tSum.qMC ,"qMC/I");
  treeOut->Branch("ptMC", &tSum.ptMC ,"ptMC/F");
  treeOut->Branch("etaMC", &tSum.etaMC ,"etaMC/F");
  treeOut->Branch("phiMC", &tSum.phiMC ,"phiMC/F");
  treeOut->Branch("drMC", &tSum.drMC ,"drMC/F");
  treeOut->Branch("dzMC", &tSum.dzMC ,"dzMC/F");
  //
  treeOut->Branch("nrec", &tSum.nrec ,"nrec/I");
  //
  treeOut->Branch("q", &tSum.q ,"q/I");
  treeOut->Branch("zvSPD", &tSum.zvSPD ,"zvSPD/F");
  treeOut->Branch("zvTrc", &tSum.zvTrc ,"zvTrc/F");
  treeOut->Branch("pt", &tSum.pt ,"pt/F");
  treeOut->Branch("eta", &tSum.eta ,"eta/F");
  treeOut->Branch("phi", &tSum.phi ,"phi/F");
  treeOut->Branch("dca", &tSum.dca ,"dca[2]/F");
  treeOut->Branch("dcaErr", &tSum.dcaErr ,"dcaErr[2]/F");
  treeOut->Branch("sigpti", &tSum.sigpti ,"sigpti/F");
  treeOut->Branch("tpcChi2", &tSum.tpcChi2 ,"tpcChi2/F");
  treeOut->Branch("trdChi2", &tSum.trdChi2 ,"trdChi2/F");
  //
  treeOut->Branch("nclITS", &tSum.nclITS ,"nclITS/b");
  treeOut->Branch("nclITSF", &tSum.nclITSF ,"nclITSF/b");
  treeOut->Branch("nclTPC", &tSum.nclTPC ,"nclTPC/b");
  treeOut->Branch("nclTRD", &tSum.nclTRD ,"nclTRD/b");
  //
  treeOut->Branch("clITS", &tSum.clITS ,"clITS/b");
  treeOut->Branch("clITSF", &tSum.clITSF ,"clITSF/b");
  treeOut->Branch("clTRD", &tSum.clTRD ,"clTRD/b");
  //
  treeOut->Branch("mlt05", &tSum.mlt05 ,"mlt05/I");
  //
  treeOut->Branch("lbl", &tSum.lbl ,"lbl/I");
  treeOut->Branch("lblTPC", &tSum.lblTPC ,"lblTPC/I");  
  treeOut->Branch("lblITS", &tSum.lblITS ,"lblITS/I");
  //
  treeOut->Branch("fFlags", &tSum.fFlags ,"fFlags/l");
  treeOut->Branch("fMyBits",&tSum.fMyBits ,"fMyBits/l");
}

//_______________________________________________
void CloseTree()
{
  treeOutFile->cd();
  treeOut->Write();
  delete treeOut;
  treeOutFile->Close();
  delete treeOutFile;
}
