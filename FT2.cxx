#include "FT2.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliGRPManager.h"
#include "AliLog.h"
#include "AliITSUReconstructor.h"
#include "AliITSURecoDet.h"
#include <TGeoGlobalMagField.h>
#include "AliMagF.h"
#include "AliGeomManager.h"
#include "AliITSUGeomTGeo.h"
#include "AliTrackerBase.h"
#include "AliITSUTrackerGlo.h"
#include "AliITSURecoSens.h"
#include "AliITSURecoLayer.h"
#include "AliITSsegmentation.h"
#include "AliESDVertex.h"
#include <TParticle.h>
#include <TDatabasePDG.h>
#include <TRandom.h>
#include <../TPC/Base/AliTPCcalibDB.h>
#include <../TPC/Base/AliTPCParam.h>
#include "AliPIDResponse.h"
#include "AliDetectorPID.h"


using namespace AliITSUAux;


ClassImp(FT2)

float FT2::fgMaxStepTGeo = 1.0;

//________________________________________________
FT2::FT2() :
fITSRec(0)
  ,fITS(0)
  ,fUsePIDForTracking(0)
  ,fIsTPC(kFALSE)
  ,fPIDResponse(0)
  ,fTPCSignalElectron(0)
  ,fTPCSignalMuon(0)
  ,fTPCSignalPion(0)
  ,fTPCSignalKaon(0)
  ,fTPCSignalProton(0)
  ,fTPCSectorEdge(0)
  ,fMaxSnpTPC(0.95)
  ,fTPCLayers()
  ,fTPCHitLr()
  ,fNTPCHits(0)
  ,fMinTPCHits(60)
  ,fMinITSLrHit(4)
  ,fProbe()
  ,fKalmanOutward(0)
  ,fUseKalmanOut(kTRUE)
  //,fProbeMass(0.14)
  ,fBz(0.5)
  ,fSimMat(kFALSE)
  ,fCurrITSLr(-1)
  ,fNClTPC(0)
  ,fNClITS(0)
  ,fNClITSFakes(0)
  ,fITSPatternFake(0)  
  ,fITSPattern(0)
  ,fChi2TPC(0)
  ,fChi2ITS(0)
  ,fSigYITS(3.14e-4) // 5e-4
  ,fSigZITS(3.38e-4) // 5e-4
  ,fNITSLrHit(0)
  ,fdNdY(-1)
{
  //
  AliInfo("Default Constructor");

  memset(fITSSensHit,0,kMaxITSLr*2*sizeof(AliITSURecoSens*));
}

//________________________________________________
FT2::~FT2()
{
  //destroy
  AliInfo("Destroy");
  delete[] fKalmanOutward;
  delete fITSRec;
}

//________________________________________________
void FT2::InitTPCSignalFiles(const char *TPCsignalElectrons,const char *TPCsignalMuons,const char *TPCsignalPions,const char *TPCsignalKaons,const char *TPCsignalProtons)
{
  AliInfo("Setting Files");
	
  fTPCSignalElectron = TFile::Open(TPCsignalElectrons);
  if(fTPCSignalElectron->IsZombie()) AliFatal("Problem with opening TPC signal file for electrons - File not available!");
	
  fTPCSignalMuon = TFile::Open(TPCsignalMuons);
  if(fTPCSignalMuon->IsZombie()) AliFatal("Problem with opening TPC signal file for muons - File not available!");
	
  fTPCSignalPion = TFile::Open(TPCsignalPions);
  if(fTPCSignalPion->IsZombie()) AliFatal("Problem with opening TPC signal file for pions - File not available!");
	
  fTPCSignalKaon = TFile::Open(TPCsignalKaons);
  if(fTPCSignalKaon->IsZombie()) AliFatal("Problem with opening TPC signal file for kaons - File not available!");
	
  fTPCSignalProton = TFile::Open(TPCsignalProtons);
  if(fTPCSignalProton->IsZombie()) AliFatal("Problem with opening TPC signal file for protons - File not available!");
}
//________________________________________________
void FT2::InitDetector(Bool_t addTPC, Float_t sigYTPC,Float_t sigZTPC,Float_t effTPC,Float_t scEdge)
{
  AliInfo("Initializing Detector");

  // read full geometry and initialize qauzi-layers for TPC if requested
  AliCDBManager* man = AliCDBManager::Instance();
  if (!man->IsDefaultStorageSet() || man->GetRun()<0) AliFatal("CDB manager is not configured");
  //
  if ( !TGeoGlobalMagField::Instance()->GetField() ) {
    AliWarning("Magnetic field is not initialized, loading from GRP");
    AliGRPManager grpMan;
    if(!grpMan.ReadGRPEntry()) AliFatal("Cannot get GRP entry");
    if(!grpMan.SetMagField()) AliFatal("Problem with magnetic field setup");
  }
	
  fBz = AliTrackerBase::GetBz();
  //
  AliGeomManager::LoadGeometry("geometry.root");
  AliCDBEntry* ent = man->Get("ITS/Calib/RecoParam");
  AliITSURecoParam* par = (AliITSURecoParam*)((TObjArray*)ent->GetObject())->At(1); // need just to initialize the interface
  fITSRec = new AliITSUReconstructor();
  fITSRec->SetRecoParam(par);
  fITSRec->Init();
  //
  fITS = fITSRec->GetITSInterface();
  if (addTPC) AddTPC(sigYTPC,sigZTPC,effTPC,scEdge);
  int nLrITS = fITS->GetNLayersActive();
  fKalmanOutward = new AliExternalTrackParam[nLrITS];
  //
  AliTPCcalibDB * calib = AliTPCcalibDB::Instance();
  const AliMagF * field = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  calib->SetExBField(field);
  calib->SetRun(138871);
	
  //
}
//____________________________________________________
void FT2::InitTPCPIDResponse()
{
  AliInfo("Initializing TPC PID Response");
	
  fPIDResponse = new AliPIDResponse(kTRUE);
  fPIDResponse->GetTPCResponse().SetUseDatabase(kFALSE);
  AliTPCParam* param = AliTPCcalibDB::Instance()->GetParameters();
  if (param){
    TVectorD *paramBB = param->GetBetheBlochParameters();
    if (paramBB){
      fPIDResponse->GetTPCResponse().SetBetheBlochParameters((*paramBB)(0),(*paramBB)(1),(*paramBB)(2),(*paramBB)(3),(*paramBB)(4));
      AliInfo(Form("Setting BB parameters from OCDB (AliTPCParam): %.2g, %.2g, %.2g, %.2g, %.2g",(*paramBB)(0),(*paramBB)(1),(*paramBB)(2),(*paramBB)(3),(*paramBB)(4)));
    }
    else {
      AliError("Couldn't get BB parameters from OCDB, the old default values will be used instead");
    }
  }
  else {
    AliError("Couldn't get TPC parameters");
  }
}
//____________________________________________________
void FT2::InitEnvLocal()
{
  AliInfo("Initializing Local Environment");

  // init with environment of local ITSU generation
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  man->SetSpecificStorage("GRP/*/*","alien://folder=/alice/data/2010/OCDB");
  man->SetSpecificStorage("ITS/*/*", "alien://folder=/alice/simulation/LS1_upgrade/Ideal");
  man->SetSpecificStorage("TPC/*/*", "alien://folder=/alice/simulation/2008/v4-15-Release/Ideal/");
  //	man->SetSpecificStorage("GRP/GRP/Data","alien://folder=/alice/data/2010/OCDB");
  //	man->SetSpecificStorage("ITS/Align/Data", "alien://folder=/alice/simulation/LS1_upgrade/Ideal");
  //	man->SetSpecificStorage("ITS/Calib/RecoParam", "alien://folder=/alice/simulation/LS1_upgrade/Ideal");
  man->SetRun(138871);
  man->Print("");
  /*
    AliCDBManager* man = AliCDBManager::Instance();
    man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    man->SetSpecificStorage("GRP/GRP/Data",
    Form("local://%s",gSystem->pwd()));
    man->SetSpecificStorage("ITS/Align/Data",
    Form("local://%s",gSystem->pwd()));
    man->SetSpecificStorage("ITS/Calib/RecoParam",
    Form("local://%s",gSystem->pwd()));
    man->SetRun(0);*/
  //
}

//____________________________________________________
void FT2::AddTPC(Float_t sigY, Float_t sigZ, Float_t eff,Float_t scEdge)
{
  // add TPC mock-up
  //
  const Float_t kRadLPerRow = 0.000036;
  //
  const Float_t kTPCInnerRadialPitch  =    0.75 ;    // cm
  const Float_t kTPCMiddleRadialPitch =    1.0  ;    // cm
  const Float_t kTPCOuterRadialPitch  =    1.5  ;    // cm
  const Int_t   kTPCInnerRows            =   63 ;
  const Int_t   kTPCMiddleRows           =   64  ;
  const Int_t   kTPCOuterRows            =   32  ;
  const Int_t   kTPCRows       =   (kTPCInnerRows + kTPCMiddleRows + kTPCOuterRows) ;
  const Float_t kTPCRowOneRadius         =   85.2  ;    // cm
  const Float_t kTPCRow64Radius          =  135.1  ;    // cm
  const Float_t kTPCRow128Radius         =  199.2  ;    // cm
  //
  if (fIsTPC) {
    printf("TPC was already added\n");
    return;
  }
  fIsTPC = kTRUE;
  fTPCSectorEdge = scEdge;
  //
  for ( Int_t k = 0 ; k < kTPCRows ; k++ ) {
    //
    Float_t rowRadius =0;
    if (k<kTPCInnerRows) rowRadius =  kTPCRowOneRadius + k*kTPCInnerRadialPitch ;
    else if ( k>=kTPCInnerRows && k<(kTPCInnerRows+kTPCMiddleRows) )
      rowRadius =  kTPCRow64Radius + (k-kTPCInnerRows+1)*kTPCMiddleRadialPitch ;
    else if (k>=(kTPCInnerRows+kTPCMiddleRows) && k<kTPCRows )
      rowRadius = kTPCRow128Radius + (k-kTPCInnerRows-kTPCMiddleRows+1)*kTPCOuterRadialPitch ;
    AddTPCLayer(rowRadius,kRadLPerRow,sigY,sigZ,eff);
  }
  //
  fTPCHitLr.resize(fTPCLayers.size());
}

//____________________________________________________
void FT2::AddTPCLayer(Float_t x, Float_t x2x0,Float_t sigY, Float_t sigZ, Float_t eff)
{
  // add single TPC layer
  fTPCLayers.push_back(FT2TPCLayer(x,x2x0,sigY,sigZ,eff));
}

//____________________________________________________
void FT2::PrintLayout()
{
  // print setup
  if (fITS) fITS->Print("lr");
  if (fIsTPC) {
    printf("TPC inactive sector edge: %.3f\n",fTPCSectorEdge);
    printf("TPC \t  R   \tx2x0\tsgRPhi\t sigZ  \t Eff  \n");
    for (int ilr=0;ilr<(int)fTPCLayers.size();ilr++) {
      FT2TPCLayer& lr = fTPCLayers[ilr];
      printf("%3d\t%6.2f\t%1.4f\t",ilr,lr.x,lr.x2x0);
      //
      if (lr.isDead) printf("      \t");
      else printf("%6.4f\t",lr.rphiRes);
      //
      if (lr.isDead) printf("      \t");
      else printf("%6.4f\t",lr.zRes);
      //
      if (lr.isDead) printf("      \t");
      else printf("%6.4f\t",lr.eff);
      //
      printf("\n");
    }
  }
}

//____________________________________________________
Bool_t FT2::ProcessTrack(TParticle* part, AliVertex* vtx)
{
  // track the particle throught the setup; it must be provided in the production point
  // then relate to vtx
  //
  static AliESDVertex vtx0;
  static Bool_t first = kTRUE;
  if (first) {
    double vd[6] = {0};
    vtx0.SetXYZ(&vd[0]);
    vtx0.SetCovarianceMatrix(&vd[0]);
    first = kFALSE;
  }
  if (!InitProbe(part)) {
#if DEBUG
    printf("Initialization failed\n");
    part->Print();
    fProbe.Print();
#endif
    return kFALSE;
  }
  PrepareProbe();
  //
  // check if there is something to reconstruct
  if ( (fIsTPC && fNTPCHits<fMinTPCHits) || (fNITSLrHit<fMinITSLrHit)) {
#if DEBUG
    printf("Track hit requirements are not satisfied: Hit ITS Layers: %d (need %d) TPCHits: %d (need %d)\n",
	   fNITSLrHit,fMinITSLrHit, fNTPCHits, fIsTPC ? fMinTPCHits:0);
    fProbe.Print();
#endif
    return kFALSE;
  }
	
  ResetCovMat(&fProbe);
  //
  if (!ReconstructProbe()) {    
#if DEBUG
    printf("Track reconstruction failed\n");
    fProbe.Print();
#endif
    return kFALSE;
  }
  //
  // go to innermost radius of ITS (including beam pipe)
  if (!PropagateToR(fITS->GetRMin(),-1, kTRUE, kFALSE, kTRUE)) {
#if DEBUG
    printf("Track propagatation to ITS RMin=%f failed\n",fITS->GetRMin());
    fProbe.Print();
#endif  
    return kFALSE; // don't exit on faile propagation, may simply not reach this point
  }
  //
  AliVertex* vtuse = vtx ? vtx : (AliVertex*)&vtx0; // if no vertex provided, relate to 0,0,0
  //
  fDCACov[0] = fDCACov[1] = fDCACov[2] = 0;
  Bool_t res = fProbe.PropagateToDCA(vtuse,fBz, fITS->GetRMin(), fDCA, fDCACov);

#if DEBUG
  if (!res) {
    printf("Track propagation to vertex failed\n");
    vtuse->Print();
    fProbe.Print();
  }
#endif  
  //
  /*
    printf("Tracking done: NclITS: %d (chi2ITS=%.3f) NclTPC: %3d (chi2TPC=%.3f)\n",fNClITS,fChi2ITS,fNClTPC,fChi2TPC);
    fProbe.Print();
    //
    if (res) {
    printf("DCA: {%e %e} DCACov: {%e %e %e}\n",fDCA[0],fDCA[1],fDCACov[0],fDCACov[1],fDCACov[2]);
    }
    else {
    printf("Failed to propagate to vertex\n");
    }
  */
  //
  return res;
}

//____________________________________________________
Bool_t FT2::InitProbe(TParticle* part)
{
  // create probe (code copied from AliESDtrack::AliESDtrack(TParticle * part) :
  //
  Double_t xref;
  Double_t alpha;
  Double_t param[5];
  Double_t covar[15] = {0};
  //
  fNTPCHits = fNClTPC = fNClITS = fNClITSFakes = fNITSLrHit = fITSPatternFake = fITSPattern = 0;
  fChi2TPC = fChi2ITS = 0;
  fDCA[0] = fDCA[1] = fDCACov[0] = fDCACov[1] = fDCACov[2] = 0.;
  //
  // Calculate alpha: the rotation angle of the corresponding local system (TPC sector)
  alpha = part->Phi()*180./TMath::Pi();
  if (alpha<0) alpha+= 360.;
  if (alpha>360) alpha -= 360.;
  //
  Int_t sector = (Int_t)(alpha/20.);
  alpha = 10. + 20.*sector;
  alpha /= 180;
  alpha *= TMath::Pi();
  //
  // Get the vertex of origin and the momentum
  TVector3 ver(part->Vx(),part->Vy(),part->Vz());
  TVector3 mom(part->Px(),part->Py(),part->Pz());
	
  // Rotate to the local coordinate system (TPC sector)
  ver.RotateZ(-alpha);
  mom.RotateZ(-alpha);
	
  // X of the referense plane
  xref = ver.X();
	
  Int_t pdgCode = part->GetPdgCode();
  TParticlePDG* pdgp = TDatabasePDG::Instance()->GetParticle(pdgCode);
  Double_t charge = pdgp->Charge();
  fProbe.fProbeMass = pdgp->Mass();
  fProbe.fAbsPdgCode = TMath::Abs(pdgCode);
  //
  param[0] = ver.Y();
  param[1] = ver.Z();
  param[2] = TMath::Sin(mom.Phi());
  param[3] = mom.Pz()/mom.Pt();
  param[4] = TMath::Sign(1/mom.Pt(),charge);
  //	AliInfo(Form("Probe Mass is: %f ",fProbe.fProbeMass));
  // Set AliExternalTrackParam
  fProbe.Set(xref, alpha, param, covar);
  ResetCovMat(&fProbe);
  return kTRUE;
}

//____________________________________________________
Bool_t FT2::PrepareProbe()
{
  // propagate the probe to max allowed R of the setup
  int nlrITS = fITS->GetNLayers(); // including passive layers
  //
  double xyzTmp[3];
  //
  for (int ilr=0;ilr<nlrITS;ilr++) {
    AliITSURecoLayer* lr = fITS->GetLayer(ilr);
    if (lr->IsPassive()) { // cylindric passive layer, just go to this layer
      if (!PropagateToR(lr->GetRMax(),1,kFALSE,fSimMat,fSimMat)) return kFALSE;
    }
    else { // active layer, need to simulate the hit positions
      if (!PassActiveITSLayer(lr)) return kFALSE;
    }
  }
  //
  fNTPCHits = 0;
  if (fIsTPC) {
    const double kTanSectH = 1.76326980708464975e-01; // tangent of  10 degree (sector half opening)
    double xyz[3];
    fProbe.GetXYZ(xyz);
    Double_t alphan = TMath::ATan2(xyz[1], xyz[0]);
    if (!fProbe.RotateParamOnly(alphan)) {
      return kFALSE; // start in the frame with X pointing to track position
    }
    if (!PropagateToR(fTPCLayers[0].x-0.1,1,kFALSE,fSimMat,fSimMat)) {
      return kFALSE; // reach inner layer
    }
    //
    int sector = -1;
    for (int ilr=0;ilr<(int)fTPCLayers.size();ilr++) {
#if DEBUG>5
      printf("At TPC lr %d\n",ilr);
      fProbe.Print();
#endif
      FT2TPCLayer_t &tpcLr = fTPCLayers[ilr];
      tpcLr.hitSect = -1;
      fProbe.GetXYZ(xyzTmp);
      double phi = TMath::ATan2(xyzTmp[1],xyzTmp[0]);
      BringTo02Pi(phi);
      Int_t sectNew = (Int_t)(phi*TMath::RadToDeg()/20.);
      if (sector!=sectNew) {
	sector = sectNew;
	phi = 10. + 20.*sector;
	phi /= 180;
	phi *= TMath::Pi();
	if (!fProbe.RotateParamOnly(phi)) {
#if DEBUG
	  printf("TPC:Failed to rotate to alpha %.4f (sc %d)\n",phi,sector);
	  fProbe.Print();
#endif
	  return kFALSE; // go to sector frame
	}
      }
      if (!fProbe.PropagateParamOnlyTo(tpcLr.x, fBz)) {
#if DEBUG
	printf("TPC:Failed to go to X: %.4f\n",tpcLr.x);
	fProbe.Print();
#endif
	return kFALSE; // no materials inside TPC
      }
      if (TMath::Abs(fProbe.GetZ())>kMaxZTPC) break; // exit from the TPC
      double maxY = kTanSectH*tpcLr.x - fTPCSectorEdge; // max allowed Y in the sector
      if (TMath::Abs(fProbe.GetY())>maxY) {
#if DEBUG>3
	printf("TPC: No hit in dead zone: %f\n",maxY);
	fProbe.Print();
#endif
	continue; // in dead zone
      }
      //
      // if track Snp > limit, stop propagation
      if (TMath::Abs(fProbe.GetSnp())>fMaxSnpTPC) {
#if DEBUG>2
	printf("TPC: Max snp %f reached\n",fMaxSnpTPC);
	fProbe.Print();
#endif
	return kFALSE;
      }
			
      // register hit
#if DEBUG>5
      printf("Register TPC hit at sect:%d, Y:%.4f Z:%.4f\n",sector,fProbe.GetY(),fProbe.GetZ());
      fProbe.Print();
#endif
      if (tpcLr.isDead || (tpcLr.eff<1 && gRandom->Rndm()>tpcLr.eff)) continue;
      //
      tpcLr.hitSect = sector;
      tpcLr.hitY = fProbe.GetY();
      tpcLr.hitZ = fProbe.GetZ();
      fTPCHitLr[fNTPCHits++] = ilr;
      //
    }
  }
  //
  return kTRUE;
}


//________________________________________________________________________
Bool_t FT2::ReconstructProbe()
{
  // reconstruct the probe
	
  //
  int sect = -1;
  fChi2TPC = 0;
  fChi2ITS = 0;
  fNClITS = 0;
  fNClITSFakes = 0;
  fITSPatternFake = 0;
  fITSPattern = 0;
  fNClTPC = 0;
  fCurrITSLr = -1;	
  //
  if (fNTPCHits) {
    // TPC tracking
#if DEBUG>5
    AliInfo(Form("Number of clusters generated on the pad rows: %i",fNTPCHits));
#endif
    for (int ih=fNTPCHits;ih--;) {
#if DEBUG>5
      AliInfo(Form("Entering TPC cluster loop: %i",ih));
#endif
      if(ih==0 && fUsePIDForTracking){
	double signalTPC, p = fProbe.P(); fProbe.fTPCmomentum = p;

	TH1* hdedx;
	AliPID::EParticleType typePID=AliPID::EParticleType(2);

	if(fProbe.fAbsPdgCode==11)
	  {
	    int pbin = ((TH1*)fTPCSignalElectron->Get("hMomentumAxis"))->FindBin(p);
	    hdedx = (TH1*)fTPCSignalElectron->Get(Form("electronDEDX%i",pbin));
	    typePID=AliPID::EParticleType(0);
	  }
	else if(fProbe.fAbsPdgCode==13)
	  {
	    int pbin = ((TH1*)fTPCSignalMuon->Get("hMomentumAxis"))->FindBin(p);
	    hdedx = (TH1*)fTPCSignalMuon->Get(Form("muonDEDX%i",pbin));
	    typePID=AliPID::EParticleType(1);
	  }
	else if(fProbe.fAbsPdgCode==211)
	  {
	    int pbin = ((TH1*)fTPCSignalPion->Get("hMomentumAxis"))->FindBin(p);
	    hdedx = (TH1*)fTPCSignalPion->Get(Form("pionDEDX%i",pbin));
	    typePID=AliPID::EParticleType(2);
	  }
	else if(fProbe.fAbsPdgCode==321)
	  {
	    int pbin = ((TH1*)fTPCSignalKaon->Get("hMomentumAxis"))->FindBin(p);
	    hdedx = (TH1*)fTPCSignalKaon->Get(Form("kaonDEDX%i",pbin));
	    typePID=AliPID::EParticleType(3);
	  }
	else if(fProbe.fAbsPdgCode==2212)
	  {
	    int pbin = ((TH1*)fTPCSignalProton->Get("hMomentumAxis"))->FindBin(p);
	    hdedx = (TH1*)fTPCSignalProton->Get(Form("protonDEDX%i",pbin));
	    typePID=AliPID::EParticleType(4);
	  }
	else {hdedx=0;AliFatal("PDC Code %4.f cannot be treated in the code - This case should not be possible here!");}
#if DEBUG>5
	AliInfo(Form("PDG code of Probe: %4.f",fProbe.fAbsPdgCode));
	AliInfo(Form("Momentum of Probe: %f",fProbe.fTPCmomentum));
#endif
	if(hdedx->Integral()>0){ signalTPC = ((TH1*)hdedx)->GetRandom();} //30
	else{
					
	  Double_t mean = fPIDResponse->GetTPCResponse().GetExpectedSignal(&fProbe,typePID,AliTPCPIDResponse::kdEdxDefault,kFALSE,kFALSE);
	  fProbe.fTPCSignalN = (UShort_t)mean;
	  Double_t sigma = fPIDResponse->GetTPCResponse().GetExpectedSigma(&fProbe,typePID,AliTPCPIDResponse::kdEdxDefault,kFALSE,kFALSE);
	  signalTPC = (40./50)*gRandom->Gaus(mean,sigma); // 40. or 33. ?
#if DEBUG>5
	  AliInfo(Form("TPC Signal for %i from scaling %f",typePID,signalTPC));
#endif
	}
	fProbe.fTPCSignal = signalTPC;
	fProbe.fTPCSignalN = (UShort_t)signalTPC;

	Double_t prob[AliPID::kSPECIES]={0.};
	GetComputeTPCProbability(&fProbe,AliPID::kSPECIES,prob);
				
	Float_t max=0.,min=1.e9;
	Int_t pid=-1;
	for (Int_t i=0; i<AliPID::kSPECIES; ++i) {
	  if (prob[i]>max) {pid=i; max=prob[i];}
	  if (prob[i]<min) min=prob[i];
	}
				
	if (pid>AliPID::kSPECIES-1 || (min>=max)) pid = AliPID::kPion;
				
#if DEBUG>5
	for(Int_t i=0;i<AliPID::kSPECIES;i++){
	  AliInfo(Form("Pid Prob : %f",prob[i]));
	}
#endif
	Int_t pdgCode = 0;
	if(pid==0) pdgCode = 11;
	if(pid==1) pdgCode = 13;
	if(pid==2) pdgCode = 211;
	if(pid==3) pdgCode = 321;
	if(pid==4) pdgCode = 2212;
	if(pdgCode==0) AliFatal("This should not have happened");
#if DEBUG>5
	AliInfo(Form("Probe mass for tracking: %f",fProbe.fProbeMass));
#endif
	if(pid==0) fProbe.fAbsPdgCodeForTracking = 211;			// e(11)	--> track as pion due to bug in fMC
	if(pid==1) fProbe.fAbsPdgCodeForTracking = 211;			// mu(13)	--> track as pion due to bug in fMC
	if(pid==2) fProbe.fAbsPdgCodeForTracking = 211;			// pi(211)
	if(pid==3) fProbe.fAbsPdgCodeForTracking = 321;			// ka(321)
	if(pid==4) fProbe.fAbsPdgCodeForTracking = 2212;		// pr(2212)
	if(fProbe.fAbsPdgCodeForTracking==0) AliFatal("This should not have happened");
	TParticlePDG* pdgp = TDatabasePDG::Instance()->GetParticle(fProbe.fAbsPdgCodeForTracking);
	fProbe.fProbeMass = pdgp->Mass();
#if DEBUG>5
	AliInfo(Form("PID for tracking: %i with probe mass %f",pid,fProbe.fProbeMass));
#endif
      }
			
      FT2TPCLayer_t &tpcLr = fTPCLayers[fTPCHitLr[ih]];
      if (tpcLr.hitSect!=sect) { // rotate to hit sector
	sect = tpcLr.hitSect;
	if (!fProbe.Rotate( TMath::DegToRad()*(10.+20.*sect))) return kFALSE;
      }
      // go to the padrow
      if (!fProbe.PropagateTo(tpcLr.x,fBz)) return kFALSE;
      if (tpcLr.x2x0>0 && !fProbe.CorrectForMeanMaterial(tpcLr.x2x0,0,fProbe.fProbeMass) ) return kFALSE;
      //
      if (tpcLr.isDead) continue;
      double chi = UpdateKalman(&fProbe,tpcLr.hitY,tpcLr.hitZ,tpcLr.rphiRes,tpcLr.zRes,kTRUE,kFALSE);
      if (chi<0) return kFALSE;
      fChi2TPC += chi;
      fNClTPC++;
    }
    // go to ITS/TPC matching R, accounting for TGeo materials
    if (!PropagateToR(fITS->GetRITSTPCRef(),-1,kTRUE, kFALSE, kTRUE)) return kFALSE;
  }
  //
  // ITS tracking
  for (int ilr=fITS->GetNLayersActive();ilr--;) {
    fCurrITSLr = ilr;
#if DEBUG>5
    AliInfoF("Entering ITS cluster loop on lr%d: %d clusters",ilr,fNITSHits[ilr]);
#endif
    if (!fNITSHits[ilr]) continue;
    // for 2 hits (overlap) chose only 1
    int iht = fNITSHits[ilr]==1 ? 0 : (gRandom->Rndm()>0.5 ? 0:1);
    AliITSURecoSens* sens = fITSSensHit[ilr][iht];
    if (!fProbe.Rotate(sens->GetPhiTF())) {
#if DEBUG
      printf("Failed to rotate to sensor phi %f at lr %d\n",sens->GetPhiTF(),ilr); 
      fProbe.Print();
      return kFALSE;
#endif
    }
    if (!PropagateToX(sens->GetXTF(),-1,kTRUE, kFALSE, kTRUE)) return kFALSE; // account for materials
    double chi = UpdateKalman(&fProbe,fITSHitYZ[ilr][iht][0],fITSHitYZ[ilr][iht][1],fSigYITS,fSigZITS,kTRUE,kTRUE);
 #if DEBUG>5
    AliInfoF("Updated at ITS Lr%d, chi2:%f NclITS:%d Fakes: %d",ilr,chi,fNClITS,fNClITSFakes);
#endif   
    if (chi<0) return kFALSE;
    fChi2ITS += chi;
    fNClITS++;
    fITSPattern |= 0x1<<ilr;
  }
  //
  return kTRUE;
}

//________________________________________________________________________
Double_t FT2::UpdateKalman(AliExternalTrackParam* trc, double y,double z,double sigY,double sigZ,Bool_t randomize, Bool_t fake)
{
  // calman update, return chi2 increment or -1 on failure
  //
  if (randomize) {
    double ry,rz;
    gRandom->Rannor(ry,rz);
    y += sigY*ry;
    z += sigZ*rz;
  }
  double meas[2] = {y,z};
  double measErr2[3] = {sigY*sigY,0,sigZ*sigZ};
  //
  // do we want to simulate fakes?
  if (fake && fdNdY>0) {
    double err2a,err2b;   
    double trCov[3];
    double* covInw = (double*)trc->GetCovariance();
    if (fCurrITSLr>=0 && fKalmanOutward[fCurrITSLr].GetSigmaY2()>0) { // weight with outward calman track
      double* covOut = (double*)fKalmanOutward[fCurrITSLr].GetCovariance();
      double detInw = covInw[0]*covInw[2]-covInw[1]*covInw[1];
      double detOut = covOut[0]*covOut[2]-covOut[1]*covOut[1];
      if (detInw<1e-16 || detOut<1e-16) {
#if DEBUG
	AliInfoF("Failed Update for fakes: detInw: %e detOutw: %e",detInw,detOut);
	printf("Inward:  ");trc->Print();
	printf("Outward: ");fKalmanOutward[fCurrITSLr].Print();	
#endif
	return -1;
      }
      detInw = 1./detInw;
      detOut = 1./detOut;
      double si00 = covInw[2]*detInw+covOut[2]*detOut;
      double si11 = covInw[0]*detInw+covOut[0]*detOut;
      double si01 =-covInw[1]*detInw-covOut[1]*detOut;
      double det = si00*si11-si01*si01;
      if (det<1e-1) {
#if DEBUG
	AliInfoF("Failed Update for fakes: detFinal: %e",det);
	printf("Inward:  ");trc->Print();
	printf("Outward: ");fKalmanOutward[fCurrITSLr].Print();	
#endif
	return -1;
      }
      trCov[0] = si11/det;
      trCov[2] = si00/det;
      trCov[1] = -si01/det;
#if DEBUG//>5
      printf("Errors on Lr %d: \n",fCurrITSLr);
      printf("Inward  :  %+e %+e %+e\n",covInw[0],covInw[1],covInw[2]);
      printf("Outward :  %+e %+e %+e\n",covOut[0],covOut[1],covOut[2]);
      printf("Weighted:  %+e %+e %+e\n",trCov[0],trCov[1],trCov[2]);
#endif
    }
    else for (int i=3;i--;) trCov[i] = covInw[i];
    //
    if (!DiagonalizeErrors(trCov,err2a,err2b)) {
#if DEBUG
      printf("Failed to diagonalize erros\n"); trc->Print();
#endif
      return -1;
    }
    double xyz[3]; 
    trc->GetXYZ(xyz);
    double rho = HitDensity(xyz[0]*xyz[0]+xyz[1]*xyz[1],trc->GetTgl());
    //
    // prob. of good hit (http://rnc.lbl.gov/~wieman/HitFinding2D.htm)
    Double_t sx2, sy2, probFake;
    sx2 = 2 * TMath::Pi()*err2a*rho;
    sy2 = 2 * TMath::Pi()*err2b*rho;
    probFake = 1.-TMath::Sqrt(1./((1+sx2)*(1+sy2)));
#if DEBUG
    printf("DiagErr: %e %e Rho:%e PFake: %e\n",err2a,err2b, rho, probFake);
#endif
    if (gRandom->Rndm()<probFake) { // need to fake the hit
      BiasAsFake(meas, trc->GetParameter(), trCov);
      if (fCurrITSLr>=0) {
	fITSPatternFake |= 0x1<<fCurrITSLr;
	fNClITSFakes++;
      }
      /*
	printf("rho:%f errs: %f %f -> %f -> %d r=%f\n",rho,err2a,err2b,probFake,fNClITSFakes, 
	TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]));
      */
    }
  }
  //
  double chi2 = trc->GetPredictedChi2(meas,measErr2);
  if (fake) {
    printf("Updating {%e %e}/{%e %e} {%e %e %e} Ncl:%d NclF:%d -> %f\n", meas[0],meas[1],y,z,
	   measErr2[0],measErr2[1],measErr2[2],fNClITS,fNClITSFakes,chi2); 
    trc->Print();
  }
  if (!trc->Update(meas,measErr2)) {
#if DEBUG
    printf("Failed to Update {%e %e}/{%e %e} {%e %e %e} Ncl:%d NclF:%d\n", meas[0],meas[1],y,z,
	   measErr2[0],measErr2[1],measErr2[2],fNClITS,fNClITSFakes); 
    trc->Print();
#endif
     return -1;
  }
  return chi2;
}

//________________________________________________________________________
void FT2::ResetCovMat(AliExternalTrackParam* trc)
{
  // assign huge errors to probe at max.radius
  enum {kY,kZ,kSnp,kTgl,kPtI};              // track parameter aliases
  enum {kY2,kYZ,kZ2,kYSnp,kZSnp,kSnp2,kYTgl,kZTgl,kSnpTgl,kTgl2,kYPtI,kZPtI,kSnpPtI,kTglPtI,kPtI2}; // cov.matrix aliases
  const double kLargeErr2Coord = 5*5;
  const double kLargeErr2Dir = 0.7*0.7;
  const double kLargeErr2PtI = 30.5*30.5;
  double *trPars = (double*)trc->GetParameter();
  double *trCov  = (double*)trc->GetCovariance();
  //
  for (int ic=15;ic--;) trCov[ic] = 0.;
  trCov[kY2]   = trCov[kZ2]   = kLargeErr2Coord;
  trCov[kSnp2] = trCov[kTgl2] = kLargeErr2Dir;
  trCov[kPtI2] = kLargeErr2PtI*trPars[kPtI]*trPars[kPtI];
  trc->CheckCovariance();
  //
}


//________________________________________________________________________
Bool_t FT2::PropagateToR(double r, int dir,
			 Bool_t propErr, // propagate errors
			 Bool_t simMat,  // simulate material effects
			 Bool_t useTGeo)
{
  // propagate track to given R (not X!)
  double xTgt;
  if (!fProbe.GetXatLabR(r, xTgt, fBz, dir)) {
#if DEBUG
    printf("PropagateToR failed in GetXatLabR(%.3f,%.3f, %.3f, %d)\n", r,xTgt, fBz, dir);
    fProbe.Print();
#endif
    return kFALSE;
  }
  double dx = xTgt-fProbe.GetX();
  if (dir*dx<0) return kTRUE; // probe is already above given R
  return PropagateToX(xTgt,dir,propErr,simMat,useTGeo);
}

//________________________________________________________________________
Bool_t FT2::PropagateToX(double xTgt, int dir,
			 Bool_t propErr, // propagate errors
			 Bool_t simMat,  // simulate material effects
			 Bool_t useTGeo)
{
  // propagate track to given X
  //
  Bool_t doKalmanOut = kFALSE;
  if (simMat && dir>0 && fUseKalmanOut && xTgt<fITS->GetRMax()) {
    doKalmanOut = kTRUE;
  }
  double maxStep = useTGeo ? fgMaxStepTGeo : 1e6;
  //
  Double_t xyz0[3],xyz1[3],param[7];
  //
  if (useTGeo) {
    fProbe.GetXYZ(xyz0);
  }
  while ( dir*(xTgt-fProbe.GetX()) > kTrackToler) {
    Double_t step = dir*TMath::Min(TMath::Abs(xTgt-fProbe.GetX()), maxStep);
    Double_t x    = fProbe.GetX()+step;
    Bool_t res = (propErr||doKalmanOut) ? fProbe.PropagateTo(x,fBz) : fProbe.PropagateParamOnlyTo(x,fBz);
    if (!res)  {
#if DEBUG
      printf("Failed Prop PropagateToX(%.1f,%d,%d,%d,%d) doKalmanOut=%d",xTgt,dir,propErr,simMat,useTGeo,doKalmanOut);
      fProbe.Print();
#endif
      return kFALSE;
    }
    if (useTGeo) {
      fProbe.GetXYZ(xyz1);
      AliTrackerBase::MeanMaterialBudget(xyz0,xyz1,param);
      Double_t xrho=param[0]*param[4], xx0=param[1];
      if (dir>0) xrho = -xrho;
      if (simMat) {
	if (!ApplyMSEloss(xx0,xrho)) {
#if DEBUG
	  printf("Failed ApplyMSEloss(%f,%f) in PropagateToX(%.1f,%d,%d,%d,%d)  doKalmanOut=%d",
		 xx0,xrho,xTgt,dir,propErr,simMat,useTGeo,doKalmanOut);
	  fProbe.Print();
#endif
	  return kFALSE;
	}
	if (doKalmanOut) { // special mode for outward kalman: we don't want to affect params by materials since
	  // it is already done, onle cov.matrix need to be updated
	  AliExternalTrackParam trTmp = fProbe;
	  if (!trTmp.CorrectForMeanMaterial(xx0,xrho,fProbe.fProbeMass,kTRUE)) {
#if DEBUG
	    printf("Failed MatCorr(%f %f) in PropagateToX(%.1f,%d,%d,%d,%d) for doKalmanOut",xx0,xrho,xTgt,dir,propErr,simMat,useTGeo);
	      trTmp.Print();
#endif
	  }
	  memcpy((double*)fProbe.GetCovariance(),trTmp.GetCovariance(),15*sizeof(double)); // restore parameter
	}
      }
      else if (!fProbe.CorrectForMeanMaterial(xx0,xrho,fProbe.fProbeMass,kTRUE)) {
#if DEBUG
	printf("Failed CorrMeanMat(%f,%f) in PropagateToX(%.1f,%d,%d,%d,%d)",xx0,xrho,xTgt,dir,propErr,simMat,useTGeo);
	fProbe.Print();
#endif
	return kFALSE;
      }
    }
    memcpy(xyz0,xyz1,3*sizeof(double));
    //
  }
  //
  return kTRUE;
}

//______________________________________________________________
Bool_t FT2::ApplyMSEloss(double x2X0, double xrho)
{
  // simulate random modification of track params due to the MS
  if (x2X0<=0) return kTRUE;
  //  printf("BeforeMAT: "); fProbe.Print();
  double alpha = fProbe.GetAlpha(); // store original alpha
  double mass2 = fProbe.fProbeMass*fProbe.fProbeMass;
  //
  double snp = fProbe.GetSnp();
  double dip = fProbe.GetTgl();
  Double_t angle=TMath::Sqrt((1.+ dip*dip)/((1-snp)*(1.+snp)));
  x2X0 *= angle;
  //
  static double covCorr[15],covDum[21]={0};
  static double mom[3],pos[3];
  double *cov = (double*) fProbe.GetCovariance();
  memcpy(covCorr,cov,15*sizeof(double));
  fProbe.GetXYZ(pos);
  fProbe.GetPxPyPz(mom);
  double pt2 = mom[0]*mom[0]+mom[1]*mom[1];
  double pt = TMath::Sqrt(pt2);
  double ptot2 = pt2 + mom[2]*mom[2];
  double ptot  = TMath::Sqrt(ptot2);
  double beta = ptot/TMath::Sqrt(ptot2 + mass2);
  double sigth = TMath::Sqrt(x2X0)*0.014/(ptot*beta);
  //
  // a la geant
  double phiSC = gRandom->Rndm()*TMath::Pi();
  double thtSC = gRandom->Gaus(0,1.4142*sigth);
  //  printf("MS phi: %+.5f tht: %+.5f\n",phiSC,thtSC);
  double sn = TMath::Sin(thtSC);
  double dx = sn*TMath::Sin(phiSC);
  double dy = sn*TMath::Cos(phiSC);
  double dz = TMath::Cos(thtSC);
  double v[3];
  //  printf("Before: %+.3e %+.3e %+.3e | MS: %+.3e %+.3e\n",mom[0],mom[1],mom[2],thtSC,phiSC);
  for (int i=3;i--;) mom[i] /= ptot;
  double vmm = TMath::Sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
  if (!IsZero(pt)) {
    double pd1 = mom[0]/vmm;
    double pd2 = mom[1]/vmm;
    v[0] = pd1*mom[2]*dx - pd2*dy + mom[0]*dz;
    v[1] = pd2*mom[2]*dx + pd1*dy + mom[1]*dz;
    v[2] = -vmm*dx                + mom[2]*dz;
  }
  else {
    v[0] = dx;
    v[1] = dy;
    v[2] = dz*TMath::Sign(1.,mom[2]);
  }
  //
  // account for eloss
  Double_t etot=TMath::Sqrt(ptot2 + mass2);
  Double_t dE=fProbe.BetheBlochSolid(ptot/fProbe.fProbeMass)*xrho;
  //  printf("X:%e E=%e dE=%e (xrho:%e)-> %e | ptot=%e M2 = %e\n",fProbe.GetX(),etot,dE,xrho,etot+dE,ptot,mass2);
  if ( TMath::Abs(dE) > 0.3*ptot ) {
#if DEBUG>1
    printf("StopEloss ");fProbe.Print();
#endif
    return kFALSE;
  } //30% energy loss is too much!
  etot += dE;
  if (etot<fProbe.fProbeMass+1./kRidiculous) {
#if DEBUG>1
    printf("StopEloss ");fProbe.Print();
#endif
    return kFALSE;
  } // stopped
  ptot = TMath::Sqrt(etot*etot - mass2);
  //
  double nrm = TMath::Sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  //  printf("before :%+e %+e %+e  || %+e %+e %+e %+e\n",mom[0],mom[1],mom[2],  sigth, x2X0, pt, beta);
  //  fProbe.Print();
  // direction cosines -> p
  for (int i=3;i--;) mom[i] = ptot*v[i]/nrm;
  //  printf("After : %+.3e %+.3e %+.3e\n",mom[0],mom[1],mom[2]);
  fProbe.Set(pos,mom,covDum,fProbe.Charge());
  //
  fProbe.RotateParamOnly(alpha);
  memcpy(cov,covCorr,15*sizeof(double));
	
  //  printf("AfterMAT: "); fProbe.Print();
  return kTRUE;
  //
}

//__________________________________________________________________
Bool_t FT2::GetRoadWidth(AliITSURecoLayer* lrA,double *pr)
{
  static AliExternalTrackParam sc;   // seed copy for manipulations
#if DEBUG>5
  printf(">>RW at Lr%d\n",lrA->GetActiveID());
#endif
  sc = fProbe;
  double xt;
  if (!sc.GetXatLabR(lrA->GetRMin(),xt,fBz,1)) return kFALSE;
  if (!sc.PropagateParamOnlyTo(xt,fBz)) return kFALSE;
  sc.GetXYZ(&pr[AliITSUTrackerGlo::kTrXIn]);
  pr[AliITSUTrackerGlo::kTrPhiIn] = TMath::ATan2(pr[AliITSUTrackerGlo::kTrYIn],pr[AliITSUTrackerGlo::kTrXIn]);
  if (!sc.RotateParamOnly(pr[AliITSUTrackerGlo::kTrPhiIn])) return kFALSE; // go to the frame of the entry point into the layer
  BringTo02Pi(pr[AliITSUTrackerGlo::kTrPhiIn]);
  double dr  = lrA->GetDR();                              // approximate X dist at the inner radius
  if (!sc.GetXYZAt(sc.GetX()-dr, fBz, pr + AliITSUTrackerGlo::kTrXOut)) {
    // special case: track does not reach inner radius, might be tangential
    double r = sc.GetD(0,0,fBz);
    double x;
    if (!sc.GetXatLabR(r,x,fBz,-1)) return kFALSE;
    dr = Abs(sc.GetX() - x);
    if (!sc.GetXYZAt(x, fBz, pr + AliITSUTrackerGlo::kTrXOut)) return kFALSE;
  }
  //
  pr[AliITSUTrackerGlo::kTrPhiOut] = ATan2(pr[AliITSUTrackerGlo::kTrYOut],pr[AliITSUTrackerGlo::kTrXOut]);
  BringTo02Pi(pr[AliITSUTrackerGlo::kTrPhiOut]);
  double sgy = 1e-6; // dummy spread
  double sgz = 1e-6; // dummy spread
  double phi0  = MeanPhiSmall(pr[AliITSUTrackerGlo::kTrPhiOut],pr[AliITSUTrackerGlo::kTrPhiIn]);
  double dphi0 = DeltaPhiSmall(pr[AliITSUTrackerGlo::kTrPhiOut],pr[AliITSUTrackerGlo::kTrPhiIn]);
  //
  pr[AliITSUTrackerGlo::kTrPhi0] = phi0;
  pr[AliITSUTrackerGlo::kTrZ0]   = 0.5*(pr[AliITSUTrackerGlo::kTrZOut]+pr[AliITSUTrackerGlo::kTrZIn]);
  dphi0 += sgy/lrA->GetR();
  pr[AliITSUTrackerGlo::kTrDPhi] =  dphi0<PiOver2() ? dphi0 : PiOver2();
  pr[AliITSUTrackerGlo::kTrDZ]   = 0.5*Abs(pr[AliITSUTrackerGlo::kTrZOut]-pr[AliITSUTrackerGlo::kTrZIn])   + sgz;
  //
#if DEBUG>5
  printf("<<RW at Lr%d: phi0: %e+-%e z0: %e+-%e\n",lrA->GetActiveID(),pr[AliITSUTrackerGlo::kTrPhi0],
	 pr[AliITSUTrackerGlo::kTrDPhi],pr[AliITSUTrackerGlo::kTrZ0],pr[AliITSUTrackerGlo::kTrDZ]);
#endif
  return kTRUE;
  //
}

//__________________________________________________________
Bool_t FT2::PassActiveITSLayer(AliITSURecoLayer* lr)
{
  // pass the layer and register the hits
  Double_t trImpData[AliITSUTrackerGlo::kNTrImpData];
  AliITSURecoSens *hitSens[AliITSURecoLayer::kMaxSensMatching];
  AliITSUGeomTGeo* gm = fITS->GetGeom();
  //
  int lrAID = lr->GetActiveID();
  ((double*)fKalmanOutward[lrAID].GetCovariance())[0] = -1; // invalidate outward track
  //
  if (!GetRoadWidth(lr, trImpData)) return kFALSE;
  fNITSHits[lrAID] = 0;
  const AliITSsegmentation* segm = gm->GetSegmentation(lrAID);
  int nsens = lr->FindSensors(&trImpData[AliITSUTrackerGlo::kTrPhi0], hitSens);
  int idxHit[2] = {0,1};
  //
#if DEBUG>5
  printf("Will test %d sensors on lr %d\n",nsens,lr->GetActiveID());
  for (int isn=0;isn<nsens;isn++) hitSens[isn]->Print();
  fProbe.Print();
#endif
  //
  if (nsens>1) { // in case of overlaping sensors need to traverse them in increasing R order
    double r2hit[AliITSURecoLayer::kMaxSensMatching];
    AliExternalTrackParam tmpt;
    for (int isn=0;isn<nsens;isn++) {
      AliITSURecoSens* sens =  hitSens[isn]; // estimate hit point, we need them ordered in R
      tmpt = fProbe;
      if (!tmpt.RotateParamOnly(sens->GetPhiTF())) {
#if DEBUG
	printf("testFailed to rotate for sens %d ",isn);sens->Print();
	tmpt.Print();
#endif
	return kFALSE;
      }
      double x = sens->GetXTF(), y;
      if (!tmpt.GetYAt(x,fBz,y)) {
#if DEBUG
	printf("tesdFailed to do GetYAt for sensor %d ",isn);sens->Print();
	tmpt.Print();
#endif
	return kFALSE;
      }
      r2hit[isn] = x*x+y*y;
    }
    if (r2hit[0]>r2hit[1]) {idxHit[0]=1;idxHit[1]=0;} // only 2 hits per layer are possible
    nsens=2;
  }
  //
  for (int isn=0;isn<nsens;isn++) {
    AliITSURecoSens* sens =  hitSens[idxHit[isn]];
    Bool_t res = fUseKalmanOut ? fProbe.Rotate(sens->GetPhiTF()) : fProbe.RotateParamOnly(sens->GetPhiTF());
    if (!res) {
#if DEBUG
      printf("Failed to rotate for sens %d ",isn); sens->Print();
      fProbe.Print();
#endif
      return kFALSE;
    }
    if (!PropagateToX(sens->GetXTF(),1,kFALSE,fSimMat,fSimMat)) {
      return kFALSE;
    }
    const TGeoHMatrix* mt = gm->GetMatrixT2L(sens->GetID());
    double xyzL[3],xyzT[3] = {fProbe.GetX(),fProbe.GetY(),fProbe.GetZ()};
    mt->LocalToMaster(xyzT,xyzL);
    int ix,iz;
    if (!segm->LocalToDet(xyzL[0],xyzL[2],ix,iz)) {
#if DEBUG>5
      double xyzt[3];
      double ph = TMath::ATan2(xyzt[1],xyzt[0]);
      BringTo02Pi(ph);
      printf("NoHit with XZloc %+.4f %+.4f, phi %.4f\n",xyzL[0],xyzL[2],ph);
      fProbe.Print();
      sens->Print();
      segm->Print();
#endif
      continue; // not in the active zone
    }
    // register hit
    int nh = fNITSHits[lrAID];
    fITSSensHit[lrAID][nh] = sens;
    fITSHitYZ[lrAID][nh][0] = fProbe.GetY();
    fITSHitYZ[lrAID][nh][1] = fProbe.GetZ();
    fNITSHits[lrAID]++;
    if (fNITSHits[lrAID]) fNITSLrHit++;
#if DEBUG>5
    printf("Register hit %d at lr%d\n",fNITSHits[lrAID],lr->GetActiveID());
    fProbe.Print();
#endif
    // emulate outward Kalman track
    if (fUseKalmanOut) {
      AliExternalTrackParam& trKO = fKalmanOutward[lrAID];
      if (((double*)trKO.GetCovariance())[0]<0) {
	trKO = fProbe;
	AliExternalTrackParam trKOUp = fProbe;
	double chi = UpdateKalman(&trKOUp,fITSHitYZ[lrAID][nh][0],fITSHitYZ[lrAID][nh][1],
				  fSigYITS,fSigZITS,kTRUE,kFALSE);
	if (chi<0) {
	  printf("Failed on emulating outward Kalman on lr %d\n",lrAID);
	  trKOUp.Print();
	  return kFALSE;
	}
	// assign updated error matrix to fProbe
	memcpy((double*)fProbe.GetCovariance(),trKOUp.GetCovariance(),15*sizeof(double));
      } // update for outward Kalman
    }
  }
  //
  return kTRUE;
}
//__________________________________________________________
AliPIDResponse::EDetPidStatus FT2::GetComputeTPCProbability (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const
{
  // Similar as AliPIDResponse
  // Compute PID response for the TPC
  //
  // set flat distribution (no decision)
  for (Int_t k=0; k<nSpecies; k++) p[k]=1./nSpecies;
	
	
  Double_t dedx=track->GetTPCsignal();
  Bool_t mismatch=kTRUE/*, heavy=kTRUE*/;
	
  Double_t bethe = 0.;
  Double_t sigma = 0.;
	
  for (Int_t j=0; j<nSpecies; j++) {
    AliPID::EParticleType type=AliPID::EParticleType(j);
	 
    bethe=fPIDResponse->GetTPCResponse().GetExpectedSignal(track,type,AliTPCPIDResponse::kdEdxDefault,kFALSE, kFALSE);
    sigma=fPIDResponse->GetTPCResponse().GetExpectedSigma(track, type,AliTPCPIDResponse::kdEdxDefault,kFALSE, kFALSE);
#if DEBUG>5
    AliInfo(Form("Species Hypothesis %i - dEdx: %f - Bethe :%f +/- %f",j,dedx,bethe,sigma));
#endif
    if (TMath::Abs(dedx-bethe) > 5.*sigma) {
      p[j]=TMath::Exp(-0.5*5.*5.)/sigma;
    } else {
      p[j]=TMath::Exp(-0.5*(dedx-bethe)*(dedx-bethe)/(sigma*sigma))/sigma;
      mismatch=kFALSE;
    }
  }
  if (mismatch){
    for (Int_t j=0; j<nSpecies; j++) p[j]=1./nSpecies;
  }
  return AliPIDResponse::kDetPidOk;
}

//__________________________________________________________
Bool_t FT2::DiagonalizeErrors(const double* cov, double &sy2d, double &sz2d) 
{
  // diagonalize cov. matrix
  double dd = cov[0]-cov[2]; 
  dd = TMath::Sqrt(dd*dd + 4.*cov[1]*cov[1]);
  double sd = cov[0]+cov[2];
  sy2d = 0.5*(sd - dd);
  sz2d = 0.5*(sd + dd);
  if (sy2d<=0 || sz2d<=0) return kFALSE;
  return kTRUE;
}

//_____________________________________________
Double_t FT2::HitDensity(double r2, double tgl) const
{
  // hit density from single collision at radius r
  double den = fdNdY/(2.*TMath::Pi()*r2);
  return den/TMath::Sqrt(1 + tgl*tgl);
}

//_____________________________________________
Bool_t FT2::BiasAsFake(double yz[2], const double* extyz, const double *cov) const
{
  // assign instead of the "true" yz hit position a fake position such that it
  // will have better chi2 wrt the extrapolation point
  Double_t r00=cov[0], r01=cov[1], r11=cov[2];
  Double_t det=r00*r11 - r01*r01;
  if (TMath::Abs(det) < 1e-16) return kFALSE;
  Double_t tmp=r00; r00=r11/det; r11=tmp/det; r01=-r01/det;
  double dy0 = yz[0]-extyz[0];
  double dz0 = yz[1]-extyz[1];
  double d2max = dy0*dy0*r00+dz0*dz0*r11+2.*dy0*dz0*r01;
  dy0 = TMath::Sqrt(d2max*cov[0]);
  dz0 = TMath::Sqrt(d2max*cov[2]);
  double d2 = 1e9;
  //  printf("d2max: %f dy:%f dz:%f\n",d2max,dy0,dz0);
  double dy(dy0),dz(dz0);
  while(d2>=d2max) {
    dy = (gRandom->Rndm()-0.5)*2*dy0;
    dz = (gRandom->Rndm()-0.5)*2*dz0;
    d2 = dy*dy*r00+dz*dz*r11+2.*dy*dz*r01;
  }
  //
  printf("AddFake DY:%f DZ:%f Chi2Max:%f Chi2F:%f\n",dy,dz, d2max,d2);
  yz[0] = extyz[0]+dy;
  yz[1] = extyz[1]+dz;
  return kTRUE;
}
