#ifndef FT2_H
#define FT2_H
#include <TObject.h>
#include "AliExternalTrackParam.h"
#include "AliITSUAux.h"
#include "AliESDpid.h"
#include "AliPIDResponse.h"

//#define DEBUG

class AliITSUReconstructor;
class AliITSURecoDet;
class AliITSURecoLayer;
class AliITSURecoSens;
class TParticle;
class AliVertex;

const float kRidiculous = 999999;
const float kMaxZTPC = 250;
const double kTrackToler = 1.e-6; // tracking tolerance

class FTProbe : public AliExternalTrackParam
{
 public:
  FTProbe(){};
  virtual ~FTProbe(){};
  virtual Double_t GetTPCmomentum() const {return fTPCmomentum;}
  virtual Double_t GetTPCsignal() const {return fTPCSignal;}
  virtual UShort_t GetTPCsignalN() const {return fTPCSignalN;}
  virtual Int_t GetTPCTrackingPID() const {return fAbsPdgCodeForTracking;}
	
 protected:
  Double_t fProbeMass;	// true mass
  Double_t fTPCmomentum;	// momentum after TPC reconstruction
  Double_t fTPCSignal;	// TPC signal
  UShort_t fTPCSignalN;
  Double_t fAbsPdgCode;	// pdg code of particle
  Int_t fAbsPdgCodeForTracking;
  ClassDef(FTProbe,1)
    };


class FT2 : public TObject
{
 public:
  enum {kMaxITSLr=7, kMaxHitPerLr=2};
  struct FT2TPCLayer {
  FT2TPCLayer(float xr=0,float x2x=0,float resRPhi=0,float resZ=0,float effL=0) :
    x(xr),x2x0(x2x),rphiRes(resRPhi),zRes(resZ),eff(effL),isDead(resRPhi>kRidiculous||resZ>kRidiculous),
      hitY(0),hitZ(0),hitSect(-1) {if (isDead) effL=-1;}
    Float_t  x;
    Float_t  x2x0;
    Float_t  rphiRes;
    Float_t  zRes;
    Float_t  eff;
    Bool_t   isDead;
    //
    Float_t  hitY;
    Float_t  hitZ;
    Float_t  hitSect;
  };
  typedef struct FT2TPCLayer FT2TPCLayer_t;
	
	
  FT2();
  virtual ~FT2();
	
  void InitEnvLocal();
  void InitTPCSignalFiles(const char *TPCsignalElectrons,const char *TPCsignalMuons,const char *TPCsignalPions,const char *TPCsignalKaons,const char *TPCsignalProtons);
  void InitTPCPIDResponse();
  void InitDetector(Bool_t addTPC=kTRUE, Float_t sigYTPC=0.1, Float_t sigZTPC=0.1,
		    Float_t effTPC=0.99, Float_t scEdge=2.6);
  //
  void PrintLayout();
  //
  Bool_t ProcessTrack(TParticle* trc, AliVertex* vtx);
  void   SetSimMat(Bool_t v=kTRUE) {fSimMat = v;}
  void   SetMinTPCHits(int v=60)  {fMinTPCHits=v;}
  void   SetMinITSLrHit(int v=4)  {fMinITSLrHit=v;}
  void   SetUsePIDForTracking(Bool_t usePID) {fUsePIDForTracking=usePID;}
  //
  Int_t    GetNClITSFakes()  const {return fNClITSFakes;}
  Int_t    GetNClITS()  const {return fNClITS;}
  Int_t    GetNClTPC()  const {return fNClTPC;}
  Double_t GetChi2ITS() const {return fChi2ITS;}
  Double_t GetChi2TPC() const {return fChi2TPC;}
  FTProbe& GetProbe() const {return (FTProbe&)fProbe;}
  AliExternalTrackParam& GetKalmanOut(int i) {return (AliExternalTrackParam&)fKalmanOutward[i];}
  //
  const Double_t* GetDCA()    const {return &fDCA[0];}
  const Double_t* GetDCACov() const {return &fDCACov[0];}
  static float GetMaxStepTGeo() {return fgMaxStepTGeo;}
  static void  SetMaxStepTGeo(float stp=1.) {fgMaxStepTGeo = stp;}
  //
  void     SetdNdY(double v=-1) {fdNdY = v;}
  Double_t GetdNdY() const {return fdNdY;}
  Double_t HitDensity(double r2, double tgl) const;
  Bool_t   BiasAsFake(double yz[2], const double* extyz, const double *cov) const;
  Bool_t DiagonalizeErrors(const double *cov, double &sy2d, double &sz2d);
  //
 protected:
  void AddTPC(Float_t sigY=0.1, Float_t sigZ=0.1, Float_t eff=0.99, Float_t scEdge=2.6);
  void AddTPCLayer(Float_t x, Float_t x2x0,Float_t sigY, Float_t sigZ, Float_t eff);
  Bool_t InitProbe(TParticle* trc);
  Bool_t PrepareProbe();
  Bool_t ApplyMSEloss(double x2X0, double xrho);
  Bool_t PropagateToX(double xTgt, int dir,Bool_t propErr,Bool_t simMat,Bool_t useTGeo);
  Bool_t PropagateToR(double xTgt, int dir,Bool_t propErr,Bool_t simMat,Bool_t useTGeo);
  Bool_t IsZero(double val, double tol=1e-9) const {return TMath::Abs(val)<tol;}
  Bool_t PassActiveITSLayer(AliITSURecoLayer* lr);
  Bool_t GetRoadWidth(AliITSURecoLayer* lrA,double *pr);
  void   ResetCovMat(AliExternalTrackParam* trc);
  Double_t UpdateKalman(AliExternalTrackParam* trc, double y,double z,double sigY,double sigZ,Bool_t randomize=kTRUE,Bool_t fake=kFALSE);
  Bool_t ReconstructProbe();
  AliPIDResponse::EDetPidStatus GetComputeTPCProbability (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const;
  //
 protected:
  AliITSUReconstructor* fITSRec; // interface for reconstructor
  AliITSURecoDet* fITS;          // interface for ITS
  Bool_t fUsePIDForTracking;
  Bool_t fIsTPC;                 // TPC added
  AliPIDResponse* fPIDResponse;
  TFile *fTPCSignalElectron;
  TFile *fTPCSignalMuon;
  TFile *fTPCSignalPion;
  TFile *fTPCSignalKaon;
  TFile *fTPCSignalProton;
  Float_t fTPCSectorEdge;        // cut in cm on sector edge
  Double_t fMaxSnpTPC;           // stop particle if snp>fMaxSnpTPC
  std::vector<FT2TPCLayer_t> fTPCLayers;
  std::vector<int>           fTPCHitLr;
  Int_t                      fNTPCHits; //!
  Int_t                      fMinTPCHits; // require min amount of TPC hits
  Int_t                      fMinITSLrHit; // require min amount of ITS lr hit
  Double_t fDCA[2],fDCACov[3];   //! dca to vertex and its covariance
  //
  FTProbe fProbe;  // track
  AliExternalTrackParam* fKalmanOutward; //! parameters of outward kalman 
  //Double_t              fProbeMass; // probe mass
  Double_t              fBz;        // bz
  Bool_t                fSimMat;    // simulate material effects in probe preparation
  //
  Int_t                 fNClTPC;    //! N used TPC clusters
  Int_t                 fNClITS;    //! N used ITS clusters
  Int_t                 fNClITSFakes;    //! N used ITS Fake clusters
  Double_t              fChi2TPC;   //! total chi2 in TPC
  Double_t              fChi2ITS;   //! total chi2 in ITS
  //
  // hit info in the ITS
  Double_t fSigYITS,fSigZITS;       // nominal ITS later resolution
  Int_t fNITSHits[kMaxITSLr]; //! n hits per ITS layer
  Int_t fNITSLrHit;           //! n of ITS layers whith hit
  AliITSURecoSens* fITSSensHit[kMaxITSLr][2]; //! hit sensors
  Double_t fITSHitYZ[kMaxITSLr][kMaxHitPerLr][2]; //! tracking Y,Z of each hit
  //
  Double_t fdNdY;             //! if positive, use it for fakes simulation
  //
  static float fgMaxStepTGeo; // max step for tracking accounting for TGeo materials
	
  ClassDef(FT2,1)
};

#endif
