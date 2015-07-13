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
#include "AliTPCcalibDB.h"
#include "AliTPCParam.h"
#include "AliTPCClusterParam.h"
#include "AliTPCRecoParam.h"
#include "AliPIDResponse.h"
#include "AliDetectorPID.h"
#include "TH2F.h"
#include "TF3.h"
#include "AliNDLocalRegression.h"
#include "AliMathBase.h"

ClassImp(FTProbe)
ClassImp(FT2)

TH2F* hfake = 0;

using namespace AliITSUAux;

float FT2::fgMaxStepTGeo = 1.0;

//________________________________________________
FTProbe::FTProbe()
:fProbeMass(0.14)
,fTPCmomentum(0)
,fTPCSignal(0)
,fTPCSignalN(0)
,fAbsPdgCode(0)
,fPdgCode(0)
,fAbsPdgCodeForTracking(211)
,fTrueMass(0)
,fProbeNClTPC(0)
,fProbeNClITS(0)
,fProbeNClITSFakes(0)
,fProbeITSPatternFake(0)
,fProbeITSPattern(0)
,fProbeChi2TPC(0)
,fProbeChi2ITS(0)
,fIsDecayed(0)
,fIsAbsorbed(0)
,fDecayRadius(0)
,fAbsorbtionRadius(0)
,fLostInItsTpcMatching(0)
,fTrackToClusterChi2CutITS(0)
,fProbeZAtCutOffCheck(0)
{
	// def. c-tor
}

//________________________________________________
FT2::FT2() :
fITSRec(0)
,fITS(0)
,fUsePIDForTracking(0)
,fRunNumber("")
,fIsTPC(kFALSE)
,fPIDResponse(0)
,fTPCParaFile(0)
,fXSectionFile(0)
,fTPCSectorEdge(0)
,fMaxSnpTPC(0.95)
,fTPCInnerRows(63)
,fTPCMiddleRows(64)
,fTPCOuterRows(32)
,fTPCLayers()
,fTPCHitLr()
,fNTPCHits(0)
,fMinTPCHits(0)
,fMinITSLrHit(0)
,fTrueMCtrackMult(-1)
,fProbe()
,fProbeIni()
,fKalmanOutward(0)
,fUseKalmanOut(kTRUE)
,fBz(0.5)
,fSimMat(kFALSE)
,fAllowDecay(kFALSE)
,fCurrITSLr(-1)
,fNClTPC(0)
,fNClITS(0)
,fNClITSFakes(0)
,fITSPatternFake(0)
,fITSPattern(0)
,fChi2TPC(0)
,fChi2ITS(0)
,fTPCMap()
,fSigYITS(3.14e-4) // 5e-4
,fSigZITS(3.38e-4) // 5e-4
,fITSerrSclY(2.0)
,fITSerrSclZ(1.0)
,fNITSLrHit(0)
,fdNdY(-1)
,fTPCClsLossProbIROC(0)
,fTPCClsLossProbOROCmedium(0)
,fTPCClsLossProbOROClong(0)
,fTPCDistortionRPhi(0)
,fTPCDistortionR(0)
,fTPCfMCChi2(0)
,fTPCft2Chi2(0)
,fItsTpcMatchingfMC(0)
,fItsTpcMatchingft2(0)
,fTPCClusterErrorParamY(0)
,fTPCClusterErrorParamZ(0)
,fTPCParam(0)
,fTPCRecoParam(0)
,fAllocCorrelatedITSFakes(kTRUE)
{
	//
	AliInfo("Default Constructor");
	//
	memset(fITSSensHit,0,kMaxITSLr*2*sizeof(AliITSURecoSens*));
	//
	const double accAmp[kMaxITSLr]  = {0.017,0.013,0.013,0.014,0.016,0.015,0.015};
	const double accSigY[kMaxITSLr] = { 9e-2, 9e-2, 9e-2, 9e-2, 9e-2, 9e-2, 9e-2};
	const double accSigZ[kMaxITSLr] = { 9e-2, 9e-2, 9e-2, 9e-2, 9e-2, 9e-2, 9e-2};
	//
	float c0tr2clChi2[7]	= {20,25,30,40,45,45,70};	// cut on cluster to track chi2
	float c0gloChi2[7]		= {6,10,20,30,60,60,70};	// cut on seed global norm chi2
	float c0missPen[7]		= {2.,2.,2.,2.,2.,2.,2.};	// missing cluster penalty
	
	for (int i=0;i<kMaxITSLr;i++) {
		fNCorrelITSFakes[i] = accAmp[i];		// integral of accompanying hits
		fCorrelITSFakesSigY[i] = accSigY[i];	// width in rphi
		fCorrelITSFakesSigZ[i] = accSigZ[i];	// width in z
		fC0tr2clChi2[i] = c0tr2clChi2[i];		// cut on cluster to track chi2
		fC0gloChi2[i]	= c0gloChi2[i];		// cut on seed global norm chi2
		fC0missPen[i]	= c0missPen[i];		// missing cluster penalty
		
	}
	fTpcPidSignal[0]=fTpcPidSignal[1]=fTpcPidSignal[2]=fTpcPidSignal[3]=fTpcPidSignal[4] = 0;
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
void FT2::InitTPCParaFile(const char *TPCParaFile)
{
	AliInfo("Setting TPC Files");
	
	fTPCParaFile = TFile::Open(TPCParaFile);
	if(fTPCParaFile->IsZombie()) AliFatal("Problem with opening TPC Parameterization File - File not available!");
	//
	fTPCClsLossProbIROC				= (AliNDLocalRegression*)fTPCParaFile->Get("tpcNclProb_IROC");
	fTPCClsLossProbOROCmedium = (AliNDLocalRegression*)fTPCParaFile->Get("tpcNclProb_OROCmedium");
	fTPCClsLossProbOROClong		= (AliNDLocalRegression*)fTPCParaFile->Get("tpcNclProb_OROClong");
	//
	fTpcPidSignal[0] = (AliNDLocalRegression*)fTPCParaFile->Get("tpcPidSignal_electrons");
	fTpcPidSignal[1] = (AliNDLocalRegression*)fTPCParaFile->Get("tpcPidSignal_muons");
	fTpcPidSignal[2] = (AliNDLocalRegression*)fTPCParaFile->Get("tpcPidSignal_pions");
	fTpcPidSignal[3] = (AliNDLocalRegression*)fTPCParaFile->Get("tpcPidSignal_kaons");
	fTpcPidSignal[4] = (AliNDLocalRegression*)fTPCParaFile->Get("tpcPidSignal_protons");
	//
	if(AliTrackerBase::GetBz()<0){
		// B--: parameterization of 137544
		fTPCDistortionRPhi	= (AliNDLocalRegression*)fTPCParaFile->Get("tpcExBDeltaRPhi_Bminus");
		fTPCDistortionR			= (AliNDLocalRegression*)fTPCParaFile->Get("tpcExBDeltaR_Bminus");
	}
	else{
		// B++: parameterization of 138871
		fTPCDistortionRPhi	= (AliNDLocalRegression*)fTPCParaFile->Get("tpcExBDeltaRPhi_Bplus");
		fTPCDistortionR			= (AliNDLocalRegression*)fTPCParaFile->Get("tpcExBDeltaR_Bplus");
	}
	
	fTPCfMCChi2 = (AliNDLocalRegression*)fTPCParaFile->Get("fmcTpcChi2");
	fTPCft2Chi2 = (AliNDLocalRegression*)fTPCParaFile->Get("ft2TpcChi2");

	fItsTpcMatchingfMC = (AliNDLocalRegression*)fTPCParaFile->Get("fmcITSTPCmatching");
	fItsTpcMatchingft2 = (AliNDLocalRegression*)fTPCParaFile->Get("ft2ITSTPCmatching");

	fTpcClusterAcc = (TF1*)fTPCParaFile->Get("tpcClusterAcc");
	
	//	fTpcClusterAcc = new TF1("fTpcClusterAcc","[0]+[1]/x",80,260);
	//	fTpcClusterAcc->SetParameters(1.05,-25.1827);
}
//________________________________________________
void FT2::InitXSectionFile(const char *XSectionFile)
{
	AliInfo("Setting X-Section Files");
	
	fXSectionFile = TFile::Open(XSectionFile);
	if(fXSectionFile->IsZombie()) AliFatal("Problem with opening X-Section File - File not available!");
	
	fXSectionHp[0] = (TH1F*)fXSectionFile->Get("hPiplusPXSection");
	fXSectionHp[1] = (TH1F*)fXSectionFile->Get("hPiminusPXSection");
	fXSectionHp[2] = (TH1F*)fXSectionFile->Get("hKplusPXSection");
	fXSectionHp[3] = (TH1F*)fXSectionFile->Get("hKminusPXSection");
	fXSectionHp[4] = (TH1F*)fXSectionFile->Get("hPplusPXSection");
	fXSectionHp[5] = (TH1F*)fXSectionFile->Get("hPminusPXSection");
}
//________________________________________________
void FT2::InitDetector(Bool_t addTPC,Float_t scEdge)
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
	if (addTPC) AddTPC(scEdge);
	int nLrITS = fITS->GetNLayersActive();
	fKalmanOutward = new AliExternalTrackParam[nLrITS];
	//
	AliTPCcalibDB * calib = AliTPCcalibDB::Instance();
	const AliMagF * field = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
	calib->SetExBField(field);
	calib->SetRun(fRunNumber.Atoi());
	// get cluster param from calibDB
	AliTPCClusterParam *clsParam = calib->GetClusterParam();
	clsParam->Instance();
	clsParam->SetInstance(clsParam);
//	fTPCClusterErrorParam = new TF1("fTPCClusterErrorParam","AliTPCClusterParam::SGetError0Par([0],[1],x,[2])",0.,250.);
	fTPCClusterErrorParamY = new TF1("fTPCClusterErrorParamY","clusterResolutionYAliRoot(x,[0],[1])",0,250.);
	fTPCClusterErrorParamZ = new TF1("fTPCClusterErrorParamZ","clusterResolutionZAliRoot(x,[0],[1])",0,250.);

	// get tpc param from calibDB
	fTPCParam = calib->GetParameters();
	// get tpc reco param from calibDB
	fTPCRecoParam = calib->GetRecoParam(1); // to be checked with marian
}
//____________________________________________________
void FT2::InitTPCPIDResponse()
{
	AliInfo("Initializing TPC PID Response");
	
	fPIDResponse = new AliPIDResponse(kTRUE);
	fPIDResponse->GetTPCResponse().SetUseDatabase(kFALSE);
	//AliTPCParam* param = AliTPCcalibDB::Instance()->GetParameters();
	if (fTPCParam){
		TVectorD *paramBB = fTPCParam->GetBetheBlochParameters();
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
//	man->SetSpecificStorage("ITS/*/*", "alien://folder=/alice/simulation/LS1_upgrade/Ideal");
//	man->SetSpecificStorage("TPC/*/*", "alien://folder=/alice/simulation/2008/v4-15-Release/Ideal/");
	
	man->SetSpecificStorage("ITS/Align/Data", "alien://folder=/alice/simulation/LS1_upgrade/Ideal");
	man->SetSpecificStorage("ITS/Calib/RecoParam", "alien://folder=/alice/simulation/LS1_upgrade/Ideal");
	man->SetDefaultStorage("alien://folder=/alice/data/2010/OCDB");
	man->SetSpecificStorage("MUON/Align/Data","alien://folder=/alice/simulation/2008/v4-15-Release/Residual/");
	man->SetSpecificStorage("TPC/Align/Data", "alien://folder=/alice/simulation/2008/v4-15-Release/Residual/");
	man->SetSpecificStorage("TPC/Calib/ClusterParam", "alien://folder=/alice/simulation/2008/v4-15-Release/Residual/");
	man->SetSpecificStorage("TPC/Calib/RecoParam", "alien://folder=/alice/simulation/2008/v4-15-Release/Residual/");
	man->SetSpecificStorage("TPC/Calib/TimeGain", "alien://folder=/alice/simulation/2008/v4-15-Release/Residual/");
	man->SetSpecificStorage("TPC/Calib/AltroConfig", "alien://folder=/alice/simulation/2008/v4-15-Release/Residual/");
	man->SetSpecificStorage("TPC/Calib/TimeDrift", "alien://folder=/alice/simulation/2008/v4-15-Release/Residual/");
	man->SetSpecificStorage("TPC/Calib/Correction", "alien://folder=/alice/simulation/2008/v4-15-Release/Residual/");

	//	man->SetSpecificStorage("GRP/GRP/Data","alien://folder=/alice/data/2010/OCDB");
	//	man->SetSpecificStorage("ITS/Align/Data", "alien://folder=/alice/simulation/LS1_upgrade/Ideal");
	//	man->SetSpecificStorage("ITS/Calib/RecoParam", "alien://folder=/alice/simulation/LS1_upgrade/Ideal");
	man->SetRun(fRunNumber.Atoll());
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
void FT2::AddTPC(Float_t scEdge)
{
	// add TPC mock-up
	//
	const Float_t kRadLPerRow = 0.000036;
	//
	const Float_t kTPCInnerRadialPitch  =    0.75 ;    // cm
	const Float_t kTPCMiddleRadialPitch =    1.0  ;    // cm
	const Float_t kTPCOuterRadialPitch  =    1.5  ;    // cm
	const Int_t   kTPCInnerRows            =   fTPCInnerRows;		// 63
	const Int_t   kTPCMiddleRows           =   fTPCMiddleRows;	// 64
	const Int_t   kTPCOuterRows            =   fTPCOuterRows;		// 32
	const Int_t   kTPCRows       =   (kTPCInnerRows + kTPCMiddleRows + kTPCOuterRows) ;
	const Float_t kTPCRowOneRadius         =   85.2  ;    // cm
	const Float_t kTPCRow64Radius          =  135.1  ;    // cm
	const Float_t kTPCRow128Radius         =  199.2  ;    // cm
	Float_t pitch;
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
		if (k<kTPCInnerRows) {
			rowRadius =  kTPCRowOneRadius + k*kTPCInnerRadialPitch ;
			pitch = kTPCInnerRadialPitch;
		}
		else if ( k>=kTPCInnerRows && k<(kTPCInnerRows+kTPCMiddleRows) ){
			rowRadius =  kTPCRow64Radius + (k-kTPCInnerRows+1)*kTPCMiddleRadialPitch ;
			pitch = kTPCMiddleRadialPitch;
		}
		else if (k>=(kTPCInnerRows+kTPCMiddleRows) && k<kTPCRows ){
			rowRadius = kTPCRow128Radius + (k-kTPCInnerRows-kTPCMiddleRows+1)*kTPCOuterRadialPitch ;
			pitch = kTPCOuterRadialPitch;
		}
		AddTPCLayer(k,rowRadius,kRadLPerRow,pitch);
	}
	//
	fTPCHitLr.resize(fTPCLayers.size());
}

//____________________________________________________
void FT2::AddTPCLayer(Int_t rowId,Float_t x, Float_t x2x0, Float_t pitch)
{
	// add single TPC layer
	fTPCLayers.push_back(FT2TPCLayer(rowId,x,x2x0,pitch));
}

//____________________________________________________
void FT2::PrintLayout()
{
	// print setup
	if (fITS) fITS->Print("lr");
	if (fIsTPC) {
		printf("TPC inactive sector edge: %.3f\n",fTPCSectorEdge);
		printf("TPC \t  R   \tx2x0\t Pitch  \n");
		for (int ilr=0;ilr<(int)fTPCLayers.size();ilr++) {
			FT2TPCLayer& lr = fTPCLayers[ilr];
			printf("%3d\t%6.2f\t%1.6f\t",ilr,lr.x,lr.x2x0);
			//
			if (lr.isDead) printf("      \t");
			else printf("%6.4f\t",lr.pitch);
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
	if(fTrueMCtrackMult==-1){AliFatal("True MC Multiplicity is 0. This must not happen! Set multiplicity before FT2::InitEnvLocal()!\n");}
	
	
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
	if (!ReconstructProbe(part)) {
#if DEBUG
		printf("Track reconstruction failed\n");
		fProbe.Print();
#endif
		return kFALSE;
	}
	
	// check if tracks are rejected by ITS chi2 per layer
	fProbe.fTrackToClusterChi2CutITS = CutOnTrackToClusterChi2ITS();
	
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
	
	if(fTuneOnDataOrMC){
#if DEBUG
		printf("\n\n\n################################");
		printf("### TuneOnDataOrMC is Active ###\n# TPC min. cluster? %d < %d \n# ITS kAny? %d \n",fNClTPC,30,(fITSPattern&(0x1<<0) || (fITSPattern&(0x1<<1) && fITSPattern&(0x1<<2))));
#endif
		//
		if(fNClTPC<30 || (!(fITSPattern&(0x1<<0)) && !(fITSPattern&(0x1<<1) && fITSPattern&(0x1<<2)))){
#if DEBUG
			printf("################################\n\n\n");
#endif
			fProbe.fLostInItsTpcMatching=kTRUE;
		}
		//
		Double_t pointChi[3] = {1./AliMathBase::BetheBlochAleph(0.0001+part->P()/TDatabasePDG::Instance()->GetParticle(fProbe.fAbsPdgCode)->Mass()),
			abs(part->Eta()),fTrueMCtrackMult};
		//
		Double_t pointMatch[3] = {fNClTPC,TMath::Abs(part->Eta()),(1./part->Pt())};
		//
		Double_t templateFmcMatchingEff = fItsTpcMatchingfMC->Eval(pointMatch);
		//
#if DEBUG
		Double_t templateFt2MatchingEff = fItsTpcMatchingft2->Eval(pointMatch);
		printf("# fMC ITS+TPC match.? %0.2f\n# FT2 ITS+TPC match.? %0.2f\n",templateFmcMatchingEff,templateFt2MatchingEff);
#endif
		if(gRandom->Rndm()>templateFmcMatchingEff){
#if DEBUG
			printf("################################\n\n\n");
#endif
			fProbe.fLostInItsTpcMatching=kTRUE;
		}
		//
		Double_t templateFmcChi2stand = fTPCfMCChi2->Eval(pointChi); // TPC standardized chi2 from either full MC or data
		Double_t templateFt2Chi2stand = fTPCft2Chi2->Eval(pointChi); // TPC standardized chi2 from ft2 first iteration
		if(templateFt2Chi2stand==0)templateFt2Chi2stand=templateFmcChi2stand;
		fChi2TPC *=(templateFmcChi2stand/templateFt2Chi2stand);
#if DEBUG
		Double_t oldChi2 = fChi2TPC;
		printf("# fMC chi2? %0.2f\n# FT2 chi2? %0.2f\n# Old Chi2? %0.2f\n# New Chi2? %0.2f\n",templateFmcChi2stand,templateFt2Chi2stand,oldChi2,fChi2TPC);
		printf("################################\n\n\n");
#endif
	}
	
	fProbe.fProbeNClTPC					= fNClTPC;
	fProbe.fProbeNClITS					= fNClITS;
	fProbe.fProbeNClITSFakes		= fNClITSFakes;
	fProbe.fProbeITSPatternFake	= fITSPatternFake;
	fProbe.fProbeITSPattern			= fITSPattern;
	fProbe.fProbeChi2TPC				= fChi2TPC;
	fProbe.fProbeChi2ITS				= fChi2ITS;
	
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
	Int_t pdgCode				= part->GetPdgCode();
	TParticlePDG* pdgp	= TDatabasePDG::Instance()->GetParticle(pdgCode);
	Double_t charge			= pdgp->Charge();
	fProbe.fProbeMass		= pdgp->Mass();
	fProbe.fAbsPdgCode	= TMath::Abs(pdgCode);
	fProbe.fPdgCode			= pdgCode;
	fProbe.fAbsPdgCodeForTracking = 0;
	fProbe.fTrueMass		= fProbe.fProbeMass;
	//
	fProbe.fProbeNClTPC					= fNClTPC;
	fProbe.fProbeNClITS					= fNClITS;
	fProbe.fProbeNClITSFakes		= fNClITSFakes;
	fProbe.fProbeITSPatternFake	= fITSPatternFake;
	fProbe.fProbeITSPattern			= fITSPattern;
	fProbe.fProbeChi2TPC				= fChi2TPC;
	fProbe.fProbeChi2ITS				= fChi2ITS;
	fProbe.fIsDecayed = fProbe.fIsAbsorbed = 0;
	fProbe.fAbsorbtionRadius = fProbe.fDecayRadius = 0.;
	fProbe.fTrackToClusterChi2CutITS = 0;
	fProbe.fProbeZAtCutOffCheck = 0.;
	fProbe.fLostInItsTpcMatching = 0;
	//
	fProbe.fTPCSignal = 0.;
	fProbe.fTPCSignalN = 0;
	//
/*
	fProbe.fInnerTrackParameters[7] = {0.};
	fProbe.fProbeITSClusterIsCls[8] = {0.};
	fProbe.fProbeITSClusterX[8] = {0.};
	fProbe.fProbeITSClusterY[8] = {0.};
	fProbe.fProbeITSClusterZ[8] = {0.};
	fProbe.fProbeTPCClusterIsCls[160] = {0.};
	fProbe.fProbeTPCClusterX[160] = {0.};
	fProbe.fProbeTPCClusterY[160] = {0.};
	fProbe.fProbeTPCClusterZ[160] = {0.};
	chiwITS[7] = {0.};
*/
	
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
	fTPCMap.ResetAllBits();
	fProbeIni = fProbe;
	return kTRUE;
}

//____________________________________________________
Int_t FT2::ProbeDecayAbsorb(double* posIni)
{
	// check if the probe has decayed (return 1) or absorbed (return -1).
	// if survived, return 0
	// Before returning, assign to posIni current position
	double posCurr[3];
	fProbe.GetXYZ(posCurr);
	double params[8];
	AliTrackerBase::MeanMaterialBudget(posIni,posCurr,params);
	Double_t dist=params[4];
	//
#if DEBUG>5
	printf("New step length %.2f from XYZ= %+.2f %+.2f %+.2f\n",dist,posIni[0],posIni[1],posIni[2]);
#endif
	//
	if(gRandom->Rndm()<ParticleDecayProbability(dist)) {
		fProbe.fIsDecayed=kTRUE;
		fProbe.fIsAbsorbed=kFALSE;
		fProbe.fDecayRadius=TMath::Sqrt(posIni[0]*posIni[0]+posIni[1]*posIni[1]);
  return 1;
	}
	if(gRandom->Rndm()<ParticleAbsorptionProbability(params[4],params[0],params[2],params[3])){
		fProbe.fIsDecayed=kFALSE;
		fProbe.fIsAbsorbed=kTRUE;
		fProbe.fAbsorbtionRadius=TMath::Sqrt(posIni[0]*posIni[0]+posIni[1]*posIni[1]);
		return -1;
	}
	for (int j=3;j--;) posIni[j] = posCurr[j];
	fProbe.fIsAbsorbed = fProbe.fIsDecayed = 0;
	fProbe.fAbsorbtionRadius = fProbe.fDecayRadius = 0;
	return 0;
	//
}

//____________________________________________________
Bool_t FT2::PrepareProbe()
{
	// propagate the probe to max allowed R of the setup
	int nlrITS = fITS->GetNLayers(); // including passive layers
	//
	double xyzIni[3];
	if(fAllowDecay){fProbe.GetXYZ(xyzIni);}
	//
	for (int ilr=0;ilr<nlrITS;ilr++) {
		AliITSURecoLayer* lr = fITS->GetLayer(ilr);
		//
		if (lr->IsPassive()) { // cylindric passive layer, just go to this layer
			if (!PropagateToR(lr->GetRMax(),1,kFALSE,fSimMat,fSimMat)) return kFALSE;
		}
		else { // active layer, need to simulate the hit positions
			if (!PassActiveITSLayer(lr)) return kFALSE;
			fProbe.fProbeITSClusterIsCls[ilr]=1;
			Double_t posClsITS[3];
			fProbe.GetXYZ(posClsITS);
			fProbe.fProbeITSClusterX[ilr] = posClsITS[0];
			fProbe.fProbeITSClusterY[ilr] = posClsITS[1];
			fProbe.fProbeITSClusterZ[ilr] = posClsITS[2];
  }
		if(fAllowDecay && ProbeDecayAbsorb(xyzIni)) return kFALSE;
	}
	//
	fNTPCHits = 0;
	if (fIsTPC) {
		const double kTanSectH = 1.76326980708464975e-01; // tangent of  10 degree (sector half opening)
		if (!fProbe.RotateParamOnly( fProbe.PhiPos() )) {
			return kFALSE; // start in the frame with X pointing to track position
		}
		if (!PropagateToR(fTPCLayers[0].x-0.1,1,kFALSE,fSimMat,fSimMat)) {
			return kFALSE; // reach inner layer
		}
		if(fAllowDecay && ProbeDecayAbsorb(xyzIni)) return kFALSE;
		
		//
		int sector = -1;
		for (int ilrReset=0;ilrReset<(int)fTPCLayers.size();ilrReset++) {
			FT2TPCLayer_t &tpcLrReset = fTPCLayers[ilrReset];
			tpcLrReset.hitSect = -1;
		}
		for (int ilr=0;ilr<(int)fTPCLayers.size();ilr++) {
#if DEBUG>5
			printf("At TPC lr %d\n",ilr);
			fProbe.Print();
#endif
			
			FT2TPCLayer_t &tpcLr = fTPCLayers[ilr];
			
			double phi = fProbe.PhiPos();
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
			
			fProbe.fProbeZAtCutOffCheck = fProbe.GetZ();
			
			// cut taken from AliTPCclusterer, line 1179
			/*
			//Int_t crtime = Int_t((fParam->GetZLength(fSector)-fRecoParam->GetCtgRange()*fRx)/fZWidth-fParam->GetNTBinsL1()-5);
			//if (i%fMaxTime<=crtime) continue;
			Int_t fMaxTime = 1006;
			Float_t fZWidth = fTPCParam->GetZWidth();
			Float_t fRx			= fTPCParam->GetPadRowRadii(sector,ilr);
			Int_t crtime = Int_t((fTPCParam->GetZLength(sector)-fTPCRecoParam->GetCtgRange()*fRx)/fZWidth-fTPCParam->GetNTBinsL1()-5);
			cout << "This is the crtime " << crtime << endl;
			printf("fZWidth:  %f\n",fZWidth);
			printf("fRx:      %f\n",fRx);
			printf("ZLength:  %f\n",fTPCParam->GetZLength(sector));
			printf("CtgRange: %f\n",fTPCRecoParam->GetCtgRange());
			printf("NTBinsL1: %f\n",fTPCParam->GetNTBinsL1());
			printf("fMaxTime: %f\n",fMaxTime);
			printf("Cut?:     %d\n\n",ilr%fMaxTime);
			
	//		    if (ilr%fMaxTime<=crtime) continue;
	//		if (i%fMaxTime<=crtime) continue;
			*/
			if(TMath::Abs(fProbe.GetZ()/fProbe.GetX())>fTpcClusterAcc->Eval(fProbe.GetX())) continue;
			
			// cut during cluster seeding, taken from Ruben
			Double_t kDeltaZ = 10;
			if(TMath::Abs(fProbe.GetZ())>(/*AliTPCReconstructor::GetCtgRange()*/1.05*fProbe.GetX()+kDeltaZ)) continue ;
			//if (TMath::Abs(fProbe.GetZ())>kMaxZTPC) break; // exit from the TPC
// was here before
			
			if(fAllowDecay && ProbeDecayAbsorb(xyzIni)) return kFALSE;
			
			double maxY = kTanSectH*tpcLr.x - fTPCSectorEdge; // max allowed Y in the sector
			
			Double_t z = fProbe.GetZ();
			if(z>245.) z=245.;
			
			Double_t pointExB[3] = {TMath::Sqrt(fProbe.GetX()*fProbe.GetX()+fProbe.GetY()*fProbe.GetY())/250.,fProbe.GetZ()/250.,fProbe.Phi()-TMath::Pi()};
			Double_t TPCdistortionRPhi = fTPCDistortionRPhi->Eval(pointExB); // include tpc field distortion in rphi, assuming rphi==Y
			
			if((fProbe.GetY()<0 && TPCdistortionRPhi<0) || (fProbe.GetY()>0 && TPCdistortionRPhi>0)){
				if ((TMath::Abs(fProbe.GetY())+TMath::Abs(TPCdistortionRPhi))>maxY) {
#if DEBUG>3
					printf("No hit in dead zone because\n");
					printf("R: %f \t Z: %f \t phi: %f\n",pointExB[0],pointExB[1],pointExB[2]);
					printf("Probe.Y: %f with \n",fProbe.GetY());
					printf("Rphi distortion |%f| larger than pad maxY %0.4f\n",TPCdistortionRPhi,maxY);
					fProbe.Print();
#endif
					continue; // in dead zone
				}
			}
			else{
				if (TMath::Abs(fProbe.GetY()+TPCdistortionRPhi)>maxY) {
#if DEBUG>3
					printf("No hit in dead zone because\n");
					printf("R: %f \t Z: %f \t phi: %f\n",pointExB[0],pointExB[1],pointExB[2]);
					printf("Probe.Y: %f with \n",fProbe.GetY());
					printf("Rphi distortion |%f| larger than pad maxY %0.4f\n",TPCdistortionRPhi,maxY);
					fProbe.Print();
#endif
					continue; // in dead zone
				}
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
			cout << fNTPCHits << endl;
			//		fProbe.Print();
			#endif
			//
			Double_t TPCdistortionR = fTPCDistortionR->Eval(pointExB); // include tpc field distortion in R
			//
			if(TMath::Abs(TPCdistortionR)>tpcLr.pitch){		// distortion is larger than pad pitch => change in R
				Int_t ilrDis =ilr;
				if(TPCdistortionR<0){												// distortion to pad at lower radius
					ilrDis=ilr-1;
				}
				else if(TPCdistortionR>0){									// distortion to pad at higher radius
					ilrDis=ilr+1;
				}
				FT2TPCLayer_t &tpcLrDisR = fTPCLayers[ilrDis];
				//if (tpcLrDisR.isDead) continue;
				#if DEBUG>5
				printf("Reassigning cluster position to pad row because\n");
				printf("R distortion |%0.4f| larger than pad pitch %0.4f\n",TPCdistortionR,tpcLr.pitch);
				printf("Was there a hit before? %0.1f\n",tpcLrDisR.hitSect); // sector number set if hit before
				printf("Updating # of TPC hits? %i\n",fNTPCHits);
				#endif
				if(tpcLrDisR.hitSect!=-1){ // if hit before, do not increase TPC hits
					fTPCHitLr[fNTPCHits] = ilrDis;
				}
				else{	// if not hit before, increase TPC hits
					fTPCHitLr[fNTPCHits++] = ilrDis;
				}
				tpcLrDisR.hitSect = sector;
				tpcLrDisR.hitY = fProbe.GetY();
				tpcLrDisR.hitZ = fProbe.GetZ();
				//
				fProbe.fProbeTPCClusterIsCls[fNTPCHits]=1;
				Double_t posCls[3];
				fProbe.GetXYZ(posCls);
				fProbe.fProbeTPCClusterX[fNTPCHits] = posCls[0];
				fProbe.fProbeTPCClusterY[fNTPCHits] = posCls[1];
				fProbe.fProbeTPCClusterZ[fNTPCHits] = posCls[2];
				
				#if DEBUG>5
				printf("-----> %i\n",fNTPCHits);
				#endif
			}
			else{
				//
				//if (tpcLr.isDead) continue;
				//
				//printf("Track not distorted in R\n");
				if(tpcLr.hitSect!=-1){ // if hit before, do not increase TPC hits
					//printf("was hit before\n");
					fTPCHitLr[fNTPCHits] = ilr;
				}
				else{
					//printf("was not hit before\n");
					fTPCHitLr[fNTPCHits++] = ilr;
				}
				//printf("TPC hits: %i\n",fNTPCHits);
				tpcLr.hitSect = sector;
				tpcLr.hitY = fProbe.GetY();
				tpcLr.hitZ = fProbe.GetZ();
				//
				fProbe.fProbeTPCClusterIsCls[fNTPCHits]=1;
				Double_t posCls[3];
				fProbe.GetXYZ(posCls);
				fProbe.fProbeTPCClusterX[fNTPCHits] = posCls[0];
				fProbe.fProbeTPCClusterY[fNTPCHits] = posCls[1];
				fProbe.fProbeTPCClusterZ[fNTPCHits] = posCls[2];
			}
		}
	}
	// RS: temporary fix: the errors assigned in clusterizer is sqrt(pixel_extent/12)
	double eta = fProbe.Eta();
	fITSerrSclZ = 20e-4*(1.4+0.61*eta*eta)*TMath::Sqrt(1./12.)/fSigZITS;
	//
	return kTRUE;
}


//________________________________________________________________________
Bool_t FT2::ReconstructProbe(TParticle* part)
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
	fProbe.fAbsPdgCodeForTracking = 211; // default assumption is pion
	fProbe.fProbeMass = TDatabasePDG::Instance()->GetParticle(fProbe.fAbsPdgCodeForTracking)->Mass();
	//
	// needed for parameterizations
	Double_t point[3] = {
		1./AliMathBase::BetheBlochAleph(0.0001+part->P()/TDatabasePDG::Instance()->GetParticle(fProbe.fAbsPdgCode)->Mass()),
		abs(part->Eta()),
		fTrueMCtrackMult};
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
				Int_t typePID=-1;
				Double_t p = fProbe.P(); fProbe.fTPCmomentum = p;
		  if(fProbe.fAbsPdgCode==11){
				typePID=0;
			}
			else if(fProbe.fAbsPdgCode==13){
				typePID=1;
			}
			else if(fProbe.fAbsPdgCode==211){
				typePID=2;
			}
			else if(fProbe.fAbsPdgCode==321){
				typePID=3;
			}
			else if(fProbe.fAbsPdgCode==2212){
				typePID=4;
			}
			else {AliFatal("PDC Code cannot be treated in the code - This case should not be possible here!");}
#if DEBUG>5
				AliInfo(Form("PDG code of Probe: %4.f",fProbe.fAbsPdgCode));
				AliInfo(Form("Momentum of Probe: %f",fProbe.fTPCmomentum));
#endif
				Double_t signalTPC = fTpcPidSignal[typePID]->Eval(point);
				
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
			//
			// cluster errror parameterization
			Int_t typeROC;
			Double_t angleProbe = fProbe.GetSnp()*fProbe.GetSnp();
			angleProbe = TMath::Sqrt(TMath::Abs(angleProbe/(1.-angleProbe)));
			//
			Double_t xyzForDrift[3];
			fProbe.GetXYZ(xyzForDrift);
			Double_t zDriftLength = TMath::Abs(TMath::Abs(xyzForDrift[2])-250.);
			//
			// TPC cluter pickup probability (1/dEdx,eta,multiplicity)
			if(ih<fTPCInnerRows){
				// IROC
				typeROC = 0;
				if(gRandom->Rndm()>fTPCClsLossProbIROC->Eval(point)){continue;}
#if DEBUG>5
				AliInfo(Form("TPC Cluster Loss Probability IROC: %f for Pdg %f - 1/dEdx: %f - #eta: %f - mult: %f",fTPCClsLossProbIROC->Eval(point),fProbe.fAbsPdgCode,point[0],point[1],point[2]));
#endif
			}
			else if(ih>=fTPCInnerRows && ih<(fTPCInnerRows+fTPCMiddleRows)){
				// OROC (medium)
				typeROC = 1;
				if(gRandom->Rndm()>fTPCClsLossProbOROCmedium->Eval(point)){continue;}
#if DEBUG>5
				AliInfo(Form("TPC Cluster Loss Probability OROC (m): %f for Pdg %f - 1/dEdx: %f - #eta: %f - mult: %f",fTPCClsLossProbOROCmedium->Eval(point),fProbe.fAbsPdgCode,point[0],point[1],point[2]));
#endif
			}
			else if(ih>=(fTPCInnerRows+fTPCMiddleRows)){
				// OROC (long)
				typeROC = 2;
				if(gRandom->Rndm()>fTPCClsLossProbOROClong->Eval(point)){continue;}
#if DEBUG>5
				AliInfo(Form("TPC Cluster Loss Probability OROC (l): %f for Pdg %f - 1/dEdx: %f - #eta: %f - mult: %f",fTPCClsLossProbOROClong->Eval(point),fProbe.fAbsPdgCode,point[0],point[1],point[2]));
#endif
			}
			//
			fTPCClusterErrorParamY->SetParameters(typeROC,angleProbe); // (y fixed; TPC ROC; angle)
			Double_t tpcSigY = fTPCClusterErrorParamY->Eval(zDriftLength);
			//
			fTPCClusterErrorParamZ->SetParameters(typeROC,angleProbe); // (z fixed; TPC ROC; angle)
			Double_t tpcSigZ = fTPCClusterErrorParamZ->Eval(zDriftLength);
			//
#if DEBUG>5
			printf("TPC Cluster error parameterization\n");
			printf("Layer?\t\t\t %i\n",ih);
			printf("ROC?\t\t\t %i\n",typeROC);
			printf("Angle?\t\t\t %f\n",angleProbe);
			printf("L_drift?\t\t %f\n",zDriftLength);
			printf("sig_y --> %f\n",tpcSigY);
			printf("sig_z --> %f\n",tpcSigZ);
			cout << "------------------------------------------" << endl;
#endif
			double chi = UpdateKalman(&fProbe,tpcLr.hitY,tpcLr.hitZ,tpcSigY,tpcSigZ,kTRUE);
			if (chi<0) return kFALSE;
			fChi2TPC += chi;
			fNClTPC++;
			
			if(ih==0){ // set inner track parameters
				fProbe.fInnerTrackParameters[0] = fProbe.GetAlpha();
				fProbe.fInnerTrackParameters[1] = fProbe.GetX();
				fProbe.fInnerTrackParameters[2] = fProbe.GetY();
				fProbe.fInnerTrackParameters[3] = fProbe.GetZ();
				fProbe.fInnerTrackParameters[4] = fProbe.GetSnp();
				fProbe.fInnerTrackParameters[5] = fProbe.GetTgl();
				fProbe.fInnerTrackParameters[6] = fProbe.GetSigned1Pt();
			}
			
		}
		// go to ITS/TPC matching R, accounting for TGeo materials
		if (!PropagateToR(fITS->GetRITSTPCRef(),-1,kTRUE, kFALSE, kTRUE)) return kFALSE;
	}
	//
	//
#if DEBUG
	printf("Outward kalman results:\n");
	for (int ilr=0;ilr<fITS->GetNLayersActive();ilr++) {
		if (GetKalmanOut(ilr).GetSigmaY2()<=0) continue;
		printf("Extrap to Lr %d : ",ilr); GetKalmanOut(ilr).Print();
	}
#endif
	//
	// ITS tracking
	if (fUseKalmanOut) MakeITSKalmanOut();
	int maxMiss = fITS->GetNLayersActive() - fMinITSLrHit;
	int nMiss = 0;
	for (int ilr=fITS->GetNLayersActive();ilr--;) {
		int res = ReconstructOnITSLayer(ilr);
		if (res<1) nMiss++;
		if (res<0) break;
		if (nMiss>maxMiss) return kFALSE;
	} // true hit is present
	//
	return fNClTPC>=fMinTPCHits && fNClITS>=fMinITSLrHit;
}

//________________________________________________________________________
Double_t FT2::UpdateKalman(AliExternalTrackParam* trc, double y,double z,double sigY,double sigZ,Bool_t randomize,double sclY,double sclZ)
{
	// calman update, return chi2 increment or -1 on failure
	//
	if (randomize) {
		double ry,rz;
		gRandom->Rannor(ry,rz);
		y += sigY*ry;
		z += sigZ*rz;
	}
	sigY *= sclY;
	sigZ *= sclZ;
	double meas[2] = {y,z};
	double measErr2[3] = {sigY*sigY,0,sigZ*sigZ};
	double chi2 = trc->GetPredictedChi2(meas,measErr2);
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
		Bool_t res = propErr ? fProbe.PropagateTo(x,fBz) : fProbe.PropagateParamOnlyTo(x,fBz);
		if (!res)  {
#if DEBUG
			printf("Failed Prop PropagateToX(%.1f,%d,%d,%d,%d)",xTgt,dir,propErr,simMat,useTGeo);
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
					printf("Failed ApplyMSEloss(%f,%f) in PropagateToX(%.1f,%d,%d,%d,%d)",
						   xx0,xrho,xTgt,dir,propErr,simMat,useTGeo);
					fProbe.Print();
#endif
					return kFALSE;
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
Bool_t FT2::GetRoadWidth(AliITSURecoLayer* lrA,double *pr,Int_t nstd)
{
	static AliExternalTrackParam sc;   // seed copy for manipulations
#if DEBUG>5
	printf(">>RW at Lr%d\n",lrA->GetActiveID());
#endif
	sc = fProbe;
	double xt;
	if (!sc.GetXatLabR(lrA->GetRMax(),xt,fBz,0/*1*/)) return kFALSE;
	Bool_t res = nstd>0 ? sc.PropagateTo(xt,fBz) : sc.PropagateParamOnlyTo(xt,fBz);
	if (!res) return kFALSE;
	sc.GetXYZ(&pr[AliITSUTrackerGlo::kTrXIn]);
	pr[AliITSUTrackerGlo::kTrPhiIn] = TMath::ATan2(pr[AliITSUTrackerGlo::kTrYIn],pr[AliITSUTrackerGlo::kTrXIn]);
	res = nstd>0 ? sc.Rotate(pr[AliITSUTrackerGlo::kTrPhiIn]) : sc.RotateParamOnly(pr[AliITSUTrackerGlo::kTrPhiIn]);
	if (!res) return kFALSE; // go to the frame of the entry point into the layer
	BringTo02Pi(pr[AliITSUTrackerGlo::kTrPhiIn]);
	double dr  = lrA->GetDR();                              // approximate X dist at the inner radius
	if (!sc.GetXYZAt(sc.GetX()-dr, fBz, pr + AliITSUTrackerGlo::kTrXOut)) {
		// special case: track does not reach inner radius, might be tangential
		double r = sc.GetD(0,0,fBz);
		double x;
		if (!sc.GetXatLabR(r,x,fBz,0/*-1*/)) return kFALSE;
		dr = Abs(sc.GetX() - x);
		if (!sc.GetXYZAt(x, fBz, pr + AliITSUTrackerGlo::kTrXOut)) return kFALSE;
	}
	//
	pr[AliITSUTrackerGlo::kTrPhiOut] = ATan2(pr[AliITSUTrackerGlo::kTrYOut],pr[AliITSUTrackerGlo::kTrXOut]);
	BringTo02Pi(pr[AliITSUTrackerGlo::kTrPhiOut]);
	double sgy = 1e-6; // dummy spread
	double sgz = 1e-6; // dummy spread
	if (nstd>0) {
		sgy = nstd*TMath::Sqrt(sc.GetSigmaY2());
		sgz = nstd*TMath::Sqrt(sc.GetSigmaZ2());
	}
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
	AliITSUGeomTGeo* gm = fITS->GetGeom();
	//
	int lrAID = lr->GetActiveID();
	AliITSURecoSens **hitSens = &fITSSensCand[lrAID][0];
	//
	if (!GetRoadWidth(lr, trImpData, -1)) return kFALSE;
	fNITSHits[lrAID] = 0;
	fNITSSensCand[lrAID] = 0;
	//
	const AliITSsegmentation* segm = gm->GetSegmentation(lrAID);
	int nsens = fNITSSensCand[lrAID] = lr->FindSensors(&trImpData[AliITSUTrackerGlo::kTrPhi0], hitSens);
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
		Bool_t res = fProbe.RotateParamOnly(sens->GetPhiTF());
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
		if (!fNITSHits[lrAID]) fNITSLrHit++;
		fNITSHits[lrAID]++;
#if DEBUG>5
		printf("Register hit %d at lr%d\n",fNITSHits[lrAID],lr->GetActiveID());
		fProbe.Print();
#endif
	}
	//
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
	if (!hfake) hfake = new TH2F("hfake","hf",100,1,-1,100,1,-1);
	hfake->Fill(dy,dz);
	//
#if DEBUG>1
	printf("AddFake DY:%f DZ:%f Chi2Max:%f Chi2F:%f\n",dy,dz, d2max,d2);
#endif
	yz[0] = extyz[0]+dy;
	yz[1] = extyz[1]+dz;
	return kTRUE;
}

//_____________________________________________
Int_t FT2::ReconstructOnITSLayer(int ilr, double chi2Cut)
{
	// find track prolongation on its layer ilr
	//
	fCurrITSLr = ilr;
	AliITSURecoLayer* lrA = fITS->GetLayerActive(ilr);
#if DEBUG>5
	AliInfoF("Entering ITS cluster loop on lr%d: %d clusters",ilr,fNITSHits[ilr]);
#endif
	AliITSURecoSens* sens = 0;
	int iht = -1;
	if (fNITSHits[ilr]) {  // for >1 hits (overlap) chose only 1st ons
		iht = 0;
		sens = fITSSensHit[ilr][iht];
	}
	else { // there was no hit but fake may be picked
		Double_t trImpData[AliITSUTrackerGlo::kNTrImpData];
		if (!GetRoadWidth(lrA, trImpData)) return -1; // track cannot be continued
		AliITSURecoSens **hitSens = &fITSSensCand[ilr][0];
		int nsens = lrA->FindSensors(&trImpData[AliITSUTrackerGlo::kTrPhi0], hitSens);
		if (nsens) sens = hitSens[0];
		else return 0; // no candidate -> no update on this layer
	}
	if (!fProbe.Rotate(sens->GetPhiTF())) {
#if DEBUG
		printf("Failed to rotate to sensor phi %f at lr %d\n",sens->GetPhiTF(),ilr);
		fProbe.Print();
		return -1; //  track cannot be continued
#endif
	}
	//
	if (!PropagateToX(sens->GetXTF(),-1,kTRUE, kFALSE, kTRUE)) return -1; // account for materials
	//
	double trCov[3],trPos[2]; // get inward/outward smoothed position
	const AliExternalTrackParam* trSmooth = GetSmoothedEstimate(ilr,&fProbe,trPos,trCov);
	if (!trSmooth) return -1;
	//  if (!GetSmoothedEstimate(ilr,&fProbe,trPos,trCov)) return -1;
	double xyz[3];
	fProbe.GetXYZ(xyz);
	double rho = HitDensity(xyz[0]*xyz[0]+xyz[1]*xyz[1],fProbe.GetTgl());  // get hit density at given r, eta
	//
	double chi2Best = chi2Cut;
	double sgy = fSigYITS*fITSerrSclY;
	double sgz = fSigZITS*fITSerrSclZ;
	double yzw[2],yzf[2],yzcov[3]={sgy*sgy,0,sgz*sgz};
	Int_t winOK = -1; // -1: no winner, 0: fake winner, 1: correct winner
	if (iht>=0) {
		double ry,rz;
		gRandom->Rannor(ry,rz);
		yzw[0] = fITSHitYZ[ilr][iht][0] + fSigYITS*ry;
		yzw[1] = fITSHitYZ[ilr][iht][1] + fSigZITS*rz;
		//    double chi2 = fProbe.GetPredictedChi2(yzw,yzcov);
		double chi2 = trSmooth->GetPredictedChi2(yzw,yzcov);
		if (chi2<chi2Cut) {
			chi2Best = chi2;
			winOK = 1;
			yzw[0] -= fSigYITS*ry; // restore true positions,
			yzw[1] -= fSigZITS*rz; // randomization will be done in the update
		}
	}
	//
	//  double nstd = TMath::Sqrt(chi2Cut);
	double dYsearch = 2*TMath::Max(0.00001,TMath::Sqrt((trCov[0]+yzcov[0])*chi2Cut));
	double dZsearch = 2*TMath::Max(0.00001,TMath::Sqrt((trCov[2]+yzcov[2])*chi2Cut));
	int nFakeCandRnd = gRandom->Poisson(rho*dYsearch*dZsearch); // expected number of surrounding hits
	int nFakeCandCorr = 0;
	if (fAllocCorrelatedITSFakes) nFakeCandCorr = gRandom->Poisson(GetNCorrelITSFakes(ilr));
	int nFakeCandTot = nFakeCandCorr + nFakeCandRnd;
	//
#if DEBUG>1
	printf("Lr:%d NF:%3d NFC:%3d rho:%f S:%.3fx%.3f\n",ilr,nFakeCandRnd,nFakeCandCorr,rho,dYsearch,dZsearch);
#endif
	if (nFakeCandTot) {
		AliITSUGeomTGeo* gm = fITS->GetGeom();
		const AliITSsegmentation* segm = gm->GetSegmentation(ilr);
		for (int ifc=nFakeCandTot;ifc--;) {
			if (ifc<nFakeCandRnd) { // flat component
				yzf[0] = trPos[0] + (gRandom->Rndm()-0.5)*dYsearch;
				yzf[1] = trPos[1] + (gRandom->Rndm()-0.5)*dZsearch;
			}
			else { // correlated component
				double sy,sz;
				gRandom->Rannor(sy,sz);
				if (iht>=0) { // there was a true hit, put the accompanying hits around
					yzf[0] = fITSHitYZ[ilr][iht][0] + sy*GetCorrelITSFakesSigY(ilr);
					yzf[1] = fITSHitYZ[ilr][iht][1] + sz*GetCorrelITSFakesSigZ(ilr);
				}
				else {
					yzf[0] = trPos[0] + sy*GetCorrelITSFakesSigY(ilr);
					yzf[1] = trPos[1] + sz*GetCorrelITSFakesSigZ(ilr);
				}
			}
			//      double chi2 = fProbe.GetPredictedChi2(yzf,yzcov);
			double chi2 = trSmooth->GetPredictedChi2(yzf,yzcov);
			if (chi2>chi2Best) continue;
			//
			// check if the hit is in valid area
			double xyzL[3],xyzT[3] = {fProbe.GetX(),yzf[0],yzf[1]};
			gm->GetMatrixT2L(sens->GetID())->LocalToMaster(xyzT,xyzL);
			int ix,iz;
			if (!segm->LocalToDet(xyzL[0],xyzL[2],ix,iz)) continue;
			//
			// fake wins
			yzw[0] = yzf[0];
			yzw[1] = yzf[1];
			winOK = 0;
			chi2Best = chi2;
		}
	}
	//
	if (winOK<0) return 0; // no update on this layer
	double chiw = UpdateKalman(&fProbe,yzw[0],yzw[1],fSigYITS,fSigZITS,kTRUE,fITSerrSclY,fITSerrSclZ);
	fProbe.chiwITS[ilr]=chiw;
#if DEBUG>5
	AliInfoF("Updated at ITS Lr%d, chi2:%f NclITS:%d Fakes: %d",ilr,chiw,fNClITS,fNClITSFakes);
#endif
	if (chiw>0) {
		fChi2ITS += chiw;
		fNClITS++;
		fITSPattern |= 0x1<<ilr;
		if (!winOK) { // register fake
			fITSPatternFake |= 0x1<<ilr;
			fNClITSFakes++;
		}
		return 1;
	}
	return 0;
}


//_____________________________________________
const AliExternalTrackParam* FT2::GetSmoothedEstimate(int ilr,const AliExternalTrackParam* trcInw, double* trPos,double* trCov)
{
	static AliExternalTrackParam trJoint;
	if (fUseKalmanOut && fKalmanOutward[ilr].GetSigmaY2()>0) {
  //
		trJoint = fKalmanOutward[ilr];
		if (TMath::Abs(trJoint.GetAlpha()-trcInw->GetAlpha())>1e-6 && !trJoint.Rotate(trcInw->GetAlpha())) {
#if DEBUG
			AliInfoF("Failed to rotate outward Kalman track to alpha=%f",trcInw->GetAlpha());
			trJoint.Print();
			return 0;
#endif
		}
		if (TMath::Abs(trJoint.GetX()-trcInw->GetX())>1e-4 && !trJoint.PropagateTo(trcInw->GetX(),fBz)) {
#if DEBUG
			AliInfoF("Failed to propagate outward Kalman track to X=%f",trcInw->GetX());
			trJoint.Print();
			return 0;
#endif
		}
		double inwPos[2] = {trcInw->GetY(),trcInw->GetZ()}; // inward track extrapolation
		double inwCov[3] = {trcInw->GetSigmaY2(),trcInw->GetSigmaZY(),trcInw->GetSigmaZ2()};
		if (!trJoint.Update(inwPos,inwCov)) { // updated it with outward one
#if DEBUG
			AliInfo("Failed to update outward Kalman track by inward one");
			trJoint.Print();
			trcInw->Print();
			return 0;
#endif
		}
	}
	else {
		trJoint = *trcInw;
	}
#if DEBUG//>5
	printf("Errors on Lr %d: \n",ilr);
	printf("Inward  :  {%+e,%+e,%+e} {%+e,%+e}\n",trcInw->GetSigmaY2(),trcInw->GetSigmaZY(),
		   trcInw->GetSigmaZ2(), trcInw->GetY(),trcInw->GetZ());
	if (fUseKalmanOut && fKalmanOutward[ilr].GetSigmaY2()>0) {
		printf("Outward :  {%+e,%+e,%+e} {%+e,%+e}\n",fKalmanOutward[ilr].GetSigmaY2(),
			   fKalmanOutward[ilr].GetSigmaZY(),fKalmanOutward[ilr].GetSigmaZ2(),
			   fKalmanOutward[ilr].GetY(),fKalmanOutward[ilr].GetZ());
	}
	printf("Weighted:  {%+e,%+e,%+e} {%+e,%+e}\n",trJoint.GetSigmaY2(),trJoint.GetSigmaZY(),trJoint.GetSigmaZ2(),
		   trJoint.GetY(),trJoint.GetZ());
#endif
	trPos[0] = trJoint.GetY();
	trPos[1] = trJoint.GetZ();
	trCov[0] = trJoint.GetSigmaY2();
	trCov[1] = trJoint.GetSigmaZY();
	trCov[2] = trJoint.GetSigmaZ2();
	//
	return &trJoint;
}

//_____________________________________________
Bool_t FT2::MakeITSKalmanOut()
{
	// prepare kalman outward track
	// save probe kinematics
	//
	int lr0 = 0;
	for (lr0=0;lr0<fITS->GetNLayersActive();lr0++) {
		((double*)fKalmanOutward[lr0].GetCovariance())[0] = -1; // invalidate outward track
		if (fNITSHits[lr0]<1) continue;
		break;
	}
	if (lr0>=fITS->GetNLayersActive()) return kFALSE; // 1st hit layer
	//
	AliExternalTrackParam prbSav = fProbe;
	ResetCovMat(&fProbeIni);
	fProbe.AliExternalTrackParam::operator=(fProbeIni);
	//
	AliITSURecoSens* sens = 0;
	int nclAdd = 0;
	for (int ilA=lr0;ilA<fITS->GetNLayersActive();ilA++) {
		AliITSURecoLayer* lrA = fITS->GetLayerActive(ilA);
		if (fNITSHits[ilA]<1 && nclAdd) { // no hits
			if (!PropagateToR(lrA->GetRMax(),1,kTRUE,kFALSE,fSimMat)) {
				fProbe.AliExternalTrackParam::operator=(prbSav); return kFALSE;
			}
			fKalmanOutward[ilA] = fProbe;
			continue;
		}
		sens = fITSSensHit[ilA][0];
		if (!fProbe.Rotate(sens->GetPhiTF()) || !PropagateToX(sens->GetXTF(),1,kTRUE, kFALSE, kTRUE)) {
			fProbe.AliExternalTrackParam::operator=(prbSav); return kFALSE;
		}
		fKalmanOutward[ilA] = fProbe;
		//
		double chi2l = UpdateKalman(&fProbe,fITSHitYZ[ilA][0][0],fITSHitYZ[ilA][0][1],fSigYITS,fSigZITS,kTRUE,fITSerrSclY,fITSerrSclZ);
		if (chi2l<0) {fProbe.AliExternalTrackParam::operator=(prbSav); return kFALSE;}
		nclAdd++;
	}
	//
	fProbe.AliExternalTrackParam::operator=(prbSav);  // restore
	return kTRUE;
}

//_________________________________________________________
Double_t FT2::GetITSRMin() const
{
	return fITS->GetRMin();
}
//_________________________________________________________
Double_t FT2::ParticleDecayProbability(Double_t step){
	
	TParticlePDG* pdgp = TDatabasePDG::Instance()->GetParticle(fProbe.fAbsPdgCode);
	
	Double_t mass		= fProbe.fTrueMass;
	Double_t energy		= TMath::Sqrt(mass*mass+fProbe.GetP()*fProbe.GetP());
	Double_t lifetime	= pdgp->Lifetime();
	Double_t sol		= 3.E10; // cm/s
	
	if(lifetime==0){return -1;}
	else{
		return (1.-TMath::Exp(-step/((energy/mass)*sol*lifetime)));
	}
}
//_________________________________________________________
Double_t FT2::ParticleAbsorptionProbability(Double_t length,Double_t rho, Double_t A, Double_t Z){
	
	Double_t Navo = 6.022E23; // Avogadro
	
	Double_t sigma0 = 0.;
	Int_t pdg = fProbe.fPdgCode;
	Double_t mom = fProbe.P();
	
	
	if(pdg==+211){		sigma0 = fXSectionHp[0]->GetBinContent(fXSectionHp[0]->FindBin(mom));}
	else if(pdg==-211){	sigma0 = fXSectionHp[1]->GetBinContent(fXSectionHp[1]->FindBin(mom));}
	else if(pdg==+321){	sigma0 = fXSectionHp[2]->GetBinContent(fXSectionHp[2]->FindBin(mom));}
	else if(pdg==-321){	sigma0 = fXSectionHp[3]->GetBinContent(fXSectionHp[3]->FindBin(mom));}
	else if(pdg==+2212){sigma0 = fXSectionHp[4]->GetBinContent(fXSectionHp[4]->FindBin(mom));}
	else if(pdg==-2212){sigma0 = fXSectionHp[5]->GetBinContent(fXSectionHp[5]->FindBin(mom));}
	else if(pdg==+11){
		Double_t mass		= fProbe.fTrueMass;
		//Double_t energy		= TMath::Sqrt(mass*mass+mom*mom);
		//Double_t gamma		= energy/mass;
		Double_t radLength	= (1432.8*A)/(Z*(Z+1)*(11.319-TMath::Log(Z))); // g/cm2
		return (1.-TMath::Exp(-length*rho/radLength));
	}
	else if(pdg==-11){
		Double_t mass		= fProbe.fTrueMass;
		Double_t energy		= TMath::Sqrt(mass*mass+mom*mom);
		Double_t gamma		= energy/mass;
		Double_t radLength	= (1432.8*A)/(Z*(Z+1)*(11.319-TMath::Log(Z))); // g/cm2
		Double_t radiusEl	= 2.8179403267E-13; // cm
		Double_t sigmaAni	= Z*TMath::Pi()*radiusEl*radiusEl/(gamma+1)*((gamma*gamma+4.*gamma+1)/(gamma*gamma-1)*TMath::Log(gamma+TMath::Sqrt(gamma*gamma-1))-(gamma+3)/TMath::Sqrt(gamma*gamma-1));
		
		Double_t lambdaAni = A/(Navo*sigmaAni);
		
		return (1.-TMath::Exp(-length*rho/radLength))*(1-TMath::Exp(-length*rho/lambdaAni));
	}
	else return -1;
	sigma0*=1E-27; // X-Section from mb to cm2
	
	Double_t lambda = TMath::Power(A,1./3.)/(Navo*sigma0);
#if DEBUG>5
	AliInfo(Form("\n### %i with p = %f\n### Rho: %f\n### A: %f\n### Z: %f\n### Navo: %E\n### sigma0: %E\n### x: %E\n### Lambda: %f\n### xrho/La: %f\n",pdg,mom,rho,A,Z,Navo,sigma0,length,lambda,length*rho/lambda));
#endif
	return (1.-TMath::Exp(-length*rho/lambda));
	
}
//_________________________________________________________
Int_t FT2::CutOnTrackToClusterChi2ITS(){
	
	float gloChi2 = 0, penalty=0;
	int ncl=0;
	
	Int_t accept = 8;
	for (int ilr=7;ilr--;) {
		if (fProbe.chiwITS[ilr]<0){penalty+=fC0missPen[ilr];}
		else {
			if (fProbe.chiwITS[ilr]>fC0tr2clChi2[ilr]) {accept = ilr; break;}
			gloChi2 += fProbe.chiwITS[ilr];
			ncl++;
		}
		int ndf = TMath::Max(2*ncl-5,1);
		float gloChi2norm = gloChi2/ndf + penalty;
		if (gloChi2norm>fC0gloChi2[ilr]) {accept = ilr; break;}
	}
	return accept;
}
