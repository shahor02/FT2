#include <TMath.h>
#include <RooEllipse.h>
#include <TH2F.h>
#include <TPad.h>
#include <TMarker.h>

TH2* hh=0;
RooEllipse* elp=0;
TMarker* pl=0;

void DrawCov(double *posCov, double *errCov, double* posHit=0,Bool_t elpOnly=kFALSE)
{
  double axX = TMath::Sqrt(errCov[0]);
  double axY = TMath::Sqrt(errCov[2]);
  double corr = errCov[1]/(axX*axY);
  //
  double xmn=posCov[0]-2*axX;
  double xmx=posCov[0]+2*axX;
  double ymn=posCov[1]-2*axY;
  double ymx=posCov[1]+2*axY;
  if (posHit) {
    if (xmn>posHit[0]) xmn = posHit[0]-0.1*(xmx-posHit[0]);
    else if (xmx<posHit[0]) xmx = posHit[0]+0.1*(posHit[0]-xmn);
    if (ymn>posHit[1]) ymn = posHit[1]-0.1*(ymx-posHit[1]);
    else if (ymx<posHit[1]) ymx = posHit[1]+0.1*(posHit[1]-ymn);
  }
  if (!elpOnly) {
    hh = new TH2F("err","",20,xmn,xmx,20,ymn,ymx);
    hh->Draw();
  }
  elp = new RooEllipse("elp",posCov[0],posCov[1],axX,axY,corr);
  elp->Draw();
  gPad->SetGrid(1,1);
  //
  if (posHit) {
    pl = new TMarker(posHit[0],posHit[1],20);
    pl->Draw();
  }
}
