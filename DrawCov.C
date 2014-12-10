#include <TMath.h>
#include <TPolyLine.h>
#include <TVector2.h>
#include <TH2F.h>

TH2* hh;
TPolyLine* pl;

Bool_t DiagonalizeErrors(const double* cov, double &sy2d, double &sz2d);

void DrawCov(double sig00, double sig01,double sig11)
{
  double cov[3] = {sig00,sig01,sig11};
  double sA2,sB2;
  if (!DiagonalizeErrors(cov,sA2,sB2)) return;
  sA2 = TMath::Sqrt(sA2);
  sB2 = TMath::Sqrt(sB2);
  double tht = 0.5*TMath::ATan2(2*sig01,(sig00-sig11));
  printf("Axes: %f %f, tht: %f\n",sA2,sB2,tht);
  TVector2 vs;
  int np=100;
  double dph = TMath::TwoPi()/np;
  pl = new TPolyLine(np);
  //  
  double axX = TMath::Sqrt(sig00);
  double axY = TMath::Sqrt(sig11);
  hh = new TH2F("err","",20,-2*axX,2*axX, 20,-2*axY,2*axY);
  hh->Draw();
  for (int i=0;i<np;i++) {
    vs.Set(sA2*TMath::Cos(dph*i), sB2*TMath::Sin(dph*i));
    vs = vs.Rotate(tht);
    pl->SetPoint(i,vs.X(),vs.Y());
  }
  pl->Draw();
}

//__________________________________________________________
Bool_t DiagonalizeErrors(const double* cov, double &sy2d, double &sz2d) 
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
