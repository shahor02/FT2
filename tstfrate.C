{
  hft = new TH1F("hft","",7,-0.6,6.5);
  hftn = new TH1F("hftn","",7,-0.6,6.5);
  hf0 = new TH1F("hf0","",7,-0.6,6.5);
  hf0n = new TH1F("hf0n","",7,-0.6,6.5);
  //
  TCut ctgen = "nclTPC>60&&nclITS>4";
  for (int i=0;i<7;i++) tt->Draw(Form("%d >>+ hft",i),ctgen*Form("(clITSF&(0x1<<%d))>0",i));
  for (int i=0;i<7;i++) tt->Draw(Form("%d >>+ hftn",i),ctgen*Form("(clITS&(0x1<<%d))>0",i));
  hft->Sumw2();
  hft->Divide(hftn);

  for (int i=0;i<7;i++) t0->Draw(Form("%d >>+ hf0",i),ctgen*Form("(clITSF&(0x1<<%d))>0",i));
  for (int i=0;i<7;i++) t0->Draw(Form("%d >>+ hf0n",i),ctgen*Form("(clITS&(0x1<<%d))>0",i));
  hf0->Sumw2();
  hf0->Divide(hf0n);
   //
}
