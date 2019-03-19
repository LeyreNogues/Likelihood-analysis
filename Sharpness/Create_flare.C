void Create_flare(Double_t a, Double_t b, Double_t c, Double_t d, Double_t e){
  
  Double_t A=a;
  Double_t B=b;
  Double_t C=c;
  Double_t D=d;
  Double_t E=e;
  
  TF2* f = new TF2("f", "[0]+[1]*exp(-0.5*((x-[2]-[4]*y)/[3])**2)", 0., 100., 0.1, 10.);
  f->SetParameters(A, B, C, D, E);
  FILE* fout;
  fout = fopen(Form("simflare_%1.1f_%1.1f_%1.1f_%1.1f_%1.1f.dat",a , b, c, d, e),"wb");
  for(int i = 0; i<1e3; i++)
    {
    double ax, ay;
    f->GetRandom2(ax,ay);
    fprintf(fout, "%lf %lf \n", ax, ay);
    }
}
