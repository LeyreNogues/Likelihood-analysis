void Read()
{

 MEnergyEst         *Eest=0;
 MTime              *Time;
 MHadronness        *Hadron=0;
 MStereoPar         *Theta=0;

//TString path = "/home/leyre/MAGIC.../*.root";

//TFile* outputfile = new TFile("List");
TChain* filechain = new TChain("Events");  //Tree
filechain->Add("*Crab*.root");

filechain->SetBranchStatus("*", 0); //TBrowser
filechain->SetBranchStatus("MEnergyEst.fEnergy", 1);
filechain->SetBranchStatus("MTime_1.*", 1);
filechain->SetBranchStatus("MHadronness.fHadronness", 1); 
filechain->SetBranchStatus("MStereoPar.fTheta2", 1); 

filechain->SetBranchAddress("MEnergyEst.", &Eest);
filechain->SetBranchAddress("MTime_1.", &Time);
filechain->SetBranchAddress("MHadronness.", &Hadron);
filechain->SetBranchAddress("MStereoPar.", &Theta);

Float_t nEvents;
FILE *fout;

fout = fopen("Listafinal_Crab.dat","wb");

nEvents  = filechain->GetEntries();
 fprintf(fout, "Hadronness        Time(MJD)          Energy(GeV) \n");

for (Int_t event = 0; event < nEvents; event++)
  {
    filechain->GetEvent(event);
    Double_t Energy = Eest->GetEnergy();
    Double_t  Tiempo = Time->GetMjd();
    Double_t Hadronn = Hadron->GetHadronness();
    Double_t Teta = Theta->GetTheta2();
    Byte_t h, m, s;
    UShort_t ms;
    Time->GetTime(h,m,s,ms);
    // if( Energy > 0)cout << Energy << "     " << Form("%0.8f",Tiempo) << "     "  << (Int_t)h << ":"  << (Int_t)m << ":"  << (Int_t)s << "."  << ms << "     "  << Form("%0.3f",Hadronn) << "     " << Form("%0.3f",Teta) << endl;
   

     if( Energy > 0)
    {
      fprintf(fout,"  %0.3lf        %0.8lf         %0.3lf     \n", Hadronn, Tiempo, Energy);
    }
   }
}
