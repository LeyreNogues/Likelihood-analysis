void take_plots_statusdisplay(){

  const TString stringStatusFluxlc   = "./Status_Output_fluxlc.root";
  //const TString stringStatusFluxlc   = "/nfs/magic-ifae01/afernandez/CrabNebula/flute/cycleX/Status_flute.root";

TFile* fileStatusFluxlc      = TFile::Open(stringStatusFluxlc,"READ");
MStatusArray* StatusOutputFluxlc = fileStatusFluxlc->Get("MStatusDisplay");

TCanvas * MigMatrix = (TCanvas *) StatusOutputFluxlc->FindCanvas("Migration Matrix");

TObject *obj;
TIter next(MigMatrix->GetListOfPrimitives());
TPad *pad;

TH2D* Mig_Matrix = new TH2D();

while ((obj = next())) {
  cout << "Reading: " << obj->GetName() << endl;
  if (obj->InheritsFrom("TH2D")) {
    cout << "TH2D: " << obj->GetName() << endl;
    // obj->Print("all");
    Mig_Matrix = (TH2D*) obj;
    
  }
 }
 TCanvas * c2 = new TCanvas();
 c2->SetLogx();
 c2->SetLogy();
 // Mig_Matrix->Print("all");
 Mig_Matrix->DrawClone("colz");
 
}
