{
//=========Macro generated from canvas: c1/c1
//=========  (Tue Feb 23 15:26:24 2016) by ROOT version5.34/14
   TCanvas *c1 = new TCanvas("c1", "c1",292,120,700,502);
//   c1->SetHighLightColor(2);
//   c1->Range(-0.08109689,1.164741,5.567437,5.842741);
//   c1->SetFillColor(0);
//   c1->SetBorderMode(0);
//   c1->SetBorderSize(2);
   c1->SetLogx();
   c1->SetLogy();
   c1->SetGridx();
   c1->SetGridy();
//   c1->SetLeftMargin(0.13);
//   c1->SetRightMargin(0.15);
//   c1->SetBottomMargin(0.13);
//   c1->SetFrameBorderMode(0);
//   c1->SetFrameBorderMode(0);
   
  TGraph* gre = new TGraph(28);
   gre->SetName("Graph_from_mga1");
   gre->SetTitle("Effective Collection Area Before Cuts");

   gre->SetMarkerStyle(22);
   gre->SetMarkerSize(0.7);
   gre->SetPoint(0,5.973739,0);
   gre->SetPoint(1,8.300483,0);
   gre->SetPoint(2,11.53348,3.24533);
   gre->SetPoint(3,16.02572,20.51407);
   gre->SetPoint(4,22.26767,123.8702);
   gre->SetPoint(5,30.94083,761.136);
   gre->SetPoint(6,42.99214,3727.596);
   gre->SetPoint(7,59.73739,13194.41);
   gre->SetPoint(8,83.00483,32261.97);
   gre->SetPoint(9,115.3348,52874.8);
   gre->SetPoint(10,160.2572,68198);
   gre->SetPoint(11,222.6767,77421.93);
   gre->SetPoint(12,309.4083,85609.01);
   gre->SetPoint(13,429.9214,94620.26);
   gre->SetPoint(14,597.3739,102467.7);
   gre->SetPoint(15,830.0483,105256.5);
   gre->SetPoint(16,1153.348,106260.6);
   gre->SetPoint(17,1602.572,107973);
   gre->SetPoint(18,2226.767,107017.8);
   gre->SetPoint(19,3094.083,120304.3);
   gre->SetPoint(20,4299.214,140517.1);
   gre->SetPoint(21,5973.739,0);
   gre->SetPoint(22,8300.483,0);
   gre->SetPoint(23,11533.48,0);
   gre->SetPoint(24,16025.72,0);
   gre->SetPoint(25,22267.67,0);
   gre->SetPoint(26,30940.83,0);
   gre->SetPoint(27,42992.14,0);
   
   gre->Draw();


}
