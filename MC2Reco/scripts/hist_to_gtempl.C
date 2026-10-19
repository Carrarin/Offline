// Template histogram of laser_template_stream.C (templates.root, hLaserTemplate, 0.5 ns bins) -> file with TGraph "gtempl"
// at the bin centres (no resampling), maximum = 1, peak shifted to x = 0: the input of MC2TemplateFit / MC2CaloRecoDigiMaker.
void hist_to_gtempl(TString in, TString out, TString hname = "hLaserTemplate") {
  TFile fi(in); TH1* h = (TH1*)fi.Get(hname); if (!h) { printf("ERROR: no %s in %s\n", hname.Data(), in.Data()); return; }
  const double xp = h->GetBinCenter(h->GetMaximumBin()), m = h->GetMaximum();
  TGraphErrors* g = new TGraphErrors(); g->SetName("gtempl");
  for (int i = 1; i <= h->GetNbinsX(); ++i) g->AddPoint(h->GetBinCenter(i) - xp, h->GetBinContent(i) / m);
  TFile fo(out, "RECREATE"); g->Write("gtempl"); fo.Close();
  printf("written %s: %d points, step %.2f ns, x %.2f..%.2f ns\n", out.Data(), g->GetN(), h->GetBinWidth(1), g->GetX()[0], g->GetX()[g->GetN() - 1]); }
