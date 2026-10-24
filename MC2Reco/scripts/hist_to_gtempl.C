// hist_to_gtempl.C: laser pulse template histogram -> TGraph "gtempl", the template file of MC2TemplateFit /
// MC2CaloRecoDigiMaker.
//
// The input is the template histogram of laser_template_stream.C (templates.root, hLaserTemplate, 0.5 ns bins). The
// graph has one point per bin at the bin centre (no resampling), normalised to maximum = 1 and with the peak shifted
// to x = 0.
//
// Usage (any ROOT):
//
//   root -l -b -q 'Offline/MC2Reco/scripts/hist_to_gtempl.C("templates.root","gtempl.root")'
//   root -l -b -q 'Offline/MC2Reco/scripts/hist_to_gtempl.C("templates.root","gtempl.root","hLaserTemplate")'
//
// Arguments: in = file with the histogram; out = file to write (RECREATE); hname = histogram name (default
//            hLaserTemplate).
//
// Output: <out> with the TGraphErrors "gtempl" (errors 0). Prints the number of points, the step and the x range.

#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TString.h>

#include <cstdio>

void hist_to_gtempl(TString in, TString out, TString hname = "hLaserTemplate") {
  TFile inputFile(in);
  TH1* templateHist = (TH1*)inputFile.Get(hname);
  if (!templateHist) {
    printf("ERROR: no %s in %s\n", hname.Data(), in.Data());
    return;
  }

  // peak position and height of the histogram
  const double peakX = templateHist->GetBinCenter(templateHist->GetMaximumBin());
  const double peakHeight = templateHist->GetMaximum();

  // one point per bin, peak at x = 0, maximum = 1
  TGraphErrors* graph = new TGraphErrors();
  graph->SetName("gtempl");
  for (int i = 1; i <= templateHist->GetNbinsX(); ++i) {
    graph->AddPoint(templateHist->GetBinCenter(i) - peakX, templateHist->GetBinContent(i) / peakHeight);
  }

  TFile outputFile(out, "RECREATE");
  graph->Write("gtempl");
  outputFile.Close();

  printf("written %s: %d points, step %.2f ns, x %.2f..%.2f ns\n", out.Data(), graph->GetN(),
         templateHist->GetBinWidth(1), graph->GetX()[0], graph->GetX()[graph->GetN() - 1]);
}
