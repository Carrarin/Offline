//
// MC2TemplateFit: template fit of one calorimeter waveform (see the description in MC2TemplateFit.hh).
//
#include "Offline/MC2Reco/inc/MC2TemplateFit.hh"

#include "TFile.h"
#include "TGraph.h"
#include "TParameter.h"
#include "TSpline.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>

namespace mu2e {

  namespace {

    constexpr double kSamplePeriod = 5.0;   // ns, digitizer period

    // Determinant of a 3x3 matrix (expansion along the first row).
    double determinant3(const double m[3][3]) {
      return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
           - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
           + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    }

    // Determinant of the matrix m with its column k replaced by v (Cramer's rule).
    double determinantWithColumn(const double m[3][3], const double v[3], int k) {
      double q[3][3];
      for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
          q[a][b] = (b == k) ? v[a] : m[a][b];
        }
      }
      return determinant3(q);
    }

  }

  //--------------------------------------------------------------------------------------------------------------------
  // Construction: read the template and tabulate it
  //--------------------------------------------------------------------------------------------------------------------

  MC2TemplateFit::MC2TemplateFit(const Config& cfg) : cfg_(cfg) {
    std::unique_ptr<TFile> file(TFile::Open(cfg_.templateFile.c_str(), "READ"));
    if (!file || file->IsZombie()) {
      throw std::runtime_error("MC2TemplateFit: cannot open " + cfg_.templateFile);
    }

    if (cfg_.shapeMode > 0) {
      // Amplitude model: T0, T1, T2 and the reference amplitude A0
      auto* graph0 = dynamic_cast<TGraph*>(file->Get("T0"));
      auto* graph1 = dynamic_cast<TGraph*>(file->Get("T1"));
      auto* graph2 = dynamic_cast<TGraph*>(file->Get("T2"));
      auto* referenceAmplitude = dynamic_cast<TParameter<double>*>(file->Get("A0"));
      if (!graph0 || !graph1 || !graph2 || !referenceAmplitude) {
        throw std::runtime_error("MC2TemplateFit: no T0, T1, T2, A0 (amplitude model) in " + cfg_.templateFile);
      }
      a0_ = referenceAmplitude->GetVal();
      const double normalisation = buildTemplate(graph0);
      buildShapeTerms(graph1, graph2, normalisation);
      return;
    }

    auto* graph = dynamic_cast<TGraph*>(file->Get(cfg_.graphName.c_str()));
    if (!graph || graph->GetN() < 10) {
      throw std::runtime_error("MC2TemplateFit: no graph " + cfg_.graphName + " in " + cfg_.templateFile);
    }
    buildTemplate(graph);
  }

  MC2TemplateFit::MC2TemplateFit(const Config& cfg, const TGraph& graph) : cfg_(cfg) {
    if (graph.GetN() < 10) {
      throw std::runtime_error("MC2TemplateFit: template graph with fewer than 10 points");
    }
    cfg_.shapeMode = 0;   // the amplitude model is only read from a file
    TGraph copy(graph);
    buildTemplate(&copy);
  }

  // Cubic spline through the points of the graph, evaluated at nPoints points every step from xStart.
  std::vector<double> MC2TemplateFit::tabulate(TGraph* graph, double xStart, double step, int nPoints) {
    graph->Sort();
    TSpline3 spline("sp3", graph);
    std::vector<double> table(nPoints);
    for (int i = 0; i < nPoints; ++i) {
      table[i] = spline.Eval(xStart + i * step);
    }
    return table;
  }

  // Tables of the amplitude terms T1 and T2, on the same points as T0 and divided by the normalisation of T0.
  void MC2TemplateFit::buildShapeTerms(TGraph* graph1, TGraph* graph2, double normalisation) {
    table1_ = tabulate(graph1, xStart_, step_, table0_.size());
    table2_ = tabulate(graph2, xStart_, step_, table0_.size());
    for (auto& value : table1_) {
      value /= normalisation;
    }
    for (auto& value : table2_) {
      value /= normalisation;
    }
  }

  // Table of the template, normalised to 1 at its maximum; returns the normalisation.
  // The cubic spline (as caloreco: TSpline3 of gtempl) is tabulated every 0.05 ns: the production files have 1 ns
  // points and a linear interpolation between them is too coarse on the rising edge.
  double MC2TemplateFit::buildTemplate(TGraph* graph) {
    graph->Sort();
    xStart_ = graph->GetX()[0];
    step_ = 0.05;
    const double xEnd = graph->GetX()[graph->GetN() - 1];
    const int nPoints = int(std::floor((xEnd - xStart_) / step_)) + 1;
    table0_ = tabulate(graph, xStart_, step_, nPoints);

    // x_peak = the highest point of the input graph (not of the spline): the fit window is defined from t0 + x_peak,
    // and the window end of the MC2 studies (-31 ns) was optimised with this definition
    double highestY = -1e30;
    for (int i = 0; i < graph->GetN(); ++i) {
      if (graph->GetY()[i] > highestY) {
        highestY = graph->GetY()[i];
        xPeak_ = graph->GetX()[i];
      }
    }

    double maximum = 0;
    for (int i = 0; i < nPoints; ++i) {
      if (table0_[i] > maximum) {
        maximum = table0_[i];
      }
    }
    if (maximum <= 0) {
      throw std::runtime_error("MC2TemplateFit: template without a positive maximum");
    }
    for (auto& value : table0_) {
      value /= maximum;
    }
    return maximum;
  }

  // Linear interpolation in a table; 0 outside it.
  double MC2TemplateFit::interpolate(const std::vector<double>& table, double x) const {
    const double position = (x - xStart_) / step_;
    if (position <= 0 || position >= table.size() - 1) {
      return 0.;
    }
    const int i = int(position);
    return table[i] + (position - i) * (table[i + 1] - table[i]);
  }

  // Template at x: T0, or T0 + u T1 + u^2 T2 when the shape is tied to the amplitude (shapeMode 1 and 3).
  double MC2TemplateFit::templateAt(double x) const {
    if (cfg_.shapeMode != 1 && cfg_.shapeMode != 3) {
      return interpolate(table0_, x);
    }
    return interpolate(table0_, x) + u_ * interpolate(table1_, x) + u_ * u_ * interpolate(table2_, x);
  }

  //--------------------------------------------------------------------------------------------------------------------
  // chi2 at fixed t0
  //--------------------------------------------------------------------------------------------------------------------

  // Residual sum of squares at fixed t0 on the samples [iFirst, iEnd), with A and B solved linearly; B is clamped to
  // +-maxBaseline (then A is solved again with B fixed).
  double MC2TemplateFit::chi2AtTime(const std::vector<double>& y, int iFirst, int iEnd, double t0,
                                    double& A, double& B) const {
    if (cfg_.shapeMode >= 2) {
      return chi2WithShapeTerm(y, iFirst, iEnd, t0, A, B);
    }

    // Sums for the normal equations of y = A s + B
    double n = 0, sumS = 0, sumSS = 0, sumY = 0, sumSY = 0;
    for (int i = iFirst; i < iEnd; ++i) {
      if (!isUsed(i)) {
        continue;
      }
      const double s = templateAt(i * kSamplePeriod - t0);
      n += 1;
      sumS += s;
      sumSS += s * s;
      sumY += y[i];
      sumSY += s * y[i];
    }
    if (sumSS <= 0) {
      A = B = 0;
      return 1e30;
    }
    const double det = n * sumSS - sumS * sumS;
    if (det <= 0) {
      A = B = 0;
      return 1e30;
    }
    A = (n * sumSY - sumS * sumY) / det;
    B = (sumSS * sumY - sumS * sumSY) / det;

    if (B > cfg_.maxBaseline || B < -cfg_.maxBaseline) {
      B = B > 0 ? cfg_.maxBaseline : -cfg_.maxBaseline;
      A = (sumSY - B * sumS) / sumSS;
    }

    double chi2 = 0;
    for (int i = iFirst; i < iEnd; ++i) {
      if (!isUsed(i)) {
        continue;
      }
      const double residual = y[i] - A * templateAt(i * kSamplePeriod - t0) - B;
      chi2 += residual * residual;
    }
    return chi2;
  }

  // shapeMode 2 and 3: y = A T + C T1 + B (T = T0 for 2, T(u) for 3), normal equations of (T, T1, 1) solved with
  // Cramer's rule; shapeMode 3 adds the penalty penaltyWeight * C^2. B is clamped as in chi2AtTime.
  double MC2TemplateFit::chi2WithShapeTerm(const std::vector<double>& y, int iFirst, int iEnd, double t0,
                                           double& A, double& B) const {
    double m[3][3] = {};
    double v[3] = {};
    for (int i = iFirst; i < iEnd; ++i) {
      if (!isUsed(i)) {
        continue;
      }
      const double x = i * kSamplePeriod - t0;
      const double f[3] = { templateAt(x), interpolate(table1_, x), 1. };
      for (int a = 0; a < 3; ++a) {
        v[a] += f[a] * y[i];
        for (int b = 0; b < 3; ++b) {
          m[a][b] += f[a] * f[b];
        }
      }
    }
    m[1][1] += penaltyWeight_;   // 0 for shapeMode 2

    const double det = determinant3(m);
    if (m[0][0] <= 0 || std::fabs(det) < 1e-12) {
      A = B = shapeCoefficient_ = 0;
      return 1e30;
    }
    A = determinantWithColumn(m, v, 0) / det;
    shapeCoefficient_ = determinantWithColumn(m, v, 1) / det;
    B = determinantWithColumn(m, v, 2) / det;

    if (B > cfg_.maxBaseline || B < -cfg_.maxBaseline) {
      // B fixed at the limit, A and C solved again from y - B
      B = B > 0 ? cfg_.maxBaseline : -cfg_.maxBaseline;
      const double w0 = v[0] - B * m[0][2];
      const double w1 = v[1] - B * m[1][2];
      const double det2 = m[0][0] * m[1][1] - m[0][1] * m[1][0];
      if (std::fabs(det2) < 1e-12) {
        A = shapeCoefficient_ = 0;
        return 1e30;
      }
      A = (w0 * m[1][1] - m[0][1] * w1) / det2;
      shapeCoefficient_ = (m[0][0] * w1 - m[1][0] * w0) / det2;
    }

    double chi2 = 0;
    for (int i = iFirst; i < iEnd; ++i) {
      if (!isUsed(i)) {
        continue;
      }
      const double x = i * kSamplePeriod - t0;
      const double residual = y[i] - A * templateAt(x) - shapeCoefficient_ * interpolate(table1_, x) - B;
      chi2 += residual * residual;
    }
    return chi2 + penaltyWeight_ * shapeCoefficient_ * shapeCoefficient_;
  }

  // t0 of the minimum chi2 in [tLow, tHigh], golden-section search (at most 40 steps, down to 1e-4 ns).
  double MC2TemplateFit::goldenSectionMinimum(const std::vector<double>& y, int iFirst, int iEnd,
                                              double tLow, double tHigh) const {
    const double goldenRatio = 0.6180339887;
    double A, B;
    double t1 = tHigh - goldenRatio * (tHigh - tLow);
    double t2 = tLow + goldenRatio * (tHigh - tLow);
    double chi2At1 = chi2AtTime(y, iFirst, iEnd, t1, A, B);
    double chi2At2 = chi2AtTime(y, iFirst, iEnd, t2, A, B);
    for (int step = 0; step < 40 && tHigh - tLow > 1e-4; ++step) {
      if (chi2At1 < chi2At2) {
        tHigh = t2;
        t2 = t1;
        chi2At2 = chi2At1;
        t1 = tHigh - goldenRatio * (tHigh - tLow);
        chi2At1 = chi2AtTime(y, iFirst, iEnd, t1, A, B);
      } else {
        tLow = t1;
        t1 = t2;
        chi2At1 = chi2At2;
        t2 = tLow + goldenRatio * (tHigh - tLow);
        chi2At2 = chi2AtTime(y, iFirst, iEnd, t2, A, B);
      }
    }
    return 0.5 * (tLow + tHigh);
  }

  //--------------------------------------------------------------------------------------------------------------------
  // Fit
  //--------------------------------------------------------------------------------------------------------------------

  MC2TemplateFit::Result MC2TemplateFit::fit(const std::vector<int>& waveform) const {
    Result result;
    if (int(waveform.size()) <= cfg_.nPedSamples + 4) {
      return result;
    }

    // Pedestal: mean of the first samples
    double pedestal = 0;
    for (int i = 0; i < cfg_.nPedSamples; ++i) {
      pedestal += waveform[i];
    }
    pedestal /= cfg_.nPedSamples;

    // Pedestal-subtracted samples; saturated samples are left out
    std::vector<double> y(waveform.size());
    std::vector<char> use(waveform.size(), 1);
    int nSaturated = 0;
    for (size_t i = 0; i < waveform.size(); ++i) {
      y[i] = waveform[i] - pedestal;
      if (waveform[i] >= cfg_.satADC) {
        use[i] = 0;
        ++nSaturated;
      }
    }

    result = nSaturated ? fit(y, use) : fit(y);
    result.pedestal = pedestal;
    result.nSat = nSaturated;
    return result;
  }

  MC2TemplateFit::Result MC2TemplateFit::fit(const std::vector<double>& y, const std::vector<char>& use) const {
    // The sample mask is used during this fit only
    use_ = use.size() == y.size() ? &use : nullptr;
    struct ClearMaskAtExit {
      const std::vector<char>*& mask;
      ~ClearMaskAtExit() { mask = nullptr; }
    } clearMask{use_};

    if (cfg_.shapeMode != 1 && cfg_.shapeMode != 3) {
      return fitWithFixedShape(y);
    }

    // Shape tied to the amplitude: u from the largest sample, then twice from the fitted amplitude. u is kept inside
    // the range of the model (300 - 2000 ADC: below, the quadratic term is an extrapolation; above, every pulse is
    // saturated).
    auto setShapeFromAmplitude = [&](double amplitude) {
      u_ = std::clamp(std::log(std::max(amplitude, 1.) / a0_), std::log(300. / a0_), std::log(2000. / a0_));
    };
    // shapeMode 3: penalty weight noise^2 / (shapeSigma A)^2, with A from the previous iteration
    // (the largest sample first)
    auto setPenaltyFromAmplitude = [&](double amplitude) {
      if (cfg_.shapeMode == 3) {
        const double sigma = cfg_.shapeSigma * std::max(amplitude, 1.);
        penaltyWeight_ = cfg_.noise * cfg_.noise / (sigma * sigma);
      }
    };

    const double largestSample = y.empty() ? 0. : *std::max_element(y.begin(), y.end());
    setShapeFromAmplitude(largestSample);
    setPenaltyFromAmplitude(largestSample);
    Result result = fitWithFixedShape(y);
    for (int iteration = 0; iteration < 2 && result.ok; ++iteration) {
      setShapeFromAmplitude(result.amplitude);
      setPenaltyFromAmplitude(result.amplitude);
      result = fitWithFixedShape(y);
    }
    result.shape = cfg_.shapeMode == 3 && result.ok ? u_ + result.shape : u_;
    return result;
  }

  // One fit with the current shape: steps 3 and 4 of the description in the header.
  MC2TemplateFit::Result MC2TemplateFit::fitWithFixedShape(const std::vector<double>& y) const {
    Result result;
    const int nSamples = y.size();
    if (nSamples < 5) {
      return result;
    }
    const int iMax = int(std::max_element(y.begin(), y.end()) - y.begin());
    result.iMax = iMax;
    const int iAllowed = cfg_.maxPre < 0 ? 0 : std::max(0, iMax - cfg_.maxPre);   // first sample allowed in the fit

    // Step 3: window from the largest sample, coarse scan of t0 in 1 ns steps, then refinement
    int iFirst = std::max(iAllowed, iMax + int(std::lround(-70. / kSamplePeriod)));
    int iEnd = std::min(nSamples, iMax + int(std::lround(-20. / kSamplePeriod)) + 1);
    if (iEnd - iFirst < 5) {
      return result;
    }
    const double tMax = iMax * kSamplePeriod;
    double bestChi2 = 1e30;
    double bestT0 = 0;
    double A = 0;
    double B = 0;
    for (double t0 = tMax - xPeak_ - 8; t0 <= tMax - xPeak_ + 8.01; t0 += 1.) {
      const double chi2 = chi2AtTime(y, iFirst, iEnd, t0, A, B);
      if (chi2 < bestChi2) {
        bestChi2 = chi2;
        bestT0 = t0;
      }
    }
    double t0 = goldenSectionMinimum(y, iFirst, iEnd, bestT0 - 1, bestT0 + 1);
    double chi2 = chi2AtTime(y, iFirst, iEnd, t0, A, B);
    if (A <= 0) {
      return result;
    }

    // Step 4: window [peak + fitLo, peak + fitHi] from the fitted peak, two iterations
    for (int iteration = 0; iteration < 2; ++iteration) {
      const double tPeak = t0 + xPeak_;
      iFirst = std::max(iAllowed, int(std::ceil((tPeak + cfg_.fitLo) / kSamplePeriod)));
      iEnd = std::min(nSamples, int(std::floor((tPeak + cfg_.fitHi) / kSamplePeriod)) + 1);
      if (iEnd - iFirst < 4) {
        return result;
      }
      t0 = goldenSectionMinimum(y, iFirst, iEnd, t0 - 1.5, t0 + 1.5);
      chi2 = chi2AtTime(y, iFirst, iEnd, t0, A, B);
      if (A <= 0) {
        return result;
      }
    }

    int nUsed = 0;
    for (int i = iFirst; i < iEnd; ++i) {
      if (isUsed(i)) {
        ++nUsed;
      }
    }
    result.ok = true;
    result.time = t0 + xPeak_;
    result.amplitude = A;
    result.baseline = B;
    if (cfg_.shapeMode >= 2) {
      result.shape = shapeCoefficient_ / A;
    }
    const int nFreeParameters = (cfg_.shapeMode == 2) ? 4 : 3;
    result.ndf = std::max(1, nUsed - nFreeParameters);
    result.chi2ndf = chi2 / (cfg_.noise * cfg_.noise) / result.ndf;
    return result;
  }

}
