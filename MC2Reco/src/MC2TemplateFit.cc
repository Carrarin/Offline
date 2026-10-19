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

  namespace { constexpr double kNs = 5.0; }   // digitizer period

  MC2TemplateFit::MC2TemplateFit(const Config& cfg) : cfg_(cfg) {
    std::unique_ptr<TFile> f(TFile::Open(cfg_.templateFile.c_str(), "READ"));
    if (!f || f->IsZombie()) throw std::runtime_error("MC2TemplateFit: cannot open " + cfg_.templateFile);
    if (cfg_.shapeMode > 0) {
      auto* g0 = dynamic_cast<TGraph*>(f->Get("T0")); auto* g1 = dynamic_cast<TGraph*>(f->Get("T1")); auto* g2 = dynamic_cast<TGraph*>(f->Get("T2"));
      auto* a0 = dynamic_cast<TParameter<double>*>(f->Get("A0"));
      if (!g0 || !g1 || !g2 || !a0) throw std::runtime_error("MC2TemplateFit: no T0, T1, T2, A0 (amplitude model) in " + cfg_.templateFile);
      a0_ = a0->GetVal(); const double norm = build(g0); buildShape(g1, g2, norm);
      return;
    }
    auto* g = dynamic_cast<TGraph*>(f->Get(cfg_.graphName.c_str()));
    if (!g || g->GetN() < 10) throw std::runtime_error("MC2TemplateFit: no graph " + cfg_.graphName + " in " + cfg_.templateFile);
    build(g);
  }

  MC2TemplateFit::MC2TemplateFit(const Config& cfg, const TGraph& g) : cfg_(cfg) {
    if (g.GetN() < 10) throw std::runtime_error("MC2TemplateFit: template graph with fewer than 10 points");
    cfg_.shapeMode = 0;   // the amplitude model is only read from a file
    TGraph c(g); build(&c);
  }

  std::vector<double> MC2TemplateFit::table(TGraph* g, double x0, double step, int n) {
    g->Sort(); TSpline3 sp("sp3", g); std::vector<double> t(n);
    for (int i = 0; i < n; ++i) t[i] = sp.Eval(x0 + i * step);
    return t;
  }

  void MC2TemplateFit::buildShape(TGraph* g1, TGraph* g2, double norm) {
    tab1_ = table(g1, x0_, step_, tab_.size()); tab2_ = table(g2, x0_, step_, tab_.size());
    for (auto& v : tab1_) v /= norm;
    for (auto& v : tab2_) v /= norm;
  }

  double MC2TemplateFit::build(TGraph* g) {
    // cubic spline through the points (as caloreco: TSpline3 of gtempl), tabulated every 0.05 ns: the production files have
    // 1 ns points and a linear interpolation between them is too coarse on the rising edge
    g->Sort();
    x0_ = g->GetX()[0]; step_ = 0.05;
    const int n = int(std::floor((g->GetX()[g->GetN() - 1] - x0_) / step_)) + 1;
    tab_ = table(g, x0_, step_, n);
    // x_peak = the highest point of the input graph (not of the spline): the fit window is defined from t0 + x_peak, and
    // the window end of the MC2 studies (-31 ns) was optimised with this definition
    double yg = -1e30; for (int i = 0; i < g->GetN(); ++i) if (g->GetY()[i] > yg) { yg = g->GetY()[i]; xpeak_ = g->GetX()[i]; }
    double ym = 0;
    for (int i = 0; i < n; ++i) if (tab_[i] > ym) ym = tab_[i];
    if (ym <= 0) throw std::runtime_error("MC2TemplateFit: template without a positive maximum");
    for (auto& v : tab_) v /= ym;
    return ym;
  }

  double MC2TemplateFit::look(const std::vector<double>& t, double x) const {
    const double u = (x - x0_) / step_;
    if (u <= 0 || u >= t.size() - 1) return 0.;
    const int i = int(u);
    return t[i] + (u - i) * (t[i + 1] - t[i]);
  }

  double MC2TemplateFit::eval(double x) const {
    if (cfg_.shapeMode != 1 && cfg_.shapeMode != 3) return look(tab_, x);
    return look(tab_, x) + u_ * look(tab1_, x) + u_ * u_ * look(tab2_, x);
  }

  // chi2 at fixed t0 with A, B solved linearly; B clamped to +-maxBaseline
  double MC2TemplateFit::chi2At(const std::vector<double>& y, int i0, int i1, double t0, double& A, double& B) const {
    if (cfg_.shapeMode >= 2) {   // y = A T + C T1 + B (T = T0 for 2, T(u) for 3): normal equations of (T, T1, 1), B clamped as below
      double m[3][3] = {}, v[3] = {};
      for (int i = i0; i < i1; ++i) {
        if (use_ && !(*use_)[i]) continue;
        const double f[3] = { eval(i * kNs - t0), look(tab1_, i * kNs - t0), 1. };
        for (int a = 0; a < 3; ++a) { v[a] += f[a] * y[i]; for (int b = 0; b < 3; ++b) m[a][b] += f[a] * f[b]; }
      }
      m[1][1] += lam_;   // shapeMode 3: penalty lam C^2 (0 for 2)
      const double d = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
      if (m[0][0] <= 0 || std::fabs(d) < 1e-12) { A = B = c_ = 0; return 1e30; }
      auto det3 = [&](int k) { double q[3][3]; for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) q[a][b] = b == k ? v[a] : m[a][b];
        return q[0][0] * (q[1][1] * q[2][2] - q[1][2] * q[2][1]) - q[0][1] * (q[1][0] * q[2][2] - q[1][2] * q[2][0]) + q[0][2] * (q[1][0] * q[2][1] - q[1][1] * q[2][0]); };
      A = det3(0) / d; c_ = det3(1) / d; B = det3(2) / d;
      if (B > cfg_.maxBaseline || B < -cfg_.maxBaseline) {   // B fixed at the limit, A and C from y - B
        B = B > 0 ? cfg_.maxBaseline : -cfg_.maxBaseline;
        const double w0 = v[0] - B * m[0][2], w1 = v[1] - B * m[1][2], d2 = m[0][0] * m[1][1] - m[0][1] * m[1][0];
        if (std::fabs(d2) < 1e-12) { A = c_ = 0; return 1e30; }
        A = (w0 * m[1][1] - m[0][1] * w1) / d2; c_ = (m[0][0] * w1 - m[1][0] * w0) / d2;
      }
      double c2 = 0;
      for (int i = i0; i < i1; ++i) { if (use_ && !(*use_)[i]) continue; const double r = y[i] - A * eval(i * kNs - t0) - c_ * look(tab1_, i * kNs - t0) - B; c2 += r * r; }
      return c2 + lam_ * c_ * c_;
    }
    double s1 = 0, ss = 0, sss = 0, sy = 0, ssy = 0;
    for (int i = i0; i < i1; ++i) { if (use_ && !(*use_)[i]) continue; const double s = eval(i * kNs - t0); s1 += 1; ss += s; sss += s * s; sy += y[i]; ssy += s * y[i]; }
    if (sss <= 0) { A = B = 0; return 1e30; }
    const double det = s1 * sss - ss * ss;
    if (det <= 0) { A = B = 0; return 1e30; }
    A = (s1 * ssy - ss * sy) / det; B = (sss * sy - ss * ssy) / det;
    if (B > cfg_.maxBaseline || B < -cfg_.maxBaseline) { B = B > 0 ? cfg_.maxBaseline : -cfg_.maxBaseline; A = (ssy - B * ss) / sss; }
    double c2 = 0;
    for (int i = i0; i < i1; ++i) { if (use_ && !(*use_)[i]) continue; const double r = y[i] - A * eval(i * kNs - t0) - B; c2 += r * r; }
    return c2;
  }

  double MC2TemplateFit::golden(const std::vector<double>& y, int i0, int i1, double a, double b) const {
    const double gr = 0.6180339887; double A, B;
    double x1 = b - gr * (b - a), x2 = a + gr * (b - a), f1 = chi2At(y, i0, i1, x1, A, B), f2 = chi2At(y, i0, i1, x2, A, B);
    for (int it = 0; it < 40 && b - a > 1e-4; ++it) {
      if (f1 < f2) { b = x2; x2 = x1; f2 = f1; x1 = b - gr * (b - a); f1 = chi2At(y, i0, i1, x1, A, B); }
      else         { a = x1; x1 = x2; f1 = f2; x2 = a + gr * (b - a); f2 = chi2At(y, i0, i1, x2, A, B); }
    }
    return 0.5 * (a + b);
  }

  MC2TemplateFit::Result MC2TemplateFit::fit(const std::vector<int>& wf) const {
    Result r;
    if (int(wf.size()) <= cfg_.nPedSamples + 4) return r;
    double ped = 0; for (int q = 0; q < cfg_.nPedSamples; ++q) ped += wf[q]; ped /= cfg_.nPedSamples;
    std::vector<double> y(wf.size()); std::vector<char> use(wf.size(), 1); int nSat = 0;
    for (size_t q = 0; q < wf.size(); ++q) { y[q] = wf[q] - ped; if (wf[q] >= cfg_.satADC) { use[q] = 0; ++nSat; } }
    r = nSat ? fit(y, use) : fit(y); r.pedestal = ped; r.nSat = nSat; return r;
  }

  MC2TemplateFit::Result MC2TemplateFit::fit(const std::vector<double>& y, const std::vector<char>& use) const {
    use_ = use.size() == y.size() ? &use : nullptr;
    struct Reset { const std::vector<char>*& p; ~Reset() { p = nullptr; } } reset{use_};
    if (cfg_.shapeMode != 1 && cfg_.shapeMode != 3) return fitOnce(y);
    // shape tied to the amplitude: u from the largest sample, then twice from the fitted amplitude; u kept inside the bands of the
    // model (300 - 2000 ADC: below, the quadratic term is an extrapolation; above, every pulse is saturated)
    auto setU = [&](double a) { u_ = std::clamp(std::log(std::max(a, 1.) / a0_), std::log(300. / a0_), std::log(2000. / a0_)); };
    // shapeMode 3: penalty weight noise^2 / (shapeSigma A)^2 with A from the previous iteration (the largest sample first)
    auto setLam = [&](double a) { if (cfg_.shapeMode == 3) { const double s = cfg_.shapeSigma * std::max(a, 1.); lam_ = cfg_.noise * cfg_.noise / (s * s); } };
    const double a0 = y.empty() ? 0. : *std::max_element(y.begin(), y.end());
    setU(a0); setLam(a0);
    Result r = fitOnce(y);
    for (int it = 0; it < 2 && r.ok; ++it) { setU(r.amplitude); setLam(r.amplitude); r = fitOnce(y); }
    r.shape = cfg_.shapeMode == 3 && r.ok ? u_ + r.shape : u_;
    return r;
  }

  MC2TemplateFit::Result MC2TemplateFit::fitOnce(const std::vector<double>& y) const {
    Result r;
    const int n = y.size(); if (n < 5) return r;
    const int im = int(std::max_element(y.begin(), y.end()) - y.begin()); r.iMax = im;
    const int iFirst = cfg_.maxPre < 0 ? 0 : std::max(0, im - cfg_.maxPre);   // first sample allowed in the fit
    // 1. window from the largest sample, coarse scan and refinement
    int i0 = std::max(iFirst, im + int(std::lround(-70. / kNs))), i1 = std::min(n, im + int(std::lround(-20. / kNs)) + 1);
    if (i1 - i0 < 5) return r;
    const double tpk = im * kNs; double best = 1e30, tb = 0, A = 0, B = 0;
    for (double t0 = tpk - xpeak_ - 8; t0 <= tpk - xpeak_ + 8.01; t0 += 1.) { const double c = chi2At(y, i0, i1, t0, A, B); if (c < best) { best = c; tb = t0; } }
    double t0 = golden(y, i0, i1, tb - 1, tb + 1); double c2 = chi2At(y, i0, i1, t0, A, B);
    if (A <= 0) return r;
    // 2. window from the fitted peak, two iterations
    for (int it = 0; it < 2; ++it) {
      const double tp = t0 + xpeak_;
      i0 = std::max(iFirst, int(std::ceil((tp + cfg_.fitLo) / kNs))); i1 = std::min(n, int(std::floor((tp + cfg_.fitHi) / kNs)) + 1);
      if (i1 - i0 < 4) return r;
      t0 = golden(y, i0, i1, t0 - 1.5, t0 + 1.5); c2 = chi2At(y, i0, i1, t0, A, B);
      if (A <= 0) return r;
    }
    int nUsed = 0; for (int i = i0; i < i1; ++i) if (!use_ || (*use_)[i]) ++nUsed;
    r.ok = true; r.time = t0 + xpeak_; r.amplitude = A; r.baseline = B;
    if (cfg_.shapeMode >= 2) r.shape = c_ / A;
    r.ndf = std::max(1, nUsed - (cfg_.shapeMode == 2 ? 4 : 3)); r.chi2ndf = c2 / (cfg_.noise * cfg_.noise) / r.ndf;
    return r;
  }

}
