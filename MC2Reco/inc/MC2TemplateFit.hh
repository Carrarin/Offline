#ifndef MC2Reco_MC2TemplateFit_hh
#define MC2Reco_MC2TemplateFit_hh
//
// Template fit of one calorimeter waveform, the method of the MC2 timing studies (mc2-ana scripts/laser_template_range_opt.C,
// function fitPeak, pedestal mode "clamp"), without art and without Minuit:
//
//   model      y_i = A * T(t_i - t0) + B,  t_i = 5 ns * i (i = sample index from the first sample of the hit)
//   pedestal   mean of the first nPedSamples samples, subtracted before the fit; B is a residual baseline, |B| <= maxBaseline
//   A, B       solved linearly at fixed t0 (least squares), B clamped
//   t0         1. coarse scan in 1 ns steps over +-8 ns around t_max - x_peak, samples [t_max - 70, t_max - 20] ns
//              (t_max = time of the largest sample, x_peak = x of the highest point of the template graph), golden-section refinement
//              in +-1 ns; 2. two iterations with the window [peak + fitLo, peak + fitHi] ns from the fitted peak and a
//              golden-section refinement in +-1.5 ns
//   time       t0 + x_peak, i.e. the time of the fitted peak, in ns from the first sample of the hit
//   chi2/ndf   residual sum of squares / noise^2 / (N - 3) on the last window
//   saturated  raw samples >= satADC are left out of every sum (the largest sample is then the first saturated one)
//
// The template is a TGraph (default name "gtempl") of a production-format file (splines_*.root, as splines4_run39.root),
// interpolated with a cubic spline (TSpline3, as caloreco) tabulated every 0.05 ns, then linearly; it is 0 outside its points.
//
// Amplitude-dependent template (shapeMode > 0; templateFile = the output of mc2-ana scripts/template_perf/template_amp_model.py,
// graphs T0, T1, T2 and TParameter A0):  T(x; A) = T0 + u T1 + u^2 T2,  u = ln(A / A0)
//   shapeMode 1  shape tied to the amplitude: still A, B, t0 free; u from the largest sample, then twice from the fitted A
//   shapeMode 2  free shape: y = A T0 + C T1 + B, C solved linearly with A and B (shape = C / A)
//   shapeMode 3  constrained shape: y = A T(u) + D T1 + B, u tied to the amplitude as in 1, D the deviation of this hit from the tied
//                shape, with the penalty (D / (shapeSigma A))^2 noise^2 added to the residual sum (A from the previous iteration);
//                still linear in A, D, B at fixed t0 (shape = u + D / A)
//
#include <string>
#include <vector>

class TGraph;

namespace mu2e {

  class MC2TemplateFit {
    public:
      struct Config {
        std::string templateFile;
        std::string graphName   = "gtempl";
        double      fitLo       = -80.;   // ns from the fitted peak
        double      fitHi       = -31.;   // ns from the fitted peak
        double      noise       = 3.6;    // ADC per sample, for chi2/ndf only
        int         nPedSamples = 4;
        double      maxBaseline = 10.;    // ADC
        int         satADC      = 4090;   // raw samples >= satADC are saturated and left out of the fit
        int         maxPre      = 16;     // samples before the largest one that may enter the fit (as in the MC2 studies); < 0: all
        int         shapeMode   = 0;      // 0 fixed template; 1 shape tied to the fitted amplitude; 2 free shape coefficient; 3 constrained
        double      shapeSigma  = 0.5;    // shapeMode 3: spread of the shape between hits, in units of u (prior width of D / A)
      };
      struct Result {
        bool   ok        = false;
        double time      = 0.;   // ns from the first sample
        double amplitude = 0.;   // ADC (template normalised to 1 at its maximum)
        double baseline  = 0.;   // residual baseline B, ADC
        double pedestal  = 0.;   // ADC, mean of the first samples
        double chi2ndf   = 0.;
        int    ndf       = 0;
        int    iMax      = -1;   // index of the largest sample
        int    nSat      = 0;    // saturated samples (left out of the fit)
        double shape     = 0.;   // shapeMode 1: u = ln(A / A0) of the last fit; 2: C / A; 3: u + D / A
      };

      explicit MC2TemplateFit(const Config& cfg);
      MC2TemplateFit(const Config& cfg, const TGraph& g);   // template from a graph in memory (cfg.templateFile not read)
      Result fit(const std::vector<int>& waveform) const;
      Result fit(const std::vector<double>& y, const std::vector<char>& use = {}) const;   // pedestal-subtracted samples; use[i] = 0 leaves sample i out
      double xPeak() const { return xpeak_; }

    private:
      double build(TGraph* g);   // sorted copy of the template graph -> spline table, x_peak; returns the normalisation
      void buildShape(TGraph* g1, TGraph* g2, double norm);   // T1, T2 on the table of T0, divided by the normalisation of T0
      static std::vector<double> table(TGraph* g, double x0, double step, int n);
      double look(const std::vector<double>& t, double x) const;
      double eval(double x) const;    // T0 + u T1 + u^2 T2 at the current u (T0 for shapeMode 0 and 2)
      Result fitOnce(const std::vector<double>& y) const;
      double chi2At(const std::vector<double>& y, int i0, int i1, double t0, double& A, double& B) const;
      double golden(const std::vector<double>& y, int i0, int i1, double a, double b) const;
      mutable const std::vector<char>* use_ = nullptr;   // sample mask of the current fit: one fit at a time per instance (not thread safe)

      Config              cfg_;
      std::vector<double> tab_, tab1_, tab2_;
      double              a0_ = 700.;
      mutable double      u_ = 0., c_ = 0., lam_ = 0.;   // current u (shapeMode 1, 3); C (2) or D (3) of the last chi2At; penalty weight (3)
      double              x0_ = 0., step_ = 0.05, xpeak_ = 0.;
  };

}
#endif
