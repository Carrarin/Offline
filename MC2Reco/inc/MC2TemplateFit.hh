#ifndef MC2Reco_MC2TemplateFit_hh
#define MC2Reco_MC2TemplateFit_hh
//
// MC2TemplateFit: template fit of one calorimeter waveform.
//
// It is the method of the MC2 timing studies (mc2-ana scripts/laser_template_range_opt.C, function fitPeak, pedestal
// mode "clamp"), written without art and without Minuit so that it can be used both in the art modules and in plain
// programs.
//
// Model of the waveform:
//
//   y_i = A * T(t_i - t0) + B        t_i = 5 ns * i, i = sample index from the first sample of the hit
//
//   T        the template, normalised to 1 at its maximum
//   A        amplitude (ADC)
//   t0       time shift of the template (ns)
//   B        residual baseline (ADC), limited to |B| <= maxBaseline
//
// Steps of a fit:
//
//   1. pedestal = mean of the first nPedSamples samples, subtracted from every sample;
//      raw samples >= satADC are saturated and left out of every sum
//   2. at fixed t0, A and B are the linear least-squares solution (B clamped to the limit)
//   3. first estimate of t0: scan in 1 ns steps over +-8 ns around t_max - x_peak, using the samples between 70 ns
//      and 20 ns before the largest sample, then a golden-section refinement in +-1 ns
//      (t_max = time of the largest sample, x_peak = x of the highest point of the template graph)
//   4. two more iterations with the window [peak + fitLo, peak + fitHi] ns around the fitted peak, each refined with a
//      golden-section search in +-1.5 ns
//
// Result:
//
//   time       t0 + x_peak, the time of the fitted peak, in ns from the first sample of the hit
//   chi2/ndf   residual sum of squares / noise^2 / (N - 3) on the last window
//
// The template is a TGraph (default name "gtempl") of a production-format file (splines_*.root, as
// splines4_run39.root), interpolated with a cubic spline (TSpline3, as caloreco), tabulated every 0.05 ns and then
// interpolated linearly; it is 0 outside its points.
//
// Amplitude-dependent template (shapeMode > 0). The template file is then the output of mc2-ana
// scripts/template_perf/template_amp_model.py, with the graphs T0, T1, T2 and the TParameter A0:
//
//   T(x; A) = T0 + u T1 + u^2 T2,   u = ln(A / A0)
//
//   shapeMode 1  shape tied to the amplitude: A, B, t0 free; u from the largest sample, then twice from the fitted A
//   shapeMode 2  free shape: y = A T0 + C T1 + B, with C solved linearly together with A and B (shape = C / A)
//   shapeMode 3  constrained shape: y = A T(u) + D T1 + B, u tied to the amplitude as in 1, D the deviation of this hit
//                from the tied shape, with the penalty (D / (shapeSigma A))^2 noise^2 added to the residual sum
//                (A from the previous iteration); still linear in A, D, B at fixed t0 (shape = u + D / A)
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
        double      noise       = 3.6;    // ADC per sample, used for chi2/ndf only
        int         nPedSamples = 4;
        double      maxBaseline = 10.;    // ADC
        int         satADC      = 4090;   // raw samples >= satADC are saturated and left out of the fit
        int         maxPre      = 16;     // samples before the largest one that may enter the fit; < 0: all
        int         shapeMode   = 0;      // 0 fixed, 1 tied to the amplitude, 2 free shape, 3 constrained shape
        double      shapeSigma  = 0.5;    // shapeMode 3: spread of the shape between hits, in units of u
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

      // Template read from cfg.templateFile.
      explicit MC2TemplateFit(const Config& cfg);

      // Template from a graph in memory (cfg.templateFile is not read).
      MC2TemplateFit(const Config& cfg, const TGraph& graph);

      // Fit of a raw waveform (ADC samples): pedestal subtraction, saturated samples left out, then the fit.
      Result fit(const std::vector<int>& waveform) const;

      // Fit of pedestal-subtracted samples; use[i] = 0 leaves sample i out (empty: all samples are used).
      Result fit(const std::vector<double>& y, const std::vector<char>& use = {}) const;

      double xPeak() const { return xPeak_; }

    private:
      // Template tables
      double buildTemplate(TGraph* graph);
      void buildShapeTerms(TGraph* graph1, TGraph* graph2, double normalisation);
      static std::vector<double> tabulate(TGraph* graph, double xStart, double step, int nPoints);
      double interpolate(const std::vector<double>& table, double x) const;
      double templateAt(double x) const;

      // Fit steps
      Result fitWithFixedShape(const std::vector<double>& y) const;
      double chi2AtTime(const std::vector<double>& y, int iFirst, int iEnd, double t0, double& A, double& B) const;
      double chi2WithShapeTerm(const std::vector<double>& y, int iFirst, int iEnd, double t0,
                               double& A, double& B) const;
      double goldenSectionMinimum(const std::vector<double>& y, int iFirst, int iEnd, double tLow, double tHigh) const;
      bool isUsed(int i) const { return use_ == nullptr || (*use_)[i]; }

      Config cfg_;

      std::vector<double> table0_;   // T0 (or the only template), normalised to 1 at its maximum
      std::vector<double> table1_;   // T1 / normalisation of T0 (amplitude model)
      std::vector<double> table2_;   // T2 / normalisation of T0 (amplitude model)
      double xStart_ = 0.;           // x of the first point of the tables
      double step_   = 0.05;         // ns between the points of the tables
      double xPeak_  = 0.;           // x of the highest point of the template graph
      double a0_     = 700.;         // amplitude model: reference amplitude A0

      // State of the fit in progress: one fit at a time per instance (not thread safe).
      mutable const std::vector<char>* use_ = nullptr;   // sample mask
      mutable double u_ = 0.;                            // shape parameter u (shapeMode 1, 3)
      mutable double shapeCoefficient_ = 0.;             // C (shapeMode 2) or D (shapeMode 3), last chi2 evaluation
      mutable double penaltyWeight_ = 0.;                // shapeMode 3: weight of the penalty on D
  };

}
#endif
