#include <vector>

using std::vector;
namespace LBFGS_ns{
  class LBFGS{
    private : 
      int M_;
      int N_;
      int k_;
      int maxiter_;
      double tol_;
      double maxstep_;
      int nfev_;
      double max_f_rise_;

      // variables representing the state of the system
      std::vector<double> x_;
      double f_;
      std::vector<double> g_;
      double rms_;

      // places to store the lbfgs memory
      std::vector<vector<double> > s_;
      std::vector<vector<double> > y_;
      std::vector<double> rho_;
      double H0_;

      // 
      std::vector<double> step_;

    public :
      LBFGS(double const * x0, int N, int M);
      ~LBFGS();
      void one_iteration();
      void run();

    private :
      void update_memory(
          std::vector<double> & xold,
          std::vector<double> & gold,
          std::vector<double> & xnew,
          std::vector<double> & gnew
          );
      void compute_lbfgs_step();
      void backtracking_linesearch();
      int stop_criterion_satisfied();
      void compute_func_gradient(std::vector<double> & x, double & energy,
          std::vector<double> & gradient);

  };
}
