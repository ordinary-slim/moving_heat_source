#ifndef MATLABSESSION
#define MATLABSESSION
#include "Solver.h"
#include <stdexcept>
#ifdef hasMATLAB
#include <MatlabDataArray.hpp>
#include <MatlabEngine.hpp>

using namespace matlab::engine;

class MatLabSession {
  public:
    std::string name = "movingHeatSourceSolver";
    std::unique_ptr<MATLABEngine> ptr;
    MatLabSession() {
      ptr = startMATLAB();
    }
};

class MatLabSolver : public Solver {
  public:
    void solve(LinearSystem *ls);
    void setSession(MatLabSession* ms);
  private:
    MatLabSession* matlab;
};

#else

class MatLabSolver : public Solver {
  public:
    MatLabSolver() {
      throw std::invalid_argument("This library is compiled without MATLAB\n");
    }
};

#endif
#endif
