#include "LinearSystem.h"
#include "Solver.h"
#include <MatlabDataArray.hpp>
#include <MatlabEngine.hpp>

using namespace matlab::engine;

class MatLabSession {
  public:
    std::string name = "movingHeatSourceSolver";
    std::unique_ptr<MATLABEngine> matlabPtr;
    MatLabSession() {
      matlabPtr = startMATLAB();
    }
};

void MatLabSolver::solve(LinearSystem *ls) {
}
