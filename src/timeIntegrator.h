#ifndef TIMEINTEGRATOR
#include <Eigen/Core>

using namespace std;

class TimeIntegratorHandler {
  public:
      double lhsCoeff;
      Eigen::VectorXd rhsCoeff;

      // integrator
      int currentIntegrator = 1;
      int desiredIntegrator = 1;
      int nstepsRequired = 1;
      int nstepsStored   = 0;

      void setRequiredSteps(int icase) {
        if (icase == 0) {//ForwardEuler
          desiredIntegrator = 0;
          nstepsRequired = 1;
        } else if (icase == 1) {//BDF1
          desiredIntegrator = 1;
          nstepsRequired = 1;
        } else if (icase == 2) {//BDF2
          desiredIntegrator = 2;
          nstepsRequired = 2;
        }
      }

      void setCoeffs(int icase) {
        if (icase == 0) {//ForwardEuler
        } else if (icase == 1) {//BDF1
          lhsCoeff = 1;
          rhsCoeff.resize(1);
          rhsCoeff << 1;
        } else if (icase == 2) {//BDF2
          lhsCoeff = 1.5;
          rhsCoeff.resize(2);
          rhsCoeff << 2, -0.5;
        }
      }
    };
#define TIMEINTEGRATOR
#endif
