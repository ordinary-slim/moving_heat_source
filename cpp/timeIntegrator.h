#ifndef TIMEINTEGRATOR
#include <Eigen/Core>

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
        if (icase == 1) {//BDF1
          desiredIntegrator = 1;
          nstepsRequired = 1;
        } else if (icase == 2) {//BDF2
          desiredIntegrator = 2;
          nstepsRequired = 2;
        } else if (icase == 3) {//BDF3
          desiredIntegrator = 3;
          nstepsRequired = 3;
        } else if (icase == 4) {//BDF4
          desiredIntegrator = 4;
          nstepsRequired = 4;
        }
      }

      void setCurrentIntegrator(int icase) {
        if (icase == 1) {//BDF1
          currentIntegrator = 1;
          lhsCoeff = 1;
          rhsCoeff.resize(1);
          rhsCoeff << 1;
        } else if (icase == 2) {//BDF2
          currentIntegrator = 2;
          lhsCoeff = 1.5;
          rhsCoeff.resize(2);
          rhsCoeff << 2, -0.5;
        } else if (icase == 3) {//BDF3
          currentIntegrator = 3;
          lhsCoeff = 11.0/6.0;
          rhsCoeff.resize(3);
          rhsCoeff << 3.0, -1.5, 1.0/3.0;
        } else if (icase == 4) {//BDF4
          currentIntegrator = 4;
          lhsCoeff = 25.0/12.0;
          rhsCoeff.resize(4);
          rhsCoeff << 4.0, -3.0, 4.0/3.0, -1.0/4.0;
        }
      }
    };
#define TIMEINTEGRATOR
#endif
