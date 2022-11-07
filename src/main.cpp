#include <iostream>
#include <vector>
#include "line.h"
#include "mesh.h"
#include <numeric>
#include <Eigen/Core>
#include <Eigen/Sparse>
using namespace std;

//plotting
#include <matplotlib-cpp/matplotlibcpp.h>
namespace plt = matplotlibcpp;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

void compute_pulse(Eigen::VectorXd &pulse, double x0, double radius, double P, Mesh &m);

int main() {
  // params
  double L=10.0, radius=2.0, P=100.0;
  double x0 = 0.0, speed=10, t=0.0, Tfinal=10.0;
  double rho = 4000, k = 40, c_p = 10;
  // numerical params
  double dt;
  double CFL;
  double m_ij, k_ij, ip;
  int nels = 100;
  int iter=0, maxIter=100;

  Mesh mesh;
  mesh.initialize1DMesh( 0.0, L, nels);

  // set timestep
  /*
  CFL = 5.0;
  dt = CFL * (L / nels) / speed;
  */
  dt = 0.05;
  maxIter = L / (speed*dt) + 5;

  // Assembly
  SpMat M(mesh.nnodes, mesh.nnodes); // mass mat
  vector<T> M_coeffs;
  M_coeffs.reserve( 3*mesh.nnodes );
  SpMat K(mesh.nnodes, mesh.nnodes); // stiffness mat
  vector<T> K_coeffs;
  K_coeffs.reserve( 3*mesh.nnodes );
  Eigen::VectorXd pulse( mesh.nnodes ); // source term

  Line l;

  while ((t < Tfinal)&&(iter < maxIter)) {
      // matrices assembly
      M.setZero();
      K.setZero();
      for (int ielem = 0; ielem < mesh.nels; ielem++ ) {
        l = mesh.getElement( ielem );
        for (int inode = 0; inode < l.nnodes; inode++) {
          for (int jnode = 0; jnode < l.nnodes; jnode++) {
            m_ij = 0;
            k_ij = 0;
            for (int igp = 0; igp < l.nnodes; igp++) {
              // mass matrix
              m_ij += l.gpweight[igp] * l.baseFunGpVals[inode][igp]*l.baseFunGpVals[jnode][igp]*l.vol;
              // stiffness matrix
              ip = inner_product(l.baseFunGradGpVals[inode][igp].begin(),
                    l.baseFunGradGpVals[inode][igp].end(),
                    l.baseFunGradGpVals[jnode][igp].begin(),
                    0.0);
              k_ij += l.gpweight[igp] * ip * l.vol;
            }
            m_ij *= rho*c_p;
            k_ij *= k;
            M_coeffs.push_back( T( l.con[inode], l.con[jnode], m_ij ) );
            K_coeffs.push_back( T( l.con[inode], l.con[jnode], k_ij ) );
          }
        }
      }
    M.setFromTriplets( M_coeffs.begin(), M_coeffs.end() );
    K.setFromTriplets( K_coeffs.begin(), K_coeffs.end() );

    //update laser position
    x0 = t * speed;
    //load vector computation/assembly
    pulse.setZero();
    compute_pulse(pulse, x0, radius, P, mesh);

    //solve
    
    //update plot
    //convert from Eigen::Vector to std::vector
    vector<double> vpulse( pulse.data(), pulse.data() + pulse.size());
    plt::clf();
    cout << "current positoin:" << x0 << endl;
    plt::plot( std::vector<double>({x0}), std::vector<double>({0.0}),
        {{"marker", "x"}, {"color", "red"}, {"markersize", "20"}});
    //plt::plot( mesh.pos, vpulse );
    plt::xlim( 0.0, L );
    double maxPower = max( maxPower, *max_element( pulse.begin(), pulse.end() ) );
    //plt::ylim( 0.0, 1.1*maxPower );

    plt::pause( 0.04 );

    ++iter;
    t += dt;
  }


}
