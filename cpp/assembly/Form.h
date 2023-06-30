#ifndef FORM
#define FORM
#include "../Problem.h"
#include "../mesh/Element.h"
#include <Eigen/Core>

class Form {
  public:
    Form( const Problem *problem ) {
      p = problem;
    }
    // Initializations at element level
    virtual void preGauss(const mesh::Element *e){};
    virtual void inGauss(int igp, const mesh::Element *e){};
  protected:
    const Problem *p;
};

class LinearForm : public Form {
  public:
    LinearForm( const Problem *problem )
        : Form( problem ) {
    }
    virtual double contribute( int igp, int inode, const mesh::Element *e ) {
      return 0.0;
  }
};

class BilinearForm : public Form {
  public:
    BilinearForm( const Problem *problem )
      : Form( problem ) {
    }
    virtual double contribute( int igp, int inode, int jnode, const mesh::Element *e ) {
      return 0.0;
    }
};

class MassForm : public BilinearForm {
  public:
    MassForm( const Problem *problem )
      : BilinearForm( problem ) {
    }
    double contribute( int igp, int inode, int jnode, const mesh::Element *e ) {
      return p->density * p->specificHeat *
        e->BaseGpVals[inode][igp]*e->BaseGpVals[jnode][igp]*
        e->gpweight[igp]  * e->vol;
    }
};

class DiffusionForm : public BilinearForm {
  public:
    DiffusionForm( const Problem *problem )
      : BilinearForm( problem ) {
    }
    double contribute( int igp, int inode, int jnode, const mesh::Element *e ) {
      double ip = e->GradBaseGpVals[inode][igp].dot( e->GradBaseGpVals[jnode][igp] );
      return p->conductivity * e->gpweight[igp] * ip * e->vol;
    }
};

class AdvectionForm : public BilinearForm {
  public:
    AdvectionForm( const Problem *problem )
      : BilinearForm( problem ) {
    }
    double contribute( int igp, int inode, int jnode, const mesh::Element *e ) {
      double ip = e->GradBaseGpVals[jnode][igp].dot(p->advectionSpeed);
      return p->density * p->specificHeat * e->gpweight[igp] * (ip * e->BaseGpVals[inode][igp]) * e->vol;
    }
};

class SourceForm : public LinearForm {
  private:
    Eigen::Vector3d xgp;
  public:
    SourceForm( const Problem *problem )
      : LinearForm( problem ) {
    }
    void inGauss(int igp, const mesh::Element *e){
      xgp = e->gpos.row( igp );
    }
    double contribute( int igp, int inode, const mesh::Element *e ) {
      return e->gpweight[igp] * e->BaseGpVals[inode][igp] * e->vol *
        p->mhs->operator()(xgp, p->time);
    }
};

class ASSSBilinearForm : public BilinearForm {
  private:
    double h = -1, tau = -1;
    double SCA = 2, SCD = 4;//stabilization cte advection / diffusion
    Eigen::Vector3d advectionSpeed;
    double norm_advectionSpeed = 0.0;
  public:
    ASSSBilinearForm( const Problem *problem )
      : BilinearForm( problem ) {
        advectionSpeed = p->advectionSpeed;
        norm_advectionSpeed = advectionSpeed.norm();
    }
    void preGauss(const mesh::Element *e){
      //Compute tau
      h = e->getSizeAlongVector( p->advectionSpeed );
      double advectionEstimate = h / SCA / (p->density*p->specificHeat*norm_advectionSpeed);
      tau = advectionEstimate;
      if (p->conductivity != 0) {
        double diffusionEstimate = pow(h, 2) / (SCD * p->conductivity);
        tau = 1 / ( 1/advectionEstimate + 1/diffusionEstimate );
      } else {
        tau = advectionEstimate;
      }
    }
    double contribute( int igp, int inode, int jnode, const mesh::Element *e ) {
      return (e->gpweight[igp] * e->vol)* p->density * p->specificHeat * tau *
        e->GradBaseGpVals[inode][igp].dot( advectionSpeed ) *
        e->GradBaseGpVals[jnode][igp].dot( advectionSpeed );
    }
};

class ASSSLinearForm : public LinearForm {
  private:
    Eigen::Vector3d xgp;
    double h = -1, tau = -1;
    double SCA = 2, SCD = 4;//stabilization cte advection / diffusion
    Eigen::Vector3d advectionSpeed;
    double norm_advectionSpeed = 0.0;
  public:
    ASSSLinearForm( const Problem *problem )
      : LinearForm( problem ) {
        advectionSpeed = p->advectionSpeed;
        norm_advectionSpeed = advectionSpeed.norm();
    }
    void preGauss(const mesh::Element *e){
      //Compute tau
      h = e->getSizeAlongVector( p->advectionSpeed );
      double advectionEstimate = h / SCA / (p->density*p->specificHeat*norm_advectionSpeed);
      tau = advectionEstimate;
      if (p->conductivity != 0) {
        double diffusionEstimate = pow(h, 2) / (SCD * p->conductivity);
        tau = 1 / ( 1/advectionEstimate + 1/diffusionEstimate );
      } else {
        tau = advectionEstimate;
      }
    }
    void inGauss(int igp, const mesh::Element *e){
      xgp = e->gpos.row( igp );
    }
    double contribute( int igp, int inode, const mesh::Element *e ) {
        double f_xgp = p->mhs->operator()(xgp, p->time);
        return (e->gpweight[igp] * e->vol) * tau * f_xgp *
          e->GradBaseGpVals[inode][igp].dot( advectionSpeed );;
    }
};

class ConvectionBilinearForm : public BilinearForm {
  public:
    ConvectionBilinearForm( const Problem *problem )
      : BilinearForm( problem ) {
    }
    double contribute( int igp, int inode, int jnode, const mesh::Element *e ) {
      return e->gpweight[igp] * e->vol * e->BaseGpVals[inode][igp] * e->BaseGpVals[jnode][igp] *
        p->convectionCoeff / p->conductivity;
    }
};

class ConvectionLinearForm : public LinearForm {
  public:
    ConvectionLinearForm( const Problem *problem )
      : LinearForm( problem ) {
    }
    double contribute( int igp, int inode, const mesh::Element *e ) {
        return  e->gpweight[igp] * e->vol * e->BaseGpVals[inode][igp] * p->Tenv *
        p->convectionCoeff / p->conductivity;
    }
};

class NeumannLinearForm : public LinearForm {
  public:
    NeumannLinearForm( const Problem *problem )
      : LinearForm( problem ) {
    }
    void inGauss(int igp, const mesh::Element *e) {
      normalDerivative = p->neumannFluxes[e->ient][igp] / p->conductivity;
    }
    double contribute( int igp, int inode, const mesh::Element *e ) {
        return e->gpweight[igp] * e->vol * e->BaseGpVals[inode][igp] * normalDerivative;
    }
  private:
    double normalDerivative = 0.0;
};
#endif
