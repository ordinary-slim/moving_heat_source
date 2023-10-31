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
    const Problem *p = nullptr;
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
    MassForm()
      : BilinearForm( nullptr ) {
      }
    double contribute( int igp, int inode, int jnode, const mesh::Element *e ) {
      return 
        e->BaseGpVals[inode][igp]*e->BaseGpVals[jnode][igp]*
        e->gpweight[igp]  * e->vol;
    }
};

class TimeMassForm : public BilinearForm {
  public:
    TimeMassForm( const Problem *problem )
      : BilinearForm( problem ) {
      }
    double contribute( int igp, int inode, int jnode, const mesh::Element *e ) {
      return 
        p->materials[e->imat].density*p->materials[e->imat].specificHeat*
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
      return p->materials[e->imat].conductivity * e->gpweight[igp] * ip * e->vol;
    }
};

class AdvectionForm : public BilinearForm {
  public:
    AdvectionForm( const Problem *problem )
      : BilinearForm( problem ) {
    }
    double contribute( int igp, int inode, int jnode, const mesh::Element *e ) {
      double ip = e->GradBaseGpVals[jnode][igp].dot(p->advectionSpeed);
      return p->materials[e->imat].density * p->materials[e->imat].specificHeat * e->gpweight[igp] * (ip * e->BaseGpVals[inode][igp]) * e->vol;
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

class LumpedSourceForm : public LinearForm {
  private:
    bool isElementHeated = false;
    const heat::LumpedHeatSource* lumpedHs;
  public:
    LumpedSourceForm( const Problem *problem )
      : LinearForm( problem ) {
        // Casting a unique pointer!
        // Safe because LumpedSourceForm lifetime is short compared to problem
        // https://stackoverflow.com/a/36120483/12948600
        lumpedHs = static_cast<heat::LumpedHeatSource*>(problem->mhs.get());
    }
    void preGauss(const mesh::Element *e){
      isElementHeated = bool( lumpedHs->heatedElements[e->ient] );
    }
    double contribute( int igp, int inode, const mesh::Element *e ) {
      if (isElementHeated) {
        return e->gpweight[igp] * e->BaseGpVals[inode][igp] * e->vol * lumpedHs->pd;
      } else {
        return 0.0;
      }
    }
};

class SUPG {
  public:
    double h = -1, tau = -1;
    Eigen::Vector3d advectionSpeed;
    double norm_advectionSpeed = 0.0;

    SUPG( const Problem *problem ) {
      advectionSpeed = problem->advectionSpeed;
      norm_advectionSpeed = advectionSpeed.norm();
    }

    void setTau(const mesh::Element *e, const Problem *p) {
      /*
       * (Codina, 2000)
       */
      h = e->getSizeAlongVector( p->advectionSpeed );
      tau = pow(h, 2) / (2 * h * p->materials[e->imat].density * p->materials[e->imat].specificHeat * norm_advectionSpeed +
          4 * p->materials[e->imat].conductivity);
    }
};


class SUPGTimeBilinearForm : public BilinearForm, public SUPG {
  public:
    SUPGTimeBilinearForm( const Problem *problem )
      : BilinearForm( problem ),
        SUPG( problem ) {
          // pass
    }
    void preGauss(const mesh::Element *e){
      setTau(e, p);
    }
    double contribute( int igp, int inode, int jnode, const mesh::Element *e ) {
      return (e->gpweight[igp] * e->vol)*tau*
        p->materials[e->imat].density * p->materials[e->imat].specificHeat * e->BaseGpVals[jnode][igp] *
        p->materials[e->imat].density * p->materials[e->imat].specificHeat * e->GradBaseGpVals[inode][igp].dot( advectionSpeed );
    }
};

class SUPGBilinearForm : public BilinearForm, public SUPG {
  public:
    SUPGBilinearForm( const Problem *problem )
      : BilinearForm( problem ),
        SUPG( problem ) {
          // pass
    }
    void preGauss(const mesh::Element *e){
      setTau(e, p);
    }
    double contribute( int igp, int inode, int jnode, const mesh::Element *e ) {
      return (e->gpweight[igp] * e->vol)*tau*
        p->materials[e->imat].density * p->materials[e->imat].specificHeat * e->GradBaseGpVals[jnode][igp].dot( advectionSpeed ) *
        p->materials[e->imat].density * p->materials[e->imat].specificHeat * e->GradBaseGpVals[inode][igp].dot( advectionSpeed );
    }
};

class SUPGLinearForm : public LinearForm, public SUPG {
  private:
    Eigen::Vector3d xgp;
  public:
    SUPGLinearForm( const Problem *problem )
      : LinearForm( problem ),
        SUPG( problem )
    {
          //pass
    }
    void preGauss(const mesh::Element *e){
      setTau(e, p);
    }
    void inGauss(int igp, const mesh::Element *e){
      xgp = e->gpos.row( igp );
    }
    double contribute( int igp, int inode, const mesh::Element *e ) {
        double f_xgp = p->mhs->operator()(xgp, p->time);
        return (e->gpweight[igp] * e->vol)*tau*
          f_xgp *
          p->materials[e->imat].density * p->materials[e->imat].specificHeat * e->GradBaseGpVals[inode][igp].dot( advectionSpeed );
    }
};

class ConvectionBilinearForm : public BilinearForm {
  public:
    ConvectionBilinearForm( const Problem *problem )
      : BilinearForm( problem ) {
    }
    double contribute( int igp, int inode, int jnode, const mesh::Element *e ) {
      return e->gpweight[igp] * e->vol * e->BaseGpVals[inode][igp] * e->BaseGpVals[jnode][igp] *
        p->materials[e->imat].convectionCoeff;
    }
};

class ConvectionLinearForm : public LinearForm {
  public:
    ConvectionLinearForm( const Problem *problem )
      : LinearForm( problem ) {
    }
    double contribute( int igp, int inode, const mesh::Element *e ) {
        return  e->gpweight[igp] * e->vol * e->BaseGpVals[inode][igp] * p->Tenv *
        p->materials[e->imat].convectionCoeff;
    }
};

class NeumannLinearForm : public LinearForm {
  public:
    NeumannLinearForm( const Problem *problem )
      : LinearForm( problem ) {
    }
    void inGauss(int igp, const mesh::Element *e) {
      normalDerivative = p->neumannFluxes[e->ient][igp] / p->materials[e->imat].conductivity;
    }
    double contribute( int igp, int inode, const mesh::Element *e ) {
        return e->gpweight[igp] * e->vol * e->BaseGpVals[inode][igp] * normalDerivative;
    }
  private:
    double normalDerivative = 0.0;
};
#endif
