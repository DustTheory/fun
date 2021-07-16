#include <cmath>
#include <limits>
#include <string.h>
#include "Debug.H"
#include "Level2D.H"
#include "Fort2D.H"
#include "MemoryDebugger.H"
#include "Tracer.H"

struct const_common {
  double e,pi,rt2, half,one,two,zero;
  int debug_on;
};
extern const_common F77NAME(const);

int Level::smoother_iterations=1;
Level::SMOOTHER Level::smoother=Level::GAUSS_SEIDEL_SMOOTHER;
Level::RESTRICTION_PROLONGATION Level::restriction_prolongation=
  Level::FINITE_ELEMENT_RESTRICTION_PROLONGATION;

const char* Level::smoother_name[3]={
  "gauss_seidel_smoother","gauss_seidel_red_black_smoother",
  "richardson_smoother"
};
const char* Level::restriction_prolongation_name[2]={
  "algebraic_multigrid_restriction_prolongation",
  "finite_element_restriction_prolongation"
};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::setArrayBounds() {
  for (int id=0;id<2;id++) {
    fc[id]=0;     // corner indices
    lc[id]=n[id]; // corner indices
    ifirst[id]=1;      // interior corner indices
    ilast[id]=n[id]-1; // interior corner indices
    fe[id]=0;       // cell indices
    le[id]=n[id]-1; // cell indices
  }
  number_corners=(n[0]+1)*(n[1]+1); // consistent with runOnce
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//constructs the hierarchy of Level's, beginning with the finest
Level::Level(const int *nc,const double *m,int si,SMOOTHER s,
RESTRICTION_PROLONGATION rp) : coarser(0),finer(0),level_number(0),
mu(0.),matrix(0),Ax_b(0),first_residual(0),solution_increment(0),
prolongation_vector(0),prolongcell(0),prolongside0(0),prolongside1(0) {
//TRACER_CALL(t,"Level::Level");
  smoother_iterations=si;
  smoother=s;
  restriction_prolongation=rp;
//restriction_prolongation=FINITE_ELEMENT_RESTRICTION_PROLONGATION;
#ifdef DEBUG
//cout << "\tsmoother = " << static_cast<int>(smoother) << " "
//     << smoother_name[smoother] << endl;
//cout << "\trestriction_prolongation = "
//     << static_cast<int>(restriction_prolongation) << " "
//     << restriction_prolongation_name[restriction_prolongation]
//     << endl;
#endif
  for (int id=0;id<2;id++) n[id]=nc[id];
  setArrayBounds();
#ifdef DEBUG
//cout << "\tn = " << n[0] << " " << n[1] << endl;
#endif
  if (n[0]%2==0 && n[1]%2==0 && min(n[0],n[1])>2) {
    coarser=OPERATOR_NEW Level(this);
  }
  setup(m);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//called from public Level::Level
Level::Level(Level *f) : coarser(0),finer(f),
level_number(f->level_number+1),mu(0.),matrix(0),Ax_b(0),
first_residual(0),solution_increment(0),prolongation_vector(0),
prolongcell(0),prolongside0(0),prolongside1(0) {
//TRACER_CALL(t,"Level::Level");
  ASSERT(finer->n[0]%2==0 && finer->n[1]%2==0);
  for (int id=0;id<2;id++) n[id]=finer->n[id]/2;
  setArrayBounds();
#ifdef DEBUG
//cout << "\tn = " << n[0] << " " << n[1] << endl;
#endif
  if (n[0]%2==0 && n[1]%2==0 && min(n[0],n[1])>2) {
    coarser=OPERATOR_NEW Level(this);
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//called from public Level::Level after all Level's are constructed
void Level::setup(const double *m) {
//TRACER_CALL(t,"Level::setup");
#ifdef DEBUG
//cout << "\tlevel_number = " << level_number << endl;
#endif

//allocate space for matrix, residuals, solution increment
  matrix=OPERATOR_NEW_BRACKET(double,number_corners*9);
  Ax_b=OPERATOR_NEW_BRACKET(double,number_corners);
  first_residual=OPERATOR_NEW_BRACKET(double,number_corners);
  solution_increment=OPERATOR_NEW_BRACKET(double,number_corners);
  if (coarser!=0) {
    prolongation_vector=OPERATOR_NEW_BRACKET(double,number_corners);
    for (int i=0;i<number_corners;i++) prolongation_vector[i]=0.;
  }

//
  if (finer==0) {
//  TRACER_CALL(t,"Level::setup finer==0");
//  no finer level exists ==> copy matrix
    memcpy(matrix,m,number_corners*9*sizeof(*matrix));
#ifdef DEBUG
//  for (int i=-4;i<=4;i++) {
//    cout << "\tmatrix[" << i << "] = " << endl;
//    F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//      ifirst[0],ifirst[1],ilast[0],ilast[1],matrix+(i+4)*number_corners);
//  }
#endif
  } else {
//  TRACER_CALL(t,"Level::setup finer!=0");
//  compute matrix on this level from matrix on finer level
    for (int i=0;i<finer->number_corners;i++) {
      finer->first_residual[i]=0.;
    }
    for (int START1=1;START1<=min(3,n[1]-1);START1++) {
      for (int START0=1;START0<=min(3,n[0]-1);START0++) {
        for (int I=0;I<number_corners;I++) solution_increment[I]=0.;
        for (int i=0;i<finer->number_corners;i++) {
          finer->solution_increment[i]=0.;
        }
#ifdef DEBUG
//      cout << "\n\tSTART = " << START0 << " " << START1 << endl;
//      cout << "\tfine solution_increment:" << endl;
//      F77NAME(rbugcell)(finer->fc[0],finer->fc[1],finer->lc[0],finer->lc[1],
//        finer->fc[0],finer->fc[1],finer->lc[0],finer->lc[1],finer->solution_increment);
#endif

        F77_NAME(vector_setup)(fc[0],fc[1],lc[0],lc[1],
          START0,START1,solution_increment);
#ifdef DEBUG
//      cout << "\tcoarse solution_increment:" << endl;
//      F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//        fc[0],fc[1],lc[0],lc[1],solution_increment);
#endif
        finer->prolongation(); // coarse soln_inc->fine soln_inc
#ifdef DEBUG
//      cout << "\tfine solution_increment:" << endl;
//      F77NAME(rbugcell)(finer->fc[0],finer->fc[1],finer->lc[0],finer->lc[1],
//        finer->fc[0],finer->fc[1],finer->lc[0],finer->lc[1],
//        finer->solution_increment);
#endif
#ifdef INDEF
        for (int i=0;i<finer->number_corners;i++) {
          finer->Ax_b[i]=numeric_limits<double>::infinity();
        }
#endif
        finer->computeResidual(finer->solution_increment,
          finer->first_residual,finer->Ax_b);
#ifdef DEBUG
//      cout << "\tfine Ax:" << endl;
//      F77NAME(rbugcell)(finer->fc[0],finer->fc[1],finer->lc[0],finer->lc[1],
//        finer->fc[0],finer->fc[1],finer->lc[0],finer->lc[1],finer->Ax_b);
#endif
        for (int i=0;i<number_corners;i++) first_residual[i]=0.;
        finer->restriction(); // fine Ax_b -> coarse first_residual
#ifdef DEBUG
//      cout << "\tcoarse Ax:" << endl;
//      F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//        fc[0],fc[1],lc[0],lc[1],first_residual);
#endif
        F77_NAME(matrix_setup)(fc[0],fc[1],lc[0],lc[1],START0,START1,
          first_residual,matrix);
      }
    }
#ifdef DEBUG
//  for (int i=-4;i<=4;i++) {
//    cout << "\tmatrix[" << i << "] = " << endl;
//    F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//      ifirst[0],ifirst[1],ilast[0],ilast[1],matrix+(i+4)*number_corners);
//  }
#endif
  }
  if (smoother==RICHARDSON_SMOOTHER) {
//  Gerschgorin circle theorem to get Richardson mu
    mu=0.;
    for (int i=0;i<number_corners;i++) {
      int i1=i/(n[0]+1);
      int i0=i%(n[0]+1);
      if (i0>1 && i0<n[0]-1 && i1>1 && i1<n[1]-1) {
        double sum=abs(matrix[i])+abs(matrix[i+number_corners])
          +abs(matrix[i+2*number_corners])+abs(matrix[i+3*number_corners])
          +abs(matrix[i+4*number_corners]);
        mu=max(mu,sum);
      }
    }
  }
  for (int i=0;i<number_corners;i++) solution_increment[i]=0.;
//recurse
  if (coarser!=0) {
    if (restriction_prolongation==
    ALGEBRAIC_MULTIGRID_RESTRICTION_PROLONGATION) {
      int (&cn)[2]=coarser->n;
      prolongcell=OPERATOR_NEW_BRACKET(double,cn[0]*cn[1]*4);
      prolongside0=OPERATOR_NEW_BRACKET(double,(cn[0]+1)*cn[1]*2);
      prolongside1=OPERATOR_NEW_BRACKET(double,(cn[1]+1)*cn[0]*2);
      F77_NAME(setup_am)(coarser->fe[0],coarser->fe[1],coarser->le[0],coarser->le[1],
        fe[0],fe[1],le[0],le[1],matrix,
        prolongcell,prolongside0,prolongside1);
    }
    coarser->setup(0);
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int Level::numberLevels() const {
  if (coarser==0) return level_number+1;
  else return coarser->numberLevels();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Level::~Level() {
//TRACER_CALL(t,"Level::~Level");
  if (coarser!=0) delete coarser;
  delete [] Ax_b;
  delete [] first_residual;
  delete [] matrix;
  delete [] solution_increment;
  if (prolongation_vector!=0) delete [] prolongation_vector;
  if (prolongcell!=0) delete [] prolongcell;
  if (prolongside0!=0) delete [] prolongside0;
  if (prolongside1!=0) delete [] prolongside1;
  coarser=0;
  finer=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//applies smoother to residual to obtain solution_increment
void Level::preSmooth() const {
//TRACER_CALL(t,"Level::preSmooth");
#ifdef DEBUG
//cout << "\tsmoother = " << static_cast<int>(smoother) << " "
//     << smoother_name[smoother] << endl;
//cout << "\tfirst_residual:" << endl;
//F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//  ifirst[0],ifirst[1],ilast[0],ilast[1],first_residual);
//cout << "\tsolution_increment:" << endl;
//F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//  ifirst[0],ifirst[1],ilast[0],ilast[1],solution_increment);
#endif
//
  double residual_max;
  for (int it=0;it<smoother_iterations;it++) {
    switch (smoother) {
      case GAUSS_SEIDEL_RED_BLACK_SMOOTHER: {
//      TRACER_CALL(t,"Level::preSmooth red-black");
//      coarse nodes, then fine nodes
        F77_NAME(gauss_seidel_red_black)(fc[0],fc[1],lc[0],lc[1],
          ifirst[0],ifirst[1],ilast[0],ilast[1],0,residual_max,matrix,
          first_residual, solution_increment);
        break;
      }
      case GAUSS_SEIDEL_SMOOTHER: {
//      TRACER_CALL(t,"Level::preSmooth to-fro");
        F77_NAME(gauss_seidel_to_fro)(fc[0],fc[1],lc[0],lc[1],
          ifirst[0],ifirst[1],ilast[0],ilast[1],0,residual_max,matrix,
          first_residual, solution_increment);
        break;
      }
      case RICHARDSON_SMOOTHER:
      default: {
//      TRACER_CALL(t,"Level::preSmooth Richardson");
        double *r=OPERATOR_NEW_BRACKET(double,number_corners);
        F77NAME(richardson)(fc[0],fc[1],lc[0],lc[1],
          ifirst[0],ifirst[1],ilast[0],ilast[1],mu,residual_max,matrix,
          first_residual, solution_increment,r);
        delete [] r;
      }
    }
  }
//
/*
//identity preconditioner, for debugging
  for (int i=0;i<number_corners;i++) {
    solution_increment[i]+=first_residual[i];
  }
*/
#ifdef DEBUG
//cout << "\tsolution_increment:" << endl;
//F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//  ifirst[0],ifirst[1],ilast[0],ilast[1],solution_increment);
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//transpose of preSmooth
void Level::postSmooth() const {
//TRACER_CALL(t,"Level::postSmooth");
#ifdef DEBUG
//cout << "\tsmoother = " << static_cast<int>(smoother) << " "
//     << smoother_name[smoother] << endl;
//cout << "\tfirst_residual:" << endl;
//F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//  ifirst[0],ifirst[1],ilast[0],ilast[1],first_residual);
//cout << "\tsolution_increment:" << endl;
//F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//  ifirst[0],ifirst[1],ilast[0],ilast[1],solution_increment);
#endif
//
  double residual_max=numeric_limits<double>::infinity();
  for (int it=0;it<smoother_iterations;it++) {
    switch (smoother) {
//    fine nodes first, then coarse nodes
      case GAUSS_SEIDEL_RED_BLACK_SMOOTHER: {
//      TRACER_CALL(t,"Level::postSmooth red-black");
        F77_NAME(gauss_seidel_red_black)(fc[0],fc[1],lc[0],lc[1],
          ifirst[0],ifirst[1],ilast[0],ilast[1],1,residual_max,matrix,
          first_residual, solution_increment);
        break;
      }
      case GAUSS_SEIDEL_SMOOTHER: {
//      TRACER_CALL(t,"Level::postSmooth to-fro");
//      last to first
        F77_NAME(gauss_seidel_to_fro)(fc[0],fc[1],lc[0],lc[1],
          ifirst[0],ifirst[1],ilast[0],ilast[1],1,residual_max,matrix,
          first_residual, solution_increment);
        break;
      }
      case RICHARDSON_SMOOTHER:
      default: {
//      TRACER_CALL(t,"Level::postSmooth Richardson");
        double *r=OPERATOR_NEW_BRACKET(double,number_corners);
        F77NAME(richardson)(fc[0],fc[1],lc[0],lc[1],
          ifirst[0],ifirst[1],ilast[0],ilast[1],mu,residual_max,matrix,
          first_residual, solution_increment,r);
        delete [] r;
      }
    }
  }
//
/*
//identity preconditioner, for debugging
  for (int i=0;i<number_corners;i++) {
    solution_increment[i]+=Ax_b[i];
  }
*/
#ifdef DEBUG
//cout << "\tsolution_increment:" << endl;
//F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//  ifirst[0],ifirst[1],ilast[0],ilast[1],solution_increment);
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//applies restriction to compute coarser residual from finer
void Level::restriction() const {
//TRACER_CALL(t,"Level::restriction");
  ASSERT(coarser!=0);
#ifdef DEBUG
//cout << "\tAx_b:" << endl;
//F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],ifirst[0],ifirst[1],
//  ilast[0],ilast[1], Ax_b);
#endif
  switch (restriction_prolongation) {
    case FINITE_ELEMENT_RESTRICTION_PROLONGATION: {
//    TRACER_CALL(t,"Level::restriction fem");
      F77_NAME(restrict_fem)(coarser->fc[0],coarser->fc[1],coarser->lc[0],coarser->lc[1],
        fc[0],fc[1],lc[0],lc[1],Ax_b,coarser->first_residual);
      break;
    }
    case ALGEBRAIC_MULTIGRID_RESTRICTION_PROLONGATION:
    default: {
//    TRACER_CALL(t,"Level::restriction am");
      F77_NAME(restrict_am)(coarser->fe[0],coarser->fe[1],coarser->le[0],coarser->le[1],
        fe[0],fe[1],le[0],le[1],prolongcell,prolongside0,prolongside1,
        Ax_b, coarser->first_residual);
      break;
    }
  }
#ifdef DEBUG
//cout << "\tcoarser first_residual:" << endl;
//F77NAME(rbugcell)(coarser->fc[0],coarser->fc[1],coarser->lc[0],coarser->lc[1],
//  coarser->ifirst[0],coarser->ifirst[1],
//  coarser->ilast[0],coarser->ilast[1],
//  coarser->first_residual);
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//applies prolongation to compute finer solution_increment from coarser
//must be minus adjoint of Level::restriction
void Level::prolongation() const {
//TRACER_CALL(t,"Level::prolongation");
  ASSERT(coarser!=0);
#ifdef DEBUG
//cout << "\tcoarser solution_increment:" << endl;
//F77NAME(rbugcell)(coarser->fc[0],coarser->fc[1],coarser->lc[0],coarser->lc[1],
//  coarser->fc[0],coarser->fc[1],coarser->lc[0],coarser->lc[1],
//  coarser->solution_increment);
#endif
  switch (restriction_prolongation) {
    case FINITE_ELEMENT_RESTRICTION_PROLONGATION: {
//    injection of coarse space into fine
//    TRACER_CALL(t,"Level::prolongation fem");
      F77_NAME(prolong_fem)(coarser->fc[0],coarser->fc[1],coarser->lc[0],coarser->lc[1],
        fc[0],fc[1],lc[0],lc[1],coarser->solution_increment,
        prolongation_vector);
      break;
    }
    case ALGEBRAIC_MULTIGRID_RESTRICTION_PROLONGATION:
    default: {
//    TRACER_CALL(t,"Level::prolongation am");
      F77_NAME(prolong_am)(coarser->fe[0],coarser->fe[1],coarser->le[0],coarser->le[1],
        fe[0],fe[1],le[0],le[1],prolongcell,prolongside0,prolongside1,
        coarser->solution_increment, prolongation_vector);
      break;
    }
  }
  for (int i=0;i<number_corners;i++) {
    solution_increment[i]+=prolongation_vector[i];
  }
#ifdef DEBUG
//cout << "\tsolution_increment:" << endl;
//F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//  ifirst[0],ifirst[1],ilast[0],ilast[1], solution_increment);
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double Level::computeResidual(const double *x,const double *b,double *r)
const {
//TRACER_CALL(t,"Level::computeResidual");
  double residual_max;
#ifdef DEBUG
//cout << "\tx:" << endl;
//F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//  fc[0],fc[1],lc[0],lc[1],x);
#endif
//compute Ax-b:
  F77_NAME(compute_residual)(fc[0],fc[1],lc[0],lc[1],
    ifirst[0],ifirst[1],ilast[0],ilast[1], residual_max,matrix,b,x,r);
#ifdef DEBUG
//cout << "\tr:" << endl;
//F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//  fc[0],fc[1],lc[0],lc[1],r);
#endif
  return residual_max;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::updateResidual() const {
//TRACER_CALL(t,"Level::updateResidual");
#ifdef DEBUG
//cout << "\tfirst_residual:" << endl;
//F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//  ifirst[0],ifirst[1],ilast[0],ilast[1],first_residual);
//cout << "\tsolution_increment:" << endl;
//F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//  ifirst[0],ifirst[1],ilast[0],ilast[1],solution_increment);
#endif
  double residual_max;
  F77_NAME(compute_residual)(fc[0],fc[1],lc[0],lc[1],
    ifirst[0],ifirst[1],ilast[0],ilast[1], residual_max,
    matrix,first_residual,solution_increment,Ax_b);
  for (int i=0;i<number_corners;i++) Ax_b[i]=-Ax_b[i];
#ifdef DEBUG
//cout << "\tAx_b:" << endl;
//F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//  ifirst[0],ifirst[1],ilast[0],ilast[1],Ax_b);
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//the following is usually called on the coarsest level only
void Level::solve() const {
//TRACER_CALL(t,"Level::solve");
  double residual_max=numeric_limits<double>::infinity();
  double *ap=OPERATOR_NEW_BRACKET(double,number_corners);
  double *cg_direction=OPERATOR_NEW_BRACKET(double,number_corners);
  double gamma=0.;
  double rhs_max=0.;
#ifdef DEBUG
//cout << "\tfirst_residual:" << endl;
//F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//  ifirst[0],ifirst[1],ilast[0],ilast[1],first_residual);
//for (int i=-4;i<=4;i++) {
//  cout << "\tmatrix[" << i << "] = " << endl;
//  F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//    ifirst[0],ifirst[1],ilast[0],ilast[1],matrix+(i+4)*number_corners);
//}
#endif
  for (int i=0;i<number_corners;i++) {
    ap[i]=0.;
    solution_increment[i]=0.;
    Ax_b[i]=first_residual[i];
    cg_direction[i]=Ax_b[i];
    gamma+=cg_direction[i]*Ax_b[i];
    rhs_max=max(rhs_max,abs(Ax_b[i]));
  }
  if (gamma>0.) {
    for (int it=0;it<number_corners && residual_max>1.e-16*rhs_max;it++)
    {
      F77_NAME(matrix_multiply)(fc[0],fc[1],lc[0],lc[1],
        ifirst[0],ifirst[1],ilast[0],ilast[1], matrix,cg_direction, ap);
#ifdef DEBUG
//    cout << "\tap:" << endl;
//    F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//      ifirst[0],ifirst[1],ilast[0],ilast[1],ap);
#endif
      double dot=F77NAME(ddot)(number_corners,cg_direction,1,ap,1);
      ASSERT(dot>0.);
      double alpha=gamma/dot;
      residual_max=0.;
      for (int i=0;i<number_corners;i++) {
        solution_increment[i]+=cg_direction[i]*alpha;
        Ax_b[i]-=ap[i]*alpha;
        residual_max=max(residual_max,abs(Ax_b[i]));
      }
#ifdef DEBUG
//    cout << "\tsolution_increment:" << endl;
//    F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//      ifirst[0],ifirst[1],ilast[0],ilast[1],solution_increment);
//    cout << "\tAx_b:" << endl;
//    F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//      fc[0],fc[1],lc[0],lc[1],Ax_b);
//    cout << "\tresidual_max = " << residual_max << endl;
#endif
      double delta=F77NAME(ddot)(number_corners,Ax_b,1,Ax_b,1);
      double beta=delta/gamma;
      for (int i=0;i<number_corners;i++) {
        cg_direction[i]=Ax_b[i]+cg_direction[i]*beta;
      }
#ifdef DEBUG
//    cout << "\tcg_direction:" << endl;
//    F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//      ifirst[0],ifirst[1],ilast[0],ilast[1],cg_direction);
#endif
      gamma=delta;
    }
  }
  for (int i=0;i<number_corners;i++) Ax_b[i]=0.;
  delete [] ap;
  delete [] cg_direction;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//one step of the multigrid algorithm
void Level::multigridStep(double *ax_b, double *d) const {
//TRACER_CALL(t,"Level::multigridStep");
  if (finer==0) {
    memcpy(first_residual,ax_b,number_corners*sizeof(double));
#ifdef DEBUG
//  cout << "\tfirst_residual:" << endl;
//  F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//    ifirst[0],ifirst[1],ilast[0],ilast[1],first_residual);
#endif
  }
  for (int i=0;i<number_corners;i++) solution_increment[i]=0;
  if (coarser==0) {
    solve();
#ifdef DEBUG
//  cout << "\tsolution_increment:" << endl;
//  F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//    ifirst[0],ifirst[1],ilast[0],ilast[1],solution_increment);
#endif
  } else {
    preSmooth();
#ifdef DEBUG
//  cout << "\tsolution_increment:" << endl;
//  F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//    ifirst[0],ifirst[1],ilast[0],ilast[1],solution_increment);
#endif
    updateResidual();
#ifdef DEBUG
//  cout << "\tAx_b:" << endl;
//  F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//    ifirst[0],ifirst[1],ilast[0],ilast[1],Ax_b);
#endif
    restriction();
#ifdef DEBUG
//  cout << "\tcoarser->first_residual:" << endl;
//  F77NAME(rbugcell)(coarser->fc[0],coarser->fc[1],coarser->lc[0],coarser->lc[1],
//    coarser->ifirst[0],coarser->ifirst[1],coarser->ilast[0],coarser->ilast[1],
//    coarser->first_residual);
#endif
    coarser->multigridStep(coarser->Ax_b,coarser->solution_increment);
    prolongation();
#ifdef DEBUG
//  cout << "\tsolution_increment:" << endl;
//  F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//    ifirst[0],ifirst[1],ilast[0],ilast[1],solution_increment);
#endif
    updateResidual();
#ifdef DEBUG
//  cout << "\tAx_b:" << endl;
//  F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//    ifirst[0],ifirst[1],ilast[0],ilast[1],Ax_b);
#endif
    postSmooth();
#ifdef DEBUG
//  cout << "\tsolution_increment:" << endl;
//  F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//    ifirst[0],ifirst[1],ilast[0],ilast[1],solution_increment);
#endif
  }
  if (finer==0) {
    memcpy(d,solution_increment,number_corners*sizeof(double));
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::checkProlongationAndRestrictionRandom() const {
  TRACER_CALL(t,"Level::checkProlongationAndRestrictionRandom");
  for (int i=0;i<number_corners;i++) {
    Ax_b[i]=drand48();
    solution_increment[i]=0.;
  }
  for (int i=0;i<coarser->number_corners;i++) {
    coarser->solution_increment[i]=drand48();
    coarser->first_residual[i]=0.;
  }
  restriction(); // fine Ax_b -> coarse first_residual
  prolongation(); // coarse solution_inc -> fine solution_inc
  double fine_sum=0.;
  for (int i=0;i<number_corners;i++) {
    fine_sum+=Ax_b[i]*solution_increment[i];
  }
  double coarse_sum=0.;
  for (int i=0;i<coarser->number_corners;i++) {
    coarse_sum+=
      coarser->first_residual[i]*coarser->solution_increment[i];
  }
  cout << "\tcoarse_sum,fine_sum = " << coarse_sum << " " 
       << fine_sum << " " 
       << abs(coarse_sum-fine_sum)
         /max(abs(coarse_sum),abs(fine_sum)) << endl;
  ASSERT(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::checkProlongationAndRestriction() const {
  TRACER_CALL(t,"Level::checkProlongationAndRestriction");
#ifdef DEBUG
//
//for (int i=-4;i<=4;i++) {
//  cout << "\tmatrix[" << i << "] = " << endl;
//  F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//    ifirst[0],ifirst[1],ilast[0],ilast[1],matrix+(i+4)*number_corners);
//}
//
#endif
  for (int i=0;i<number_corners;i++) 
//int i=7;
//F77NAME(const).debug_on=1;
  {
    for (int j=0;j<number_corners;j++) Ax_b[j]=0.;
    Ax_b[i]=1.;
    for (int J=0;J<coarser->number_corners;J++) {
      coarser->first_residual[J]=0.;
    }
    restriction();
#ifdef DEBUG
/*
    cout << "\tAx_b:" << endl;
    F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
      ifirst[0],ifirst[1],ilast[0],ilast[1],Ax_b);
    cout << "\tcoarser->first_residual:" << endl;
    F77NAME(rbugcell)(coarser->fc[0],coarser->fc[1],coarser->lc[0],coarser->lc[1],
      coarser->ifirst[0],coarser->ifirst[1],coarser->ilast[0],coarser->ilast[1],
      coarser->first_residual);
*/
#endif
    for (int I=0;I<coarser->number_corners;I++) 
//  int I=4;
    {
      for (int J=0;J<coarser->number_corners;J++) {
        coarser->solution_increment[J]=0.;
      }
      coarser->solution_increment[I]=1.;
      for (int j=0;j<number_corners;j++) solution_increment[j]=0.;
      prolongation();
#ifdef DEBUG
/*
      cout << "\tcoarser->solution_increment:" << endl;
      F77NAME(rbugcell)(coarser->fc[0],coarser->fc[1],coarser->lc[0],coarser->lc[1],
        coarser->ifirst[0],coarser->ifirst[1],coarser->ilast[0],coarser->ilast[1],
        coarser->solution_increment);
      cout << "\tsolution_increment:" << endl;
      F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
        ifirst[0],ifirst[1],ilast[0],ilast[1],solution_increment);
*/
#endif
      if (abs(coarser->first_residual[I]-solution_increment[i])>
      1.e-12*max(abs(coarser->first_residual[I]),
      abs(solution_increment[i]))) {
        cout << "\tR,P[" << I << "," << i << "] = " 
             << coarser->first_residual[I] << " " 
             << solution_increment[i] << " "
             << abs(coarser->first_residual[I]-solution_increment[i])
               /max(abs(coarser->first_residual[I]),
                    abs(solution_increment[i]))
             << endl;
      }
    }
  }
  ASSERT(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::checkMatrixSymmetryRandom() const {
  TRACER_CALL(t,"Level::checkMatrixSymmetryRandom");
  double *b=OPERATOR_NEW_BRACKET(double,number_corners);
  double *x=OPERATOR_NEW_BRACKET(double,number_corners);
  double *y=OPERATOR_NEW_BRACKET(double,number_corners);
  double *ax=OPERATOR_NEW_BRACKET(double,number_corners);
  double *ay=OPERATOR_NEW_BRACKET(double,number_corners);
  for (int i=0;i<number_corners;i++) {
    b[i]=0.;
    int i1=i/(n[0]+1);
    int i0=i%(n[0]+1);
    if (i0>0 && i0<n[0] && i1>0 && i1<n[1]) {
      x[i]=drand48();
      y[i]=drand48();
    } else {
      x[i]=0.;
      y[i]=0.;
    }
  }
  computeResidual(x,b,ax); 
  computeResidual(y,b,ay); 
  double sum1=0.,sum2=0.;
  for (int i=0;i<number_corners;i++) {
    sum1+=y[i]*ax[i];
    sum2+=ay[i]*x[i];
  }
  cout << "\tsum1,sum2 = " << sum1 << " " << sum2 << " "
       << abs(sum1-sum2)/max(abs(sum1),abs(sum2)) << endl;
  delete [] b;
  delete [] x;
  delete [] y;
  delete [] ax;
  delete [] ay;
  ASSERT(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::checkMatrixSymmetry() const {
  TRACER_CALL(t,"Level::checkMatrixSymmetry");
  double *b=OPERATOR_NEW_BRACKET(double,number_corners);
  double *x=OPERATOR_NEW_BRACKET(double,number_corners);
  double *y=OPERATOR_NEW_BRACKET(double,number_corners);
  double *ax=OPERATOR_NEW_BRACKET(double,number_corners);
  double *ay=OPERATOR_NEW_BRACKET(double,number_corners);
  for (int i=0;i<number_corners;i++) {
    b[i]=0.;
  }
  for (int j=0;j<number_corners;j++) 
//int j=16;
  {
    int j1=j/(n[0]+1);
    int j0=j%(n[0]+1);
    if (j0>0 && j0<n[0] && j1>0 && j1<n[1]) {
      for (int k=0;k<number_corners;k++) y[k]=0;
      y[j]=1.;
      computeResidual(y,b,ay); 
      for (int i=0;i<number_corners;i++) 
//    int i=17;
      {
        int i1=i/(n[0]+1);
        int i0=i%(n[0]+1);
        if (i0>0 && i0<n[0] && i1>0 && i1<n[1]) {
          for (int k=0;k<number_corners;k++) x[k]=0;
          x[i]=1.;
          computeResidual(x,b,ax); 
          if (abs(ay[i]-ax[j])>1.e-12*max(abs(ay[i]),abs(ax[j]))) {
            cout << "\tA[" << i << "," << j << "] = "
                 << ay[i] << " " << ax[j] << endl;
          }
        }
      }
    }
  }
  delete [] b;
  delete [] x;
  delete [] y;
  delete [] ax;
  delete [] ay;
  ASSERT(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::checkVCycleSymmetryRandom() const {
  TRACER_CALL(t,"Level::checkVCycleSymmetryRandom");
  double *r=OPERATOR_NEW_BRACKET(double,number_corners);
  double *s=OPERATOR_NEW_BRACKET(double,number_corners);
  double *vr=OPERATOR_NEW_BRACKET(double,number_corners);
  double *vs=OPERATOR_NEW_BRACKET(double,number_corners);
  for (int i=0;i<number_corners;i++) {
    int i1=i/(n[0]+1);
    int i0=i%(n[0]+1);
    if (i0>0 && i0<n[0] && i1>0 && i1<n[1]) {
      r[i]=drand48();
      s[i]=drand48();
    } else {
      r[i]=0.;
      s[i]=0.;
    }
  }
#ifdef DEBUG
  cout << "\tr:" << endl;
  F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
    ifirst[0],ifirst[1],ilast[0],ilast[1],r);
  cout << "\ts:" << endl;
  F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
    ifirst[0],ifirst[1],ilast[0],ilast[1],s);
#endif

  memcpy(first_residual,r,number_corners*sizeof(double));
  for (int i=0;i<number_corners;i++) solution_increment[i]=0;
  if (coarser==0) {
    solve();
  } else {
    preSmooth();
    updateResidual();
    restriction();
    coarser->multigridStep(coarser->Ax_b,coarser->solution_increment);
    prolongation();
    updateResidual();
    postSmooth();
  }
  memcpy(vr,solution_increment,number_corners*sizeof(double));

  memcpy(first_residual,s,number_corners*sizeof(double));
  for (int i=0;i<number_corners;i++) solution_increment[i]=0;
  if (coarser==0) {
    solve();
  } else {
    preSmooth();
    updateResidual();
    restriction();
    coarser->multigridStep(coarser->Ax_b,coarser->solution_increment);
    prolongation();
    updateResidual();
    postSmooth();
  }
  memcpy(vs,solution_increment,number_corners*sizeof(double));

#ifdef DEBUG
  cout << "\tvr:" << endl;
  F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
    ifirst[0],ifirst[1],ilast[0],ilast[1],vr);
  cout << "\tvs:" << endl;
  F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
    ifirst[0],ifirst[1],ilast[0],ilast[1],vs);
#endif
  double sum1=0.,sum2=0.;
  for (int i=0;i<number_corners;i++) {
    sum1+=s[i]*vr[i];
    sum2+=vs[i]*r[i];
  }
  cout << "\tsum1,sum2 = " << sum1 << " " << sum2 << " "
       << abs(sum1-sum2)/max(abs(sum1),abs(sum2)) << endl;
  delete [] r;
  delete [] s;
  delete [] vr;
  delete [] vs;
  ASSERT(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::checkVCycleSymmetry() const {
  TRACER_CALL(t,"Level::checkVCycleSymmetry");
  double *r=OPERATOR_NEW_BRACKET(double,number_corners);
  double *s=OPERATOR_NEW_BRACKET(double,number_corners);
  double *vr=OPERATOR_NEW_BRACKET(double,number_corners);
  double *vs=OPERATOR_NEW_BRACKET(double,number_corners);
  for (int j=0;j<number_corners;j++) {
    int j1=j/(n[0]+1);
    int j0=j%(n[0]+1);
    if (j0>0 && j0<n[0] && j1>0 && j1<n[1]) {
      for (int k=0;k<number_corners;k++) s[k]=0;
      s[j]=1.;
      multigridStep(s,vs);
      for (int i=0;i<number_corners;i++) {
        int i1=i/(n[0]+1);
        int i0=i%(n[0]+1);
        if (i0>0 && i0<n[0] && i1>0 && i1<n[1]) {
          for (int k=0;k<number_corners;k++) r[k]=0;
          r[i]=1.;
          multigridStep(r,vr);
          if (abs(vs[i]-vr[j])>1.e-12*max(abs(vs[i]),abs(vr[j]))) {
            cout << "\tV[" << i << "," << j << "] = "
                 << vs[i] << " " << vr[j] << endl;
          }
        }
      }
    }
  }
  delete [] r;
  delete [] s;
  delete [] vr;
  delete [] vs;
  ASSERT(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::checkSmoother() const {
  TRACER_CALL(t,"Level::checkSmoother");
  double *si=OPERATOR_NEW_BRACKET(double,number_corners);
  for (int j=0;j<number_corners;j++) 
//int j=6;
  {
    for (int k=0;k<number_corners;k++) {
      first_residual[k]=0.;
      solution_increment[k]=0.;
    }
    first_residual[j]=1.;
    preSmooth();
#ifdef DEBUG
//  cout << "\tfirst_residual:" << endl;
//  F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//    ifirst[0],ifirst[1],ilast[0],ilast[1],first_residual);
//  cout << "\tsolution_increment:" << endl;
//  F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//    ifirst[0],ifirst[1],ilast[0],ilast[1],solution_increment);
#endif
    memcpy(si,solution_increment,number_corners*sizeof(double));
    for (int i=0;i<number_corners;i++) 
//  int i=7;
    {
      for (int k=0;k<number_corners;k++) {
        first_residual[k]=0.;
        solution_increment[k]=0.;
      }
      first_residual[i]=1.;
      postSmooth();
#ifdef DEBUG
//    cout << "\tfirst_residual:" << endl;
//    F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//      ifirst[0],ifirst[1],ilast[0],ilast[1],first_residual);
//    cout << "\tsolution_increment:" << endl;
//    F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//      ifirst[0],ifirst[1],ilast[0],ilast[1],solution_increment);
#endif
      if (abs(si[i]-solution_increment[j])>
      1.e-12*max(abs(si[i]),abs(solution_increment[j]))) {
        cout << "\tsi[" << i << "," << j << "] = " 
             << si[i] << " " << solution_increment[j] << endl;
      }
    }
  }
  delete [] si;
  ASSERT(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::checkCoarseGridProjection() const {
  TRACER_CALL(t,"Level::checkCoarseGridProjection");
  for (int I=0;I<coarser->number_corners;I++) 
//int I=5;
  {
    for (int k=0;k<number_corners;k++) {
      first_residual[k]=0.;
      solution_increment[k]=0.;
    }
    for (int K=0;K<coarser->number_corners;K++) {
      coarser->solution_increment[K]=0.;
    }
    coarser->solution_increment[I]=1.;
#ifdef DEBUG
//  cout << "\tcoarser solution_increment:" << endl;
//  F77NAME(rbugcell)(coarser->fc[0],coarser->fc[1],coarser->lc[0],coarser->lc[1],
//    coarser->ifirst[0],coarser->ifirst[1],coarser->ilast[0],coarser->ilast[1],
//    coarser->solution_increment);
#endif
    prolongation();
#ifdef DEBUG
//  cout << "\tsolution_increment:" << endl;
//  F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//    ifirst[0],ifirst[1],ilast[0],ilast[1],solution_increment);
#endif
    updateResidual();
#ifdef DEBUG
//  cout << "\tAx_b:" << endl;
//  F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//    ifirst[0],ifirst[1],ilast[0],ilast[1],Ax_b);
#endif
    restriction();
#ifdef DEBUG
//  cout << "\tcoarser first_residual:" << endl;
//  F77NAME(rbugcell)(coarser->fc[0],coarser->fc[1],coarser->lc[0],coarser->lc[1],
//    coarser->ifirst[0],coarser->ifirst[1],coarser->ilast[0],coarser->ilast[1],
//    coarser->first_residual);
#endif
    coarser->multigridStep(coarser->Ax_b,coarser->solution_increment);
#ifdef DEBUG
//  cout << "\tcoarser solution_increment:" << endl;
//  F77NAME(rbugcell)(coarser->fc[0],coarser->fc[1],coarser->lc[0],coarser->lc[1],
//    coarser->ifirst[0],coarser->ifirst[1],coarser->ilast[0],coarser->ilast[1],
//    coarser->solution_increment);
#endif
    prolongation();
#ifdef DEBUG
//  cout << "\tsolution_increment:" << endl;
//  F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//    ifirst[0],ifirst[1],ilast[0],ilast[1],solution_increment);
#endif
    double si_max,si_min;
    F77NAME(rcmpcell)(fc[0],fc[1],lc[0],lc[1],
      ifirst[0],ifirst[1],ilast[0],ilast[1],solution_increment,si_max,si_min);
    cout << "\tnullspace[" << I << "] = " << si_min << " " << si_max
         << endl;
  }

/*
  double *restriction_matrix=
    OPERATOR_NEW_BRACKET(double,coarser->number_corners*number_corners);
  double *singular_values=
    OPERATOR_NEW_BRACKET(double,coarser->number_corners);
  double *right_singular_vectors=
    OPERATOR_NEW_BRACKET(double,number_corners*number_corners);
  int lwork=max(3*coarser->number_corners+number_corners,
                5*coarser->number_corners-4);
  double *work=OPERATOR_NEW_BRACKET(double,lwork);
  for (int j=0;j<number_corners;j++) {
    for (int i=0;i<number_corners;i++) Ax_b[i]=0.;
    Ax_b[j]=1.;
    restriction();
    for (int I=0;I<coarser->number_corners;I++) {
      restriction_matrix[I+j*coarser->number_corners]=
        coarser->first_residual[I];
    }
  }
  int info=0;
  F77NAME(dgesvd)('N','A',coarser->number_corners,number_corners,
    restriction_matrix,coarser->number_corners,singular_values,
    0,1,right_singular_vectors,number_corners,work,lwork,info);
//for (int j=coarser->number_corners;j<number_corners;j++) 
  int j=11;
  {
    for (int i=0;i<number_corners;i++) {
      first_residual[i]=right_singular_vectors[j+i*number_corners];
    }
#ifdef DEBUG
//  cout << "\tfirst_residual:" << endl;
//  F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//    ifirst[0],ifirst[1],ilast[0],ilast[1],first_residual);
#endif
    solve();
#ifdef DEBUG
//  cout << "\tsolution_increment:" << endl;
//  F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//    ifirst[0],ifirst[1],ilast[0],ilast[1],solution_increment);
#endif
    for (int i=0;i<number_corners;i++) first_residual[i]=0.;
    updateResidual();
#ifdef DEBUG
//  cout << "\tAx_b:" << endl;
//  F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//    ifirst[0],ifirst[1],ilast[0],ilast[1],Ax_b);
#endif
    for (int i=0;i<number_corners;i++) solution_increment[i]=0.;
    restriction();
#ifdef DEBUG
//  cout << "\tcoarser first_residual:" << endl;
//  F77NAME(rbugcell)(coarser->fc[0],coarser->fc[1],coarser->lc[0],coarser->lc[1],
//    coarser->ifirst[0],coarser->ifirst[1],coarser->ilast[0],coarser->ilast[1],
//    coarser->first_residual);
#endif
    coarser->multigridStep(coarser->Ax_b,coarser->solution_increment);
#ifdef DEBUG
//  cout << "\tcoarser solution_increment:" << endl;
//  F77NAME(rbugcell)(coarser->fc[0],coarser->fc[1],coarser->lc[0],coarser->lc[1],
//    coarser->ifirst[0],coarser->ifirst[1],coarser->ilast[0],coarser->ilast[1],
//    coarser->solution_increment);
#endif
    prolongation();
#ifdef DEBUG
//  cout << "\tsolution_increment:" << endl;
//  F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//    ifirst[0],ifirst[1],ilast[0],ilast[1],solution_increment);
#endif
    double si_max,si_min;
    F77NAME(rcmpcell)(fc[0],fc[1],lc[0],lc[1],
      ifirst[0],ifirst[1],ilast[0],ilast[1],solution_increment,si_max,si_min);
    cout << "\trange[" << j << "] = " << si_min << " " << si_max
         << endl;
  }
  delete [] restriction_matrix;
  delete [] singular_values;
  delete [] right_singular_vectors;
*/
  ASSERT(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::checkPointSourceRandom() const {
  TRACER_CALL(t,"Level::checkPointSourceRandom");
  double *b=OPERATOR_NEW_BRACKET(double,number_corners);
  double *d=OPERATOR_NEW_BRACKET(double,number_corners);
  double *r=OPERATOR_NEW_BRACKET(double,number_corners);
  double *x=OPERATOR_NEW_BRACKET(double,number_corners);
  for (int k=0;k<100;k++) {
    for (int i=0;i<number_corners;i++) {
      int i0=i%(n[0]+1);
      int i1=i/(n[0]+1);
      x[i]=(i0>0 && i0<n[0] && i1>0 && i1<n[1] ? 
        static_cast<double>(rand())/static_cast<double>(RAND_MAX) : 0.);
      b[i]=0.;
      d[i]=0.;
      r[i]=0.;
    }
    computeResidual(x,b,r); 
    multigridStep(r,d);
    for (int i=0;i<number_corners;i++) x[i]-=d[i];
    double si_max,si_min;
    F77NAME(rcmpcell)(fc[0],fc[1],lc[0],lc[1],
      ifirst[0],ifirst[1],ilast[0],ilast[1],x,si_max,si_min);
    cout << "\tpoint error[" << k << "] = " << max(-si_min,si_max)
         << endl;
  }
  delete [] b;
  delete [] d;
  delete [] r;
  delete [] x;
  ASSERT(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::checkPointSource() const {
  TRACER_CALL(t,"Level::checkPointSource");
  double *b=OPERATOR_NEW_BRACKET(double,number_corners);
  double *d=OPERATOR_NEW_BRACKET(double,number_corners);
  double *r=OPERATOR_NEW_BRACKET(double,number_corners);
  double *x=OPERATOR_NEW_BRACKET(double,number_corners);
  for (int i=0;i<number_corners;i++) 
//int i=4;
  {
    int i0=i%(n[0]+1);
    int i1=i/(n[0]+1);
    if (i0>0 && i0<n[0] && i1>0 && i1<n[1]) {
      for (int k=0;k<number_corners;k++) {
        b[k]=0.;
        r[k]=0.;
        x[k]=0.;
      }
      x[i]=1.;
      computeResidual(x,b,r); 
      multigridStep(r,d);
      for (int k=0;k<number_corners;k++) x[k]-=d[k];
      double si_max,si_min;
      F77NAME(rcmpcell)(fc[0],fc[1],lc[0],lc[1],
        ifirst[0],ifirst[1],ilast[0],ilast[1],x,si_max,si_min);
      cout << "\tpoint error[" << i << "] = " << max(-si_min,si_max)
           << endl;
    }
  }
  delete [] b;
  delete [] d;
  delete [] r;
  delete [] x;
  ASSERT(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::checkSmootherUpdate() const {
  TRACER_CALL(t,"Level::checkSmootherUpdate");
  double *si=OPERATOR_NEW_BRACKET(double,number_corners);
  for (int j=0;j<number_corners;j++) 
//int j=6;
  {
    int j1=j/(n[0]+1);
    int j0=j%(n[0]+1);
    if (j0>0 && j0<n[0] && j1>0 && j1<n[1]) {
      for (int k=0;k<number_corners;k++) {
        first_residual[k]=0.;
        solution_increment[k]=0.;
      }
      first_residual[j]=1.;
      preSmooth();
      updateResidual();
      memcpy(si,Ax_b,number_corners*sizeof(double));
      for (int i=0;i<number_corners;i++) 
//    int i=6;
      {
        int i1=i/(n[0]+1);
        int i0=i%(n[0]+1);
        if (i0>0 && i0<n[0] && i1>0 && i1<n[1]) {
          for (int k=0;k<number_corners;k++) {
            first_residual[k]=0.;
            solution_increment[k]=0.;
          }
          solution_increment[i]=1.;
          updateResidual();
          postSmooth();
          if (abs(si[i]-solution_increment[j])>
          1.e-12*max(abs(si[i]),abs(solution_increment[j]))) {
            cout << "\tsi[" << i << "," << j << "] = " 
                 << si[i] << " " << solution_increment[j] << endl;
          }
        }
      }
    }
  }
  delete [] si;
  ASSERT(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::checkSmootherUpdateRestrict() const {
  TRACER_CALL(t,"Level::checkSmootherUpdateRestrict");
  ASSERT(coarser!=0);
  double *coarse_answer=
    OPERATOR_NEW_BRACKET(double,coarser->number_corners);
  for (int i=0;i<number_corners;i++) 
//int i=12;
  {
    int i1=i/(n[0]+1);
    int i0=i%(n[0]+1);
    if (i0>0 && i0<n[0] && i1>0 && i1<n[1]) {
      for (int k=0;k<number_corners;k++) solution_increment[k]=0.;
      for (int k=0;k<number_corners;k++) first_residual[k]=0.;
      first_residual[i]=1.;
#ifdef DEBUG
//    cout << "\tfirst_residual:" << endl;
//    F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//      ifirst[0],ifirst[1],ilast[0],ilast[1],first_residual);
#endif
//    preSmooth(); // d_f = d_f + S_f r_f^0 = S_f r_f^0
#ifdef DEBUG
//    cout << "\tsolution_increment:" << endl;
//    F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//      ifirst[0],ifirst[1],ilast[0],ilast[1],solution_increment);
#endif
      updateResidual(); // r_f = ( I - A_f S_f ) r_f^0
#ifdef DEBUG
//    cout << "\tAx_b:" << endl;
//    F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//      ifirst[0],ifirst[1],ilast[0],ilast[1],Ax_b);
#endif
      restriction(); // r_c = R ( I - A_f S_f ) r_f^0
      memcpy(coarse_answer,coarser->first_residual,
             coarser->number_corners*sizeof(double));
#ifdef DEBUG
//    cout << "\tcoarser first_residual:" << endl;
//    F77NAME(rbugcell)(coarser->fc[0],coarser->fc[1],coarser->lc[0],coarser->lc[1],
//      coarser->ifirst[0],coarser->ifirst[1],coarser->ilast[0],coarser->ilast[1],
//      coarser->first_residual);
#endif
//    coarser->multigridStep(coarser->Ax_b,coarser->solution_increment);
//    memcpy(coarse_answer,coarser->solution_increment,
//           coarser->number_corners*sizeof(double));
        // d_c = V_c R ( I - A_f S_f ) r_f^0
#ifdef DEBUG
//    cout << "\tcoarser solution_increment:" << endl;
//    F77NAME(rbugcell)(coarser->fc[0],coarser->fc[1],coarser->lc[0],coarser->lc[1],
//      coarser->ifirst[0],coarser->ifirst[1],coarser->ilast[0],coarser->ilast[1],
//      coarser->solution_increment);
#endif

      for (int I=0;I<coarser->number_corners;I++) 
//    int I=4;
      {
        int I1=I/(coarser->n[0]+1);
        int I0=I%(coarser->n[0]+1);
        if (I0>0 && I0<coarser->n[0] && I1>0 && I1<coarser->n[1]) {
          for (int k=0;k<number_corners;k++) {
            first_residual[k]=0.;
            solution_increment[k]=0.;
          }
          for (int K=0;K<coarser->number_corners;K++) {
            coarser->solution_increment[K]=0.; //start with prolongation
//          coarser->first_residual[K]=0.; //start with multigridStep
          }
          coarser->solution_increment[I]=1.;
//        coarser->first_residual[I]=1.;
#ifdef DEBUG
//        cout << "\tcoarser first_residual:" << endl;
//        F77NAME(rbugcell)(coarser->fc[0],coarser->fc[1],coarser->lc[0],coarser->lc[1],
//          coarser->ifirst[0],coarser->ifirst[1],coarser->ilast[0],coarser->ilast[1],
//          coarser->first_residual);
#endif
//        coarser->
//          multigridStep(coarser->Ax_b,coarser->solution_increment);
#ifdef DEBUG
//        cout << "\tcoarser solution_increment:" << endl;
//        F77NAME(rbugcell)(coarser->fc[0]coarser->fc[1],coarser->lc[0],coarser->lc[1],
//          coarser->ifirst[0],coarser->ifirst[1],coarser->ilast[0],coarser->ilast[1],
//          coarser->solution_increment);
#endif
          prolongation(); // d_f = d_f + P d_c = P d_c
#ifdef DEBUG
//        cout << "\tsolution_increment:" << endl;
//        F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//          ifirst[0],ifirst[1],ilast[0],ilast[1],solution_increment);
#endif
          updateResidual(); // r_f = r_f^0 - A_f d_f = - A_f P d_c
#ifdef DEBUG
//        cout << "\tAx_b:" << endl;
//        F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//          ifirst[0],ifirst[1],ilast[0],ilast[1],Ax_b);
#endif
//        postSmooth();//d_f = d_f + S_f^T r_f = ( I - S_f^T A_f ) P d_c
#ifdef DEBUG
//        cout << "\tsolution_increment:" << endl;
//        F77NAME(rbugcell)(fc[0],fc[1],lc[0],lc[1],
//          ifirst[0],ifirst[1],ilast[0],ilast[1],solution_increment);
#endif

          if (abs(4.*coarse_answer[I]-solution_increment[i])>1.e-12*
          max(abs(4.*coarse_answer[I]),abs(solution_increment[i]))) {
            cout << "\tcoarse_answer,solution_increment[" 
                 << I << "," << i << "] = " << 4.*coarse_answer[I] 
                 << " " << solution_increment[i] << endl;
          }
        }
      }
    }
  }
  delete [] coarse_answer;
  ASSERT(0);
}
