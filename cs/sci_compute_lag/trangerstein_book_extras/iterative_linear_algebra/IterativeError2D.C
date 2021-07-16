#include <float.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdlib.h>
#include <math.h>
#include "Arch.H"
#include "AreaGraphTool.H"
#include "DDecl.H"
#include "Debug.H"
#include "Fort2D.H"
#include "GUIEnumInputParameter.H"
#include "GUIVirtualInput.H"
#include "InputParameter.H"
#include "Level2D.H"
#include "Palette.H"
#include "SetTraps.H"
#include "TimedObject.H"
#include "Tracer.H"
#include "VGT.H"
#include "XYGraphTool.H"
#include <unistd.h>
#ifdef USE_GTK
#include <gtk/gtkmain.h>
#include "GTKGUI.H"
#include "GTKGUIVirtualInput.H"
#include "GTKWindow.H"
#else
#include "GUI.H"
#include "GUIVirtualInput.H"
#endif

struct const_common {
  double e,pi,rt2, half,one,two,zero;
  int debug_on;
};
extern const_common F77NAME(const);

struct machine_common {
  double roundoff,small,huge,undefind;
};
extern machine_common F77NAME(machine);

clock_t start;
struct tms usage;

GUI_INPUT_PARAMETER_LIST_TYPE *param_list=0;
bool skip_gui=FALSE;
int ncells[2]={8,8};
int nsteps=10;
//int nbands=2;
int smoother_iterations=1;
double decay=1.;
double mu=4.;
double omega=1.;
double tol=1.e-12;
double winsize=0.5;
int ncontours=5;
char* display_name=0;

bool run_richardson=false;
bool run_jacobi=false;
bool run_gauss_seidel_red_black=false;
bool run_gauss_seidel_to_fro=false;
bool run_sor=false;
//bool run_incomplete_factorization=false;
bool run_conjugate_gradients=false;
bool run_multigrid=true;

enum EQUATION{LAPLACE_EQUATION,HEAT_EQUATION};
EQUATION equation=LAPLACE_EQUATION;
int iequation=equation;
static const char *equation_name[2]={"laplace_equation","heat_equation"};

enum DIFFUSION{CONSTANT_DIFFUSION,RANDOM_DIFFUSION};
//DIFFUSION diffusion=CONSTANT_DIFFUSION;
DIFFUSION diffusion=RANDOM_DIFFUSION;
int idiffusion=diffusion;
static const char *diffusion_name[2]={
  "constant_diffusion","random_diffusion"
};

enum ITERATION{RICHARDSON_ITERATION,JACOBI_OMEGA_ITERATION,
  GAUSS_SEIDEL_RED_BLACK_ITERATION,GAUSS_SEIDEL_TO_FRO_ITERATION,
  SOR_ITERATION, // INCOMPLETE_FACTORIZATION_ITERATION,
  CONJUGATE_GRADIENTS_ITERATION,MULTIGRID_ITERATION};
ITERATION iteration=RICHARDSON_ITERATION;
int iiteration=iteration;
static const char *iteration_name[7 /* 8 */]={
  "richardson","jacobi_omega","gauss_seidel_red_black",
  "gauss_seidel_to_fro","sor", // "incomplete_factorization",
  "conjugate_gradients","multigrid"
};

enum PRECONDITIONER{IDENTITY_PRECONDITIONER,JACOBI_PRECONDITIONER,
  BLOCK_JACOBI_PRECONDITIONER, // INCOMPLETE_CHOLESKY_PRECONDITIONER,
  MULTIGRID_PRECONDITIONER};
PRECONDITIONER preconditioner=IDENTITY_PRECONDITIONER;
int ipreconditioner=preconditioner;
static const char *preconditioner_name[4 /* 5 */ ]={
  "identity","jacobi","block_jacobi", // "incomplete cholesky",
  "multigrid"
};

Level::SMOOTHER smoother=Level::GAUSS_SEIDEL_RED_BLACK_SMOOTHER;
int ismoother=smoother;

//algebraic multigrid prolongation not debugged
Level::RESTRICTION_PROLONGATION restriction_prolongation=
  Level::ALGEBRAIC_MULTIGRID_RESTRICTION_PROLONGATION; 
//Level::FINITE_ELEMENT_RESTRICTION_PROLONGATION; 
int irestriction_prolongation=restriction_prolongation;

static const char *color_name[11]={
  "red","green","blue","cyan","yellow","magenta","brown","purple",
  "grey","pink","tan"
};

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void checkMainInput() {
//TRACER_CALL(tr,"checkMainInput");
//iteration=static_cast<ITERATION>(iiteration);
  if (run_richardson) iteration=RICHARDSON_ITERATION;
  else if (run_jacobi) iteration=JACOBI_OMEGA_ITERATION;
  else if (run_gauss_seidel_red_black) {
    iteration=GAUSS_SEIDEL_RED_BLACK_ITERATION;
  } else if (run_gauss_seidel_to_fro) {
    iteration=GAUSS_SEIDEL_TO_FRO_ITERATION;
  } else if (run_sor) iteration=SOR_ITERATION;
  else if (run_conjugate_gradients) {
    iteration=CONJUGATE_GRADIENTS_ITERATION;
  } else iteration=MULTIGRID_ITERATION;
  iiteration=static_cast<int>(iteration);
  diffusion=static_cast<DIFFUSION>(idiffusion);
  equation=static_cast<EQUATION>(iequation);
  preconditioner=static_cast<PRECONDITIONER>(ipreconditioner);
  smoother=static_cast<Level::SMOOTHER>(ismoother);
  restriction_prolongation=
   static_cast<Level::RESTRICTION_PROLONGATION>(irestriction_prolongation);
  if (iteration==RICHARDSON_ITERATION || 
  smoother==Level::RICHARDSON_SMOOTHER){
    mu=max(mu,4.);
  }
  if (iteration==JACOBI_OMEGA_ITERATION) {
    omega=min(omega,1.);
  }
#ifdef DEBUG
//cout << "\titeration = " << iiteration << " "
//     << iteration_name[iiteration] << endl;;
//cout << "\tdiffusion = " << idiffusion << " "
//     << diffusion_name[idiffusion] << endl; 
//cout << "\tequation = " << iequation << " " << equation_name[iequation]
//     << endl;
//cout << "\tpreconditioner = " << ipreconditioner << " "
//     << preconditioner_name[ipreconditioner] << endl;
//cout << "\tsmoother = " << ismoother << " "
//     << Level::smoother_name[ismoother] << endl;
//cout << "\trestriction_prolongation = "
//     << irestriction_prolongation << " "
//     << Level::restriction_prolongation_name[irestriction_prolongation]
//     << endl;
//cout << "\tmu,omega = " << mu << " " << omega << endl;
#endif
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void contour_corner(AreaGraphTool &g,const int *fc,const int *lc,
const int *ifirst,const int *ilast,double val,const double *w) {
//TRACER_CALL(tr,"contour_corner");
#ifdef DEBUG
//cout << "\tw = " << endl;
//F77NAME(rbugcell)(DIM_ARG(fc),DIM_ARG(lc),
//  DIM_ARG(ifirst),DIM_ARG(ilast),w);
#endif
  double dif[5],xc[5][2],pos,x[4][2];
  int cstride[2],ie[2];
  for (int id=0;id<SPACEDIM;id++) cstride[id]=lc[id]-fc[id]+1;
  for (ie[1]=fc[1];ie[1]<lc[1];ie[1]++) {
    for (ie[0]=fc[0];ie[0]<lc[0];ie[0]++) {
      dif[0]=w[ie[0]  -fc[0]+cstride[0]*(ie[1]  -fc[1])]-val;
      dif[1]=w[ie[0]+1-fc[0]+cstride[0]*(ie[1]  -fc[1])]-val;
      dif[2]=w[ie[0]+1-fc[0]+cstride[0]*(ie[1]+1-fc[1])]-val;
      dif[3]=w[ie[0]  -fc[0]+cstride[0]*(ie[1]+1-fc[1])]-val;
      xc[0][0]=static_cast<double>(ie[0]  -fc[0]);
      xc[1][0]=static_cast<double>(ie[0]+1-fc[0]);
      xc[2][0]=static_cast<double>(ie[0]+1-fc[0]);
      xc[3][0]=static_cast<double>(ie[0]  -fc[0]);
      xc[0][1]=static_cast<double>(ie[1]  -fc[1]);
      xc[1][1]=static_cast<double>(ie[1]  -fc[1]);
      xc[2][1]=static_cast<double>(ie[1]+1-fc[1]);
      xc[3][1]=static_cast<double>(ie[1]+1-fc[1]);

      dif[4]=dif[0];
      xc[4][0]=xc[0][0];
      xc[4][1]=xc[0][1];
#ifdef DEBUG
//    cout << "\n\tdif = " << dif[0] << " " << dif[1] << " " << dif[2]
//         << " " << dif[3] << " " << dif[4] << endl;
//    cout << "\txc[0] = " << xc[0][0] << " " << xc[1][0] << " " 
//         << xc[2][0] << " " << xc[3][0] << " " << xc[4][0] << endl;
//    cout << "\txc[1] = " << xc[0][1] << " " << xc[1][1] << " " 
//         << xc[2][1] << " " << xc[3][1] << " " << xc[4][1] << endl;
#endif
      int count=0;
      int first,last;
      for (int i=0;i<4;i++) {
        if (dif[i]*dif[i+1]<0.) {
          if (count==0) first=i;
          else last=i;
          count++;
          pos=dif[i]/(dif[i+1]-dif[i]);
          for (int id=0;id<2;id++) {
            x[i][id]=xc[i][id]-pos*(xc[i+1][id]-xc[i][id]);
          }
        }
      }
#ifdef DEBUG
//    cout << "\tcount,first,last = " << count << " " << first << " "
//         << last << endl;
#endif
      if (count==4) {
        g.movePen(x[0][0],x[0][1]);
        g.drawLine(x[2][0],x[2][1]);
        g.movePen(x[1][0],x[1][1]);
        g.drawLine(x[3][0],x[3][1]);
      } else if (count==2) {
        g.movePen(x[first][0],x[first][1]);
        g.drawLine(x[last][0],x[last][1]);
      }
    }
  }
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void runOnce(int *nc) {
//TRACER_CALL(tr,"runOnce");
  if (run_richardson) iteration=RICHARDSON_ITERATION;
  else if (run_jacobi) iteration=JACOBI_OMEGA_ITERATION;
  else if (run_gauss_seidel_red_black) {
    iteration=GAUSS_SEIDEL_RED_BLACK_ITERATION;
  }
  else if (run_gauss_seidel_to_fro) {
    iteration=GAUSS_SEIDEL_TO_FRO_ITERATION;
  }
  else if (run_sor) iteration=SOR_ITERATION;
//else if (run_incomplete_factorization) {
//  iteration=INCOMPLETE_FACTORIZATION_ITERATION;
//}
  else if (run_conjugate_gradients) {
    iteration=CONJUGATE_GRADIENTS_ITERATION;
  }
  else if (run_multigrid) {
    iteration=MULTIGRID_ITERATION;
  }
#ifdef DEBUG
//int iiteration=static_cast<int>(iteration);
//cout << "\titeration = " << iiteration << " "
//     << iteration_name[iiteration] << endl;;
#endif
//array bounds for Fortran calls
  int fc[2],lc[2],ifirst[2],ilast[2];
  double rlo[2],rhi[2];
  int ncorners=1;
  int total_ncells=1;
  for (int id=0;id<2;id++) {
    fc[id]=0;
    lc[id]=nc[id];
    ifirst[id]=1;
    ilast[id]=nc[id]-1;
    ncorners*=lc[id]-fc[id]+1;
    total_ncells*=lc[id]-fc[id];
    rlo[id]=static_cast<double>(fc[id]);
    rhi[id]=static_cast<double>(lc[id]);
  }
#ifdef DEBUG
//cout << "\ttotal_ncells = " << total_ncells << endl;
//cout << "\tncorners = " << ncorners << endl;
#endif
  double log10=log(10.);

  double *diffusion_coef=OPERATOR_NEW_BRACKET(double,total_ncells);
  double *ap=OPERATOR_NEW_BRACKET(double,ncorners);
  double *cg_direction=OPERATOR_NEW_BRACKET(double,ncorners);
  double *cg_preconditioner_solve=OPERATOR_NEW_BRACKET(double,ncorners);
  double *cg_residual=OPERATOR_NEW_BRACKET(double,ncorners);
  double *diagonal_block=OPERATOR_NEW_BRACKET(double,(lc[0]-fc[0]+1)*2);
  double *matrix=OPERATOR_NEW_BRACKET(double,ncorners*9);
  double *residual=OPERATOR_NEW_BRACKET(double,ncorners);
  double *rhs=OPERATOR_NEW_BRACKET(double,ncorners);
  double *solution=OPERATOR_NEW_BRACKET(double,ncorners);
  double *solution_increment=OPERATOR_NEW_BRACKET(double,ncorners);
  double *true_solution=OPERATOR_NEW_BRACKET(double,ncorners);
  double *log_error=OPERATOR_NEW_BRACKET(double,nsteps);
  Level *finest_level=0;

//int n=(ilast[0]-ifirst[0]+1)*(ilast[1]-ifirst[1]+1);
//int nnz=4*(ilast[0]-ifirst[0]+1)*(ilast[1]-ifirst[1]+1)
//       -3*(ilast[0]-ifirst[0]+1)-3*(ilast[1]-ifirst[1]+1)+2;
//double *a=OPERATOR_NEW_BRACKET(double,nnz+n*nbands);
//double *diag=OPERATOR_NEW_BRACKET(double,n);
//int *col_ptr=OPERATOR_NEW_BRACKET(int,n+1);
//int *row_ind=OPERATOR_NEW_BRACKET(int,nnz+n*nbands);

  double gamma=numeric_limits<double>::infinity();
  double residual_max=numeric_limits<double>::infinity();
  double rhs_max=numeric_limits<double>::infinity();

  {
//  TRACER_CALL(tr,"runOnce initialization");
//  initialize arrays to help IEEE exception handling catch unassigned
    for (int i=0;i<ncorners;i++) {
      ap[i]=0.;
      cg_direction[i]=0.;
      cg_preconditioner_solve[i]=0.;
      cg_residual[i]=0.;
      residual[i]=0.;
#ifdef INDEF
      rhs[i]=numeric_limits<double>::infinity();
      solution[i]=numeric_limits<double>::infinity();
      true_solution[i]=numeric_limits<double>::infinity();
      solution_increment[i]=numeric_limits<double>::infinity();
#endif
    }
#ifdef INDEF
    for (int i=0;i<9*ncorners;i++) {
      matrix[i]=numeric_limits<double>::infinity();
    }
#endif
    double c=1./static_cast<double>(RAND_MAX);
    switch (diffusion) {
      case RANDOM_DIFFUSION: {
//      TRACER_CALL(tr,"runOnce RANDOM_DIFFUSION");
        for (int i=0;i<total_ncells;i++) {
          diffusion_coef[i]=static_cast<double>(rand())*c;
        }
        break;
      }
      case CONSTANT_DIFFUSION:
      default: {
//      TRACER_CALL(tr,"runOnce CONSTANT_DIFFUSION");
        for (int i=0;i<total_ncells;i++) diffusion_coef[i]=1.;
        break;
      }
    }
    double dx0=M_PI/static_cast<double>(nc[0]);
    double dx1=M_PI/static_cast<double>(nc[1]);
    for (int i=0;i<ncorners;i++) {
      int i0=i%(nc[0]+1);
      int i1=i/(nc[0]+1);
//    true_solution[i]=
//       static_cast<double>(i0)/static_cast<double>(nc[0])
//      +static_cast<double>(i1)/static_cast<double>(nc[1]);
      true_solution[i]=sin(static_cast<double>(i0)*dx0)
                      *sin(static_cast<double>(i1)*dx1);
      solution[i]=static_cast<double>(rand())*c;
    }
    switch (equation) {
      case HEAT_EQUATION: {
//      TRACER_CALL(tr,"runOnce HEAT_EQUATION");
        F77_NAME(setup_heat)(DIM_ARG(fc),DIM_ARG(lc),
          DIM_ARG(ifirst),DIM_ARG(ilast),decay,
          diffusion_coef,true_solution,matrix,solution);
        break;
      }
      case LAPLACE_EQUATION:
      default: {
//      TRACER_CALL(tr,"runOnce LAPLACE_EQUATION");
        F77_NAME(setup_laplace)(DIM_ARG(fc),DIM_ARG(lc),
          DIM_ARG(ifirst),DIM_ARG(ilast),diffusion_coef,true_solution,
          matrix,solution);
        break;
      }
    }
//  since residual=0, the following sets rhs=matrix * true_solution
    F77_NAME(compute_residual)(DIM_ARG(fc),DIM_ARG(lc),
      DIM_ARG(ifirst),DIM_ARG(ilast),rhs_max,
      matrix,residual,true_solution,rhs);
//  compute matrix * solution - rhs, and store in cg_residual
    F77_NAME(compute_residual)(DIM_ARG(fc),DIM_ARG(lc),
      DIM_ARG(ifirst),DIM_ARG(ilast),residual_max,
      matrix,rhs,solution, cg_residual);
#ifdef DEBUG
//  cout << "\ttrue_solution = " << endl;
//  F77NAME(rbugcell)(DIM_ARG(fc),DIM_ARG(lc),
//    DIM_ARG(ifirst),DIM_ARG(ilast),true_solution);
//  cout << "\trhs = " << endl;
//  F77NAME(rbugcell)(DIM_ARG(fc),DIM_ARG(lc),
//    DIM_ARG(ifirst),DIM_ARG(ilast),rhs);
//  cout << "\tsolution = " << endl;
//  F77NAME(rbugcell)(DIM_ARG(fc),DIM_ARG(lc),
//    DIM_ARG(ifirst),DIM_ARG(ilast),solution);
//  cout << "\tcg_residual = " << endl;
//  F77NAME(rbugcell)(DIM_ARG(fc),DIM_ARG(lc),
//    DIM_ARG(ifirst),DIM_ARG(ilast),cg_residual);
//  for (int i=-2;i<=2;i++) {
//    cout << "\tmatrix[" << i << "] = " << endl;
//    F77NAME(rbugcell)(DIM_ARG(fc),DIM_ARG(lc),
//      DIM_ARG(ifirst),DIM_ARG(ilast),matrix+(i+2)*ncorners);
//  }
#endif
    switch (iteration) {
/*
      case INCOMPLETE_FACTORIZATION_ITERATION: {
//      TRACER_CALL(tr,"runOnce INCOMPLETE_FACTORIZATION_ITERATION");
        int nnz_dicf;
        F77_NAME(iterative_improvement_dicf_setup)(
          DIM_ARG(fc),DIM_ARG(lc),DIM_ARG(ifirst),DIM_ARG(ilast),matrix,
          n,nnz_dicf,a,diag,col_ptr,row_ind);
#ifdef DEBUG
//      cout << "\tnnz = " << nnz << " " << nnz_dicf << endl;
//      for (int i=0;i<n;i++) {
//        cout << "\tdiag[" << i+1 << "] = " << diag[i] << endl;
//      }
//      for (int j=0;j<n;j++) {
//        for (int i=col_ptr[j];i<col_ptr[j+1];i++) {
//          cout << "\ta[" << row_ind[i-1] << "," << j+1 << "] = "
//               << a[i-1] << endl;
//        }
//      }
#endif
        ASSERT(nnz==nnz_dicf);
//      int info;
//      int *indr=OPERATOR_NEW_BRACKET(int,n);
//      int *indf=OPERATOR_NEW_BRACKET(int,n);
//      int *list=OPERATOR_NEW_BRACKET(int,n);
//      double *w=OPERATOR_NEW_BRACKET(double,n);
//      F77NAME(dicf)(n,nnz,a,diag,col_ptr,row_ind,nbands,info,
//        indr,indf,list,w);
//      delete [] indr;
//      delete [] indf;
//      delete [] list;
//      delete [] w;

        double *ta=OPERATOR_NEW_BRACKET(double,n);
        int *ifirst_istdic=OPERATOR_NEW_BRACKET(int,n);
        int *list=OPERATOR_NEW_BRACKET(int,n);
        int info=F77NAME(istdic)(n,diag,a,col_ptr,row_ind,ta,
          ifirst_istdic,list);
        delete [] ta;
        delete [] ifirst_istdic;
        delete [] list;
        break;
      } 
*/
      case CONJUGATE_GRADIENTS_ITERATION: {
//      TRACER_CALL(tr,"runOnce initialization conjugate gradients");
        for (int i=0;i<ncorners;i++) {
          cg_residual[i]=-cg_residual[i];
        }
        switch (preconditioner) {
          case IDENTITY_PRECONDITIONER: {
            for (int i=0;i<ncorners;i++) {
              cg_direction[i]=cg_residual[i];
            }
            break;
          }
          case JACOBI_PRECONDITIONER:
            F77_NAME(jacobi_preconditioner)(DIM_ARG(fc),DIM_ARG(lc),
              DIM_ARG(ifirst),DIM_ARG(ilast),
              matrix,cg_residual, cg_direction);
            break;
          case BLOCK_JACOBI_PRECONDITIONER:
            F77_NAME(block_jacobi_preconditioner)(
              DIM_ARG(fc),DIM_ARG(lc),DIM_ARG(ifirst),DIM_ARG(ilast),
              matrix,cg_residual, diagonal_block,cg_direction);
            break;
/*
          case INCOMPLETE_CHOLESKY_PRECONDITIONER: {
            F77_NAME(iterative_improvement_dicf_setup)(
              DIM_ARG(fc),DIM_ARG(lc),DIM_ARG(ifirst),DIM_ARG(ilast),
              matrix, n,nnz,a,diag,col_ptr,row_ind);
//          int info;
//          int *indr=OPERATOR_NEW_BRACKET(int,n);
//          int *indf=OPERATOR_NEW_BRACKET(int,n);
//          int *list=OPERATOR_NEW_BRACKET(int,n);
//          double *w=OPERATOR_NEW_BRACKET(double,n);
//          F77NAME(dicf)(n,nnz,a,diag,col_ptr,row_ind,nbands,info,
//            indr,indf,list,w);
//          delete [] indr;
//          delete [] indf;
//          delete [] list;
//          delete [] w;

            double *ta=OPERATOR_NEW_BRACKET(double,n);
            int *ifirst_istdic=OPERATOR_NEW_BRACKET(int,n);
            int *list=OPERATOR_NEW_BRACKET(int,n);
            int info=F77NAME(istdic)(n,diag,a,col_ptr,row_ind,ta,
              ifirst_istdic,list);
            F77_NAME(dicf_preconditioner)(
              DIM_ARG(fc),DIM_ARG(lc),DIM_ARG(ifirst),DIM_ARG(ilast),
              n,nnz, a,diag,col_ptr,cg_residual,row_ind, 
              cg_residual,residual);
            delete [] ta;
            delete [] ifirst_istdic;
            delete [] list;
            break;
          }
*/
          case MULTIGRID_PRECONDITIONER: {
          default:
            finest_level=OPERATOR_NEW Level(nc,matrix,smoother_iterations,
              smoother,restriction_prolongation);
            break;
          }
        }
        gamma=0.;
        for (int i=0;i<ncorners;i++) {
          gamma+=cg_residual[i]*cg_direction[i];
        }
#ifdef DEBUG
//      cout << "\tmatrix:" << endl;
//      for (int i=0;i<ncorners;i++) {
//        cout << "\tmatrix[" << i << "] = " << matrix[i] << " "
//             << matrix[i+ncorners] << " "
//             << matrix[i+2*ncorners] << " "
//             << matrix[i+3*ncorners] << " "
//             << matrix[i+4*ncorners] << endl;
//      }
//      for (int i=0;i<ncorners;i++) {
//        cout << "\trhs[" << i << "] = " << rhs[i] << endl;
//      }
//      for (int i=0;i<ncorners;i++) {
//        cout << "\tsolution[" << i << "] = " << solution[i] << endl;
//      }
//      for (int i=0;i<ncorners;i++) {
//        cout << "\tcg_residual[" << i << "] = " << cg_residual[i] 
//             << endl;
//      }
//      for (int i=0;i<ncorners;i++) {
//        cout << "\tcg_direction[" << i << "] = " << cg_direction[i] 
//             << endl;
//      }
//      cout << "\tgamma = " << gamma << endl;
#endif
        break;
      } 
      case MULTIGRID_ITERATION: {
//      TRACER_CALL(tr,"runOnce initialization multigrid");
        finest_level=OPERATOR_NEW Level(nc,matrix,smoother_iterations,
          smoother,restriction_prolongation);
#ifdef DEBUG
//      finest_level->checkMatrixSymmetryRandom();
//      finest_level->coarserLevel()->checkMatrixSymmetryRandom();

//      finest_level->checkProlongationAndRestrictionRandom();
//      finest_level->checkProlongationAndRestriction();
//      finest_level->coarserLevel()->checkProlongationAndRestriction();

//      finest_level->checkSmoother();
//      finest_level->checkSmootherUpdate();
//      finest_level->checkSmootherUpdateRestrict();

//      finest_level->checkCoarseGridProjection();
//      finest_level->checkPointSource();
//      finest_level->checkPointSourceRandom();

//      finest_level->checkVCycleSymmetryRandom();
//      finest_level->coarserLevel()->checkVCycleSymmetryRandom();
//      finest_level->coarserLevel()->coarserLevel()->checkVCycleSymmetry();
//      finest_level->checkVCycleSymmetry();
#endif
        break;
      }
    }
  }

//find min,max data values
  double ulo,uhi;
  F77NAME(rcmpcell)(DIM_ARG(fc),DIM_ARG(lc),
    DIM_ARG(ifirst),DIM_ARG(ilast),solution,uhi,ulo);

  double elo=numeric_limits<double>::infinity();
  double ehi=-numeric_limits<double>::infinity();

//setup interactive graphics
  Palette pal;
  AreaGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
  char label[LENGTH_NAME];
  snprintf(label,LENGTH_NAME,"%s : solution",iteration_name[iteration]);
  AreaGraphTool gt(label,"x0","x1",rlo,rhi,&cmap,0,winsize);

  {
//  TRACER_CALL(tr0,"runOnce initialize graphics");
//  initialize graphics display
    gt.setbgColor("white");
    gt.setfgColor("black");
    gt.drawAxes();
    gt.movePen(rlo[0],rlo[1]);
    gt.drawLine(rhi[0],rlo[1]);
    gt.drawLine(rhi[0],rhi[1]);
    gt.drawLine(rlo[0],rhi[1]);
    gt.drawLine(rlo[0],rlo[1]);

//  draw numerical solution
    double valdif=(uhi-ulo)/static_cast<double>(ncontours+1);
#ifdef DEBUG
//  cout << "\tuhi,ulo,valdif = " << uhi << " " << ulo << " " << valdif
//       << endl;
#endif
    if (valdif>0.) {
      double val=ulo+0.5*valdif;
      for (int il=0;il<ncontours;il++,val+=valdif) {
        gt.setfgColor(il,ncontours);
        contour_corner(gt,fc,lc,ifirst,ilast,val,solution);
      }
    }
//  force X to perform the requests
    gt.flush();
  }

  int stop_step=1;
  int step=0;
  for (;step<nsteps;step++) {
//  TRACER_CALL(tr,"runOnce step");
    switch (iteration) {
      case RICHARDSON_ITERATION: {
//      TRACER_CALL(tr,"runOnce step richardson");
        F77NAME(richardson)(DIM_ARG(fc),DIM_ARG(lc),
          DIM_ARG(ifirst),DIM_ARG(ilast), mu,residual_max, 
          matrix,rhs, solution, residual);
        break;
      }
      case JACOBI_OMEGA_ITERATION: {
//      TRACER_CALL(tr,"runOnce step jacobi");
        F77_NAME(jacobi_omega)(DIM_ARG(fc),DIM_ARG(lc),
          DIM_ARG(ifirst),DIM_ARG(ilast), omega,residual_max,
          matrix,rhs, solution, residual);
        break;
      }
      case GAUSS_SEIDEL_RED_BLACK_ITERATION: {
//      TRACER_CALL(tr,"runOnce step gauss seidel red black");
        F77_NAME(gauss_seidel_red_black)(DIM_ARG(fc),DIM_ARG(lc),
          DIM_ARG(ifirst),DIM_ARG(ilast), step,residual_max,
          matrix,rhs, solution);
        break;
      }
      case GAUSS_SEIDEL_TO_FRO_ITERATION: {
//      TRACER_CALL(tr,"runOnce step gauss seidel to fro");
        F77_NAME(gauss_seidel_to_fro)(DIM_ARG(fc),DIM_ARG(lc),
          DIM_ARG(ifirst),DIM_ARG(ilast), step,residual_max,
        matrix,rhs, solution);
        break;
      }
      case SOR_ITERATION: {
//      TRACER_CALL(tr,"runOnce step sor");
        F77NAME(sor)(DIM_ARG(fc),DIM_ARG(lc),
          DIM_ARG(ifirst),DIM_ARG(ilast), omega,residual_max,
          matrix,rhs, solution);
        break;
      }
/*
      case INCOMPLETE_FACTORIZATION_ITERATION: {
//      TRACER_CALL(tr,"runOnce step iterative_improvement_dicf");
        F77_NAME(iterative_improvement_dicf)(DIM_ARG(fc),DIM_ARG(lc),
          DIM_ARG(ifirst),DIM_ARG(ilast),n,nnz,residual_max,
          a,diag,col_ptr, matrix,rhs,row_ind, solution,residual);
        break;
      }
*/
      case CONJUGATE_GRADIENTS_ITERATION: {
//      TRACER_CALL(tr,"runOnce step conjugate gradients");
        F77_NAME(matrix_multiply)(DIM_ARG(fc),DIM_ARG(lc),
          DIM_ARG(ifirst),DIM_ARG(ilast), matrix,cg_direction, ap);
        double dot=F77NAME(ddot)(ncorners,cg_direction,1,ap,1);
        double alpha=gamma/dot;
        residual_max=0.;
        for (int i=0;i<ncorners;i++) {
          solution[i]+=cg_direction[i]*alpha; // daxpy
          cg_residual[i]-=ap[i]*alpha; // daxpy
          residual_max=max(residual_max,abs(cg_residual[i]));
        }
        switch (preconditioner) {
          case IDENTITY_PRECONDITIONER: {
            memcpy(cg_preconditioner_solve,cg_residual,
              ncorners*sizeof(double));
            break;
          }
          case JACOBI_PRECONDITIONER:
            F77_NAME(jacobi_preconditioner)(DIM_ARG(fc),DIM_ARG(lc),
              DIM_ARG(ifirst),DIM_ARG(ilast),
              matrix,cg_residual, cg_preconditioner_solve);
            break;
          case BLOCK_JACOBI_PRECONDITIONER:
            F77_NAME(block_jacobi_preconditioner)(
              DIM_ARG(fc),DIM_ARG(lc),DIM_ARG(ifirst),DIM_ARG(ilast),
              matrix,cg_residual, 
              diagonal_block,cg_preconditioner_solve);
            break;
/*
          case INCOMPLETE_CHOLESKY_PRECONDITIONER:
            F77_NAME(dicf_preconditioner)(
              DIM_ARG(fc),DIM_ARG(lc),DIM_ARG(ifirst),DIM_ARG(ilast),
              n,nnz, a,diag,col_ptr,cg_residual,row_ind, 
              cg_preconditioner_solve,residual);
            break;
*/
          case MULTIGRID_PRECONDITIONER:
          default:
            finest_level->
              multigridStep(cg_residual,cg_preconditioner_solve);
            break;
        }
        double delta=F77NAME(ddot)(ncorners,cg_preconditioner_solve,1,
          cg_residual,1);
        double beta=delta/gamma;
        for (int i=0;i<ncorners;i++) {
          cg_direction[i]=cg_preconditioner_solve[i]
                         +cg_direction[i]*beta;
        }
        gamma=delta;
        break;
      }
      case MULTIGRID_ITERATION: {
//      TRACER_CALL(tr,"runOnce step multigrid");
        finest_level->multigridStep(cg_residual,solution_increment);
        for (int i=0;i<ncorners;i++) solution[i]-=solution_increment[i];
        residual_max=
          finest_level->computeResidual(solution,rhs,cg_residual);
#ifdef DEBUG
//      cout << "\tsolution_increment = " << endl;
//      F77NAME(rbugcell)(DIM_ARG(fc),DIM_ARG(lc),
//        DIM_ARG(ifirst),DIM_ARG(ilast),solution_increment);
//      cout << "\tsolution = " << endl;
//      F77NAME(rbugcell)(DIM_ARG(fc),DIM_ARG(lc),
//        DIM_ARG(ifirst),DIM_ARG(ilast),solution);
//      cout << "\tcg_residual = " << endl;
//      F77NAME(rbugcell)(DIM_ARG(fc),DIM_ARG(lc),
//        DIM_ARG(ifirst),DIM_ARG(ilast),cg_residual);
#endif
        break;
      }
      default:
        break;
    }
    log_error[step]=
      (residual_max>0. ? log(residual_max/rhs_max)/log10 : -16.);
#ifdef DEBUG
//  cout << "\tresidual_max,log_error[" << step << "] = " 
//       << residual_max << " " << log_error[step] << endl;
#endif
    elo=min(elo,log_error[step]);
    ehi=max(ehi,log_error[step]);

//  initialize graphics display
    gt.newPage();
    gt.setbgColor("white");
    gt.setfgColor("black");
    gt.drawAxes();
    gt.movePen(rlo[0],rlo[1]);
    gt.drawLine(rhi[0],rlo[1]);
    gt.drawLine(rhi[0],rhi[1]);
    gt.drawLine(rlo[0],rhi[1]);
    gt.drawLine(rlo[0],rlo[1]);

//  draw numerical solution
    F77NAME(rcmpcell)(DIM_ARG(fc),DIM_ARG(lc),
      DIM_ARG(ifirst),DIM_ARG(ilast),solution,uhi,ulo);
    double valdif=(uhi-ulo)/static_cast<double>(ncontours+1);
    if (valdif>0.) {
      double val=ulo+0.5*valdif;
      for (int il=0;il<ncontours;il++,val+=valdif) {
        gt.setfgColor(il,ncontours);
        contour_corner(gt,fc,lc,ifirst,ilast,val,solution);
      }
    }
    gt.flush();
        
#ifdef DEBUG
    if (step+1==stop_step) {
      cout << "\tstep number = " << stop_step << endl;
      AreaGraphTool::WINDOW_TYPE::QuitButton qb;
        stop_step *= 10;
    }
#endif
    if (residual_max<=tol*rhs_max) {
      step++;
      break;
    }
  } // end of time step loop
  step--;
#ifdef DEBUG
  cout << "\titeration terminated at " << step << " steps" << endl;
#endif
  snprintf(label,LENGTH_NAME,"%s : error",iteration_name[iteration]);
  XYGraphTool gte(label,"iteration number",
    "log10( error )",0.,step,elo,ehi,&cmap,0,winsize);
  gte.setbgColor("white");
  gte.setfgColor("black");
  gte.drawAxes();
  int s=0;
  gte.setfgColor("blue");
  gte.movePen(static_cast<double>(s),log_error[s]);
  for (s=1;s<step;s++) {
    gte.drawLine(static_cast<double>(s),log_error[s]);
  }
  gte.flush();
  XYGraphTool::WINDOW_TYPE::QuitButton qb;

//since x, u, flux were created with operator new, must delete
  delete [] ap;
  delete [] cg_direction;
  delete [] cg_preconditioner_solve;
  delete [] cg_residual;
  delete [] diagonal_block;
  delete [] diffusion_coef;
  delete [] log_error;
  delete [] solution;
  delete [] solution_increment;
  delete [] true_solution;
  delete [] rhs;
  delete [] residual;
  delete [] matrix;
//delete [] a;
//delete [] diag;
//delete [] col_ptr;
//delete [] row_ind;
  if (finest_level!=0) delete finest_level;
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double runSimulation(int *nc,double &log_error,double &log_it,
double &log_time) {
//TRACER_CALL(tr,"runSimulation");
#ifdef DEBUG
//cout << "\titeration_name = " << iteration_name[iteration] << endl;
#endif
  if (iteration==RICHARDSON_ITERATION) {
    mu=2.*static_cast<double>(nc[0]*nc[0]+nc[1]*nc[1]);
  }
//array bounds for Fortran calls
  int fc[2],lc[2],ifirst[2],ilast[2];
  int ncorners=1;
  int total_ncells=1;
  for (int id=0;id<2;id++) {
    fc[id]=0;
    lc[id]=nc[id];
    ifirst[id]=1;
    ilast[id]=nc[id]-1;
    ncorners*=lc[id]-fc[id]+1;
    total_ncells*=lc[id]-fc[id];
  }
  double log10=log(10.);

  double *diffusion_coef=OPERATOR_NEW_BRACKET(double,total_ncells);
  double *ap=OPERATOR_NEW_BRACKET(double,ncorners);
  double *cg_direction=OPERATOR_NEW_BRACKET(double,ncorners);
  double *cg_preconditioner_solve=OPERATOR_NEW_BRACKET(double,ncorners);
  double *cg_residual=OPERATOR_NEW_BRACKET(double,ncorners);
  double *diagonal_block=OPERATOR_NEW_BRACKET(double,lc[0]-fc[0]+1);
  double *matrix=OPERATOR_NEW_BRACKET(double,ncorners*9);
  double *residual=OPERATOR_NEW_BRACKET(double,ncorners);
  double *rhs=OPERATOR_NEW_BRACKET(double,ncorners);
  double *solution=OPERATOR_NEW_BRACKET(double,ncorners);
  double *solution_increment=OPERATOR_NEW_BRACKET(double,ncorners);
  double *true_solution=OPERATOR_NEW_BRACKET(double,ncorners);
  Level *finest_level=0;

//int n=(ilast[0]-ifirst[0]+1)*(ilast[1]-ifirst[1]+1);
//int nnz=4*(ilast[0]-ifirst[0]+1)*(ilast[1]-ifirst[1]+1)
//       -3*(ilast[0]-ifirst[0]+1)-3*(ilast[1]-ifirst[1]+1)+2;
//double *a=OPERATOR_NEW_BRACKET(double,nnz+n*nbands);
//double *diag=OPERATOR_NEW_BRACKET(double,n);
//int *col_ptr=OPERATOR_NEW_BRACKET(int,n+1);
//int *row_ind=OPERATOR_NEW_BRACKET(int,nnz+n*nbands);

  double gamma=numeric_limits<double>::infinity();
  double residual_max=numeric_limits<double>::infinity();
  double rhs_max=numeric_limits<double>::infinity();

  {
//  TRACER_CALL(tr,"runSimulation initialization");
#ifdef INDEF
//  initialize arrays to help IEEE exception handling catch unassigned
    for (int i=0;i<ncorners;i++) {
      ap[i]=0.;
      cg_direction[i]=0.;
      cg_preconditioner_solve[i]=0.;
      cg_residual[i]=0.;
      residual[i]=0.;
      rhs[i]=numeric_limits<double>::infinity();
      solution[i]=numeric_limits<double>::infinity();
      true_solution[i]=numeric_limits<double>::infinity();
      solution_increment[i]=numeric_limits<double>::infinity();
    }
    for (int i=0;i<9*ncorners;i++) {
      matrix[i]=numeric_limits<double>::infinity();
    }
#endif
    double c=1./static_cast<double>(RAND_MAX);
    switch (diffusion) {
      case RANDOM_DIFFUSION: {
//      TRACER_CALL(tr,"runOnce RANDOM_DIFFUSION");
        for (int i=0;i<total_ncells;i++) {
          diffusion_coef[i]=static_cast<double>(rand())*c;
        }
        break;
      }
      case CONSTANT_DIFFUSION:
      default: {
//      TRACER_CALL(tr,"runOnce CONSTANT_DIFFUSION");
        for (int i=0;i<total_ncells;i++) diffusion_coef[i]=1.;
        break;
      }
    }
    double dx0=M_PI/static_cast<double>(nc[0]);
    double dx1=M_PI/static_cast<double>(nc[1]);
    for (int i=0;i<ncorners;i++) {
      int i0=i%(nc[0]+1);
      int i1=i/(nc[0]+1);
//    true_solution[i]=
//       static_cast<double>(i0)/static_cast<double>(nc[0])
//      +static_cast<double>(i1)/static_cast<double>(nc[1]);
      true_solution[i]=sin(static_cast<double>(i0)*dx0)
                      *sin(static_cast<double>(i1)*dx1);
      solution[i]=static_cast<double>(rand())*c;
    }
    switch(equation) {
      case HEAT_EQUATION:
        F77_NAME(setup_heat)(DIM_ARG(fc),DIM_ARG(lc),
          DIM_ARG(ifirst),DIM_ARG(ilast),decay,
          diffusion_coef,true_solution, matrix,solution);
        break;
      case LAPLACE_EQUATION:
      default:
        F77_NAME(setup_laplace)(DIM_ARG(fc),DIM_ARG(lc),
          DIM_ARG(ifirst),DIM_ARG(ilast),diffusion_coef,true_solution,
          matrix,solution);
        break;
    }
//  since residual=0, the following sets rhs=matrix * true_solution
    F77_NAME(compute_residual)(DIM_ARG(fc),DIM_ARG(lc),
      DIM_ARG(ifirst),DIM_ARG(ilast),rhs_max,
      matrix,residual,true_solution,rhs);
//  compute matrix * solution - rhs, and store in cg_residual
    F77_NAME(compute_residual)(DIM_ARG(fc),DIM_ARG(lc),
      DIM_ARG(ifirst),DIM_ARG(ilast),residual_max,
      matrix,rhs,solution, cg_residual);
#ifdef DEBUG
//  cout << "\ttrue_solution = " << endl;
//  F77NAME(rbugcell)(DIM_ARG(fc),DIM_ARG(lc),
//    DIM_ARG(ifirst),DIM_ARG(ilast),true_solution);
//  cout << "\trhs = " << endl;
//  F77NAME(rbugcell)(DIM_ARG(fc),DIM_ARG(lc),
//    DIM_ARG(ifirst),DIM_ARG(ilast),rhs);
//  cout << "\tsolution = " << endl;
//  F77NAME(rbugcell)(DIM_ARG(fc),DIM_ARG(lc),
//    DIM_ARG(ifirst),DIM_ARG(ilast),solution);
//  cout << "\tcg_residual = " << endl;
//  F77NAME(rbugcell)(DIM_ARG(fc),DIM_ARG(lc),
//    DIM_ARG(ifirst),DIM_ARG(ilast),cg_residual);
//  for (int i=-2;i<=2;i++) {
//    cout << "\tmatrix[" << i << "] = " << endl;
//    F77NAME(rbugcell)(DIM_ARG(fc),DIM_ARG(lc),
//      DIM_ARG(ifirst),DIM_ARG(ilast),matrix+(i+2)*ncorners);
//  }
#endif
    switch (iteration) {
/*
      case INCOMPLETE_FACTORIZATION_ITERATION: {
//      TRACER_CALL(tr,"runSimulation initialization dicf");
        int nnz_dicf;
        F77_NAME(iterative_improvement_dicf_setup)(
          DIM_ARG(fc),DIM_ARG(lc),DIM_ARG(ifirst),DIM_ARG(ilast),matrix,
          n,nnz_dicf,a,diag,col_ptr,row_ind);
        ASSERT(nnz==nnz_dicf);
//      int info;
//      int *indr=OPERATOR_NEW_BRACKET(int,n);
//      int *indf=OPERATOR_NEW_BRACKET(int,n);
//      int *list=OPERATOR_NEW_BRACKET(int,n);
//      double *w=OPERATOR_NEW_BRACKET(double,n);
//      F77NAME(dicf)(n,nnz,a,diag,col_ptr,row_ind,nbands,info,
//        indr,indf,list,w);
//      delete [] indr;
//      delete [] indf;
//      delete [] list;
//      delete [] w;

        double *ta=OPERATOR_NEW_BRACKET(double,n);
        int *ifirst_istdic=OPERATOR_NEW_BRACKET(int,n);
        int *list=OPERATOR_NEW_BRACKET(int,n);
        int info=F77NAME(istdic)(n,diag,a,col_ptr,row_ind,ta,
          ifirst_istdic,list);
        delete [] ta;
        delete [] ifirst_istdic;
        delete [] list;
        break;
      }
*/
      case CONJUGATE_GRADIENTS_ITERATION: {
//      TRACER_CALL(tr,"runSimulation initialization conjugate gradients");
        F77_NAME(compute_residual)(DIM_ARG(fc),DIM_ARG(lc),
          DIM_ARG(ifirst),DIM_ARG(ilast),residual_max,
          matrix,rhs,solution, cg_residual);
        for (int i=0;i<ncorners;i++) {
          cg_residual[i]=-cg_residual[i];
        }
        switch (preconditioner) {
          case IDENTITY_PRECONDITIONER: {
            for (int i=0;i<ncorners;i++) {
              cg_direction[i]=cg_residual[i];
            }
            break;
          }
          case JACOBI_PRECONDITIONER:
            F77_NAME(jacobi_preconditioner)(DIM_ARG(fc),DIM_ARG(lc),
              DIM_ARG(ifirst),DIM_ARG(ilast),
              matrix,cg_residual, cg_direction);
            break;
          case BLOCK_JACOBI_PRECONDITIONER:
            F77_NAME(block_jacobi_preconditioner)(
              DIM_ARG(fc),DIM_ARG(lc),DIM_ARG(ifirst),DIM_ARG(ilast),
              matrix,cg_residual, diagonal_block,cg_direction);
            break;
/*
          case INCOMPLETE_CHOLESKY_PRECONDITIONER: {
            F77_NAME(iterative_improvement_dicf_setup)(
              DIM_ARG(fc),DIM_ARG(lc),DIM_ARG(ifirst),DIM_ARG(ilast),
              matrix, n,nnz,a,diag,col_ptr,row_ind);
//          int info;
//          int *indr=OPERATOR_NEW_BRACKET(int,n);
//          int *indf=OPERATOR_NEW_BRACKET(int,n);
//          int *list=OPERATOR_NEW_BRACKET(int,n);
//          double *w=OPERATOR_NEW_BRACKET(double,n);
//          F77NAME(dicf)(n,nnz,a,diag,col_ptr,row_ind,nbands,info,
//            indr,indf,list,w);
//          delete [] indr;
//          delete [] indf;
//          delete [] list;
//          delete [] w;

            double *ta=OPERATOR_NEW_BRACKET(double,n);
            int *ifirst_istdic=OPERATOR_NEW_BRACKET(int,n);
            int *list=OPERATOR_NEW_BRACKET(int,n);
            int info=F77NAME(istdic)(n,diag,a,col_ptr,row_ind,ta,
              ifirst_istdic,list);
            delete [] ta;
            delete [] ifirst_istdic;
            delete [] list;
            break;
          }
*/
          default:
            break;
        }
        gamma=0.;
        for (int i=0;i<ncorners;i++) {
          gamma+=cg_residual[i]*cg_direction[i];
        }
        break;
      }
      case MULTIGRID_ITERATION: {
//      TRACER_CALL(tr,"runOnce initialization multigrid");
        finest_level=OPERATOR_NEW Level(nc,matrix,smoother_iterations,
          smoother,restriction_prolongation);
        break;
      }
    }
  }

#ifdef DEBUG
//cout << "\tresidual_max,log_error[0] = " << residual_max << " "
//     << log_error[0] << endl;
#endif

  TimedObject iterate_timing("iterate");
  int step=0;
  double initial_residual_norm=numeric_limits<double>::infinity();
  {
    Timer timer(&iterate_timing);
    for (;step<nsteps;step++) {
//    TRACER_CALL(tr,"runSimulation loop");
      switch (iteration) {
        case RICHARDSON_ITERATION:
          F77NAME(richardson)(DIM_ARG(fc),DIM_ARG(lc),
            DIM_ARG(ifirst),DIM_ARG(ilast), mu,residual_max, 
            matrix,rhs, solution, residual);
          break;
        case JACOBI_OMEGA_ITERATION:
          F77_NAME(jacobi_omega)(DIM_ARG(fc),DIM_ARG(lc),
            DIM_ARG(ifirst),DIM_ARG(ilast), omega,residual_max,
            matrix,rhs, solution, residual);
          break;
        case GAUSS_SEIDEL_RED_BLACK_ITERATION:
          F77_NAME(gauss_seidel_red_black)(DIM_ARG(fc),DIM_ARG(lc),
            DIM_ARG(ifirst),DIM_ARG(ilast), step,residual_max,
            matrix,rhs, solution);
          break;
        case GAUSS_SEIDEL_TO_FRO_ITERATION:
            F77_NAME(gauss_seidel_to_fro)(DIM_ARG(fc),DIM_ARG(lc),
              DIM_ARG(ifirst),DIM_ARG(ilast), step,residual_max,
            matrix,rhs, solution);
          break;
        case SOR_ITERATION:
          F77NAME(sor)(DIM_ARG(fc),DIM_ARG(lc),
            DIM_ARG(ifirst),DIM_ARG(ilast), omega,residual_max,
            matrix,rhs, solution);
          break;
/*
        case INCOMPLETE_FACTORIZATION_ITERATION:
          F77_NAME(iterative_improvement_dicf)(DIM_ARG(fc),DIM_ARG(lc),
            DIM_ARG(ifirst),DIM_ARG(ilast),n,nnz,residual_max,
            a,diag,col_ptr, matrix,rhs,row_ind, solution,residual);
          break;
*/
        case CONJUGATE_GRADIENTS_ITERATION: {
          F77_NAME(matrix_multiply)(DIM_ARG(fc),DIM_ARG(lc),
            DIM_ARG(ifirst),DIM_ARG(ilast), matrix,cg_direction, ap);
          double dot=F77NAME(ddot)(ncorners,cg_direction,1,ap,1);
          double alpha=gamma/dot;
          residual_max=0.;
          for (int i=0;i<ncorners;i++) {
            solution[i]+=cg_direction[i]*alpha; // daxpy
            cg_residual[i]-=ap[i]*alpha; // daxpy
            residual_max=max(residual_max,abs(cg_residual[i]));
          }
          switch (preconditioner) {
            case IDENTITY_PRECONDITIONER: {
              memcpy(cg_preconditioner_solve,cg_residual,
                ncorners*sizeof(double));
              break;
            }
            case JACOBI_PRECONDITIONER:
              F77_NAME(jacobi_preconditioner)(DIM_ARG(fc),DIM_ARG(lc),
                DIM_ARG(ifirst),DIM_ARG(ilast),
                matrix,cg_residual, cg_preconditioner_solve);
              break;
            case BLOCK_JACOBI_PRECONDITIONER:
              F77_NAME(block_jacobi_preconditioner)(
                DIM_ARG(fc),DIM_ARG(lc),DIM_ARG(ifirst),DIM_ARG(ilast),
                matrix,cg_residual, 
                diagonal_block,cg_preconditioner_solve);
              break;
/*
            case INCOMPLETE_CHOLESKY_PRECONDITIONER:
              F77_NAME(dicf_preconditioner)(
                DIM_ARG(fc),DIM_ARG(lc),DIM_ARG(ifirst),DIM_ARG(ilast),
                n,nnz, a,diag,col_ptr,cg_residual,row_ind, 
                cg_preconditioner_solve,residual);
              break;
*/
            default:
              break;
          }
          double delta=F77NAME(ddot)(ncorners,cg_preconditioner_solve,1,
            cg_residual,1);
          double beta=delta/gamma;
          for (int i=0;i<ncorners;i++) {
            cg_direction[i]=cg_preconditioner_solve[i]
                           +cg_direction[i]*beta;
          }
          gamma=delta;
          break;
        }
        case MULTIGRID_ITERATION: {
//        TRACER_CALL(tr,"runOnce step multigrid");
          finest_level->multigridStep(cg_residual,solution_increment);
          for (int i=0;i<ncorners;i++) solution[i]-=solution_increment[i];
          residual_max=
            finest_level->computeResidual(solution,rhs,cg_residual);
#ifdef DEBUG
//        cout << "\tsolution_increment = " << endl;
//        F77NAME(rbugcell)(DIM_ARG(fc),DIM_ARG(lc),
//          DIM_ARG(ifirst),DIM_ARG(ilast),solution_increment);
//        cout << "\tsolution = " << endl;
//        F77NAME(rbugcell)(DIM_ARG(fc),DIM_ARG(lc),
//          DIM_ARG(ifirst),DIM_ARG(ilast),solution);
//        cout << "\tcg_residual = " << endl;
//        F77NAME(rbugcell)(DIM_ARG(fc),DIM_ARG(lc),
//          DIM_ARG(ifirst),DIM_ARG(ilast),cg_residual);
#endif
          break;
        }
        default:
          break;
      }
      if (step==0) initial_residual_norm=residual_max;
      if (residual_max<=tol*rhs_max) break;
#ifdef DEBUG
//  cout << "\tresidual_max[" << step << "] = " << residual_max << endl;
#endif
    }
  }
  log_error=(residual_max>0. ? log(residual_max/rhs_max)/log10 : -16.);
  log_it=log(static_cast<double>(step))/log10;
  double elapsed=iterate_timing.totalRunTime();
  log_time=(elapsed>0. ? log(elapsed)/log10 : -3. );
  double log_spectral_radius=
    (log_error-log(initial_residual_norm/rhs_max)/log10)/step;

  delete [] ap;
  delete [] cg_direction;
  delete [] cg_preconditioner_solve;
  delete [] cg_residual;
  delete [] diagonal_block;
  delete [] solution;
  delete [] rhs;
  delete [] residual;
  delete [] matrix;
//delete [] a;
//delete [] diag;
//delete [] col_ptr;
//delete [] row_ind;
  return log_spectral_radius;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void runMultipleSimulations(unsigned int number_simulations,int *nc,
unsigned int runs,double &elo,double &ehi,double &tlo,double &thi,
double &itlo,double &ithi,double *log_error,double *log_time,
double *log_it,double *log_size) {
//TRACER_CALL(z,"runMultipleSimulations");
  switch (runs) {
    case 1:
      iteration=RICHARDSON_ITERATION;
      break;
    case 2:
      iteration=JACOBI_OMEGA_ITERATION;
      break;
    case 4:
      iteration=GAUSS_SEIDEL_RED_BLACK_ITERATION;
      break;
    case 8:
      iteration=GAUSS_SEIDEL_TO_FRO_ITERATION;
      break;
    case 16:
      iteration=SOR_ITERATION;
      break;
    case 32:
/*
      iteration=INCOMPLETE_FACTORIZATION_ITERATION;
      break;
    case 64:
*/
      iteration=CONJUGATE_GRADIENTS_ITERATION;
      break;
  }
    
#ifdef DEBUG
  cout << "\n" << iteration_name[iteration] << " : " 
       << color_name[iteration] << endl;
#endif
  double log10=log(10.);
  for (int nsim=0;nsim<number_simulations;nsim++) {
    log_size[nsim]=log(static_cast<double>((nc[0]-1)*(nc[1]-1)))/log10;
    double rho=
      runSimulation(nc,log_error[nsim],log_it[nsim],log_time[nsim]);
#ifdef DEBUG
    cout << "\tlog time,its,error,rho(ncells = " << nc[0] << ") = "
         << log_time[nsim] << " " << log_it[nsim] << " "
         << log_error[nsim] << " " << rho << endl;
//  if (nsim>0) {
//    cout << " slope = "
//         << (log_error[nsim]-log_error[nsim-1])
//           /(log_it[nsim]-log_it[nsim-1]);
//  }
    cout << endl;
#endif
    elo=min(elo,log_error[nsim]);
    ehi=max(ehi,log_error[nsim]);
    tlo=min(tlo,log_time[nsim]);
    thi=max(thi,log_time[nsim]);
    itlo=0.;
    ithi=max(ithi,log_it[nsim]);
    for (int id=0;id<SPACEDIM;id++) nc[id]*=2;
  }
#ifdef DEBUG
//cout << "\telo,ehi = " << elo << " " << ehi << endl;
//cout << "\ttlo,thi = " << tlo << " " << thi << endl;
//cout << "\titlo,ithi = " << itlo << " " << ithi << endl;
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void markPoint(XYGraphTool *g,double x,double y,double dx,int s) {
  switch (s) {
    case 0:
      g->drawBoxGivenCenter(x,y,dx);
      break;
    case 1:
      g->drawDiamondGivenCenter(x,y,dx);
      break;
    case 2:
      g->drawBoxGivenCenter(x,y,dx);
      g->drawPlus(x,y,dx);
      break;
    case 3:
      g->drawDiamondGivenCenter(x,y,dx);
      g->drawCross(x,y,dx/sqrt(2.));
      break;
    case 4:
      g->drawBoxGivenCenter(x,y,dx);
      g->drawCross(x,y,dx);
      break;
    case 5:
      g->drawDiamondGivenCenter(x,y,dx);
      g->drawPlus(x,y,dx);
      break;
    case 6:
      g->drawBoxGivenCenter(x,y,dx);
      g->drawPlus(x,y,dx);
      g->drawCross(x,y,dx);
      break;
    case 7:
      g->drawDiamondGivenCenter(x,y,dx);
      g->drawPlus(x,y,dx);
      g->drawCross(x,y,dx/sqrt(2.));
      break;
    case 8:
      g->drawBoxGivenCenter(x,y,dx);
      g->drawDiamondGivenCenter(x,y,dx);
      break;
    case 9:
      g->drawBoxGivenCenter(x,y,dx);
      g->drawPlus(x,y,dx);
      g->drawDiamondGivenCenter(x,y,dx);
      break;
    case 10:
      g->drawBoxGivenCenter(x,y,dx);
      g->drawCross(x,y,dx);
      g->drawDiamondGivenCenter(x,y,dx);
      break;
    case 11:
      g->drawBoxGivenCenter(x,y,dx);
      g->drawPlus(x,y,dx);
      g->drawCross(x,y,dx);
      g->drawDiamondGivenCenter(x,y,dx);
      break;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void plotErrors(XYGraphTool *gt1,XYGraphTool *gt2,ITERATION s,
int number_simulations,double itlo,double ithi,double tlo,double thi,
const double *log_size,const double *log_time,const double *log_it) {
//TRACER_CALL(tr,"plotErrors");
#ifdef DEBUG
//cout << "\tcolor_name[" << s << "] = " << color_name[s] << endl;
//cout << "\titlo,ithi = " << itlo << " " << ithi << endl;
//cout << "\ttlo,thi = " << tlo << " " << thi << endl;
//for (int i=0;i<number_simulations;i++) {
//  cout << "\tlogs = " << log_it[i] << " " << log_time[i] << " "
//       << log_size[i] << endl;
//}
#endif
  gt1->setfgColor(color_name[s]);
  gt1->movePen(log_size[0],log_it[0]);
  for (int i=1;i<number_simulations;i++) {
    gt1->drawLine(log_size[i],log_it[i]);
  }
  double db=(ithi-itlo)/static_cast<double>(number_simulations*10);
  for (int i=0;i<number_simulations;i++) {
    markPoint(gt1,log_size[i],log_it[i],db,s);
  }
  gt2->setfgColor(color_name[s]);
  gt2->movePen(log_size[0],log_time[0]);
  for (int i=1;i<number_simulations;i++) {
    gt2->drawLine(log_size[i],log_time[i]);
  }
  db=(thi-tlo)/static_cast<double>(number_simulations*10);
  for (int i=0;i<number_simulations;i++) {
    markPoint(gt2,log_size[i],log_time[i],db,s);
  }
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void runMain(bool /*called_before*/) {
//TRACER_CALL(tr,"runMain");
  if (ncells[0]>0 && ncells[1]>0) {
    runOnce(ncells);
    return;
  }

  unsigned int number_simulations=4;
  double *log_error_r=OPERATOR_NEW_BRACKET(double,number_simulations);
  double *log_nits_r=OPERATOR_NEW_BRACKET(double,number_simulations);
  double *log_time_r=OPERATOR_NEW_BRACKET(double,number_simulations);

  double *log_error_j=OPERATOR_NEW_BRACKET(double,number_simulations);
  double *log_nits_j=OPERATOR_NEW_BRACKET(double,number_simulations);
  double *log_time_j=OPERATOR_NEW_BRACKET(double,number_simulations);

  double *log_error_gsrb=
    OPERATOR_NEW_BRACKET(double,number_simulations);
  double *log_nits_gsrb=OPERATOR_NEW_BRACKET(double,number_simulations);
  double *log_time_gsrb=OPERATOR_NEW_BRACKET(double,number_simulations);

  double *log_error_gstf=
    OPERATOR_NEW_BRACKET(double,number_simulations);
  double *log_nits_gstf=OPERATOR_NEW_BRACKET(double,number_simulations);
  double *log_time_gstf=OPERATOR_NEW_BRACKET(double,number_simulations);

  double *log_error_sor=OPERATOR_NEW_BRACKET(double,number_simulations);
  double *log_nits_sor=OPERATOR_NEW_BRACKET(double,number_simulations);
  double *log_time_sor=OPERATOR_NEW_BRACKET(double,number_simulations);

  double *log_error_if=OPERATOR_NEW_BRACKET(double,number_simulations);
  double *log_nits_if=OPERATOR_NEW_BRACKET(double,number_simulations);
  double *log_time_if=OPERATOR_NEW_BRACKET(double,number_simulations);

  double *log_error_cg=OPERATOR_NEW_BRACKET(double,number_simulations);
  double *log_nits_cg=OPERATOR_NEW_BRACKET(double,number_simulations);
  double *log_time_cg=OPERATOR_NEW_BRACKET(double,number_simulations);

  double *log_error_m=OPERATOR_NEW_BRACKET(double,number_simulations);
  double *log_nits_m=OPERATOR_NEW_BRACKET(double,number_simulations);
  double *log_time_m=OPERATOR_NEW_BRACKET(double,number_simulations);

  double *log_size=OPERATOR_NEW_BRACKET(double,number_simulations);

  int nc[2]={10,10};
  double elo=numeric_limits<double>::infinity();
  double ehi=-numeric_limits<double>::infinity();
  double tlo=numeric_limits<double>::infinity();
  double thi=-numeric_limits<double>::infinity();
  double itlo=numeric_limits<double>::infinity();
  double ithi=-numeric_limits<double>::infinity();

  unsigned int runs=1;
  if (run_richardson) {
    runMultipleSimulations(number_simulations,nc,runs,elo,ehi,tlo,thi,
      itlo,ithi,log_error_r,log_time_r,log_nits_r,log_size);
  }
  runs*=2;
  nc[0]=nc[1]=10;
  if (run_jacobi) {
    runMultipleSimulations(number_simulations,nc,runs,elo,ehi,tlo,thi,
      itlo,ithi,log_error_j,log_time_j,log_nits_j,log_size);
  }
  runs*=2;
  nc[0]=nc[1]=10;
  if (run_gauss_seidel_red_black) {
    runMultipleSimulations(number_simulations,nc,runs,elo,ehi,tlo,thi,
      itlo,ithi,log_error_gsrb,log_time_gsrb,log_nits_gsrb,log_size);
  }
  runs*=2;
  nc[0]=nc[1]=10;
  if (run_gauss_seidel_to_fro) {
    runMultipleSimulations(number_simulations,nc,runs,elo,ehi,tlo,thi,
      itlo,ithi,log_error_gstf,log_time_gstf,log_nits_gstf,log_size);
  }
  runs*=2;
  nc[0]=nc[1]=10;
  if (run_sor) {
    runMultipleSimulations(number_simulations,nc,runs,elo,ehi,tlo,thi,
      itlo,ithi,log_error_sor,log_time_sor,log_nits_sor,log_size);
  }
/*
  runs*=2;
  nc[0]=nc[1]=10;
  if (run_incomplete_factorization) {
    runMultipleSimulations(number_simulations,nc,runs,elo,ehi,tlo,thi,
      itlo,ithi,log_error_if,log_time_if,log_nits_if,log_size);
  }
*/
  runs*=2;
  nc[0]=nc[1]=10;
  if (run_conjugate_gradients) {
    runMultipleSimulations(number_simulations,nc,runs,elo,ehi,tlo,thi,
      itlo,ithi,log_error_cg,log_time_cg,log_nits_cg,log_size);
  }
  runs*=2;
  nc[0]=nc[1]=10;
  if (run_multigrid) {
    runMultipleSimulations(number_simulations,nc,runs,elo,ehi,tlo,thi,
      itlo,ithi,log_error_cg,log_time_cg,log_nits_cg,log_size);
  }

//setup interactive graphics
  Palette pal;
  XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
  char label[LENGTH_NAME];
  snprintf(label,LENGTH_NAME,"%s : iterations vs unknowns",
    iteration_name[iteration]);
  XYGraphTool gt1(label,"log(no. unknowns)","log(iterations)",
    log_size[0],log_size[number_simulations-1],itlo,ithi,
    &cmap,0,winsize);
  snprintf(label,LENGTH_NAME,"%s : time vs unknowns",
    iteration_name[iteration]);
  XYGraphTool gt2(label,"log(no. unknowns)","log(time)",
    log_size[0],log_size[number_simulations-1],tlo,thi,
    &cmap,0,winsize);

  gt1.newPage();
  gt1.setbgColor("white");
  gt1.setfgColor("black");
  gt1.drawAxes();
  gt2.newPage();
  gt2.setbgColor("white");
  gt2.setfgColor("black");
  gt2.drawAxes();
  if (run_richardson) {
    plotErrors(&gt1,&gt2,RICHARDSON_ITERATION,number_simulations,
      itlo,ithi,tlo,thi,log_size,log_time_r,log_nits_r);
  }
  if (run_jacobi) {
    plotErrors(&gt1,&gt2,JACOBI_OMEGA_ITERATION,number_simulations,
      itlo,ithi,tlo,thi,log_size,log_time_j,log_nits_j);
  }
  if (run_gauss_seidel_red_black) {
    plotErrors(&gt1,&gt2,GAUSS_SEIDEL_RED_BLACK_ITERATION,
      number_simulations,itlo,ithi,tlo,thi,
      log_size,log_time_gsrb,log_nits_gsrb);
  }
  if (run_gauss_seidel_to_fro) {
    plotErrors(&gt1,&gt2,GAUSS_SEIDEL_TO_FRO_ITERATION,
      number_simulations,itlo,ithi,tlo,thi,
      log_size,log_time_gstf,log_nits_gstf);
  }
  if (run_sor) {
    plotErrors(&gt1,&gt2,SOR_ITERATION,
      number_simulations,itlo,ithi,tlo,thi,
      log_size,log_time_sor,log_nits_sor);
  }
/*
  if (run_incomplete_factorization) {
    plotErrors(&gt1,&gt2,INCOMPLETE_FACTORIZATION_ITERATION,
      number_simulations,itlo,ithi,tlo,thi,
      log_size,log_time_if,log_nits_if);
  }
*/
  if (run_conjugate_gradients) {
    plotErrors(&gt1,&gt2,CONJUGATE_GRADIENTS_ITERATION,
      number_simulations,itlo,ithi,tlo,thi,
      log_size,log_time_cg,log_nits_cg);
  }
  if (run_multigrid) {
    plotErrors(&gt1,&gt2,CONJUGATE_GRADIENTS_ITERATION,
      number_simulations,itlo,ithi,tlo,thi,
      log_size,log_time_m,log_nits_m);
  }
  XYGraphTool::WINDOW_TYPE::QuitButton qb;

  delete [] log_error_r;
  delete [] log_nits_r;
  delete [] log_time_r;

  delete [] log_error_j;
  delete [] log_nits_j;
  delete [] log_time_j;

  delete [] log_error_gsrb;
  delete [] log_nits_gsrb;
  delete [] log_time_gsrb;

  delete [] log_error_gstf;
  delete [] log_nits_gstf;
  delete [] log_time_gstf;

  delete [] log_error_sor;
  delete [] log_nits_sor;
  delete [] log_time_sor;

  delete [] log_error_if;
  delete [] log_nits_if;
  delete [] log_time_if;

  delete [] log_error_cg;
  delete [] log_nits_cg;
  delete [] log_time_cg;

  delete [] log_error_m;
  delete [] log_nits_m;
  delete [] log_time_m;

  delete [] log_size;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void cleanup() {
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void shutdown() {
  while (param_list->notEmpty()) delete param_list->delAfter(0);
  delete param_list; param_list=0;
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int main(int argc,char* argv[]) {
#ifdef DEBUG
//setTraps();
#endif
#ifdef MEM_DEBUG
  MemoryDebugger md(1);
#endif
#ifdef USE_GTK
  GTKWindow::gtkInit(argc,argv);
#endif
  {
//  set machine-dependent constants for Fortran
    F77NAME(machine).roundoff=DBL_EPSILON;
    F77NAME(machine).small=DBL_MIN;
    F77NAME(machine).huge=DBL_MAX;
    F77NAME(machine).undefind=numeric_limits<double>::infinity();
    F77NAME(const).e=M_E;
    F77NAME(const).pi=M_PI;
    F77NAME(const).rt2=M_SQRT2;
    F77NAME(const).half=0.5;
    F77NAME(const).one=1.;
    F77NAME(const).two=2.;
    F77NAME(const).zero=0.;
    F77NAME(const).debug_on=0;

//  define input parameters and bounds
    param_list=OPERATOR_NEW GUI_INPUT_PARAMETER_LIST_TYPE("Main List");
    { const char *group="Numerical Method Parameters";
      param_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
        iequation,"equation",equation_name[0],equation_name[1]));
      param_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
        idiffusion,"diffusion",diffusion_name[0],diffusion_name[1]));
    }
    { const char *group="Numerical Method Parameters";
      param_list->append(OPERATOR_NEW GUIInputParameter<int>(ncells[0],
        "ncells[0]",0,INT_MAX,group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<int>(ncells[1],
        "ncells[1]",0,INT_MAX,group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<int>(nsteps,
        "nsteps",1,INT_MAX,group) );
//    param_list->append(OPERATOR_NEW GUIInputParameter<int>(nbands,
//      "nbands",1,INT_MAX,group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<int>(
        smoother_iterations,"smoother_iterations",1,INT_MAX,group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<double>(decay,
        "heat equation decay number",0.,DBL_MAX,group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<double>(mu,
        "Richardson mu",2.,DBL_MAX,group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<double>(omega,
        "Jacobi/SOR omega",0.,2.,group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<double>(tol,
        "convergence tolerance",0.,1.,group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<bool>(
        run_richardson,"richardson",false,true,group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<bool>(
        run_jacobi,"jacobi",false,true,group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<bool>(
        run_gauss_seidel_red_black,"gauss-seidel red-black",false,true,
        group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<bool>(
        run_gauss_seidel_to_fro,"gauss-seidel to-fro",false,true,
        group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<bool>(
        run_sor,"sor",false,true, group) );
//    param_list->append(OPERATOR_NEW GUIInputParameter<bool>(
//      run_incomplete_factorization,"incomplete factorization",
//      false,true, group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<bool>(
       run_conjugate_gradients,"conjugate gradients",false,true,group));
      param_list->append(OPERATOR_NEW GUIInputParameter<bool>(
       run_multigrid,"multigrid",false,true,group));
//    param_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
//      iiteration,"iteration",iteration_name[0],iteration_name[1],
//      iteration_name[2],iteration_name[3],iteration_name[4],
//      iteration_name[5],iteration_name[6],iteration_name[7]));
      param_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
        ipreconditioner,"preconditioner",preconditioner_name[0],
        preconditioner_name[1],preconditioner_name[2],
        preconditioner_name[3] /* ,preconditioner_name[4] */ ));
      param_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
        ismoother,"smoother",Level::smoother_name[0],
        Level::smoother_name[1],Level::smoother_name[2]));
      param_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
        irestriction_prolongation,"restriction_prolongation",
        Level::restriction_prolongation_name[0],
        Level::restriction_prolongation_name[1]));
    }
    { const char *group="Graphics Parameters";
      param_list->append(OPERATOR_NEW GUIInputParameter<int>(ncontours,
        "ncontours",1,30,group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<double>(winsize,
        "winsize",0.,1.,group) );
    }

//  read input parameters 
    if (argc>1) {
      ifstream in_file;
      in_file.open(argv[1],ios::in);
      if (in_file) {
        istream &infile(in_file);
        infile.clear(ios::goodbit);
        infile.seekg(0,ios::beg);
        char name[LENGTH_NAME], comment[LENGTH_COMMENT];
        bool found_main=FALSE;
        while ( in_file >> setw(LENGTH_NAME) >> name ) {
          if ( strcmp(name,"Main") == 0 ) {
            in_file.getline( comment, LENGTH_COMMENT);
            found_main=TRUE;
            while ( in_file >> setw(LENGTH_NAME) >> name ) {
              if ( strcmp(name,"end") == 0 ) break;
              else if (strcmp(name,"skip_gui") == 0) {
                int iskip_gui; // type bool not read correctly
                in_file >> iskip_gui;
                skip_gui= iskip_gui!=0;
              } else param_list->formattedRead(in_file,name);
              in_file.getline( comment, LENGTH_COMMENT);
            }
          } else in_file.getline(comment,LENGTH_COMMENT);
        }
        if ( !found_main ) skip_gui=FALSE;
      }
      if (skip_gui) checkMainInput();
    } else skip_gui=FALSE;

    if (skip_gui) {
      runMain(FALSE);
      cleanup();
      shutdown();
    } else {
#ifdef USE_GTK
      GTKGUI gui(argv[0],display_name,param_list,&runMain,
      &checkMainInput,&cleanup,&shutdown,TRUE);
#else
      GUI gui(argv[0],display_name,param_list,&runMain,&checkMainInput,
      &cleanup,&shutdown,TRUE);
#endif
      gui.createFileMenu();
      gui.createViewMenu();
      gui.createHelpMenu();
      gui.createMainWindow(argc,argv);
      gui.eventLoop();
    }
  } // scope of MemoryDebugger
  return EXIT_SUCCESS;
}

template class InputParameter<int>;
template class InputParameter<double>;
template class InputParameter<bool>;
