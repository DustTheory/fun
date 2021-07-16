#include <float.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdlib.h>
#include <math.h>
#include "Arch.H"
#include "Debug.H"
#include "Fort1D.H"
#include "GUIEnumInputParameter.H"
#include "GUIVirtualInput.H"
#include "InputParameter.H"
#include "Level1D.H"
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

struct machine_common {
  double roundoff,small,huge,undefind;
};
extern machine_common F77NAME(machine);

clock_t start;
struct tms usage;

GUI_INPUT_PARAMETER_LIST_TYPE *param_list=0;
bool skip_gui=FALSE;
Level *finest_level=0;
int ncells=8;
int nsteps=10;
int smoother_iterations=1;
double decay_number=1.;
double mu_in=4.;
double mu=numeric_limits<double>::infinity();
double omega=1.;
double winsize=0.5;
char* display_name=0;

enum EQUATION{LAPLACE_EQUATION,HEAT_EQUATION /* ,SINH */ };
EQUATION equation=LAPLACE_EQUATION;
int iequation=equation;
static const char *equation_name[2]={
  "laplace_equation","heat_equation" /* ,"sinh" */
};

enum DIFFUSION{CONSTANT_DIFFUSION,RANDOM_DIFFUSION};
DIFFUSION diffusion=RANDOM_DIFFUSION;
int idiffusion=diffusion;
static const char *diffusion_name[2]={
  "constant_diffusion","random_diffusion"
};

enum ITERATION{
  RICHARDSON_ITERATION,//JACOBI,
  JACOBI_OMEGA_ITERATION,
  GAUSS_SEIDEL_ITERATION,
  GAUSS_SEIDEL_RED_BLACK_ITERATION,
  GAUSS_SEIDEL_TO_FRO_ITERATION,
  SOR_ITERATION,
  CONJUGATE_GRADIENTS_ITERATION,
  MULTIGRID_ITERATION};
ITERATION iteration=MULTIGRID_ITERATION;
int iiteration=iteration;
static const char *iteration_name[8]={
  "richardson",//"jacobi",
  "jacobi_omega","gauss_seidel", "gauss_seidel_red_black",
  "gauss_seidel_to_fro","sor","conjugate_gradients","multigrid"
};

enum PRECONDITIONER{IDENTITY_PRECONDITIONER,JACOBI_PRECONDITIONER,
  MULTIGRID_PRECONDITIONER};
PRECONDITIONER preconditioner=MULTIGRID_PRECONDITIONER;
int ipreconditioner=preconditioner;
static const char *preconditioner_name[3]={
  "identity","jacobi","multigrid"
};

//Level::SMOOTHER smoother=Level::GAUSS_SEIDEL_RED_BLACK_SMOOTHER;
Level::SMOOTHER smoother=Level::GAUSS_SEIDEL_SMOOTHER;
int ismoother=smoother;
//Level::RESTRICTION_PROLONGATION restriction_prolongation=
//  Level::ALGEBRAIC_MULTIGRID_RESTRICTION_PROLONGATION;
Level::RESTRICTION_PROLONGATION restriction_prolongation=
  Level::FINITE_ELEMENT_RESTRICTION_PROLONGATION;
int irestriction_prolongation=restriction_prolongation;

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void checkMainInput() {
//TRACER_CALL(tr,"checkMainInput");
  equation=static_cast<EQUATION>(iequation);
  diffusion=static_cast<DIFFUSION>(idiffusion);
  iteration=static_cast<ITERATION>(iiteration);
  if (iteration==JACOBI_OMEGA_ITERATION) {
    omega=min(omega,2./3.); // 2 / max no. nonzeros in any row of matrix
  }
  if (iteration==RICHARDSON_ITERATION) {
    mu=mu_in*static_cast<double>(ncells); // matrix is proportional to 1/dx
//  cout << "\tmu = " << mu << endl;
  }
  preconditioner=static_cast<PRECONDITIONER>(ipreconditioner);
  smoother=static_cast<Level::SMOOTHER>(ismoother);
  restriction_prolongation=static_cast<Level::RESTRICTION_PROLONGATION>(
    irestriction_prolongation);
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void plotSolution(XYGraphTool &gt,const int &ifirst,const int &ilast,
const double *solution,const double *true_solution) {
//initialize graphics display
  double slo=numeric_limits<double>::infinity();
  double shi=-numeric_limits<double>::infinity();
  for (int ic=ifirst-1;ic<=ilast+1;ic++) {
    slo=min(slo,true_solution[ic]);
    shi=max(shi,true_solution[ic]);
  }
  gt.rescale(static_cast<double>(ifirst),static_cast<double>(ilast),
    slo,shi);

  gt.setbgColor("white");
  gt.setfgColor("black");
  gt.drawAxes();
//draw numerical solution
  gt.setfgColor("blue");
  int ic=ifirst;
  gt.movePen(static_cast<double>(ic),solution[ic]);
#ifdef DEBUG
  if (abs(solution[ic])>10.) abort();
#endif
  for (ic=ifirst+1;ic<=ilast;ic++) {
    gt.drawLine(static_cast<double>(ic),solution[ic]);
#ifdef DEBUG
    if (abs(solution[ic])>10.) abort();
#endif
  }
  gt.setfgColor("red");
  ic=ifirst;
  gt.movePen(static_cast<double>(ic),true_solution[ic]);
  for (ic=ifirst+1;ic<=ilast;ic++) {
    gt.drawLine(static_cast<double>(ic),true_solution[ic]);
  }
//force X to perform the requests
  gt.flush();
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void preconditionedConjugateGradients(const int &fi,const int &la,
const int &ifirst,const int &ilast,double &gamma,
const double *matrix,
double *cg_direction,double *cg_residual,double *solution,
double *ap,double *cg_preconditioner_solve) {
//TRACER_CALL(tr,"preconditionedConjugateGradients");
  F77_NAME(matrix_multiply)(fi,la,ifirst,ilast, matrix,cg_direction, ap);
  int nc=ilast-ifirst+1;
  double dot=F77NAME(ddot)(nc,&cg_direction[ifirst],1,&ap[ifirst],1);
  double alpha=gamma/dot;
  F77NAME(daxpy)(nc,alpha,&cg_direction[ifirst],1,&solution[ifirst],1);

  F77NAME(daxpy)(nc,-alpha,&ap[ifirst],1,&cg_residual[ifirst],1);
//F77_NAME(compute_residual)(fi,la,ifirst,ilast,matrix,rhs,solution,
//  cg_residual);
//for (int i=ifirst;i<=ilast;i++) cg_residual[i]=-cg_residual[i];

//int i=F77NAME(idamax)(nc,&cg_residual[ifirst],1)-1;
//double residual_max=abs(cg_residual[i]);
  switch (preconditioner) {
    case IDENTITY_PRECONDITIONER: {
      memcpy(&cg_preconditioner_solve[ifirst],&cg_residual[ifirst],
        nc*sizeof(double));
      break;
    }
    case JACOBI_PRECONDITIONER:
      F77_NAME(jacobi_preconditioner)(fi,la,ifirst,ilast,
        matrix,cg_residual, cg_preconditioner_solve);
      break;
    case MULTIGRID_PRECONDITIONER:
    default:
      finest_level->multigridStep(cg_residual,cg_preconditioner_solve);
      break;
  }
  double delta=F77NAME(ddot)(nc,&cg_preconditioner_solve[ifirst],1,
    &cg_residual[ifirst],1);
  double beta=delta/gamma;
  for (int i=ifirst;i<=ilast;i++) {
    cg_direction[i]=cg_preconditioner_solve[i]
                   +cg_direction[i]*beta;
  }
  gamma=delta;
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void computeError(const int &ifirst,const int &ilast, 
const double *solution,const double *true_solution, double &log_err) {
  log_err=0.;
  for (int i=ifirst;i<=ilast;i++) {
    log_err=max(log_err,abs(solution[i]-true_solution[i]));
  }
  if (log_err<=0.) log_err=-16.*log(10.);
  else log_err=log(log_err);
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void runMain(bool /*called_before*/) {
//TRACER_CALL(tr,"runMain");
#ifdef DEBUG
//cout << "\tncells = " << ncells << endl;
#endif
//array bounds for Fortran calls
  double log10=log(10.);
  double dx=1./static_cast<double>(ncells);
  double dt=decay_number*dx*dx;

  int n=ncells+1; // interior unknowns plus boundary values
  int n2=2*n;
  int n3=3*n;
  double *matrix_c=OPERATOR_NEW_BRACKET(double,n3);
  double *residual_c=OPERATOR_NEW_BRACKET(double,n);
  double *rhs_c=OPERATOR_NEW_BRACKET(double,n);
  double *solution_c=OPERATOR_NEW_BRACKET(double,n);
  double *ap_c=OPERATOR_NEW_BRACKET(double,n);
  double *cg_direction_c=OPERATOR_NEW_BRACKET(double,n);
  double *cg_preconditioner_solve_c=OPERATOR_NEW_BRACKET(double,n);
  double *cg_residual_c=OPERATOR_NEW_BRACKET(double,n);
//double *mg_residual_c=OPERATOR_NEW_BRACKET(double,n);
  double *true_solution_c=OPERATOR_NEW_BRACKET(double,n);
  double *solution_increment_c=OPERATOR_NEW_BRACKET(double,n);
  double *matrix_copy_c=OPERATOR_NEW_BRACKET(double,n3);
  double *log_error=OPERATOR_NEW_BRACKET(double,nsteps+1);
#ifdef INDEF
//initialize arrays to help IEEE exception handling catch unassigned
  for (int i=0;i<n;i++) {
    ap_c[i]=numeric_limits<double>::infinity();
    cg_direction_c[i]=numeric_limits<double>::infinity();
    cg_preconditioner_solve_c[i]=numeric_limits<double>::infinity();
    cg_residual_c[i]=numeric_limits<double>::infinity();
//  mg_residual_c[i]=numeric_limits<double>::infinity();
    residual_c[i]=numeric_limits<double>::infinity();
    rhs_c[i]=numeric_limits<double>::infinity();
    solution_c[i]=numeric_limits<double>::infinity();
    solution_increment_c[i]=numeric_limits<double>::infinity();
    true_solution_c[i]=numeric_limits<double>::infinity();
  }
  for (int i=0;i<3*n;i++) {
    matrix_c[i]=numeric_limits<double>::infinity();
  }
  for (int i=0;i<nsteps;i++) {
    log_error[i]=numeric_limits<double>::infinity();
  }
#endif

  int fi=0; //loop address of first stored array index
  int la=ncells;
  int ifirst=1; // loop address of first unknown array index
  int ilast=ncells-1;

  double *matrix=matrix_c-fi;
  double *residual=residual_c-fi;
  double *rhs=rhs_c-fi;
  double *solution=solution_c-fi;
  double *ap=ap_c-fi;
  double *cg_direction=cg_direction_c-fi;
  double *cg_preconditioner_solve=cg_preconditioner_solve_c-fi;
  double *cg_residual=cg_residual_c-fi;
//double *mg_residual=mg_residual_c-fi;
  double *true_solution=true_solution_c-fi;
  double *solution_increment=solution_increment_c-fi;
  double *matrix_copy=matrix_copy_c-fi;

  double gamma=numeric_limits<double>::infinity();

  {
//  TRACER_CALL(tr,"initialization");
//  steady-state heat equation (ie Laplace) with constant diffusion
    double *diffusion_coef=OPERATOR_NEW_BRACKET(double,ncells);
    switch (diffusion) {
      case RANDOM_DIFFUSION: {
        double c=1./(dx*static_cast<double>(RAND_MAX));
        for (int i=fi;i<la;i++) {
          diffusion_coef[i]=static_cast<double>(rand())*c;
        }
        break;
      }
      case CONSTANT_DIFFUSION:
      default: {
        double c=1./dx;
        for (int i=fi;i<la;i++) diffusion_coef[i]=c;
        break;
      }
    }
    switch (equation) {
      case LAPLACE_EQUATION: {
        for (int i=ifirst;i<=ilast;i++) {
          matrix[i]=-diffusion_coef[i-1];
          matrix[i+n]=diffusion_coef[i-1]+diffusion_coef[i];
          matrix[i+n2]=-diffusion_coef[i];
        }
        break;
      }
      case HEAT_EQUATION: { // backward Euler matrix
        F77NAME(dscal)(ilast-ifirst+1,dt,&diffusion_coef[ifirst],1);
        for (int i=ifirst;i<=ilast;i++) {
          matrix[i]=-diffusion_coef[i-1];
          matrix[i+n]=dx+diffusion_coef[i-1]+diffusion_coef[i];
          matrix[i+n2]=-diffusion_coef[i];
        }
        break;
      }
/*
      case SINH: {
        double off_diag=0.5*(-static_cast<double>(ncells)+dx);
        double diag=static_cast<double>(ncells)+dx;
        for (int i=fi;i<=la;i++) {
          matrix[i]=off_diag;
          matrix[i+n]=diag;
          matrix[i+n2]=off_diag;
          rhs[i]=0.;
        }

        int i=0;
        double bv=sinh(-1.);
        matrix[i+n]=0.5*diag;
        rhs[0]=matrix[i+n]*bv;
        rhs[1]=-matrix[i+1]*bv;
        matrix[i+1]=0.;
        matrix[i+n2]=0.;

        i=la;
        bv=sinh(1.);
        matrix[i+n]=0.5*diag;
        rhs[i]=matrix[i+n]*bv;
        rhs[i-1]=-matrix[i-1+n2]*bv;
        matrix[i-1+n2]=0.;
        matrix[i]=0.;
        break;
      }
*/
    }
    delete [] diffusion_coef; diffusion_coef=0;
#ifdef DEBUG
//  for (int i=fi;i<=la;i++) {
//    cout << "\tmatrix[" << i << "] = " << matrix[i] << " "
//         << matrix[i+n] << " " << matrix[i+n2] <<endl;
//  }
#endif
    memcpy(matrix_copy_c,matrix_c,3*n*sizeof(double));

//  Dirichlet bc at left and right
/*
    if (equation==SINH) {
      for (int i=fi;i<=la;i++) {
        true_solution[i]=sinh(-1.+dx*static_cast<double>(2*i));
      }
    } else {
*/
      for (int i=fi;i<=la;i++) {
        true_solution[i]=1.-static_cast<double>(i)*dx;
//      true_solution[i]=0.;
//      true_solution[i]=sin(M_PI*static_cast<double>(i)*dx);
      }
//  }
#ifdef DEBUG
//for (int i=fi;i<=la;i++) {
//  cout << "\ttrue_solution[" << i << "] = " << true_solution[i] <<endl;
//}
#endif

//  if (equation!=SINH) {
      F77_NAME(matrix_multiply)(fi,la,ifirst,ilast,matrix,true_solution,
        rhs);
//  }
#ifdef DEBUG
//  for (int i=fi;i<=la;i++) {
//    cout << "\trhs[" << i << "] = " << rhs[i] <<endl;
//  }
#endif

/*
    if (equation==SINH) {
      solution[fi]=sinh(-1.);
      for (int i=ifirst;i<=ilast;i++) solution[i]=0.;
      solution[la]=sinh(1.);
    } else {
*/
//    random initial guess for iteration
      double c=1./static_cast<double>(RAND_MAX);
      for (int i=ifirst;i<=ilast;i++) {
        solution[i]=static_cast<double>(rand())*c;
      }
//  }
#ifdef DEBUG
//  for (int i=fi;i<=la;i++) {
//    cout << "\tsolution[" << i << "] = " << solution[i] <<endl;
//  }
#endif
  }

  computeError(ifirst,ilast,solution,true_solution,log_error[0]);
  log_error[0]/=log10;
#ifdef DEBUG
//cout << "\tlog_error[0] = " << log_error[0] << endl;
#endif

  switch (iteration) {
    case CONJUGATE_GRADIENTS_ITERATION: {
//    TRACER_CALL(tr,"CONJUGATE_GRADIENTS_ITERATION");
      F77_NAME(compute_residual)(fi,la,ifirst,ilast,matrix,rhs,solution,
        cg_residual);
      for (int i=ifirst;i<=ilast;i++) cg_residual[i]=-cg_residual[i];
      switch (preconditioner) {
        case IDENTITY_PRECONDITIONER:
          memcpy(cg_direction_c,cg_residual_c,n*sizeof(double));
          break;
        case JACOBI_PRECONDITIONER:
          F77_NAME(jacobi_preconditioner)(fi,la,ifirst,ilast,
            matrix,cg_residual, cg_direction);
          break;
        case MULTIGRID_PRECONDITIONER:
        default: {
          finest_level=OPERATOR_NEW Level(ncells,matrix,smoother_iterations,
            smoother,restriction_prolongation);
          finest_level->multigridStep(cg_residual,cg_direction);
          for (int i=ifirst;i<=ilast;i++) cg_direction[i]=-cg_direction[i];
          break;
        }
      }
#ifdef DEBUG
//    for (int i=ifirst;i<=ilast;i++) {
//      cout << "\tcg_residual[" << i << "] = " << cg_residual[i] <<endl;
//    }
//    for (int i=ifirst;i<=ilast;i++) {
//      cout << "\tcg_direction[" << i << "] = " << cg_direction[i] <<endl;
//    }
#endif
      gamma=F77NAME(ddot)(ilast-ifirst+1,&cg_residual[ifirst],1,
        &cg_direction[ifirst],1);
      break;
    }
    case MULTIGRID_ITERATION: {
      F77_NAME(compute_residual)(fi,la,ifirst,ilast,matrix,rhs,solution,
        cg_residual);
      finest_level=OPERATOR_NEW Level(ncells,matrix,smoother_iterations,
        smoother,restriction_prolongation);
#ifdef DEBUG
//    for (int i=ifirst;i<=ilast;i++) {
//      cout << "\tcg_residual[" << i << "] = " << cg_residual[i] <<endl;
//    }
//
//    finest_level->checkProlongationAndRestrictionRandom();
//    finest_level->checkProlongationAndRestriction();
//    finest_level->checkSmootherRandom();
//    finest_level->checkSmoother();

//    finest_level->checkCoarseGridProjectionRandom();
//    finest_level->checkCoarseGridProjection();
//    finest_level->checkMatrixSymmetryRandom();
//    finest_level->coarserLevel()->checkMatrixSymmetryRandom();
//    finest_level->checkMatrixSymmetry();
//    finest_level->checkVCycleSymmetryRandom();
//    finest_level->checkVCycleSymmetry();
#endif
      break;
    }
  }

//find min,max data values
  double xlo=static_cast<double>(fi);
  double xhi=static_cast<double>(la);
  double ulo=0.;
  double uhi=1.;
  double elo=log_error[0];
  double ehi=log_error[0];

  {
//  TRACER_CALL(tr0,"initialize graphics");
//  setup interactive graphics
    Palette pal;
    XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
    double winsize=0.5;
    char label[LENGTH_NAME];
    snprintf(label,LENGTH_NAME,"%s : solution",iteration_name[iteration]);
    XYGraphTool gt(label,"x","u",xlo,xhi,ulo,uhi,&cmap,0,winsize);

    plotSolution(gt,ifirst,ilast,solution,true_solution);

//  TimedObject integrate_timing("solve");
    times(&usage);
    start=usage.tms_utime;
    int stop_step=1;
    for (int step=0;step<nsteps;step++) {
//    TRACER_CALL(tr,"step");
      {
//      Timer timer(&solve_timing);

        switch (iteration) {
          case RICHARDSON_ITERATION:
            F77NAME(richardson)(fi,la,ifirst,ilast, mu, 
              matrix,rhs, solution, residual);
            break;
//        case JACOBI:
//          F77NAME(jacobi)(fi,la,ifirst,ilast,
//            matrix,rhs, solution, residual);
//          break;
          case JACOBI_OMEGA_ITERATION:
            F77_NAME(jacobi_omega)(fi,la,ifirst,ilast, omega,
              matrix,rhs, solution, residual);
            break;
          case GAUSS_SEIDEL_ITERATION:
            F77_NAME(gauss_seidel_to_fro)(fi,la,ifirst,ilast,0,
              matrix,rhs, solution);
            break;
          case GAUSS_SEIDEL_RED_BLACK_ITERATION:
            F77_NAME(gauss_seidel_red_black)(fi,la,ifirst,ilast, step,
              matrix,rhs, solution);
            break;
          case GAUSS_SEIDEL_TO_FRO_ITERATION:
              F77_NAME(gauss_seidel_to_fro)(fi,la,ifirst,ilast, step,
              matrix,rhs, solution);
            break;
          case SOR_ITERATION:
            F77NAME(sor)(fi,la,ifirst,ilast, omega,
              matrix,rhs, solution);
            break;
          case CONJUGATE_GRADIENTS_ITERATION: {
            preconditionedConjugateGradients(fi,la,ifirst,ilast,gamma,
              matrix, cg_direction,cg_residual,solution,
              ap,cg_preconditioner_solve);
            break;
          }
          case MULTIGRID_ITERATION:
          default: {
            finest_level->multigridStep(cg_residual,solution_increment);
            for (int i=ifirst;i<=ilast;i++) {
              solution[i]-=solution_increment[i];
            }
            finest_level->computeResidual(solution,rhs,cg_residual);
            break;
          }
        }
      }
#ifdef DEBUG
//    cout << "\n\tstep = " << step << endl;
//    for (int i=ifirst;i<=ilast;i++) {
//      cout << "\tsolution[" << i << "] = " << solution[i] <<endl;
//    }
#endif
      computeError(ifirst,ilast,solution,true_solution,log_error[step+1]);
      log_error[step+1]/=log10;
#ifdef DEBUG
//    cout << "\tlog_error[" << step+1 << "] = " << log_error[step+1]
//         << endl;
#endif

      elo=min(elo,log_error[step+1]);
      ehi=max(ehi,log_error[step+1]);

      plotSolution(gt,ifirst,ilast,solution,true_solution);

//    REAL elapsed=static_cast<REAL>(usage.tms_utime-start)
//                /static_cast<REAL>(sysconf(_SC_CLK_TCK));

      if (step+1==stop_step) {
        cout << "\tstep number = " << stop_step << endl;
        XYGraphTool::WINDOW_TYPE::QuitButton qb;
        stop_step *= 10;
      }
    } // end of time step loop
    snprintf(label,LENGTH_NAME,"%s : error",iteration_name[iteration]);
    XYGraphTool gte(label,"iteration number",
      "log10( error )",0.,nsteps,elo,ehi,&cmap,0,winsize);
    gte.setbgColor("white");
    gte.setfgColor("black");
    gte.drawAxes();
    int step=0;
    gte.setfgColor("blue");
    gte.movePen(static_cast<double>(step),log_error[step]);
    for (step=1;step<=nsteps;step++) {
      gte.drawLine(static_cast<double>(step),log_error[step]);
    }
    gte.flush();
    XYGraphTool::WINDOW_TYPE::QuitButton qb;
  } // pal,cmap,gt go out of scope here 
//integrate_timing.printOn(cout);

//since x, u, flux were created with operator new, must delete
  delete [] log_error;
  delete [] matrix_copy_c;
  delete [] solution_increment_c;
  delete [] true_solution_c;
  delete [] solution_c;
  delete [] rhs_c;
  delete [] residual_c;
  delete [] matrix_c;
  delete [] ap_c;
  delete [] cg_direction_c;
  delete [] cg_preconditioner_solve_c;
  delete [] cg_residual_c;
//delete [] mg_residual_c;
  if (finest_level!=0) delete finest_level; finest_level=0;
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
    cout << setprecision(16);
//  set machine-dependent constants for Fortran
    F77NAME(machine).roundoff=DBL_EPSILON;
    F77NAME(machine).small=DBL_MIN;
    F77NAME(machine).huge=DBL_MAX;
    F77NAME(machine).undefind=numeric_limits<double>::infinity();

//  define input parameters and bounds
    param_list=OPERATOR_NEW GUI_INPUT_PARAMETER_LIST_TYPE("Main List");
    { const char *group="Numerical Method Parameters";
      param_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
        iequation,"equation",equation_name[0],equation_name[1]
        /* equation_name[2] */ ));
      param_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
        idiffusion,"diffusion",diffusion_name[0],diffusion_name[1]));
      param_list->append(OPERATOR_NEW GUIInputParameter<double>(
        decay_number,"heat equation decay_number",0.,DBL_MAX,group));
    }
    { const char *group="Numerical Method Parameters";
      param_list->append(OPERATOR_NEW GUIInputParameter<int>(ncells,
        "ncells",1,INT_MAX,group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<int>(nsteps,
        "nsteps",1,INT_MAX,group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<double>(mu_in,
        "Richardson mu",2.,DBL_MAX,group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<double>(omega,
        "Jacobi/SOR omega",0.,2.,group) );
      param_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
        iiteration,"iteration",iteration_name[0],iteration_name[1],
        iteration_name[2],iteration_name[3],iteration_name[4],
        iteration_name[5],iteration_name[6],iteration_name[7]));
      param_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
        ipreconditioner,"CG preconditioner",preconditioner_name[0],
        preconditioner_name[1],preconditioner_name[2]));
      param_list->append(OPERATOR_NEW GUIInputParameter<int>(
        smoother_iterations,"MG smoother_iterations",1,INT_MAX,group) );
      param_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
        ismoother,"MG smoother",Level::smoother_name[0],
        Level::smoother_name[1],Level::smoother_name[2]));
      param_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
        irestriction_prolongation,"MG restriction/prolongation",
        Level::restriction_prolongation_name[0],
        Level::restriction_prolongation_name[1]));
    }
    { const char *group="Graphics Parameters";
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
