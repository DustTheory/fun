#include <float.h>
#include <fstream>
#include <iostream>
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
BOOLEAN skip_gui=FALSE;
int ncells[2]={8,8};
int nsteps=10;
int nbands=2;
int smoother_iterations=1;
double decay=1.;
double peclet=1.;
double mu=4.;
double omega=1.;
double tol=1.e-2;
double winsize=0.5;
int ncontours=5;
char* display_name=0;

int nelt=0;
int fc[2],lc[2],ifirst[2],ilast[2];
double *x=0;
double *ax=0;
double *matrix=0;

bool run_jacobi=false;
bool run_gauss_seidel_red_black=false;
bool run_gauss_seidel_to_fro=false;
bool run_sor=false;
bool run_gmres=false;
bool run_orthomin=false;

enum ITERATION{JACOBI_OMEGA_ITERATION,GAUSS_SEIDEL_RED_BLACK_ITERATION,
  GAUSS_SEIDEL_TO_FRO_ITERATION,SOR_ITERATION,
  GMRES_ITERATION,ORTHOMIN_ITERATION};
ITERATION iteration=ORTHOMIN_ITERATION;
int iiteration=iteration;
static char *iteration_name[6]={
  "jacobi_omega","gauss_seidel_red_black","gauss_seidel_to_fro","sor",
  "gmres","orthomin"
};

static char *color_name[11]={
  "red","green","blue","cyan","yellow","magenta","brown","purple",
  "grey","pink","tan"
};

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void checkMainInput() {
  iteration=static_cast<ITERATION>(iiteration);
  if (iteration==JACOBI_OMEGA_ITERATION) {
    omega=min(omega,1.);
  }
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void contour_corner(AreaGraphTool &g,const int *fc,const int *lc,
const int *ifirst,const int *ilast,double val,const double *w) {
//TRACER_CALL(tr,"contour_corner");
#ifdef DEBUG
//cout << "\tw = " << endl;
//F77NAME(rbugcorn)(DIM_ARG(fc),DIM_ARG(lc),
//  DIM_ARG(ifirst),DIM_ARG(ilast),w);
#endif
  double dif[5],xc[5][2],pos,x4[4][2];
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
            x4[i][id]=xc[i][id]-pos*(xc[i+1][id]-xc[i][id]);
          }
        }
      }
#ifdef DEBUG
//    cout << "\tcount,first,last = " << count << " " << first << " "
//         << last << endl;
#endif
      if (count==4) {
        g.movePen(x4[0][0],x4[0][1]);
        g.drawLine(x4[2][0],x4[2][1]);
        g.movePen(x4[1][0],x4[1][1]);
        g.drawLine(x4[3][0],x4[3][1]);
      } else if (count==2) {
        g.movePen(x4[first][0],x4[first][1]);
        g.drawLine(x4[last][0],x4[last][1]);
      }
    }
  }
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void F77NAME(mymatvec)(const int &n,const double *dlap_x,
double *dlap_ax,const int &nelt,const int *ia,const int *ja,
const double *a,const int &isym) {
//TRACER_CALL(tr,"mymatvec");
  F77_NAME(copy_dlap_solution)(DIM_ARG(fc),DIM_ARG(lc),
    DIM_ARG(ifirst),DIM_ARG(ilast), dlap_x,x);
  F77_NAME(matrix_multiply)(DIM_ARG(fc),DIM_ARG(lc),
    DIM_ARG(ifirst),DIM_ARG(ilast),matrix,x, ax);
  F77_NAME(setup_dlap_system)(DIM_ARG(fc),DIM_ARG(lc),
    DIM_ARG(ifirst),DIM_ARG(ilast), ax, dlap_ax);
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void F77NAME(mymsolve)(const int &n,const double *r,double *z,
const int &nelt,const int *ia,const int *ja,const double *a,
const int &isym,double *rwork,int *iwork) {
//TRACER_CALL(tr,"mymsolve");
  memcpy(z,r,n*sizeof(double)); // identity preconditioner
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void runOnce(int *nc) {
  TRACER_CALL(tr,"runOnce");
//if (run_richardson) iteration=RICHARDSON_ITERATION;
  if (run_jacobi) iteration=JACOBI_OMEGA_ITERATION;
  else if (run_gauss_seidel_red_black) {
    iteration=GAUSS_SEIDEL_RED_BLACK_ITERATION;
  }
  else if (run_gauss_seidel_to_fro) {
    iteration=GAUSS_SEIDEL_TO_FRO_ITERATION;
  }
  else if (run_sor) iteration=SOR_ITERATION;
  else if (run_gmres) iteration=GMRES_ITERATION;
  else if (run_orthomin) iteration=ORTHOMIN_ITERATION;
//array bounds for Fortran calls
  double rlo[2],rhi[2];
  int ncorners=1;
  int total_ncells=1;
  for (int id=0;id<2;id++) {
    fc[id]=0; // includes boundary values
    lc[id]=nc[id];
    ifirst[id]=1; // excludes Dirichlet boundary values
    ilast[id]=nc[id]-1;
    ncorners*=lc[id]-fc[id]+1;
    total_ncells*=lc[id]-fc[id];
    rlo[id]=static_cast<double>(fc[id]);
    rhi[id]=static_cast<double>(lc[id]);
  }
  double log10=log(10.);
#ifdef DEBUG
//cout << "\tfc = " << fc[0] << " " << fc[1] << endl;
//cout << "\tlc = " << lc[0] << " " << lc[1] << endl;
//cout << "\tifirst = " << ifirst[0] << " " << ifirst[1] << endl;
//cout << "\tilast = " << ilast[0] << " " << ilast[1] << endl;
//cout << "\tncorners = " << ncorners << endl;
//cout << "\ttotal_ncells = " << total_ncells << endl;
#endif

  double *diagonal_block=OPERATOR_NEW_BRACKET(double,lc[0]-fc[0]+1);
  matrix=OPERATOR_NEW_BRACKET(double,ncorners*9);
  double *residual=OPERATOR_NEW_BRACKET(double,ncorners);
  double *rhs=OPERATOR_NEW_BRACKET(double,ncorners);
  double *solution=OPERATOR_NEW_BRACKET(double,ncorners);
  x=OPERATOR_NEW_BRACKET(double,ncorners);
  ax=OPERATOR_NEW_BRACKET(double,ncorners);
  double *true_solution=OPERATOR_NEW_BRACKET(double,ncorners);
  double *dlap_b=OPERATOR_NEW_BRACKET(double,ncorners);
  double *dlap_x=OPERATOR_NEW_BRACKET(double,ncorners);
  double *dlap_a=OPERATOR_NEW_BRACKET(double,ncorners*9);
  int *dlap_ia=OPERATOR_NEW_BRACKET(int,ncorners*9);
  int *dlap_ja=OPERATOR_NEW_BRACKET(int,ncorners*9);
  double *log_error=OPERATOR_NEW_BRACKET(double,nsteps);

  int n=(ilast[0]-ifirst[0]+1)*(ilast[1]-ifirst[1]+1);
  int nnz=4*(ilast[0]-ifirst[0]+1)*(ilast[1]-ifirst[1]+1)
         -3*(ilast[0]-ifirst[0]+1)-3*(ilast[1]-ifirst[1]+1)+2;
  double *a=OPERATOR_NEW_BRACKET(double,nnz+n*nbands);
  double *diag=OPERATOR_NEW_BRACKET(double,n);
  int *col_ptr=OPERATOR_NEW_BRACKET(int,n+1);
  int *row_ind=OPERATOR_NEW_BRACKET(int,nnz+n*nbands);

  double gamma=numeric_limits<double>::infinity();
  double residual_max=numeric_limits<double>::infinity();
  double rhs_max=numeric_limits<double>::infinity();

  {
//  TRACER_CALL(tr,"runOnce initialization");
#ifdef INDEF
//  initialize arrays to help IEEE exception handling catch unassigned
    for (int i=0;i<ncorners;i++) {
      residual[i]=0.;
      x[i]=0.;
      ax[i]=0.;
      rhs[i]=numeric_limits<double>::infinity();
      solution[i]=numeric_limits<double>::infinity();
      true_solution[i]=numeric_limits<double>::infinity();
    }
    for (int i=0;i<9*ncorners;i++) {
      matrix[i]=numeric_limits<double>::infinity();
      dlap_a[i]=numeric_limits<double>::infinity();
      dlap_ia[i]=INT_MAX;
      dlap_ja[i]=INT_MAX;
    }
    for (int i=0;i<nsteps;i++) log_error[i]=numeric_limits<double>::infinity();
#endif
    double c=1./static_cast<double>(RAND_MAX);
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
#ifdef DEBUG
//  cout << "\ttrue_solution = " << endl;
//  F77NAME(rbugcell)(DIM_ARG(fc),DIM_ARG(lc),
//    DIM_ARG(ifirst),DIM_ARG(ilast),true_solution);
#endif
    F77_NAME(setup_convection_diffusion)(DIM_ARG(fc),DIM_ARG(lc),
      DIM_ARG(ifirst),DIM_ARG(ilast),decay,peclet,
      true_solution,matrix,solution);
    if (iteration==GMRES_ITERATION || iteration==ORTHOMIN_ITERATION) {
      F77_NAME(setup_dlap_matrix)(DIM_ARG(fc),DIM_ARG(lc),
        DIM_ARG(ifirst),DIM_ARG(ilast), nelt,
        matrix, dlap_a,dlap_ia,dlap_ja);
    } 
#ifdef DEBUG
//  for (int i=-2;i<=2;i++) {
//    cout << "\tmatrix[" << i << "] = " << endl;
//    F77NAME(rbugcell)(DIM_ARG(fc),DIM_ARG(lc),
//      DIM_ARG(ifirst),DIM_ARG(ilast),matrix+(i+2)*ncorners);
//  }
#endif
//  since residual=0, the following sets rhs=matrix * true_solution
    F77_NAME(compute_residual)(DIM_ARG(fc),DIM_ARG(lc),
      DIM_ARG(ifirst),DIM_ARG(ilast),rhs_max,
      matrix,residual,true_solution,rhs);
#ifdef DEBUG
//  cout << "\trhs = " << endl;
//  F77NAME(rbugcell)(DIM_ARG(fc),DIM_ARG(lc),
//    DIM_ARG(ifirst),DIM_ARG(ilast),rhs);
//  cout << "\tsolution = " << endl;
//  F77NAME(rbugcell)(DIM_ARG(fc),DIM_ARG(lc),
//    DIM_ARG(ifirst),DIM_ARG(ilast),solution);
#endif
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
      case ORTHOMIN_ITERATION: {
//      TRACER_CALL(tr,"runOnce step domn");
        F77_NAME(setup_dlap_system)(DIM_ARG(fc),DIM_ARG(lc),
          DIM_ARG(ifirst),DIM_ARG(ilast), rhs, dlap_b);
        F77_NAME(setup_dlap_system)(DIM_ARG(fc),DIM_ARG(lc),
          DIM_ARG(ifirst),DIM_ARG(ilast), solution, dlap_x);
        int dlap_n=(ilast[0]-ifirst[0]+1)*(ilast[1]-ifirst[1]+1);
#ifdef DEBUG
//      for (int i=0;i<dlap_n;i++) {
//        cout << "\tdlap_b[" << i << "] = " << dlap_b[i] << endl;
//      }
//      for (int i=0;i<dlap_n;i++) {
//        cout << "\tdlap_x[" << i << "] = " << dlap_x[i] << endl;
//      }
//      for (int i=0;i<nelt;i++) {
//        cout << "\tdlap_a[" << dlap_ia[i] << "," << dlap_ja[i]
//             << "] = " << dlap_a[i] << endl;
//      }
#endif
        int isym=0;
        int itol=1;
        double tol=1.e3*F77NAME(d1mach)(3);
        int itmax=step+1;
        int iter=INT_MAX;
        double err=numeric_limits<double>::infinity();
        int ierr=INT_MAX;
        int iunit=0;
        int nsave=step+1;
        double *r=OPERATOR_NEW_BRACKET(double,dlap_n);
        double *z=OPERATOR_NEW_BRACKET(double,dlap_n);
        double *p=OPERATOR_NEW_BRACKET(double,dlap_n*(nsave+1));
        double *ap=OPERATOR_NEW_BRACKET(double,dlap_n*(nsave+1));
        double *emap=OPERATOR_NEW_BRACKET(double,dlap_n*(nsave+1));
        double *dz=OPERATOR_NEW_BRACKET(double,dlap_n);
        double *csav=OPERATOR_NEW_BRACKET(double,nsave);
        double *rwork=0;
        int *iwork=0;
        F77NAME(domn)(dlap_n,dlap_b,dlap_x,nelt,
          dlap_ia,dlap_ja,dlap_a,isym,F77NAME(mymatvec),
          F77NAME(mymsolve),nsave,itol,tol,itmax,iter,err,ierr,iunit,
          r,z,p,ap,emap,dz,csav,rwork,iwork);
        residual_max=err;
#ifdef DEBUG
//      cout << "\tstep = " << step << endl;
//      cout << "\titer = " << iter << endl;
//      cout << "\terr = " << err << endl;
//      cout << "\tierr = " << ierr << endl;
#endif
        delete [] r; r=0;
        delete [] z; z=0;
        delete [] p; p=0;
        delete [] ap; ap=0;
        delete [] emap; emap=0;
        delete [] dz; dz=0;
        delete [] csav; csav=0;
        F77_NAME(copy_dlap_solution)(DIM_ARG(fc),DIM_ARG(lc),
          DIM_ARG(ifirst),DIM_ARG(ilast), dlap_x, solution);
        break;
      }
      case GMRES_ITERATION: 
      default: {
//      TRACER_CALL(tr,"runOnce step gmres");
// this is really inefficient: we start from the beginning each step
//  but dgmres doesn't let us see intermediate solution values, just errors
        F77_NAME(setup_dlap_system)(DIM_ARG(fc),DIM_ARG(lc),
          DIM_ARG(ifirst),DIM_ARG(ilast), rhs, dlap_b);
        F77_NAME(setup_dlap_system)(DIM_ARG(fc),DIM_ARG(lc),
          DIM_ARG(ifirst),DIM_ARG(ilast), solution, dlap_x);
        int dlap_n=(ilast[0]-ifirst[0]+1)*(ilast[1]-ifirst[1]+1);
#ifdef DEBUG
//      for (int i=0;i<dlap_n;i++) {
//        cout << "\tdlap_b[" << i << "] = " << dlap_b[i] << endl;
//      }
//      for (int i=0;i<dlap_n;i++) {
//        cout << "\tdlap_x[" << i << "] = " << dlap_x[i] << endl;
//      }
//      for (int i=0;i<nelt;i++) {
//        cout << "\tdlap_a[" << dlap_ia[i] << "," << dlap_ja[i]
//             << "] = " << dlap_a[i] << endl;
//      }
#endif
        int isym=0;
        int itol=0;
        double tol=0.;
        int itmax=step+1;
        int iter=INT_MAX;
        double err=numeric_limits<double>::infinity();
        int ierr=INT_MAX;
        int iunit=0;
        double *sb=0;
        double *sx=0;
        int maxl=step+1;
        int lrgw=max(131+16*dlap_n,2+dlap_n*(maxl+6)+maxl*(maxl+3));
        double *rgwk=OPERATOR_NEW_BRACKET(double,lrgw);
        int ligw=20;
        int *igwk=OPERATOR_NEW_BRACKET(int,ligw);
        igwk[0]=maxl;
        igwk[1]=maxl;
        igwk[2]=0;
        igwk[3]=0;
        igwk[4]=-1;
        double *rwork=0;
        int *iwork=0;
        F77NAME(dgmres)(dlap_n,dlap_b,dlap_x,
          nelt,dlap_ia,dlap_ja,dlap_a,isym,
          F77NAME(mymatvec),F77NAME(mymsolve),itol,tol,itmax,iter,
          err,ierr,iunit,sb,sx,rgwk,lrgw,igwk,ligw,rwork,iwork);
        residual_max=err;
        delete [] rgwk; rgwk=0;
        delete [] igwk; igwk=0;
        F77_NAME(copy_dlap_solution)(DIM_ARG(fc),DIM_ARG(lc),
          DIM_ARG(ifirst),DIM_ARG(ilast), dlap_x, solution);
        break;
      }
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
  cout << "\titeration converged in " << step << " steps" << endl;
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
  delete [] diagonal_block;
  delete [] log_error;
  delete [] solution;
  delete [] true_solution;
  delete [] rhs;
  delete [] residual;
  delete [] a;
  delete [] diag;
  delete [] col_ptr;
  delete [] row_ind;

  delete [] matrix; matrix=0;
  delete [] x; x=0;
  delete [] ax; ax=0;
  delete [] dlap_b; dlap_b=0;
  delete [] dlap_x; dlap_x=0;
  delete [] dlap_a; dlap_a=0;
  delete [] dlap_ia; dlap_ia=0;
  delete [] dlap_ja; dlap_ja=0;
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double runSimulation(int *nc,double &log_error,double &log_it,
double &log_time) {
  TRACER_CALL(tr,"runSimulation");
#ifdef DEBUG
  cout << "\titeration_name = " << iteration_name[iteration] << endl;
#endif
//array bounds for Fortran calls
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

  double *diagonal_block=OPERATOR_NEW_BRACKET(double,lc[0]-fc[0]+1);
  matrix=OPERATOR_NEW_BRACKET(double,ncorners*9);
  double *residual=OPERATOR_NEW_BRACKET(double,ncorners);
  double *rhs=OPERATOR_NEW_BRACKET(double,ncorners);
  double *solution=OPERATOR_NEW_BRACKET(double,ncorners);
  double *true_solution=OPERATOR_NEW_BRACKET(double,ncorners);
  double *dlap_b=OPERATOR_NEW_BRACKET(double,ncorners);
  double *dlap_x=OPERATOR_NEW_BRACKET(double,ncorners);
  double *dlap_a=OPERATOR_NEW_BRACKET(double,ncorners*9);
  int *dlap_ia=OPERATOR_NEW_BRACKET(int,ncorners*9);
  int *dlap_ja=OPERATOR_NEW_BRACKET(int,ncorners*9);

  double *x=OPERATOR_NEW_BRACKET(double,ncorners);
  double *ax=OPERATOR_NEW_BRACKET(double,ncorners);

  int n=(ilast[0]-ifirst[0]+1)*(ilast[1]-ifirst[1]+1);
  int nnz=4*(ilast[0]-ifirst[0]+1)*(ilast[1]-ifirst[1]+1)
         -3*(ilast[0]-ifirst[0]+1)-3*(ilast[1]-ifirst[1]+1)+2;
  double *a=OPERATOR_NEW_BRACKET(double,nnz+n*nbands);
  double *diag=OPERATOR_NEW_BRACKET(double,n);
  int *col_ptr=OPERATOR_NEW_BRACKET(int,n+1);
  int *row_ind=OPERATOR_NEW_BRACKET(int,nnz+n*nbands);

  double gamma=numeric_limits<double>::infinity();
  double residual_max=numeric_limits<double>::infinity();
  double rhs_max=numeric_limits<double>::infinity();

  {
//  TRACER_CALL(tr,"runSimulation initialization");
#ifdef INDEF
//  initialize arrays to help IEEE exception handling catch unassigned
    for (int i=0;i<ncorners;i++) {
      residual[i]=0.;
      rhs[i]=numeric_limits<double>::infinity();
      solution[i]=numeric_limits<double>::infinity();
      true_solution[i]=numeric_limits<double>::infinity();
    }
    for (int i=0;i<9*ncorners;i++) {
      matrix[i]=numeric_limits<double>::infinity();
      dlap_a[i]=numeric_limits<double>::infinity();
      dlap_ia[i]=INT_MAX;
      dlap_ja[i]=INT_MAX;
    }
#endif
    for (int i=0;i<ncorners;i++) {
      int i0=i%(nc[0]+1);
      int i1=i/(nc[0]+1);
//    true_solution[i]=
//       static_cast<double>(i0)/static_cast<double>(nc[0])
//      +static_cast<double>(i1)/static_cast<double>(nc[1]);
      true_solution[i]=
        sin(M_PI*static_cast<double>(i0)/static_cast<double>(nc[0]))
       *sin(M_PI*static_cast<double>(i1)/static_cast<double>(nc[1]));
      solution[i]=static_cast<double>(rand())
                 /static_cast<double>(RAND_MAX);
    }
    F77_NAME(setup_convection_diffusion)(DIM_ARG(fc),DIM_ARG(lc),
      DIM_ARG(ifirst),DIM_ARG(ilast),decay,peclet,
      true_solution,matrix,solution);
//  since residual=0, the following sets rhs=matrix * true_solution
    F77_NAME(compute_residual)(DIM_ARG(fc),DIM_ARG(lc),
      DIM_ARG(ifirst),DIM_ARG(ilast),residual_max,
      matrix,residual,true_solution,rhs);
    switch (iteration) {
      case ORTHOMIN_ITERATION:
      case GMRES_ITERATION:
        F77_NAME(setup_dlap_matrix)(DIM_ARG(fc),DIM_ARG(lc),
          DIM_ARG(ifirst),DIM_ARG(ilast), nelt,
          matrix, dlap_a,dlap_ia,dlap_ja);
        break;
    } 
  }

//find min,max data values

#ifdef DEBUG
//cout << "\tresidual_max,log_error[0] = " << residual_max << " "
//     << log_error[0] << endl;
#endif

  TimedObject iterate_timing("iterate");
  int step=0;
  double initial_residual_norm=numeric_limits<double>::infinity();
  {
    Timer timer(&iterate_timing);
//  TRACER_CALL(tr,"runSimulation loop");
    switch (iteration) {
      case JACOBI_OMEGA_ITERATION:
        for (;step<nsteps;step++) {
          F77_NAME(jacobi_omega)(DIM_ARG(fc),DIM_ARG(lc),
            DIM_ARG(ifirst),DIM_ARG(ilast), omega,residual_max,
            matrix,rhs, solution, residual);
          if (step==0) initial_residual_norm=residual_max;
          if (residual_max<=tol*rhs_max) break;
        }
        break;
      case GAUSS_SEIDEL_RED_BLACK_ITERATION:
        for (;step<nsteps;step++) {
          F77_NAME(gauss_seidel_red_black)(DIM_ARG(fc),DIM_ARG(lc),
            DIM_ARG(ifirst),DIM_ARG(ilast), step,residual_max,
            matrix,rhs, solution);
          if (step==0) initial_residual_norm=residual_max;
          if (residual_max<=tol*rhs_max) break;
        }
        break;
      case GAUSS_SEIDEL_TO_FRO_ITERATION:
        for (;step<nsteps;step++) {
          F77_NAME(gauss_seidel_to_fro)(DIM_ARG(fc),DIM_ARG(lc),
            DIM_ARG(ifirst),DIM_ARG(ilast), step,residual_max,
            matrix,rhs, solution);
          if (step==0) initial_residual_norm=residual_max;
          if (residual_max<=tol*rhs_max) break;
        }
        break;
      case SOR_ITERATION:
        for (;step<nsteps;step++) {
          F77NAME(sor)(DIM_ARG(fc),DIM_ARG(lc),
            DIM_ARG(ifirst),DIM_ARG(ilast), omega,residual_max,
            matrix,rhs, solution);
          if (step==0) initial_residual_norm=residual_max;
          if (residual_max<=tol*rhs_max) break;
        }
        break;
      case ORTHOMIN_ITERATION: {
        F77_NAME(setup_dlap_system)(DIM_ARG(fc),DIM_ARG(lc),
          DIM_ARG(ifirst),DIM_ARG(ilast), rhs, dlap_b);
        F77_NAME(setup_dlap_system)(DIM_ARG(fc),DIM_ARG(lc),
          DIM_ARG(ifirst),DIM_ARG(ilast), solution, dlap_x);
        int dlap_n=(ilast[0]-ifirst[0]+1)*(ilast[1]-ifirst[1]+1);
        int isym=0;
        int itol=1;
//      double tol=1.e3*F77NAME(d1mach)(3);
        int itmax=nsteps;
        int iter=INT_MAX;
        double err=numeric_limits<double>::infinity();
        int ierr=INT_MAX;
        int iunit=0;
        int nsave=10;
        double *r=OPERATOR_NEW_BRACKET(double,dlap_n);
        double *z=OPERATOR_NEW_BRACKET(double,dlap_n);
        double *p=OPERATOR_NEW_BRACKET(double,dlap_n*(nsave+1));
        double *ap=OPERATOR_NEW_BRACKET(double,dlap_n*(nsave+1));
        double *emap=OPERATOR_NEW_BRACKET(double,dlap_n*(nsave+1));
        double *dz=OPERATOR_NEW_BRACKET(double,dlap_n);
        double *csav=OPERATOR_NEW_BRACKET(double,nsave);
        double *rwork=0;
        int *iwork=0;
        F77NAME(domn)(dlap_n,dlap_b,dlap_x,nelt,dlap_ia,dlap_ja,dlap_a,
          isym,F77NAME(mymatvec),F77NAME(mymsolve),nsave,itol,tol,itmax,
          iter,err,ierr,iunit,r,z,p,ap,emap,dz,csav,rwork,iwork);
#ifdef DEBUG
        cout << "\titer = " << iter << endl;
        cout << "\terr = " << err << endl;
        cout << "\tierr = " << ierr << endl;
#endif
        delete [] r; r=0;
        delete [] z; z=0;
        delete [] p; p=0;
        delete [] ap; ap=0;
        delete [] emap; emap=0;
        delete [] dz; dz=0;
        delete [] csav; csav=0;
        F77_NAME(copy_dlap_solution)(DIM_ARG(fc),DIM_ARG(lc),
          DIM_ARG(ifirst),DIM_ARG(ilast), dlap_x, solution);
        break;
      }
      case GMRES_ITERATION:
      default: {
        F77_NAME(setup_dlap_system)(DIM_ARG(fc),DIM_ARG(lc),
          DIM_ARG(ifirst),DIM_ARG(ilast), rhs, dlap_b);
        F77_NAME(setup_dlap_system)(DIM_ARG(fc),DIM_ARG(lc),
          DIM_ARG(ifirst),DIM_ARG(ilast), solution, dlap_x);
        int dlap_n=ilast-ifirst+1;
        int isym=0;
        int itol=0;
//      double tol=0.;
        int itmax=nsteps;
        int iter=INT_MAX;
        double err=numeric_limits<double>::infinity();
        int ierr=INT_MAX;
        int iunit=0;
        double *sb=0;
        double *sx=0;
        int maxl=10;
        int nrmax=10;
        int lrgw=1+dlap_n*(maxl+6)+maxl*(maxl+3);
        double *rgwk=OPERATOR_NEW_BRACKET(double,lrgw);
        int ligw=20;
        int *igwk=OPERATOR_NEW_BRACKET(int,ligw);
        igwk[0]=maxl;
        igwk[1]=maxl;
        igwk[2]=0;
        igwk[3]=0;
        igwk[4]=nrmax;
        double *rwork=0;
        int *iwork=0;
        F77NAME(dgmres)(dlap_n,dlap_b,dlap_x,
          nelt,dlap_ia,dlap_ja,dlap_a,isym,
          F77NAME(mymatvec),F77NAME(mymsolve),itol,tol,itmax,iter,
          err,ierr,iunit,sb,sx,rgwk,lrgw,igwk,ligw,rwork,iwork);
        delete [] rgwk; rgwk=0;
        delete [] igwk; igwk=0;
        F77_NAME(copy_dlap_solution)(DIM_ARG(fc),DIM_ARG(lc),
          DIM_ARG(ifirst),DIM_ARG(ilast), dlap_x, solution);
        break;
      }
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

  delete [] diagonal_block;
  delete [] solution;
  delete [] rhs;
  delete [] residual;
  delete [] a;
  delete [] diag;
  delete [] col_ptr;
  delete [] row_ind;

  delete [] matrix; matrix=0;
  delete [] x; x = 0;
  delete [] ax; ax = 0;
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
      iteration=JACOBI_OMEGA_ITERATION;
      break;
    case 2:
      iteration=GAUSS_SEIDEL_RED_BLACK_ITERATION;
      break;
    case 4:
      iteration=GAUSS_SEIDEL_TO_FRO_ITERATION;
      break;
    case 8:
      iteration=SOR_ITERATION;
      break;
    case 16:
      iteration=ORTHOMIN_ITERATION;
      break;
    case 32:
      iteration=GMRES_ITERATION;
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
void runMain(BOOLEAN /*called_before*/) {
//TRACER_CALL(tr,"runMain");
  if (ncells[0]>0 && ncells[1]>0) {
    runOnce(ncells);
    return;
  }

  unsigned int number_simulations=4;
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

  double *log_error_domn=OPERATOR_NEW_BRACKET(double,number_simulations);
  double *log_nits_domn=OPERATOR_NEW_BRACKET(double,number_simulations);
  double *log_time_domn=OPERATOR_NEW_BRACKET(double,number_simulations);

  double *log_error_dgmres=OPERATOR_NEW_BRACKET(double,number_simulations);
  double *log_nits_dgmres=OPERATOR_NEW_BRACKET(double,number_simulations);
  double *log_time_dgmres=OPERATOR_NEW_BRACKET(double,number_simulations);

  double *log_size=OPERATOR_NEW_BRACKET(double,number_simulations);

  int nc[2]={10,10};
  double elo=numeric_limits<double>::infinity();
  double ehi=-numeric_limits<double>::infinity();
  double tlo=numeric_limits<double>::infinity();
  double thi=-numeric_limits<double>::infinity();
  double itlo=numeric_limits<double>::infinity();
  double ithi=-numeric_limits<double>::infinity();

  unsigned int runs=1;
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
  runs*=2;
  nc[0]=nc[1]=10;
  if (run_orthomin) {
    runMultipleSimulations(number_simulations,nc,runs,elo,ehi,tlo,thi,
      itlo,ithi,log_error_domn,log_time_domn,log_nits_domn,log_size);
  }
  runs*=2;
  nc[0]=nc[1]=10;
  if (run_gmres) {
    runMultipleSimulations(number_simulations,nc,runs,elo,ehi,tlo,thi,
      itlo,ithi,log_error_dgmres,log_time_dgmres,log_nits_dgmres,log_size);
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
  if (run_orthomin) {
    plotErrors(&gt1,&gt2,ORTHOMIN_ITERATION,
      number_simulations,itlo,ithi,tlo,thi,
      log_size,log_time_domn,log_nits_domn);
  }
  if (run_gmres) {
    plotErrors(&gt1,&gt2,GMRES_ITERATION,
      number_simulations,itlo,ithi,tlo,thi,
      log_size,log_time_dgmres,log_nits_dgmres);
  }
  XYGraphTool::WINDOW_TYPE::QuitButton qb;

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

  delete [] log_error_domn;
  delete [] log_nits_domn;
  delete [] log_time_domn;

  delete [] log_error_dgmres;
  delete [] log_nits_dgmres;
  delete [] log_time_dgmres;

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
    { char *group="Problem Parameters";
      param_list->append(OPERATOR_NEW GUIInputParameter<double>(decay,
        "decay number",0.,DBL_MAX,group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<double>(peclet,
        "cell peclet number",0.,DBL_MAX,group) );
    }
    { char *group="Numerical Method Parameters";
      param_list->append(OPERATOR_NEW GUIInputParameter<int>(ncells[0],
        "ncells[0]",0,INT_MAX,group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<int>(ncells[1],
        "ncells[1]",0,INT_MAX,group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<int>(nsteps,
        "nsteps",1,INT_MAX,group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<int>(nbands,
        "nbands",1,INT_MAX,group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<int>(
        smoother_iterations,"smoother_iterations",1,INT_MAX,group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<double>(omega,
        "omega",0.,2.,group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<double>(tol,
        "convergence tolerance",0.,1.,group) );
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
      param_list->append(OPERATOR_NEW GUIInputParameter<bool>(
       run_orthomin,"orthomin",false,true,group));
      param_list->append(OPERATOR_NEW GUIInputParameter<bool>(
       run_gmres,"gmres",false,true,group));
    }
    { char *group="Graphics Parameters";
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
        BOOLEAN found_main=FALSE;
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
