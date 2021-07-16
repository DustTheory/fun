#include <float.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdlib.h>
#include <math.h>
#include "Arch.H"
#include "Debug.H"
#include "GUIInputParameter.H"
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
#include "XColormap.H"
#include "GUI.H"
#include "GUIVirtualInput.H"
#endif

//#define LENGTH_NAME 80

extern "C" {
  void F77NAME(sor)(const int &fi,const int &la,
    const int &ifirst,const int &ilast, const double &omega,
    double *matrix,double *rhs, double *solution);
//void F77NAME(print_loc)(void *p) {
//  cout << p << endl;
//  cout << "\tdouble = " << *(double*) p << endl;
//  cout << "\tfloat = " << *(float*) p << endl;
//  cout << "\tint = " << *(int*) p << endl;
//}
}
struct machine_common {
  double roundoff,small,huge,undefind;
};
extern machine_common F77NAME(machine);

GUI_INPUT_PARAMETER_LIST_TYPE *param_list=0;
bool skip_gui=FALSE;
int ncells=10;
int nsteps=10;
int nomegas=1;
double omega=1.;
double winsize=0.5;
char* display_name=0;

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void checkMainInput() {
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void runMain(bool /*called_before*/) {
//array bounds for Fortran calls
  int fi=0;       // for cell-centered arrays, ie u
  int la=ncells-1;
  int ifirst=fi;
  int ilast=la;
  double log10=log(10.);

//since ncells is determined dynamically, 
//must use dynamic memory allocation
  double *matrix=OPERATOR_NEW_BRACKET(double,ncells*3);
  double *residual=OPERATOR_NEW_BRACKET(double,ncells);
  double *rhs=OPERATOR_NEW_BRACKET(double,ncells);
  double *solution=OPERATOR_NEW_BRACKET(double,ncells);
  double *true_solution=OPERATOR_NEW_BRACKET(double,ncells);
  double *matrix_copy=OPERATOR_NEW_BRACKET(double,ncells*3);
  double *log_error=OPERATOR_NEW_BRACKET(double,nsteps+1);
  double *log_spectral_radius=OPERATOR_NEW_BRACKET(double,nomegas);

  {
//  TRACER_CALL(tr,"initialization");
#ifdef INDEF
//  initialize arrays to help IEEE exception handling catch unassigned
    for (int i=0;i<ncells;i++) {
      residual[i]=numeric_limits<double>::infinity();
      rhs[i]=numeric_limits<double>::infinity();
      solution[i]=numeric_limits<double>::infinity();
    }
    for (int i=0;i<3*ncells;i++) {
      matrix[i]=numeric_limits<double>::infinity();
    }
#endif
    for (int i=0;i<ncells;i++) {
      matrix[i]=-1.;
      rhs[i]=0.;
      solution[i]=static_cast<double>(rand())
                 /static_cast<double>(RAND_MAX);
    }
    rhs[0]=1.;
    for (int i=ncells;i<2*ncells;i++) {
      matrix[i]=2.;
    }
    for (int i=2*ncells;i<3*ncells;i++) {
      matrix[i]=-1.;
    }
  }
  {
#ifdef INDEF
    for (int i=0;i<ncells;i++) {
      true_solution[i]=numeric_limits<double>::infinity();
    }
    for (int i=0;i<nsteps;i++) {
      log_error[i]=numeric_limits<double>::infinity();
    }
#endif
    for (int i=0;i<3*ncells;i++) matrix_copy[i]=matrix[i];

//  factor and forward-solve
    true_solution[0]=rhs[0];
    for (int i=1;i<ncells;i++) {
      matrix_copy[i]/=matrix_copy[i+ncells-1];
      matrix_copy[i+ncells]-=
        matrix_copy[i]*matrix_copy[i+2*ncells-1];
      true_solution[i]=rhs[i]-matrix_copy[i]*true_solution[i-1];
    }
//  back-solve
    true_solution[ncells-1]/=matrix_copy[2*ncells-1];
    for (int i=ncells-2;i>=0;i--) {
      true_solution[i]=
        (true_solution[i]-matrix_copy[i+2*ncells]*true_solution[i+1])
        /matrix_copy[i+ncells];
    }
//  for (int i=0;i<ncells;i++) {
//    cout << "\tsolution[" << i << "] = " << true_solution[i] <<endl;
//  }

    log_error[0]=0.;
    for (int i=0;i<ncells;i++) {
      log_error[0]=
        max(log_error[0],abs(solution[i]-true_solution[i]));
    }
    log_error[0]=log(log_error[0])/log10;
  }
  {
#ifdef INDEF
    for (int i=0;i<nomegas;i++) {
      log_spectral_radius[i]=numeric_limits<double>::infinity();
    }
#endif
  }

//find min,max data values
  double xlo=static_cast<double>(fi);
  double xhi=static_cast<double>(la);
  double ulo=0.;
  double uhi=1.;
  double elo=log_error[0];
  double ehi=log_error[0];

  if (nomegas<=1) {
//  TRACER_CALL(tr0,"initialize graphics");
//  setup interactive graphics
    Palette pal;
    XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
    double winsize=0.5;
    XYGraphTool gt("sor : solution","x","u",xlo,xhi,ulo,uhi,
                   &cmap,0,winsize);

//  initialize graphics display
    gt.setbgColor("white");
    gt.setfgColor("black");
    gt.drawAxes();

//  draw numerical solution
    gt.setfgColor("blue");
    int ic=ifirst;
    gt.movePen(static_cast<double>(ic),solution[ic]);
    for (ic=ifirst+1;ic<=ilast;ic++) {
      gt.drawLine(static_cast<double>(ic),solution[ic]);
    }
//  force X to perform the requests
    gt.flush();

//  TimedObject integrate_timing("solve");
    for (int step=0;step<nsteps;step++) {
//    TRACER_CALL(tr,"step");
      {
//      TRACER_CALL(tr,"solve");
//      Timer timer(&solve_timing);

//      F77NAME(richardson)(fi,la,ifirst,ilast, mu, 
//        matrix,rhs, solution, residual);
//      F77NAME(jacobi)(fi,la,ifirst,ilast,
//        matrix,rhs, solution, residual);
        F77NAME(sor)(fi,la,ifirst,ilast, omega,
          matrix,rhs, solution);
      }
      log_error[step+1]=0.;
      for (int i=0;i<ncells;i++) {
        log_error[step+1]=
          max(log_error[step+1],abs(solution[i]-true_solution[i]));
      }
      log_error[step+1]=log(log_error[step+1])/log10;
      elo=min(elo,log_error[step+1]);
      ehi=max(ehi,log_error[step+1]);

//    initialize graphics display
      gt.newPage();
      gt.setbgColor("white");
      gt.setfgColor("black");
      gt.drawAxes();

//    draw numerical solution
      int ic=ifirst;
      gt.setfgColor("blue");
      gt.movePen(static_cast<double>(ic),solution[ic]);
      for (ic=ifirst+1;ic<=ilast;ic++) {
        gt.drawLine(static_cast<double>(ic),solution[ic]);
      }

//    force X to perform the requests
      gt.flush();
      if (step==0 || step==10 || step==100) { 
        cout << "step = " << step << endl;
        XYGraphTool::WINDOW_TYPE::QuitButton qb;
      }

    } // end of time step loop
    XYGraphTool gte("error","step","log10(error)",0.,nsteps,elo,ehi,
                   &cmap,0,winsize);
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
  else {
//  TRACER_CALL(tr,"nomegas > 1");
    double om=omega;
    double domega=(2.-om)/static_cast<double>(nomegas);
    double rholo=0.;
    double rhohi=0.;
    double optimal_omega=om;
    for (int k=0;k<nomegas;k++,om+=domega) {
      log_error[0]=0.;
      for (int i=0;i<ncells;i++) {
          solution[i]=static_cast<double>(rand())
                   /static_cast<double>(RAND_MAX);
//      log_error[0]=
//        max(log_error[0],abs(solution[i]-true_solution[i]));
        double e=solution[i]-true_solution[i];
        log_error[0]+=e*e;
      }
//    log_error[0]=log(log_error[0]);
      log_error[0]=log(sqrt(log_error[0]))/log10;
      double spectral_radius=0.;
      int count=0;
#ifdef DEBUG
//    cout << "\tnsteps = " << nsteps << endl;
#endif
      for (int step=0;step<nsteps;step++) {
        F77NAME(sor)(fi,la,ifirst,ilast, om,
          matrix,rhs, solution);

        log_error[step+1]=0.;
        for (int i=0;i<ncells;i++) {
//        log_error[step+1]=
//          max(log_error[step+1],abs(solution[i]-true_solution[i]));
          double e=solution[i]-true_solution[i];
          log_error[step+1]+=e*e;
        }
//      log_error[step+1]=log(log_error[step+1]);
        log_error[step+1]=log(sqrt(log_error[step+1]))/log10;
#ifdef DEBUG
//      cout << "\tlog_error[" << step+1 << "] = "
//           << log_error[step+1] << endl;
#endif
        if (log_error[step+1]<0. && log_error[step+1]>-14. 
        && step>10) {
          spectral_radius+=
             (log_error[step]-log_error[step-10])*log10
          /static_cast<double>(10);
          count++;
        }
#ifdef DEBUG
//      cout << "\tlog_error[" << step+1 << "] = "
//           << log_error[step+1] << endl;
//      cout << "\tcount = " << count << endl;
#endif
      }
//    cout << "\tcount[" << omega << "] = " << count << endl;
      if (count>0) {
        log_spectral_radius[k]=
          spectral_radius/static_cast<double>(count);
      } else log_spectral_radius[k]=rhohi;
      if (log_spectral_radius[k]<rholo) {
        rholo=log_spectral_radius[k];
        optimal_omega=om;
      }
      rhohi=max(rhohi,log_spectral_radius[k]);
    }
    cout << "\toptimal omega = " << optimal_omega << endl;

    Palette pal;
    XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
    double winsize=0.5;
    XYGraphTool gt("log spectral radius vs omega","omega",
      "log(rho)",omega,2.,rholo,rhohi,&cmap,0,winsize);
    gt.setbgColor("white");
    gt.setfgColor("black");
    gt.drawAxes();

    int ic=ifirst;
    gt.setfgColor("blue");
    om=omega;
    gt.movePen(om,log_spectral_radius[0]);
    for (int k=1;k<nomegas;k++,om+=domega) {
      gt.drawLine(om,log_spectral_radius[k]);
    }
    gt.flush();
    XYGraphTool::WINDOW_TYPE::QuitButton qb;
  }
//integrate_timing.printOn(cout);

//since x, u, flux were created with operator new, must delete
  delete [] log_spectral_radius;
  delete [] log_error;
  delete [] true_solution;
  delete [] solution;
  delete [] rhs;
  delete [] residual;
  delete [] matrix;
  delete [] matrix_copy;
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

//  define input parameters and bounds
    param_list=OPERATOR_NEW GUI_INPUT_PARAMETER_LIST_TYPE("Main List");
    { const char *group="Numerical Method Parameters";
      param_list->append(OPERATOR_NEW GUIInputParameter<int>(ncells,
        "ncells",1,INT_MAX,group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<int>(nsteps,
        "nsteps",1,INT_MAX,group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<int>(nomegas,
        "nomegas",1,INT_MAX,group) );
      param_list->append(OPERATOR_NEW GUIInputParameter<double>(omega,
        "omega",0,2.,group) );
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
