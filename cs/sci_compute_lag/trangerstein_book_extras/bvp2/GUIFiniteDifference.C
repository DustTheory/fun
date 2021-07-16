#include <iostream>
#include <limits>
//#include <math.h>

#include "Debug.H"
#include "GUIEnumInputParameter.H"
#include "GUIInputs.H"
#include "MemoryDebugger.H"
#include "NumPtr.H"
#include "SetTraps.H"
//#include "TimedObject.H"
//#include "Tracer.H"
#include "XYGraphTool.H"

#ifdef USE_GTK
#include <gtk/gtkmain.h>
#include "GTKGUI.H"
#include "GTKGUIVirtualInput.H"
#include "GTKWindow.H"
#else
#include "GUI.H"
#include "GUIVirtualInput.H"
#endif

GUI_INPUT_PARAMETER_LIST_TYPE *main_list=0;
bool skip_gui=FALSE;
char *display_name=0;
double winsize=0.5;

enum MESH{UNIFORM_MESH,RANDOM_MESH};
MESH mesh=UNIFORM_MESH;
int imesh=static_cast<int>(mesh);
const char *mesh_name[2]={"uniform mesh","random mesh"};

enum BOUNDARY_CONDITION{DIRICHLET_BOUNDARY_CONDITION,
  NEUMANN_BOUNDARY_CONDITION};
BOUNDARY_CONDITION boundary_condition=DIRICHLET_BOUNDARY_CONDITION;
int iboundary_condition=static_cast<int>(boundary_condition);
const char *boundary_condition_name[2]={"Dirichlet","Neumann"};

int ncells=1; 
int max_refine=2;
int nplot=512;

//#define LENGTH_NAME 80

extern "C" {
  void F77_NAME(analytical_solution)(const int &fi,const int &la,
    const double *mesh, double *solution);
  void F77_NAME(initialize_dirichlet)(const int &fi,const int &la,
    const double *mesh, double *matrix,double *rhs);
  void F77_NAME(initialize_neumann)(const int &fi,const int &la,
    const double *mesh, double *matrix,double *rhs);
  void F77NAME(solve)(const int &fi,const int &la,
    const double *matrix,const double *rhs_in_soln_out);
  void F77NAME(print_loc)(void *p) {
    cout << p << endl;
    cout << "\tdouble = " << *(double*) p << endl;
    cout << "\tfloat = " << *(float*) p << endl;
    cout << "\tint = " << *(int*) p << endl;
  }
}
struct machine_common {
  double roundoff,small,huge,undefind,pi;
};
extern machine_common F77NAME(machine);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void processCommandLine(int argc,char *argv[],char *display_name,
bool &skip_gui,ifstream &in_file) {
//TRACER_CALL(t,"processCommandLine");
  if (argc<2) skip_gui=FALSE;
  else {
    in_file.open(argv[1],ios::in);
    CHECK_TEST(!in_file.fail());
    if (in_file) {
      int i=2;
      while (i<argc) {
        if (strcmp(argv[i],"-d")==0) {
          display_name=OPERATOR_NEW_BRACKET(char,strlen(argv[++i])+1);
          strcpy(display_name,argv[i]);
        } else {
          cout << " >>>> error...invalid command line option:"
               << argv[i] << endl;
          abort();
        }
        i++;
      }   
    } else {
      cerr << "\ncannot open " << argv[1] << " for input" << endl;
      exit(1);
    }   
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void readMainInput(ifstream &in_file,InputParameterList *main_list,
bool &skip_gui) {
//TRACER_CALL(t,"readMainInput");
  bool found_main=FALSE;
  char name[LENGTH_NAME], comment[LENGTH_COMMENT];

  in_file.clear(ios::goodbit);
  in_file.seekg(0,ios::beg);
  while ( in_file >> setw(LENGTH_NAME) >> name ) {
    if ( strcmp(name,"Main") == 0 ) {
      in_file.getline( comment, LENGTH_COMMENT);
      found_main=TRUE;
      while ( in_file >> setw(LENGTH_NAME) >> name ) {
#ifdef DEBUG
//      cout << "\tname = " << name << endl;
#endif
        if ( strcmp(name,"end") == 0 ) break;
        else if (strcmp(name,"skip_gui") == 0) {
          int iskip_gui; // type bool not read correctly, so use int
          in_file >> iskip_gui;
          skip_gui= iskip_gui!=0;
        }
        else main_list->formattedRead(in_file,name);
        in_file.getline( comment, LENGTH_COMMENT);
#ifdef DEBUG
//      cout << "\tcomment = " << comment << endl;
#endif
      }
    } else in_file.getline(comment,LENGTH_COMMENT);
  }
#ifdef DEBUG
//main_list->printOn(cout);
#endif
  if ( !found_main ) skip_gui=FALSE;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void makeMainList(GUI_INPUT_PARAMETER_LIST_TYPE *&main_list) {
//TRACER_CALL(t,"makeMainList");
  main_list=OPERATOR_NEW GUI_INPUT_PARAMETER_LIST_TYPE("Main");

  { const char *group="Numerical Method Parameters";
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,imesh,
      "mesh",mesh_name[0],mesh_name[1]));
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
      iboundary_condition,"boundary_condition",boundary_condition_name[0],
      boundary_condition_name[1]));
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(ncells,
      "ncells",0,INT_MAX,group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(max_refine,
      "max_refine",2,INT_MAX,group) );
  }

  { const char *group="Graphics";
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(winsize,
      "winsize",0.,1.,group));
  }
#ifdef DEBUG
//main_list->printOn(cout);
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void checkMainInput() {
//TRACER_CALL(t,"checkMainInput");
  mesh=static_cast<MESH>(imesh);
  boundary_condition=static_cast<BOUNDARY_CONDITION>(iboundary_condition);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void runMain(bool /*called_before*/) {
  if (ncells>1) {
//  TRACER_CALL(tr0,"solve once");
//    array bounds for Fortran calls
    int fi=1;
    int la=ncells-1;

//  since ncells is determined dynamically, 
//  must use dynamic memory allocation
    double *solution=OPERATOR_NEW_BRACKET(double,nplot+1);
    double *x_plot=OPERATOR_NEW_BRACKET(double,nplot+1);
    double *matrix=OPERATOR_NEW_BRACKET(double,3*ncells);
    double *rhs=OPERATOR_NEW_BRACKET(double,ncells);
    double *x=OPERATOR_NEW_BRACKET(double,ncells+1);
//  double *matrix_copy=OPERATOR_NEW_BRACKET(double,3*ncells);
//  double *rhs_copy=OPERATOR_NEW_BRACKET(double,ncells);

#ifdef INDEF
//  initialize arrays to help IEEE exception handling catch unassigned
    for (int i=0;i<ncells;i++) {
      rhs[i]=numeric_limits<double>::infinity();
    }
    for (int i=0;i<3*ncells;i++) {
      matrix[i]=numeric_limits<double>::infinity();
    }
    for (int i=0;i<=ncells;i++) {
      x[i]=numeric_limits<double>::infinity();
    }
    for (int i=0;i<=nplot;i++) {
      solution[i]=numeric_limits<double>::infinity();
      x_plot[i]=numeric_limits<double>::infinity();
    }
#endif
    x_plot[0]=0.;
    double dx=1./static_cast<double>(nplot);
    for (int i=1;i<nplot;i++) {
      x_plot[i]=dx*static_cast<double>(i);
    }
    x_plot[nplot]=1.;

    x[0]=0.;
    switch (mesh) {
      case UNIFORM_MESH: {
        double dx=(boundary_condition==DIRICHLET_BOUNDARY_CONDITION ?
          1./static_cast<double>(ncells) :
          1./(static_cast<double>(ncells)-0.5) );
        for (int i=1;i<ncells;i++) {
          x[i]=dx*static_cast<double>(i);
        }
        break;
      }
      case RANDOM_MESH: {
        NumPtrC<double> mesh_points(ncells-1);
        for (int i=0;i<ncells-1;i++) {
          mesh_points[i]=drand48();
        }
        mesh_points.heapSort();
        for (int i=1;i<ncells;i++) {
          x[i]=mesh_points[i-1];
        }
        break;
      }
      default:
        OBSOLETE("unknown mesh type");
    }
    x[ncells]=1.;
    switch (boundary_condition) {
      case DIRICHLET_BOUNDARY_CONDITION:
        F77_NAME(initialize_dirichlet)(fi,la, x, matrix,rhs);
        break;
      case NEUMANN_BOUNDARY_CONDITION:
        F77_NAME(initialize_neumann)(fi,la, x, matrix,rhs);
        break;
    }
    F77_NAME(analytical_solution)(1,nplot-1, x_plot, solution);
#ifdef DEBUG
//  for (int i=0;i<=ncells;i++) {
//    cout << "\tx[" << i << "] = " << x[i] << endl;
//  }
//  for (int i=0;i<ncells;i++) {
//    cout << "\trhs[" << i << "] = " << rhs[i] << endl;
//  }
//  for (int i=0;i<=nplot;i++) {
//    cout << "\tx_plot,solution[" << i << "] = " << x_plot[i] << " "
//         << solution[i] << endl;
//  }
//  memcpy(matrix_copy,matrix,3*ncells*sizeof(double));
//  memcpy(rhs_copy,rhs,ncells*sizeof(double));
#endif

    F77NAME(solve)(fi,la, matrix,rhs);
#ifdef DEBUG
//  for (int i=0;i<ncells;i++) {
//    cout << "\trhs[" << i << "] = " << rhs[i] << endl;
//  }

//  cout << "\terror[1] = "
//       << rhs_copy[0]-matrix_copy[ncells]*rhs[0]
//                     -matrix_copy[2*ncells]*rhs[1]
//       << endl;
//  for (int i=1;i<ncells-1;i++) {
//    cout << "\terror[" << i+1 << "] = "
//         << rhs_copy[i]-matrix_copy[i]*rhs[i-1]
//                       -matrix_copy[i+ncells]*rhs[i]
//                       -matrix_copy[i+2*ncells]*rhs[i+1]
//         << endl;
//  }
//  int i=ncells-1;
//  cout << "\terror[" << i+1 << "] = "
//       << rhs_copy[i]-matrix_copy[i]*rhs[i-1]
//                     -matrix_copy[i+ncells]*rhs[i] << endl;
#endif

//  find min,max data values
    double xlo=0.;
    double xhi=1.;
    double ulo=solution[0];
    double uhi=solution[0];
    for (int i=1;i<nplot;i++) {
      double ui=solution[i];
      ulo=min(ulo,ui);
      uhi=max(uhi,ui);
//    double Ui=rhs[i-1];
//    ulo=min(ulo,min(ui,Ui));
//    uhi=max(uhi,max(ui,Ui));
    }
    double ui=solution[nplot];
    ulo=min(ulo,ui);
    uhi=max(uhi,ui);

//  setup interactive graphics
    Palette pal;
    XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
    double winsize=0.5;
    XYGraphTool gt("solution","x","u",xlo,xhi,ulo,uhi,&cmap,0,winsize);

//  initialize graphics display
    gt.setbgColor("white");
    gt.setfgColor("black");
    gt.drawAxes();

//  draw numerical solution
    gt.setfgColor("blue");
    gt.movePen(x[0],0.);
    for (int i=1;i<=ncells;i++) {
      gt.drawLine(x[i],rhs[i-1]);
    }

//  draw analytical solution
    gt.setfgColor("red");
    gt.movePen(x_plot[0],solution[0]);
    for (int i=1;i<=nplot;i++) {
      gt.drawLine(x_plot[i],solution[i]);
    }

//  force X to perform the requests
    gt.flush();
    XYGraphTool::WINDOW_TYPE::QuitButton qb;

    delete [] solution;
    delete [] x_plot;
    delete [] rhs;
    delete [] matrix;
    delete [] x;
  } // qb,pal,cmap,gt go out of scope here 
  else {
    TRACER_CALL(tr0,"mesh refinement study");
    double *log_h=OPERATOR_NEW_BRACKET(double,max_refine);
    double *error_infinity=OPERATOR_NEW_BRACKET(double,max_refine);
    double *error_1=OPERATOR_NEW_BRACKET(double,max_refine);
    double *error_2=OPERATOR_NEW_BRACKET(double,max_refine);

    int nc=4;
    double log10=log(10.);
    for (int k=0;k<max_refine;nc*=2,k++) {
      int fi=1;
      int la=nc-1;
      double *matrix=OPERATOR_NEW_BRACKET(double,3*nc);
      double *rhs=OPERATOR_NEW_BRACKET(double,nc);
      double *solution=OPERATOR_NEW_BRACKET(double,nc+1);
      double *x=OPERATOR_NEW_BRACKET(double,nc+1);
#ifdef DEBUG
//    F77_NAME(print_loc)(rhs);
#endif
      x[0]=0.;
      switch (mesh) {
        case UNIFORM_MESH: {
          double dx=(boundary_condition==DIRICHLET_BOUNDARY_CONDITION ?
            1./static_cast<double>(nc) :
            1./(static_cast<double>(nc)-0.5) );
          for (int i=1;i<nc;i++) {
            x[i]=dx*static_cast<double>(i);
          }
          break;
        }
        case RANDOM_MESH: {
          NumPtrC<double> mesh_points(nc-1);
          for (int i=0;i<nc-1;i++) {
            mesh_points[i]=drand48();
          }
          mesh_points.heapSort();
          for (int i=1;i<nc;i++) {
            x[i]=mesh_points[i-1];
          }
          break;
        }
        default:
          OBSOLETE("unknown mesh type");
      }
      x[nc]=1.;
      switch (boundary_condition) {
        case DIRICHLET_BOUNDARY_CONDITION:
          F77_NAME(initialize_dirichlet)(fi,la, x, matrix,rhs);
          break;
        case NEUMANN_BOUNDARY_CONDITION:
          F77_NAME(initialize_neumann)(fi,la, x, matrix,rhs);
          break;
      }
      F77_NAME(analytical_solution)(fi,la, x, solution);
      F77NAME(solve)(fi,la, matrix,rhs);

      double hmax=0.;
      error_infinity[k]=0.;
      error_1[k]=0.;
      error_2[k]=0.;
      for (int i=1;i<nc;i++) {
        double h=0.5*(x[i+1]-x[i-1]); 
        double ei=abs(rhs[i-1]-solution[i]);
        error_infinity[k]=max(error_infinity[k],ei);
        error_1[k]+=ei*h;
        error_2[k]+=ei*ei*h;
        hmax=max(hmax,h);
      }
      log_h[k]=-log(hmax)/log10;
      error_infinity[k]=log(error_infinity[k])/log10;
      error_1[k]=log(error_1[k])/log10;
      error_2[k]=log(sqrt(error_2[k]))/log10;
#ifdef DEBUG
//    cout << "\terror[" << hmax << "] = " << error_1[k] << " "
//         << error_2[k] << " " << error_infinity[k] << endl;
#endif
      delete [] solution;
      delete [] rhs;
      delete [] matrix;
      delete [] x;
    }
    nc/=2;
    Palette pal;
    XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
    {
      double elo=0.;
      double ehi=0.;
      for (int k=0;k<max_refine;k++) {
        elo=min(elo,error_infinity[k]);
        ehi=max(ehi,error_infinity[k]);
      }
      XYGraphTool gt("L^infinity error","log_10(h)","log_10(error)",
        log_h[0],log_h[max_refine-1],elo,ehi,&cmap,0,winsize);
      gt.setbgColor("white");
      gt.setfgColor("black");
      gt.drawAxes();
      gt.setfgColor("blue");
      gt.movePen(log_h[0],error_infinity[0]);
      cout << "\n infinity-norm errors:" << endl;
      for (int k=1;k<max_refine;k++) {
        gt.drawLine(log_h[k],error_infinity[k]);
        cout << "\tslope[" << k << "] = " 
             << (error_infinity[k]-error_infinity[k-1])
               /(log_h[k]-log_h[k-1]) << endl;
      }
      gt.flush();
      XYGraphTool::WINDOW_TYPE::QuitButton qb;
    }
    {
      double elo=0.;
      double ehi=0.;
      for (int k=0;k<max_refine;k++) {
        elo=min(elo,error_1[k]);
        ehi=max(ehi,error_1[k]);
      }
      XYGraphTool gt("L^1 error","log_10(h)","log_10(error)",
        log_h[0],log_h[max_refine-1],elo,ehi,&cmap,0,winsize);
      gt.setbgColor("white");
      gt.setfgColor("black");
      gt.drawAxes();
      gt.setfgColor("blue");
      gt.movePen(log_h[0],error_1[0]);
      cout << "\n 1-norm errors:" << endl;
      for (int k=1;k<max_refine;nc*=2,k++) {
        gt.drawLine(log_h[k],error_1[k]);
        cout << "\tslope[" << k << "] = " 
             << (error_1[k]-error_1[k-1])/(log_h[k]-log_h[k-1]) 
             << endl;
      }
      gt.flush();
      XYGraphTool::WINDOW_TYPE::QuitButton qb;
    }
    {
      double elo=0.;
      double ehi=0.;
      for (int k=0;k<max_refine;k++) {
        elo=min(elo,error_2[k]);
        ehi=max(ehi,error_2[k]);
      }
      XYGraphTool gt("L^2 error","log_10(h)","log_10(error)",
        log_h[0],log_h[max_refine-1],elo,ehi,&cmap,0,winsize);
      gt.setbgColor("white");
      gt.setfgColor("black");
      gt.drawAxes();
      gt.setfgColor("blue");
      gt.movePen(log_h[0],error_2[0]);
      cout << "\n 2-norm errors:" << endl;
      for (int k=1;k<max_refine;nc*=2,k++) {
        gt.drawLine(log_h[k],error_2[k]);
        cout << "\tslope[" << k << "] = " 
             << (error_2[k]-error_2[k-1])/(log_h[k]-log_h[k-1]) 
             << endl;
      }
      gt.flush();
      XYGraphTool::WINDOW_TYPE::QuitButton qb;
    }
    delete [] log_h;
    delete [] error_infinity;
    delete [] error_1;
    delete [] error_2;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void cleanup() {
//TRACER_CALL(t,"cleanup");
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//use this for things that need to be done at the end of GUI operations
void shutdown() {
//TRACER_CALL(t,"shutdown");
//things allocated in processCommandLine:
  if (display_name) delete display_name; display_name=0;

//things allocated in makeMainList:
  while (main_list->notEmpty()) {
    delete main_list->delAfter(0);
  }
  delete main_list; main_list=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int main(int argc,char* argv[]) {
//cout << "\tin main" << endl;
#ifdef DEBUG
//setTraps();
#endif
#if MEM_DEBUG
  MemoryDebugger md(1);
#endif
#ifdef USE_GTK
  GTKWindow::gtkInit(argc,argv);
#endif

//set machine-dependent constants for Fortran
  F77NAME(machine).roundoff=numeric_limits<double>::epsilon();
  F77NAME(machine).small=numeric_limits<double>::min();
  F77NAME(machine).huge=numeric_limits<double>::max();
  F77NAME(machine).undefind=numeric_limits<double>::infinity();
  F77NAME(machine).pi=M_PI;

  ifstream in_file;
  processCommandLine(argc,argv,display_name,skip_gui,in_file);
  makeMainList(main_list);
  if (in_file) {
    readMainInput(in_file,main_list,skip_gui);
    if (skip_gui) checkMainInput();
  }

  if (skip_gui) {
    runMain(0);
    cleanup();
    shutdown();
  } else {
#ifdef USE_GTK
    GTKGUI gui(argv[0],display_name,main_list,&runMain,
      &checkMainInput,&cleanup,&shutdown,TRUE);
#else
    GUI gui(argv[0],display_name,main_list,&runMain,&checkMainInput,
      &cleanup,&shutdown,TRUE);
#endif
    gui.createFileMenu();
    gui.createViewMenu();
    gui.createHelpMenu();
    gui.createMainWindow(argc,argv);
    gui.eventLoop();
  }
  return EXIT_SUCCESS;
}
#include "NumPtr.C"
INSTANTIATE_NUMPTRC(double);
//template int std::__cmath_power<int>(int,unsigned);
template class InputParameter<bool>;

//template class InputParameter<int>;
//template class InputParameter<double>;
//template class InputParameter<bool>;
