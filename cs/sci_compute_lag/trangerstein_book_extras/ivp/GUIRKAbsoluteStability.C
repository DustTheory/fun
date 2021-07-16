#include <cmath>
#include <iostream>
#include <limits>
#include <stdlib.h>

//#include "Debug.H"
#include "GUIEnumInputParameter.H"
#include "GUIInputs.H"
#include "MemoryDebugger.H"
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

enum SCHEME{FIRST_ORDER,SECOND_ORDER,THIRD_ORDER,FOURTH_ORDER};
SCHEME scheme=FIRST_ORDER;
int ischeme=static_cast<int>(scheme);

int nplot=128;
//int order=1;
complex<double> one_(1.,0.);
complex<double> half_(0.5,0.);
complex<double> sixth_(1./6.,0.);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double fcn3(double r, double theta,double val) {
  complex<double> hl=polar(r,theta)-one_;
  complex<double> rho=one_+hl*(one_+hl*(half_+hl*sixth_));
  return norm(rho)-val*val;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double fcn3deriv(double r, double theta) {
  complex<double> dhldr=polar(1.,theta);
  complex<double> hl=dhldr*r-one_;
  complex<double> rho=one_+hl*(one_+hl*(half_+hl*sixth_));
  complex<double> drhodr=dhldr*(one_+hl*(one_+hl*half_));
  return 2.*real(drhodr*conj(rho));
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double solveFcn3(double r,double theta,double val) {
  double eps=sqrt(numeric_limits<double>::epsilon());
  double f3=fcn3(r,theta,val);
  double f3eps=fcn3(r*(1.+eps),theta,val);
  for (int i=0;i<20;i++) {
    double f=fcn3(r,theta,val);
    double fprime=fcn3deriv(r,theta);
    r-=f/fprime;
    r=max(numeric_limits<double>::epsilon(),r);
    if (abs(f)<=1.e-14*r*abs(fprime)) break;
  }
  return r;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double fcn4(double r, double theta,double val) {
  complex<double> hl=polar(r,theta)-one_;
  return norm(one_+hl*(one_+hl*(half_+hl*(sixth_+hl/24.))))-val*val;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double fcn4deriv(double r, double theta) {
  complex<double> dhldr=polar(1.,theta);
  complex<double> hl=dhldr*r-one_;
  complex<double> rho=one_+hl*(one_+hl*(half_+hl*(sixth_+hl/24.)));
  complex<double> drhodr=dhldr*(one_+hl*(one_+hl*(half_+hl*sixth_)));
  return 2.*real(drhodr*conj(rho));
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double solveFcn4(double r,double theta,double val) {
  for (int i=0;i<20;i++) {
    double f=fcn4(r,theta,val);
    double fprime=fcn4deriv(r,theta);
    r-=f/fprime;
    r=max(numeric_limits<double>::epsilon(),r);
    if (abs(f)<=1.e-14*r*abs(fprime)) break;
  }
  return r;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void makeMainList(GUI_INPUT_PARAMETER_LIST_TYPE *&main_list) {
//TRACER_CALL(t,"makeMainList");
  main_list=OPERATOR_NEW GUI_INPUT_PARAMETER_LIST_TYPE("Main");

  { const char *group="Numerical Method Parameters";
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,ischeme,
      "scheme","forward Euler","modified Euler","Heun 3rd order",
      "classical 4th order") );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(nplot,
      "nplot",2,INT_MAX,group) );
  }

  { const char *group="Graphics";
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(winsize,
      "winsize",0.,1.,group));
  }
#ifdef DEBUG
//main_list->printOn(cout);
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void checkMainInput() {
//TRACER_CALL(t,"checkMainInput");
  scheme=static_cast<SCHEME>(ischeme);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void runMain(bool /*called_before*/) {
//TRACER_CALL(z,"runMain");
  double dtheta=2.*M_PI/static_cast<double>(nplot);
  complex<double> *point=OPERATOR_NEW_BRACKET(complex<double>,nplot+1);
  complex<double> *point1=OPERATOR_NEW_BRACKET(complex<double>,nplot+1);
  complex<double> *point9=OPERATOR_NEW_BRACKET(complex<double>,nplot+1);
  double theta=0.;
  switch (scheme) {
    case FIRST_ORDER: {
      for (int j=0;j<=nplot;j++,theta+=dtheta) {
        complex<double> z=polar(1.,theta);
        point[j]=z-one_;
        point9[j]=z*0.9-one_;
        point1[j]=z*1.1-one_;
      }
      break;
    }
    case SECOND_ORDER: {
      for (int j=0;j<=nplot;j++,theta+=dtheta) {
        double cs=cos(2.*theta);
        double r=sqrt(sqrt(3.+cs*cs)-cs);
        point[j]=polar(r,theta)-one_;
        r=sqrt(sqrt(2.24+cs*cs)-cs);
        point9[j]=polar(r,theta)-one_;
        r=sqrt(sqrt(3.84+cs*cs)-cs);
        point1[j]=polar(r,theta)-one_;
      }
      break;
    }
    case THIRD_ORDER: {
      double r=1.;
      for (int j=0;j<=nplot;j++,theta+=dtheta) {
        double rnew=solveFcn3(r,theta,1.);
        point[j]=polar(rnew,theta)-one_;
        rnew=solveFcn3(r,theta,0.9);
        point9[j]=polar(rnew,theta)-one_;
        rnew=solveFcn3(r,theta,1.1);
        point1[j]=polar(rnew,theta)-one_;
        r=rnew;
      }
      break;
    }
    case FOURTH_ORDER: {
      double r=1.;
      for (int j=0;j<=nplot;j++,theta+=dtheta) {
        double rnew=solveFcn4(r,theta,1.);
        point[j]=polar(rnew,theta)-one_;
        rnew=solveFcn4(r,theta,0.9);
        point9[j]=polar(rnew,theta)-one_;
        rnew=solveFcn4(r,theta,1.1);
        point1[j]=polar(rnew,theta)-one_;
        r=rnew;
      }
      break;
    }
  }
  double xmin=numeric_limits<double>::max();
  double xmax=-numeric_limits<double>::max();
  double ymin=numeric_limits<double>::max();
  double ymax=-numeric_limits<double>::max();
  theta=0.;
  for (int j=0;j<=nplot;j++,theta+=dtheta) {
    xmin=min(xmin,point[j].real());
    xmax=max(xmax,point[j].real());
    ymin=min(ymin,point[j].imag());
    ymax=max(ymax,point[j].imag());

    xmin=min(xmin,point9[j].real());
    xmax=max(xmax,point9[j].real());
    ymin=min(ymin,point9[j].imag());
    ymax=max(ymax,point9[j].imag());

    xmin=min(xmin,point1[j].real());
    xmax=max(xmax,point1[j].real());
    ymin=min(ymin,point1[j].imag());
    ymax=max(ymax,point1[j].imag());
  }
  Palette pal;
  XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
  XYGraphTool gt("absolute stability boundary","Real Part(lambda)",
    "Imaginary Part(lambda)",xmin,xmax,ymin,ymax,&cmap,0,winsize);
  gt.setbgColor("white");

  gt.newPage();
  gt.setfgColor("black");
  gt.drawAxes();
  gt.setfgColor("blue");
  gt.movePen(point[0].real(),point[0].imag());
  for (int j=1;j<=nplot;j++) {
    gt.drawLine(point[j].real(),point[j].imag());
  }
  gt.setfgColor("red");
  gt.movePen(point1[0].real(),point1[0].imag());
  for (int j=1;j<=nplot;j++) {
    gt.drawLine(point1[j].real(),point1[j].imag());
  }
  gt.setfgColor("green");
  gt.movePen(point9[0].real(),point9[0].imag());
  for (int j=1;j<=nplot;j++) {
    gt.drawLine(point9[j].real(),point9[j].imag());
  }
  gt.flush();
  XYGraphTool::WINDOW_TYPE::QuitButton qb;
  delete [] point;
  delete [] point9;
  delete [] point1;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void cleanup() {
//TRACER_CALL(t,"cleanup");
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
template int std::__cmath_power<int>(int,unsigned);
template class InputParameter<bool>;
