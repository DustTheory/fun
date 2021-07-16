#include <cmath>
#include <iostream>
#include <limits>
#include <math.h>
#include <stdlib.h>

#include "BinomialCoefficient.H"
#include "Debug.H"
#include "GUIEnumInputParameter.H"
#include "GUIInputs.H"
#include "MemoryDebugger.H"
#include "SetTraps.H"
#include "SpecializedMatrix.H"
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

enum SCHEME{ADAMS_BASHFORTH,ADAMS_MOULTON,BACKWARD_DIFFERENTIATION_FORMULA};
SCHEME scheme=ADAMS_BASHFORTH;
int ischeme=static_cast<int>(scheme);

int nplot=128;
int order=1;
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
      "scheme","adams-bashforth","adams-moulton","bdf") );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(nplot,
      "nplot",2,INT_MAX,group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(order,
      "order",1,10,group) );
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
  cout << "order = " << order << endl;
  double dtheta=2.*M_PI/static_cast<double>(nplot);
  complex<double> *point=OPERATOR_NEW_BRACKET(complex<double>,nplot+1);
  bool *ok=OPERATOR_NEW_BRACKET(bool,nplot+1);
  double eps2=sqrt(numeric_limits<double>::epsilon());
  switch (scheme) {
    case ADAMS_BASHFORTH: {
      double *gamma=OPERATOR_NEW_BRACKET(double,order);
      for (int i=0;i<order;i++) {
        gamma[i]=1.;
        for (int j=0;j<i;j++) gamma[i]-=gamma[j]/static_cast<double>(i-j+1);
      }
//    cout << "gamma =";
//    for (int i=0;i<order;i++) cout << " " << gamma[i];
//    cout << endl;
      BinomialCoefficient bc(order-1);
      double *companion=OPERATOR_NEW_BRACKET(double,order);
      double si=1.;
      for (int i=0;i<order;i++,si=-si) {
        double sum=0.;
        for (int j=i;j<order;j++) sum+=gamma[j]*bc.value(j,i);
        companion[order-1-i]=sum*si;
      }
//    cout << "companion =";
//    for (int i=0;i<order;i++) cout << " " << companion[i];
//    cout << endl;
      double theta=0.;
      for (int j=0;j<=nplot;j++,theta+=dtheta) {
        complex<double> z(cos(theta),sin(theta));
        complex<double> denom=gamma[0];
        for (int k=1;k<order;k++) {
          denom+=gamma[k]*pow(1.-1./z,k);
        }
        point[j]=(z-1.)/denom;
        Vector<double,complex<double> > p(order+1);
        p[0]=1.;
        for (int i=0;i<order;i++) p[i+1]=-companion[order-1-i]*point[j];
        p[1]-=1.;
//      cout << "p = " << endl;
//      p.printOn(cout);
        CompanionMatrix<double,complex<double> > C(p);
//      cout << "C = " << endl;
//      C.printOn(cout);
        SquareMatrix<double,complex<double> > *V=0;
        SquareMatrix<double,complex<double> > *U=0;
        Vector<double,complex<double> > *eig=C.eigenvalues(V,U);
        double eigmax=0.;
        for (int i=0;i<eig->size();i++) eigmax=max(eigmax,abs((*eig)[i]));
        ok[j]=(eigmax<=1.+eps2);
//      cout << "\tpoint,ok,eigs = " << point[j] << " " << ok[j];
//      for (int i=0;i<order;i++) cout << " " << abs((*eig)[i]);
//      cout << endl;
        delete eig; eig=0;
      }
      delete [] companion; companion=0;
      delete [] gamma; gamma=0;
      break;
    }
    case ADAMS_MOULTON: {
      double *gamma=OPERATOR_NEW_BRACKET(double,order);
      gamma[0]=1.;
      for (int i=1;i<order;i++) {
        gamma[i]=0.;
        for (int j=0;j<i;j++) gamma[i]-=gamma[j]/static_cast<double>(i-j+1);
      }
//    cout << "gamma =";
//    for (int i=0;i<order;i++) cout << " " << gamma[i];
//    cout << endl;
      BinomialCoefficient bc(order-1);
      double *companion=OPERATOR_NEW_BRACKET(double,order);
      double si=1.;
      for (int i=0;i<order;i++,si=-si) {
        double sum=0.;
        for (int j=i;j<order;j++) sum+=gamma[j]*bc.value(j,i);
        companion[order-1-i]=sum*si;
      }
//    cout << "companion =";
//    for (int i=0;i<order;i++) cout << " " << companion[i];
//    cout << endl;
      double theta=0.;
      bool print_now=true;
      for (int j=0;j<=nplot;j++,theta+=dtheta) {
        complex<double> z(cos(theta),sin(theta));
        complex<double> denom=gamma[0];
        for (int k=1;k<order;k++) {
          denom+=gamma[k]*pow(1.-1./z,k);
        }
        point[j]=(1.-1./z)/denom;
        if (order==1) {
          Vector<double,complex<double> > p(2);
          p[0]=1.-companion[order-1]*point[j];
          p[1]=-1.;
//        cout << "p = " << endl;
//        p.printOn(cout);
          CompanionMatrix<double,complex<double> > C(p);
//        cout << "C = " << endl;
//        C.printOn(cout);
          SquareMatrix<double,complex<double> > *V=0;
          SquareMatrix<double,complex<double> > *U=0;
          Vector<double,complex<double> > *eig=C.eigenvalues(V,U);
          double eigmax=0.;
          for (int i=0;i<eig->size();i++) eigmax=max(eigmax,abs((*eig)[i]));
          ok[j]=(eigmax<=1.+eps2);
//        cout << "\tpoint,ok,eigs = " << point[j] << " " << ok[j];
//        for (int i=0;i<order;i++) cout << " " << abs((*eig)[i]);
//        cout << endl;
          delete eig; eig=0;
        } else {
          Vector<double,complex<double> > p(order);
          p[0]=1.-companion[order-1]*point[j];
          for (int i=1;i<order;i++) p[i]=-companion[order-1-i]*point[j];
          p[1]-=1.;
//        cout << "p = " << endl;
//        p.printOn(cout);
          CompanionMatrix<double,complex<double> > C(p);
//        cout << "C = " << endl;
//        C.printOn(cout);
          SquareMatrix<double,complex<double> > *V=0;
          SquareMatrix<double,complex<double> > *U=0;
          Vector<double,complex<double> > *eig=C.eigenvalues(V,U);
          double eigmax=0.;
          for (int i=0;i<eig->size();i++) eigmax=max(eigmax,abs((*eig)[i]));
          ok[j]=(eigmax<=1.+eps2);
//        cout << "\tpoint,ok,eigs = " << point[j] << " " << ok[j];
//        for (int i=0;i<eig->size();i++) cout << " " << abs((*eig)[i]);
//        cout << endl;
          delete eig; eig=0;
        }
      }
      delete [] companion; companion=0;
      delete [] gamma; gamma=0;
      break;
    }
    case BACKWARD_DIFFERENTIATION_FORMULA: {
      BinomialCoefficient bc(order);
      double *companion=OPERATOR_NEW_BRACKET(double,order+1);
      double si=1.;
      for (int i=0;i<=order;i++,si=-si) {
        double sum=0.;
        for (int j=max(1,i);j<=order;j++) {
          sum+=bc.value(j,i)/static_cast<double>(j);
        }
        companion[i]=sum*si;
      }
//    cout << "companion =";
//    for (int i=0;i<=order;i++) cout << " " << companion[i];
//    cout << endl;
      double theta=0.;
      for (int j=0;j<=nplot;j++,theta+=dtheta) {
        complex<double> z(cos(theta),sin(theta));
        point[j]=0.;
        for (int k=1;k<=order;k++) {
          point[j]+=pow(1.-1./z,k)/static_cast<double>(k);
        }
        Vector<double,complex<double> > p(order+1);
        p[0]=companion[0]-point[j];
        for (int i=1;i<=order;i++) p[i]=companion[i];
//      cout << "p = " << endl;
//      p.printOn(cout);
        CompanionMatrix<double,complex<double> > C(p);
//      cout << "C = " << endl;
//      C.printOn(cout);
        SquareMatrix<double,complex<double> > *V=0;
        SquareMatrix<double,complex<double> > *U=0;
        Vector<double,complex<double> > *eig=C.eigenvalues(V,U);
        double eigmax=0.;
        for (int i=0;i<eig->size();i++) eigmax=max(eigmax,abs((*eig)[i]));
        ok[j]=(eigmax<=1.+eps2);
//      cout << "\tpoint,ok,eigs = " << point[j] << " " << ok[j];
//      for (int i=0;i<eig->size();i++) cout << " " << abs((*eig)[i]);
//      cout << endl;
        delete eig; eig=0;
      }
      delete [] companion; companion=0;
      break;
    }
  }
  double xmin=numeric_limits<double>::max();
  double xmax=-numeric_limits<double>::max();
  double ymin=numeric_limits<double>::max();
  double ymax=-numeric_limits<double>::max();
  double theta=0.;
  for (int j=0;j<=nplot;j++,theta+=dtheta) {
    if (ok[j]) {
      xmin=min(xmin,point[j].real());
      xmax=max(xmax,point[j].real());
      ymin=min(ymin,point[j].imag());
      ymax=max(ymax,point[j].imag());
    }
  }
  Palette pal;
  XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
  XYGraphTool gt("absolute stability boundary","Real Part(lambda)",
    "Imaginary Part(lambda)",0.,1.,0.,1.,&cmap,0,winsize);
  gt.setbgColor("white");

  gt.newPage();
  gt.rescale(xmin,xmax,ymin,ymax);
  gt.setfgColor("black");
  gt.drawAxes();
  gt.setfgColor("red");
  gt.movePen(point[0].real(),point[0].imag());
  bool last_ok=true;
  for (int j=1;j<=nplot;j++) {
    if (last_ok) gt.drawLine(point[j].real(),point[j].imag());
    else gt.movePen(point[j].real(),point[j].imag());
    last_ok=ok[j];
  }
  gt.flush();
  XYGraphTool::WINDOW_TYPE::QuitButton qb;
  delete [] ok; ok=0;
  delete [] point; point=0;
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
template class InputParameter<bool>;
