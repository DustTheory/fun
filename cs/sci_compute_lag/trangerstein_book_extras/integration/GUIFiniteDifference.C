#include <cmath>
#include <iostream>
#include <math.h> // for HUGE_VAL,M_PI
#include <rfftw.h>
#include <stdlib.h>

#include "Debug.H"
#include "GUIEnumInputParameter.H"
#include "GUIInputs.H"
#include "MemoryDebugger.H"
#include "SetTraps.H"
#include "TimedObject.H"
#include "Tracer.H"
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

enum SCHEME{ONE_SIDED_DIFS,CENTERED_DIFS,RICHARDSON_ONE_SIDED,
  RICHARDSON_CENTERED};
SCHEME scheme=ONE_SIDED_DIFS;
int ischeme=static_cast<int>(scheme);

int nmax=2;
double x=3.;

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
      "scheme","one sided","centered","richardson one sided",
      "richardson centered") );
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(x,
      "x",0.,DBL_MAX,group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(nmax,
      "max power of 2",2,48,group) );
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

  double *log_error=OPERATOR_NEW_BRACKET(double,nmax);
  double *log_h=OPERATOR_NEW_BRACKET(double,nmax);

  double log10=log(10.);
  double f=log(x);
  double fp=1./x;
  double h=1.;
  switch (scheme) {
    case ONE_SIDED_DIFS: {
      for (int i=0;i<nmax;i++,h*=0.5) {
        double error=abs(fp-(log(x+h)-f)/h);
        if (error<=0.) error=abs(f)*DBL_EPSILON;
#ifdef DEBUG
//      cout << "\terror(" << h << ") = " << error << endl;
#endif
        log_error[i]=log(error)/log10;
        log_h[i]=-log(h)/log10;
      }
      break;
    }
    case CENTERED_DIFS: {
      for (int i=0;i<nmax;i++,h*=0.5) {
        double error=abs(fp-(log(x+h)-log(x-h))/(2.*h));
        if (error<=0.) error=abs(f)*DBL_EPSILON;
        log_error[i]=log(error)/log10;
        log_h[i]=-log(h)/log10;
      }
      break;
    }
    case RICHARDSON_ONE_SIDED: {
//    TRACER_CALL(z,"runMain RICHARDSON_ONE_SIDED");
      double *table=OPERATOR_NEW_BRACKET(double,nmax);
      for (int i=1;i<=nmax;i++,h*=0.5) {
        table[nmax-i]=(log(x+h)-f)/h;
#ifdef DEBUG
//      cout << "\n\ti = " << i << endl;
//      cout << "\ttable[" << nmax-i << "] = " << table[nmax-i] << endl;
#endif
        double p=2.;
        for (int j=nmax+1-i;j<nmax;j++,p*=2.) {
          table[j]=table[j-1]+(table[j-1]-table[j])/(p-1.);
#ifdef DEBUG
//        cout << "\ttable[" << j << "] = " << table[j] << endl;
#endif
        }
        double error=abs(fp-table[nmax-1]);
        if (error<=0.) error=abs(f)*DBL_EPSILON;
        log_error[i-1]=log(error)/log10;
        log_h[i-1]=-log(h)/log10;
      }
      delete [] table;
      break;
    }
    case RICHARDSON_CENTERED: {
//    TRACER_CALL(z,"runMain RICHARDSON_CENTERED");
      double *table=OPERATOR_NEW_BRACKET(double,nmax);
      for (int i=1;i<=nmax;i++,h*=0.5) {
        table[nmax-i]=(log(x+h)-log(x-h))/(2.*h);
#ifdef DEBUG
//      cout << "\n\ti = " << i << endl;
//      cout << "\ttable[" << nmax-i << "] = " << table[nmax-i] << endl;
#endif
        double p=4.;
        for (int j=nmax+1-i;j<nmax;j++,p*=4.) {
          table[j]=table[j-1]+(table[j-1]-table[j])/(p-1.);
#ifdef DEBUG
//        cout << "\ttable[" << j << "] = " << table[j] << endl;
#endif
        }
        double error=abs(fp-table[nmax-1]);
        if (error<=0.) error=abs(f)*DBL_EPSILON;
        log_error[i-1]=log(error)/log10;
        log_h[i-1]=-log(h)/log10;
      }
      delete [] table;
      break;
    }
    default:
      abort();
  }
  double xmin=HUGE_VAL;
  double xmax=-HUGE_VAL;
  double ymin=HUGE_VAL;
  double ymax=-HUGE_VAL;
  for (int i=0;i<nmax;i++) {
    ymin=min(ymin,log_error[i]);
    ymax=max(ymax,log_error[i]);
    xmin=min(xmin,log_h[i]);
    xmax=max(xmax,log_h[i]);
  }

  Palette pal;
  XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
  XYGraphTool gt("log_10(error) vs -log_10(h)","log_10(h)","log_10(error)",
    xmin,xmax,ymin,ymax,&cmap,0,winsize);
  gt.setbgColor("white");
  gt.setfgColor("black");
  gt.drawAxes();
  gt.setfgColor("blue");
  for (int i=0;i<nmax;i++) {
    if (i==0) gt.movePen(log_h[i],log_error[i]);
    else gt.drawLine(log_h[i],log_error[i]);
  }
  gt.setfgColor("green");
  double dx=(xmax-xmin)/static_cast<double>(2*nmax);
  for (int i=0;i<nmax;i++) {
    gt.drawPlus(log_h[i],log_error[i],dx);
  }
  XYGraphTool::WINDOW_TYPE::QuitButton qb;
  delete [] log_error;
  delete [] log_h;
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
