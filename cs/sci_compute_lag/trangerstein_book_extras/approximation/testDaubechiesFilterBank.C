#include <iostream>
#include <stdlib.h>
#include "BinomialCoefficient.H"
#include "DaubechiesFilterBank.H"
#include "GUIEnumInputParameter.H"
#include "GUIInputs.H"
#include "MemoryDebugger.H"
#include "NumPtr.C"
#include "SetTraps.H"
#include "Tracer.H"
//#include "XColormap.H"
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

int nlevels=1;
int scaling_function_order=1;

double fcn(double t) { return exp(t); }
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

  { const char *group="Filter Bank Parameters";
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(nlevels,
      "nlevels",1,INT_MAX,group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(
      scaling_function_order,"scaling_function_order",1,INT_MAX,group) );
  }

  { const char *group="Graphics";
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(winsize,
      "winsize",0.,1.,group));
//  main_list->append(OPERATOR_NEW GUIInputParameter<int>(
//    nplot,"nplot",2,INT_MAX,group) );
  }
#ifdef DEBUG
//main_list->printOn(cout);
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void checkMainInput() {
//TRACER_CALL(t,"checkMainInput");
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void runMain(bool /*called_before*/) {
//TRACER_CALL(t,"runMain");
  int twotonlevels=pow(2,nlevels);
  Signal s(0,twotonlevels-1);
  double h=1./static_cast<double>(twotonlevels);
  double x=0;
  double ymax=-HUGE_VAL;
  double ymin=HUGE_VAL;
  for (int n=0;n<twotonlevels;n++,x+=h) {
    s.value(n)=fcn(x);
    ymax=max(ymax,s.value(n));
    ymin=min(ymin,s.value(n));
  }

  Palette pal;
  XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
  XYGraphTool gtlow("low signal","n","signal",s.firstIndex(),s.lastIndex(),
    ymin,ymax,&cmap,0,0.5);
  gtlow.newPage();
  gtlow.setfgColor("black");
  gtlow.drawAxes();
  gtlow.setfgColor("blue");
  for (int n=s.firstIndex();n<=s.lastIndex();n++) {
    gtlow.drawPlus(n,s.value(n),0.5);
  }
  gtlow.flush();
  {
  XYGraphTool::WINDOW_TYPE::QuitButton qb;
  }

  XYGraphTool gthigh("high signal","n","signal",0.,1.,0.,1.,&cmap,0,0.5);
  
  DaubechiesFilterBank sfb(scaling_function_order);
  int af=sfb.analysisFilterStart();
  int al=sfb.analysisFilterFinish();
  int sf=sfb.synthesisFilterStart();
  int sl=sfb.synthesisFilterFinish();

  NumPtr<Signal*> low_signal(nlevels+1);
  NumPtr<Signal*> high_signal(nlevels+1);
  low_signal[nlevels]=&s;
  for (int j=nlevels-1;j>=0;j--) {
#ifdef DEBUG
    cout << "analysis on level " << j << endl;
#endif
    int f=low_signal[j+1]->firstIndex();
    int l=low_signal[j+1]->lastIndex();
    int sjf=(f-al)/2;
    int sjl=(l-af)/2;
    low_signal[j]=OPERATOR_NEW Signal(sjf,sjl);
    high_signal[j]=OPERATOR_NEW Signal(sjf,sjl);
    low_signal[j]->operator=(0.);
    high_signal[j]->operator=(0.);
    sfb.analyze(*low_signal[j+1],*low_signal[j],*high_signal[j]);

    ymax=-HUGE_VAL;
    ymin=HUGE_VAL;
    for (int n=low_signal[j]->firstIndex();n<=low_signal[j]->lastIndex();
    n++) {
      ymax=max(ymax,low_signal[j]->value(n));
      ymin=min(ymin,low_signal[j]->value(n));
    }
    gtlow.rescale(low_signal[j]->firstIndex(),low_signal[j]->lastIndex(),
      ymin,ymax);
    gtlow.newPage();
    gtlow.setfgColor("black");
    gtlow.drawAxes();
    gtlow.setfgColor("blue");
    for (int n=low_signal[j]->firstIndex();n<=low_signal[j]->lastIndex();
    n++) {
      gtlow.drawPlus(n,low_signal[j]->value(n),0.5);
    }
    gtlow.flush();

    ymax=-HUGE_VAL;
    ymin=HUGE_VAL;
    for (int n=high_signal[j]->firstIndex();n<=high_signal[j]->lastIndex();
    n++) {
      ymax=max(ymax,high_signal[j]->value(n));
      ymin=min(ymin,high_signal[j]->value(n));
    }
    gthigh.rescale(high_signal[j]->firstIndex(),
      high_signal[j]->lastIndex(),ymin,ymax);
    gthigh.newPage();
    gthigh.setfgColor("black");
    gthigh.drawAxes();
    gthigh.setfgColor("blue");
    for (int n=high_signal[j]->firstIndex();n<=high_signal[j]->lastIndex();
    n++) {
      gthigh.drawPlus(n,high_signal[j]->value(n),0.5);
    }
    gthigh.flush();
    XYGraphTool::WINDOW_TYPE::QuitButton qb;
  }
  NumPtr<Signal*> output(nlevels+1);
  for (int j=1;j<=nlevels;j++) {
#ifdef DEBUG
    cout << "synthesis on level " << j << endl;
#endif
    int f=min(low_signal[j-1]->firstIndex(),high_signal[j-1]->firstIndex());
    int l=max(low_signal[j-1]->lastIndex(),high_signal[j-1]->lastIndex());
    output[j]=OPERATOR_NEW Signal(sf+2*f,sl+2*l);
    sfb.synthesize(*low_signal[j-1],*high_signal[j-1],*output[j]);

    ymax=-HUGE_VAL;
    ymin=HUGE_VAL;
    for (int n=output[j]->firstIndex();n<=output[j]->lastIndex();
    n++) {
      ymax=max(ymax,output[j]->value(n));
      ymin=min(ymin,output[j]->value(n));
    }
    gtlow.rescale(output[j]->firstIndex(),output[j]->lastIndex(),
      ymin,ymax);
    gtlow.newPage();
    gtlow.setfgColor("black");
    gtlow.drawAxes();
    gtlow.setfgColor("blue");
    for (int n=output[j]->firstIndex();n<=output[j]->lastIndex();
    n++) {
      gtlow.drawPlus(n,output[j]->value(n),0.5);
    }
    gtlow.flush();
    XYGraphTool::WINDOW_TYPE::QuitButton qb;
  }
  double error=0.;
  for (int n=s.firstIndex();n<=s.lastIndex();n++) {
    error=max(error,abs(s.value(n)-output[nlevels]->value(n)));
  }
  cout << "\tmax error in final synthesis = " << error << endl;
  for (int j=nlevels;j>=1;j--) {
    delete output[j]; output[j]=0;
  }
  for (int j=0;j<nlevels;j++) {
    delete low_signal[j]; low_signal[j]=0;
    delete high_signal[j]; high_signal[j]=0;
  }
  low_signal[nlevels]=0;
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
INSTANTIATE_NUMPTR(Signal)
