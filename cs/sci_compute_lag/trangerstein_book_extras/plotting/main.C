#include <fstream.h>
#include <iomanip.h>
#include <ostream.h>
#include <strstream.h>
#include <stdlib.h>
#include <Xm/ToggleB.h>
#include "setTraps.H"
#include "GUI.H"
#include "GUIVirtualInput.H"
#include "MemoryDebugger.H"
//#include "Palette.H"
#include "TimedObject.H"
#include "Tracer.H"
//#include "XColormap.H"
//#include "XYGraphTool.H"
#include <unistd.h>

#ifdef __GNUC__
//#pragma implementation "GUIInputParameter.C"
#endif
//#include "GUIInputParameter.C"

#define LENGTH_SUFFIX 6
#define LENGTH_NAME 80

BOOLEAN skip_gui=FALSE;
char* display_name=0;
char* sim_name=0;
double winsize=0.;
GUIInputParameterList *main_list=0;

extern void makeMainList(GUIInputParameterList*&);
extern void checkMainInput(BOOLEAN);
extern void runMain(BOOLEAN /*called_before*/);
extern void cleanup();
extern void shutdown();
extern void wait();
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void processCommandLine(int argc,char *argv[],char *display_name,
BOOLEAN &skip_gui,ifstream &in_file) {
//Tracer t("processCommandLine");
  if (argc<2) skip_gui=FALSE;
  else {
    in_file.open(argv[1],ios::in | ios::nocreate);
    if (in_file) {
      INTEGER i=2;
      while (i<argc) {
	if (strcmp(argv[i],"-d")==0) {
	  display_name=new char[strlen(argv[++i])+1];
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
#include "main.H"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int main(int argc,char* argv[]) {
#ifdef DEBUG
  setTraps();
#endif
  { MemoryDebugger md(1);
  {
    ifstream in_file;
    processCommandLine(argc,argv,display_name,skip_gui,in_file);
    makeMainList(main_list);
    if (in_file) {
      readMainInput(in_file,main_list,skip_gui);
      if (skip_gui) checkMainInput(FALSE);
    }

//  open graphical user interface; turn control of events to it
    if (skip_gui) {
      runMain(0);
      wait();
      cleanup();
      shutdown();
    } else {
      GUI gui(argv[0],display_name,main_list,&runMain,&checkMainInput,
	&cleanup,&shutdown);
      gui.createFileMenu();
      gui.createViewMenu();
      gui.createHelpMenu();
      gui.createMainWindow(argc,argv);
      gui.eventLoop();
    }

  }
  }
  return EXIT_SUCCESS;
}

#ifdef __GNUC__
extern "C" {
void MAIN__(void) {;}
}
#endif
