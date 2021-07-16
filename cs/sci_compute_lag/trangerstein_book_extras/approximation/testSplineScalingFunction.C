#include <iostream>
#include <math.h> // for HUGE_VAL,M_PI
#include <stdlib.h>

using namespace std;

#include "Debug.H"
#include "GUIEnumInputParameter.H"
#include "GUIInputs.H"
#include "MemoryDebugger.H"
#include "ScalingFunction.H"
#include "SetTraps.H"
#include "SplineFilterBank.H"
#include "TimedObject.H"
#include "Wavelet.H"
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

int spline_filter_bank_N=2; // for spline lowpass filter
int spline_filter_bank_L=4; // for alpha; beta will use M=(N+L)/2
int number_scales=2;

char* display_name=0;
bool skip_gui=false;
GUI_INPUT_PARAMETER_LIST_TYPE *param_list=0;
Palette *pal=0;
XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE *cmap=0;
XYGraphTool *gt=0;
XYGraphTool *gt2=0;
XYGraphTool *gt3=0;
XYGraphTool *gt4=0;

void setField(char* field,int field_width,int val) {
  field[field_width]='\0';
  int i=field_width-1;
  int n=val;
  for (;i>=0;i--,n/=10) field[i]='0'+n %10;
}

void checkMainInput() { }

void runMain(bool /*called_before*/) {
//TRACER_CALL(z,"runMain");
  if (pal==0) pal=OPERATOR_NEW Palette();
  if (cmap==0) {
    cmap=OPERATOR_NEW XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE(pal);
  }
  if (gt==0) {
    gt=OPERATOR_NEW XYGraphTool("scaling function","t","phi",
      0.,1.,0.,1.,cmap,0,0.5);
  }
  if (gt2==0) {
    gt2=OPERATOR_NEW XYGraphTool("dual scaling function","t","phi tilde",
      0.,1.,0.,1.,cmap,0,0.5);
  }
  if (gt3==0) {
    gt3=OPERATOR_NEW XYGraphTool("wavelet","t","psi",
      0.,1.,0.,1.,cmap,0,0.5);
  }
  if (gt4==0) {
    gt4=OPERATOR_NEW XYGraphTool("dual wavelet","t","psi tilde",
      0.,1.,0.,1.,cmap,0,0.5);
  }
  gt->setbgColor("white");
  gt2->setbgColor("white");
  gt3->setbgColor("white");
  gt4->setbgColor("white");

  SplineFilterBank sfb(spline_filter_bank_N,spline_filter_bank_L);
  ScalingFunction phi(sfb.lowpassFilter().impulseResponse());
  ScalingFunction phi_tilde(sfb.dualLowpassFilter().impulseResponse());
  Wavelet psi(sfb.highpassFilter().impulseResponse(),phi);
  Wavelet psi_tilde(sfb.dualHighpassFilter().impulseResponse(),phi_tilde);

  Signal *phi_values=phi.values(number_scales);
  double phi_min=DBL_MAX;
  double phi_max=-DBL_MAX;
  for (int n=phi_values->firstIndex();n<=phi_values->lastIndex();n++) {
    double pn=phi_values->value(n);
    phi_min=min(phi_min,pn);
    phi_max=max(phi_max,pn);
#ifdef DEBUG
//  cout << "\tphi_value[ " << n << " ] = " << pn << endl;
#endif
  }
  if (phi_max<=phi_min) {
    phi_max=max(0.,phi_max);
    phi_min=min(0.,phi_min);
    if (phi_max<=phi_min) {
      phi_max=1.;
      phi_min=0.;
    }
  }

  double tlo=phi.supportLow();
  double thi=phi.supportHigh();
  double dt=1./pow(2.,number_scales);
  gt->newPage();
  gt->rescale(tlo,thi,phi_min,phi_max);
  gt->setfgColor("black");
  gt->drawAxes();
  gt->setfgColor("blue");
  double t=tlo;
  gt->movePen(t,phi_values->value(phi_values->firstIndex()));
  for (int n=phi_values->firstIndex()+1;n<=phi_values->lastIndex();n++) {
    t+=dt;
    gt->drawLine(t,phi_values->value(n));
  }
  gt->flush();

  Signal *phi_tilde_values=phi_tilde.values(number_scales);
  double phi_tilde_min=DBL_MAX;
  double phi_tilde_max=-DBL_MAX;
  for (int n=phi_tilde_values->firstIndex();
  n<=phi_tilde_values->lastIndex();n++) {
    double pn=phi_tilde_values->value(n);
    phi_tilde_min=min(phi_tilde_min,pn);
    phi_tilde_max=max(phi_tilde_max,pn);
  }
  if (phi_tilde_max<=phi_tilde_min) {
    phi_tilde_max=max(0.,phi_tilde_max);
    phi_tilde_min=min(0.,phi_tilde_min);
    if (phi_tilde_max<=phi_tilde_min) {
      phi_tilde_max=1.;
      phi_tilde_min=0.;
    }
  }

  tlo=phi_tilde.supportLow();
  thi=phi_tilde.supportHigh();
  dt=1./pow(2.,number_scales);
  gt2->newPage();
  gt2->rescale(tlo,thi,phi_tilde_min,phi_tilde_max);
  gt2->setfgColor("black");
  gt2->drawAxes();
  gt2->setfgColor("blue");
  t=tlo;
  gt2->movePen(t,phi_tilde_values->value(phi_tilde_values->firstIndex()));
  for (int n=phi_tilde_values->firstIndex()+1;
  n<=phi_tilde_values->lastIndex();n++) {
    t+=dt;
    gt2->drawLine(t,phi_tilde_values->value(n));
  }
  gt2->flush();

  Signal *psi_values=psi.values(number_scales);
  double psi_min=DBL_MAX;
  double psi_max=-DBL_MAX;
  for (int n=psi_values->firstIndex();n<=psi_values->lastIndex();n++) {
    double pn=psi_values->value(n);
    psi_min=min(psi_min,pn);
    psi_max=max(psi_max,pn);
#ifdef DEBUG
//  cout << "\tpsi_value[ " << n << " ] = " << pn << endl;
#endif
  }
  if (psi_max<=psi_min) {
    psi_max=max(0.,psi_max);
    psi_min=min(0.,psi_min);
    if (psi_max<=psi_min) {
      psi_max=1.;
      psi_min=0.;
    }
  }

  tlo=psi.supportLow();
  thi=psi.supportHigh();
  dt=1./pow(2.,number_scales);
  gt3->newPage();
  gt3->rescale(tlo,thi,psi_min,psi_max);
  gt3->setfgColor("black");
  gt3->drawAxes();
  gt3->setfgColor("blue");
  t=tlo;
  gt3->movePen(t,psi_values->value(psi_values->firstIndex()));
  for (int n=psi_values->firstIndex()+1;n<=psi_values->lastIndex();n++) {
    t+=dt;
    gt3->drawLine(t,psi_values->value(n));
  }
  gt3->flush();
  delete psi_values; psi_values=0;

  Signal *psi_tilde_values=psi_tilde.values(number_scales);
  double psi_tilde_min=DBL_MAX;
  double psi_tilde_max=-DBL_MAX;
  for (int n=psi_tilde_values->firstIndex();
  n<=psi_tilde_values->lastIndex();n++) {
    double pn=psi_tilde_values->value(n);
    psi_tilde_min=min(psi_tilde_min,pn);
    psi_tilde_max=max(psi_tilde_max,pn);
  }
  if (psi_tilde_max<=psi_tilde_min) {
    psi_tilde_max=max(0.,psi_tilde_max);
    psi_tilde_min=min(0.,psi_tilde_min);
    if (psi_tilde_max<=psi_tilde_min) {
      psi_tilde_max=1.;
      psi_tilde_min=0.;
    }
  }

  tlo=psi_tilde.supportLow();
  thi=psi_tilde.supportHigh();
  dt=1./pow(2.,number_scales);
  gt4->newPage();
  gt4->rescale(tlo,thi,psi_tilde_min,psi_tilde_max);
  gt4->setfgColor("black");
  gt4->drawAxes();
  gt4->setfgColor("blue");
  t=tlo;
  gt4->movePen(t,psi_tilde_values->value(psi_tilde_values->firstIndex()));
  for (int n=psi_tilde_values->firstIndex()+1;
  n<=psi_tilde_values->lastIndex();n++) {
    t+=dt;
    gt4->drawLine(t,psi_tilde_values->value(n));
  }
  gt4->flush();
  delete psi_tilde_values; psi_tilde_values=0;

  XYGraphTool::WINDOW_TYPE::QuitButton qb;

  delete phi_values; phi_values=0;
  delete phi_tilde_values; phi_tilde_values=0;
}

void cleanup() {
  if (gt) delete gt; gt=0;
  if (gt2) delete gt2; gt2=0;
  if (gt3) delete gt3; gt3=0;
  if (gt4) delete gt4; gt4=0;
  if (cmap) delete cmap; cmap=0;
  if (pal) delete pal; pal=0;
}

void shutdown() {
  while (param_list->notEmpty()) delete param_list->delAfter(0);
  delete param_list; param_list=0;
}

int main(int argc,char* argv[]) {
#ifdef DEBUG
//setTraps();
#endif
  {
#ifdef MEM_DEBUG
    MemoryDebugger md(1);
#endif
#ifdef USE_GTK
    GTKWindow::gtkInit(argc,argv);
#endif

  param_list=OPERATOR_NEW GUI_INPUT_PARAMETER_LIST_TYPE("Param List");
  { const char *group="Problem Parameters";
    param_list->append(OPERATOR_NEW GUIInputParameter<int>(
      spline_filter_bank_N,"spline filter bank index N",1,INT_MAX,group));
    param_list->append(OPERATOR_NEW GUIInputParameter<int>(
      spline_filter_bank_L,"spline filter bank index L",1,INT_MAX,group));
    param_list->append(OPERATOR_NEW GUIInputParameter<int>(
      number_scales,"number of scales for scaling fcn",1,INT_MAX,group));
  }

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
              }
              else param_list->formattedRead(in_file,name);
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
      GUI gui(argv[0],display_name,param_list,&runMain,
        &checkMainInput,&cleanup,&shutdown,TRUE);
#endif
      gui.createFileMenu();
      gui.createViewMenu();
      gui.createHelpMenu();
      gui.createMainWindow(argc,argv);
      gui.eventLoop();
    }

  }
  return EXIT_SUCCESS;
}

template class InputParameter<bool>;
template class InputParameter<double>;
template class InputParameter<int>;
