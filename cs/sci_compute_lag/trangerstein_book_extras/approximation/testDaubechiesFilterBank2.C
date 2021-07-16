#include <iostream>
#include <stdlib.h>
#include "BinomialCoefficient.H"
#include "DaubechiesFilterBank.H"
#include "GUIEnumInputParameter.H"
#include "GUIInputs.H"
#include "MemoryDebugger.H"
#include "NumPtr.C"
#include "ScalingFunction.H"
#include "SetTraps.H"
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

int nlevels=1;
int scaling_function_order=1;
int quadrature_scale=1;
int nerr=1;

double a=0.;
double b=1.; // approximate f(x) on [a,b]
double f(double t) { return exp(t); }
//double f(double x) { return x; }
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
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(
      scaling_function_order,"scaling_function_order",1,INT_MAX,group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(nlevels,
      "approximation scale",1,INT_MAX,group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(quadrature_scale,
      "quadrature scale",1,INT_MAX,group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(nerr,
      "error sampling scale",1,INT_MAX,group) );
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
  if (quadrature_scale<floor(scaling_function_order*0.5+0.5)) {
    cerr << "quadrature scale too small for scaling function order"
         << endl;
    quadrature_scale=
      static_cast<int>(floor(scaling_function_order*0.5+0.5));
  }
  if (nerr<nlevels) {
    cerr << "error sampling scale smaller than approximation scale"
         << endl;
    nerr=nlevels;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void runMain(bool /*called_before*/) {
//TRACER_CALL(t,"runMain");
  DaubechiesFilterBank dfb(scaling_function_order);
  const Filter &lowpass=dfb.lowpassFilter();
  const Filter &dual_lowpass=dfb.dualLowpassFilter();
  ScalingFunction phi(lowpass.impulseResponse());
  ScalingFunction phi_tilde(dual_lowpass.impulseResponse());

  int twotoqs=pow(2,quadrature_scale);
  NumPtr<double> errors(nlevels+1);
  double errors_max=-HUGE_VAL;
  double errors_min=HUGE_VAL;

  Palette pal;
  XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
  XYGraphTool gt("functions","t","function",a,b,0.,1.,&cmap,0,0.5);

  double log2=log(2.);
  for (int j=0;j<=nlevels;j++) {
    int twotonlevels=pow(2,j);
    int twotonlpqs=twotonlevels*twotoqs;

//  f(t) approximated in [a,b] by sum_n gamma_n phi(2^N t - n)
//    if support(phi) = [alpha,beta]
//    then 2^N a - beta <= n <= 2^N b - alpha
//
//  gamma_n = integral f( 2^{-N} [ tau + n ] ) phi_tilde( tau ) d tau
//    if support(phi_tilde) = [alpha_tilde,beta_tilde]
//    then
//     a-2^{-N}(beta-alpha_tilde)<=2^{-N}[tau+n]<=b+2^{-N}(beta_tilde-alpha)
//   to compute integrals, need values of f
//     at points tau_q = 2^{-Q} q between these bounds
//     ==> 2^{N+Q} [ a - 2^{-N}(beta-alpha_tilde) ]
//         <= 2^Q[tau+n]
//         <= 2^{N+Q} [ b + 2^{-N}(beta_tilde-alpha) ]
//  floor(x+.5) <==> round x to nearest integer
    int fv_lo=static_cast<int>(floor(twotonlpqs*a
      -twotoqs*(phi.supportHigh()-phi_tilde.supportLow())+0.5));
    int fv_hi=static_cast<int>(floor(twotonlpqs*b
      +twotoqs*(phi_tilde.supportHigh()-phi.supportLow())+0.5));
    Signal fvalues(fv_lo,fv_hi);
    double h=1./static_cast<double>(twotonlpqs);
    double x=static_cast<double>(fv_lo)*h;
    double ymax=-HUGE_VAL;
    double ymin=HUGE_VAL;
    for (int i=fv_lo;i<=fv_hi;i++,x+=h) {
      fvalues.value(i)=f(x); 
      if (x>=a && x<=b) {
        ymax=max(ymax,fvalues.value(i));
        ymin=min(ymin,fvalues.value(i));
      }
    }
  
//  gamma_n = integral_alpha^beta f(2^{-N}[tau+n]) phi_tilde(tau) d tau
//          approximated by sum_q f(2^{-N}[tau_q+n]) phi_tilde(tau_q) w_q
//    where tau_q = 2^{-Q} q <==> 2^Q alpha <= q <= 2^Q beta
    Signal *phi_tilde_values=phi_tilde.values(quadrature_scale);
    int phi_tildev_lo=phi_tilde_values->firstIndex();
    int phi_tildev_hi=phi_tilde_values->lastIndex();

//  f(t) approximated by sum_n gamma_n phi_tilde(2^N t - n)
//    if support(phi) = [alpha,beta]
//    then alpha <= 2^N t - n <= beta
//    ==> 2^N a - beta <= n <= 2^N b - alpha
    int gamma_lo=
      static_cast<int>(floor(twotonlevels*a-phi.supportHigh()+0.5));
    int gamma_hi=
      static_cast<int>(floor(twotonlevels*b-phi.supportLow()+0.5));
    Signal gamma(gamma_lo,gamma_hi);
    NumPtr<double> table(quadrature_scale+1);//quad extrapolation table
    h=1./static_cast<double>(twotonlevels);
    for (int n=gamma_lo;n<=gamma_hi;n++) {
//    to compute integral for gamma_n,
//      sample phi_tilde on scales 2^{-j}, 0 <= j <= Q
//      and extrapolate trapezoidal rule Q times
//    if we only computed trapezoidal rule with h=support(phi_tilde)
//      the quadrature would produce value zero for all except Haar wavelet
      double stride=twotoqs;
      double hl=1.;
      double trapezoidal=0.;
      int foffset=twotoqs*n;
      for (int i=phi_tildev_lo+stride;i<=phi_tildev_hi;i+=stride) {
//      fvalues computed on scale 2^{-N-Q}
//      phi_tilde_values computed on scale 2^{-Q}
//        ==> alpha_tilde <= 2^{-Q} i <= beta_tilde
        trapezoidal+=fvalues.value(i+foffset)
                    *phi_tilde_values->value(i);
      }
      if (scaling_function_order==1) {
//      only case where phi_tilde is nonzero at an endpoint
        trapezoidal+=0.5*(fvalues.value(phi_tildev_lo+foffset)
                         +fvalues.value((phi_tildev_hi+1)+foffset))
                        *phi_tilde_values->value(phi_tildev_lo);
      }
      trapezoidal*=hl;
//    extrapolation
      for (int j=1;j<=quadrature_scale;j++) {
        double old_table=trapezoidal;
        int stride2=stride/2;
        double midpoint=0.;
        for (int i=phi_tildev_lo+stride2;i<=phi_tildev_hi;i+=stride) {
          midpoint+=fvalues.value(i+foffset)*phi_tilde_values->value(i);
        }
        midpoint*=hl;
        trapezoidal=(trapezoidal+midpoint)*0.5;
        table[0]=trapezoidal;
        int fourtok=4;
        for (int k=1;k<=j;k++) {
          double extrapolant=table[k-1]
            +(table[k-1]-old_table)/static_cast<double>(fourtok-1);
          old_table=table[k];
          table[k]=extrapolant;
          fourtok*=4;
        }
        hl*=0.5;
        stride=stride2;
      }
      gamma.value(n)=table[quadrature_scale];
    }

    int jerr=j+nerr-nlevels;
    int twotonerr=pow(2,jerr);
    h=1./static_cast<double>(twotonerr);

    gt.rescale(a,b,ymin,ymax);
    gt.newPage();
    gt.setfgColor("black");
    gt.drawAxes();
    gt.setfgColor("red");
    double t=a;
    gt.movePen(t,f(t));
    t+=h;
    for (int i=static_cast<int>(floor(twotonerr*a+0.5)+1);
    i<=static_cast<int>(floor(twotonerr*b+0.5));i++,t+=h) {
      gt.drawLine(t,f(t));
    }
    gt.setfgColor("blue");

//  err_m = f(2^{-M} i) - sum_n gamma_n phi(2^N 2^{-M} i - n)
//    where a <= 2^{-M} i <= b <==> 2^M a <= i <= 2^M b
//  and
//    alpha <= 2^{N-M} i - n <= beta
//      <==> 2^N a - beta <= n <= 2^N b - alpha
    Signal *phi_values=phi.values(nerr-nlevels);
    int phiv_lo=phi_values->firstIndex();
    int phiv_hi=phi_values->lastIndex();
    int offset=static_cast<int>(a*static_cast<double>(twotonlevels));
    int stride=twotonerr/twotonlevels;
    double errmax=-HUGE_VAL;
    for (int i=static_cast<int>(floor(twotonerr*a+0.5));
    i<=static_cast<int>(floor(twotonerr*b+0.5));i++) {
      double error=0.;
      for (int n=
      static_cast<int>(floor(a*twotonlevels-phi_tilde.supportHigh()+0.5));
      n<=static_cast<int>(floor(b*twotonlevels-phi_tilde.supportLow()+0.5));
      n++) {
        int m=i-n*stride;
        if (m>=phiv_lo && m <= phiv_hi) {
          error+=gamma.value(n)*phi_values->value(m);
        }
      }
      gt.drawPlus(i*h,error,h*0.5);
      error=abs(f(i*h)-error);
      errmax=max(errmax,error);
    }
#ifdef DEBUG
    cout << "j,errmax = " << j << " " << errmax << endl;
#endif
    errors[j]=log(errmax)/log2;
    errors_min=min(errors_min,errors[j]);
    errors_max=max(errors_max,errors[j]);
    gt.flush();
    XYGraphTool::WINDOW_TYPE::QuitButton qb;
    delete phi_tilde_values; phi_tilde_values=0;
    delete phi_values; phi_values=0;
  }
  XYGraphTool gte("log_2 error vs scale","scale","error",0,nlevels,
    errors_min,errors_max,&cmap,0,0.5);
  gte.newPage();
  gte.setfgColor("black");
  gte.drawAxes();
  gte.setfgColor("blue");
  gt.movePen(0,errors[0]);
  for (int j=1;j<=nlevels;j++) {
    gte.drawLine(j,errors[j]);
  }
  gt.flush();
  XYGraphTool::WINDOW_TYPE::QuitButton qb;
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
