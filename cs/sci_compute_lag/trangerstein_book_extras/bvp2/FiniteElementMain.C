#include <float.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include "Arch.H"
#include "Debug.H"
#include "SetTraps.H"
#include "InputParameter.H"
#include "Palette.H"
#include "TimedObject.H"
#include "Tracer.H"
#include "VGT.H"
#include "XColormap.H"
#include "XYGraphTool.H"
#include <unistd.h>

using namespace std;

//#define LENGTH_NAME 80
#define NPTS 1000

extern "C" {
  double F77NAME(approximation)(const int &continuity,
    const int &nelements,const int &nnodes,const int &order,
    const double &x,
    const int *element_to_node,const double *mesh,const double* soln);
  void F77NAME(canonical)(const int &continuity,const int &order,
    double *basis,double *basis_deriv,double *gauss,double *weight);
  void F77NAME(grid)(const int &continuity,const int &nelements,
    const int &nnodes,const int &order,
    int *element_to_node,double *mesh);
  void F77NAME(initialize)(const int &continuity,const int &nelements,
    const int &nnodes,const int &order,
    const double *basis,const double *basis_deriv,
    const int *element_to_node,const double *gauss,const double *mesh,
    const double *weight,
    bool *dirichlet,double *pgauss,double *rgauss,double *residual,
    double *soln);
  double F77NAME(solution)(const double &x);
  void F77NAME(precg)(const int &limit,const int &ndigit,
    const int &nnelements,const int &nnodes,const int &order,
    const double *basis,const double *basis_deriv,const bool *dirichlet,
    const int *element_to_node,const double *pgauss,
    const double *rgauss,
    double *b,double *x,
    double *w);
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

int main(int argc,char* argv[]) {
  ASSERT(argc>1);
#ifdef DEBUG
  setTraps();
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
    F77NAME(machine).undefind=HUGE_VAL;
    F77NAME(machine).pi=M_PI;

    int continuity=0;
    int limit=1;
    int nelements=1;
    int ndigit=1;
    int order=1;

//  define input parameters and bounds
    InputParameterList param_list("Main List");
    { char *group="Numerical Method Parameters";
      param_list.append(OPERATOR_NEW InputParameter<int>(continuity,
        "continuity",0,2,group) );
      param_list.append(OPERATOR_NEW InputParameter<int>(limit,
        "limit",1,INT_MAX,group) );
      param_list.append(OPERATOR_NEW InputParameter<int>(ndigit,
        "ndigit",1,DBL_DIG,group) );
      param_list.append(OPERATOR_NEW InputParameter<int>(nelements,
        "nelements",1,INT_MAX,group) );
      param_list.append(OPERATOR_NEW InputParameter<int>(order,
        "order",1,3,group) );
    }

//  read input parameters 
    ifstream in_file;
    in_file.open(argv[1],ios::in);
    ASSERT(in_file);
    istream &infile(in_file);
    infile.clear(ios::goodbit);
    infile.seekg(0,ios::beg);
    char name[LENGTH_NAME], comment[LENGTH_COMMENT];
    while ( in_file >> setw(LENGTH_NAME) >> name ) {
      if ( strcmp(name,"Main") == 0 ) {
        in_file.getline( comment, LENGTH_COMMENT);
        while ( in_file >> setw(LENGTH_NAME) >> name ) {
          if ( strcmp(name,"end") == 0 ) break;
          else param_list.formattedRead(in_file,name);
          in_file.getline( comment, LENGTH_COMMENT);
        }
      } else in_file.getline(comment,LENGTH_COMMENT);
    }
    continuity=min(continuity,order-1);
    if (order==3) continuity=min(continuity,1);
    int nnodes=(order-continuity)*nelements+continuity+1;
//  limit=min(limit,nnodes);
//  cout << "\tnelements = " << nelements << endl;
//  cout << "\tnnodes = " << nnodes << endl;
//  cout << "\torder = " << order << endl;
//  cout << "\tcontinuity = " << continuity << endl;
//  cout << "\tlimit = " << limit << endl;
//  cout << "\tndigit = " << ndigit << endl;

    if (nelements>1) {
      TRACER_CALL(tr,"solve once");
//    since ncells is determined dynamically, 
//    must use dynamic memory allocation
      double *basis=OPERATOR_NEW_BRACKET(double,(order+1)*order);
      double *basis_deriv=OPERATOR_NEW_BRACKET(double,(order+1)*order);
      double *gauss=OPERATOR_NEW_BRACKET(double,order);
      double *weight=OPERATOR_NEW_BRACKET(double,order);

      bool *dirichlet=OPERATOR_NEW_BRACKET(bool,nnodes);
      int *element_to_node=
        OPERATOR_NEW_BRACKET(int,nelements*(order+1));
      double *pgauss=OPERATOR_NEW_BRACKET(double,nelements*order);
      double *rgauss=OPERATOR_NEW_BRACKET(double,nelements*order);
      double *residual=OPERATOR_NEW_BRACKET(double,nnodes);
      double *soln=OPERATOR_NEW_BRACKET(double,nnodes);
      double *mesh=OPERATOR_NEW_BRACKET(double,nelements+1);

      double *work=OPERATOR_NEW_BRACKET(double,2*nnodes);

#ifdef INDEF
//    initialize arrays to help IEEE exception handling catch unassigned
      {
//      TRACER_CALL(tr0,"default");
        for (int i=0;i<order;i++) {
          gauss[i]=HUGE_VAL;
          weight[i]=HUGE_VAL;
        }
        for (int i=0;i<(order+1)*order;i++) {
          basis[i]=HUGE_VAL;
          basis_deriv[i]=HUGE_VAL;
        }
        for (int i=0;i<nnodes;i++) {
          dirichlet[i]=false;
          residual[i]=HUGE_VAL;
          soln[i]=HUGE_VAL;
        }
        for (int i=0;i<nelements*order;i++) {
          pgauss[i]=HUGE_VAL;
          rgauss[i]=HUGE_VAL;
        }
        for (int i=0;i<nelements*(order+1);i++) element_to_node[i]=-1;
        for (int i=0;i<=nelements;i++) mesh[i]=HUGE_VAL;
        for (int i=0;i<2*nnodes;i++) work[i]=HUGE_VAL;
      }
#endif
      F77NAME(canonical)(continuity,order,
        basis,basis_deriv,gauss,weight);
      F77NAME(grid)(continuity,nelements,nnodes,order,
        element_to_node,mesh);
      F77NAME(initialize)(continuity,nelements,nnodes,order,
        basis,basis_deriv,element_to_node,gauss,mesh,weight,
        dirichlet,pgauss,rgauss,residual,soln);
      F77NAME(precg)(limit,ndigit,nelements,nnodes,order,
        basis,basis_deriv,dirichlet,element_to_node,pgauss,rgauss,
        residual,soln,
        work);
//    for (int i=0;i<nnodes;i++) {
//      cout << "\tsoln[" << i << "] = " << soln[i] << endl;
//    }

//    find min,max data values
      double xlo=mesh[0];
      double xhi=mesh[nelements];
      double dx=(xhi-xlo)/static_cast<double>(NPTS);
      double ulo=F77NAME(solution)(0.);
      double uhi=ulo;
      for (int i=1;i<NPTS;i++) {
        double x=xlo+static_cast<double>(i)*dx;
        double ui=F77NAME(solution)(x);
        ulo=min(ulo,ui);
        uhi=max(uhi,ui);
      }
//    cout << "\tulo,uhi = " << ulo << " " << uhi << endl;

//    setup interactive graphics
      Palette pal;
      XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
      double winsize=0.5;
      XYGraphTool gt("solution","x","u",xlo,xhi,ulo,uhi,
                     &cmap,NULL,winsize);

//    initialize graphics display
      gt.setbgColor("white");
      gt.setfgColor("black");
      gt.drawAxes();

//    draw analytical solution
      gt.setfgColor("red");
      gt.movePen(0.,F77NAME(solution)(0.));
      for (int i=1;i<=NPTS;i++) {
        double x=static_cast<double>(i)*dx;
        gt.drawLine(x,F77NAME(solution)(x));
      }

//    draw numerical solution
      gt.setfgColor("blue");
      double x=0.;
      double y=
        F77NAME(approximation)(continuity,nelements,nnodes,order, x, 
          element_to_node,mesh,soln);
      double err=y-F77NAME(solution)(x);
      double err_min=err;
      double err_max=err;
      gt.movePen(x,y);
      for (int i=1;i<=NPTS;i++) {
        double x=static_cast<double>(i)*dx;
        double y=
          F77NAME(approximation)(continuity,nelements,nnodes,order, x, 
            element_to_node,mesh,soln);
        gt.drawLine(x,y);
        err=y-F77NAME(solution)(x);
        err_min=min(err_min,err);
        err_max=max(err_max,err);
//      cout << "\terr[" << x << "] = " << err_min << " " << err_max
//           << endl;
      }
      gt.flush();

      XYGraphTool gte("error","e","x",xlo,xhi,err_min,err_max,&cmap,NULL,
        winsize);
      gte.setbgColor("white");
      gte.setfgColor("black");
      gte.drawAxes();
      gte.setfgColor("green");

      x=0.;
      y=F77NAME(approximation)(continuity,nelements,nnodes,order, x, 
          element_to_node,mesh,soln)
       -F77NAME(solution)(x);
      double z=F77NAME(solution)(x);
//    cout << "\terr[" << x << "] = " << y << endl;
      gte.movePen(x,y);
      for (int i=1;i<=NPTS;i++) {
        double x=static_cast<double>(i)*dx;
        double y=
          F77NAME(approximation)(continuity,nelements,nnodes,order, x, 
            element_to_node,mesh,soln)
         -F77NAME(solution)(x);
//      cout << "\terr[" << x << "] = " << y << endl;
        gte.drawLine(x,y);
      }
      gte.flush();

/*
      XYGraphTool gtec("error mapped to canonical element","e","xi",0.,1.,
        -1.,1.,&cmap,NULL,winsize);
      gtec.setbgColor("white");
      gtec.setfgColor("black");
      gtec.drawAxes();
      gtec.setfgColor("green");

//    int e=static_cast<int>(static_cast<double>(nelements)
//      *static_cast<double>(rand())/static_cast<double>(RAND_MAX));
//    cout << "\te = " << e << endl;
      for (int e=0;e<nelements;e++) 
      {
        x=mesh[e];
        double dx=(mesh[e+1]-mesh[e]);
        double xi=0.;

        y=F77NAME(approximation)(continuity,nelements,nnodes,order, x, 
            element_to_node,mesh,soln)
         -F77NAME(solution)(x);
        double err_max=abs(y);
        for (int i=1;i<=NPTS;i++) {
          xi=static_cast<double>(i)/static_cast<double>(NPTS);
          x=mesh[e]+xi*dx;
          double y=
            F77NAME(approximation)(continuity,nelements,nnodes,order, x,
              element_to_node,mesh,soln)
           -F77NAME(solution)(x);
          err_max=max(err_max,abs(y));
        }

        x=mesh[e];
        xi=0.;
        y=F77NAME(approximation)(continuity,nelements,nnodes,order, x, 
            element_to_node,mesh,soln)
         -F77NAME(solution)(x);
        gtec.movePen(xi,y/err_max);
        for (int i=1;i<=NPTS;i++) {
          xi=static_cast<double>(i)/static_cast<double>(NPTS);
          x=mesh[e]+xi*dx;
          double y=
            F77NAME(approximation)(continuity,nelements,nnodes,order, x,
              element_to_node,mesh,soln)
           -F77NAME(solution)(x);
          gtec.drawLine(xi,y/err_max);
        }
      }
      gtec.flush();
*/

      wait();

      delete [] work;
      delete [] mesh;
      delete [] soln;
      delete [] residual;
      delete [] rgauss;
      delete [] pgauss;
      delete [] element_to_node;
      delete [] dirichlet;
      delete [] weight;
      delete [] gauss;
      delete [] basis_deriv;
      delete [] basis;
    } // pal,cmap,gt go out of scope here 
    else {
      TRACER_CALL(tr0,"error refinement study");
      int max_refine=10;
      double *basis=OPERATOR_NEW_BRACKET(double,(order+1)*order);
      double *basis_deriv=OPERATOR_NEW_BRACKET(double,(order+1)*order);
      double *gauss=OPERATOR_NEW_BRACKET(double,order);
      double *weight=OPERATOR_NEW_BRACKET(double,order);

      double *log_nnodes=OPERATOR_NEW_BRACKET(double,max_refine);
      double *error_infinity=OPERATOR_NEW_BRACKET(double,max_refine);
      double *error_1=OPERATOR_NEW_BRACKET(double,max_refine);
      double *error_2=OPERATOR_NEW_BRACKET(double,max_refine);
      double *error_nodes=OPERATOR_NEW_BRACKET(double,max_refine);

      F77NAME(canonical)(continuity,order,
        basis,basis_deriv,gauss,weight);

      int nelements=2;
      double log10=log(10.);
      for (int k=0;k<max_refine;nelements*=2,k++) {
        nnodes=(order-continuity)*nelements+continuity+1;
//      limit=nnodes;

        bool *dirichlet=OPERATOR_NEW_BRACKET(bool,nnodes);
        int *element_to_node=
          OPERATOR_NEW_BRACKET(int,nelements*(order+1));
        double *pgauss=OPERATOR_NEW_BRACKET(double,nelements*order);
        double *rgauss=OPERATOR_NEW_BRACKET(double,nelements*order);
        double *residual=OPERATOR_NEW_BRACKET(double,nnodes);
        double *soln=OPERATOR_NEW_BRACKET(double,nnodes);
        double *mesh=OPERATOR_NEW_BRACKET(double,nelements+1);
        double *work=OPERATOR_NEW_BRACKET(double,2*nnodes);
        F77NAME(grid)(continuity,nelements,nnodes,order,
          element_to_node,mesh);
        F77NAME(initialize)(continuity,nelements,nnodes,order,
          basis,basis_deriv,element_to_node,gauss,mesh,weight,
          dirichlet,pgauss,rgauss,residual,soln);
        F77NAME(precg)(limit,ndigit,nelements,nnodes,order,
          basis,basis_deriv,dirichlet,element_to_node,pgauss,rgauss,
          residual,soln,
          work);

        error_infinity[k]=0.;
        error_1[k]=0.;
        error_2[k]=0.;
        for (int element=0;element<nelements;element++) {
          double dx=mesh[element+1]-mesh[element];
          for (int ngauss=0;ngauss<order;ngauss++) {
            double xg=mesh[element]+gauss[ngauss]*dx;
            double ei=abs(F77NAME(solution)(xg)
                         -F77NAME(approximation)(continuity,nelements,
                           nnodes,order,xg,element_to_node,mesh,soln));
            error_infinity[k]=max(error_infinity[k],ei);
            error_1[k]+=ei*weight[ngauss]*dx;
            error_2[k]+=ei*ei*weight[ngauss]*dx;
          }
        }
        log_nnodes[k]=log(static_cast<double>(nnodes))/log10;
        error_infinity[k]=log(error_infinity[k])/log10;
        error_1[k]=log(error_1[k])/log10;
        error_2[k]=log(sqrt(error_2[k]))/log10;

        if (continuity==0) {
          error_nodes[k]=0.;
          for (int node=0;node<=nelements;node++) {
            double xg=mesh[node];
            double ei=abs(F77NAME(solution)(xg)
                         -F77NAME(approximation)(continuity,nelements,
                           nnodes,order,xg,element_to_node,mesh,soln));
            error_nodes[k]=max(error_nodes[k],ei);
          }
          if (order==2) {
            double frac=0.5;
            for (int element=0;element<=nelements;element++) {
              double xg=mesh[element]
                       +frac*(mesh[element+1]-mesh[element]);
              double ei=abs(F77NAME(solution)(xg)
                           -F77NAME(approximation)(continuity,nelements,
                            nnodes,order,xg,element_to_node,mesh,soln));
              error_nodes[k]=max(error_nodes[k],ei);
            }
          } else if (order==3) {
            double frac=0.5*(1.-sqrt(0.2));
//          for (int element=0;element<=nelements;element++) {
//            double xg=mesh[element]
//                     +frac*(mesh[element+1]-mesh[element]);
//            double ei=abs(F77NAME(solution)(xg)
//                         -F77NAME(approximation)(continuity,nelements,
//                          nnodes,order,xg,element_to_node,mesh,soln));
//            error_nodes[k]=max(error_nodes[k],ei);
//          }
//          frac=0.5*(1.+sqrt(0.2));
//          for (int element=0;element<=nelements;element++) {
//            double xg=mesh[element]
//                     +frac*(mesh[element+1]-mesh[element]);
//            double ei=abs(F77NAME(solution)(xg)
//                         -F77NAME(approximation)(continuity,nelements,
//                          nnodes,order,xg,element_to_node,mesh,soln));
//            error_nodes[k]=max(error_nodes[k],ei);
//          }
          }
        } else {
          if (order==3) {
            double frac=0.5*(1.-sqrt(1./3.));
            for (int element=0;element<=nelements;element++) {
              double xg=mesh[element]
                       +frac*(mesh[element+1]-mesh[element]);
              double ei=abs(F77NAME(solution)(xg)
                           -F77NAME(approximation)(continuity,nelements,
                            nnodes,order,xg,element_to_node,mesh,soln));
              error_nodes[k]=max(error_nodes[k],ei);
            }
            frac=0.5*(1.+sqrt(1./3.));
            for (int element=0;element<=nelements;element++) {
              double xg=mesh[element]
                       +frac*(mesh[element+1]-mesh[element]);
              double ei=abs(F77NAME(solution)(xg)
                           -F77NAME(approximation)(continuity,nelements,
                            nnodes,order,xg,element_to_node,mesh,soln));
              error_nodes[k]=max(error_nodes[k],ei);
            }
          }
        }
        error_nodes[k]=log(error_nodes[k])/log10;
//      cout << "\terror[" << k << "] = " << error_1[k] << " "
//           << error_2[k] << " " << error_infinity[k] << endl;
        cout << "\terror_nodes[" << k << "] = " << error_nodes[k] 
             << endl;
        delete [] dirichlet;
        delete [] element_to_node;
        delete [] pgauss;
        delete [] rgauss;
        delete [] residual;
        delete [] soln;
        delete [] mesh;
        delete [] work;
      }
      nelements/=2;
      Palette pal;
      XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
      double winsize=0.5;
      double bs=(log_nnodes[max_refine-1]-log_nnodes[0])
               /static_cast<double>(10*max_refine);
/*
      {
        double elo=0.;
        double ehi=0.;
        for (int k=0;k<max_refine;k++) {
          elo=min(elo,error_infinity[k]);
          ehi=max(ehi,error_infinity[k]);
        }
        XYGraphTool gt("L^infinity error","e","h",log_nnodes[0],
          log_nnodes[max_refine-1],elo,ehi,&cmap,NULL,winsize);
        gt.setbgColor("white");
        gt.setfgColor("black");
        gt.drawAxes();
        gt.setfgColor("blue");
        gt.movePen(log_nnodes[0],error_infinity[0]);
        for (int k=1;k<max_refine;k++) {
          gt.drawLine(log_nnodes[k],error_infinity[k]);
          cout << "\tslope[" << k << "] = " 
               << (error_infinity[k]-error_infinity[k-1])
                 /(log_nnodes[k]-log_nnodes[k-1]) << endl;
        }

        gt.drawBoxGivenCenter(log_nnodes[0],error_infinity[0],bs);
        for (int k=1;k<max_refine;k++) {
          gt.drawBoxGivenCenter(log_nnodes[k],error_infinity[k],bs);
        }
        gt.flush();
        wait();
      }
*/
/*
      {
        double elo=0.;
        double ehi=0.;
        for (int k=0;k<max_refine;k++) {
          elo=min(elo,error_1[k]);
          ehi=max(ehi,error_1[k]);
        }
        XYGraphTool gt("L^1 error","e","h",log_nnodes[0],
          log_nnodes[max_refine-1],elo,ehi,&cmap,NULL,winsize);
        gt.setbgColor("white");
        gt.setfgColor("black");
        gt.drawAxes();
        gt.setfgColor("blue");
        gt.movePen(log_nnodes[0],error_1[0]);
        for (int k=1;k<max_refine;nelements*=2,k++) {
          gt.drawLine(log_nnodes[k],error_1[k]);
          cout << "\tslope[" << k << "] = " 
               << (error_1[k]-error_1[k-1])
                 /(log_nnodes[k]-log_nnodes[k-1]) 
               << endl;
        }

        gt.drawBoxGivenCenter(log_nnodes[0],error_1[0],bs);
        for (int k=1;k<max_refine;k++) {
          gt.drawBoxGivenCenter(log_nnodes[k],error_1[k],bs);
        }
        gt.flush();
        gt.flush();
        wait();
      }
*/
//
      {
        double elo=0.;
        double ehi=0.;
        for (int k=0;k<max_refine;k++) {
          elo=min(elo,error_2[k]);
          ehi=max(ehi,error_2[k]);
        }
        XYGraphTool gt("L^2 error","e","h",log_nnodes[0],
          log_nnodes[max_refine-1],elo,ehi,&cmap,NULL,winsize);
        gt.setbgColor("white");
        gt.setfgColor("black");
        gt.drawAxes();
        gt.setfgColor("blue");
        gt.movePen(log_nnodes[0],error_2[0]);
        for (int k=1;k<max_refine;nelements*=2,k++) {
          gt.drawLine(log_nnodes[k],error_2[k]);
          cout << "\tslope[" << k << "] = " 
               << (error_2[k]-error_2[k-1])
                 /(log_nnodes[k]-log_nnodes[k-1]) 
               << endl;
        }

        gt.drawBoxGivenCenter(log_nnodes[0],error_2[0],bs);
        for (int k=1;k<max_refine;k++) {
          gt.drawBoxGivenCenter(log_nnodes[k],error_2[k],bs);
        }
        gt.flush();
        gt.flush();
        wait();
      }
//
      {
        double elo=0.;
        double ehi=0.;
        for (int k=0;k<max_refine;k++) {
          elo=min(elo,error_nodes[k]);
          ehi=max(ehi,error_nodes[k]);
        }
        XYGraphTool gt("L^infinity error at nodes","e","h",log_nnodes[0],
          log_nnodes[max_refine-1],elo,ehi,&cmap,NULL,winsize);
        gt.setbgColor("white");
        gt.setfgColor("black");
        gt.drawAxes();
        gt.setfgColor("blue");
        gt.movePen(log_nnodes[0],error_nodes[0]);
        for (int k=1;k<max_refine;nelements*=2,k++) {
          gt.drawLine(log_nnodes[k],error_nodes[k]);
          cout << "\tslope[" << k << "] = " 
               << (error_nodes[k]-error_nodes[k-1])
                 /(log_nnodes[k]-log_nnodes[k-1]) 
               << endl;
        }

        gt.drawBoxGivenCenter(log_nnodes[0],error_nodes[0],bs);
        for (int k=1;k<max_refine;k++) {
          gt.drawBoxGivenCenter(log_nnodes[k],error_nodes[k],bs);
        }
        gt.flush();
        gt.flush();
        wait();
      }
      delete [] basis;
      delete [] basis_deriv;
      delete [] gauss;
      delete [] weight;

      delete [] log_nnodes;
      delete [] error_infinity;
      delete [] error_1;
      delete [] error_2;
      delete [] error_nodes;
    }

  } // scope of MemoryDebugger
  return EXIT_SUCCESS;
}

template class InputParameter<int>;
template class InputParameter<double>;
template class InputParameter<bool>;

extern "C" {
void MAIN__(void) {;}
}
