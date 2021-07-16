#include <iostream>
#include <limits>
//#include <math.h>
#include "BandMatrix.H"
#include "Debug.H"
#include "GUIEnumInputParameter.H"
#include "GUIInputParameter.H"
#include "GUIInputs.H"
#include "MemoryDebugger.H"
#include "Polynomial.H"
#include "Quadrature.H"
#include "SetTraps.H"
#include "TimedObject.H"
#include "Tracer.H"
#include "Vector.H"
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

enum POLYNOMIAL_FAMILY{LAGRANGE_POLYNOMIAL,HIERARCHICAL_POLYNOMIAL};
POLYNOMIAL_FAMILY poly_family=LAGRANGE_POLYNOMIAL;
int ipoly_family=static_cast<int>(poly_family);
const char *poly_family_names[2]={"lagrange","hierarchical"};

enum QUADRATURE_RULE{GAUSSIAN_QUADRATURE,LOBATTO_QUADRATURE,
  NEWTON_COTES_QUADRATURE,CLENSHAW_CURTIS_QUADRATURE};
QUADRATURE_RULE quad_rule=GAUSSIAN_QUADRATURE;
int iquad_rule=static_cast<int>(quad_rule);
const char *quad_rule_names[4]=
  {"gaussian","lobatto","newton-cotes","clenshaw-curtis"};

enum BOUNDARY_CONDITION{ESSENTIAL_BOUNDARY_CONDITION,
  NATURAL_BOUNDARY_CONDITION};
BOUNDARY_CONDITION left_bc=ESSENTIAL_BOUNDARY_CONDITION;
int ileft_bc=static_cast<int>(left_bc);
double left_bc_val=0.;
BOUNDARY_CONDITION right_bc=ESSENTIAL_BOUNDARY_CONDITION;
int iright_bc=static_cast<int>(right_bc);
double right_bc_val=0.;
const char *bc_names[2]={"essential","natural"};

enum LAGRANGE_POLYNOMIAL_NODES{EQUALLY_SPACED_NODES,CHEBYSHEV_NODES};
LAGRANGE_POLYNOMIAL_NODES poly_nodes=EQUALLY_SPACED_NODES;
int ipoly_nodes=poly_nodes;
const char *poly_node_names[2]={"equally spaced","chebyshev"};

enum MESH_TYPE{UNIFORM_MESH,RANDOM_MESH};
MESH_TYPE mesh_type=UNIFORM_MESH;
int imesh_type=static_cast<int>(mesh_type);
const char *mesh_name[2]={"uniform","random"};

TimedObject *assembly_timing=0;
TimedObject *solver_timing=0;

int nelements=1;
int nquad_pts=1;
int order=1;
int max_refine=10;

double p(double x) { return 1.; }
//
double r(double x) { return 1.; }
double f(double x) { return 1.; }
double solution(double x) {
//solve -(p u')' + r u = f:
  switch(left_bc) {
    case ESSENTIAL_BOUNDARY_CONDITION: {
      switch (right_bc) {
        case ESSENTIAL_BOUNDARY_CONDITION: {
          double A=left_bc_val-1.;
          double B=(right_bc_val-1.-A*cosh(1.))/sinh(1.);
          return 1.+A*cosh(x)+B*sinh(x);
        }
        case NATURAL_BOUNDARY_CONDITION: {
          double A=left_bc_val-1.;
          double B=(right_bc_val-A*sinh(1.))/cosh(1.);
          return 1.+A*cosh(x)+B*sinh(x);
        }
      }
    }
    case NATURAL_BOUNDARY_CONDITION: {
      switch (right_bc) {
        case ESSENTIAL_BOUNDARY_CONDITION: {
          double B=left_bc_val;
          double A=(right_bc_val-1.-B*sinh(1.))/cosh(1.);
          return 1.+A*cosh(x)+B*sinh(x);
        }
        case NATURAL_BOUNDARY_CONDITION: {
          double B=left_bc_val;
          double A=(right_bc_val-B*cosh(1.))/sinh(1.);
          return 1.+A*cosh(x)+B*sinh(x);
        }
      }
    }
  }
}
//
/*
double r(double x) { return 1.; }
double f(double x) {
  if (order<=1) return x;
  return pow(x,order-2)*(x*x-static_cast<double>(order*(order-1)));
}
double solution(double x) { return pow(x,order); }
*/
/*
double r(double x) { return 0.; }
double f(double x) {
  return  -pow(x,order-2)*static_cast<double>(order*(order-1));
}
double solution(double x) { return pow(x,order); }
*/

#define NPTS 1000

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

  { const char *group="Polynomial Parameters";
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
      ipoly_family,"polynomial family",poly_family_names[0],
        poly_family_names[1]));
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
      ipoly_nodes,"Lagrange poly nodes",poly_node_names[0],
      poly_node_names[1]));
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(order,
      "order",1,INT_MAX,group) );
  }
  { const char *group="Quadrature Parameters";
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
      iquad_rule,"quadrature rule",quad_rule_names[0],
      quad_rule_names[1],quad_rule_names[2],quad_rule_names[3]));
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(nquad_pts,
      "number quadrature points",1,INT_MAX,group) );
  }
  { const char *group="Mesh Parameters";
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
      imesh_type,"mesh type",mesh_name[0],mesh_name[1]));
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(nelements,
      "nelements (<=1 ==> refinement study)",0,INT_MAX,group) );
  }
  { const char *group="Boundary Conditions";
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
      ileft_bc,"left type",bc_names[0],bc_names[1]));
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(left_bc_val,
      "left value",-DBL_MAX,DBL_MAX,group) );
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
      iright_bc,"right type",bc_names[0],bc_names[1]));
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(right_bc_val,
      "right value",-DBL_MAX,DBL_MAX,group) );
  }
  { const char *group="Refinement Study";
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(max_refine,
      "number refinement meshes",2,INT_MAX,group) );
  }
  { const char *group="Graphics";
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(winsize,
      "winsize",0.,1.,group));
  }
#ifdef DEBUG
//main_list->printOn(cout);
#endif
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void checkMainInput() {
//TRACER_CALL(t,"checkMainInput");
  poly_family=static_cast<POLYNOMIAL_FAMILY>(ipoly_family);
  poly_nodes=static_cast<LAGRANGE_POLYNOMIAL_NODES>(ipoly_nodes);
  quad_rule=static_cast<QUADRATURE_RULE>(iquad_rule);
  mesh_type=static_cast<MESH_TYPE>(imesh_type);
  left_bc=static_cast<BOUNDARY_CONDITION>(ileft_bc);
  right_bc=static_cast<BOUNDARY_CONDITION>(iright_bc);
  switch (quad_rule ) {
    case GAUSSIAN_QUADRATURE:
      nquad_pts=max(1,nquad_pts);
      break;
    case LOBATTO_QUADRATURE:
    case CLENSHAW_CURTIS_QUADRATURE:
      nquad_pts=max(2,nquad_pts);
      break;
    case NEWTON_COTES_QUADRATURE:
      nquad_pts=max(2,min(9,nquad_pts));
      break;
  }
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double approximation(int e,double xi,const Polynomial &poly,
const Vector<double,double> &soln) {
//TRACER_CALL(t,"approximation");
#ifdef DEBUG
//for (int i=0;i<soln.size();i++) {
//  cout << "\tsoln[" << i << "] = " << soln[i] << endl;
//}
#endif
  NumPtr<double> basis_poly_values(order+1);
  poly.values(xi,basis_poly_values);
  int offset=e*order;
  double sum=0.;
  for (int dof=0;dof<=order;dof++) {
    sum+=basis_poly_values[dof]*soln[offset+dof];
  }
  return sum;
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double finiteElementSolution(const Polynomial &poly,
const Quadrature<1> &quad,double *&mesh,Vector<double,double> *&soln) {
//TRACER_CALL(t,"finiteElementSolution");
  mesh=OPERATOR_NEW_BRACKET(double,nelements+1);
  mesh[0]=0.;
  mesh[nelements]=1.;
  switch (mesh_type) {
    case UNIFORM_MESH: {
      double dx=1./static_cast<double>(nelements);
      double xi=dx;
      for (int i=1;i<nelements;i++,xi+=dx) mesh[i]=xi;
      break;
    }
    case RANDOM_MESH: {
      NumPtrC<double> mesh_points(nelements-1);
      for (int i=0;i<nelements-1;i++) mesh_points[i]=drand48();
      mesh_points.heapSort();
      for (int i=1;i<nelements;i++) mesh[i]=mesh_points[i-1];
      break;
    }
  }
#ifdef DEBUG
//for (int e=0;e<=nelements;e++) {
//  cout << "\tmesh[" << e << "] = " << mesh[e] << endl;
//}
#endif

  int nnodes=order*nelements+1;
#ifdef DEBUG
//cout << "\tnelements,nnodes = " << nelements << " " << nnodes << endl;
#endif
  SymmetricPositiveBandMatrix<double,double> matrix(nnodes,order+1,0.);
  Vector<double,double> rhs(nnodes,0.);

  NumPtr<double> basis_poly_values(order+1);
  NumPtr<double> basis_poly_slopes(order+1);
  {
    Timer tim(assembly_timing);
    for (int q=0;q<quad.numberPoints();q++) {
      double xi=quad.point(q)[0];
      poly.values(xi,basis_poly_values);
      poly.slopes(xi,basis_poly_slopes);
      double wq=quad.weight(q);
#ifdef DEBUG
//    cout << "\txi[" << q << "] = " << xi << endl;
//    cout << "\twq[" << q << "] = " << wq << endl;
//    cout << "\tbasis_poly_values = ";
//    for (int dof=0;dof<=order;dof++) {
//      cout << basis_poly_values[dof] << " ";
//    }
//    cout << endl;
//    cout << "\tbasis_poly_slopes = ";
//    for (int dof=0;dof<=order;dof++) {
//      cout << basis_poly_slopes[dof] << " ";
//    }
//    cout << endl;
#endif
      for (int e=0;e<nelements;e++) {
        double he=mesh[e+1]-mesh[e];
        double xq=mesh[e]+xi*he;
        double pq=p(xq)*wq/he;
        double rq=r(xq)*wq*he;
        double fq=f(xq)*wq*he;
        int offset=e*order;
        for (int dofj=0;dofj<=order;dofj++) {
          int column=offset+dofj;
          double ps=pq*basis_poly_slopes[dofj];
          double rv=rq*basis_poly_values[dofj];
          for (int dofi=dofj;dofi<=order;dofi++) {
            int row=offset+dofi;
            matrix(row,column)+=//SymmetricBandMatrix stores on & below diag
              basis_poly_slopes[dofi]*ps+basis_poly_values[dofi]*rv;
          }
          rhs[column]+=basis_poly_values[dofj]*fq;
        }
      }
    }
#ifdef DEBUG
//  for (int row=0;row<matrix.size(0);row++) {
//    cout << "matrix(" << row << ",*) = ";
//    for (int col=max(0,row-order);col<=row;col++) {
//      cout << matrix(row,col) << " ";
//    }
//    cout << endl;
//  }
//  for (int row=0;row<matrix.size(0);row++) {
//    cout << "rhs[" << row << "] = " << rhs[row] << endl;
//  }
#endif
    switch (left_bc) {
      case ESSENTIAL_BOUNDARY_CONDITION: {
//      TRACER_CALL(tle,"left bc");
#ifdef DEBUG
//      cout << "\tleft_bc_val = " << left_bc_val << endl;
#endif
        int essential_dof=0;
        rhs[0]=matrix(essential_dof,essential_dof)*left_bc_val;
        for (int dofi=0;dofi<=order;dofi++) {
          if (dofi!=essential_dof) {
            rhs[dofi]-=matrix(dofi,essential_dof)*left_bc_val;
            matrix(dofi,essential_dof)=0.;
          }
        }
        break;
      }
      case NATURAL_BOUNDARY_CONDITION: {
//      TRACER_CALL(tle,"left bc");
#ifdef DEBUG
//      cout << "\tleft_bc_val = " << left_bc_val << endl;
#endif
        poly.values(0.,basis_poly_values);
        for (int dofi=0;dofi<=order;dofi++) {
          rhs[dofi]+=basis_poly_values[dofi]*left_bc_val;
        }
        break;
      }
    }
    switch (right_bc) {
      case ESSENTIAL_BOUNDARY_CONDITION: {
//      TRACER_CALL(tre,"right bc");
#ifdef DEBUG
//      cout << "\tright_bc_val = " << right_bc_val << endl;
#endif
        int essential_dof=(poly_family==LAGRANGE_POLYNOMIAL ? order : 1);
        int offset=order*(nelements-1);
        int rowe=offset+essential_dof;
        rhs[rowe]=matrix(rowe,rowe)*right_bc_val;
        for (int dofj=0;dofj<=order;dofj++) {
          int col=offset+dofj;
          if (dofj<essential_dof) {
            rhs[col]-=matrix(rowe,col)*right_bc_val;
            matrix(rowe,col)=0.;
          } else if (dofj>essential_dof) {
            rhs[col]-=matrix(col,rowe)*right_bc_val;
            matrix(col,rowe)=0.;
          }
        }
        break;
      }
      case NATURAL_BOUNDARY_CONDITION: {
//      TRACER_CALL(tre,"right bc");
#ifdef DEBUG
//      cout << "\tright_bc_val = " << right_bc_val << endl;
#endif
        poly.values(1.,basis_poly_values);
        for (int dofi=0;dofi<=order;dofi++) {
          rhs[dofi]+=basis_poly_values[dofi]*right_bc_val;
        }
        break;
      }
    }
  }
#ifdef DEBUG
//for (int row=0;row<matrix.size(0);row++) {
//  cout << "matrix(" << row << ",*) = ";
//  for (int col=max(0,row-order);col<=row;col++) {
//    cout << matrix(row,col) << " ";
//  }
//  cout << endl;
//}
//for (int row=0;row<matrix.size(0);row++) {
//  cout << "rhs[" << row << "] = " << rhs[row] << endl;
//}
#endif
  soln=OPERATOR_NEW Vector<double,double>(nnodes);
#ifdef DEBUG
//double dx=1./static_cast<double>(nnodes-1);
//double x=0.;
//for (int i=0;i<nnodes;i++,x+=dx) (*soln)[i]=pow(x,order);
//Vector<double,double> *err=matrix*(*soln);
//(*err)-=rhs;
//for (int i=0;i<nnodes;i++) {
//  cout << "\terr[" << i << "] = " << (*err)[i] << endl;
//}
//delete err; err=0;
#endif
  {
    Timer tim(solver_timing);
    matrix.solve(rhs,*soln);
  }
#ifdef DEBUG
//for (int i=0;i<nnodes;i++) {
//  cout << "\tsoln[" << i << "] = " << (*soln)[i] << endl;
//}
#endif
  BandMatrix<double,double> *band_matrix=matrix.makeBandMatrix();
  double cni=band_matrix->reciprocalConditionNumber('1');
  delete band_matrix; band_matrix=0;
  return cni;
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void runMain(bool /*called_before*/) {
//TimedObjectList tol;
  assembly_timing=OPERATOR_NEW TimedObject("assembly timing");
  solver_timing=OPERATOR_NEW TimedObject("solver timing");
  Polynomial *poly=0;
  switch (poly_family) {
    case LAGRANGE_POLYNOMIAL: {
      switch (poly_nodes) {
        case EQUALLY_SPACED_NODES: {
          poly=OPERATOR_NEW C0LagrangePolynomial(order);
#ifdef DEBUG
//        double dxi=1./static_cast<double>(order);
//        double xi=0.;
//        for (int i=0;i<=order;i++,xi+=dxi) {
//          cout << "\n\txi = " << xi << endl;
//          NumPtr<double> v(order+1);
//          poly->values(xi,v);
//          for (int j=0;j<=order;j++) {
//            cout << "v[" << j << "] = " << v[j] << endl; 
//          }
//          NumPtr<double> vv(order+1);
//          NumPtr<double> s(order+1);
//          poly->values(xi+1.e-7,vv);
//          poly->slopes(xi,s);
//          for (int j=0;j<=order;j++) {
//            cout << "s[" << j << "] = " << s[j] << " "
//                 << (vv[j]-v[j])*1.e7 << endl; 
//          }
//        }
#endif
          break;
        }
        case CHEBYSHEV_NODES: {
//        TRACER_CALL(tr,"poly_family");
          NumPtr<double> nodes(order+1);
          nodes[0]=0.;
          nodes[order]=1.;
          if (order>1) {
            double dtheta=M_PI/static_cast<double>(order-1.);
            double theta=0.5*dtheta;
            for (int j=1;j<order;j++,theta+=dtheta) {
              nodes[order-j]=.5*cos(theta)+.5;
            }
          }
#ifdef DEBUG
//        for (int i=0;i<nodes.getNumber();i++) {
//          cout << "\tnodes[" << i << "] = " << nodes[i] << endl;
//        }
#endif
          poly=OPERATOR_NEW C0LagrangePolynomial(nodes);
          break;
        }
      }
      break;
    }
    case HIERARCHICAL_POLYNOMIAL:
      poly=OPERATOR_NEW HierarchicalPolynomial();
      break;
  }

  Quadrature<1> *quad=0;
  switch (quad_rule) {
    case GAUSSIAN_QUADRATURE :
      quad=OPERATOR_NEW GaussianQuadrature<1>(nquad_pts);
#ifdef DEBUG
//    for (int i=0;i<2*nquad_pts;i++) {
//      double sum=0.;
//      for (int q=0;q<nquad_pts;q++) {
//        sum+=pow(quad->point(q)[0],i)*quad->weight(q);
//      }
//      cout << "\t1 / integral of x^" << i << " = " << 1./sum << endl;
//    }
#endif
      break;
    case LOBATTO_QUADRATURE :
      quad=OPERATOR_NEW LobattoQuadrature<1>(nquad_pts);
      break;
    case NEWTON_COTES_QUADRATURE :
      quad=OPERATOR_NEW NewtonCotesQuadrature<1>(nquad_pts);
      break;
    case CLENSHAW_CURTIS_QUADRATURE :
      quad=OPERATOR_NEW ClenshawCurtisQuadrature<1>(nquad_pts);
      break;
  }
#ifdef DEBUG
//quad->printOn(cout);
#endif
//cout << "\tnelements = " << nelements << endl;
//cout << "\tnnodes = " << nnodes << endl;
//cout << "\torder = " << order << endl;
//cout << "\tcontinuity = " << continuity << endl;

  if (nelements>1) { // computation for a single mesh
//  TRACER_CALL(tr,"solve once");
    Vector<double,double> *soln=0;
    double *mesh=0;
    double condition_number_inverse=
      finiteElementSolution(*poly,*quad,mesh,soln);
//  find min,max data values
    double xlo=mesh[0];
    double xhi=mesh[nelements];
    double dx=(xhi-xlo)/static_cast<double>(NPTS);
    double ulo=solution(0.);
    double uhi=ulo;
    double x=dx;
    for (int i=1;i<NPTS;i++,x+=dx) {
      double ui=solution(x);
      ulo=min(ulo,ui);
      uhi=max(uhi,ui);
    }
//  cout << "\tulo,uhi = " << ulo << " " << uhi << endl;

//  setup interactive graphics
    Palette pal;
    XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
    double winsize=0.5;
    XYGraphTool gt("solution","x","u",xlo,xhi,ulo,uhi,
                   &cmap,NULL,winsize);

//  initialize graphics display
    gt.setbgColor("white");
    gt.setfgColor("black");
    gt.drawAxes();

//  draw analytical solution
    gt.setfgColor("red");
    gt.movePen(0.,solution(0.));
    x=dx;
    for (int i=1;i<=NPTS;i++,x+=dx) {
      gt.drawLine(x,solution(x));
    }

//  draw numerical solution
    gt.setfgColor("blue");
    int e=0;
    double xi=0.;
    double y=approximation(e,xi,*poly,*soln);
    x=0.;
    double err=y-solution(x);
    double err_min=err;
    double err_max=err;
    gt.movePen(x,y);
    x=dx;
    for (int i=1;i<=NPTS;i++,x+=dx) {
      while (e<nelements-1 && x>mesh[e+1]) e++;
      xi=(x-mesh[e])/(mesh[e+1]-mesh[e]);
      y=approximation(e,xi,*poly,*soln);
      gt.drawLine(x,y);
      err=y-solution(x);
      err_min=min(err_min,err);
      err_max=max(err_max,err);
//    cout << "\terr[" << x << "] = " << err_min << " " << err_max
//         << endl;
    }
    gt.flush();

    XYGraphTool gte("error","x","e",xlo,xhi,err_min,err_max,&cmap,NULL,
      winsize);
    gte.setbgColor("white");
    gte.setfgColor("black");
    gte.drawAxes();
    gte.setfgColor("green");

    e=0;
    x=0.;
    xi=0.;
    y=approximation(e,xi,*poly,*soln)-solution(x);
//  cout << "\terr[" << x << "] = " << y << endl;
    gte.movePen(x,y);
    x=dx;
    for (int i=1;i<=NPTS;i++,x+=dx) {
      while (e<nelements-1 && x>mesh[e+1]) e++;
      xi=(x-mesh[e])/(mesh[e+1]-mesh[e]);
      y=approximation(e,xi,*poly,*soln)-solution(x);
//    cout << "\terr[" << x << "] = " << y << endl;
      gte.drawLine(x,y);
    }
    gte.flush();
    XYGraphTool::WINDOW_TYPE::QuitButton qb;

    delete [] mesh; mesh=0;
    delete soln; soln=0;
  } else {
//  TRACER_CALL(tr0,"error refinement study");
    double *log_nnodes=OPERATOR_NEW_BRACKET(double,max_refine);
//  double *error_infinity=OPERATOR_NEW_BRACKET(double,max_refine);
//  double *error_1=OPERATOR_NEW_BRACKET(double,max_refine);
    double *error_2=OPERATOR_NEW_BRACKET(double,max_refine);
    double *conno=OPERATOR_NEW_BRACKET(double,max_refine);

    nelements=2;
    double log10=log(10.);
    for (int k=0;k<max_refine;nelements*=2,k++) {
#ifdef DEBUG
      cout << "\tk,nelements = " << k << " " << nelements << endl;
#endif
      double *mesh=0;
      Vector<double,double> *soln=0;
      conno[k]=finiteElementSolution(*poly,*quad,mesh,soln);
      conno[k]=-log(conno[k])/log10;

//    error_infinity[k]=0.;
//    error_1[k]=0.;
      error_2[k]=0.;
      for (int e=0;e<nelements;e++) {
        double dx=mesh[e+1]-mesh[e];
        for (int q=0;q<quad->numberPoints();q++) {
          double xi=quad->point(q)[0];
          double wt=quad->weight(q)*dx;
          double x=mesh[e]+xi*dx;
          double ei=abs(solution(x)-approximation(e,xi,*poly,*soln));
//        error_infinity[k]=max(error_infinity[k],ei);
//        error_1[k]+=ei*wt;
          error_2[k]+=ei*ei*wt;
        }
      }
      log_nnodes[k]=log(static_cast<double>(soln->size()))/log10;
//    error_infinity[k]=log(error_infinity[k])/log10;
//    error_1[k]=log(error_1[k])/log10;
      error_2[k]=log(sqrt(error_2[k]))/log10;
      delete [] mesh; mesh=0;
      delete soln; soln=0;
    }
    nelements/=2;
    Palette pal;
    XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
    double winsize=0.5;
    double bs=(log_nnodes[max_refine-1]-log_nnodes[0])
             /static_cast<double>(10*max_refine);
    {
      double elo=0.;
      double ehi=0.;
      for (int k=0;k<max_refine;k++) {
        elo=min(elo,error_2[k]);
        ehi=max(ehi,error_2[k]);
      }
      XYGraphTool gt("L^2 error versus mesh width","h","e",log_nnodes[0],
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

      double clo=0.;
      double chi=0.;
      for (int k=0;k<max_refine;k++) {
        clo=min(clo,conno[k]);
        chi=max(chi,conno[k]);
      }
      XYGraphTool gtc("condition number vs mesh width","h","kappa",
        log_nnodes[0],log_nnodes[max_refine-1],clo,chi,&cmap,NULL,winsize);
      gtc.setbgColor("white");
      gtc.setfgColor("black");
      gtc.drawAxes();
      gtc.setfgColor("green");
      gtc.movePen(log_nnodes[0],conno[0]);
      for (int k=1;k<max_refine;nelements*=2,k++) {
        gtc.drawLine(log_nnodes[k],conno[k]);
      }
      gtc.drawBoxGivenCenter(log_nnodes[0],conno[0],bs);
      for (int k=1;k<max_refine;k++) {
        gtc.drawBoxGivenCenter(log_nnodes[k],conno[k],bs);
      }
      gtc.flush();
      XYGraphTool::WINDOW_TYPE::QuitButton qb;
    }
/*
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
      XYGraphTool::WINDOW_TYPE::QuitButton qb;
    }
*/
    delete [] log_nnodes; log_nnodes=0;
//  delete [] error_infinity; error_infinity=0;
//  delete [] error_1; error_1=0;
    delete [] error_2; error_2=0;
    delete [] conno; conno=0;
  }
  delete poly; poly=0;
  delete quad; quad=0;
  assembly_timing->printOn(cout);
  solver_timing->printOn(cout);
  delete assembly_timing; assembly_timing=0;
  delete solver_timing; solver_timing=0;
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
//F77NAME(machine).roundoff=DBL_EPSILON;
//F77NAME(machine).small=DBL_MIN;
//F77NAME(machine).huge=DBL_MAX;
//F77NAME(machine).undefind=HUGE_VAL;
//F77NAME(machine).pi=M_PI;

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
//#include "InputParameter.C"
//template class InputParameter<int>;
//template class InputParameter<double>;
template class InputParameter<bool>;
