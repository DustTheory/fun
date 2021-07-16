// "$Header:$"
//**********************************************************************
// Copyright 2009 John A. Trangenstein
//
// This software is made available for research and instructional use 
// only. 
// You may copy and use this software without charge for these 
// non-commercial purposes, provided that the copyright notice and 
// associated text is reproduced on all copies.  
// For all other uses (including distribution of modified versions), 
// please contact the author at
//   John A. Trangenstein
//   Department of Mathematics
//   Duke University
//   Durham, NC 27708-0320
//   USA
// or
//   johnt@math.duke.edu
// 
// This software is made available "as is" without any assurance that it
// is completely correct, or that it will work for your purposes.  
// Use the software at your own risk.
//**********************************************************************

#include <cmath>
#include <iostream>
#include <math.h> // for HUGE_VAL,M_PI
#include <stdlib.h>

#include "Debug.H"
#include "GUIEnumInputParameter.H"
#include "GUIInputs.H"
#include "MemoryDebugger.H"
#include "Quadrature.H"
#include "SetTraps.H"
#include "TimedObject.H"
#include "Tracer.H"
#include "Types.H"
#include "XColormap.H"
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

#if (SPACEDIM==1)
enum QUADRATURE_RULE{GAUSSIAN,LOBATTO,NEWTON_COTES,CLENSHAW_CURTIS};
QUADRATURE_RULE quadrature_rule=GAUSSIAN;
const char *quadrature_rule_name[4]={
  "Gaussian","Lobatto","Newton-Cotes","Clenshaw-Curtis"
};
#elif (SPACEDIM==2)
enum QUADRATURE_RULE
  {GAUSSIAN,LOBATTO,NEWTON_COTES,WANDZURA,DUNAVANT,FEKETE};
QUADRATURE_RULE quadrature_rule=GAUSSIAN;
const char *quadrature_rule_name[6]={
  "Gaussian","Lobatto","Newton-Cotes","Wandzura","Dunavant","Fekete"
};
#else
enum QUADRATURE_RULE{GAUSSIAN,LOBATTO,NEWTON_COTES,KEAST,FELIPPA,
  GRUNDMANN_MOELLER};
QUADRATURE_RULE quadrature_rule=GAUSSIAN;
const char *quadrature_rule_name[6]={
  "Gaussian","Lobatto","Newton-Cotes","Keast","Felippa",
  "Grundmann-Moeller"
};
#endif
int iquadrature_rule=quadrature_rule;

GUI_INPUT_PARAMETER_LIST_TYPE *main_list=0;
BOOLEAN skip_gui=FALSE;
char *display_name=0;
REAL winsize=0.5;

INTEGER npoints=1;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void processCommandLine(int argc,char *argv[],char *display_name,
ifstream &in_file) {
//TRACER_CALL(t,"processCommandLine");
  if (argc<2) skip_gui=FALSE;
  else {
    in_file.open(argv[1],ios::in);
    CHECK_TEST(!in_file.fail());
    if (in_file) {
      INTEGER i=2;
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
void readMainInput(ifstream &in_file) {
//TRACER_CALL(t,"readMainInput");
  BOOLEAN found_main=FALSE;
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
template double __cmath_power<double>(double,unsigned);
template class InputParameter<BOOLEAN>;
template class InputParameter<INTEGER>;
template class InputParameter<REAL>;
#include "NumPtr.C"
INSTANTIATE_NUMPTR(BOOLEAN)
INSTANTIATE_NUMPTR(NumPtr<REAL>)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void makeMainList() {
//TRACER_CALL(t,"makeMainList");
  main_list=OPERATOR_NEW GUI_INPUT_PARAMETER_LIST_TYPE("Main");

  { const char *group="Quadrature Parameters";
#if (SPACEDIM==1)
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
      iquadrature_rule,"quadrature rule",quadrature_rule_name[0],
      quadrature_rule_name[1],quadrature_rule_name[2],
      quadrature_rule_name[3]));
#elif (SPACEDIM==2)
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
      iquadrature_rule,"quadrature rule",quadrature_rule_name[0],
      quadrature_rule_name[1],quadrature_rule_name[2],
      quadrature_rule_name[3],quadrature_rule_name[4],
      quadrature_rule_name[5]));
#else
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
      iquadrature_rule,"quadrature rule",quadrature_rule_name[0],
      quadrature_rule_name[1],quadrature_rule_name[2],
      quadrature_rule_name[3],quadrature_rule_name[4],
      quadrature_rule_name[5]));
#endif
    main_list->append(OPERATOR_NEW GUIInputParameter<INTEGER>(
      npoints,"npoints",1,32,group));
  }

  { const char *group="Graphics";
    main_list->append(OPERATOR_NEW GUIInputParameter<REAL>(winsize,
      "winsize",0.,1.,group));
  }
#ifdef DEBUG
//main_list->printOn(cout);
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void checkMainInput() {
//TRACER_CALL(t,"checkMainInput");
  quadrature_rule=static_cast<QUADRATURE_RULE>(iquadrature_rule);
  if (quadrature_rule==LOBATTO) {
    npoints=max(npoints,2);
  }
  if (quadrature_rule==NEWTON_COTES) {
    npoints=max(min(npoints,8),2);
  }
#if (SPACEDIM==1)
  if (quadrature_rule==CLENSHAW_CURTIS) {
    npoints=max(npoints,2);
  }
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==1)
void runCase(const Quadrature<1> &q) {
//TRACER_CALL(t,"runCase");
#ifdef DEBUG
//cout << "\tnpoints = " << q.numberPoints() << endl;
//cout << "\torder = " << q.order() << endl;
#endif
  const NumPtr< Point<1> > &qp=q.points();
  const NumPtr<REAL> &qw=q.weights();
#ifdef DEBUG
//for (int p=0;p<q.numberPoints();p++) {
//  cout << "\t" << p << " : quadrature_point[" << qp[p] 
//       << "], weight[" << qw[p] << "]" << endl;
//}
#endif
  for (int i=0;i<q.order();i++) {
    REAL sum=0.;
    for (int p=0;p<q.numberPoints();p++) {
      sum+=pow(qp[p][0],i)*qw[p];
    }
    REAL exact=1./static_cast<REAL>(i+1);
    REAL error=(abs(sum-exact)/abs(exact))/DBL_EPSILON;
    if (error>128.) {
      cout << "\t" << i << " : error = " << error*DBL_EPSILON << " "
           << error << endl;
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==2)
void runQuadrilateralCase(const Quadrature<2> &q) {
//TRACER_CALL(t,"runQuadrilateralCase");
#ifdef DEBUG
//cout << "\tnpoints = " << q.numberPoints() << endl;
//cout << "\torder = " << q.order() << endl;
#endif
  const NumPtr< Point<2> > &qp=q.points();
  const NumPtr<REAL> &qw=q.weights();
#ifdef DEBUG
//for (int p=0;p<q.numberPoints();p++) {
//  cout << "\t" << p << " : quadrature_point[" << qp[p] 
//       << "], weight[" << qw[p] << "]" << endl;
//}
#endif
  for (int i=0;i<q.order();i++) {
    for (int j=0;j<q.order();j++) {
      REAL sum=0.;
      for (int p=0;p<q.numberPoints();p++) {
        sum+=pow(qp[p][0],i)*pow(qp[p][1],j)*qw[p];
      }
      REAL exact=1./static_cast<REAL>((i+1)*(j+1));
      REAL error=(abs(sum-exact)/abs(exact))/DBL_EPSILON;
      if (error>128.) {
        cout << "\t" << i << " " << j << " : error = " 
             << error*DBL_EPSILON << " " << error << endl;
      }
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void runHexahedronCase(const Quadrature<3> &q) {
//TRACER_CALL(t,"runHexahedronCase");
#ifdef DEBUG
//cout << "\tnpoints = " << q.numberPoints() << endl;
//cout << "\torder = " << q.order() << endl;
#endif
  const NumPtr< Point<3> > &qp=q.points();
  const NumPtr<REAL> &qw=q.weights();
#ifdef DEBUG
//for (int p=0;p<q.numberPoints();p++) {
//  cout << "\t" << p << " : quadrature_point[" << qp[p] 
//       << "], weight[" << qw[p] << "]" << endl;
//}
#endif
  for (int i=0;i<q.order();i++) {
    for (int j=0;j<q.order();j++) {
      for (int k=0;k<q.order();k++) {
        REAL sum=0.;
        for (int p=0;p<q.numberPoints();p++) {
          sum+=pow(qp[p][0],i)*pow(qp[p][1],j)*pow(qp[p][2],k)*qw[p];
        }
        REAL exact=1./static_cast<REAL>((i+1)*(j+1)*(k+1));
        REAL error=(abs(sum-exact)/abs(exact))/DBL_EPSILON;
        if (error>128.) {
          cout << "\t" << i << " " << j << " " << k << " : error = " 
               << error*DBL_EPSILON << " " << error << endl;
        }
      }
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double factorial(int n) {
  if (n<=1) return 1.;
  return static_cast<double>(n)*factorial(n-1);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==2)
void runTriangleCase(const Quadrature<2> &q) {
//TRACER_CALL(t,"runTriangleCase");
#ifdef DEBUG
//cout << "\tnpoints = " << q.numberPoints() << endl;
//cout << "\torder = " << q.order() << endl;
#endif
  const NumPtr< Point<2> > &qp=q.points();
  const NumPtr<REAL> &qw=q.weights();
#ifdef DEBUG
//for (int p=0;p<q.numberPoints();p++) {
//  cout << "\t" << p << " : quadrature_point[" << qp[p] 
//       << "], weight[" << qw[p] << "]" << endl;
//}
#endif
  for (int i=0;i<q.order();i++) {
    for (int j=0;i+j<q.order();j++) {
      REAL sum=0.;
      for (int p=0;p<q.numberPoints();p++) {
        sum+=pow(qp[p][0],i)*pow(qp[p][1],j)*qw[p];
      }
      REAL exact=factorial(i)*factorial(j)/factorial(2+i+j);
      REAL error=(abs(sum-exact)/exact)/DBL_EPSILON;
      if (error>128.) {
        cout << "\t" << i << " " << j << " : error = " 
           << error*DBL_EPSILON << " " << " " << sum
           << " " << factorial(i)*factorial(j)/factorial(2+i+j)
           << " " << error << endl;
      }
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void runTetrahedronCase(const Quadrature<3> &q) {
//TRACER_CALL(t,"runTetrahedronCase");
#ifdef DEBUG
//cout << "\tnpoints = " << q.numberPoints() << endl;
//cout << "\torder = " << q.order() << endl;
#endif
  const NumPtr< Point<3> > &qp=q.points();
  const NumPtr<REAL> &qw=q.weights();
#ifdef DEBUG
//for (int p=0;p<q.numberPoints();p++) {
//  cout << "\t" << p << " : quadrature_point[" << qp[p] 
//       << "], weight[" << qw[p] << "]" << endl;
//}
#endif
  for (int i=0;i<q.order();i++) {
    for (int j=0;i+j<q.order();j++) {
      for (int k=0;i+j+k<q.order();k++) {
        REAL sum=0.;
        for (int p=0;p<q.numberPoints();p++) {
          sum+=pow(qp[p][0],i)*pow(qp[p][1],j)*pow(qp[p][2],k)*qw[p];
        }
        REAL exact=factorial(i)*factorial(j)*factorial(k)
                  /factorial(3+i+j+k);
        REAL error=(abs(sum-exact)/exact)/DBL_EPSILON;
        if (error>128.) {
          cout << "\t" << i << " " << j << " " << k << " : error = " 
          << error*DBL_EPSILON << " " << error << endl;
        }
      }
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM==3)
void runPrismCase(const Quadrature<3> &q) {
//TRACER_CALL(t,"runPrismCase");
#ifdef DEBUG
//cout << "\tnpoints = " << q.numberPoints() << endl;
//cout << "\torder = " << q.order() << endl;
#endif
  const NumPtr< Point<3> > &qp=q.points();
  const NumPtr<REAL> &qw=q.weights();
#ifdef DEBUG
//for (int p=0;p<q.numberPoints();p++) {
//  cout << "\t" << p << " : quadrature_point[" << qp[p] 
//       << "], weight[" << qw[p] << "]" << endl;
//}
#endif
  for (int i=0;i<q.order();i++) {
    for (int j=0;i+j<q.order();j++) {
      for (int k=0;i+j+k<q.order();k++) {
        REAL sum=0.;
        for (int p=0;p<q.numberPoints();p++) {
          sum+=pow(qp[p][0],i)*pow(qp[p][1],j)*pow(qp[p][2],k)*qw[p];
        }
        REAL exact=factorial(i)*factorial(j)
                  /(factorial(2+i+j)*static_cast<REAL>(k+1));
        REAL error=(abs(sum-exact)/exact)/DBL_EPSILON;
        if (error>128.) {
          cout << "\t" << i << " " << j << " " << k << " : error = " 
            << error*DBL_EPSILON << " " << error << endl;
        }
      }
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void runMain(BOOLEAN /*called_before*/) {
//TRACER_CALL(t,"runMain");
  for (int n=1;n<=10;n++) {
    int m=dimensionalRoot(SPACEDIM,__cmath_power(n,SPACEDIM));
    CHECK_SAME(m,n)
  }
#if (SPACEDIM>=2)
  for (int n=1;n<=10;n++) {
    int m=triangleRoot((n*(n+1))/2);
    CHECK_SAME(m,n)
  }
#endif
#if (SPACEDIM==3)
  for (int n=1;n<=10;n++) {
    int m=tetrahedronRoot((n*(n+1)*(n+2))/6);
    CHECK_SAME(m,n)
  }
#endif
  switch (quadrature_rule) {
    case GAUSSIAN: {
//    TRACER_CALL(t,"runMain GAUSSIAN");
#if (SPACEDIM==1)
      {
        TRACER_CALL(t,"GaussianQuadrature");
        for (int n=1;n<=npoints;n++) {
          cout << "\tn = " << n << endl;
          GaussianQuadrature<1> gq(n);
          int ord=gq.order();
          CHECK_SAME(gq.numberPoints(),
                     GaussianQuadrature<1>::nPointsGivenOrder(ord));
//        cout << "\torder = " << ord << endl;
//        cout << "\targGivenOrder = " << gq.argGivenOrder(ord) << endl;
          runCase(gq);
        }
      }
#elif (SPACEDIM==2)
//
      {
        TRACER_CALL(t,"GaussianQuadrature");
        for (int n=1;n<=npoints;n++) {
          cout << "\tn = " << n << endl;
          GaussianQuadrature<2> gq(n);
          int ord=gq.order();
          CHECK_SAME(gq.numberPoints(),
                     GaussianQuadrature<2>::nPointsGivenOrder(ord));
//        cout << "\torder = " << ord << endl;
//        cout << "\targGivenOrder = " << gq.argGivenOrder(ord) << endl;
          runQuadrilateralCase(gq);
        }
      }
//
//
      {
        TRACER_CALL(t,"TriangleGaussianQuadrature");
//      constructor argument is base_order, not npoints
        for (int n=2;n<=min(npoints,31);n++) {
          cout << "\tn = " << n << endl;
          TriangleGaussianQuadrature tgq(n);
          int ord=tgq.order();
          CHECK_SAME(tgq.numberPoints(),
                     TriangleGaussianQuadrature::nPointsGivenOrder(ord));
//        cout << "\torder = " << ord << endl;
//        cout << "\targGivenOrder = " << tgq.argGivenOrder(ord) 
//             << endl;
          runTriangleCase(tgq);
        }
      }
//
      {
        TRACER_CALL(t,"TriangleProductGaussianQuadrature");
        for (int n=1;n<=npoints;n++) {
          cout << "\tn = " << n << endl;
          TriangleProductGaussianQuadrature tpgq(n);
          int ord=tpgq.order();
          CHECK_SAME(tpgq.numberPoints(),
               TriangleProductGaussianQuadrature::nPointsGivenOrder(ord));
          cout << "\torder = " << ord << endl;
          cout << "\targGivenOrder = " << tpgq.argGivenOrder(ord) << endl;
          runTriangleCase(tpgq);
        }
      }
//
#else
//
      {
        TRACER_CALL(t,"GaussianQuadrature");
        for (int n=1;n<=npoints;n++) {
          cout << "\tn = " << n << endl;
          GaussianQuadrature<3> gq(n);
          int ord=gq.order();
          CHECK_SAME(gq.numberPoints(),
               GaussianQuadrature<3>::nPointsGivenOrder(ord));
          cout << "\torder = " << ord << endl;
          cout << "\targGivenOrder = " << gq.argGivenOrder(ord) << endl;
          runHexahedronCase(gq);
        }
      }
//
//
      {
        TRACER_CALL(t,"TetrahedronGaussianQuadrature");
        for (int n=2;n<=min(npoints,4);n++) {
          cout << "\tn = " << n << endl;
          TetrahedronGaussianQuadrature tgq(n);
          int ord=tgq.order();
          CHECK_SAME(tgq.numberPoints(),
               TetrahedronGaussianQuadrature::nPointsGivenOrder(ord));
          cout << "\torder = " << ord << endl;
          cout << "\targGivenOrder = " << tgq.argGivenOrder(ord) << endl;
          runTetrahedronCase(tgq);
        }
      }
//
//
      {
        TRACER_CALL(t,"TetrahedronProductGaussianQuadrature");
        for (int n=2;n<=npoints;n++) {
          cout << "\tn = " << n << endl;
          TetrahedronProductGaussianQuadrature tpgq(n);
          int ord=tpgq.order();
          CHECK_SAME(tpgq.numberPoints(),
            TetrahedronProductGaussianQuadrature::nPointsGivenOrder(ord));
          cout << "\torder = " << ord << endl;
          cout << "\targGivenOrder = " << tpgq.argGivenOrder(ord) << endl;
          runTetrahedronCase(tpgq);
        }
      }
//
      {
        TRACER_CALL(t,"ProductQuadrature");
        for (int n=1;n<=npoints;n++) {
          cout << "\tn = " << n << endl;
          TriangleProductGaussianQuadrature 
            tpgq2(TriangleProductGaussianQuadrature::argGivenOrder(n));
          GaussianQuadrature<1> 
            gq1(GaussianQuadrature<1>::argGivenOrder(n));
          cout << "\ttpgq order = " << tpgq2.order() << endl;
          cout << "\tgq order = " << gq1.order() << endl;
          ProductQuadrature pq(tpgq2,gq1);
          cout << "\torder = " << pq.order() << endl;
          runPrismCase(pq);
        }
      }
//
#endif
      break;
    }
    case LOBATTO: {
//    TRACER_CALL(t,"runMain LOBATTO");
#if (SPACEDIM==1)
      {
        TRACER_CALL(t,"LobattoQuadrature");
        for (int n=2;n<=npoints;n++) {
          cout << "\tn = " << n << endl;
          LobattoQuadrature<1> lq(n);
          int ord=lq.order();
          CHECK_SAME(lq.numberPoints(),
                     LobattoQuadrature<1>::nPointsGivenOrder(ord));
//        cout << "\torder = " << ord << endl;
//        cout << "\targGivenOrder = " << lq.argGivenOrder(ord) << endl;
          runCase(lq);
        }
      }
#elif (SPACEDIM==2)
//
      {
        TRACER_CALL(t,"LobattoQuadrature");
        for (int n=2;n<=npoints;n++) {
          cout << "\tn = " << n << endl;
          LobattoQuadrature<2> lq(n);
          int ord=lq.order();
          CHECK_SAME(lq.numberPoints(),
                     LobattoQuadrature<2>::nPointsGivenOrder(ord));
          cout << "\torder = " << ord << endl;
          cout << "\targGivenOrder = " << lq.argGivenOrder(ord) << endl;
          runQuadrilateralCase(lq);
        }
      }
//
      {
        TRACER_CALL(t,"TriangleLobattoQuadrature");
        for (int n=2;n<=min(19,npoints);n++) {
          cout << "\tn = " << n << endl;
          TriangleLobattoQuadrature tlq(n);
          int ord=tlq.order();
          CHECK_SAME(tlq.numberPoints(),
                     TriangleLobattoQuadrature::nPointsGivenOrder(ord));
          cout << "\torder = " << ord << endl;
          cout << "\targGivenOrder = " << tlq.argGivenOrder(ord) << endl;
          runTriangleCase(tlq);
        }
      }
#else
//
      {
        TRACER_CALL(t,"LobattoQuadrature");
        for (int n=2;n<=npoints;n++) {
          cout << "\tn = " << n << endl;
          LobattoQuadrature<3> lq(n);
          int ord=lq.order();
          CHECK_SAME(lq.numberPoints(),
                     LobattoQuadrature<3>::nPointsGivenOrder(ord));
          cout << "\torder = " << ord << endl;
          cout << "\targGivenOrder = " << lq.argGivenOrder(ord) << endl;
          runHexahedronCase(lq);
        }
      }
//
//
      if (npoints>=2) {
        TRACER_CALL(t,"TetrahedronLobattoQuadrature");
        for (int n=2;n<=min(2,npoints);n++) {
          cout << "\tn = " << n << endl;
          TetrahedronLobattoQuadrature tlq(n);
          int ord=tlq.order();
          CHECK_SAME(tlq.numberPoints(),
                     TetrahedronLobattoQuadrature::nPointsGivenOrder(ord));
          cout << "\torder = " << ord << endl;
          cout << "\targGivenOrder = " << tlq.argGivenOrder(ord) << endl;
          runTetrahedronCase(tlq);
        }
      }
//
//
      {
        TRACER_CALL(t,"ProductQuadrature");
        for (int n=1;n<=min(19,npoints);n++) {
          cout << "\tn = " << n << endl;
          TriangleLobattoQuadrature 
            tlq2(TriangleLobattoQuadrature::argGivenOrder(n));
          LobattoQuadrature<1> 
            lq1(LobattoQuadrature<1>::argGivenOrder(n));
          ProductQuadrature pq(tlq2,lq1);
          cout << "\torder = " << pq.order() << endl;
          runPrismCase(pq);
        }
      }
//
#endif
      break;
    }
    case NEWTON_COTES: {
//    TRACER_CALL(t,"runMain NEWTON_COTES");
#if (SPACEDIM==1)
      {
        TRACER_CALL(t,"NewtonCotesQuadrature");
        for (int n=2;n<=min(npoints,8);n++) {
          cout << "\tn = " << n << endl;
          NewtonCotesQuadrature<1> ncq(n);
          int ord=ncq.order();
          if (n==NewtonCotesQuadrature<1>::argGivenOrder(ord)) {
            CHECK_SAME(ncq.numberPoints(),
                       NewtonCotesQuadrature<1>::nPointsGivenOrder(ord));
          }
//        cout << "\torder = " << ord << endl;
//        cout << "\targGivenOrder = " << ncq.argGivenOrder(ord) << endl;
          runCase(ncq);
        }
      }
#elif (SPACEDIM==2)
//
      {
        TRACER_CALL(t,"NewtonCotesQuadrature");
        for (int n=2;n<=min(npoints,8);n++) {
          cout << "\tn = " << n << endl;
          NewtonCotesQuadrature<2> ncq(n);
          int ord=ncq.order();
          if (n==ncq.argGivenOrder(ord)) {
            CHECK_SAME(ncq.numberPoints(),
                       NewtonCotesQuadrature<2>::nPointsGivenOrder(ord));
          }
          cout << "\torder = " << ord << endl;
          cout << "\targGivenOrder = " << ncq.argGivenOrder(ord) << endl;
          runQuadrilateralCase(ncq);
        }
      }
//
//
      {
        TRACER_CALL(t,"TriangleNewtonCotesQuadrature");
        for (int n=2;n<=min(npoints,6);n++) {
          cout << "\tn = " << n << endl;
          TriangleNewtonCotesQuadrature tncq(n);
          int ord=tncq.order();
          if (n==tncq.argGivenOrder(ord)) {
            CHECK_SAME(tncq.numberPoints(),
                   TriangleNewtonCotesQuadrature::nPointsGivenOrder(ord));
          }
          cout << "\torder = " << ord << endl;
          cout << "\targGivenOrder = " << tncq.argGivenOrder(ord) << endl;
          runTriangleCase(tncq);
        }
      }
//
      {
        TRACER_CALL(t,"TriangleFullNewtonCotesQuadrature");
        for (int n=2;n<=min(npoints,9);n++) {
          cout << "\tn = " << n << endl;
          TriangleFullNewtonCotesQuadrature tfncq(n);
          int ord=tfncq.order();
          if (n==tfncq.argGivenOrder(ord)) {
            CHECK_SAME(tfncq.numberPoints(),
               TriangleFullNewtonCotesQuadrature::nPointsGivenOrder(ord));
          }
          cout << "\torder = " << ord << endl;
          cout << "\targGivenOrder = " << tfncq.argGivenOrder(ord) << endl;
          runTriangleCase(tfncq);
        }
      }
#else
//
      {
        TRACER_CALL(t,"NewtonCotesQuadrature");
        for (int n=2;n<=min(npoints,8);n++) {
          cout << "\tn = " << n << endl;
          NewtonCotesQuadrature<3> ncq(n);
          int ord=ncq.order();
          if (n==ncq.argGivenOrder(ord)) {
            CHECK_SAME(ncq.numberPoints(),
                       NewtonCotesQuadrature<3>::nPointsGivenOrder(ord));
          }
          cout << "\torder = " << ord << endl;
          cout << "\targGivenOrder = " << ncq.argGivenOrder(ord) << endl;
          runHexahedronCase(ncq);
        }
      }
//
//
      {
        TRACER_CALL(t,"TetrahedronNewtonCotesQuadrature");
        for (int n=2;n<=min(npoints,4);n++) {
          cout << "\tn = " << n << endl;
          TetrahedronNewtonCotesQuadrature tncq(n);
          int ord=tncq.order();
          if (n==tncq.argGivenOrder(ord)) {
            CHECK_SAME(tncq.numberPoints(),
              TetrahedronNewtonCotesQuadrature::nPointsGivenOrder(ord));
          }
          cout << "\torder = " << ord << endl;
          cout << "\targGivenOrder = " << tncq.argGivenOrder(ord) << endl;
          runTetrahedronCase(tncq);
        }
      }
//
//
      {
        TRACER_CALL(t,"TetrahedronFullNewtonCotesQuadrature");
        for (int n=2;n<=min(npoints,7);n++) {
          cout << "\tn = " << n << endl;
          TetrahedronFullNewtonCotesQuadrature tfncq(n);
          int ord=tfncq.order();
          if (n==tfncq.argGivenOrder(ord)) {
            CHECK_SAME(tfncq.numberPoints(),
              TetrahedronFullNewtonCotesQuadrature::nPointsGivenOrder(ord));
          }
          cout << "\torder = " << ord << endl;
          cout << "\targGivenOrder = " << tfncq.argGivenOrder(ord) << endl;
          runTetrahedronCase(tfncq);
        }
      }
//
//
      {
        TRACER_CALL(t,"ProductQuadrature");
        for (int n=2;n<=min(npoints,4);n++) {
          cout << "\tn = " << n << endl;
          TriangleNewtonCotesQuadrature // 2 <= order <= 4 
            tncq2(TriangleNewtonCotesQuadrature::argGivenOrder(n));
          NewtonCotesQuadrature<1> // 2 <= order <= 8
            ncq1(NewtonCotesQuadrature<1>::argGivenOrder(n));
          ProductQuadrature pq(tncq2,ncq1);
//        cout << "\torder = " << pq.order() << endl;
          runPrismCase(pq);
        }
      }
//
#endif
      break;
    }
#if (SPACEDIM==1)
    case CLENSHAW_CURTIS: {
      {
        TRACER_CALL(t,"ClenshawCurtisQuadrature");
        for (int n=2;n<=npoints;n++) {
          cout << "\tn = " << n << endl;
          ClenshawCurtisQuadrature<1> ccq(n);
//        int ord=ccq.order();
//        if (n==ClenshawCurtisQuadrature<1>::argGivenOrder(ord)) {
//          CHECK_SAME(ccq.numberPoints(),
//                     NewtonCotesQuadrature<1>::nPointsGivenOrder(ord));
//        }
//        cout << "\torder = " << ord << endl;
//        cout << "\targGivenOrder = " << ccq.argGivenOrder(ord) << endl;
          runCase(ccq);
        }
      }
      break;
    }
#endif
#if (SPACEDIM==2)
    case WANDZURA: {
      {
        TRACER_CALL(t,"TriangleWandzuraQuadrature");
//      constructor argument is base_order, not npoints
        for (int n=1;n<=min(6,npoints);n++) {
          cout << "\tn = " << n << endl;
          TriangleWandzuraQuadrature twq(n);
          int ord=twq.order();
          CHECK_SAME(twq.numberPoints(),
                 TriangleWandzuraQuadrature::nPointsGivenOrder(ord));
          cout << "\torder = " << ord << endl;
          cout << "\targGivenOrder = " << twq.argGivenOrder(ord) << endl;
          runTriangleCase(twq);
        }
      }
      break;
    }
    case DUNAVANT: {
      {
        TRACER_CALL(t,"TriangleDunavantQuadrature");
//      constructor argument is base_order, not npoints
        for (int n=1;n<=min(19,npoints);n++) {
          cout << "\tn = " << n << endl;
          TriangleDunavantQuadrature tdq(n);
          int ord=tdq.order();
          CHECK_SAME(tdq.numberPoints(),
                 TriangleDunavantQuadrature::nPointsGivenOrder(ord));
          cout << "\torder = " << ord << endl;
          cout << "\targGivenOrder = " << tdq.argGivenOrder(ord) << endl;
          runTriangleCase(tdq);
        }
      }
      break;
    }
    case FEKETE: {
      {
        TRACER_CALL(t,"TriangleFeketeQuadrature");
//      constructor argument is base_order, not npoints
        for (int n=1;n<=min(7,npoints);n++) {
          cout << "\tn = " << n << endl;
          TriangleFeketeQuadrature tfq(n);
          int ord=tfq.order();
          CHECK_SAME(tfq.numberPoints(),
                 TriangleFeketeQuadrature::nPointsGivenOrder(ord));
          cout << "\torder = " << ord << endl;
          cout << "\targGivenOrder = " << tfq.argGivenOrder(ord) << endl;
          runTriangleCase(tfq);
        }
      }
      break;
    }
#endif
#if (SPACEDIM==3)
    case KEAST: {
      {
        TRACER_CALL(t,"TetrahedronKeastQuadrature");
//      constructor argument is base_order, not npoints
        for (int n=1;n<=min(8,npoints);n++) {
          cout << "\tn = " << n << endl;
          TetrahedronKeastQuadrature tkq(n);
          int ord=tkq.order();
          if (n==tkq.argGivenOrder(ord)) {
            CHECK_SAME(tkq.numberPoints(),
                   TetrahedronKeastQuadrature::nPointsGivenOrder(ord));
          }
          cout << "\torder = " << ord << endl;
          cout << "\targGivenOrder = " << tkq.argGivenOrder(ord) << endl;
          runTetrahedronCase(tkq);
        }
      }
      break;
    }
    case FELIPPA: {
      {
        TRACER_CALL(t,"TetrahedronFelippaQuadrature");
//      constructor argument is base_order, not npoints
        for (int n=1;n<=min(9,npoints);n++) {
          cout << "\tn = " << n << endl;
          TetrahedronFelippaQuadrature tfq(n);
          int ord=tfq.order();
          if (n==tfq.argGivenOrder(ord)) {
            CHECK_SAME(tfq.numberPoints(),
                   TetrahedronFelippaQuadrature::nPointsGivenOrder(ord));
          }
          cout << "\torder = " << ord << endl;
          cout << "\targGivenOrder = " << tfq.argGivenOrder(ord) << endl;
          runTetrahedronCase(tfq);
        }
      }
      break;
    }
    case GRUNDMANN_MOELLER: {
      {
        TRACER_CALL(t,"TetrahedronGrundmannMoellerQuadrature");
//      constructor argument is base_order, not npoints
        for (int n=1;n<=min(4,npoints);n++) {
          cout << "\tn = " << n << endl;
          TetrahedronGrundmannMoellerQuadrature tgmq(n);
          int ord=tgmq.order();
          if (n==tgmq.argGivenOrder(ord)) {
            CHECK_SAME(tgmq.numberPoints(),
             TetrahedronGrundmannMoellerQuadrature::nPointsGivenOrder(ord));
          }
          cout << "\torder = " << ord << endl;
          cout << "\targGivenOrder = " << tgmq.argGivenOrder(ord) << endl;
          runTetrahedronCase(tgmq);
        }
      }
      break;
    }
#endif
    default:
      OBSOLETE("unknown quadrature_rule");
  }
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
  processCommandLine(argc,argv,display_name,in_file);
  makeMainList();
  if (in_file) {
    readMainInput(in_file);
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
