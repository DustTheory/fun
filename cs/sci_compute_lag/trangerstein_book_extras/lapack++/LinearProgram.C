#include "LinearProgram.H"
//#include "MyInline.H"
//#include "Tracer.H"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> VirtualLinearProgram<T>::VirtualLinearProgram(
const Matrix<T,T> &Ain,const Vector<T,T> &bin,const Vector<T,T> &cin) : 
A(0),b(0),c(0),current_status(UNKNOWN) {
  int M=Ain.size(0),N=Ain.size(1);
  CHECK_SAME(M,bin.size());
  CHECK_SAME(N,cin.size());
  A=OPERATOR_NEW Matrix<T,T>(M,N);
  A->copyFrom('A',M,N,Ain);
  b=OPERATOR_NEW Vector<T,T>(M);
  b->copyFrom(M,bin);
  c=OPERATOR_NEW Vector<T,T>(N);
  c->copyFrom(N,cin);
  x=OPERATOR_NEW Vector<T,T>(N);
  y=OPERATOR_NEW Vector<T,T>(M);
}

template<class T> T VirtualLinearProgram<T>::currentValue() const {
  switch (current_status) {
    case UNKNOWN:
    case INFEASIBLE:
    case UNBOUNDED:
      return huge_;
    default:
      break;
  }
  if (current_status==PRIMAL_FEASIBLE) return c->dot(*x);
  else return y->dot(*b);
}

template<class T> void VirtualLinearProgram<T>::printOn(ostream &s) 
const {
  s << "VirtualLinearProgram: current_status = " << current_status 
    << " \n";
  s << "A: \n" << endl;
  A->printOn(s);
  s << "b: \n" << endl;
  b->printOn(s);
  s << "c: \n" << endl;
  c->printOn(s);
  s << "x: \n" << endl;
  x->printOn(s);
  s << "y: \n" << endl;
  y->printOn(s);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> SFLinearProgram<T>::SFLinearProgram(
const Matrix<T,T> &Ain,const Vector<T,T> &bin,
const Vector<T,T> &cin) : VirtualLinearProgram<T>(Ain,bin,cin),r(0),
h(0),QR(0),Abasic(0),xbasic(0),cbasic(0),resid(0),column(0) {
//TRACER_CALL(t,"SFLinearProgram<T>::SFLinearProgram");
  int M=this->A->size(0),N=this->A->size(1);
  CHECK_TEST(M<=N);

  for (int i=0;i<M;i++) {
    if ((*this->b)[i]<this->zero_) {
      for (int j=0;j<N;j++) (*this->A)(i,j)=-(*this->A)(i,j);
      (*this->b)[i]=-(*this->b)[i];
    }
  }

  r=OPERATOR_NEW Vector<T,T>(N-M);
  h=OPERATOR_NEW Vector<T,T>(M);
  column=OPERATOR_NEW_BRACKET(int,N);
  for (int j=0;j<N;j++) column[j]=j;
  new_basic=N; //initial column index invalid
#ifdef DEBUG
//printOn(cout);
#endif
}

/*
template<class T> T SFLinearProgram<T>::currentValue() const {
  if (this->current_status==PRIMAL_FEASIBLE ||
  this->current_status==OPTIMAL) {
    return this->c->dot(*this->x);
  } else return this->huge_;
}
*/

template<class T> void SFLinearProgram<T>::specifyBasicVariables(
int *basic_index) {
//TRACER_CALL(t,"SFLinearProgram::specifyBasicVariables");
  CHECK_TEST(this->current_status==UNKNOWN);
  int M=this->A->size(0);
#ifdef DEBUG
  for (int j=1;j<M;j++) {
    CHECK_TEST(basic_index[j-1]<basic_index[j]);
  }
#endif
  Abasic=OPERATOR_NEW SquareMatrix<T,T>(M);
  for (int j=0;j<M;j++) {
    int t=column[j];
    column[j]=basic_index[j];
    column[basic_index[j]]=t;
    memcpy(Abasic->addr(0,j),this->A->addr(0,column[j]),M*sizeof(T));
  }
  QR=OPERATOR_NEW GramSchmidtQRFactorization<T,T>(*Abasic);
  xbasic=OPERATOR_NEW Vector<T,T>(M);
  resid=OPERATOR_NEW Vector<T,T>(M);
  QR->solveOverdetermined(*this->b,*xbasic,*resid);

  *this->x=this->zero_;
  cbasic=OPERATOR_NEW Vector<T,T>(M);
  for (int j=0;j<M;j++) {
    (*this->x)[column[j]]=(*xbasic)[j];
    (*cbasic)[j]=(*this->c)[column[j]];
  }
  QR->solveUnderdetermined(*cbasic,*this->y);

//find cost reduction vector
  int jmin=smallestDualSlack();
  this->current_status=
    ((*this->r)[jmin]<this->zero_ ? PRIMAL_FEASIBLE : OPTIMAL);
  new_basic=jmin+M;
#ifdef DEBUG
//printOn(cout);
#endif
}

template<class T> STATUS_OPTION 
SFLinearProgram<T> ::findBasicFeasibleGuess() {
//TRACER_CALL(t,"SFLinearProgram::findBasicFeasibleGuess");
//if we have already dinked with this program, get out now
  switch (this->current_status) {
    case PRIMAL_FEASIBLE:
    case OPTIMAL:
      return this->current_status;
    case DUAL_FEASIBLE:
    case INFEASIBLE:
    case UNBOUNDED:
      CHECK_TEST(0);
    default:
      break;
  }

//augment the program and find a basic feasible solution
  int M=this->A->size(0),N=this->A->size(1);
  Matrix<T,T> augmented_A(M,M+N,this->zero_);
  Vector<T,T> augmented_c(M+N,this->zero_);
  for (int j=0;j<M;j++) {
    augmented_A(j,j)=this->one_;
    augmented_c[j]=this->one_;
  }
  memcpy(augmented_A.addr(0,M),this->A->addr(),M*N*sizeof(T));
#ifdef DEBUG
//cout << "\taugmented_c = " << endl;
//augmented_c.printOn(cout);
//cout << "\taugmented_A = " << endl;
//augmented_A.printOn(cout);
#endif

//initialize and solve the augmented program
  SFLinearProgram<T> LP(augmented_A,*this->b,augmented_c);
  int *basic_indices=OPERATOR_NEW_BRACKET(int,M);
  for (int i=0;i<M;i++) basic_indices[i]=i;
  LP.specifyBasicVariables(basic_indices);
  delete basic_indices; basic_indices=0;
  while (LP.current_status == PRIMAL_FEASIBLE) LP.simplexStep();
  CHECK_TEST(LP.currentStatus()==OPTIMAL);

  *this->x=this->zero_;
  cbasic=OPERATOR_NEW Vector<T,T>(M);
  for (int j=0;j<M;j++) {
    int lp_colj=LP.column[j];
    if (lp_colj<M) return this->current_status=INFEASIBLE;
    int colj=lp_colj-M;
    column[j]=colj;
    (*this->x)[colj]=(*LP.x)[lp_colj];
    (*cbasic)[j]=(*this->c)[colj];
  }
  int i=M;
  for (int j=M;j<M+N;j++) {
    int lp_colj=LP.column[j];
    if (lp_colj>=M) {
      column[i]=lp_colj-M;
      i=i+1;
    }
  }
  QR=OPERATOR_NEW GramSchmidtQRFactorization<T,T>(*LP.QR);
  QR->solveUnderdetermined(*cbasic,*this->y);

//find cost reduction vector
  int jmin=smallestDualSlack();
  this->current_status=
    ((*r)[jmin]<this->zero_ ? PRIMAL_FEASIBLE : OPTIMAL);
  new_basic=jmin+M;
#ifdef DEBUG
//printOn(cout);
#endif
  return this->current_status;
}

template<class T> STATUS_OPTION SFLinearProgram<T>::simplexStep() {
//TRACER_CALL(t,"SFLinearProgram::simplexStep");
  switch (this->current_status) {
    case OPTIMAL:
      return this->current_status;
    case INFEASIBLE:
    case UNBOUNDED:
      CHECK_TEST(0);
    default:
      break;
  }
  int M=this->A->size(0),N=this->A->size(1);

//find largest feasible simplex step
  int imin=N; 
  int colmin=N;
  T epsmin=this->zero_;
  int col_new_basic=column[new_basic];
#ifdef DEBUG
//cout << "\tM,N = " << M << " " << N << endl;
//cout << "\tnew_basic = " << new_basic << endl;
//cout << "\tcol_new_basic = " << col_new_basic << endl;
//cout << "\tx = " << endl;
//x->printOn(cout);
#endif
  Vector<T,T> An(M);
  memcpy(An.addr(),this->A->addr(0,col_new_basic),M*sizeof(T));
  Vector<T,T> resid(M);
  QR->solveOverdetermined(An,*h,resid);
  for (int i=0;i<M;i++) {
    T hi=(*h)[i];
    int coli=column[i];
    if (hi>this->zero_) {
      T xi=(*this->x)[coli];
      if (imin>=N || xi < epsmin*hi) { 
	imin=i; colmin=coli; epsmin=xi/hi;
      } else if (xi<=this->zero_ && coli<colmin) { 
	imin=i; colmin=coli; epsmin=this->zero_; 
      }
    }
  }
#ifdef DEBUG
//cout << "\th = A_B^inv * A_N * e_jmin" << endl;
//h->printOn(cout);
//cout << "\timin = " << imin << endl;
//cout << "\tepsmin = " << epsmin << endl;
#endif
  if (imin>=N) return this->current_status=UNBOUNDED;
  QR->exchangeColumn(imin,col_new_basic,column[imin],*this->A);
#ifdef DEBUG
//QR->printOn(cout);
#endif

//update solution
//QR->solveOverdetermined(*b,*xbasic); // then permute xbasic to get x
  for (int i=0;i<M;i++) (*this->x)[column[i]] -= (*h)[i] * epsmin;
  (*this->x)[column[imin]]=this->zero_;
  (*this->x)[col_new_basic]=epsmin;
#ifdef DEBUG
//cout << "\tupdated x" << endl;
//x->printOn(cout);
#endif

//update dual solution
  (*cbasic)[imin]=(*this->c)[col_new_basic];
  QR->solveUnderdetermined(*cbasic,*this->y);
#ifdef DEBUG
//cout << "\tupdated y" << endl;
//y->printOn(cout);
#endif

//switch columns
  int col=column[imin];
  column[imin]=column[new_basic];
  column[new_basic]=col;
#ifdef DEBUG
//cout << "\tupdated column = " << endl;
//cout << *column << endl;
#endif

//find cost reduction vector
  int jmin=smallestDualSlack();
  this->current_status=
    ((*r)[jmin]<this->zero_ ? PRIMAL_FEASIBLE : OPTIMAL);
  new_basic=jmin+M;
#ifdef DEBUG
//cout << "\tupdated r= " << endl;
//r->printOn(cout);
//cout << "new_basic = " << column[new_basic] << endl;
//printOn(cout);
#endif
  return this->current_status;
}

template<class T> T SFLinearProgram<T>::costSensitivity(int j,
Vector<T,T> &dydcj) const {
  int M=this->A->size(0);
  CHECK_SAME(M,dydcj.size())
  for (int i=0;i<M;i++) {
    if (column[i]==j) {
      Vector<T,T> axisi(M,this->zero_);
      axisi[i]=this->one_;
      QR->solveUnderdetermined(axisi,dydcj);
      return (*this->x)[j]; // partial objective / partial cj = xj
    }
  }
//nonbasic cost
  dydcj=this->zero_;
  return this->zero_;
}

//if perturbation in single component b_i between lower[i] and upper[i],
//then perturbed optimal x nonnegative
template<class T> void SFLinearProgram<T>::constraintBounds(
Vector<T,T> &lower,Vector<T,T> &upper) const {
  int M=this->A->size(0);
  assert(lower.size()==M);
  assert(upper.size()==M);
  for (int i=0;i<M;i++) {
    T li=-this->huge_;
    T ui=this->huge_;
    Vector<T,T> axisi(M,this->zero_);
    axisi[i]=this->one_;
    Vector<T,T> gi(M);
    QR->solveOverdetermined(axisi,gi,*resid);
    for (int j=0;j<M;j++) {
      T denom=gi[j];
      T xj=(*this->x)[column[j]];
      if (denom>this->zero_ && -xj/denom>li) li=max(li,-xj/denom);
      else if (denom<this->zero_ && -xj/denom<ui) ui=min(ui,-xj/denom);
    }
    lower[i]=li;
    upper[i]=ui;
  }
}

template<class T> T SFLinearProgram<T>::constraintSensitivity(int i,
Vector<T,T> &dxdbi) const {
  int M=this->A->size(0),N=this->A->size(1);
  assert(dxdbi.size()==N);
  Vector<T,T> axisi(M,this->zero_);
  axisi[i]=this->one_;
  QR->solveOverdetermined(axisi,*xbasic,*resid);
  dxdbi=this->zero_;
  for (int j=0;j<M;j++) {
    dxdbi[column[j]]=(*xbasic)[j];
  }
  return (*this->y)[i]; // partial objective / partial bi = yi
}

/*
template<class T> T SFLinearProgram<T>::arraySensitivity(int i,int j,
Vector<T,T> &dxdA,Vector<T,T> &dydA) const {
  int M=this->A->size(0),N=this->A->size(1);
  assert(dxdA.size()==N);
  assert(dydA.size()==M);
//perturbation in basic column of A
  T yi=(*this->y)[i];
  T xj=(*this->x)[j];
  for (int jj=0;jj<M;jj++) {
    if (column[jj]==j) {
      for (int k=0;k<M;k++) {
	dxdA[column[k]]=-(*basic_inverse)(k,i)*xj;
	dydA[k]=-yi*(*basic_inverse)(jj,k);
      }
      return yi*xj;
    }
  }

//perturbation in non-basic column of A
  dxdA=this->zero_;
  dydA=this->zero_;
  return this->zero_;
}
*/

template<class T> void SFLinearProgram<T>::printOn(ostream &s) const {
  VirtualLinearProgram<T>::printOn(s);
  s << "SFLinearProgram: new_basic = " << new_basic << endl;
  if (r != 0) {
    s << "r = " << endl;
    r->printOn(s);
  }
  if (h != 0) {
    s << "h = " << endl;
    h->printOn(s);
  }
  if (Abasic != 0) {
    s << "Abasic = " << endl;
    Abasic->printOn(s);
  }
  if (xbasic != 0) {
    s << "xbasic = " << endl;
    xbasic->printOn(s);
  }
  if (cbasic != 0) {
    s << "cbasic = " << endl;
    cbasic->printOn(s);
  }
  if (resid != 0) {
    s << "resid = " << endl;
    resid->printOn(s);
  }
  if (column != 0 ) {
    s << "column= ";
    int N=this->A->size(1);
    for (int j=0;j<N;j++) s << column[j] << " ";
    s << endl;
  }
  if (QR != 0) {
    s << "QR = " << endl;
    QR->printOn(s);
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> void testSFLinearProgram(T scalar) {
  int *basic_index=0;

/*
//Giapetto problem, p. 48 of Winston 
//x=[0,0,20,20,60]
//y=[-3,-1,-2]
  Matrix<T,T> A(3,5); 
  Vector<T,T> b(3); 
  Vector<T,T> c(5); 
  c[0]=0.;   c[1]=0.;   c[2]=0.;   c[3]=-3.;   c[4]=-2.;
  A(0,0)=1.; A(0,1)=0.; A(0,2)=0.; A(0,3)= 2.; A(0,4)= 1.; b[0]=100.;
  A(1,0)=0.; A(1,1)=1.; A(1,2)=0.; A(1,3)= 1.; A(1,4)= 1.; b[1]= 80.;
  A(2,0)=0.; A(2,1)=0.; A(2,2)=1.; A(2,3)= 1.; A(2,4)= 0.; b[2]= 40.;
  basic_index=OPERATOR_NEW_BRACKET(int,3);
  basic_index[0]=0; basic_index[1]=1; basic_index[2]=2;
*/

/*
//Dorian problem, p. 130 of Winston 
  Matrix<T,T> A(4,7); 
  Vector<T,T> b(4); 
  Vector<T,T> c(7); 
  c[0]=0.;   c[1]=0.;   c[2]=0.;   c[3]=0.;   c[4]=-60.;   c[5]=-30. ;   c[6]=-20. ;
  A(0,0)=1.; A(0,1)=0.; A(0,2)=0.; A(0,3)=0.; A(0,4)=  8.; A(0,5)=  6. ; A(0,6)= 1. ; b[0]=48.;
  A(1,0)=0.; A(1,1)=1.; A(1,2)=0.; A(1,3)=0.; A(1,4)=  4.; A(1,5)=  2. ; A(1,6)= 1.5; b[1]=20.;
  A(2,0)=0.; A(2,1)=0.; A(2,2)=1.; A(2,3)=0.; A(2,4)=  2.; A(2,5)=  1.5; A(2,6)= 0.5; b[2]= 8.;
  A(3,0)=0.; A(3,1)=0.; A(3,2)=0.; A(3,3)=1.; A(3,4)=  0.; A(3,5)=  1. ; A(3,6)= 0. ; b[3]= 5.;
  basic_index=OPERATOR_NEW_BRACKET(int,4);
  basic_index[0]=0; basic_index[1]=1; basic_index[2]=2; basic_index[3]=3;
*/

/*
//first two unknowns not basic, feasible
  Matrix<T,T> A(2,5); 
  Vector<T,T> b(2); 
  Vector<T,T> c(5); 
  c[0]= 5.;   c[1]=0.;   c[2]=0.;   c[3]= 3.;   c[4]=-2.;
  A(0,0)=-6.; A(0,1)=0.; A(0,2)=1.; A(0,3)=-2.; A(0,4)= 2.; b[0]= 6.;
  A(1,0)=-3.; A(1,1)=1.; A(1,2)=0.; A(1,3)= 6.; A(1,4)= 3.; b[1]=15.;
*/

/*
//problem unbounded, p. 146
  Matrix<T,T> A(2,6); 
  Vector<T,T> b(2); 
  Vector<T,T> c(6); 
  c[0]=0.;   c[1]=0.;   c[2]=-36.;   c[3]=-30.;   c[4]= 3.;   c[5]= 4.;
  A(0,0)=1.; A(0,1)=0.; A(0,2)=  1.; A(0,3)=  1.; A(0,4)=-1.; A(0,5)= 0.; b[0]= 5.;
  A(1,0)=0.; A(1,1)=1.; A(1,2)=  6.; A(1,3)=  5.; A(1,4)= 0.; A(1,5)=-1.; b[1]=10.;
  basic_index=OPERATOR_NEW_BRACKET(int,2);
  basic_index[0]=0; basic_index[1]=1;
*/

/*
//basic feasible solution is degenerate, p. 159
  Matrix<T,T> A(2,4);
  Vector<T,T> b(2); 
  Vector<T,T> c(4); 
  c[0]=0.;   c[1]=0.;   c[2]=-5.;   c[3]=-2.;
  A(0,0)=1.; A(0,1)=0.; A(0,2)= 1.; A(0,3)= 1.; b[0]=6.;
  A(1,0)=0.; A(1,1)=1.; A(1,2)= 1.; A(1,3)=-1.; b[1]=0.;
  basic_index=OPERATOR_NEW_BRACKET(int,2);
  basic_index[0]=0; basic_index[1]=1;
*/

//problem cycles without lexico-graphic ordering, p. 161
  Matrix<T,T> A(2,6);
  Vector<T,T> b(2); 
  Vector<T,T> c(6); 
  c[0]=0.;   c[1]=0.;   c[2]=-2.  ;   c[3]=-3.;   c[4]= 1.   ;   c[5]=12.;
  A(0,0)=1.; A(0,1)=0.; A(0,2)=-2.  ; A(0,3)=-9.; A(0,4)= 1.   ; A(0,5)= 9.; b[0]=0.;
  A(1,0)=0.; A(1,1)=1.; A(1,2)=1./3.; A(1,3)= 1.; A(1,4)=-1./3.; A(1,5)=-2.; b[1]=0.;
  basic_index=OPERATOR_NEW_BRACKET(int,2);
  basic_index[0]=0; basic_index[1]=1;
//

  int m=A.size(0),n=A.size(1);
  SFLinearProgram<T> LP(A,b,c);
  if (basic_index!=0) {
    LP.specifyBasicVariables(basic_index);
    delete basic_index; basic_index=0;
  } else {
    LP.findBasicFeasibleGuess();
  }

  while (LP.currentStatus()==PRIMAL_FEASIBLE) {
    LP.simplexStep();
//  cout << "current value = " << LP.currentValue() << endl;
//  cout << "current solution = " << endl;
//  LP.currentPrimalSolution().printOn(cout);
//  char dummy[80];
//  cout << "hit RETURN to continue" << endl;
//  cin.getline(dummy,80);
  }

  cout << "\n\nafter iteration, status = " << LP.currentStatus()
       << endl;

  if (LP.currentStatus()==OPTIMAL) {
    cout << "final value = " << LP.currentValue() << endl;
    cout << "final primal solution = " << endl;
    LP.currentPrimalSolution().printOn(cout);
    cout << "final dual solution = " << endl;
    LP.currentDualSolution().printOn(cout);
    cout << "\tc . x = " << c.dot(LP.currentPrimalSolution()) << endl;
    cout << "\ty . b = " << b.dot(LP.currentDualSolution()) << endl;
  }
//if (LP.currentStatus()==UNBOUNDED) {
//  Vector<T,T> r(n,0.);
//  r.copy(c);
//  A.gemv(-1.,LP.currentDualSolution(),1.,r,'T');
//  cout << "\tc - A^T * y = " << endl;
//  r.printOn(cout);
//}

/*
  Vector<T,T> c_lower(n);
  Vector<T,T> c_upper(n);
  LP.costBounds(c_lower,c_upper);
  cout << "lower bounds on cost vector:" << endl;
  c_lower.printOn(cout);
  cout << "upper bounds on cost vector:" << endl;
  c_upper.printOn(cout);

  Vector<T,T> dydcj(m);
  for (int j=0;j<n;j++) {
    T dzdc=LP.costSensitivity(j,dydcj);
    cout << "\nj,dzdc = " << j << " " << dzdc << endl;
    cout << "dydcj = " << endl;
    dydcj.printOn(cout);
  }

  Vector<T,T> b_lower(m);
  Vector<T,T> b_upper(m);
  LP.constraintBounds(b_lower,b_upper);
  cout << "lower bounds on constraint vector:" << endl;
  b_lower.printOn(cout);
  cout << "upper bounds on constraint vector:" << endl;
  b_upper.printOn(cout);

  Vector<T,T> dxdbi(n);
  for (int i=0;i<m;i++) {
    T dzdb=LP.constraintSensitivity(i,dxdbi);
    cout << "\ni,dzdb = " << i << " " << dzdb << endl;
    cout << "dxdbi = " << endl;
    dxdbi.printOn(cout);
  }
*/

//T A_lower,A_upper;
//cout << "lower bounds on array:" << endl;
//for (int i=0;i<m;i++) {
//  for (int j=0;j<n;j++) {
//    LP.arrayBounds(i,j,A_lower,A_upper);
//    cout << A_lower << " ";
//  }
//  cout << endl;
//}
//cout << "upper bounds on array:" << endl;
//for (int i=0;i<m;i++) {
//  for (int j=0;j<n;j++) {
//    LP.arrayBounds(i,j,A_lower,A_upper);
//    cout << A_upper << " ";
//  }
//  cout << endl;
//}

//for (int i=0;i<m;i++) {
//  for (int j=0;j<n;j++) {
//    T dzda=
//    LP.arraySensitivity(i,j,*primal_derivative,*dual_derivative);
//    cout << "\ni,j,dzda = " << i << " " << j << " " << dzda <<endl;
//    cout << "primal_derivative = " << endl;
//    cout << *primal_derivative << endl;
//    cout << "dual_derivative = " << endl;
//    cout << *dual_derivative << endl;
//  }
//}
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
template<class T> T LinearProgram<T>::currentValue() const {
  switch (this->current_status) {
    case UNKNOWN:
    case INFEASIBLE:
    case UNBOUNDED:
      return this->huge_;
    case PRIMAL_FEASIBLE: {
      T val=this->zero_;
      for (int j=0;j<basic_number;j++) {
	int colj=column[j];
	val += (*this->c)[colj] * (*this->x)[colj];
      }
      return val;
    }
    default: {
      T val=this->zero_;
      for (int i=0;i<basic_number;i++) {
	int rowi=row[i];
	val += (*this->y)[rowi] * (*this->b)[rowi];
      }
      return val;
    }
  }
}
*/

/*
template<class T> int LinearProgram<T>::smallestDualSlack() {
  TRACER_CALL(t,"LinearProgram::smallestDualSlack");
  int N=this->A->size(1);
  if (this->basic_number>=N) return N;
  int jmin=N; 
  T r_min_value=this->huge_;
  for (int j=basic_number;j<N;j++) {
    int colj=column[j];
    T rj=(*this->c)[colj];
    for (int i=0;i<basic_number;i++) {
      int rowi=row[i];
      rj -= (*this->y)[rowi]*(*this->A)(rowi,colj);
    }
#ifdef DEBUG
//  cout << "j,rj = " << j << " " << rj << endl;
#endif
    if (rj < r_min_value) { jmin=j; r_min_value=rj; }
    (*this->r)[j-basic_number]=rj;
  }
  return jmin;
}
*/

/*
template<class T> int LinearProgram<T>::smallestPrimalSlack() {
  TRACER_CALL(t,"LinearProgram::smallestPrimalSlack");
  int M=this->A->size(0);
  if (this->basic_number>=M) return M;
  int imin=M; 
  T s_min_value=this->huge_;
  for (int i=basic_number;i<M;i++) {
    int rowi=row[i];
    T si=-(*this->b)[rowi];
    for (int j=0;j<basic_number;j++) {
      int colj=column[j];
      si += (*this->A)(rowi,colj) * (*this->x)[colj];
    }
#ifdef DEBUG
//  cout << "i,slack = " << i << " " << si << endl;
#endif
    if (si < s_min_value) { imin=i; s_min_value=si; }
    (*this->s)[i-basic_number]=si;
  }
  return imin;
}
*/

template<class T> void LinearProgram<T>::computeSolution() {
//TRACER_CALL(t,"LinearProgram::computeSolution");
  int M=this->A->size(0),N=this->A->size(1);
#ifdef DEBUG
/*
  cout << "basic_number = " << basic_number << endl;
  if (column != 0 ) {
    cout << "column= ";
    for (int j=0;j<N;j++) cout << column[j] << " ";
    cout << endl;
  }
  if (row != 0 ) {
    cout << "row= ";
    for (int i=0;i<M;i++) cout << row[i] << " ";
    cout << endl;
  }
  if (QR!=0) {
    cout << "\tQR = " << endl;
    QR->printOn(cout);
    {
      cout << "\tQ * R = " << endl;
      Matrix<T,T> *prod=QR->orthogonalPart() * QR->rightTrapezoidalPart();
      prod->printOn(cout);
      delete prod; prod=0;
    }
    {
      SquareMatrix<T,T> Abasic(basic_number);
      for (int j=0;j<basic_number;j++) {
        for (int i=0;i<basic_number;i++) {
          Abasic(i,j)=(*this->A)(row[i],column[j]);
        }
      }
      cout << "\tAbasic = " << endl;
      Abasic.printOn(cout);
    }
  }
*/
#endif
  if (bbasic!=0) delete bbasic; bbasic=0;
  if (cbasic!=0) delete cbasic; cbasic=0;
  if (xbasic!=0) delete xbasic; xbasic=0;
  if (ybasic!=0) delete ybasic; ybasic=0;
  if (r!=0) delete r; r=0;
  if (s!=0) delete s; s=0;

  *this->x=this->zero_;
  *this->y=this->zero_;

  if (basic_number>0) {
    bbasic=OPERATOR_NEW Vector<T,T>(basic_number);
    cbasic=OPERATOR_NEW Vector<T,T>(basic_number);
    xbasic=OPERATOR_NEW Vector<T,T>(basic_number);
    ybasic=OPERATOR_NEW Vector<T,T>(basic_number);
    Vector<T,T> resid(basic_number);

    for (int i=0;i<basic_number;i++) {
      (*bbasic)[i]=(*this->b)[row[i]];
    }
    QR->solveOverdetermined(*bbasic,*xbasic,resid);
    for (int j=0;j<basic_number;j++) {
      (*cbasic)[j]=(*this->c)[column[j]];
    }
    QR->solveUnderdetermined(*cbasic,*ybasic);

    for (int j=0;j<basic_number;j++) (*this->x)[column[j]]=(*xbasic)[j];
    for (int i=0;i<basic_number;i++) (*this->y)[row[i]]=(*ybasic)[i];
  }

  if (basic_number<M) {
    s=OPERATOR_NEW Vector<T,T>(M-basic_number);
    for (int i=basic_number;i<M;i++) {
      int rowi=row[i];
      T si=-(*this->b)[rowi];
      for (int j=0;j<basic_number;j++) {
        si+=(*this->A)(rowi,column[j])*(*xbasic)[j];
      }
      (*s)[i-basic_number]=si;
    }
  }
  if (basic_number<N) {
    r=OPERATOR_NEW Vector<T,T>(N-basic_number);
    for (int j=basic_number;j<N;j++) {
      int colj=column[j];
      T rj=(*this->c)[colj];
      for (int i=0;i<basic_number;i++) {
        rj-=(*ybasic)[i]*(*this->A)(row[i],colj);
      }
      (*r)[j-basic_number]=rj;
    }
  }
#ifdef DEBUG
//printOn(cout);
#endif
}

/*
template<class T> T LinearProgram<T>::costSensitivity(int j,
Vector<T,T> &dydcj) const {
  TRACER_CALL(t,"LinearProgram::costSensitivity");
  int M=this->A->size(0);
  CHECK_SAME(M,dydcj.size())
  for (int i=0;i<M;i++) {
    if (column[i]==j) {
      Vector<T,T> axisi(M,this->zero_);
      axisi[i]=this->one_;
      QR->solveUnderdetermined(axisi,dydcj);
      return (*this->x)[j]; // partial objective / partial cj = xj
    }
  }
//nonbasic cost
  dydcj=this->zero_;
  return this->zero_;
}
*/

/*
template<class T> void LinearProgram<T>::constraintBounds(
Vector<T,T> &lower,Vector<T,T> &upper) const {
  TRACER_CALL(t,"LinearProgram::constraintBounds");
  CHECK_TEST(0);
//int M=this->A->size(0);
//assert(lower.size()==M);
//assert(upper.size()==M);

//perturbation in basic variable
//Vector<T,T> s(M-basic_number);
//{ for (int i=basic_number;i<M;i++) {
//  int rowi=row[i];
//  T si=-this->b[rowi];
//  for (int j=0;j<basic_number;j++) {
//    int colj=column[j];
//    si += (*this->A)(rowi,colj)*(*this->x)[colj];
//  }
//  s[i-basic_number]=si;
//} }
//for (int j=0;j<basic_number;j++) {
//  T lo=-this->huge_;
//  T hi=this->huge_;
//  { for (int i=0;i<basic_number;i++) {
//    int coli=column[i];
//    T denom=(*basic_inverse)(i,j);
//    T xi=(*this->x)[coli];
//    if (denom>this->zero_ && lo<-xi/denom) lo=-xi/denom;
//    else if (denom<this->zero_ && hi>-xi/denom) hi=-xi/denom;
//  } }
//  { for (int i=basic_number;i<M;i++) {
//    int rowi=row[i];
//    T denom=this->zero_;
//    for (int k=0;k<basic_number;k++) {
//	int colk=column[k];
//      denom+=(*this->A)(rowi,colk) * (*basic_inverse)(k,j);
//    }
//    T si=s[i-basic_number];
//    if (denom<this->zero_ && lo<si/denom) lo=si/denom;
//    else if (denom>this->zero_ && hi>si/denom) hi=si/denom;
//  } }
//  int rowj=row[j];
//  lower[rowj]=this->b[rowj]+lo;
//  upper[rowj]=this->b[rowj]+hi;
//}
//perturbation in non-basic variable
//{ for (int i=basic_number;i<M;i++) {
//  int rowi=row[i];
//  upper[rowi]=this->b[rowi]+s[i-basic_number];
//  lower[rowi]=-this->huge_;
//} }
}

template<class T> T LinearProgram<T>::constraintSensitivity(int i,
Vector<T,T> &dxdbi) const {
  TRACER_CALL(t,"LinearProgram::constraintSensitivity");
  int M=this->A->size(0),N=this->A->size(1);
  assert(dxdbi.size()==N);
  Vector<T,T> axisi(M,this->zero_);
  axisi[i]=this->one_;
  CHECK_TEST(0);
//QR->solveOverdetermined(axisi,*xbasic,*resid);
  dxdbi=this->zero_;
  for (int j=0;j<M;j++) {
    CHECK_TEST(0);
//  dxdbi[column[j]]=(*xbasic)[j];
  }
  return (*this->y)[i]; // partial objective / partial bi = yi
}
*/

template<class T> STATUS_OPTION LinearProgram<T>::simplexStep() {
//TRACER_CALL(t,"LinearProgram::simplexStep");
  switch (this->current_status) {
    case PRIMAL_FEASIBLE:
      return primalSimplexStep();
    case DUAL_FEASIBLE:
      return dualSimplexStep();
    case INFEASIBLE:
    case UNBOUNDED:
      assert(0);
    default:
      break;
  }
  return this->current_status;
}

template<class T> int LinearProgram<T>::basicPrimalPivot(
const Vector<T,T> &gh,T &min_ratio) {
//TRACER_CALL(t,"LinearProgram::basicPrimalPivot");
#ifdef DEBUG
//cout << "\tx_basic:" << endl;
//xbasic->printOn(cout);
//cout << "\tgh:" << endl;
//gh.printOn(cout);
#endif
  int j_basic=basic_number;
  min_ratio=this->huge_;
  for (int j=0;j<basic_number;j++) {
    T ghj=gh[j];
    if (ghj>this->zero_) {
      T xj=(*xbasic)[j];
      if (xj<ghj*min_ratio) {
        j_basic=j;
        min_ratio=xj/ghj;
      }
    }
  }
#ifdef DEBUG
//cout << "\tj_basic, min_ratio = " << j_basic << " " << min_ratio
//     << endl;
#endif
  return j_basic;
}

template<class T> int LinearProgram<T>::nonBasicPrimalPivot(
const Vector<T,T> &gh,T &min_ratio) {
//TRACER_CALL(t,"LinearProgram::nonBasicPrimalPivot");
  int M=this->A->size(0);
#ifdef DEBUG
//cout << "\ts_nonbasic:" << endl;
//s->printOn(cout);
//cout << "\tgh:" << endl;
//gh.printOn(cout);
#endif
  int i_non_basic=M;
  min_ratio=this->huge_;
  for (int i=0;i<M-basic_number;i++) {
    T ghi=gh[i];
    if (ghi>this->zero_) {
      T primal_slack=(*this->s)[i];
      if (primal_slack<min_ratio*ghi) {
        i_non_basic=i;
        min_ratio=primal_slack/ghi;
      }
    }
  }
#ifdef DEBUG
//cout << "\ti_non_basic, min_ratio = " << i_non_basic << " "
//     << min_ratio << endl;
#endif
  return i_non_basic+basic_number;
}

template<class T> int LinearProgram<T>::basicDualPivot(
Vector<T,T> &gh,T &min_ratio) {
//TRACER_CALL(t,"LinearProgram::basicDualPivot");
#ifdef DEBUG
//cout << "\tybasic = " << endl;
//ybasic->printOn(cout);
//cout << "\tgh = " << endl;
//gh.printOn(cout);
#endif
  int i_basic=basic_number;
  min_ratio=this->huge_;
  for (int i=0;i<basic_number;i++) {
    T ghi=gh[i];
    if (ghi>this->zero_) {
      T yi=(*ybasic)[i];
      if (yi<ghi*min_ratio) {
        i_basic=i;
        min_ratio=yi/ghi;
      }
    }
  }
#ifdef DEBUG
//cout << "\ti_basic,min_ratio = " << i_basic << " " << min_ratio
//     << endl;
#endif
  return i_basic;
}

template<class T> int LinearProgram<T>::nonBasicDualPivot(
const Vector<T,T> &gh,T &min_ratio) {
//TRACER_CALL(t,"LinearProgram::nonBasicDualPivot");
#ifdef DEBUG
//cout << "\tr = " << endl;
//r->printOn(cout);
//cout << "\tgh :\n" << endl;
//gh.printOn(cout);
#endif
  int N=this->A->size(1);
  int j_non_basic=N;
  min_ratio=this->huge_;
  for (int j=0;j<N-basic_number;j++) {
    T ghj=gh[j];
    T dual_slack=(*this->r)[j];
    if (ghj>this->zero_) {
      if (dual_slack<ghj*min_ratio) {
        j_non_basic=j;
        min_ratio=dual_slack/ghj;
      }
    }
  }
#ifdef DEBUG
//cout << "\tj_non_basic,min_ratio = " << j_non_basic
//     << " " << min_ratio << endl;
#endif
  return j_non_basic+basic_number;
}

template<class T> void LinearProgram<T>::drop(int i_basic,int j_basic) {
//TRACER_CALL(t,"LinearProgram::drop");
#ifdef DEBUG
//cout << "\ti_basic,j_basic = " << i_basic << " " << j_basic << endl;
//cout << "dropping row " << row[i_basic] << endl;
//cout << "dropping column " << column[j_basic] << endl;
//printOn(cout);
#endif
  QR->dropRowAndColumn(i_basic,j_basic);
  int jdrop=column[j_basic];
  for (int j=j_basic+1;j<basic_number;j++) column[j-1]=column[j];
  column[basic_number-1]=jdrop;
  int idrop=row[i_basic];
  for (int i=i_basic+1;i<basic_number;i++) row[i-1]=row[i];
  row[basic_number-1]=idrop;
  basic_number--;
#ifdef DEBUG
//cout << "\tcolumn = ";
//for (int jj=0;jj<this->A->size(1);jj++) cout << column[jj] << " ";
//cout << endl;
//cout << "\trow = ";
//for (int ii=0;ii<this->A->size(0);ii++) cout << row[ii] << " ";
//cout << endl;
#endif
}

template<class T> void LinearProgram<T>::add(int i_non_basic,
int j_non_basic) {
//TRACER_CALL(t,"LinearProgram::add");
#ifdef DEBUG
//cout << "\ti_non_basic,j_non_basic = " << i_non_basic << " "
//     << j_non_basic << endl;
//cout << "adding original row " << row[i_non_basic] << endl;
//cout << "adding original column " << column[j_non_basic] << endl;
//if (QR!=0) {
//  QR->printOn(cout);
//  {
//    cout << "\tQ * R = " << endl;
//    Matrix<T,T> *prod=QR->orthogonalPart() * QR->rightTrapezoidalPart();
//    prod->printOn(cout);
//    delete prod; prod=0;
//  }
//  {
//    SquareMatrix<T,T> Abasic(basic_number);
//    for (int j=0;j<basic_number;j++) {
//      for (int i=0;i<basic_number;i++) {
//        Abasic(i,j)=(*this->A)(row[i],column[j]);
//      }
//    }
//    cout << "\tAbasic = " << endl;
//    Abasic.printOn(cout);
//  }
//}
#endif
  if (QR==0) {
    CHECK_TEST(basic_number==0);
    SquareMatrix<T,T> Abasic(1);
    Abasic(0,0)=(*this->A)(row[i_non_basic],column[j_non_basic]);
    QR=OPERATOR_NEW GramSchmidtQRFactorization<T,T>(Abasic);
    int iadd=row[i_non_basic];
    row[i_non_basic]=row[0];
    row[0]=iadd;
    int jadd=column[j_non_basic];
    column[j_non_basic]=column[0];
    column[0]=jadd;
  } else {
    Matrix<T,T> Arow(1,basic_number);
    int iadd=row[i_non_basic];
    for (int j=0;j<basic_number;j++) {
      Arow(0,j)=(*this->A)(iadd,column[j]);
    }
    QR->addRow(0,Arow);
    row[i_non_basic]=row[basic_number];
    row[basic_number]=iadd;

    Matrix<T,T> Acol(basic_number+1,1);
    int jadd=column[j_non_basic];
    for (int i=0;i<basic_number;i++) {
      Acol(i,0)=(*this->A)(row[i],jadd);
    }
    Acol(basic_number,0)=(*this->A)(iadd,jadd);
    QR->addColumn(0,Acol);
    column[j_non_basic]=column[basic_number];
    column[basic_number]=jadd;
  }
  basic_number++;
#ifdef DEBUG
//printOn(cout);
//{
//  cout << "\tQ * R = " << endl;
//  Matrix<T,T> *prod=QR->orthogonalPart() * QR->rightTrapezoidalPart();
//  prod->printOn(cout);
//  delete prod; prod=0;
//}
//{
//  SquareMatrix<T,T> Abasic(basic_number);
//  for (int j=0;j<basic_number;j++) {
//    for (int i=0;i<basic_number;i++) {
//      Abasic(i,j)=(*this->A)(row[i],column[j]);
//    }
//  }
//  cout << "\tAbasic = " << endl;
//  Abasic.printOn(cout);
//}
#endif
}

template<class T> void  LinearProgram<T>::switchRows(int i_basic,
int i_non_basic) {
//TRACER_CALL(t,"LinearProgram::switchRows");
#ifdef DEBUG
//cout << "\ti_basic,i_non_basic = " << i_basic << " " << i_non_basic
//     << endl;
//cout << "switching basic row " << row[i_basic] 
//     << " with non-basic row " << row[i_non_basic] << endl;
#endif
  int rowib=row[i_basic],rowin=row[i_non_basic];
  Matrix<T,T> Abasic(2,basic_number);
  for (int j=0;j<basic_number;j++) {
    int colj=column[j];
    Abasic(0,j)=(*this->A)(rowib,colj);
    Abasic(1,j)=(*this->A)(rowin,colj);
  }
  QR->exchangeRow(i_basic,1,0,Abasic);
  int irow=row[i_basic];
  row[i_basic]=row[i_non_basic];
  row[i_non_basic]=irow;
}

template<class T> void  LinearProgram<T>::switchColumns(int j_basic,
int j_non_basic) {
//TRACER_CALL(t,"LinearProgram::switchColumns");
#ifdef DEBUG
//cout << "\tj_basic,j_non_basic = " << j_basic << " " << j_non_basic
//     << endl;
//cout << "switching basic column " << column[j_basic] 
//     << " with non-basic column " << column[j_non_basic] << endl;
#endif
  int coljb=column[j_basic],coljn=column[j_non_basic];
  Matrix<T,T> Abasic(basic_number,2);
  for (int i=0;i<basic_number;i++) {
    int rowi=row[i];
    Abasic(i,0)=(*this->A)(rowi,coljb);
    Abasic(i,1)=(*this->A)(rowi,coljn);
  }
  QR->exchangeColumn(j_basic,1,0,Abasic);
  int jcol=column[j_basic];
  column[j_basic]=column[j_non_basic];
  column[j_non_basic]=jcol;
}

/*
template<class T> void LinearProgram<T>::costBounds(
Vector<T,T> &lower,Vector<T,T> &upper) const {
  TRACER_CALL(t,"LinearProgram::costBounds");
  CHECK_TEST(0);
#ifdef DEBUG
//printOn(cout);
#endif
//int N=this->A->size(1);
//assert(lower.size(1)==N);
//assert(upper.size(1)==N);

//perturbation in basic variable
//Vector<T,T> r(N-basic_number);
//for (int j=basic_number;j<N;j++) {
//  int colj=column[j];
//  T rj=this->c[colj];
//  for (int i=0;i<basic_number;i++) {
//    int rowi=row[i];
//    rj -= (*this->y)[rowi]*(*this->A)(rowi,colj);
//  }
//  r[j-basic_number]=rj;
//}
#ifdef DEBUG
//cout << "\tr= " << r<< endl;
#endif
//for (int i=0;i<basic_number;i++) {
#ifdef DEBUG
//  cout << "\tperturbation in basic cost i = " << (*column)[i] << endl;
#endif
//  T lo=-this->huge_;
//  T hi=this->huge_;
//  { for (int j=0;j<basic_number;j++) {
//    int rowj=row[j];
//    T denom=(*basic_inverse)(i,j);
//    T yj=(*this->y)[rowj];
//    if (denom>this->zero_ && lo<-yj/denom) lo=-yj/denom;
//    else if (denom<this->zero_ && hi>-yj/denom) hi=-yj/denom;
#ifdef DEBUG
//    cout << "\tbasic j,ratio = " << rowj << " " << -yj/denom 
//         << "\tyj,denom = " << yj << " " << denom 
//         << "\n\tlo,hi = " << lo << " " << hi << endl;
#endif
//  } }
//  { for (int j=basic_number;j<N;j++) {
//    int colj=column[j];
//    T denom=this->zero_;
//    for (int k=0;k<basic_number;k++) {
//	int rowk=row[k];
//      denom+=(*basic_inverse)(i,k) * (*this->A)(rowk,colj);
//    }
//    T rj=r(0,j-basic_number);
//    if (denom<this->zero_ && lo<rj/denom) lo=rj/denom;
//    else if (denom>this->zero_ && hi>rj/denom) hi=rj/denom;
#ifdef DEBUG
//    cout << "\tnon-basic j,ratio = " << colj << " " <<  rj/denom 
//         << "\trj,denom = " << rj << " " << denom 
//         << "\n\tlo,hi = " << lo << " " << hi << endl;
#endif
//  } }
//  int coli=column[i];
//  lower[coli]=this->c[coli]+lo;
//  upper[coli]=this->c[coli]+hi;
//}
//perturbation in non-basic variable
//{ for (int j=basic_number;j<N;j++) {
//  int colj=column[j];
//  lower[colj]=this->c[colj]-r[j-basic_number];
//  upper[colj]=this->huge_;
//} }
#ifdef DEBUG
//cout << "\tlower bounds = " << lower
//     << "\n\tupper bounds = " << upper<< endl;
#endif
}
*/

template<class T> void LinearProgram<T>::printOn(ostream &os) const {
  VirtualLinearProgram<T>::printOn(os);
  os << "LinearProgram:" << endl;
  os << "basic_number = " << basic_number << endl;
  os << "currentValue = " << this->currentValue() << endl;
//
  if (r != 0) {
    os << "r = " << endl;
    r->printOn(os);
  }
  if (s != 0) {
    os << "s = " << endl;
    s->printOn(os);
  }
//
  if (xbasic != 0) {
    os << "xbasic = " << endl;
    xbasic->printOn(os);
  }
  if (ybasic != 0) {
    os << "ybasic = " << endl;
    ybasic->printOn(os);
  }
//
  if (bbasic != 0) {
    os << "bbasic = " << endl;
    bbasic->printOn(os);
  }
  if (cbasic != 0) {
    os << "cbasic = " << endl;
    cbasic->printOn(os);
  }
//
  if (column != 0 ) {
    os << "column= ";
    int N=this->A->size(1);
    for (int j=0;j<N;j++) os << column[j] << " ";
    os << endl;
  }
  if (row != 0 ) {
    os << "row= ";
    int M=this->A->size(0);
    for (int i=0;i<M;i++) os << row[i] << " ";
    os << endl;
  }
//
  if (QR != 0) {
    os << "QR = " << endl;
    QR->printOn(os);
  }
//
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> void testLinearProgram(T scalar) {
  scalar=1.;
  int basic_no=0;
  int *basic_row_index=0;
  int *basic_column_index=0;

/*
//problem 6, Winston p. 286
//solves dual
//adds one variable
//y=[1,0]
  Matrix<T,T> A(2,2); 
  Vector<T,T> b(2); 
  Vector<T,T> c(2); 
  c[0]=1.;   c[1]= 2.;
  A(0,0)=1.; A(0,1)=-1.; b[0]= 3.;
  A(1,0)=1.; A(1,1)= 1.; b[1]= 1.;
//basic_no=1;
//basic_row_index=OPERATOR_NEW_BRACKET(int,1);
//basic_column_index=OPERATOR_NEW_BRACKET(int,1);
//basic_row_index[0]=0;
//basic_column_index[0]=0;
*/

/*
//problem 6, Winston p. 286
//solves primal
//adds one variable
//x=[1,0]
  Matrix<T,T> A(2,2); 
  Vector<T,T> b(2); 
  Vector<T,T> c(2); 
  c[0]=-3.;   c[1]=-1.;
  A(0,0)=-1.; A(0,1)=-1.; b[0]=-1.;
  A(1,0)= 1.; A(1,1)=-1.; b[1]=-2.;
*/

/*
//Leather Limited problem, Winston p. 121
//solves dual
//adds row 0 and col 1; adds row 1 and col 0
//y_trans = [20,20]
  Matrix<T,T> A(2,2); 
  Vector<T,T> b(2); 
  Vector<T,T> c(2); 
  c[0]=40.;   c[1]=60.;
  A(0,0)= 1.; A(0,1)= 2.; b[0]= 4.;
  A(1,0)= 1.; A(1,1)= 1.; b[1]= 3.;
*/

/*
//Leather Limited problem, Winston p. 121
//solves primal
//adds row 1 and col 0; adds row 0 and col 1
//x = [20,20]
  Matrix<T,T> A(2,2); 
  Vector<T,T> b(2); 
  Vector<T,T> c(2); 
  c[0]=-4.;   c[1]=-3.;
  A(0,0)=-1.; A(0,1)=-1.; b[0]=-40.;
  A(1,0)=-2.; A(1,1)=-1.; b[1]=-60.;
*/

/*
//Giapetto problem, Winston p. 45
//solves dual
//adds row 0 & col 2; adds row 1 and col 0; switch cols 1 and 2
//y_trans = [20,60]
  Matrix<T,T> A(2,3); 
  Vector<T,T> b(2); 
  Vector<T,T> c(3); 
  c[0]=100.;   c[1]=80.;   c[2]=40.;
  A(0,0)=  2.; A(0,1)= 1.; A(0,2)= 1.; b[0]= 3.;
  A(1,0)=  1.; A(1,1)= 1.; A(1,2)= 0.; b[1]= 2.;
*/

/*
//Giapetto problem, Winston p. 45
//solves primal
//adds row 2 & col 0; adds row 0 and col 1; switch rows 1 and 2
//x = [20,60]
  Matrix<T,T> A(3,2); 
  Vector<T,T> b(3); 
  Vector<T,T> c(2); 
  c[0]=-3.;   c[1]=-2.;
  A(0,0)=-2.; A(0,1)=-1.; b[0]=-100.;
  A(1,0)=-1.; A(1,1)=-1.; b[1]= -80.;
  A(2,0)=-1.; A(2,1)= 0.; b[2]= -40.;
*/

/*
//Dorian problem, Winston p. 58
//solves dual
//adds row 0 & col 0; adds row 1 and col 1
//y_trans = [5,7.5]
  Matrix<T,T> A(2,2); 
  Vector<T,T> b(2); 
  Vector<T,T> c(2); 
  c[0]=50.;   c[1]=100.;
  A(0,0)= 7.; A(0,1)=  2.; b[0]= 28.;
  A(1,0)= 2.; A(1,1)= 12.; b[1]= 24.;
*/

/*
//Auto company, Winston p. 61
//multiple solutions
//solves dual
//adds row 0 and col 0
//y_trans = [40,0]
  Matrix<T,T> A(2,2); 
  Vector<T,T> b(2); 
  Vector<T,T> c(2); 
  c[0]=120.;   c[1]=50.;
  A(0,0)=  3.; A(0,1)= 1.; b[0]= 3.;
  A(1,0)=  2.; A(1,1)= 1.; b[1]= 2.;
*/

/*
//Problem 2, Winston p. 198
//solves primal
//adds row 0 and col 0
//x = [40,0]
  Matrix<T,T> A(2,2); 
  Vector<T,T> b(2); 
  Vector<T,T> c(2); 
  c[0]=-4.;   c[1]= 1.;
  A(0,0)=-3.; A(0,1)=-1.; b[0]=-6.;
  A(1,0)= 1.; A(1,1)=-2.; b[1]= 0.;
*/

/*
//SOMETHING HANGS
//Auto company, Winston p. 63
//primal infeasible, dual unbounded
//solves dual
//adds row 0 and col 0
  Matrix<T,T> A(2,4); 
  Vector<T,T> b(2); 
  Vector<T,T> c(4); 
  c[0]=120.;   c[1]=50.;   c[2]=-30.;   c[3]=-20.;
  A(0,0)=  3.; A(0,1)= 1.; A(0,2)=  1.; A(0,3)=  0.; b[0]= 3.;
  A(1,0)=  2.; A(1,1)= 1.; A(1,2)=  0.; A(1,3)=  1.; b[1]= 2.;
  basic_no=2;
  basic_row_index=OPERATOR_NEW_BRACKET(int,basic_no);
  basic_row_index[0]=0; basic_row_index[1]=1;
  basic_column_index=OPERATOR_NEW_BRACKET(int,basic_no);
  basic_column_index[0]=2; basic_column_index[1]=3;
*/

/*
//Winston p. 64
//unbounded
  Matrix<T,T> A(2,2); 
  Vector<T,T> b(2); 
  Vector<T,T> c(2); 
  c[0]= 1.;   c[1]=-6.;
  A(0,0)= 1.; A(0,1)= 2.; b[0]= 2.;
  A(1,0)=-1.; A(1,1)= 1.; b[1]=-1.;
  basic_no=2;
  basic_row_index=OPERATOR_NEW_BRACKET(int,basic_no);
  basic_row_index[0]=0; basic_row_index[1]=1;
  basic_column_index=OPERATOR_NEW_BRACKET(int,basic_no);
  basic_column_index[0]=0; basic_column_index[1]=1;
*/

/*
//Dakotta problem, Winston p. 280
//solves dual
//adds row 0 and col 2; adds row2 and col 1
//y_trans = [2,0,8]
  Matrix<T,T> A(3,3); 
  Vector<T,T> b(3); 
  Vector<T,T> c(3); 
  c[0]= 48.;   c[1]=20. ;   c[2]= 8. ;
  A(0,0)=  8.; A(0,1)= 4. ; A(0,2)= 2. ; b[0]= 60.;
  A(1,0)=  6.; A(1,1)= 2. ; A(1,2)= 1.5; b[1]= 30.;
  A(2,0)=  1.; A(2,1)= 1.5; A(2,2)= 0.5; b[2]= 20.;
*/

/*
//Diet problem, Winston p. 68
//solves dual
//adds row 0 and col 1; adds row 2 and col 2; switches row 0 with row 1
//y_trans=[0,2.5,7.5,0]
  Matrix<T,T> A(4,4); 
  Vector<T,T> b(4); 
  Vector<T,T> c(4); 
  c[0]= 50.;   c[1]= 20.;   c[2]= 30.;   c[3]= 80.;
  A(0,0)=400.; A(0,1)=200.; A(0,2)=150.; A(0,3)=500.; b[0]= 500.;
  A(1,0)=  3.; A(1,1)=  2.; A(1,2)=  0.; A(1,3)=  0.; b[1]=   6.;
  A(2,0)=  2.; A(2,1)=  2.; A(2,2)=  4.; A(2,3)=  4.; b[2]=  10.;
  A(3,0)=  2.; A(3,1)=  4.; A(3,2)=  1.; A(3,3)=  5.; b[3]=   8.;
*/

//
//Post Office, Winston p. 71
//solves dual
//adds row 3 and col 0
//adds row 5 and col 1
//adds row 0 and col 3
//adds row 2 and col 5
//adds row 6 and col 2
//y_trans = [1/3,0,1/3,1/3,0,1/3,0]
  Matrix<T,T> A(7,7); 
  Vector<T,T> b(7); 
  Vector<T,T> c(7); 
  for (int j=0;j<7;j++) {
    c[j]=1.;
    for (int k=0;k<5;k++) A((j+k)%7,j)=1.;
    A((j+5)%7,j)=A((j+6)%7,j)=0.;
  }
  b[0]=17.; b[1]=13.; b[2]=15.; b[3]=19.; b[4]=14.; b[5]=16.; b[6]=11.;
//

/*
//Chvatal p 135a
//x=[0,2,0,3,0,2,0,1,1,0,0,0]
//y=[0,6,0,15,2,1,1,0]
  Matrix<T,T> A(8,12,0.);
  Vector<T,T> b(8);
  Vector<T,T> c(12);
  A(0,0)= 1.; A(0,1)=-1.; A(0,2)= 0.; A(0,3)= 0.; A(0, 4)=1.; b[0]=-3;
  A(1,0)= 0.; A(1,1)= 0.; A(1,2)= 1.; A(1,3)=-1.; A(1, 5)=1.; b[1]=-1;
  A(2,0)= 3.; A(2,1)=-3.; A(2,2)=-2.; A(2,3)= 2.; A(2, 6)=1.; b[2]=-1;
  A(3,0)= 1.; A(3,1)=-1.; A(3,2)=-1.; A(3,3)= 1.; A(3, 7)=1.; b[3]= 2;
  A(4,0)=-5.; A(4,1)= 5.; A(4,2)= 4.; A(4,3)=-4.; A(4, 8)=1.; b[4]=-1;
  A(5,0)=-2.; A(5,1)= 2.; A(5,2)= 1.; A(5,3)=-1.; A(5, 9)=1.; b[5]= 1;
  A(6,0)= 4.; A(6,1)=-4.; A(6,2)=-3.; A(6,3)= 3.; A(6,10)=1.; b[6]= 1;
  A(7,0)=-6.; A(7,1)= 6.; A(7,2)= 5.; A(7,3)=-5.; A(7,11)=1.; b[7]=-4;
  c[0]=7;     c[1]=-7;    c[2]=-3;    c[3]=3;
  c[4]=8;     c[5]=6;     c[6]=4;     c[7]=15;
  c[8]=2;     c[9]=10;    c[10]=10;   c[11]=3;
*/

/*
//Chvatal p 135a
//x=[0,6,0,15,2,1,1,0]
//y=[0,2,0,3,0,2,0,1,1,0,0,0]
  Matrix<T,T> A(12,8,0.);
  Vector<T,T> c(8);
  Vector<T,T> b(12);
  A(0,0)=-1.; A(1,0)= 1.; A(2,0)= 0.; A(3,0)= 0.; A( 4,0)=-1.; c[0]= 3;
  A(0,1)= 0.; A(1,1)= 0.; A(2,1)=-1.; A(3,1)= 1.; A( 5,1)=-1.; c[1]= 1;
  A(0,2)=-3.; A(1,2)= 3.; A(2,2)= 2.; A(3,2)=-2.; A( 6,2)=-1.; c[2]= 1;
  A(0,3)=-1.; A(1,3)= 1.; A(2,3)= 1.; A(3,3)=-1.; A( 7,3)=-1.; c[3]=-2;
  A(0,4)= 5.; A(1,4)=-5.; A(2,4)=-4.; A(3,4)= 4.; A( 8,4)=-1.; c[4]= 1;
  A(0,5)= 2.; A(1,5)=-2.; A(2,5)=-1.; A(3,5)= 1.; A( 9,5)=-1.; c[5]=-1;
  A(0,6)=-4.; A(1,6)= 4.; A(2,6)= 3.; A(3,6)=-3.; A(10,6)=-1.; c[6]=-1;
  A(0,7)= 6.; A(1,7)=-6.; A(2,7)=-5.; A(3,7)= 5.; A(11,7)=-1.; c[7]= 4;
  b[0]=-7;    b[1]=  7;   b[2]=   3;  b[3]= -3;
  b[4]=-8;    b[5]= -6;   b[6]=  -4;  b[7]=-15;
  b[8]=-2;    b[9]=-10;   b[10]=-10;  b[11]=-3;
*/

//int m=A.size(0),n=A.size(1);
  LinearProgram<T> LP(A,b,c);
  if (LP.currentStatus()==UNKNOWN) {
    if (basic_row_index!=0 && basic_column_index!=0) {
      LP.specifyBasicVariables(basic_no,basic_row_index,
        basic_column_index);
    } else {
      LP.findBasicFeasibleGuess();
    }
  }
  if (basic_row_index!=0) delete basic_row_index; basic_row_index=0;
  if (basic_column_index!=0) delete basic_column_index;
    basic_column_index=0;

  while (LP.currentStatus()==PRIMAL_FEASIBLE || 
  LP.currentStatus()==DUAL_FEASIBLE) {
    LP.simplexStep();
/*
//  cout << LP << endl;
    cout << "\nstatus = " << LP.currentStatus()
         << endl;
    cout << "current primal value = " << LP.currentValue() <<endl;
    cout << "current primal solution = " << endl;
    LP.currentPrimalSolution().printOn(cout); 
    cout << "current dual solution = " << endl;
    LP.currentDualSolution().printOn(cout);
    char dummy[80];
    cout << "hit RETURN to continue" << endl;
    cin.getline(dummy,80);
*/
  }

  cout << "\n\nafter iteration, status = " << LP.currentStatus()
       << endl;
  cout << "final primal value = " << LP.currentValue() << endl;
  cout << "final primal solution = " << endl;
  LP.currentPrimalSolution().printOn(cout);
  cout << "final dual solution = " << endl;
  LP.currentDualSolution().printOn(cout);
}
