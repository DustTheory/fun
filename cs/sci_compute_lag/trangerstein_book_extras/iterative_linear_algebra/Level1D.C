#include <cmath>
#include <string.h>
#include "Level1D.H"
#include "Fort1D.H"
#include "MemoryDebugger.H"
#include "Tracer.H"

int Level::smoother_iterations=1;
Level::SMOOTHER Level::smoother=Level::GAUSS_SEIDEL_SMOOTHER;
Level::RESTRICTION_PROLONGATION Level::restriction_prolongation=
  Level::ALGEBRAIC_MULTIGRID_RESTRICTION_PROLONGATION;
const char* Level::smoother_name[3]={
  "gauss_seidel_smoother","gauss_seidel_red_black_smoother",
  "richardson_smoother"};
const char* Level::restriction_prolongation_name[2]={
  "algebraic_multigrid_restriction_prolongation",
   "finite_element_restriction_prolongation"};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//constructs the hierarchy of Level's, beginning with the finest
Level::Level(int nc,const double *m,int si,SMOOTHER s,
RESTRICTION_PROLONGATION rp) : coarser(0),finer(0),
ncells(nc),level_number(0),fi(0),la(ncells),stride(ncells+1),
ifirst(1),ilast(ncells-1), /*decay(d),*/ /*error(e),*/
matrix_copy_c(0) {
//TRACER_CALL(t,"Level::Level");
#ifdef DEBUG
//cout << "\tlevel_number = " << level_number << endl;
//cout << "\tncells = " << ncells << endl;
//cout << "\tfi,la,stride = " << fi << " " << la << " " << stride
//     << endl;
//cout << "\tifirst,ilast = " << ifirst << " " << ilast << endl;
#endif
  smoother_iterations=si;
  smoother=s;
  restriction_prolongation=rp;
  if (nc%2==0 && nc>2) coarser=OPERATOR_NEW Level(this);
  else {
//  matrix bandwidth assumed to be 3
    matrix_copy_c=OPERATOR_NEW_BRACKET(double,3*stride);
    matrix_copy=matrix_copy_c-fi;
  }
  setup(m);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//called from public Level::Level
Level::Level(Level *f) : coarser(0),finer(f),ncells(f->ncells/2),
level_number(f->level_number+1),fi(0),la(ncells),stride(ncells+1),
ifirst(1),ilast(ncells-1), /*error(0),*/ matrix_copy_c(0)
{
//TRACER_CALL(t,"Level::Level");
  ASSERT(finer->ncells%2==0);
//decay=finer->decay/(1.+4.*finer->decay);
  if (ncells%2==0 && ncells>2) coarser=OPERATOR_NEW Level(this);
  else {
//  matrix bandwidth assumed to be 3
    matrix_copy_c=OPERATOR_NEW_BRACKET(double,3*stride);
    matrix_copy=matrix_copy_c-fi;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Level::~Level() {
//TRACER_CALL(t,"Level::~Level");
  if (coarser!=0) delete coarser;
  delete [] matrix_c;
  delete [] Ax_b_c;
  delete [] first_residual_c;
  delete [] solution_increment_c;
  delete [] prolongation_vector_c;
  delete [] richardson_residual_c;
  if (matrix_copy_c!=0) delete [] matrix_copy_c;
  coarser=0;
  finer=0;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//called from public Level::Level after all Level's are constructed
void Level::setup(const double *m) {
//TRACER_CALL(t,"Level::setup");
#ifdef DEBUG
//cout << "\tlevel_number = " << level_number << endl;
#endif

//allocate space for matrix, residuals, solution increment
//  matrix bandwidth assumed to be 3
  matrix_c=OPERATOR_NEW_BRACKET(double,3*stride);
  Ax_b_c=OPERATOR_NEW_BRACKET(double,stride);
  first_residual_c=OPERATOR_NEW_BRACKET(double,stride);
  solution_increment_c=OPERATOR_NEW_BRACKET(double,stride);
  prolongation_vector_c=OPERATOR_NEW_BRACKET(double,stride);
  richardson_residual_c=OPERATOR_NEW_BRACKET(double,stride);

//offset memory to simplify addressing
  matrix=matrix_c-fi;
  Ax_b=Ax_b_c-fi;
  first_residual=first_residual_c-fi;
  solution_increment=solution_increment_c-fi;
  prolongation_vector=prolongation_vector_c-fi;
  richardson_residual=richardson_residual_c-fi;

  int stride2=2*stride;
  if (finer==0) {
    memcpy(matrix,m,3*stride*sizeof(double));
  } else {
//  coarser level: compute matrix from finer
    for (int i=finer->ifirst;i<=finer->ilast;i++) {
      finer->first_residual[i]=0.;
    }
    for (int START=ifirst;START<=min(3,ilast);START++) {
//    cout << "\n\tSTART = " << START << endl;
      for (int I=ifirst;I<=ilast;I++) solution_increment[I]=0.;
      for (int i=finer->ifirst;i<=finer->ilast;i++) {
        finer->solution_increment[i]=0.;
      }

      for (int I=START;I<=ilast;I+=3) solution_increment[I]=1.;
//    prolong solution_increment on this level
//      to solution_increment on finer level 
      finer->prolongation(); 
//    compute matrix * solution_increment on finer level
      finer->computeResidual(finer->solution_increment,
        finer->first_residual,finer->Ax_b);
//    restrict matrix * solution_increment from finer level to this
      finer->restriction();
//    store column of coarsened matrix in this level
      for (int I=START;I<=ilast;I+=3) {
        if (I>ifirst) matrix[I-1+stride2]=first_residual[I-1];
        matrix[I+stride]=first_residual[I];
        if (I<ilast) matrix[I+1]=first_residual[I+1];
      }
    }
  }
#ifdef DEBUG
//for (int i=ifirst;i<=ilast;i++) {
//  cout << "\tmatrix[" << i << "] = " << matrix[i] << " " 
//       << matrix[i+stride] << " " << matrix[i+stride*2] << endl;
//}
#endif
  if (smoother==RICHARDSON_SMOOTHER) {
//  Gerschgorin circle theorem to get Richardson mu
    mu=max(abs(matrix[ifirst+stride])+abs(matrix[ifirst+stride2]),
           abs(matrix[ilast])+abs(matrix[ilast+stride]));
    for (int i=ifirst+1;i<ilast;i++) {
      double sum=abs(matrix[i])+abs(matrix[i+stride])
                +abs(matrix[i+stride2]);
      mu=max(mu,sum);
    }
#ifdef DEBUG
//  cout << "\tmu = " << mu << endl;
#endif
  }

  for (int i=fi;i<=la;i++) solution_increment[i]=0.;
//recurse
  if (coarser==0) {
    memcpy(matrix_copy,matrix,3*stride*sizeof(double));
//  gaussian factorization:
    for (int i=ifirst+1;i<=ilast;i++) {
      matrix_copy[i]/=matrix_copy[i-1+stride];
      matrix_copy[i+stride]-=matrix_copy[i]*matrix_copy[i-1+2*stride];
    }
#ifdef DEBUG
//  for (int i=ifirst;i<=ilast;i++) {
//    cout << "\tmatrix_copy[" << i << "] = " << matrix_copy[i] << " " 
//         << matrix_copy[i+stride] << " " << matrix_copy[i+stride*2] 
//         << endl;
//  }
#endif
  } else coarser->setup(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int Level::numberLevels() const {
  if (coarser==0) return level_number+1;
  else return coarser->numberLevels();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//applies smoother to residual to obtain solution_increment
void Level::preSmooth() const {
//TRACER_CALL(t,"Level::preSmooth");
//int m=n-1;
  for (int it=0;it<smoother_iterations;it++) {
    switch (smoother) {
      case GAUSS_SEIDEL_SMOOTHER:
        F77_NAME(gauss_seidel)(fi,la,ifirst,ilast,
          matrix,first_residual, solution_increment);
        break;
      case GAUSS_SEIDEL_RED_BLACK_SMOOTHER:
        F77_NAME(gauss_seidel_red_black)(fi,la,ifirst,ilast, 0,
          matrix,first_residual, solution_increment);
        break;
      case RICHARDSON_SMOOTHER:
      default:
        F77NAME(richardson)(fi,la,ifirst,ilast, mu,
          matrix,first_residual, solution_increment, richardson_residual);
        break;
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//transpose of preSmooth
void Level::postSmooth() const {
//TRACER_CALL(t,"Level::postSmooth");
//ASSERT(n%2==0);
//int m=n-1;
  for (int it=0;it<smoother_iterations;it++) {
    switch (smoother) {
      case GAUSS_SEIDEL_SMOOTHER:
        F77_NAME(gauss_seidel_reverse)(fi,la,ifirst,ilast,
          matrix,first_residual, solution_increment);
        break;
      case GAUSS_SEIDEL_RED_BLACK_SMOOTHER:
        F77_NAME(gauss_seidel_red_black)(fi,la,ifirst,ilast, 1,
          matrix,first_residual, solution_increment);
        break;
      case RICHARDSON_SMOOTHER:
      default:
        F77NAME(richardson)(fi,la,ifirst,ilast, mu,
          matrix,first_residual, solution_increment, richardson_residual);
        break;
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//applies restriction to compute coarser residual from finer
void Level::restriction() const {
  ASSERT(coarser!=0);
//TRACER_CALL(t,"Level::restriction");
#ifdef DEBUG
//for (int i=ifirst;i<=ilast;i++) {
//  cout << "\tAx_b[" << i << "] = " << Ax_b[i] << endl;
//}
#endif
  switch (restriction_prolongation) {
    case ALGEBRAIC_MULTIGRID_RESTRICTION_PROLONGATION:
      F77_NAME(algebraic_multigrid_restriction)(coarser->fi,coarser->la,
        fi,la,coarser->ifirst,coarser->ilast, matrix,Ax_b, 
        coarser->first_residual);
      break;
    case FINITE_ELEMENT_RESTRICTION_PROLONGATION:
      F77_NAME(finite_element_restriction)(coarser->fi,coarser->la,
        fi,la,coarser->ifirst,coarser->ilast, Ax_b, 
        coarser->first_residual);
      break;
  }
#ifdef DEBUG
//for (int I=coarser->ifirst;I<=coarser->ilast;I++) {
//  cout << "\tfirst_residual[" << I << "] = "
//       << coarser->first_residual[I] << endl;
//}
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//applies prolongation to compute finer solution_increment from coarser
//must be adjoint of Level::restriction
void Level::prolongation() const {
  ASSERT(coarser!=0);
//TRACER_CALL(t,"Level::prolongation");
#ifdef DEBUG
//for (int I=coarser->ifirst;I<=coarser->ilast;I++) {
//  cout << "\tsolution_increment[" << I << "] = "
//       << coarser->solution_increment[I] << endl;
//}
#endif
//should prolong into prolongation_vector, then add to solution_increment
  switch (restriction_prolongation) {
    case ALGEBRAIC_MULTIGRID_RESTRICTION_PROLONGATION:
      F77_NAME(algebraic_multigrid_prolongation)(
        coarser->fi,coarser->la,fi,la,coarser->ifirst,coarser->ilast, 
        matrix,coarser->solution_increment, prolongation_vector);
      break;
    case FINITE_ELEMENT_RESTRICTION_PROLONGATION:
      F77_NAME(finite_element_prolongation)(coarser->fi,coarser->la,
        fi,la,coarser->ifirst,coarser->ilast, 
        coarser->solution_increment, prolongation_vector);
      break;
  }
  for (int i=ifirst;i<=ilast;i++) {
    solution_increment[i] += prolongation_vector[i];
  }
#ifdef DEBUG
//for (int i=ifirst;i<=ilast;i++) {
//  cout << "\tsolution_increment[" << i << "] = "
//       << solution_increment[i] << endl;
//}
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::computeResidual(const double *x,const double *b,double *ax_b) 
const {
//TRACER_CALL(t,"Level::computeResidual");
  F77_NAME(compute_residual)(fi,la,ifirst,ilast, matrix,b,x, ax_b);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::updateResidual() const {
//TRACER_CALL(t,"Level::updateResidual");
  F77_NAME(matrix_multiply)(fi,la,ifirst,ilast, 
    matrix,solution_increment, Ax_b);
  for (int i=ifirst;i<=ilast;i++) {
    Ax_b[i]=first_residual[i]-Ax_b[i];
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//the following is usually called on the coarsest level only
void Level::solve() const {
//TRACER_CALL(t,"Level::solve");
//forward solve:
  solution_increment[ifirst]=first_residual[ifirst];
  for (int i=ifirst+1;i<=ilast;i++) {
    solution_increment[i]=first_residual[i]
                         -matrix_copy[i]*solution_increment[i-1];
  }
//back-solve:
  solution_increment[ilast]/=matrix_copy[ilast+stride];
  int stride2=stride*2;
  for (int i=ilast-1;i>=ifirst;i--) {
    solution_increment[i]=(solution_increment[i]
      -matrix_copy[i+stride2]*solution_increment[i+1])
      /matrix_copy[i+stride];
  }
  for (int i=ifirst;i<=ilast;i++) Ax_b[i]=0.;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//one step of the multigrid algorithm
void Level::multigridStep(double *ax_b, double *d) const {
//TRACER_CALL(t,"Level::multigridStep");
#ifdef DEBUG
//cout << "\tlevel_number = " << level_number << endl;
#endif
  if (finer==0) {
    memcpy(&first_residual[ifirst],&ax_b[ifirst],
           (ilast-ifirst+1)*sizeof(double));
  }
  for (int i=ifirst;i<=ilast;i++) solution_increment[i]=0;
  if (coarser==0) {
    solve();
  } else {
    preSmooth();
    updateResidual();
    restriction();
    coarser->multigridStep(coarser->Ax_b,coarser->solution_increment);
    prolongation();
    updateResidual();
    postSmooth();
  }
  if (finer==0) {
    memcpy(&d[ifirst],&solution_increment[ifirst],
           (ilast-ifirst+1)*sizeof(double));
  }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::checkProlongationAndRestrictionRandom() const {
  TRACER_CALL(t,"Level::checkProlongationAndRestrictionRandom");
  for (int i=ifirst;i<=ilast;i++) {
    Ax_b[i]=drand48();
    solution_increment[i]=0.;
  }
  for (int i=coarser->ifirst;i<=coarser->ilast;i++) {
    coarser->solution_increment[i]=drand48();
    coarser->first_residual[i]=0.;
  }
  restriction(); // fine Ax_b -> coarse first_residual
  prolongation(); // coarse solution_inc -> fine solution_inc
//fine_sum=(finer->Ax_b) . (Prolongation * coarser->solution_increment)
//coarse_sum=(Restriction * finer->Ax_b) . (coarser->solution_increment)
  double fine_sum=F77NAME(ddot)(ilast-ifirst+1,&Ax_b[ifirst],1,
                                &solution_increment[ifirst],1);
  double coarse_sum=F77NAME(ddot)(coarser->ilast-coarser->ifirst+1,
    &coarser->first_residual[coarser->ifirst],1,
    &coarser->solution_increment[coarser->ifirst],1);
  cout << "\tcoarse_sum,fine_sum = " << coarse_sum << " "
       << fine_sum << " "
       << abs(coarse_sum-fine_sum)
         /max(abs(coarse_sum),abs(fine_sum)) << endl;
  ASSERT(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::checkProlongationAndRestriction() const {
  TRACER_CALL(t,"Level::checkProlongationAndRestriction");
#ifdef DEBUG
//for (int i=ifirst;i<=ilast;i++) {
//  cout << "\tmatrix[" << i << "] = " << matrix[i] << " " 
//       << matrix[i+stride] << " " << matrix[i+stride*2] << endl;
//}
#endif
  for (int i=ifirst;i<=ilast;i++) 
//int i=6;
  {
    for (int j=ifirst;j<=ilast;j++) Ax_b[j]=0.;
    Ax_b[i]=1.;
    for (int J=coarser->ifirst;J<=coarser->ilast;J++) {
      coarser->first_residual[J]=0.;
    }
    restriction();
#ifdef DEBUG
//  for (int j=ifirst;j<=ilast;j++) {
//    cout << "\tAx_b[" << j << "] = " << Ax_b[j] << endl;
//  }
//  for (int J=coarser->ifirst;J<=coarser->ilast;J++) {
//    cout << "\tfirst_residual[" << J << "] = " << first_residual[J] 
//         << endl;
//  }
#endif
    for (int I=coarser->ifirst;I<=coarser->ilast;I++) 
//  int I=4;
    {
      for (int J=coarser->ifirst;J<=coarser->ilast;J++) {
        coarser->solution_increment[J]=0.;
      }
      coarser->solution_increment[I]=1.;
      for (int j=ifirst;j<=ilast;j++) solution_increment[j]=0.;
      prolongation();
#ifdef DEBUG
//    for (int J=coarser->ifirst;J<=coarser->ilast;J++) {
//      cout << "\tsolution_increment[" << J << "] = " 
//           << solution_increment[J] << endl;
//    }
//    for (int j=ifirst;j<=ilast;j++) {
//      cout << "\tsolution_increment[" << j << "] = " 
//           << solution_increment[j] << endl;
//    }
#endif
      if (abs(coarser->first_residual[I]-solution_increment[i])>
      1.e-12*max(abs(coarser->first_residual[I]),
      abs(solution_increment[i]))) {
        cout << "\tR,P[" << I << "," << i << "] = " 
             << coarser->first_residual[I] << " " 
             << solution_increment[i] << " "
             << abs(coarser->first_residual[I]-solution_increment[i])
               /max(abs(coarser->first_residual[I]),
                    abs(solution_increment[i]))
             << endl;
      }
    }
  }
  ASSERT(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::checkMatrixSymmetryRandom() const {
  TRACER_CALL(t,"Level::checkMatrixSymmetryRandom");
  double *b=OPERATOR_NEW_BRACKET(double,stride);
  double *x=OPERATOR_NEW_BRACKET(double,stride);
  double *y=OPERATOR_NEW_BRACKET(double,stride);
  double *ax=OPERATOR_NEW_BRACKET(double,stride);
  double *ay=OPERATOR_NEW_BRACKET(double,stride);
  for (int i=ifirst;i<=ilast;i++) {
    b[i]=0.;
    x[i]=drand48();
    y[i]=drand48();
  }
  computeResidual(x,b,ax); 
  computeResidual(y,b,ay); 
  double sum1=
    F77NAME(ddot)(ilast-ifirst+1,&y[ifirst],1,&ax[ifirst],1);
  double sum2=
    F77NAME(ddot)(ilast-ifirst+1,&ay[ifirst],1,&x[ifirst],1);
  cout << "\tsum1,sum2 = " << sum1 << " " << sum2 << " "
       << abs(sum1-sum2)/max(abs(sum1),abs(sum2)) << endl;
  delete [] b;
  delete [] x;
  delete [] y;
  delete [] ax;
  delete [] ay;
  ASSERT(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::checkMatrixSymmetry() const {
  TRACER_CALL(t,"Level::checkMatrixSymmetry");
  double *b=OPERATOR_NEW_BRACKET(double,stride);
  double *x=OPERATOR_NEW_BRACKET(double,stride);
  double *y=OPERATOR_NEW_BRACKET(double,stride);
  double *ax=OPERATOR_NEW_BRACKET(double,stride);
  double *ay=OPERATOR_NEW_BRACKET(double,stride);
  for (int i=ifirst;i<=ilast;i++) {
    b[i]=0.;
  }
  for (int j=ifirst;j<=ilast;j++) 
//int j=16;
  {
    for (int k=ifirst;k<=ilast;k++) y[k]=0;
    y[j]=1.;
    computeResidual(y,b,ay); 
    for (int i=ifirst;i<=ilast;i++) 
//  int i=17;
    {
      for (int k=ifirst;k<=ilast;k++) x[k]=0;
      x[i]=1.;
      computeResidual(x,b,ax); 
      if (abs(ay[i]-ax[j])>1.e-12*max(abs(ay[i]),abs(ax[j]))) {
        cout << "\tA[" << i << "," << j << "] = "
             << ay[i] << " " << ax[j] << endl;
      }
    }
  }
  delete [] b;
  delete [] x;
  delete [] y;
  delete [] ax;
  delete [] ay;
  ASSERT(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::checkVCycleSymmetryRandom() const {
  TRACER_CALL(t,"Level::checkVCycleSymmetryRandom");
  double *r=OPERATOR_NEW_BRACKET(double,stride);
  double *s=OPERATOR_NEW_BRACKET(double,stride);
  double *vr=OPERATOR_NEW_BRACKET(double,stride);
  double *vs=OPERATOR_NEW_BRACKET(double,stride);
  for (int i=ifirst;i<=ilast;i++) {
    r[i]=drand48();
    s[i]=drand48();
  }
#ifdef DEBUG
  for (int i=ifirst;i<=ilast;i++) {
    cout << "\tr[" << i << "] = " << r[i] << endl;
  }
  for (int i=ifirst;i<=ilast;i++) {
    cout << "\ts[" << i << "] = " << s[i] << endl;
  }
#endif

  memcpy(&first_residual[ifirst],&r[ifirst],
         (ilast-ifirst+1)*sizeof(double));
  for (int i=ifirst;i<=ilast;i++) solution_increment[i]=0;
  if (coarser==0) {
    solve();
  } else {
    preSmooth();
    updateResidual();
    restriction();
    coarser->multigridStep(coarser->Ax_b,coarser->solution_increment);
    prolongation();
    updateResidual();
    postSmooth();
  }
  memcpy(&vr[ifirst],&solution_increment[ifirst],
         (ilast-ifirst+1)*sizeof(double));

  memcpy(&first_residual[ifirst],&s[ifirst],
         (ilast-ifirst+1)*sizeof(double));
  for (int i=ifirst;i<ilast;i++) solution_increment[i]=0;
  if (coarser==0) {
    solve();
  } else {
    preSmooth();
    updateResidual();
    restriction();
    coarser->multigridStep(coarser->Ax_b,coarser->solution_increment);
    prolongation();
    updateResidual();
    postSmooth();
  }
  memcpy(&vs[ifirst],&solution_increment[ifirst],
         (ilast-ifirst+1)*sizeof(double));

#ifdef DEBUG
  for (int i=ifirst;i<=ilast;i++) {
    cout << "\tvr[" << i << "] = " << vr[i] << endl;
  }
  for (int i=ifirst;i<=ilast;i++) {
    cout << "\tvs[" << i << "] = " << vs[i] << endl;
  }
#endif
  double sum1=
    F77NAME(ddot)(ilast-ifirst+1,&s[ifirst],1,&vr[ifirst],1);
  double sum2=
    F77NAME(ddot)(ilast-ifirst+1,&vs[ifirst],1,&r[ifirst],1);
  cout << "\tsum1,sum2 = " << sum1 << " " << sum2 << " "
       << abs(sum1-sum2)/max(abs(sum1),abs(sum2)) << endl;
  delete [] r;
  delete [] s;
  delete [] vr;
  delete [] vs;
  ASSERT(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::checkVCycleSymmetry() const {
  TRACER_CALL(t,"Level::checkVCycleSymmetry");
  double *r=OPERATOR_NEW_BRACKET(double,stride);
  double *s=OPERATOR_NEW_BRACKET(double,stride);
  double *vr=OPERATOR_NEW_BRACKET(double,stride);
  double *vs=OPERATOR_NEW_BRACKET(double,stride);
  for (int j=ifirst;j<=ilast;j++) {
    for (int k=ifirst;k<=ilast;k++) s[k]=0;
    s[j]=1.;
    multigridStep(s,vs);
    for (int i=ifirst;i<=ilast;i++) {
      for (int k=ifirst;k<=ilast;k++) r[k]=0;
      r[i]=1.;
      multigridStep(r,vr);
      if (abs(vs[i]-vr[j])>1.e-12*max(abs(vs[i]),abs(vr[j]))) {
        cout << "\tV[" << i << "," << j << "] = "
             << vs[i] << " " << vr[j] << endl;
      }
    }
  }
  delete [] r;
  delete [] s;
  delete [] vr;
  delete [] vs;
  ASSERT(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::checkSmootherRandom() const {
  TRACER_CALL(t,"Level::checkSmootherRandom");
  double *fr=OPERATOR_NEW_BRACKET(double,stride);
  for (int i=ifirst;i<=ilast;i++) { 
    first_residual[i]=drand48();
    solution_increment[i]=0.;
  }
  memcpy(&fr[ifirst],&first_residual[ifirst],
         (ilast-ifirst+1)*sizeof(double));
  preSmooth(); // first_residual -> solution_increment

  for (int i=ifirst;i<=ilast;i++) first_residual[i]=drand48();
  double sum1=F77NAME(ddot)(ilast-ifirst+1,&first_residual[ifirst],1,
                                &solution_increment[ifirst],1);
  for (int i=ifirst;i<=ilast;i++) solution_increment[i]=0.;

  postSmooth(); // first_residual -> solution_increment
  double sum2=F77NAME(ddot)(ilast-ifirst+1,&fr[ifirst],1,
                                &solution_increment[ifirst],1);
  cout << "\tsum1,sum2 = " << sum1 << " " << sum2 << " " 
       << abs(sum2-sum1) << endl;
  ASSERT(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::checkSmoother() const {
  TRACER_CALL(t,"Level::checkSmoother");
  double *si=OPERATOR_NEW_BRACKET(double,stride);
  for (int j=ifirst;j<=ilast;j++) 
//int j=6;
  {
    for (int k=ifirst;k<=ilast;k++) {
      first_residual[k]=0.;
      solution_increment[k]=0.;
    }
    first_residual[j]=1.;
    preSmooth();
#ifdef DEBUG
    for (int i=ifirst;i<=ilast;i++) {
      cout << "\tfirst_residual[" << i << "] = " << first_residual[i] 
           << endl;
    }
    for (int i=ifirst;i<=ilast;i++) {
      cout << "\tsolution_increment[" << i << "] = " 
           << solution_increment[i] << endl;
    }
#endif
    memcpy(&si[ifirst],&solution_increment[ifirst],
           (ilast-ifirst+1)*sizeof(double));
    for (int i=ifirst;i<=ilast;i++) 
//  int i=7;
    {
      for (int k=ifirst;k<=ilast;k++) {
        first_residual[k]=0.;
        solution_increment[k]=0.;
      }
      first_residual[i]=1.;
      postSmooth();
#ifdef DEBUG
      for (int k=ifirst;k<=ilast;k++) {
        cout << "\tfirst_residual[" << k << "] = " << first_residual[k] 
             << endl;
      }
      for (int k=ifirst;k<=ilast;k++) {
        cout << "\tsolution_increment[" << k << "] = " 
             << solution_increment[k] << endl;
      }
#endif
      if (abs(si[i]-solution_increment[j])>
      1.e-12*max(abs(si[i]),abs(solution_increment[j]))) {
        cout << "\tsi[" << i << "," << j << "] = " 
             << si[i] << " " << solution_increment[j] << endl;
      }
    }
  }
  delete [] si;
  ASSERT(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::checkCoarseGridProjectionRandom() const {
  TRACER_CALL(t,"Level::checkCoarseGridProjectionRandom");
  for (int I=coarser->ifirst;I<=coarser->ilast;I++) 
//int I=5;
  {
    for (int k=ifirst;k<=ilast;k++) {
      first_residual[k]=0.;
      solution_increment[k]=0.;
    }
    for (int K=coarser->ifirst;K<=coarser->ilast;K++) {
      coarser->solution_increment[K]=0.;
    }
    coarser->solution_increment[I]=1.;
#ifdef DEBUG
    for (int J=coarser->ifirst;J<=coarser->ilast;J++) {
      cout << "\tsolution_increment[" << J << "] = " 
           << coarser->solution_increment[J] << endl;
    }
#endif
    prolongation();
#ifdef DEBUG
    for (int j=ifirst;j<=ilast;j++) {
      cout << "\tsolution_increment[" << j << "] = " 
           << solution_increment[j] << endl;
    }
#endif
    updateResidual();
#ifdef DEBUG
    for (int j=ifirst;j<=ilast;j++) {
      cout << "\tAx_b[" << j << "] = " << Ax_b[j] << endl;
    }
#endif
    restriction();
#ifdef DEBUG
    for (int J=coarser->ifirst;J<=coarser->ilast;J++) {
      cout << "\tfirst_residual[" << J << "] = " 
           << coarser->first_residual[J] << endl;
    }
#endif
    coarser->multigridStep(coarser->Ax_b,coarser->solution_increment);
#ifdef DEBUG
    for (int J=coarser->ifirst;J<=coarser->ilast;J++) {
      cout << "\tsolution_increment[" << J << "] = " 
           << coarser->solution_increment[J] << endl;
    }
#endif
    prolongation();
#ifdef DEBUG
    for (int j=ifirst;j<=ilast;j++) {
      cout << "\tsolution_increment[" << j << "] = " 
           << solution_increment[j] << endl;
    }
#endif
    int j=F77NAME(idamax)(ilast-ifirst+1,&solution_increment[ifirst],1)
         -1+ifirst;
    cout << "\tnullspace[" << I << "] = " << solution_increment[j] << endl;
  }
  ASSERT(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Level::checkCoarseGridProjection() const {
  TRACER_CALL(t,"Level::checkCoarseGridProjection");
  for (int I=coarser->ifirst;I<=coarser->ilast;I++) 
//int I=5;
  {
    for (int k=ifirst;k<=ilast;k++) {
      first_residual[k]=0.;
      solution_increment[k]=0.;
    }
    for (int K=coarser->ifirst;K<=coarser->ilast;K++) {
      coarser->solution_increment[K]=0.;
    }
    coarser->solution_increment[I]=1.;
#ifdef DEBUG
    for (int J=coarser->ifirst;J<=coarser->ilast;J++) {
      cout << "\tsolution_increment[" << J << "] = " 
           << coarser->solution_increment[J] << endl;
    }
#endif
    prolongation();
#ifdef DEBUG
    for (int j=ifirst;j<=ilast;j++) {
      cout << "\tsolution_increment[" << j << "] = " 
           << solution_increment[j] << endl;
    }
#endif
    updateResidual();
#ifdef DEBUG
    for (int j=ifirst;j<=ilast;j++) {
      cout << "\tAx_b[" << j << "] = " << Ax_b[j] << endl;
    }
#endif
    restriction();
#ifdef DEBUG
    for (int J=coarser->ifirst;J<=coarser->ilast;J++) {
      cout << "\tfirst_residual[" << J << "] = " 
           << coarser->first_residual[J] << endl;
    }
#endif
    coarser->multigridStep(coarser->Ax_b,coarser->solution_increment);
#ifdef DEBUG
    for (int J=coarser->ifirst;J<=coarser->ilast;J++) {
      cout << "\tsolution_increment[" << J << "] = " 
           << coarser->solution_increment[J] << endl;
    }
#endif
    prolongation();
#ifdef DEBUG
    for (int j=ifirst;j<=ilast;j++) {
      cout << "\tsolution_increment[" << j << "] = " 
           << solution_increment[j] << endl;
    }
#endif
    int j=F77NAME(idamax)(ilast-ifirst+1,&solution_increment[ifirst],1)
         -1+ifirst;
    cout << "\tnullspace[" << I << "] = " << solution_increment[j] << endl;
  }
  ASSERT(0);
}
