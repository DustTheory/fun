#include <complex>
#include <float.h>
#include <stdlib.h>
#include "MemoryDebugger.H"
#include "SetTraps.H"
#include "LaEnum.C"
#include "Vector.C"
#include "BandMatrix.C"
#include "Matrix.C"
#include "OrthogonalMatrix.C"
#include "SquareMatrix.C"
#include "SymmetricMatrix.C"
#include "TrapezoidalMatrix.C"
#include "SpecializedMatrix.C"
#include "GaussianFactorization.C"
#include "CholeskyFactorization.C"
#include "MDMtFactorization.C"
#include "GramSchmidtQRFactorization.C"
#include "HouseholderQRFactorization.C"
#include "CompleteOrthogonalDecomposition.C"
#include "SingularValueDecomposition.C"
#include "LinearProgram.C"

//template<> const complex<double>
//  Vector<double,complex<double> >::undefined_(HUGE_VAL,HUGE_VAL);

int main( int /* argc */, char** /* argv */) {
#ifdef MEM_DEBUG
  MemoryDebugger md(1);
#endif
/*
  cout << "sizeof(complex<double>) = "
       << sizeof(complex<double>)/sizeof(float) << endl;
  complex<double> *zarray=OPERATOR_NEW_BRACKET(complex<double>,2);
  zarray[0].real()=1.;
  zarray[0].imag()=2.;
  zarray[1].real()=3.;
  zarray[1].imag()=4.;
  double *rarray=reinterpret_cast<double*>(zarray);
  cout << "rarray = " << rarray[0] << " " << rarray[1] << " "
    << rarray[2] << " " << rarray[3] << endl;
*/
//testVector(double(3.),double(2.));
//testVector(double(3.),complex<double>(2.,1.));

//testMatrix(double(3.),double(2.));
//testMatrix(double(3.),complex<double>(2.,1.));
//testSquareMatrix(double(3.),double(2.));
//testSquareMatrix(double(3.),complex<double>(2.,1.));

//testLowerTrapezoidalMatrix(double(3.),double(2.));
//testLowerTrapezoidalMatrix(double(3.),complex<double>(2.,1.));
//testLowerTriangularMatrix(double(3.),double(2.));
//testLowerTriangularMatrix(double(3.),complex<double>(2.,1.));
//testUnitLowerTrapezoidalMatrix(double(3.),double(2.));
//testUnitLowerTrapezoidalMatrix(double(3.),complex<double>(2.,1.));
//testUnitLowerTriangularMatrix(double(3.),double(2.));
//testUnitLowerTriangularMatrix(double(3.),complex<double>(2.,1.));
//testUpperTrapezoidalMatrix(double(3.),double(2.));
//testUpperTrapezoidalMatrix(double(3.),complex<double>(2.,1.));
//testUpperTriangularMatrix(double(3.),double(2.));
//testUpperTriangularMatrix(double(3.),complex<double>(2.,1.));
//testUnitUpperTrapezoidalMatrix(double(3.),double(2.));
//testUnitUpperTrapezoidalMatrix(double(3.),complex<double>(2.,1.));
//testUnitUpperTriangularMatrix(double(3.),double(2.));
//testUnitUpperTriangularMatrix(double(3.),complex<double>(2.,1.));

//testOrthogonalMatrix(double(3.),double(2.));
//testOrthogonalMatrix(double(3.),complex<double>(2.,1.));
  
//testSymmetricMatrix(double(3.),double(2.));
//testSymmetricMatrix(double(3.),complex<double>(2.,1.));
//testSymmetricPositiveMatrix(double(3.),double(2.));
//testSymmetricPositiveMatrix(double(3.),complex<double>(2.,1.));

//testTridiagonalMatrix(double(3.),double(2.));
//testTridiagonalMatrix(double(3.),complex<double>(2.,1.));
//testSymmetricTridiagonalMatrix(double(3.),double(2.));
//testSymmetricTridiagonalMatrix(double(3.),complex<double>(2.,1.));
//testSymmetricPositiveTridiagonalMatrix(double(3.),double(2.));
//testSymmetricPositiveTridiagonalMatrix(double(3.),complex<double>(2.,1.));
//testDiagonalMatrix(double(3.),double(2.));
//testDiagonalMatrix(double(3.),complex<double>(2.,1.));
 
//testUpperHessenbergMatrix(double(3.),double(2.));
//testUpperHessenbergMatrix(double(3.),complex<double>(2.,1.));
//testBandMatrix(double(3.),double(2.));
//testBandMatrix(double(3.),complex<double>(2.,1.));
//testSymmetricBandMatrix(double(3.),double(2.));
//testSymmetricBandMatrix(double(3.),complex<double>(2.,1.));
//testSymmetricPositiveBandMatrix(double(3.),double(2.));
//testSymmetricPositiveBandMatrix(double(3.),complex<double>(2.,1.));
//testSpecializedMatrix(double(3.),double(2.));
//testSpecializedMatrix(double(3.),complex<double>(2.,1.));

//testGaussianFactorization(double(3.),double(2.));
//testGaussianFactorization(double(3.),complex<double>(2.,1.));
//testCholeskyFactorization(double(3.),double(2.));
//testCholeskyFactorization(double(3.),complex<double>(2.,1.));
//testMDMtFactorization(double(3.),double(2.));
//testMDMtFactorization(double(3.),complex<double>(2.,1.));

//testGramSchmidtQRFactorization(double(3.),double(2.));
//testGramSchmidtQRFactorization(double(3.),complex<double>(2.,1.));
//testHouseholderQRFactorization(double(3.),double(2.));
//testHouseholderQRFactorization(double(3.),complex<double>(2.,1.));
//testCompleteOrthogonalDecomposition(double(3.),double(2.));
//testCompleteOrthogonalDecomposition(double(3.),complex<double>(2.,1.));
  testSingularValueDecomposition(double(3.),double(2.));
//testSingularValueDecomposition(double(3.),complex<double>(2.,1.));

//testSFLinearProgram(double(3.));
//testLinearProgram(double(3.));
  return EXIT_SUCCESS;
}
extern "C" {
void MAIN__(void) {;}
}
