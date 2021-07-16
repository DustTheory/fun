#include <complex>
#include <float.h>
#include <stdlib.h>
//#include "MemoryDebugger.H"
//#include "SetTraps.H"
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

//template<> const complex<float>
//  Vector<float,complex<float> >::undefined_(HUGE_VAL,HUGE_VAL);

int main( int /* argc */, char** /* argv */) {
#ifdef MEM_DEBUG
//MemoryDebugger md(1);
#endif
/*
  cout << "sizeof(complex<float>) = "
       << sizeof(complex<float>)/sizeof(float) << endl;
  complex<float> *zarray=OPERATOR_NEW_BRACKET(complex<float>,2);
  zarray[0].real()=1.;
  zarray[0].imag()=2.;
  zarray[1].real()=3.;
  zarray[1].imag()=4.;
  float *rarray=reinterpret_cast<float*>(zarray);
  cout << "rarray = " << rarray[0] << " " << rarray[1] << " "
    << rarray[2] << " " << rarray[3] << endl;
*/
//testVector(float(3.),float(2.));
//testVector(float(3.),complex<float>(2.,1.));

//testMatrix(float(3.),float(2.));
//testMatrix(float(3.),complex<float>(2.,1.));
//testSquareMatrix(float(3.),float(2.));
//testSquareMatrix(float(3.),complex<float>(2.,1.));

//testLowerTrapezoidalMatrix(float(3.),float(2.));
//testLowerTrapezoidalMatrix(float(3.),complex<float>(2.,1.));
//testLowerTriangularMatrix(float(3.),float(2.));
//testLowerTriangularMatrix(float(3.),complex<float>(2.,1.));
//testUnitLowerTrapezoidalMatrix(float(3.),float(2.));
//testUnitLowerTrapezoidalMatrix(float(3.),complex<float>(2.,1.));
//testUnitLowerTriangularMatrix(float(3.),float(2.));
//testUnitLowerTriangularMatrix(float(3.),complex<float>(2.,1.));
//testUpperTrapezoidalMatrix(float(3.),float(2.));
//testUpperTrapezoidalMatrix(float(3.),complex<float>(2.,1.));
//testUpperTriangularMatrix(float(3.),float(2.));
//testUpperTriangularMatrix(float(3.),complex<float>(2.,1.));
//testUnitUpperTrapezoidalMatrix(float(3.),float(2.));
//testUnitUpperTrapezoidalMatrix(float(3.),complex<float>(2.,1.));
//testUnitUpperTriangularMatrix(float(3.),float(2.));
//testUnitUpperTriangularMatrix(float(3.),complex<float>(2.,1.));
  
//testOrthogonalMatrix(float(3.),float(2.));
//testOrthogonalMatrix(float(3.),complex<float>(2.,1.));
  
//testSymmetricMatrix(float(3.),float(2.));
//testSymmetricMatrix(float(3.),complex<float>(2.,1.));
//testSymmetricPositiveMatrix(float(3.),float(2.));
//testSymmetricPositiveMatrix(float(3.),complex<float>(2.,1.));

//testTridiagonalMatrix(float(3.),float(2.));
//testTridiagonalMatrix(float(3.),complex<float>(2.,1.));
//testSymmetricTridiagonalMatrix(float(3.),float(2.));
//testSymmetricTridiagonalMatrix(float(3.),complex<float>(2.,1.));
//testSymmetricPositiveTridiagonalMatrix(float(3.),float(2.));
//testSymmetricPositiveTridiagonalMatrix(float(3.),complex<float>(2.,1.));
//testDiagonalMatrix(float(3.),float(2.));
//testDiagonalMatrix(float(3.),complex<float>(2.,1.));

//testUpperHessenbergMatrix(float(3.),float(2.));
//testUpperHessenbergMatrix(float(3.),complex<float>(2.,1.));
//testBandMatrix(float(3.),float(2.));
//testBandMatrix(float(3.),complex<float>(2.,1.));
//testSymmetricBandMatrix(float(3.),float(2.));
//testSymmetricBandMatrix(float(3.),complex<float>(2.,1.));
//testSymmetricPositiveBandMatrix(float(3.),float(2.));
//testSymmetricPositiveBandMatrix(float(3.),complex<float>(2.,1.));
//testSpecializedMatrix(float(3.),float(2.));
//testSpecializedMatrix(float(3.),complex<float>(2.,1.));
 
//testGaussianFactorization(float(3.),float(2.));
//testGaussianFactorization(float(3.),complex<float>(2.,1.));
//testCholeskyFactorization(float(3.),float(2.));
//testCholeskyFactorization(float(3.),complex<float>(2.,1.));
//testMDMtFactorization(float(3.),float(2.));
//testMDMtFactorization(float(3.),complex<float>(2.,1.));

//testGramSchmidtQRFactorization(float(3.),float(2.));
//testGramSchmidtQRFactorization(float(3.),complex<float>(2.,1.));
//testHouseholderQRFactorization(float(3.),float(2.));
//testHouseholderQRFactorization(float(3.),complex<float>(2.,1.));
//testCompleteOrthogonalDecomposition(float(3.),float(2.));
  testCompleteOrthogonalDecomposition(float(3.),complex<float>(2.,1.));
//testSingularValueDecomposition(float(3.),float(2.));
//testSingularValueDecomposition(float(3.),complex<float>(2.,1.));

//testSFLinearProgram(float(3.));
//testLinearProgram(float(3.));
  return EXIT_SUCCESS;
}
extern "C" {
void MAIN__(void) {;}
}
