#include <stdlib.h>
#include "LinearProgram.H"

int main(int /*argc*/, char* /*argv[]*/) {
//set up space for problem arrays
//note that vectors are arrays with 1 entry in each column (or row)
//note that c is a row array, b is a column array
  Matrix<float> A(3,5);
  Matrix<float> b(3,1);
  Matrix<float> c(1,5);

//set array values
//Giapetto problem, p. 45 of Winston
  A(0,0)=1.; A(0,1)=0.; A(0,2)=0.; A(0,3)= 2.; A(0,4)= 1.; b(0,0)=100.;
  A(1,0)=0.; A(1,1)=1.; A(1,2)=0.; A(1,3)= 1.; A(1,4)= 1.; b(1,0)= 80.;
  A(2,0)=0.; A(2,1)=0.; A(2,2)=1.; A(2,3)= 1.; A(2,4)= 0.; b(2,0)= 40.;
  c(0,0)=0.; c(0,1)=0.; c(0,2)=0.; c(0,3)=-3.; c(0,4)=-2.;

//create a standard-form linear program
//this sets up the work space to solve min c.x st Ax=b, x>0
  SFLinearProgram<float> LP(A,b,c);

//use artificial variables if necessary to get started
  LP.findBasicFeasibleGuess();

//if feasible, perform simplex steps until done
  while (LP.simplexStep()==primal_feasible) {
    cout << "\nstatus = " << LP.currentStatus() 
	 << endl;
    cout << "current value = " << LP.currentValue() << endl;
    cout << "current solution = " << LP.currentPrimalSolution() << endl;
  }

  cout << "\n\nafter iteration, status = " << LP.currentStatus() <<endl;
  cout << "final value = " << LP.currentValue() << endl;
  cout << "final solution = " << LP.currentPrimalSolution() << endl;

  return EXIT_SUCCESS;
}
