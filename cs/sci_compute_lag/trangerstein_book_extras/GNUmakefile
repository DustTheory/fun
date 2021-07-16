#	ATLAS/ DASKR/ DASPK/ blitz-0.10/ blitzwave-0.7.1/ dfftpack/
#	glpk-4.54/ lapack-3.3.1/ linpack/ mathgl-2.1.3.1/ minpack/
#	petsc-3.4.3/ stl/ tnt/ toms630/ uncmin/ vode/
all : lapackpp minpack memdebug graphics gui
	cd ode; touch depends.gnu
	cd computer; touch depends.gnu
	cd plotting; touch depends.gnu
	cd linear_algebra; touch depends.gnu
	cd nonlinear_equations; touch depends.gnu
	cd least_squares; touch depends.gnu
	cd eigenvalues; touch depends.gnu
	cd iterative_linear_algebra; touch depends.gnu
	cd nonlinear_systems; touch depends.gnu
	cd linear_programming; touch depends.gnu
	cd approximation; touch depends.gnu
	cd integration; touch depends.gnu
	cd ivp; touch depends.gnu
	cd bvp2; touch depends.gnu
	cd ode; make all
	cd computer; make all
	cd plotting; make all
	cd linear_algebra; make all
	cd nonlinear_equations; make all
	cd least_squares; make all
	cd eigenvalues; make all
	cd iterative_linear_algebra; make all
	cd nonlinear_systems; make all
	cd linear_programming; make all
	cd approximation; make all
	cd integration; make all
	cd ivp; make all
	cd bvp2; make all
atlas : FORCE
	tar jxf atlas3.11.17.tar.bz2
	mkdir ATLAS/MyObj
	cd ATLAS/MyObj; ../configure
	cd ATLAS/MyObj; make
	cd ATLAS/MyObj; make check
	cd ATLAS/MyObj; make ptcheck
	cd ATLAS/MyObj; make time
	cd ATLAS/MyObj; make install
lapackpp : FORCE
	cd lapack++; make all
mathgl : FORCE
	tar -xvf mathgl-2.1.3.1.tar
	cd mathgl-2.1.3.1; cmake .
	cd mathgl-2.1.3.1; make
	cd mathgl-2.1.3.1; make install
minpack : FORCE
	tar -xvf minpack_tarfile
	cd minpack; make libminpack.a
memdebug : FORCE
	-cd memdebug; mkdir d
	-cd memdebug; mkdir o
	cd memdebug; touch depends.gnu
	cd memdebug; make all
gui : FORCE
	-cd gui; mkdir d
	-cd gui; mkdir o
	cd gui; touch depends.gnu
	cd gui; make all
graphics : FORCE
	-cd graphics; mkdir 1d
	-cd graphics; mkdir 1o
	-cd graphics; mkdir 2d
	-cd graphics; mkdir 2o
	-cd graphics; mkdir 3d
	-cd graphics; mkdir 3o
	cd graphics; touch depends.gnu
	cd graphics; make all
clean :
	cd lapack++; make clean
	cd minpack; rm *.o *.a
	cd memdebug; make clean
	cd gui; make clean
	cd ode; make clean
	cd computer; make clean
	cd graphics; make clean
	cd plotting; make clean
	cd linear_algebra; make clean
	cd nonlinear_equations; make clean
	cd least_squares; make clean
	cd eigenvalues; make clean
	cd iterative_linear_algebra; make clean
	cd nonlinear_systems; make clean
	cd linear_programming; make clean
	cd approximation; make clean
	cd integration; make clean
	cd ivp; make clean
	cd bvp2; make clean
FORCE :
