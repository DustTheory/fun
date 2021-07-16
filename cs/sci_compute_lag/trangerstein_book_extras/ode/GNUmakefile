#a.out : euler.f
#	gfortran euler.f
all : f_exec c_exec C_exec refine efficiency
f_exec : euler.f
	gfortran -g -o f_exec euler.f
c_exec : euler.c
	gcc -g -o c_exec euler.c
C_exec : euler.C
	g++ -g -o C_exec euler.C
refine : mesh_refinement.f
	gfortran -g -o refine mesh_refinement.f
efficiency : meshRefinement.C
	g++ -g -o efficiency meshRefinement.C
clean :
	-rm *.o core*
	-rm f_exec c_exec C_exec refine efficiency
