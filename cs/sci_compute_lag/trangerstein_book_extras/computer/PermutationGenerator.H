#ifndef PERMUTATION_GENERATOR_H
#define PERMUTATION_GENERATOR_H "$Header:$"
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
#include "NumPtr.H"

//adapted from http://www.merriampark.com/perm.htm
//See Kenneth H. Rosen, Discrete Mathematics and Its Applications, 
//   2nd edition (NY: McGraw-Hill, 1991), pp. 282-284.
class PermutationGenerator {
  private:
    NumPtr<int> a;
    long num_left;
    long total;

    static long getFactorial(int);
  public:
//  for n>15, long total overflows and the following fails
    PermutationGenerator(int);
    inline bool hasMore() { return num_left>0; }
    const NumPtr<int>& getNext();
};

//adapted from http://www.fortran.com/random1.f90
//  courtesy of HDK@psuvm.psu.edu Thu Aug 25 09:11:00 MDT 1994
void  RandomPermutation(NumPtr<int>&);
#endif
