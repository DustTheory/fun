// "$Header:$"
//**********************************************************************
// Copyright 2006 John A. Trangenstein
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

#include "PermutationGenerator.H"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//for n>15, long total overflows and the following fails
long PermutationGenerator::getFactorial(int n) {
  long fact=1;
  for (int i=n;i>1;i--) fact*=i;
  return fact;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
PermutationGenerator::PermutationGenerator(int n) : a(n) {
  CHECK_TEST(n>=1);
  total=getFactorial(n);
  for (int i=0;i<n;i++) a[i]=i;
  num_left=total;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//adapted from http://www.merriampark.com/perm.htm
//See Kenneth H. Rosen, Discrete Mathematics and Its Applications, 
//  2nd edition (NY: McGraw-Hill, 1991), pp. 282-284.
const NumPtr<int>& PermutationGenerator::getNext() {
  if (num_left==total) {
    num_left--;
    return a;
  }
//find largest index j with a[j]<a[j+1]
  int n=a.getNumber();
  int j=n-2;
  while (a[j]>a[j+1]) j--;
//find index k such that a[k] is smallest integer greater than a[j]
//  to the right of a[j]
  int k=n-1;
  while (a[j]>a[k]) k--;
  int temp=a[k];
  a[k]=a[j];
  a[j]=temp;
//put tail end of permutation after jth position in increasing order
  int r=n-1;
  int s=j+1;
  while (r>s) {
    temp=a[s];
    a[s]=a[r];
    a[r]=temp;
    r--;
    s++;
  }
  num_left--;
  return a;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//adapted from http://www.fortran.com/random1.f90
//  Courtesy of HDK@psuvm.psu.edu Thu Aug 25 09:11:00 MDT 1994
void  RandomPermutation(NumPtr<int> &p) {
//Generate a random permutation, P, of the first N integers.
//  equivalent to sampling WITHOUT REPLACEMENT
//  Adaptation of Knuth Volume 2, Algorithm 3.4.2P.
  int n=p.getNumber();
  for (int i=0;i<n;i++) p[i]=i;
  for (int i=0;i<n;i++) {
    double u=drand48();
    int k=static_cast<int>(u*static_cast<double>(n-i))+i;
    int t=p[i];
    p[i]=p[k];
    p[k]=t;
  }
}
