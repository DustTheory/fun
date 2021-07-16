// "$Header:$"
//----------------------------  tensor.cc  ---------------------------
//    $Id: tensor.cc,v 1.15 2002/09/12 22:34:30 wolf Exp $
//    Version: $Name: Version-4-0-0 $
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  tensor.cc  ---------------------------
//
//modified from deal.II/base/source/tensor.cc
//  by John Trangenstein, August 2009
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

#include <Tensor.H>
#include <cmath>

template class Tensor<1,1>;
template class Tensor<1,2>;
template class Tensor<1,3>;
template class Tensor<2,1>;
template class Tensor<2,2>;
template class Tensor<2,3>;
template class Tensor<3,1>;
template class Tensor<3,2>;
template class Tensor<3,3>;
template class Tensor<4,1>;
template class Tensor<4,2>;
template class Tensor<4,3>;
