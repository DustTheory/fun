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
// "$Header: /home/faculty/johnt/cvs/deal_new/gui/GUICallback.C,v 1.1 2009/08/20 17:32:35 johnt Exp $"
#include "GUICallback.H"
#include "Tracer.H"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> 
void GUICallbackStructure<T>::printOn(ostream &os) const {
  os << "GUICallbackStructure:" << endl;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> GUICallbackStructureList<T>::
~GUICallbackStructureList() {
  while (notEmpty()) delete delAfter(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> 
void GUICallbackStructureList<T>::printOn(ostream &os) const {
  os << "GUICallbackStructureList:" << endl;
  ISLList< GUICallbackStructure<T> >::printOn(os);
}
