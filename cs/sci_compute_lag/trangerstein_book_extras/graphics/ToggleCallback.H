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
#ifndef TOGGLE_CALLBACK_H
#define TOGGLE_CALLBACK_H "$Header: /home/faculty/johnt/cvs/deal_new/graphics/ToggleCallback.H,v 1.1 2009/08/20 17:31:46 johnt Exp $"
//#include "Tracer.H"

template<class T> class ToggleButtonArrayCallbackStructure {
  private:
    T *data;
    int toggle_number;
  public:
    ToggleButtonArrayCallbackStructure(T *d,int tn) : 
    data(d),toggle_number(tn) { }
    int getToggleNumber() const { return toggle_number; }
    T* getData() const { return data; }
    void printOn(ostream &os) const {
      os << "ToggleButtonArrayCallbackStructure:" 
         << "\ttoggle_number = " << toggle_number
         << "\tdata = " << data << endl;
    }
};
#endif
