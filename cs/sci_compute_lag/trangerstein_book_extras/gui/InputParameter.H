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
#ifndef INPUT_PARAMETER_H
#define INPUT_PARAMETER_H "$Header: /home/faculty/johnt/cvs/deal_new/gui/InputParameter.H,v 1.1 2009/08/20 17:32:35 johnt Exp $"
#include <iomanip>
#include <vector>
#include "Const.H"
#include "Errors.H"
#include "MemoryDebugger.H"
#include "VirtualInput.H"

template<class T> class InputParameter : public VirtualInput {
  protected:
    T *ptr_data;
    T lower_bound,upper_bound;
    char *vstring;
  public:
    InputParameter(T &p,const char *n,T l,T u,const char* =0) : 
    VirtualInput(n), ptr_data(&p),lower_bound(l),upper_bound(u) {
      vstring=OPERATOR_NEW_BRACKET(char,LENGTH_NAME);
      CHECK_TEST(l<=u);
      if (p<l) p=l;
      if (p>u) p=u;
    }
    virtual ~InputParameter() { 
      ptr_data=0; delete vstring; vstring=0;
    }

    virtual const char* getValue() const;
    virtual ifstream& formattedRead(ifstream &ifs) {
      ifs >> *ptr_data;
      *ptr_data=min(upper_bound,max(lower_bound,*ptr_data));
      return ifs;
    }
    virtual ifstream& unformattedRead(ifstream &ifs) {
      ifs.read(reinterpret_cast<char*>(&lower_bound),
        sizeof(lower_bound));
      ifs.read(reinterpret_cast<char*>(&upper_bound),
        sizeof(upper_bound));
      ifs.read(reinterpret_cast<char*>(ptr_data),sizeof(T));
      return ifs;
    }
    virtual ofstream& unformattedWrite(ofstream &ofs) {
      ofs.write(reinterpret_cast<char*>(&lower_bound),
        sizeof(lower_bound));
      ofs.write(reinterpret_cast<char*>(&upper_bound),
        sizeof(upper_bound));
      ofs.write(reinterpret_cast<char*>(ptr_data),sizeof(T));
      return ofs;
    }
    virtual void printOn(ostream &os=cerr) const {
      os << "InputParameter: ptr_data = " << ptr_data
         << "\n\tlower_bound = " << lower_bound
         << "\n\tupper_bound = " << upper_bound << endl;
      VirtualInput::printOn(os);
    }
};
#endif
