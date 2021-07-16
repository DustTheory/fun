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
// "$Header:$"
#include "MemoryDebugger.H"
#include <cstring>
#include <cstdio>

//buffer_type MemoryNode::beforeval=0x7ff0123456789abc;
//buffer_type MemoryNode::afterval =0x7ffedcba98765432;
buffer_type MemoryNode::beforeval=0x12345678;
buffer_type MemoryNode::afterval=0x87654321;
char MemoryNode::mallocval=0xEE;
size_t MemoryNode::buffer_bytes=0U;
size_t MemoryNode::buffer_size=0U;

bool MemoryDebugger::inited = false;
int MemoryDebugger::count=0;
size_t MemoryDebugger::buffer_size=1U;
long MemoryDebugger::numalloc=0;
long MemoryDebugger::maxalloc=0;
MemoryNode* MemoryDebugger::head=0;

// The following should be selected to give maximum probability that
// pointers loaded with these values will cause an obvious crash. On
// Unix machines, a large value will cause a segment fault.
// MALLOCVAL is the value to set malloc'd data to.
char MemoryDebugger::badval=0x7A;

void MemoryNode::create(size_t n,void *d) {
  CHECK_POSITIVE(buffer_size)
//first we have storage for this
  char *dd=static_cast<char*>(d);
  size_t nd=sizeof(MemoryNode); Align(nd);

//next, we have storage for the before buffer
  buffer_type *val=static_cast<buffer_type*>(static_cast<void*>(dd+nd));
  for (size_t i=0;i<buffer_size;i++,val++) {
    *val=beforeval;
  }
  nd+=buffer_bytes;

//next, we have storage for the user data
  data=static_cast<void*>( dd+nd );
  nbytes=n; Align(nbytes);
  memset(data,mallocval,nbytes);
  nd+=nbytes;

//finally, we have storage for the after buffer
  val=static_cast<buffer_type*>(static_cast<void*>(dd+nd));
  for (size_t i=0;i<buffer_size;i++,val++) {
    *val=afterval;
  }
}

MemoryNode::MemoryNode(size_t n,void *d) : line_number(-1),bk(0),fd(0) {
  create(n,d);
  file[0]='\0';
}

MemoryNode::MemoryNode(size_t n,void *d,const char *f,int l) : 
line_number(l),bk(0),fd(0) {
  create(n,d);
  int len=strlen(f);
  len=min(FILE_LENGTH-1,len);
  strncpy(file,f,len);
  file[len]='\0';
}

bool MemoryNode::checkBuffer() const {
  CHECK_POSITIVE(buffer_size)
  char *dd=static_cast<char*>(data);
  buffer_type *val=
    static_cast<buffer_type*>(static_cast<void*>(dd-buffer_bytes));
  for (size_t i=0;i<buffer_size;i++,val++) {
    if (*val!=beforeval) {
      cout << "Pointer " << data << " underrun" << endl;
      printOn(cout);
      return false;
    }
  } 
  val=static_cast<buffer_type*>(static_cast<void*>(dd+nbytes));
  for (size_t i=0;i<buffer_size;i++,val++) {
    if (*val!=afterval) {
      cout << "Pointer " << data << " overrun" << endl;
      printOn(cout);
      return false;
    }
  } 
  return true;
}

void MemoryNode::printOrigin(ostream &os) const {
  if (line_number>=0) {
    os << "\tallocated in file " << file << " at line "
       << line_number;
  }
}

void MemoryNode::printOn(ostream &os) const {
  os << "MemoryNode:: nbytes=" << nbytes << ",data=" << data;
  printOrigin(os);
  os << endl;
}

MemoryDebugger::MemoryDebugger(size_t bs) { 
  if (!inited) {
    buffer_size=bs;
    inited=true;
#ifdef DEBUG
//  F77_NAME(test_types)();
#endif
  }
}

MemoryDebugger::~MemoryDebugger() {
  if (inited) {
    if (count!=0) {
      cout << "in MemoryDebugger::~MemoryDebugger count = " << count 
           << endl;
      register MemoryNode *mn=head;
      bool trouble=false;
      for ( ; mn!=head; mn = mn->next() ) {
        trouble |= mn->allocatedByAMR();
        cout << "Unfreed pointer: ";
        mn->printOn(cout);
      }
      CHECK_TEST(!trouble);
    }
  }
  inited = false;
  buffer_size=0;
}

void* MemoryDebugger::malloc(size_t n) {
  if (inited) {
    size_t sz=MemoryNode::allocationSize(n,buffer_size);
    void *data=::malloc(sz);
    if (data) {
      MemoryNode *mn=new(data) MemoryNode(n,data);
      append(mn);
#ifdef DEBUG
//    printf("\tMemoryDebugger::malloc data = %p\n",data);
//    mn->printOn(cout);
#endif
//    bytes allocated may be more than requested because of alignment:
      numalloc += mn->size();
      if (numalloc > maxalloc) maxalloc = numalloc;
      return mn->ptr();
    } else {
      cout << 
          "\tMemoryDebugger::malloc failed to allocate pointer of size "
           << n << endl;
      cout << "\tsz = " << sz << endl;
      ASSERT(data!=0);
    } 
    return data;
  } else {
    return ::malloc(n);
  }
}

void* MemoryDebugger::malloc(size_t n,const char *f,int l) {
//came from OPERATOR_NEW or OPERATOR_NEW_BRACKET, so watch for leaks
  if (inited) {
    size_t sz=MemoryNode::allocationSize(n,buffer_size);
    void *data=::malloc(sz);
    if (data) {
      MemoryNode *mn=new(data) MemoryNode(n,data,f,l);
      append(mn);
#ifdef DEBUG
//    printf("\tMemoryDebugger::malloc called from file %s at line %d, data = %p\n",f,l,data);
//    mn->printOn(cout);
#endif
//    bytes allocated may be more than requested because of alignment:
      numalloc += mn->size();
      if (numalloc > maxalloc) maxalloc = numalloc;
      return mn->ptr();
    } else {
      cout << "\tMemoryDebugger::malloc failed to alloc ptr of size "
           << n << endl;
      cout << "\tsz = " << sz << endl;
      ASSERT(data!=0);
    } 
    return data;
  } else return ::malloc(n);
}

bool MemoryDebugger::owns(MemoryNode *mn) {
#ifdef DEBUG
//printf("\tMemoryDebugger::owns, mn,head = %p, %p\n",mn,head);
//mn->printOn(cout);
//head->printOn(cout);
#endif
  register MemoryNode *n=head;
  if (mn==n) return true;
  for (n=n->next();n!=head;n=n->next()) {
    if (n==mn) return true;
  }
  return false;
}

void MemoryDebugger::free(void *ptr) {
#ifdef DEBUG
//printf("\tMemoryDebugger::free(void*) ptr = %p\n",ptr);
#endif
  if (!ptr) return;
//since new system code is leaky, don't watch for leaks
  if (inited) {
    MemoryNode *mn = convert(ptr);
    ASSERT(mn->checkBuffer());

    size_t bytes=mn->size();
    numalloc -= bytes;
    if (numalloc < 0) {
      cout << "freeing more bytes than allocated" << endl;
      cout << "numalloc = " << numalloc << " nbytes = " << bytes <<endl;
      CHECK_NONNEGATIVE(numalloc)
    }
    remove(mn);

//  Stomp on the freed storage to help detect references
//  after the storage was freed.
//  mn->printOn will give segmentation fault after this
    memset(mn,badval,mn->allocationSize(bytes,buffer_size));
    ::free(static_cast<void *>(mn));
//  endConvert();
  } else {
    ::free(ptr);
  }
}

void MemoryDebugger::check() {
  if (inited) {
    register MemoryNode *mn=head;
    for (;mn!=head;mn=mn->next()) ASSERT(mn->checkBuffer());
  }
}

void MemoryDebugger::checkPtr(void *p) {
  if (p!=0 && inited) {
    MemoryNode *mn=convert(p);
    mn->printOn(cout);
    ASSERT(mn->checkBuffer());
  }
}

void MemoryDebugger::printOn(ostream &os) {
  os << "MemoryDebugger: inited = " << inited << endl;
  os << "\tcount = " << count << endl;
  os << "\tbuffer_size = " << buffer_size << endl;
  os << "\tnumalloc = " << numalloc << endl;
  os << "\tmaxalloc = " << maxalloc << endl;
  os << "\tbadval = " << badval << endl;
  register MemoryNode *mn=head;
  for (;mn;mn=mn->next()) mn->printOn(os);
}

void* operator new(size_t n,const std::nothrow_t&) throw() {
  void *p=0;
  if (MemoryDebugger::active()) {
    p=MemoryDebugger::malloc(n);
  } else {
    p=malloc(n);
    ASSERT(p!=0);
  }
#ifdef DEBUG
//printf("\tin operator new(size_t,nothrow) n = %d p = %p\n",n,p);
#endif
  return p;
}

void* operator new(size_t n) throw(std::bad_alloc) {
  void *p=operator new(n,nothrow);
#ifdef DEBUG
//printf("\tin operator new(size_t) n = %d p = %p\n",n,p);
#endif
  return p;
}

void* operator new(size_t n,const char *f,int l) {
  void *p=0;
  if (MemoryDebugger::active()) {
    p=MemoryDebugger::malloc(n,f,l);
#ifdef DEBUG
//  printf("\t(active) operator new(size_t,char*,int) n = %d p = %p\n",n,p);
#endif
  } else {  
    p=malloc(n);
#ifdef DEBUG
//  printf("\t(inactive) operator new(size_t,char*,int) n = %d p = %p\n",n,p);
#endif
  }
  return p;
}

//this is called before the array entries are constructed
void* operator new[](size_t n,const std::nothrow_t&) throw() {
  void *p=0;
  if (MemoryDebugger::active()) {
    p=MemoryDebugger::malloc(n);
#ifdef DEBUG
//  printf("\t(active) operator new[](size_t,nothrow) n = %d p = %p\n",n,p);
#endif
  } else {
    p=malloc(n);
    ASSERT(p!=0);
#ifdef DEBUG
//  printf("\t(inactive) operator new[](size_t,nothrow) n = %d p = %p\n",n,p);
#endif
  }
  return p;
}

//this is called before the array entries are constructed
void* operator new[](size_t n) throw(std::bad_alloc) {
  void *p=operator new[](n,nothrow);
#ifdef DEBUG
//printf("\tin operator new[](size_t) n = %d p = %p\n",n,p);
#endif
  return p;
}

//this is called before the array entries are constructed
void* operator new[](size_t n,const char *f,int l) {
  void *p=0;
  if (MemoryDebugger::active()) {
    p=MemoryDebugger::malloc(n,f,l);
#ifdef DEBUG
//  printf("\t(active) operator new[](size_t,char*,int) n = %d p = %p\n",n,p);
#endif
  } else  {
    p=malloc(n);
#ifdef DEBUG
//  printf("\t(inactive) operator new[](size_t,char*,int) n = %d p = %p\n",n,p);
#endif
  }
  return p;
}

//this is called after the class destructor
void  operator delete(void* p,const std::nothrow_t&) throw() {
  if (MemoryDebugger::active()) {
#ifdef DEBUG
//  printf("\t(active) operator delete(void*,nothrow) p = %p\n",p);
#endif
    MemoryDebugger::free(p);
  } else {
#ifdef DEBUG
//  printf("\t(inactive) operator delete(void*,nothrow) p = %p\n",p);
#endif
    free(p);
  }
}

//this is called after the class destructor
void  operator delete(void* p) throw() {
#ifdef DEBUG
//printf("\tin operator delete(void*) p = %p\n",p);
#endif
  operator delete(p,nothrow);
}

//this is called after the array entries are destructed
void operator delete[](void *p,const std::nothrow_t&) throw() {
#ifdef DEBUG
//printf("\tin operator delete[](void*,nothrow) p = %p\n",p);
#endif
  if (MemoryDebugger::active()) {
#ifdef DEBUG
//  printf("\t(active) operator delete[](void*,nothrow) p = %p\n",p);
#endif
    MemoryDebugger::free(p);
  } else {
#ifdef DEBUG
//  printf("\t(inactive) operator delete[](void*,nothrow) p = %p\n",p);
#endif
    free(p);
  }
}

//this is called after the array entries are destructed
void operator delete[](void *p) throw() {
#ifdef DEBUG
//printf("\tin operator delete[](void*) p = %p\n",p);
#endif
  operator delete[](p,nothrow);
}

//for calls from C or Fortran:
extern "C" {
void mem_check_ptr(void *p) { MemoryDebugger::checkPtr(p); }
void mem_check() { MemoryDebugger::check(); }
void F77_NAME(mem_check)() { MemoryDebugger::check(); }
}
