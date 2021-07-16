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
#ifndef MEMORY_DEBUGGER_H
#define MEMORY_DEBUGGER_H "$Header$"

#include "limits.h"
#include "malloc.h"
#include <string>
#include "Arch.H"
#include "Errors.H"

//   o Check values are inserted before and after the alloc'ed data
//     to detect pointer underruns and overruns.
//   o Free'd pointers are checked against alloc'ed pointers.
//   o Free'd storage is cleared to smoke out references to free'd data.
//   o The routine mem_checkptr() is provided to check an alloc'ed
//     pointer.

typedef unsigned long buffer_type;
typedef unsigned long print_type;
#define FILE_LENGTH 80

//the following two defines are machine-dependent
#define ALIGN_SIZE 8
#define ALIGN_EXPONENT 3    //2^ALIGN_EXPONENT=ALIGN_SIZE
inline void Align(size_t &n) {
  CHECK_SAME(ALIGN_SIZE,CHAR_BIT); //check this machine against define's
  size_t q=n>>ALIGN_EXPONENT;
  size_t m=q*ALIGN_SIZE;
  n=  m<n ? m+ALIGN_SIZE : m;
}

class MemoryNode {
  private:
    size_t nbytes;  // size of allocation
    void* data;     // allocated data
    char file[FILE_LENGTH]; // must be of fixed length avoid calling new
    int line_number;
    MemoryNode *bk;
    MemoryNode *fd;

    static buffer_type beforeval;
    static buffer_type afterval;
    static char mallocval;
    static size_t buffer_bytes;
    static size_t buffer_size;

    void create(size_t,void*);
  public:
    MemoryNode() : nbytes(0),data(0),bk(0),fd(0) {}
    MemoryNode(size_t,void*);
    MemoryNode(size_t,void*,const char*,int);
    static size_t allocationSize(size_t n,size_t buf_sz) {
      if (buffer_size==0U) {
        buffer_bytes=sizeof(buffer_type)*buf_sz; Align(buffer_bytes);
        buffer_size=buffer_bytes/sizeof(buffer_type);
      }
      size_t nb=bytesBeforeData();
      nb+=n; Align(nb);
      nb+=sizeof(buffer_type)*buffer_size; Align(nb); 
      return nb;
    }
    virtual ~MemoryNode() {;}
    static size_t bytesBeforeData() {
      ASSERT(buffer_size>0);
      size_t nb=sizeof(MemoryNode); Align(nb);
      nb+=sizeof(buffer_type)*buffer_size; Align(nb);
      return nb;
    }

    void* ptr() const { return data; }
    size_t size() const { return nbytes; }
    bool allocatedByAMR() const { return line_number!=-1; }
    MemoryNode* prev() const { return bk; }
    MemoryNode* next() const { return fd; }
    void selfConnect() { bk=fd=this; }
    void placeAfter(MemoryNode *mn) { 
      fd=mn->fd; bk=mn;
      mn->fd->bk=this; mn->fd=this;
    }
    void disconnect() {
      bk->fd=fd;
      fd->bk=bk;
    }
    bool inBounds(const void *p) const { 
      return p>= data && p< static_cast<char*>(data)+nbytes;
    }
    bool checkBuffer() const;
    void printOn(ostream&) const;
    void printOrigin(ostream&) const;
};

class MemoryDebugger {
  private:
    static bool inited;     // TRUE if mem package is initialized
    static int count;          // # of unfreed allocs
    static size_t buffer_size; // number of beforeval's and afterval's
    static long numalloc;      // current # bytes allocated
    static long maxalloc;      // max # bytes allocated
    static MemoryNode *head;
    static char badval;

    static void append(MemoryNode *mn) {
      if (head==0) {
        head=mn;
        mn->selfConnect();
      } else mn->placeAfter(head->prev());
      count++;
    }
    static bool owns(MemoryNode*);
    static void remove(MemoryNode *mn) {
      CHECK_TEST(owns(mn));
      count--;
      if (count < 0) {
        cout << "More frees than allocs";
        ASSERT(count>=0);
      }
      if (mn==head) head=(head->next()==head ? 0 : head->next());
      mn->disconnect();
    }
  public:
  MemoryDebugger(size_t=1U);
  ~MemoryDebugger();
    static bool active() { return inited; }
    static void check() ;
    static void* malloc(size_t);
    static void* malloc(size_t,const char *,int);
    static void free(void *);
    static void checkPtr(void*);
    static MemoryNode* convert(void *p) {
      ASSERT(inited);
      char *pp=reinterpret_cast<char*>(p)-MemoryNode::bytesBeforeData();
      return reinterpret_cast<MemoryNode*>(pp);
    }
    static const MemoryNode* convert(const void *p) {
      ASSERT(inited);
      const char *pp=reinterpret_cast<const char*>(p)
                    -MemoryNode::bytesBeforeData();
      return reinterpret_cast<const MemoryNode*>(pp);
    }
    static void printOn(ostream&);
};

extern "C" {
void mem_check_ptr(void*);
void mem_check();
void F77_NAME(mem_check)();
}

//operator new with position does not involve MemoryDebugger
void* operator new(size_t,const char*,int);
void* operator new[](size_t,const char*,int);
#define OPERATOR_NEW new(__FILE__,__LINE__)
#define OPERATOR_NEW_BRACKET(T,n) new(__FILE__,__LINE__) T[n]

#endif
