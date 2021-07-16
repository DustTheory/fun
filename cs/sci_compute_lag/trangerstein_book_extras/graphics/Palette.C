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
// "$Header: /home/faculty/johnt/cvs/deal_new/graphics/Palette.C,v 1.1 2009/08/20 17:31:46 johnt Exp $"
#include "Palette.H"
#include <iomanip>
#include <string.h>
#include <stdlib.h>
#include "Const.H"
//#include "MemoryDebugger.H"
# ifndef OPERATOR_NEW
# define OPERATOR_NEW new
# endif
# ifndef OPERATOR_NEW_BRACKET
# define OPERATOR_NEW_BRACKET(T,n) new T[n]
# endif
#define DEFAULT_PALETTE_SIZE 256
#include "ISLList.C"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Palette::Palette():max_colors(DEFAULT_PALETTE_SIZE),num_named_colors(0){
  name=OPERATOR_NEW_BRACKET(char,8);
  strcpy(name,"default");
  map_index=OPERATOR_NEW_BRACKET(int,max_colors);
  color_name=OPERATOR_NEW_BRACKET(char*,max_colors);
#ifdef INDEF
  for (int i=0;i<max_colors;i++) {
    map_index[i]=0;
    color_name[i]=0;
  }
#endif
  defaultColors();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Palette::Palette(const char *na) : max_colors(DEFAULT_PALETTE_SIZE),
num_named_colors(0) {
  int length_name=strlen(na)+1;
  name=OPERATOR_NEW_BRACKET(char,length_name);
  strcpy(name,na);
  map_index=OPERATOR_NEW_BRACKET(int,max_colors);
  color_name=OPERATOR_NEW_BRACKET(char*,max_colors);
#ifdef INDEF
  for (int i=0;i<max_colors;i++) {
    map_index[i]=0;
    color_name[i]=0;
  }
#endif
  if (strcmp(name,"default")==0) defaultColors();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Palette::~Palette() {
  for (int i=0;i<num_named_colors;i++) {
    delete [] color_name[i];
  }
  delete [] color_name; color_name=0;
  delete [] map_index; map_index=0;
  delete [] name; name=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
bool Palette::isOk() const {
  if (num_named_colors<3) return false;
  if (map_index[0]!=0) return false;
  if (map_index[1]!=1) return false;
  for (int i=2;i<num_named_colors;i++) {
    if (map_index[i-1]>=map_index[i]) return false;
  }
  return true;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Palette::insertEntry(int index,const char* na) {
  CHECK_TEST(index>=0);
  CHECK_TEST(num_named_colors<max_colors);

  int i;
  if (num_named_colors==0 || map_index[num_named_colors-1]<index) { 
    i=num_named_colors;
  } else {
    for (i=0;i<num_named_colors;i++) if (map_index[i]>index) break;
    for (int j=num_named_colors;j>i;j--) {
      map_index[j]=map_index[j-1];
      color_name[j]=color_name[j-1];
    }
  }
  map_index[i]=index;
  int length_name=strlen(na)+1;
  color_name[i]=OPERATOR_NEW_BRACKET(char,length_name);
  strcpy(color_name[i],na);
  num_named_colors++;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int Palette::findEntry(const char* na) const {
  for (int i=0;i<num_named_colors;i++) {
    if (strcmp(na,color_name[i])==0) return i;
  }
  return -1;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Palette::defaultColors() {
  int i=0;
#if (SPACEDIM==1)
  insertEntry(i++,"white");
  insertEntry(i++,"red");
  insertEntry(i++,"green");
  insertEntry(i++,"blue");
  insertEntry(i++,"yellow");
  insertEntry(i++,"cyan");
  insertEntry(i++,"magenta");
  insertEntry(i++,"black");
#else
  insertEntry(i++,"white");
//insertEntry(i++,"gainsboro");
//insertEntry(i++,"LightSlateGrey");
//insertEntry(i++,"SlateGrey");
//insertEntry(i++,"NavyBlue");
//insertEntry(i++,"SteelBlue");
//insertEntry(i++,"RoyalBlue");
//insertEntry(i++,"blue");
//insertEntry(i++,"CornflowerBlue");
//insertEntry(i++,"SkyBlue");
//insertEntry(i++,"cyan");
//insertEntry(i++,"aquamarine");
//insertEntry(i++,"MediumSpringGreen");
//insertEntry(i++,"green");
//insertEntry(i++,"GreenYellow");
//insertEntry(i++,"yellow");
//insertEntry(i++,"gold");
//insertEntry(i++,"LightSalmon");
//insertEntry(i++,"orange");
//insertEntry(i++,"DarkOrange");
//insertEntry(i++,"tomato");
//insertEntry(i++,"red");
//insertEntry(i++,"DeepPink");
//insertEntry(i++,"PaleVioletRed");
//insertEntry(i++,"VioletRed");
//insertEntry(i++,"magenta");
//insertEntry(i++,"violet");
//insertEntry(i++,"DarkOrchid");
//insertEntry(i++,"BlueViolet");
//insertEntry(i++,"purple");
  insertEntry(1,"black");
  insertEntry(2,"blue");
  insertEntry(19,"cyan");
  insertEntry(39,"green");
  insertEntry(59,"yellow");
  insertEntry(79,"red");
//insertEntry(99,"magenta");
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Palette::printOn(ostream &os) const {
  os << "Palette: " << endl;
  os << "\nname = " << name 
     << "\nnum_named_colors = " << num_named_colors << endl;
  for (int i=0;i<num_named_colors;i++) {
     os << "map_index = " << map_index[i]
        << "  color_name = " << color_name[i] << endl;
  }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
PaletteList::~PaletteList() {
  while (notEmpty()) delete delAfter(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void PaletteList::makePalette(const char *palette_name,
ifstream &in_file,Palette *&palette) {
  char mycomment[LENGTH_COMMENT];
  palette=OPERATOR_NEW Palette(palette_name);
  if (palette_name!=0 && strcmp(palette_name,"default")!=0) {
    in_file.getline(mycomment,LENGTH_COMMENT);
    char name[LENGTH_NAME];
    while (in_file >> setw(sizeof(name)) >> name) {
      if (strcmp(name,"entry")==0) {
        int map_index;
        in_file >> map_index >> setw(sizeof(name)) >> name;
        palette->insertEntry(map_index,name);
        in_file.getline(mycomment,LENGTH_COMMENT);
      } else if (strcmp(name,"end")==0) break;
    }
    in_file.getline(mycomment,LENGTH_COMMENT);
  }
  append(*palette);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
const Palette* PaletteList::findByName(const char *na) const {
  for (Palette *p=first();p;p=next(p))
    if (strcmp(na,p->getName())==0) return p;
  return 0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void PaletteList::printOn(ostream &os) const {
  os << "PaletteList: this = " << this << endl;
  if (!first()) os << "list is empty" << endl;
  else {
    int i=0;
    for (Palette *p=first();p;p=next(p),i++) {
       os << i << " : Palette = " << p << endl;
    }
  }
}

template class ISLList<Palette>;
