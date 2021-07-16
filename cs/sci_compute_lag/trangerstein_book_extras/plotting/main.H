#include "VirtualInput.H"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void readMainInput(ifstream &in_file,InputParameterList *main_list,
BOOLEAN &skip_gui) {
//TRACER_CALL(t,"readMainInput");
  BOOLEAN found_main=FALSE;
  char name[LENGTH_NAME], comment[LENGTH_COMMENT];

  istream &infile(in_file);
  infile.clear(ios::goodbit);
  infile.seekg(0,ios::beg);
  while ( in_file >> setw(LENGTH_NAME) >> name ) {
    if ( strcmp(name,"Main") == 0 ) {
      in_file.getline( comment, LENGTH_COMMENT);
      found_main=TRUE;
      while ( in_file >> setw(LENGTH_NAME) >> name ) {
#ifdef DEBUG
//      cout << "\tname = " << name << endl;
#endif
        if ( strcmp(name,"end") == 0 ) break;
        else if (strcmp(name,"skip_gui") == 0) {
	  int iskip_gui; // type bool not read correctly, so use int
	  infile >> iskip_gui;
	  skip_gui= iskip_gui!=0;
	}
        else main_list->formattedRead(in_file,name);
        in_file.getline( comment, LENGTH_COMMENT);
#ifdef DEBUG
//      cout << "\tcomment = " << comment << endl;
#endif
      }
    } else in_file.getline(comment,LENGTH_COMMENT);
  }
#ifdef DEBUG
//main_list->printOn(cout);
#endif
  if ( !found_main ) skip_gui=FALSE;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#ifdef __GNUC__
#pragma implementation "InputParameter.C"
#endif
#include "InputParameter.C"
template class InputParameter<BOOLEAN>;
template class InputParameter<INTEGER>;
template class InputParameter<REAL>;
