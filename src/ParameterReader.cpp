/***********************************************************************
 * Copyright 09-09-2011 Zhi Qiu
***********************************************************************/

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <iomanip>
#include <cstdarg>
#include <vector>

#include "ParameterReader.h"

using namespace std;

//----------------------------------------------------------------------
ParameterReader::ParameterReader() {
    names = new vector<string>;
    values = new vector<double>;
}


//----------------------------------------------------------------------
ParameterReader::~ParameterReader() {
    // in principle garbage collection for vector should be automatic
    // just to be safe
    delete names; 
    delete values;
}


//----------------------------------------------------------------------
string ParameterReader::removeComments(string str, string commentSymbol) {
// Remove comments from a string "str". 
// Comments are all characters after the string "commentSymbol".
    return str.substr(0, str.find(commentSymbol));
}


//----------------------------------------------------------------------
void ParameterReader::phraseEquationWithoutComments(string equation) {
/*
  Phrase an equation like "x=1", and store the result into 
  "names" and "values". 
  The equation is first separated according to the equal sign, 
  then the left and right hand side will be trimmed, 
  after that the right hand side will be converted to double type number.
*/
    if (trim(equation).compare("")==0)
        return;
    size_t symbolPos = equation.find('=');
    if (symbolPos==string::npos) {
        cout << "ParameterReader::phraseEquationWithoutComments error:"
             << "\"=\" symbol not found in equation assignment " 
             << equation << endl;
        exit(-1);
    }
    string LHS(equation.begin(), equation.begin()+symbolPos);
    string RHS(equation.begin()+symbolPos+1, equation.end());
    setVal(LHS, stringToDouble(trim(RHS)));
}


//----------------------------------------------------------------------
long ParameterReader::find(string name) {
/*
  Check if the parameter with "name" already exists in the internal 
  "names" list. If yes, it returns its
*/
    for (unsigned int ii=0; ii<names->size(); ii++)
      if ((*names)[ii].compare(toLower(trim(name)))==0)
          return ii;
    return -1;
}


//----------------------------------------------------------------------
void ParameterReader::phraseOneLine(string str, string commentSymbol) {
/*
  Interpret a string like " x  = 1.1  #bla " to get the associated parameter 
  name and value information, and put them into the internal variables "names" 
  and "values".
*/
    if (trim(str).compare("") == 0)
        return;
    phraseEquationWithoutComments(removeComments(str, commentSymbol));
}


//----------------------------------------------------------------------
void ParameterReader::readFromFile(string filename, string commentSymbol) {
/*
  Read all lines in a file as parameter assignment list. Each line is 
  processed by the phraseOneLine function.
*/
    ifstream parameterFile(filename.c_str());
    if (!parameterFile) {
        cout << "ParameterReader::readFromFile error: file " << filename 
             << " does not exist." << endl;
        exit(-1);
    }
    char buffer[9999];
    while (!parameterFile.eof()) {
        parameterFile.getline(buffer, 9999);
        phraseOneLine(buffer);
    }
    parameterFile.close();
}


//----------------------------------------------------------------------
void ParameterReader::readFromArguments(
        long argc, char * argv[], string commentSymbol, long start_from) {
/*
  Read all strings in argv[]. Each string is processed by the phraseOneLine 
  function.
*/
    for (long ii=start_from; ii<argc; ii++) 
        phraseOneLine(argv[ii], commentSymbol);
}


//----------------------------------------------------------------------
bool ParameterReader::exist(string name) {
/*
  Return true if parameter with "name" is registered.
*/
    return find(name)==-1 ? false: true;
}


//----------------------------------------------------------------------
void ParameterReader::setVal(string name, double value) {
/*
  Set the parameter with "name" to "value". It is appended to the 
  internal "names" and "values" vector if "name" does not exist; 
  otherwise it is rewitten.
*/
    long idx = find(name);
    if (idx==-1) {
        names->push_back(toLower(trim(name))); values->push_back(value);
    } else {
        (*names)[idx]=toLower(trim(name)); (*values)[idx]=value;
    }
}


//----------------------------------------------------------------------
double ParameterReader::getVal(string name) {
/*
  Get the value for the parameter with "name".
*/
    long idx = find(name);
    if (idx!=-1) {
        return (*values)[idx];
    } else {
        cout << "ParameterReader::getVal error: parameter with name " 
             << name << " not found." << endl;
        exit(-1);
    }
}


//----------------------------------------------------------------------
void ParameterReader::echo() {
/*
  Print out all stored parameters to screen.
*/
    if (names->size() == 0) 
        return;
    cout << "=============================================================" 
         << endl;
    cout << "Parameter list: " << endl;
    cout << "=============================================================" 
         << endl;
    for (unsigned int ii = 0; ii < names->size(); ii++)
        cout << (*names)[ii] << "=" << (*values)[ii] << endl;
    cout << "Parameter list end ==========================================" 
         << endl;
}

//**********************************************************************
vector<double> ParameterReader::stringToDoubles(string str) {
    // Return a vector of doubles from the string "str". "str" should
    // be a string containing a line of data.
    // add a blank at the end so the last data will be read
    stringstream sst(str+" "); 
    vector<double> valueList;
    double val;
    sst >> val;
    while (sst.eof()==false)
    {
        valueList.push_back(val);
        sst >> val;
    }
    return valueList;
}


//**********************************************************************
double ParameterReader::stringToDouble(string str) {
    // Return the 1st doubles number read from the string "str". 
    // "str" should be a string containing a line of data.
    // add a blank at the end so the last data will be read
    stringstream sst(str+" "); 
    double val;
    sst >> val;
    return val;
}

//**********************************************************************
vector< vector<double>* >* ParameterReader::readBlockData(istream &stream_in) {
// Return a nested vector of vector<double>* object. Each column of data
// is stored in a vector<double> array and the collection is the returned
// object. Data are read from the input stream "stream_in". Each line
// of data is processed by the stringToDoubles function. Note that the
// data block is dynamicall allocated and is not release within the
// function.
// Note that all "vectors" are "new" so don't forget to delete them.
// Warning that also check if the last line is read correctly. Some files
// are not endded properly and the last line is not read.
    vector< vector<double>* >* data;
    vector<double> valuesInEachLine;
    long lineSize;
    long i; // temp variable
    char buffer[99999]; // each line should be shorter than this

    // first line:
    stream_in.getline(buffer,99999);
    valuesInEachLine = stringToDoubles(buffer);
    // see if it is empty:
    lineSize = valuesInEachLine.size();
    if (lineSize==0) {  // empty:
        cout << "readBlockData warning: input stream has empty first row; "
             << "no data read" << endl;
        return NULL;
    } else {            // not empty; allocate memory:
        data = new vector< vector<double>* >(lineSize);
        for (i=0; i<lineSize; i++) (*data)[i] = new vector<double>;
    }

    // rest of the lines:
    while (stream_in.eof()==false) {
        // set values:
        for (i=0; i<lineSize; i++)
            (*(*data)[i]).push_back(valuesInEachLine[i]);
        // next line:
        stream_in.getline(buffer,99999);
        valuesInEachLine = stringToDoubles(buffer);
    }
    return data;
}

//**********************************************************************
string ParameterReader::toLower(string str) {
// Convert all character in string to lower case
    string tmp = str;
    for (string::iterator it=tmp.begin(); it<=tmp.end(); it++) 
        *it = tolower(*it);
    return tmp;
}

//**********************************************************************
string ParameterReader::trim(string str) {
// Convert all character in string to lower case
    string tmp = str;
    long number_of_char = 0;
    for (size_t ii=0; ii<str.size(); ii++)
        if (str[ii]!=' ' && str[ii]!='\t')
        {
            tmp[number_of_char]=str[ii];
            number_of_char++;
        }
    tmp.resize(number_of_char);
    return tmp;
}

