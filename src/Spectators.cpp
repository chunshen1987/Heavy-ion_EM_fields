#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>

#include "Spectators.h"

using namespace std;

Spectators::Spectators(string filename)
{
   ostringstream infilename_stream;
   infilename_stream << results_path << filename;
   fp = new ifstream;
   fp->open(infilename_stream.str().c_str());
   if(fp->is_open() == false)
   {
      cout << "Spectators::ReadinSpectators error: data file can not be found!" << endl;
      cout << "data filename: " << infilename_stream.str().c_str() << endl;
      exit(1);

   }
   ReadinSpectatorslist();
   printinfo();
   return;
}

Spectators::~Spectators()
{
}

void Spectators::ReadinSpectatorslist()
{
   Spectator_nucleon temp;
   (*fp) >> temp.x >> temp.y >> temp.rapidity_Y ;
   while( !fp->eof() )
   {
      Spectators_list.push_back(temp);
      (*fp) >> temp.x >> temp.y >> temp.rapidity_Y;
   }
   N_spectators = Spectators_list.size();
   fp->close();
   return;
}

void Spectators::printinfo()
{
   cout << "------------------------------------------------------------" << endl;
   cout << "---Spectators information : " << endl;
   cout << "------------------------------------------------------------" << endl;
   cout << "There are totally " << N_spectators << " spectator nucleons." << endl;
   cout << "Spectators' positions are : " << endl;
   for(int i=0; i<N_spectators; i++)
   {
      cout << i << " : " 
           << " x = " 
           << setprecision(4) << setw(6) << Spectators_list[i].x << " fm,"
           << " y = " 
           << setprecision(4) << setw(6) << Spectators_list[i].y << " fm,"
           << " rapidity = "
           << setprecision(4) << setw(6) << Spectators_list[i].rapidity_Y 
           << endl;
   }
   cout << "************************************************************" << endl;
   return;
}
