/////////////////////////////////////////////////////////////////////////
//                Calculate the EM fields from Spectators 
//                         in heavy-ion collisions
//
//              author: Chun Shen <shen@mps.ohio-state.edu>
//
//  This code reads in the spectators from initial condition Monte-Carlo 
//  generator and calculates the EM fields for heavy-ion collisions.
//  Note: the output results for the EM fields are actually e*E or e*B 
//        in the unit [MeV^2]
//
//  To do in the future:
/////////////////////////////////////////////////////////////////////////

#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>

#include "Stopwatch.h"
#include "Spectators.h"
#include "EM_fields.h"

using namespace std;

int main()
{
   Stopwatch sw;

   sw.tic();

   Spectators testSp("Spectators.dat");
   EM_fields testEM;
   testEM.calculate_EM_fields(testSp);
   testEM.output_EM_fields_transverse_plane(3, 0);
   testEM.output_EM_fields_time_evolution(5, 5, 0);

   sw.toc();
   cout << "Totally takes " << sw.takeTime() << " sec." << endl;
   return(0);
}

