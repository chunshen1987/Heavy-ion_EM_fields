#ifndef SPECTATORS_H
#define SPECTATORS_H

#include<fstream>
#include<vector>

#include "parameter.h"

using namespace std;

typedef struct
{
   double x,y;
   double rapidity_Y;
}Spectator_nucleon;

class Spectators
{
   private:
      ifstream* fp;
      int N_spectators;
      vector<Spectator_nucleon> Spectators_list;
   
   public:
      Spectators(string );
      ~Spectators();

      void ReadinSpectatorslist();
      void printinfo(); //print out the information of the read in file
      int get_Spectator_list_length() {return(N_spectators);};
      double get_Spectator_position_x(int i) {return(Spectators_list[i].x);};
      double get_Spectator_position_y(int i) {return(Spectators_list[i].y);};
      double get_Spectator_rapidity(int i) {return(Spectators_list[i].rapidity_Y);};
};

#endif
