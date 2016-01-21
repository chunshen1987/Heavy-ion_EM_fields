#ifndef EM_FIELDS_H
#define EM_FIELDS_H

#include "parameter.h"
#include "Spectators.h"
#include<string>

using namespace std;

class EM_fields
{
   private:

      // matrices stored the charge density in the transverse plane from the 
      // two colliding nuclei, spectators and participants
      double **spetator_density_1, **spetator_density_2;
      double **participant_density_1, **participant_density_2;

      // The values that stored in the following arraies are e*EM_fields
      double *E_x;
      double *E_y;
      double *E_z;
      double *B_x;
      double *B_y;
      double *B_z;

      // arraies for the space-time points of the EM fields
      double *tau;
      double *x;
      double *y;
      double *eta;

   public:
      EM_fields();
      ~EM_fields();

      void calculate_EM_fields();
      void output_EM_fields(string filename);
};

#endif
