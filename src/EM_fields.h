#ifndef EM_FIELDS_H
#define EM_FIELDS_H

#include "parameter.h"
#include "Spectators.h"

using namespace std;

class EM_fields
{
   private:
      //The value that stored in the matrix are e*EM_fields
      double ****E_x;
      double ****E_y;
      double ****E_z;
      double ****B_x;
      double ****B_y;
      double ****B_z;
      double *t;
      double *x;
      double *y;
      double *z;

      double *r, *phi;
      double *r_weight, *phi_weight;
      double *cos_phi, *sin_phi;

   public:
      EM_fields();
      ~EM_fields();

      void calculate_EM_fields(Spectators );
      void output_EM_fields_transverse_plane(int itime, int iz);
      void output_EM_fields_time_evolution(int ix, int iy, int iz);

};

#endif
