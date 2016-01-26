#ifndef EM_FIELDS_H
#define EM_FIELDS_H

#include "ParameterReader.h"
#include<string>
#include<vector>

using namespace std;

class EM_fields
{
   private:

      int mode;
      int initialization_status;
      ParameterReader *paraRdr;

      int nucleon_density_grid_size;
      double nucleon_density_grid_dx;
      // matrices stored the charge density in the transverse plane from the 
      // two colliding nuclei, spectators and participants
      double *nucleon_density_grid_x_array, *nucleon_density_grid_y_array;
      int n_eta;
      double* eta_grid;
      double **spectator_density_1, **spectator_density_2;
      double **participant_density_1, **participant_density_2;

      int EM_fields_array_length;
      // The values that stored in the following arraies are e*EM_fields
      vector<double> E_x;
      vector<double> E_y;
      vector<double> E_z;
      vector<double> B_x;
      vector<double> B_y;
      vector<double> B_z;

      // arraies for the space-time points of the EM fields
      vector<double> tau;
      vector<double> x;
      vector<double> y;
      vector<double> eta;

      double charge_fraction;
      double spectator_rap;

   public:
      EM_fields(ParameterReader* paraRdr_in);
      ~EM_fields();

      void set_transverse_grid_points(double tau_local, double eta_local);
      void set_tau_grid_points(double x_local, double y_local, 
                               double eta_local);
      void read_in_densities(string path);
      void read_in_spectators_density(string filename_1, string filename_2);
      void read_in_participant_density(string filename_1, string filename_2);
      void read_in_freezeout_surface_points(string filename);
      void calculate_EM_fields();
      void output_EM_fields(string filename);
};

#endif
