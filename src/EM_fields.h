// Copyright 2016 Chun Shen
#ifndef SRC_EM_FIELDS_H_
#define SRC_EM_FIELDS_H_

#include <string>
#include <vector>

#include "./ParameterReader.h"

using namespace std;

struct vector3 {
    double x, y, z;
};

struct vector4 {
    double tau, x, y, eta;
};

struct fluidCell {
    double tau, x, y, eta;      // spatial poision of the fluid cell
    vector3 beta;               // flow velocity of the fluid cell
    vector3 E_lab, B_lab;       // E and B fields in the lab frame
    vector4 drift_u;            // drifting 4 velocity induced by EM fields
};

class EM_fields {
 private:
    int mode;
    int turn_on_bulk;
    int initialization_status;
    ParameterReader *paraRdr;

    int nucleon_density_grid_size;
    double nucleon_density_grid_dx;
    // matrices stored the charge density in the transverse plane from the
    // two colliding nuclei, spectators and participants
    double *nucleon_density_grid_x_array, *nucleon_density_grid_y_array;
    int n_eta;
    double* eta_grid;
    double *sinh_eta_array, *cosh_eta_array;
    double **spectator_density_1, **spectator_density_2;
    double **participant_density_1, **participant_density_2;

    // arraies for the space-time points of the EM fields
    int EM_fields_array_length;
    vector<fluidCell> cell_list;

    double charge_fraction;
    double spectator_rap;

 public:
    explicit EM_fields(ParameterReader* paraRdr_in);
    ~EM_fields();

    void set_transverse_grid_points(double tau_local, double eta_local);
    void set_tau_grid_points(double x_local, double y_local, double eta_local);
    void read_in_densities(string path);
    void read_in_spectators_density(string filename_1, string filename_2);
    void read_in_participant_density(string filename_1, string filename_2);
    void read_in_freezeout_surface_points_VISH2p1(string filename1,
                                                  string filename2);
    void calculate_EM_fields();
    void calculate_charge_drifting_velocity();
    void output_EM_fields(string filename);
    void output_surface_file_with_drifting_velocity(string filename);
    void lorentz_transform_vector_in_place(double *u_mu, double *v);
    void Lorentz_boost_EM_fields(double *E_lab, double *B_lab, double *beta,
                                 double *E_prime, double *B_prime);
    void cross_product(double *a, double *b, double *c);
};

#endif  // SRC_EM_FIELDS_H_
