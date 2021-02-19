// Copyright 2016 Chun Shen
#ifndef SRC_EM_FIELDS_H_
#define SRC_EM_FIELDS_H_

#include <string>
#include <vector>

#include "ParameterReader.h"
#include "data_struct.h"

using std::string;

class EM_fields {
 private:
    int debug_flag;
    int mode_;
    int verbose_level;
    int turn_on_bulk;
    int initialization_status;
    int include_participant_contributions;
    ParameterReader *paraRdr;

    int nucleon_density_grid_size;
    double nucleon_density_grid_dx;
    // matrices stored the charge density in the transverse plane from the
    // two colliding nuclei, spectators and participants
    std::vector<double> nucleon_density_grid_x_array_;
    std::vector<double> nucleon_density_grid_y_array_;
    int n_eta;
    std::vector<double> eta_grid_;
    std::vector<double> sinh_eta_array_, cosh_eta_array_;
    double **participant_density_1, **participant_density_2;

    int energy_density_grid_size_;
    double energy_density_grid_dx_;
    std::vector<chargeSource> spectators_1_, spectators_2_;
    std::vector<chargeSource> participants_1_, participants_2_;
    std::vector<double> ed_array_;

    // arraies for the space-time points of the EM fields
    int EM_fields_array_length;
    std::vector<fluidCell> cell_list;
    std::vector<fluidCellSmall> cellListSmall_;

    float hydroEvoHeader[16];

    double charge_fraction;
    double spectator_rap;

 public:
    explicit EM_fields(ParameterReader* paraRdr_in);
    ~EM_fields();

    void set_4d_grid_points();
    void set_tau_grid_points(double x_local, double y_local, double eta_local);
    void read_in_densities(string path);
    void read_in_energy_density(string filename);
    void read_in_spectators_density(string filename_1, string filename_2);
    void read_in_participant_density(string filename_1, string filename_2);
    void read_in_freezeout_surface_points_VISH2p1(string filename1,
                                                  string filename2);
    void read_in_freezeout_surface_points_Gubser(string filename);
    void read_in_freezeout_surface_points_VISH2p1_boost_invariant(
                                                            string filename);
    void read_in_freezeout_surface_points_MUSIC(string filename);
    void read_in_hydro_fluid_cells_MUSIC(string filename);
    void calculate_EM_fields();
    void calculate_EM_fields_no_electric_conductivity();
    void calculate_charge_drifting_velocity();
    void compute_averaged_EM_fields(string filename);
    void output_EM_fields(string filename);
    void output_surface_file_with_drifting_velocity(string filename);
    void output_EMfields_to_hydro_evo(string filename);
    void lorentz_transform_vector_in_place(double *u_mu, double *v);
    void lorentz_transform_vector_with_Lambda(double *u_mu, double *beta);
    void Lorentz_boost_EM_fields(double *E_lab, double *B_lab, double *beta,
                                 double *E_prime, double *B_prime);
    void Lorentz_boost_EM_fields_tensor(double *E_lab, double *B_lab,
                                        double *beta, double *E_prime,
                                        double *B_prime);
    void cross_product(double *a, double *b, double *c);
};

#endif  // SRC_EM_FIELDS_H_
