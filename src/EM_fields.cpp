// Copyright 2016 Chun Shen
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>

#include "./parameter.h"
#include "./EM_fields.h"

using namespace std;

EM_fields::EM_fields(ParameterReader* paraRdr_in) {
    initialization_status = 0;
    paraRdr = paraRdr_in;

    mode = paraRdr->getVal("mode");
    int atomic_number = paraRdr->getVal("atomic_number");
    int number_of_proton = paraRdr->getVal("number_of_proton");
    charge_fraction = (static_cast<double>(number_of_proton)
                       /static_cast<double>(atomic_number));

    double ecm = paraRdr->getVal("ecm");
    double gamma = ecm/2./0.938;  // proton mass: 0.938 GeV
    double beta = sqrt(1. - 1./(gamma*gamma));
    double beam_rapidity = atanh(beta);
    spectator_rap = beam_rapidity;
    // cout << "spectator rapidity = " << spectator_rap << endl;

    nucleon_density_grid_size = paraRdr->getVal("nucleon_density_grid_size");
    nucleon_density_grid_dx = paraRdr->getVal("nucleon_density_grid_dx");
    if (nucleon_density_grid_size <= 0) {
        cout << "EM_fields:: Error: Grid size for nucleon density profiles "
             << "needs to be larger than 0!" << endl;
        cout << "Current grid_size = " << nucleon_density_grid_size << endl;
        exit(1);
    }
    nucleon_density_grid_x_array = new double[nucleon_density_grid_size];
    nucleon_density_grid_y_array = new double[nucleon_density_grid_size];
    spectator_density_1 = new double* [nucleon_density_grid_size];
    spectator_density_2 = new double* [nucleon_density_grid_size];
    participant_density_1 = new double* [nucleon_density_grid_size];
    participant_density_2 = new double* [nucleon_density_grid_size];
    for (int i = 0; i < nucleon_density_grid_size; i++) {
        nucleon_density_grid_x_array[i] = (
            -(nucleon_density_grid_size-1)/2. + i*nucleon_density_grid_dx);
        nucleon_density_grid_y_array[i] = (
            -(nucleon_density_grid_size-1)/2. + i*nucleon_density_grid_dx);
        spectator_density_1[i] = new double[nucleon_density_grid_size];
        spectator_density_2[i] = new double[nucleon_density_grid_size];
        participant_density_1[i] = new double[nucleon_density_grid_size];
        participant_density_2[i] = new double[nucleon_density_grid_size];
        for (int j = 0; j < nucleon_density_grid_size; j++) {
            spectator_density_1[i][j] = 0.0;
            spectator_density_2[i][j] = 0.0;
            participant_density_1[i][j] = 0.0;
            participant_density_2[i][j] = 0.0;
        }
    }

    n_eta = paraRdr->getVal("n_eta");
    if (n_eta < 2) {
        cout << "EM_field: Error: n_eta needs to be at least 2" << endl;
        cout << "current n_eta = " << n_eta << endl;
        exit(1);
    }
    eta_grid = new double[n_eta];
    double deta = 2.*beam_rapidity/(n_eta - 1.);
    for (int i = 0; i < n_eta; i++) {
        eta_grid[i] = - beam_rapidity + i*deta;
    }

    read_in_densities("./results");

    if (mode == 0) {
        set_transverse_grid_points(0.2, 0.0);
    } else if (mode == 1) {
        read_in_freezeout_surface_points("./results/surface.dat");
    } else if (mode == 2) {
        set_tau_grid_points(0.0, 0.0, 0.0);
    } else {
        cout << "EM_fields:: Error: unrecognize mode! "
             << "mode = " << mode << endl;
        exit(1);
    }

    initialization_status = 1;
}

EM_fields::~EM_fields() {
    if (initialization_status == 1) {
        for (int i = 0; i < nucleon_density_grid_size; i++) {
            delete[] spectator_density_1[i];
            delete[] spectator_density_2[i];
            delete[] participant_density_1[i];
            delete[] participant_density_2[i];
        }
        delete[] spectator_density_1;
        delete[] spectator_density_2;
        delete[] participant_density_1;
        delete[] participant_density_2;
        delete[] nucleon_density_grid_x_array;
        delete[] nucleon_density_grid_y_array;

        delete[] eta_grid;

        E_x.clear();
        E_y.clear();
        E_z.clear();
        B_x.clear();
        B_y.clear();
        B_z.clear();

        tau.clear();
        x.clear();
        y.clear();
        eta.clear();
    }
    return;
}

void EM_fields::read_in_densities(string path) {
    // spectators
    ostringstream spectator_1_filename;
    spectator_1_filename << path << "/spectator_density_A_disk.dat";
    // spectator_1_filename << path
    //                      << "/spectator_density_A_fromSd_order_2.dat";
    ostringstream spectator_2_filename;
    spectator_2_filename << path << "/spectator_density_B_disk.dat";
    // spectator_2_filename << path
    //                      << "/spectator_density_B_fromSd_order_2.dat";
    read_in_spectators_density(spectator_1_filename.str(),
                               spectator_2_filename.str());
    // participants
    //ostringstream participant_1_filename;
    //participant_1_filename << path 
    //                       << "/nuclear_thickness_TA_fromSd_order_2.dat";
    //ostringstream participant_2_filename;
    //participant_2_filename << path 
    //                       << "/nuclear_thickness_TB_fromSd_order_2.dat";
    //read_in_participant_density(participant_1_filename.str(), 
    //                            participant_2_filename.str());
}

void EM_fields::read_in_spectators_density(string filename_1,
                                           string filename_2) {
    cout << "read in spectator density ...";
    ifstream spec1(filename_1.c_str());
    ifstream spec2(filename_2.c_str());

    for (int i = 0; i < nucleon_density_grid_size; i++) {
        for (int j = 0; j < nucleon_density_grid_size; j++) {
            spec1 >> spectator_density_1[i][j];
            spec2 >> spectator_density_2[i][j];
        }
    }

    spec1.close();
    spec2.close();
    cout << " done!" << endl;
}

void EM_fields::read_in_participant_density(string filename_1,
                                            string filename_2) {
    cout << "read in participant density ...";
    ifstream part1(filename_1.c_str());
    ifstream part2(filename_2.c_str());

    for (int i = 0; i < nucleon_density_grid_size; i++) {
        for (int j = 0; j < nucleon_density_grid_size; j++) {
            part1 >> participant_density_1[i][j];
            part2 >> participant_density_2[i][j];
        }
    }

    part1.close();
    part2.close();
    cout << " done!" << endl;
}

void EM_fields::set_tau_grid_points(double x_local, double y_local,
                                    double eta_local) {
    double EM_fields_grid_size = 15.0;
    double EM_fields_grid_dtau = 0.01;

    int number_of_points =
                static_cast<int>(EM_fields_grid_size/EM_fields_grid_dtau) + 1;
    for (int i = 0; i < number_of_points; i++) {
        double tau_local = 0.0 + i*EM_fields_grid_dtau;
        tau.push_back(tau_local);
        x.push_back(x_local);
        y.push_back(y_local);
        eta.push_back(eta_local);
    }
    EM_fields_array_length = tau.size();
    cout << "number of freeze-out cells: " << EM_fields_array_length << endl;
}

void EM_fields::set_transverse_grid_points(double tau_local, double eta_local) {
    double EM_fields_grid_size = 20.0;
    double EM_fields_grid_dx = 0.1;

    int number_of_points =
                static_cast<int>(EM_fields_grid_size/EM_fields_grid_dx) + 1;
    for (int i = 0; i < number_of_points; i++) {
        double x_local = - EM_fields_grid_size/2. + i*EM_fields_grid_dx;
        for (int j = 0; j < number_of_points; j++) {
            double y_local = - EM_fields_grid_size/2. + j*EM_fields_grid_dx;
            x.push_back(x_local);
            y.push_back(y_local);
            tau.push_back(tau_local);
            eta.push_back(eta_local);
        }
    }
    EM_fields_array_length = x.size();
    cout << "number of freeze-out cells: " << EM_fields_array_length << endl;
}

void EM_fields::read_in_freezeout_surface_points(string filename) {
    ifstream FOsurf(filename.c_str());

    cout << "read in freeze-out surface points ...";
    // read in freeze-out surface positions
    double dummy;
    double tau_local, x_local, y_local, eta_local;
    FOsurf >> dummy;
    while (!FOsurf.eof()) {
        FOsurf >> tau_local >> x_local >> y_local >> dummy >> dummy >> dummy;
        for (int i = 0; i < n_eta; i++) {
            eta_local = eta_grid[i];
            tau.push_back(tau_local);
            x.push_back(x_local);
            y.push_back(y_local);
            eta.push_back(eta_local);
        }
        FOsurf >> dummy;
    }
    FOsurf.close();
    cout << " done!" << endl;
    EM_fields_array_length = x.size();
    cout << "number of freeze-out cells: " << EM_fields_array_length << endl;
}

void EM_fields::calculate_EM_fields() {
    double cosh_spectator_rap = cosh(spectator_rap);
    double sinh_spectator_rap = sinh(spectator_rap);

    double dx_sq = nucleon_density_grid_dx*nucleon_density_grid_dx;
    for (int i_array = 0; i_array < EM_fields_array_length; i_array++) {
        double field_x = x[i_array];
        double field_y = y[i_array];
        double field_tau = tau[i_array];
        double field_eta = eta[i_array];
        double temp_sum_Ex_spectator = 0.0e0;
        double temp_sum_Ey_spectator = 0.0e0;
        double temp_sum_Ez_spectator = 0.0e0;
        double temp_sum_Bx_spectator = 0.0e0;
        double temp_sum_By_spectator = 0.0e0;

        double z_local_spectator_1 = field_tau*sinh(field_eta - spectator_rap);
        double z_local_spectator_2 = field_tau*sinh(field_eta + spectator_rap);
        double z_local_spectator_1_sq = z_local_spectator_1*z_local_spectator_1;
        double z_local_spectator_2_sq = z_local_spectator_2*z_local_spectator_2;

        for (int i = 0; i < nucleon_density_grid_size; i++) {
            double grid_x = nucleon_density_grid_x_array[i];
            for (int j = 0; j < nucleon_density_grid_size; j++) {
                double grid_y = nucleon_density_grid_y_array[j];

                double x_local = field_x - grid_x;
                double y_local = field_y - grid_y;
                double r_perp_local_sq = x_local*x_local + y_local*y_local;
                double r_spectator_1 = sqrt(
                                r_perp_local_sq + z_local_spectator_1_sq);
                double r_cubic_spectator_1 = (r_spectator_1*r_spectator_1
                                              *r_spectator_1);
                double r_spectator_2 = sqrt(
                                r_perp_local_sq + z_local_spectator_2_sq);
                double r_cubic_spectator_2 = (r_spectator_2*r_spectator_2
                                              *r_spectator_2);
                double spectator_integrand_1 = (
                                spectator_density_1[i][j]/r_cubic_spectator_1);
                double spectator_integrand_2 = (
                                spectator_density_2[i][j]/r_cubic_spectator_2);

                double Ex_spectator_integrand = x_local*(
                                spectator_integrand_1 + spectator_integrand_2);
                double Ey_spectator_integrand = y_local*(
                                spectator_integrand_1 + spectator_integrand_2);
                double Ez_spectator_integrand = (
                                  z_local_spectator_1*spectator_integrand_1
                                + z_local_spectator_2*spectator_integrand_2);
                double Bx_spectator_integrand = y_local*(
                                spectator_integrand_1 - spectator_integrand_2);
                double By_spectator_integrand = x_local*(
                                spectator_integrand_1 - spectator_integrand_2);
                temp_sum_Ex_spectator += Ex_spectator_integrand;
                temp_sum_Ey_spectator += Ey_spectator_integrand;
                temp_sum_Ez_spectator += Ez_spectator_integrand;
                temp_sum_Bx_spectator += Bx_spectator_integrand;
                temp_sum_By_spectator += By_spectator_integrand;
            }
        }
        E_x.push_back(unit_convert*charge_fraction*alpha_EM*(
                        cosh_spectator_rap*temp_sum_Ex_spectator)*dx_sq);
        E_y.push_back(unit_convert*charge_fraction*alpha_EM*(
                        cosh_spectator_rap*temp_sum_Ey_spectator)*dx_sq);
        E_z.push_back(unit_convert*charge_fraction*alpha_EM*(
                        temp_sum_Ez_spectator)*dx_sq);
        B_x.push_back(unit_convert*charge_fraction*alpha_EM*(
                        (-sinh_spectator_rap)*temp_sum_Bx_spectator)*dx_sq);
        B_y.push_back(unit_convert*charge_fraction*alpha_EM*(
                        sinh_spectator_rap*temp_sum_By_spectator)*dx_sq);
        B_z.push_back(0.0);
        if (i_array % static_cast<int>(EM_fields_array_length/10) == 0) {
            cout << "computing EM fields: " << setprecision(3)
                 << (static_cast<double>(i_array)
                     /static_cast<double>(EM_fields_array_length)*100)
                 << "\% done." << endl;
        }
    }
    return;
}

void EM_fields::output_EM_fields(string filename) {
    ofstream output_file(filename.c_str());

    // write a header first
    output_file << "# tau[fm]  x[fm]  y[fm]  eta  "
                << "eE_x[GeV^2]  eE_y[GeV^2]  eE_z[GeV^2]  "
                << "eB_x[GeV^2]  eB_y[GeV^2]  eB_z[GeV^2]" << endl;
    for (int i = 0; i < EM_fields_array_length; i++) {
        output_file << scientific << setprecision(8) << setw(15)
                    << tau[i] << "   " << x[i] << "   " << y[i] << "   "
                    << eta[i] << "   "
                    << E_x[i] << "   " << E_y[i] << "   " << E_z[i] << "   "
                    << B_x[i] << "   " << B_y[i] << "   " << B_z[i] << endl;
    }
    output_file.close();
    return;
}
