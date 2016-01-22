#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>

#include "parameter.h"
#include "EM_fields.h"

using namespace std;

EM_fields::EM_fields(ParameterReader* paraRdr_in)
{
    initialization_status = 0;
    paraRdr = paraRdr_in;
    int atomic_number = paraRdr->getVal("atomic_number");
    int number_of_proton = paraRdr->getVal("number_of_proton");
    charge_fraction = (double)number_of_proton/(double)atomic_number;
    
    double ecm = paraRdr->getVal("ecm");
    double gamma = ecm/2./0.938;  // proton mass: 0.938 GeV
    double beta = sqrt(1. - 1./(gamma*gamma));
    double beam_rapidity = atanh(beta);
    spectator_rap = beam_rapidity;
    //cout << "spectator rapidity = " << spectator_rap << endl;
    participant_rap = 0.0;

    nucleon_density_grid_size = paraRdr->getVal("nucleon_density_grid_size");
    nucleon_density_grid_dx = paraRdr->getVal("nucleon_density_grid_dx");
    if(nucleon_density_grid_size <= 0)
    {
        cout << "EM_fields:: Error: Grid size for nucleon density profiles " 
             << "needs to be larger than 0!" << endl;
        cout << "Current grid_size = " << nucleon_density_grid_size << endl;
        exit(1);
    }
    nucleon_density_grid_x_array = new double [nucleon_density_grid_size];
    nucleon_density_grid_y_array = new double [nucleon_density_grid_size];
    spectator_density_1 = new double* [nucleon_density_grid_size];
    spectator_density_2 = new double* [nucleon_density_grid_size];
    participant_density_1 = new double* [nucleon_density_grid_size];
    participant_density_2 = new double* [nucleon_density_grid_size];
    for(int i = 0; i < nucleon_density_grid_size; i++)
    {
        nucleon_density_grid_x_array[i] = (
            -(nucleon_density_grid_size-1)/2. + i*nucleon_density_grid_dx);
        nucleon_density_grid_y_array[i] = (
            -(nucleon_density_grid_size-1)/2. + i*nucleon_density_grid_dx);
        spectator_density_1[i] = new double [nucleon_density_grid_size];
        spectator_density_2[i] = new double [nucleon_density_grid_size];
        participant_density_1[i] = new double [nucleon_density_grid_size];
        participant_density_2[i] = new double [nucleon_density_grid_size];
        for(int j = 0; j < nucleon_density_grid_size; j++)
        {
            spectator_density_1[i][j] = 0.0;
            spectator_density_2[i][j] = 0.0;
            participant_density_1[i][j] = 0.0;
            participant_density_2[i][j] = 0.0;
        }
    }

    read_in_densities("./results");
    initialization_status = 1;
}

EM_fields::~EM_fields()
{
    if(initialization_status == 1)
    {
        for(int i = 0; i < nucleon_density_grid_size; i++)
        {
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

void EM_fields::read_in_densities(string path)
{
    // spectators
    ostringstream spectator_1_filename;
    spectator_1_filename << path 
                         << "/spectator_density_A_fromSd_order_2.dat";
    ostringstream spectator_2_filename;
    spectator_2_filename << path 
                         << "/spectator_density_B_fromSd_order_2.dat";
    read_in_spectators_density(spectator_1_filename.str(), 
                               spectator_2_filename.str());
    // participants
    ostringstream participant_1_filename;
    participant_1_filename << path 
                           << "/nuclear_thickness_TA_fromSd_order_2.dat";
    ostringstream participant_2_filename;
    participant_2_filename << path 
                           << "/nuclear_thickness_TB_fromSd_order_2.dat";
    read_in_participant_density(participant_1_filename.str(), 
                                participant_2_filename.str());
}

void EM_fields::read_in_spectators_density(string filename_1, string filename_2)
{
    ifstream spec1(filename_1.c_str());
    ifstream spec2(filename_2.c_str());

    for(int i = 0; i < nucleon_density_grid_size; i++)
    {
        for(int j = 0; j < nucleon_density_grid_size; j++)
        {
            spec1 >> spectator_density_1[i][j];
            spec2 >> spectator_density_2[i][j];
        }
    }

    spec1.close();
    spec2.close();
}

void EM_fields::read_in_participant_density(string filename_1, string filename_2)
{
    ifstream part1(filename_1.c_str());
    ifstream part2(filename_2.c_str());

    for(int i = 0; i < nucleon_density_grid_size; i++)
    {
        for(int j = 0; j < nucleon_density_grid_size; j++)
        {
            part1 >> participant_density_1[i][j];
            part2 >> participant_density_2[i][j];
        }
    }

    part1.close();
    part2.close();
}

void EM_fields::read_in_freezeout_surface_points(string filename)
{
    EM_fields_array_length = x.size();
}

void EM_fields::calculate_EM_fields()
{
    double cosh_spectator_rap = cosh(spectator_rap);
    double sinh_spectator_rap = sinh(spectator_rap);
    for(int i_array = 0; i_array < EM_fields_array_length; i_array++)
    {
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

        for(int i = 0; i < nucleon_density_grid_size; i++)
        {
            double grid_x = nucleon_density_grid_x_array[i];
            for(int j = 0; j < nucleon_density_grid_size; j++)
            {
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
        E_x[i_array] = unit_convert*charge_fraction*e_sq*(
                        cosh_spectator_rap*temp_sum_Ex_spectator);
        E_y[i_array] = unit_convert*charge_fraction*e_sq*(
                        cosh_spectator_rap*temp_sum_Ey_spectator);
        E_z[i_array] = unit_convert*charge_fraction*e_sq*(temp_sum_Ez_spectator);
        B_x[i_array] = unit_convert*charge_fraction*e_sq*(
                        sinh_spectator_rap*temp_sum_Bx_spectator);
        B_y[i_array] = unit_convert*charge_fraction*e_sq*(
                        (-sinh_spectator_rap)*temp_sum_By_spectator);
        B_z[i_array] = 0.0;
    }
    return;
}

void EM_fields::output_EM_fields(string filename)
{
   ofstream output_file(filename.c_str());
   for(int i = 0; i < EM_fields_array_length; i++)
   {
      output_file << scientific << setprecision(8) << setw(15)  
                  << tau[i] << "   " << x[i] << "   " << y[i] << "   "
                  << eta[i] << "   " 
                  << E_x[i] << "   " << E_y[i] << "   " << E_z[i] << "   " 
                  << B_x[i] << "   " << B_y[i] << "   " << B_z[i] << endl;
   }
   output_file.close();
   return;
}
