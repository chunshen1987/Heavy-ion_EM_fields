#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>

#include "parameter.h"
#include "EM_fields.h"

using namespace std;

EM_fields::EM_fields()
{
}

EM_fields::~EM_fields()
{
   delete[] E_x;
   delete[] E_y;
   delete[] E_z;
   delete[] B_x;
   delete[] B_y;
   delete[] B_z;

   delete[] tau;
   delete[] x;
   delete[] y;
   delete[] eta;

   return;
}

EM_fields::read_in_spectators_density(string filename_1, string filename_2)
{

}

EM_fields::read_in_participant_density(string filename_1, string filename_2)
{

}

void EM_fields::calculate_EM_fields(Spectators nucleon_list)
{
   for(int i=0; i<nt; i++)
      for(int j=0; j<nx; j++)
         for(int k=0; k<ny; k++)
            for(int l=0; l<nz; l++)
            {
               double temp_sum_Ex = 0.0e0;
               double temp_sum_Ey = 0.0e0;
               double temp_sum_Ez = 0.0e0;
               double temp_sum_Bx = 0.0e0;
               double temp_sum_By = 0.0e0;
               double temp_sum_Bz = 0.0e0;
               for(int list_idx=0; list_idx < list_length; list_idx++)
               {
                  double x_local, y_local, z_local;
                  x_local = x[j] - nucleon_position_x[list_idx];
                  y_local = y[k] - nucleon_position_y[list_idx];
                  z_local = z[l]*nucleon_cosh_y[list_idx] - t[i]*nucleon_sinh_y[list_idx];
                  double denominator, denominator_cubic;
                  double Ex_integrand, Ey_integrand, Ez_integrand, Bx_integrand, By_integrand, Bz_integrand;
                  double Ex_integral, Ey_integral, Ez_integral, Bx_integral, By_integral, Bz_integral;
                  switch(Nucleon_type)
                  {
                     case 0: //point-like nucleons
                       denominator = sqrt(x_local*x_local + y_local*y_local + z_local*z_local);
                       denominator_cubic = denominator*denominator*denominator;
                       Ex_integral = nucleon_cosh_y[list_idx]*x_local/denominator_cubic;
                       Ey_integral = nucleon_cosh_y[list_idx]*y_local/denominator_cubic;
                       Ez_integral = z_local/denominator_cubic;
                       Bx_integral = nucleon_sinh_y[list_idx]*y_local/denominator_cubic;
                       By_integral = - nucleon_sinh_y[list_idx]*x_local/denominator_cubic;
                       Bz_integral = 0.0e0;
                       break;
                     case 1: //disk-like nucleons
                       for(int i_r=0; i_r<n_r; i_r++)
                          for(int i_phi=0; i_phi<n_phi; i_phi++)
                          {
                             denominator = sqrt(x_local*x_local + y_local*y_local + r[i_r]*r[i_r] - 2*r[i_r]*(cos_phi[i_phi]*x_local+sin_phi[i_phi]*y_local) + z_local*z_local);
                             denominator_cubic = denominator*denominator*denominator;
                             Ex_integrand = norm_disk*nucleon_cosh_y[list_idx]*(x_local - r[i_r]*cos_phi[i_phi])/denominator_cubic;
                             Ey_integrand = norm_disk*nucleon_cosh_y[list_idx]*(y_local - r[i_r]*sin_phi[i_phi])/denominator_cubic;
                             Ez_integrand = norm_disk*z_local/denominator_cubic;
                             Bx_integrand = norm_disk*nucleon_sinh_y[list_idx]*(y_local - r[i_r]*sin_phi[i_phi])/denominator_cubic;
                             By_integrand = - norm_disk*nucleon_sinh_y[list_idx]*(x_local - r[i_r]*cos_phi[i_phi])/denominator_cubic;
                             Bz_integrand = 0.0e0;
                             
                             Ex_integral += Ex_integrand*phi_weight[i_phi]*r_weight[i_r]*r[i_r];
                             Ey_integral += Ey_integrand*phi_weight[i_phi]*r_weight[i_r]*r[i_r];
                             Ez_integral += Ez_integrand*phi_weight[i_phi]*r_weight[i_r]*r[i_r];
                             Bx_integral += Bx_integrand*phi_weight[i_phi]*r_weight[i_r]*r[i_r];
                             By_integral += By_integrand*phi_weight[i_phi]*r_weight[i_r]*r[i_r];
                             Bz_integral += Bz_integrand*phi_weight[i_phi]*r_weight[i_r]*r[i_r];
                          }
                       break;
                     case 2: //gaussian-like nucleons
                       for(int i_r=0; i_r<n_r; i_r++)
                          for(int i_phi=0; i_phi<n_phi; i_phi++)
                          {
                             double expon = exp( - r[i_r]*r[i_r]/(2*sigma_gaussian*sigma_gaussian));
                             denominator = sqrt(x_local*x_local + y_local*y_local + r[i_r]*r[i_r] - 2*r[i_r]*(cos_phi[i_phi]*x_local+sin_phi[i_phi]*y_local) + z_local*z_local);
                             denominator_cubic = denominator*denominator*denominator;
                             Ex_integrand = norm_gaussian*nucleon_cosh_y[list_idx]*expon*(x_local - r[i_r]*cos_phi[i_phi])/denominator_cubic;
                             Ey_integrand = norm_gaussian*nucleon_cosh_y[list_idx]*expon*(y_local - r[i_r]*sin_phi[i_phi])/denominator_cubic;
                             Ez_integrand = norm_gaussian*expon*z_local/denominator_cubic;
                             Bx_integrand = norm_gaussian*expon*nucleon_sinh_y[list_idx]*(y_local - r[i_r]*sin_phi[i_phi])/denominator_cubic;
                             By_integrand = - norm_gaussian*expon*nucleon_sinh_y[list_idx]*(x_local - r[i_r]*cos_phi[i_phi])/denominator_cubic;
                             Bz_integrand = 0.0e0;
                             
                             Ex_integral += Ex_integrand*phi_weight[i_phi]*r_weight[i_r]*r[i_r];
                             Ey_integral += Ey_integrand*phi_weight[i_phi]*r_weight[i_r]*r[i_r];
                             Ez_integral += Ez_integrand*phi_weight[i_phi]*r_weight[i_r]*r[i_r];
                             Bx_integral += Bx_integrand*phi_weight[i_phi]*r_weight[i_r]*r[i_r];
                             By_integral += By_integrand*phi_weight[i_phi]*r_weight[i_r]*r[i_r];
                             Bz_integral += Bz_integrand*phi_weight[i_phi]*r_weight[i_r]*r[i_r];
                          }
                       break;
                     default:
                       cout << "EM_fields::calculate_EM_fields error: wrong type of nucleons, Nucleon_type = " << Nucleon_type << endl;
                       exit(1);
                       break;
                  }
                  temp_sum_Ex += e_sq*Ex_integral;
                  temp_sum_Ey += e_sq*Ey_integral;
                  temp_sum_Ez += e_sq*Ez_integral;
                  temp_sum_Bx += e_sq*Bx_integral;
                  temp_sum_By += e_sq*By_integral;
                  temp_sum_Bz += e_sq*Bz_integral;
                  Ex_integral = 0.0;
                  Ey_integral = 0.0;
                  Ez_integral = 0.0;
                  Bx_integral = 0.0;
                  By_integral = 0.0;
                  Bz_integral = 0.0;

               }
               E_x[i][j][k][l] = temp_sum_Ex*unit_convert;   // unit of [e*E_x] is MeV^2
               E_y[i][j][k][l] = temp_sum_Ey*unit_convert;
               E_z[i][j][k][l] = temp_sum_Ez*unit_convert;
               B_x[i][j][k][l] = temp_sum_Bx*unit_convert;
               B_y[i][j][k][l] = temp_sum_By*unit_convert;
               B_z[i][j][k][l] = temp_sum_Bz*unit_convert;
            }
   return;
}

void EM_fields::output_EM_fields(string filename)
{
   ofstream output_file(filename.c_str());
   for(int i = 0; i < n_array; i++)
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
