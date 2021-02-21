/////////////////////////////////////////////////////////////////////////
//                Calculate the EM fields from Spectators
//                         in heavy-ion collisions
//
//              author: Chun Shen <chunshen@physics.mcgill.ca>
//              Copyright 2016 Chun Shen
//
//  This code reads in the spectator and participant densities from
//  initial condition Monte-Carlo generator and calculates the EM fields
//  for heavy-ion collisions at the given space-time points.
//
//  Note: the output results for the EM fields are actually e*E or e*B
//        in the unit [MeV^2]
//
/////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>

#include "./Stopwatch.h"
#include "./EM_fields.h"
#include "./ParameterReader.h"

using namespace std;

int main(int argc, char *argv[]) {
    Stopwatch sw;

    sw.tic();

    // Read-in parameters
    ParameterReader paraRdr;
    paraRdr.readFromFile("parameters.dat");
    paraRdr.readFromArguments(argc, argv);
    paraRdr.echo();
    int mode_ = paraRdr.getVal("mode");

    EM_fields testEM(&paraRdr);
    testEM.calculate_EM_fields();
    if (mode_ == 10) {
        testEM.compute_averaged_EM_fields("./results/EMavg.dat");
    } else if (mode_ == 7) {
        testEM.output_EMfields_to_hydro_evo(
                "./results/evolution_all_xyeta_withEMfields.dat");
    } else {
        testEM.output_EM_fields("./results/EM_fields.dat");
        if (mode_ != 0) {
            testEM.calculate_charge_drifting_velocity();
            testEM.output_surface_file_with_drifting_velocity(
                            "./results/surface_with_drifting_velocity.dat");
        }
    }

    sw.toc();
    cout << "Totally takes " << sw.takeTime() << " sec." << endl;
    return(0);
}

