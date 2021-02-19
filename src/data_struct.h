// Copyright 2021 Chun Shen
#ifndef DATA_STRUCT_H_
#define DATA_STRUCT_H_

struct vector3 {
    double x, y, z;
};


struct vector4 {
    double tau, x, y, eta;
};


struct fluidCell {
    double mu_m;                // the effective mass of the cell [GeV^2]
    double tau, x, y, eta;      // spatial poision of the fluid cell
    vector3 beta;               // flow velocity of the fluid cell
    vector3 E_lab, B_lab;       // E and B fields in the lab frame
    vector4 drift_u_plus;       // drifting 4 velocity induced by EM fields
    vector4 drift_u_minus;
    vector4 drift_u_plus_2;     // drifting 4 velocity induced by EM fields
    vector4 drift_u_minus_2;
};


struct chargeSource {
    double x, y;
    double weight;
    double rapidity;
};


struct fluidCellSmall {
    int itau, ix, iy, ieta;
    float tau, x, y, eta;      // spatial poision of the fluid cell
    float temperature, ed, pressure, cs2;
    vector3 u;                 // flow velocity of the fluid cell
    vector3 E_lab, B_lab;      // E and B fields in the lab frame
};


#endif  // DATA_STRUCT_H_
