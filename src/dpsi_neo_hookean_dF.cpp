#include <dpsi_neo_hookean_dF.h>
#include <cmath>
#include <iostream>
void dpsi_neo_hookean_dF(Eigen::Vector9d &dw, Eigen::Ref<const Eigen::Matrix3d> F, double C, double D) {

    double F11 = F(0, 0);
    double F12 = F(0, 1);
    double F13 = F(0, 2);
    double F21 = F(1, 0);
    double F22 = F(1, 1);
    double F23 = F(1, 2);
    double F31 = F(2, 0);
    double F32 = F(2, 1);
    double F33 = F(2, 2);


    double J = F.determinant();
    double F_square_Norm = (F.transpose() * F).trace();
    double repeatedPattern = F11 * F23 * F32 - F11 * F22 * F33 + F12 * F21 * F33 - F12 * F23 * F31 - F13 * F21 * F32 + F13 * F22 * F31;

    dw(0) = C * ((2 * F11) / pow(J, 2.0/3.0) - (2 * (F22 * F33 - F23 * F32) * (F_square_Norm)) / (3 * pow(J, 5.0/3.0))) - 2 * D * (F22 * F33 - F23 * F32) * (repeatedPattern + 1);
    dw(1) = C * ((2 * F12) / pow(J, 2.0/3.0) + (2 * (F21 * F33 - F23 * F31) * (F_square_Norm)) / (3 * pow(J, 5.0/3.0))) + 2 * D * (F21 * F33 - F23 * F31) * (repeatedPattern + 1);
    dw(2) = C * ((2 * F13) / pow(J, 2.0/3.0) - (2 * (F21 * F32 - F22 * F31) * (F_square_Norm)) / (3 * pow(J, 5.0/3.0))) - 2 * D * (F21 * F32 - F22 * F31) * (repeatedPattern + 1);
    dw(3) = C * ((2 * F21) / pow(J, 2.0/3.0) + (2 * (F12 * F33 - F13 * F32) * (F_square_Norm)) / (3 * pow(J, 5.0/3.0))) + 2 * D * (F12 * F33 - F13 * F32) * (repeatedPattern + 1);
    dw(4) = C * ((2 * F22) / pow(J, 2.0/3.0) - (2 * (F11 * F33 - F13 * F31) * (F_square_Norm)) / (3 * pow(J, 5.0/3.0))) - 2 * D * (F11 * F33 - F13 * F31) * (repeatedPattern + 1);
    dw(5) = C * ((2 * F23) / pow(J, 2.0/3.0) + (2 * (F11 * F32 - F12 * F31) * (F_square_Norm)) / (3 * pow(J, 5.0/3.0))) + 2 * D * (F11 * F32 - F12 * F31) * (repeatedPattern + 1);
    dw(6) = C * ((2 * F31) / pow(J, 2.0/3.0) - (2 * (F12 * F23 - F13 * F22) * (F_square_Norm)) / (3 * pow(J, 5.0/3.0))) - 2 * D * (F12 * F23 - F13 * F22) * (repeatedPattern + 1);
    dw(7) = C * ((2 * F32) / pow(J, 2.0/3.0) + (2 * (F11 * F23 - F13 * F21) * (F_square_Norm)) / (3 * pow(J, 5.0/3.0))) + 2 * D * (F11 * F23 - F13 * F21) * (repeatedPattern + 1);
    dw(8) = C * ((2 * F33) / pow(J, 2.0/3.0) - (2 * (F11 * F22 - F12 * F21) * (F_square_Norm)) / (3 * pow(J, 5.0/3.0))) - 2 * D * (F11 * F22 - F12 * F21) * (repeatedPattern + 1);

}