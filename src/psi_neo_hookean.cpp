#include <psi_neo_hookean.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <iostream>
#include <cmath>
void psi_neo_hookean(double &psi, 
                     Eigen::Ref<const Eigen::Matrix3d> F,
                     double C, double D, bool debug) {

    double J = F.determinant();

    if (J == 0) {
		psi = nan("1");
		return;
		//std::cout << "Determinant is zero!" << std::endl;
		//exit(1);
	}
    psi = C * (pow(J, -2.0/3.0) * (F.transpose() * F).trace() - 3) + D * (J - 1) * (J - 1);

    if (debug) {
		std::cout << "F: " << std::endl;
		std::cout << F << std::endl;
		std::cout << "J: " << J << std::endl;
		std::cout << "pow(J, -2.0/3.0)" << pow(J, -2.0/3.0) << std::endl;
		std::cout << "(F.transpose() * F).trace()" << (F.transpose() * F).trace() << std::endl;
		std::cout << "psi: " << psi << std::endl;
	}
}