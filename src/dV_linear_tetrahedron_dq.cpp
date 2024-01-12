#include <dV_linear_tetrahedron_dq.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <dpsi_neo_hookean_dF.h>
#include <quadrature_single_point.h>
#include <iostream>
#include <chrono>
void dV_linear_tetrahedron_dq(Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {
   

    // dphi/dX
    Eigen::Matrix43d dphidX;

    dphi_linear_tetrahedron_dX(dphidX, V, element, Eigen::Vector3d::Zero());

    // [x0 x1 x2 x3]
    // Can't use MatrixXd because it will create calculation error in release mode
    Eigen::Matrix34d q_temp;
    q_temp.col(0) = q.segment<3>(3 * element(0));
    q_temp.col(1) = q.segment<3>(3 * element(1));
    q_temp.col(2) = q.segment<3>(3 * element(2));
    q_temp.col(3) = q.segment<3>(3 * element(3));

    Eigen::Matrix3d vol;
    vol.col(0) = q_temp.col(1) - q_temp.col(0);
    vol.col(1) = q_temp.col(2) - q_temp.col(0);
    vol.col(2) = q_temp.col(3) - q_temp.col(0);
    
    // we are integrating over the volume of the tetrahedron in real space not reference space
    double realVol = vol.determinant() / 6.0;


    // Deformation Gradient
    Eigen::Matrix3d F = q_temp * dphidX;


    // F = Bj * qj
    Eigen::Matrix<double, 9, 12> Bj;
    Bj.setZero();

    for (int i = 0; i < 3; i++) {
        for (int col = 0; col < dphidX.cols(); col++) {
            for (int row = 0; row < dphidX.rows(); row++) {
                Bj(3 * i + col, 3 * row + i) = dphidX(row, col);
            }
        }
    }
    


   auto neohookean_linear_tet = [&](Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
        	   
	   Eigen::Vector9d dPdF;
	   dpsi_neo_hookean_dF(dPdF, F, C, D);
	   dV = Bj.transpose() * dPdF;
    };

    //quadrature_single_point(dV, q, element, std::abs(realVol), neohookean_linear_tet);
   quadrature_single_point(dV, q, element, volume, neohookean_linear_tet);
}