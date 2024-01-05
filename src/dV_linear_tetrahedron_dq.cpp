#include <dV_linear_tetrahedron_dq.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <dpsi_neo_hookean_dF.h>
#include <quadrature_single_point.h>
#include <iostream>

void dV_linear_tetrahedron_dq(Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {
   
    Eigen::Matrix3d T = Eigen::Matrix3d::Zero();
    T.col(0) = (V.row(element(1)) - V.row(element(0))).transpose();
    T.col(1) = (V.row(element(2)) - V.row(element(0))).transpose();
    T.col(2) = (V.row(element(3)) - V.row(element(0))).transpose();

    Eigen::Matrix3d T_inv = T.inverse();


    // dphi/dX
    Eigen::MatrixXd dphidX;
    dphidX.resize(4, 3);


    // -1^T * T_inv
    dphidX.row(0) = -T_inv.row(0) - T_inv.row(1) - T_inv.row(2);
    dphidX.block<3, 3>(1, 0) = T_inv;

    // [x0 x1 x2 x3]
    Eigen::MatrixXd q_temp;
    q_temp.resize(3, 4);
    q_temp.col(0) = q.segment<3>(3 * element(0));
    q_temp.col(1) = q.segment<3>(3 * element(1));
    q_temp.col(2) = q.segment<3>(3 * element(2));
    q_temp.col(3) = q.segment<3>(3 * element(3));

    // Deformation Gradient
    Eigen::Matrix3d F = q_temp * dphidX;


    // F = Bj * qj
    Eigen::MatrixXd Bj;
    Bj.resize(9, 12);
    Bj = Eigen::MatrixXd::Zero(9, 12);

    for (int i = 0; i < 3; i++) {
        for (int col = 0; col < dphidX.cols(); col++) {
            for (int row = 0; row < dphidX.rows(); row++) {
                Bj(3 * i + col, 3 * row) = dphidX(row, col);
            }
        }
    }
    

   auto neohookean_linear_tet = [&](Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
        	   
	   Eigen::Vector9d dPdF;
	   dpsi_neo_hookean_dF(dPdF, F, C, D);
	   dV = Bj.transpose() * dPdF;
    };

    quadrature_single_point(dV, q, element, volume, neohookean_linear_tet);  
    
}