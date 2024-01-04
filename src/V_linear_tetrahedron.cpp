#include <V_linear_tetrahedron.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <psi_neo_hookean.h>
#include <quadrature_single_point.h>

void V_linear_tetrahedron(double &energy, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

    Eigen::Matrix3d T = Eigen::Matrix3d::Zero();
    T.col(0) = (V.row(element(1)) - V.row(element(0))).transpose();
    T.col(1) = (V.row(element(2)) - V.row(element(0))).transpose();
    T.col(2) = (V.row(element(3)) - V.row(element(0))).transpose();

    Eigen::Matrix3d T_inv = T.inverse();
    
    
    // dphi/dX
    Eigen::MatrixXd dphidX;
    dphidX.resize(4,3);

    
    // -1^T * T_inv
    dphidX.row(0) = -T_inv.row(0) - T_inv.row(1) - T_inv.row(2);
    dphidX.block<3,3>(1,0) = T_inv;

    // [x0 x1 x2 x3]
    Eigen::MatrixXd q_temp;
    q_temp.resize(3,4);
    q_temp.col(0) = q.segment<3>(3*element(0));
    q_temp.col(1) = q.segment<3>(3*element(1));
    q_temp.col(2) = q.segment<3>(3*element(2));
    q_temp.col(3) = q.segment<3>(3*element(3));

    // Deformation Gradient
    Eigen::Matrix3d F = q_temp * dphidX;

    auto neohookean_linear_tet = [&](double &e, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
        
        psi_neo_hookean(e, F, C, D);
    };

    quadrature_single_point(energy, q, element, volume, neohookean_linear_tet);  
    
}