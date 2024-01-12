#include <V_linear_tetrahedron.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <psi_neo_hookean.h>
#include <quadrature_single_point.h>

void V_linear_tetrahedron(double &energy, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

    
    
    // dphi/dX
    Eigen::Matrix43d dphidX;
    dphi_linear_tetrahedron_dX(dphidX, V, element, Eigen::Vector3d::Zero());


    // [x0 x1 x2 x3]
    Eigen::Matrix34d q_temp;
    q_temp.col(0) = q.segment<3>(3*element(0));
    q_temp.col(1) = q.segment<3>(3*element(1));
    q_temp.col(2) = q.segment<3>(3*element(2));
    q_temp.col(3) = q.segment<3>(3*element(3));
    
    Eigen::Matrix3d vol;
    vol.col(0) = q_temp.col(1) - q_temp.col(0);
    vol.col(1) = q_temp.col(2) - q_temp.col(0);
    vol.col(2) = q_temp.col(3) - q_temp.col(0);

    // we are integrating over the volume of the tetrahedron in real space not reference space
    double realVol = vol.determinant() / 6.0;

    // Deformation Gradient
    Eigen::Matrix3d F = q_temp * dphidX;

    auto neohookean_linear_tet = [&](double &e, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
        
        psi_neo_hookean(e, F, C, D);
    };

    quadrature_single_point(energy, q, element, std::abs(realVol), neohookean_linear_tet);
    
}