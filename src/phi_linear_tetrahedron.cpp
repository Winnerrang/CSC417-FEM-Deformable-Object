#include <phi_linear_tetrahedron.h>

void phi_linear_tetrahedron(Eigen::Vector4d &phi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> x) {

    Eigen::Vector3d deltaX = x - V.row(element(0)).transpose();
    Eigen::Vector3d delta1 = V.row(element(1)).transpose() - V.row(element(0)).transpose();
    Eigen::Vector3d delta2 = V.row(element(2)).transpose() - V.row(element(0)).transpose();
    Eigen::Vector3d delta3 = V.row(element(3)).transpose() - V.row(element(0)).transpose();

    Eigen::Matrix3d T;
    T.col(0) = delta1;
    T.col(1) = delta2;
    T.col(2) = delta3;

    Eigen::Vector3d temp_phi = T.inverse() * deltaX;
    phi(0) = 1 - temp_phi.sum();
    phi(1) = temp_phi(0);
    phi(2) = temp_phi(1);
    phi(3) = temp_phi(2);
}