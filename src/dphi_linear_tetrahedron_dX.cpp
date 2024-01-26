#include <dphi_linear_tetrahedron_dX.h>
#include <phi_linear_tetrahedron.h>
#include <iostream>
void dphi_linear_tetrahedron_dX(Eigen::Matrix43d &dphi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
    Eigen::Matrix3d T;
    T.col(0) = (V.row(element(1)) - V.row(element(0))).transpose();
    T.col(1) = (V.row(element(2)) - V.row(element(0))).transpose();
    T.col(2) = (V.row(element(3)) - V.row(element(0))).transpose();

    Eigen::Matrix3d T_inv = T.inverse();

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (isnan(T_inv(i, j))) {
                std::cout << "T_inv is nan" << std::endl;
                exit(1);
            }
        }
    }

    // -1^T * T_inv
    dphi.row(0) = -T_inv.row(0) - T_inv.row(1) - T_inv.row(2);
    dphi.block<3, 3>(1, 0) = T_inv;

}