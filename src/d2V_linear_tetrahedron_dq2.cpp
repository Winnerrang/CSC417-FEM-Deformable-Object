#include <d2V_linear_tetrahedron_dq2.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <d2psi_neo_hookean_dq2.h>
#include <quadrature_single_point.h>

void d2V_linear_tetrahedron_dq2(Eigen::Matrix1212d &H, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

   auto neohookean_linear_tet = [&](Eigen::Matrix1212d &dV, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
        
      //Code to compute non-integrated hessian matrix goes here
      // dphi/dX
      Eigen::Matrix43d dphidX;

      dphi_linear_tetrahedron_dX(dphidX, V, element, Eigen::Vector3d::Zero());

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

      Eigen::Matrix99d d2psidF2;

      d2psi_neo_hookean_dF2(d2psidF2, F, C, D);

      dV = Bj.transpose() * d2psidF2 * Bj;
    };

    //integrate the non-integrated hessian across the tetrahedral element
    quadrature_single_point(H, q, element, volume, neohookean_linear_tet);  
    

    //DO NOT REMOVE THIS CODE This code ensures that the hessian matrix is symmetric postive definite by projecting all
    //negative eigenvalues to small, postive values.
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix1212d> es(H);
    
    Eigen::MatrixXd DiagEval = es.eigenvalues().real().asDiagonal();
    Eigen::MatrixXd Evec = es.eigenvectors().real();
    
    for (int i = 0; i < 12; ++i) {
        if (es.eigenvalues()[i]<1e-6) {
            DiagEval(i,i) = 1e-3;
        }
    }

    H = Evec * DiagEval * Evec.transpose();

}
