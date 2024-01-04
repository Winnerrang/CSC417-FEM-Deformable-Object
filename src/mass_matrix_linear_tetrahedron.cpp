 
 #include <mass_matrix_linear_tetrahedron.h>

 void mass_matrix_linear_tetrahedron(Eigen::Matrix1212d &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {
     M.setZero();
     for (int rowBlock = 0; rowBlock < 4; rowBlock++) {
         for (int colBlock = 0; colBlock < 4; colBlock++) {

             // got from solving integral using matlab
             double integralResult = (rowBlock == colBlock) ? 1 / 60.0 : 1 / 120.0;

             M.block<3, 3>(rowBlock * 3, colBlock * 3) = 6 * density * volume * integralResult * Eigen::Matrix3d::Identity();
         }
     }
   
 }