#include <T_linear_tetrahedron.h>
#include <mass_matrix_linear_tetrahedron.h>
#include <iostream>
void T_linear_tetrahedron(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {


    Eigen::Matrix1212d M;

	mass_matrix_linear_tetrahedron(M, qdot, element, density, volume);

    
    Eigen::Vector12d qdotj;

    for (int vertexIdx = 0; vertexIdx < 4; vertexIdx++){

        qdotj.segment<3>(3 * vertexIdx) = qdot.segment<3>(3 * element(vertexIdx));
    }
	T = 0.5 * qdotj.transpose() * M * qdotj;

}