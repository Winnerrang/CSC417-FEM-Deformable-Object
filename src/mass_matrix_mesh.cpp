#include <mass_matrix_mesh.h>
#include <mass_matrix_linear_tetrahedron.h>

void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXi> T, double density, Eigen::Ref<const Eigen::VectorXd> v0) {

    M.resize(qdot.size(), qdot.size());
    M.setZero();
    for (int tetraIdx = 0; tetraIdx < T.rows(); tetraIdx++) {
        Eigen::Matrix1212d M_tetra;
        Eigen::Vector12d qdot_tetra;
        Eigen::RowVector4i tetra_indices = T.row(tetraIdx);


        // find qdot for this tetrahedron
        qdot_tetra << qdot.segment<3>(T(tetraIdx, 0) * 3), qdot.segment<3>(T(tetraIdx, 1) * 3), qdot.segment<3>(T(tetraIdx, 2) * 3), qdot.segment<3>(T(tetraIdx, 3) * 3);

        mass_matrix_linear_tetrahedron(M_tetra, qdot_tetra, tetra_indices, density, v0(tetraIdx));

        for (int row = 0; row < M_tetra.rows(); row++) {
            for (int col = 0; col < M_tetra.cols(); col++) {

                int targetRow = tetra_indices(row / 3) * 3 + row % 3;
                int targetCol = tetra_indices(col / 3) * 3 + col % 3;
				M.coeffRef(targetRow, targetCol) += M_tetra(row, col);
			}
		}
    }

    

}