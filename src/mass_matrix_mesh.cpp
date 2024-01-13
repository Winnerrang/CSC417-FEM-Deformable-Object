#include <mass_matrix_mesh.h>
#include <mass_matrix_linear_tetrahedron.h>

void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXi> T, double density, Eigen::Ref<const Eigen::VectorXd> v0) {

    M.resize(qdot.size(), qdot.size());
    M.setZero();

    // Best way to allocate an sparse matrix according https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
    typedef Eigen::Triplet<double> Tr;
    std::vector<Tr> tripleList;

    for (int tetraIdx = 0; tetraIdx < T.rows(); tetraIdx++) {
        Eigen::Matrix1212d M_tetra;
        Eigen::Vector12d qdot_tetra;
        Eigen::RowVector4i tetra_indices = T.row(tetraIdx);


        // find qdot for this tetrahedron
        qdot_tetra << qdot.segment<3>(T(tetraIdx, 0) * 3), qdot.segment<3>(T(tetraIdx, 1) * 3), qdot.segment<3>(T(tetraIdx, 2) * 3), qdot.segment<3>(T(tetraIdx, 3) * 3);

        mass_matrix_linear_tetrahedron(M_tetra, qdot_tetra, tetra_indices, density, v0(tetraIdx));

        for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				// if row x and col y in 12x12 matrix, To find the new row,
				// We first look at the xth row in the selection matrix,
				// and then look at which column in this row is nonzero.
				// for column, we look at yth row in the selection matrix,
				// and then do the same thing.

				int row = T(tetraIdx, i);
				int col = T(tetraIdx, j);

				tripleList.push_back(Tr(row * 3 + 0, col * 3 + 0, M_tetra(i * 3 + 0, j * 3 + 0)));
				tripleList.push_back(Tr(row * 3 + 0, col * 3 + 1, M_tetra(i * 3 + 0, j * 3 + 1)));
				tripleList.push_back(Tr(row * 3 + 0, col * 3 + 2, M_tetra(i * 3 + 0, j * 3 + 2)));

				

				tripleList.push_back(Tr(row * 3 + 1, col * 3 + 0, M_tetra(i * 3 + 1, j * 3 + 0)));
				tripleList.push_back(Tr(row * 3 + 1, col * 3 + 1, M_tetra(i * 3 + 1, j * 3 + 1)));
				tripleList.push_back(Tr(row * 3 + 1, col * 3 + 2, M_tetra(i * 3 + 1, j * 3 + 2)));
				

				tripleList.push_back(Tr(row * 3 + 2, col * 3 + 0, M_tetra(i * 3 + 2, j * 3 + 0)));
				tripleList.push_back(Tr(row * 3 + 2, col * 3 + 1, M_tetra(i * 3 + 2, j * 3 + 1)));
				tripleList.push_back(Tr(row * 3 + 2, col * 3 + 2, M_tetra(i * 3 + 2, j * 3 + 2)));


			}
		}

        // for (int row = 0; row < M_tetra.rows(); row++) {
        //     for (int col = 0; col < M_tetra.cols(); col++) {

        //         int targetRow = tetra_indices(row / 3) * 3 + row % 3;
        //         int targetCol = tetra_indices(col / 3) * 3 + col % 3;
        //         tripleList.push_back(Tr(targetRow, targetCol, M_tetra(row, col)));
		// 	}
		// }
    }
    M.setFromTriplets(tripleList.begin(), tripleList.end());

    

}