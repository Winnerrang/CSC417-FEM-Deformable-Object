#include <build_skinning_matrix.h>
#include <phi_linear_tetrahedron.h>
#include <vector>
#include <iostream>
#include <cassert>
bool sameSide(Eigen::Ref<const Eigen::RowVector3d> p1,
    Eigen::Ref<const Eigen::RowVector3d> p2,
    Eigen::Ref<const Eigen::RowVector3d> p3,
    Eigen::Ref<const Eigen::RowVector3d> p4,
    Eigen::Ref<const Eigen::RowVector3d> c) {
    Eigen::RowVector3d normal = (p2 - p1).cross(p3 - p1);

    return (p4 - p1).dot(normal) * (c - p1).dot(normal) >= 0;

}
void build_skinning_matrix(Eigen::SparseMatrixd &N, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, 
                                                   Eigen::Ref<const Eigen::MatrixXd> V_skin) {

    // according to Visualize::update_vertex_positions, the skinning matrix will be times by the vertex matrix n x 3, where
    // each row is a vertex, and each column is a dimension

    N.resize(V_skin.rows(), V.rows());
    N.setZero();
    
    // Best way to allocate an sparse matrix according https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
    typedef Eigen::Triplet<double> Tr;
    std::vector<Tr> tripletList;

    for (int vertexIdx = 0; vertexIdx < V_skin.rows(); vertexIdx++) {

        Eigen::Vector4d phi;

        int tetra_idx = -1;
        for (int tetIdx = 0; tetIdx < T.rows(); tetIdx++) {
            //std::cout << "tetIdx: " << tetIdx << std::endl;
			Eigen::RowVector4i tet = T.row(tetIdx);
            Eigen::RowVector3d p1 = V.row(tet(0));
			Eigen::RowVector3d p2 = V.row(tet(1));
            Eigen::RowVector3d p3 = V.row(tet(2));
			Eigen::RowVector3d p4 = V.row(tet(3));
			Eigen::RowVector3d c = V_skin.row(vertexIdx);

			if (sameSide(p1, p2, p3, p4, c) && sameSide(p2, p3, p4, p1, c) && sameSide(p3, p4, p1, p2, c) && sameSide(p4, p1, p2, p3, c)) {
				tetra_idx = tetIdx;
                break;
            }

		}

        assert(tetra_idx != -1);
        phi_linear_tetrahedron(phi, V, T.row(tetra_idx), V_skin.row(vertexIdx).transpose());
        
        tripletList.push_back(Tr(vertexIdx, T(tetra_idx, 0), phi(0)));
        tripletList.push_back(Tr(vertexIdx, T(tetra_idx, 1), phi(1)));
        tripletList.push_back(Tr(vertexIdx, T(tetra_idx, 2), phi(2)));
        tripletList.push_back(Tr(vertexIdx, T(tetra_idx, 3), phi(3)));


    }
    N.setFromTriplets(tripletList.begin(), tripletList.end());
}