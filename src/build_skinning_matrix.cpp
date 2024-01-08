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

    // assume the order of vertex in V is the order of the generalized coordinate q
    N.resize(V_skin.rows() * 3, V.rows() * 3);
    N.setZero();


    for (int vertexIdx = 0; vertexIdx < V_skin.rows(); vertexIdx++) {
        Eigen::Vector4d phi;

        int tetra_idx = -1;
        for (int tetIdx = 0; tetIdx < T.rows(); tetIdx++) {
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
        
        for (int dim = 0; dim < 3; dim++) {
			N.coeffRef(3 * vertexIdx + dim, 3 * T(tetra_idx, 0) + dim) = phi(0);
			N.coeffRef(3 * vertexIdx + dim, 3 * T(tetra_idx, 1) + dim) = phi(1);
			N.coeffRef(3 * vertexIdx + dim, 3 * T(tetra_idx, 2) + dim) = phi(2);
			N.coeffRef(3 * vertexIdx + dim, 3 * T(tetra_idx, 3) + dim) = phi(3);
		}
    }
}