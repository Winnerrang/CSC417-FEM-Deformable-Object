#include <assemble_forces.h>
#include <iostream>
#include <dV_linear_tetrahedron_dq.h>
void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0,
                     double C, double D) { 
        
    f.resize(q.size());
    f = Eigen::VectorXd::Zero(q.size());
    for (int tetraIdx = 0; tetraIdx < T.rows(); tetraIdx++) {
        Eigen::Vector12d tetra_f;

        dV_linear_tetrahedron_dq(tetra_f, q, V, T.row(tetraIdx), v0(tetraIdx), C, D);
        //std::cout << "tetra_f:" << tetra_f << std::endl;

        for (int tetra_vertex_idx = 0; tetra_vertex_idx < 4; tetra_vertex_idx++) {

            // f = -dV/dq
            f.segment<3>(3 * T(tetraIdx, tetra_vertex_idx)) -= f.segment<3>(3 * tetra_vertex_idx);
        }
    }
    

};