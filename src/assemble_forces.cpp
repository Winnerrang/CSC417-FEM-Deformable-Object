#include <assemble_forces.h>
#include <iostream>
#include <dV_linear_tetrahedron_dq.h>
#include <chrono>
void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0,
                     double C, double D) { 
    
    //std::cout << "Init" << std::endl;
    f.resize(q.size());
    f.setZero();
    //std::cout << "Finsh Init" << std::endl;
    for (int tetraIdx = 0; tetraIdx < T.rows(); tetraIdx++) {
        Eigen::Vector12d tetra_f;

        //std::cout << "enter\n";
        dV_linear_tetrahedron_dq(tetra_f, q, V, T.row(tetraIdx), v0(tetraIdx), C, D);
        //std::cout << "finish\n";
        //std::cout << "tetra_f:" << tetra_f << std::endl;

        for (int tetra_vertex_idx = 0; tetra_vertex_idx < 4; tetra_vertex_idx++) {
            
            /*if (3 * T(tetraIdx, tetra_vertex_idx) >= f.size()) {
                std::cout << "3 * T(tetraIdx, tetra_vertex_idx):" << 3 * T(tetraIdx, tetra_vertex_idx) << std::endl;
                std::cout << "f.size():" << f.size() << std::endl;
            }*/


            // f = -dV/dq
            f.segment<3>(3 * T(tetraIdx, tetra_vertex_idx)) -= tetra_f.segment<3>(3 * tetra_vertex_idx);
        }

        //std::cout << "filling\n";
    }

    //std::cout << "finish everything\n";
};