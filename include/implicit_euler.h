#include <Eigen/Dense>
#include <EigenTypes.h>
#include <newtons_method.h>

//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  dt - the time step in seconds
//  mass - the mass matrix
//  energy(q, qdot) -  a function that computes the energy of the FEM system. This takes q and qdot as parameters, returns the energy value.
//  force(f, q, qdot) - a function that computes the force acting on the FEM system. This takes q and qdot as parameters, returns the force in f.
//  stiffness(K, q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters, returns the stiffness matrix in K.  
//  tmp_qdot - scratch space for storing velocities 
//  tmp_force - scratch space to collect forces
//  tmp_stiffness - scratch space to collect stiffness matrix
//Output:
//  q - set q to the updated generalized coordinate using linearly implicit time integration
//  qdot - set qdot to the updated generalized velocity using linearly implicit time integration
template<typename ENERGY, typename FORCE, typename STIFFNESS> 
inline void implicit_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, 
                            const Eigen::SparseMatrixd &mass,  ENERGY &energy, FORCE &force, STIFFNESS &stiffness, 
                            Eigen::VectorXd &tmp_qdot, Eigen::VectorXd &tmp_force, Eigen::SparseMatrixd &tmp_stiffness) {
    
    // initial conditions
    tmp_qdot = qdot;

    auto grad = [&](Eigen::VectorXd dVdq_dot, Eigen::VectorXd new_q_dot) {
        Eigen::VectorXd temp_force;
        force(temp_force, q + dt * new_q_dot, new_q_dot);
        dVdq_dot = mass * (new_q_dot - qdot) - dt * temp_force;
        };

    auto hessian = [&](Eigen::SparseMatrixd d2Vdq_dotdq_dot, Eigen::VectorXd new_q_dot) {
        stiffness(d2Vdq_dotdq_dot, q + dt * new_q_dot, new_q_dot);

        d2Vdq_dotdq_dot = mass - dt * dt * d2Vdq_dotdq_dot;
        };


    newtons_method(tmp_qdot, energy, grad, hessian, 10, tmp_force, tmp_stiffness);

    qdot = tmp_qdot;
    q = q + dt * qdot;

}
