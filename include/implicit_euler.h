#include <Eigen/Dense>
#include <EigenTypes.h>
#include <newtons_method.h>

#include <chrono>
#include <iostream>

using namespace std::chrono;
typedef Eigen::Triplet<double> Tr;
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
    
    //std::cout << "hello" << std::endl;
    auto start = high_resolution_clock::now();
    // initial conditions
    tmp_qdot = qdot;
    //std::cout << "1" << std::endl;
    tmp_force.resize(q.size());
    //std::cout << "2" << std::endl;
    tmp_stiffness.resize(q.size(), q.size());
    //std::cout << "world" << std::endl;
    auto grad = [&](Eigen::VectorXd &dVdq_dot, Eigen::VectorXd &new_q_dot) {
        //std::cout << "calculating gradient" << std::endl;
        force(tmp_force, q + dt * new_q_dot, new_q_dot);

        //std::cout << "Finish Calculating force" << std::endl;
        dVdq_dot = mass * (new_q_dot - qdot) - dt * tmp_force;
        
        //std::cout << "gradient calculated" << std::endl;

        };

    auto hessian = [&](Eigen::SparseMatrixd &d2Vdq_dotdq_dot, Eigen::VectorXd &new_q_dot) {
        //std::cout << "calculating hessian" << std::endl;
        stiffness(tmp_stiffness, q + dt * new_q_dot, new_q_dot);

        d2Vdq_dotdq_dot = mass - dt * dt * tmp_stiffness;
        /*std::cout << "hessian calculated" << std::endl;
         std::vector<Tr> tripletList;
         d2Vdq_dotdq_dot.setZero();
         for (int i = 0; i < q.size(); i++){
             tripletList.push_back(Tr(i, i, 1.0));
         }

        d2Vdq_dotdq_dot.setFromTriplets(tripletList.begin(), tripletList.end());*/
        };

    //std::cout << "3" << std::endl;
    newtons_method(tmp_qdot, energy, grad, hessian, 10, tmp_force, tmp_stiffness);
    //std::cout << "4" << std::endl;
    qdot = tmp_qdot;
    q = q + dt * qdot;
    //std::cout << "5" << std::endl;

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    // To get the value of duration use the count()
    // member function on the duration object
    std::cout << "FPS: " << 1000000.0 / ((double)duration.count()) << std::endl;
}
