#include <Eigen/Dense>
#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>
#include <EigenTypes.h>
#include <iostream>
#include <chrono>
using namespace std::chrono;
//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  dt - the time step in seconds
//  mass - the mass matrix
//  force(f, q, qdot) - a function that computes the force acting on the FEM system. This takes q and qdot as parameters, returns the force in f.
//  stiffness(K, q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters, returns the stiffness matrix in K.  
//  tmp_force - scratch space to collect forces
//  tmp_stiffness - scratch space to collect stiffness matrix
//Output:
//  q - set q to the updated generalized coordinate using linearly implicit time integration
//  qdot - set qdot to the updated generalized velocity using linearly implicit time integration
template<typename FORCE, typename STIFFNESS> 
inline void linearly_implicit_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, 
                            const Eigen::SparseMatrixd &mass,  FORCE &force, STIFFNESS &stiffness, 
                            Eigen::VectorXd &tmp_force, Eigen::SparseMatrixd &tmp_stiffness) {
    auto start = high_resolution_clock::now();
    //std::cout << "Implicit Euler" << std::endl;
    force(tmp_force, q, qdot);

    //std::cout << "force" << std::endl;
    stiffness(tmp_stiffness, q, qdot);

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;


    // solve the system (M-dt^2 K) qdot^(n+1) = M * qdot^n + dt f(q^n)

    solver.compute(mass - dt * dt * tmp_stiffness);
    //std::cout << "comppute" << std::endl;
    assert(solver.info() == Eigen::Success);

    qdot = solver.solve(mass * qdot + dt * tmp_force);
    //std::cout << "solve" << std::endl;
    q = q + dt * qdot;
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    // To get the value of duration use the count()
    // member function on the duration object
    std::cout << duration.count() << std::endl;
    
    std::cout << "FPS: " << 1000000.0 / ((double) duration.count()) << std::endl;
}