#include <Eigen/Dense>
#include <EigenTypes.h>
#include <iostream>
#include <chrono>
using namespace std::chrono;
//Input:
//  x0 - initial point for newtons search
//  f(x) - function that evaluates and returns the cost function at x
//  g(dx, x) - function that evaluates and returns the gradient of the cost function in dx
//  H(dH, x) - function that evaluates and returns the Hessian in dH (as a sparse matrix).
//  max steps - the maximum newton iterations to take
//  tmp_g and tmp_H are scratch space to store gradients and hessians
//Output: 
//  x0 - update x0 to new value
template<typename Objective, typename Jacobian, typename Hessian>
double newtons_method(Eigen::VectorXd &x0, Objective &f, Jacobian &g, Hessian &H, unsigned int maxSteps, Eigen::VectorXd &tmp_g, Eigen::SparseMatrixd &tmp_H) {
   
	

	double originalEnergy = f(x0);
	double currentEnergy = originalEnergy;
	double newEnergy;
	int step = 0;

	double tol = 1e-8;
	double alpha_tolerance = 1e-8;
	double scaling = 0.5;
	Eigen::VectorXd d;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

	for (int i = 0; i < maxSteps; i++) {

		g(tmp_g, x0);
		H(tmp_H, x0);


		solver.compute(tmp_H);
		d = -solver.solve(tmp_g);
		if (d.dot(tmp_g) < 0) {
			std::cout << "Yessss" << std::endl;
			exit(0);
		}
		
		double alpha = 16;
		do {
			double newEnergy = f(x0 + alpha * d);

			if (alpha < alpha_tolerance) {
				return currentEnergy;
			}
			
			if (newEnergy > currentEnergy) {
				alpha *= scaling;
			}
		} while (newEnergy > currentEnergy);
		x0 = x0 + alpha * d;

		//std::cout << "Real Energy: " << newEnergy << std::endl;

		if (currentEnergy - newEnergy < tol) {
			return newEnergy;
		}
		else {
			currentEnergy = newEnergy;
		}
	}
	return currentEnergy;

}
