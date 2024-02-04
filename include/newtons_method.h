#include <Eigen/Dense>
#include <EigenTypes.h>
#include <iostream>
#include <chrono>
#include <thread>
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
	double time_scale = 1.0;
	double right = 1.0;
	double left = 0;
	double precesion = 1e-4;

	// use binary search to find the time scale that makes the energy not nan
	// this simple method is similar to CCD method
	if (isnan(originalEnergy)) {
		while (right - left > precesion) {
			time_scale = (left + right) / 2;

			originalEnergy = f(x0 * time_scale);
			if (isnan(originalEnergy)) {
				right = time_scale;
			}
			else {
				left = time_scale;
			}
		}
	}

	// we might end up with the nan energy, so we just use the left value
	// where its energy is not nan and its precesion is good enough
	if (isnan(originalEnergy)) {
		time_scale = left;
		originalEnergy = f(x0 * time_scale);
	}

	if (isnan(originalEnergy)) {
		std::cout << "You did something wrong " << std::endl;
		exit(0);
	}
	/*do {
		originalEnergy = f(x0 * time_scale);
		if (isnan(originalEnergy)) {
			if (time_scale )
		}
			
	} while (isnan(originalEnergy));*/

	x0 = x0 * time_scale;

	double currentEnergy = originalEnergy;
	double newEnergy = currentEnergy;

	double tol = 1e-8;
	double alpha_tolerance = 1e-8;
	double scaling = 0.5;
	Eigen::VectorXd d;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	double alpha;

	for (int i = 0; i < maxSteps; i++){

		tmp_g.setZero();
		g(tmp_g, x0);

		
		tmp_H.setZero();
		H(tmp_H, x0);

		solver.compute(tmp_H);

		if (solver.info() != Eigen::Success)
		{
			std::cout <<"SHIT" << std::endl;
			exit(0);
		}

		//std::cout << "solve\n";
		d = -solver.solve(tmp_g);

		if (d.dot(tmp_g) > 0) {
			std::cout << "No!!!\n";
			exit(0);
		}
		



		 if (d.lpNorm<Eigen::Infinity>() < 1e-3){
		 	//std::cout << "Converge!!! " << currentEnergy << std::endl; 
		 	return currentEnergy;
		 }

		//line search
		alpha = 1;

		while (isnan(newEnergy) || newEnergy >= currentEnergy + 10e-8 * d.dot(tmp_g)){
			 if (alpha < alpha_tolerance){
			 	//std::cout << "No Choice...." << currentEnergy << std::endl; ;
			 	return currentEnergy;
			 } 
			
			newEnergy = f(x0 + alpha * d);

			if (isnan(newEnergy) || newEnergy >= currentEnergy + 10e-8 * d.dot(tmp_g)) alpha *= scaling;

		}

		if (isnan(newEnergy)) {
			//std::cout << x0 + alpha * d << std::endl;
			exit(0);
		}

		x0 = x0 + alpha * d;

		currentEnergy = newEnergy;
	}

	//std::cout << "I am tired.... Energy: " << currentEnergy << std::endl;

	return currentEnergy;
}
