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

		std::cout << "iteration: " << i << std::endl;
		g(tmp_g, x0);
		H(tmp_H, x0);
		
		solver.compute(tmp_H);
		d = -solver.solve(tmp_g);

		double alpha = 1;
		do {
			newEnergy = currentEnergy + alpha * d.transpose() * tmp_g;
			std::cout << "delta energy: " << alpha * d.transpose() * tmp_g << std::endl;
			std::cout << "Estimated Energy: " << newEnergy << std::endl;
			if (alpha < alpha_tolerance) {
				return currentEnergy;
			}
			
			if (newEnergy > currentEnergy) {
				alpha *= scaling;
			}
		} while (newEnergy > currentEnergy);

		
		x0 = x0 + alpha * d;
		newEnergy = f(x0);
		std::cout << "Real Energy: " << newEnergy << std::endl;

		if (currentEnergy - newEnergy < tol) {
			return newEnergy;
		}
		else {
			currentEnergy = newEnergy;
		}
	}
	
	//while (step < maxSteps) {

	//	double oldEnergy = currentEnergy;
	//	
	//	//auto start = high_resolution_clock::now();
	//	g(tmp_g, x0);
	//	//auto stop = high_resolution_clock::now();
	//	//auto duration = duration_cast<microseconds>(stop - start);

	//	// To get the value of duration use the count()
	//	// member function on the duration object
	//	//std::cout << "Caclculating g" << duration.count() << std::endl;

	//	start = high_resolution_clock::now();
	//	H(tmp_H, x0);
	//	stop = high_resolution_clock::now();
	//	duration = duration_cast<microseconds>(stop - start);
	//	//std::cout << "Caclculating H" << duration.count() << std::endl;

	//	start = high_resolution_clock::now();
	//	solver.compute(tmp_H);
	//	stop = high_resolution_clock::now();
	//	duration = duration_cast<microseconds>(stop - start);
	//	//std::cout << "Caclculating solver" << duration.count() << std::endl;
	//	
	//	assert(solver.info() == Eigen::Success);
	//	
	//	if (solver.info() != Eigen::Success) {
	//		std::cout << tmp_H << std::endl;
	//		std::cout << "Solver failed" << std::endl;
	//	}

	//	start = high_resolution_clock::now();
	//	d = -solver.solve(tmp_g);
	//	stop = high_resolution_clock::now();
	//	duration = duration_cast<microseconds>(stop - start);
	//	//std::cout << "Caclculating d" << duration.count() << std::endl;



	//	double alpha = 1;

	//	// line search
	//	start = high_resolution_clock::now();
	//	double newEnergy = f(x0 + alpha * d);
	//	stop = high_resolution_clock::now();
	//	duration = duration_cast<microseconds>(stop - start);
	//	//std::cout << "Caclculating newEnergy" << duration.count() << std::endl;

	//	start = high_resolution_clock::now();
	//	//std::cout << "Energy: " << newEnergy << std::endl;
	//	while (newEnergy > currentEnergy) {
	//		
	//		// can't find a good step size, so just return
	//		if (alpha < alpha_tolerance) {
	//			//exit(0);
	//			return currentEnergy;
	//			
	//		}
	//		alpha *= scaling;
	//		newEnergy = f(x0 + alpha * d);

	//		//std::cout << "Energy: " << newEnergy << std::endl;
	//	}

	//	stop = high_resolution_clock::now();
	//	duration = duration_cast<microseconds>(stop - start);
	//	//std::cout << "Caclculating alpha" << duration.count() << std::endl;
	//	x0 = x0 + alpha * d;
	//	currentEnergy = newEnergy;

	//	// reaches minimum early, just return
	//	if (oldEnergy - currentEnergy < tol) {
	//		//exit(0);
	//		return currentEnergy;
	//	}

	//	//std::cout << "One Iteration" << std::endl;
	//	step++;
	//}

	////exit(0);
	//return currentEnergy;


}
