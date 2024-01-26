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
   
// 	double alpha = 1.0;
// 	for (int i = 0; i < maxSteps; i++) {
// 		tmp_g.setZero();
// 		g(tmp_g, x0);
// 		if (tmp_g.norm() < 1e-8) return tmp_g.norm();

// 		tmp_H.setZero();
// 		H(tmp_H, x0);
// 		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
// 		solver.compute(tmp_H);

// 		Eigen::VectorXd dx = solver.solve(tmp_g);

// 		alpha = 1.0;
// 		while((f(x0 - alpha*dx) > f(x0) - 1e-8*alpha*tmp_g.dot(dx))&&(alpha>1e-12))
// 			   alpha *= 0.5;

// 		x0 -= alpha*dx;
//    }
//    return tmp_g.norm();

	double originalEnergy = f(x0);
	std::cout << "Original Energy: " << originalEnergy << std::endl;
	double currentEnergy = originalEnergy;
	double newEnergy = currentEnergy;
	int step = 0;

	double tol = 1e-8;
	double alpha_tolerance = 1e-8;
	double scaling = 0.5;
	Eigen::VectorXd d;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	double alpha;

	for (int i = 0; i < maxSteps; i++){
		//std::cout << "Current Energy: " << currentEnergy << std::endl;

		Eigen::VectorXd x_temp = x0;
		g(tmp_g, x_temp);

		assert(x_temp == x0);

		x_temp = x0;
		H(tmp_H, x0);
		assert(x_temp == x0);

		//tmp_H.makeCompressed();

		//std::cout << "compute\n";
		solver.compute(tmp_H);
		//std::cout << "Determinant: " << solver.determinant() << std::endl;
		//std::cout << "Heesian nonZero: " << tmp_H. << std::endl;
		//std::cout << "finish compute\n";

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
		//std::cout << "finish solve\n";
		

		//if (tmp_g.norm() < 1e-8){
		//	//std::cout << i << std::endl;
		//	//std::cout << "Converge!!!\n"; 
		//	return currentEnergy;
		//}

		//std::cout << "checking\n";
		 if (d.lpNorm<Eigen::Infinity>() < 1e-3){
		 	std::cout << "Converge!!! " << currentEnergy << std::endl; 
		 	return currentEnergy;
		 }

		//line search
		alpha = 1;

		//std::cout << "line search\n";
		while (newEnergy >= currentEnergy){
			 if (alpha < alpha_tolerance){
			 	std::cout << "No Choice...." << currentEnergy << std::endl; ;
			 	return currentEnergy;
			 } 
			
			newEnergy = f(x0 + alpha * d);

			if (newEnergy >= currentEnergy + 10e-8 * d.dot(tmp_g)) alpha *= scaling;

		}

		x0 = x0 + alpha * d;

		currentEnergy = newEnergy;
		std::cout << "finsh\n";
		std::cout << "Energy: " << currentEnergy << "Gradient Norm " << d.lpNorm<Eigen::Infinity>() << std::endl;
		if (d.lpNorm<Eigen::Infinity>() > 500) {
			//std::cout << tmp_g << std::endl;
			//std::cout << "Hessian\n" << d << std::endl;
			std::this_thread::sleep_for(std::chrono::seconds(100));
		}
	}

	// g(tmp_g, x0);
	
	// std::cout << "I am tired.... Original Energy: " << originalEnergy
	// 							<< " NewEnergy: " << newEnergy 
	// 							<< "gradient norm: " << tmp_g.norm()<< std::endl;
	/*if (tmp_g.norm() > 10e-5){
		std::cout << "Keep Going" << std::endl;
		continue;
	} */
	std::cout << "I am tired.... Energy: " << currentEnergy << std::endl;
	/*if (currentEnergy > 0.2) {
		std::cout << x0 << std::endl;
	}*/
	return currentEnergy;
}
