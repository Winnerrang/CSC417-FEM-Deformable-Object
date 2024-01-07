#include <Eigen/Dense>
#include <EigenTypes.h>

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
   
	double currentEnergy = f(x0);

	int step = 0;

	double tol = 1e-6;
	double alpha_tolerance = 1e-10;
	Eigen::VectorXd d;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

	
	while (step < maxSteps) {

		double oldEnergy = currentEnergy;

		g(tmp_g, x0);
		H(tmp_H, x0);

		solver.compute(tmp_H);
		d = -solver.solve(tmp_g);

		double alpha = 10000;

		// line search
		while (f(x0 + alpha * d) > currentEnergy) {
			
			// can't find a good step size, so just return
			if (alpha < alpha_tolerance) {
				return currentEnergy;
			}
			alpha /= 2;

		}

		x0 = x0 + alpha * d;
		currentEnergy = f(x0);

		// reaches minimum early, just return
		if (oldEnergy - currentEnergy < tol) {
			return currentEnergy;
		}
		step++;
	}
	return currentEnergy;
}
