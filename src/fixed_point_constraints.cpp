#include <fixed_point_constraints.h>
#include <algorithm>
void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {

	for (int i = 0; i < indices.size(); i++) {
		unsigned int vertexIndex = indices[i];

		for (int dimension = 0; dimension < 3; dimension++) {
			P.coeffRef(3 * vertexIndex + dimension, 3 * vertexIndex + dimension) = 0;
		}
	}

	// remove rows with all zeros
	Eigen::SparseMatrixd newP(q_size - 3 * indices.size(), q_size);

	int newRow = 0;
	for (int row = 0; row < P.rows(); row++) {
		int col = row;
		if (P.coeff(row, col) == 0) continue;

		newP.coeffRef(newRow, col) = P.coeff(row, col);

		newRow++;
	}

	P = newP;

}