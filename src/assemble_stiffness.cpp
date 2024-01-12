#include <assemble_stiffness.h>
#include <d2V_linear_tetrahedron_dq2.h>
#include <iostream>
#include <vector>
#include <thread>
#include <mutex>          // std::mutex
std::mutex mtx;           // mutex for critical section
typedef Eigen::Triplet<double> Tr;

void split_work(int startIdx, int endIdx, std::vector<Tr>* tripletList, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot,
	Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0,
	double C, double D) {
	for (int tetraIdx = startIdx; tetraIdx < endIdx; tetraIdx++) {
		Eigen::Matrix1212d d2Vdq2;

		d2V_linear_tetrahedron_dq2(d2Vdq2, q, V, T.row(tetraIdx), v0(tetraIdx), C, D);

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				// if row x and col y in 12x12 matrix, To find the new row,
				// We first look at the xth row in the selection matrix,
				// and then look at which column in this row is nonzero.
				// for column, we look at yth row in the selection matrix,
				// and then do the same thing.

				int row = T(tetraIdx, i);
				int col = T(tetraIdx, j);
				mtx.lock();
				tripletList->push_back(Tr(row * 3 + 0, col * 3 + 0, -d2Vdq2(i * 3 + 0, j * 3 + 0)));
				tripletList->push_back(Tr(row * 3 + 0, col * 3 + 1, -d2Vdq2(i * 3 + 0, j * 3 + 1)));
				tripletList->push_back(Tr(row * 3 + 0, col * 3 + 2, -d2Vdq2(i * 3 + 0, j * 3 + 2)));



				tripletList->push_back(Tr(row * 3 + 1, col * 3 + 0, -d2Vdq2(i * 3 + 1, j * 3 + 0)));
				tripletList->push_back(Tr(row * 3 + 1, col * 3 + 1, -d2Vdq2(i * 3 + 1, j * 3 + 1)));
				tripletList->push_back(Tr(row * 3 + 1, col * 3 + 2, -d2Vdq2(i * 3 + 1, j * 3 + 2)));


				tripletList->push_back(Tr(row * 3 + 2, col * 3 + 0, -d2Vdq2(i * 3 + 2, j * 3 + 0)));
				tripletList->push_back(Tr(row * 3 + 2, col * 3 + 1, -d2Vdq2(i * 3 + 2, j * 3 + 1)));
				tripletList->push_back(Tr(row * 3 + 2, col * 3 + 2, -d2Vdq2(i * 3 + 2, j * 3 + 2)));
				mtx.unlock();

			}
		}
	}
;
}
void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0, 
                     double C, double D) { 
        
    K.resize(q.size(), q.size());
	K.setZero();


	// Best way to allocate an sparse matrix according https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
	int numThread = 16;
	int delta = T.rows() / numThread;
	std::vector<Tr> tripletList;

	std::vector<std::thread> threads;
	for (int i = 0; i < numThread - 1; i++) {
		threads.push_back(std::thread(split_work, delta * i, delta * (i + 1), &tripletList, q, qdot, V, T, v0, C, D));
	}
	threads.push_back(std::thread(split_work, delta * (numThread - 1), T.rows(), &tripletList, q, qdot, V, T, v0, C, D));

	for (auto& th : threads) th.join();
	std::thread thread_object(split_work, 0, delta, &tripletList, q, qdot, V, T, v0, C, D);
	std::thread thread_object2(split_work, delta, delta * 2, &tripletList, q, qdot, V, T, v0, C, D);
	std::thread thread_object3(split_work, delta * 2, delta * 3, &tripletList, q, qdot, V, T, v0, C, D);
	std::thread thread_object4(split_work, delta * 3, T.rows(), &tripletList, q, qdot, V, T, v0, C, D);
	thread_object.join();
	thread_object2.join();
	thread_object3.join();
	thread_object4.join();

 //   for (int tetraIdx = 0; tetraIdx < T.rows(); tetraIdx++) {
 //       Eigen::Matrix1212d d2Vdq2;

 //       d2V_linear_tetrahedron_dq2(d2Vdq2, q, V, T.row(tetraIdx), v0(tetraIdx), C, D);

	//	for (int i = 0; i < 4; i++) {
	//		for (int j = 0; j < 4; j++) {
	//			// if row x and col y in 12x12 matrix, To find the new row,
	//			// We first look at the xth row in the selection matrix,
	//			// and then look at which column in this row is nonzero.
	//			// for column, we look at yth row in the selection matrix,
	//			// and then do the same thing.

	//			int row = T(tetraIdx, i);
	//			int col = T(tetraIdx, j);

	//			tripletList.push_back(Tr(row * 3 + 0, col * 3 + 0, -d2Vdq2(i * 3 + 0, j * 3 + 0)));
	//			tripletList.push_back(Tr(row * 3 + 0, col * 3 + 1, -d2Vdq2(i * 3 + 0, j * 3 + 1)));
	//			tripletList.push_back(Tr(row * 3 + 0, col * 3 + 2, -d2Vdq2(i * 3 + 0, j * 3 + 2)));

	//			

	//			tripletList.push_back(Tr(row * 3 + 1, col * 3 + 0, -d2Vdq2(i * 3 + 1, j * 3 + 0)));
	//			tripletList.push_back(Tr(row * 3 + 1, col * 3 + 1, -d2Vdq2(i * 3 + 1, j * 3 + 1)));
	//			tripletList.push_back(Tr(row * 3 + 1, col * 3 + 2, -d2Vdq2(i * 3 + 1, j * 3 + 2)));
	//			

	//			tripletList.push_back(Tr(row * 3 + 2, col * 3 + 0, -d2Vdq2(i * 3 + 2, j * 3 + 0)));
	//			tripletList.push_back(Tr(row * 3 + 2, col * 3 + 1, -d2Vdq2(i * 3 + 2, j * 3 + 1)));
	//			tripletList.push_back(Tr(row * 3 + 2, col * 3 + 2, -d2Vdq2(i * 3 + 2, j * 3 + 2)));


	//		}
	//	}
	//}    
	K.setFromTriplets(tripletList.begin(), tripletList.end());

};