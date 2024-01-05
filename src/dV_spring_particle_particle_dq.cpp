#include <dV_spring_particle_particle_dq.h>

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d>     q1, double l0, double stiffness) {

    
    Eigen::Vector3d q = q1 - q0;
	double l = q.norm();
	Eigen::Vector3d q_unit = q.normalized();
	double dV = stiffness * (l - l0);
	f.segment<3>(0) = dV * q_unit;
	f.segment<3>(3) = -dV * q_unit;
}