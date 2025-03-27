#include <iostream>
#include <string>
#include <vector>
#include "types/kalman_filter.hpp"

Kalman::Kalman() {
	xdim = 5; /*!< size of object state vector (x, y, z, xdot, ydot) */
	ydim = 3; /*!< size of measurement vector (x, y, z) */

	dt = 1.0;
	q = std::vector<double>{0.03,0.03,0.0001,0.05,0.05};
	r = std::vector<double>{0.03,0.03,0.0001};

	init_filter();
}


void Kalman::init_filter() {
	x_hat = Eigen::MatrixXd::Zero(xdim, 1);    /*!< [x, y, z, xdot, ydot] 3D position and 2D velocity */
	x_check = Eigen::MatrixXd::Zero(xdim, 1);

	A = Eigen::MatrixXd::Identity(xdim, xdim); /*!< Motion model */
	A(0, 3) = dt;
	A(1, 4) = dt;

	C = Eigen::MatrixXd::Zero(ydim, xdim); /*!< Observation model */
	C(0, 0) = 1;
	C(1, 1) = 1;
	C(2, 2) = 1;

	Q = Eigen::MatrixXd::Identity(xdim, xdim); /*!< Process Noise */
	for (int i=0; i<xdim; i++) {
		Q(i, i) = q[i];
	}

	R = Eigen::MatrixXd::Identity(ydim, ydim); /*!< Measurement Noise */
	for (int j=0; j<ydim; j++) {
		R(j, j) = r[j];
	}

	P0 = 10.0 * Q; /*!< Initial Covariance */
	P_hat = P0;  /*!< Covariance of the state */
	P_check = P_hat;

	//y_prev = Eigen::MatrixXd::Zero(ydim, 1);
}

void Kalman::init_state(double x, double y, double z) {
	x_hat = Eigen::MatrixXd::Zero(xdim, 1);
	x_hat(0) = x;
	x_hat(1) = y;
	x_hat(2) = z;

	x_check = Eigen::MatrixXd::Zero(xdim, 1);
	x_check(0) = x;
	x_check(1) = y;
	x_check(2) = z;
}

void Kalman::propagate() {
	x_check = A * x_hat;
	P_check = A * P_hat * A.transpose() + Q;
}

void Kalman::filter(double x, double y, double z) {
	propagate();

	Eigen::Vector3d y_obs = {x, y, z};
    Eigen::MatrixXd temp = C * P_check * C.transpose() + R;
    Eigen::MatrixXd K = (P_check * C.transpose()) * temp.inverse();
    Eigen::Vector3d error = y_obs - C * x_check;
    //std::cout << "error " << error.norm() << std::endl;
    x_hat = x_check + K * error;
    P_hat = (Eigen::MatrixXd::Identity(xdim, xdim) - K * C) * P_check;
}

Eigen::Vector3d Kalman::position() {
	Eigen::Vector3d pos;
	pos(0) = x_hat(0,0);
	pos(1) = x_hat(1,0);
	pos(2) = x_hat(2,0);
	return pos;
}

Eigen::Vector3d Kalman::pred_position() {
	Eigen::Vector3d pos;
	pos(0) = x_check(0,0);
	pos(1) = x_check(1,0);
	pos(2) = x_check(2,0);
	return pos;
}

Eigen::Vector2d Kalman::velocity() {
	Eigen::Vector2d vel;
	vel(0) = x_hat(3,0);
	vel(1) = x_hat(4,0);
	return vel;
}

Eigen::Vector2d Kalman::pred_velocity() {
	Eigen::Vector2d vel;
	vel(0) = x_check(3,0);
	vel(1) = x_check(4,0);
	return vel;
}

double Kalman::dT() {
	return dt;
}
