//
// A pose representing a point and direction in 2D space
//
#pragma once
#include <string>
#include <vector>
#include <Eigen/Dense>

class Kalman {
public:
	Kalman();

	void init_filter();
	void init_state(double x, double y, double z);
	void filter(double x, double y, double z);
	void propagate();

	Eigen::Vector3d position();
	Eigen::Vector2d velocity();

	Eigen::Vector3d pred_position();
	Eigen::Vector2d pred_velocity();

	double dT();

private:

	int xdim; /*!< size of object state vector (x, y, z, xdot, ydot) */
	int ydim; /*!< size of measurement vector (x, y, z) */
	double dt;
	std::vector<double> q;
	std::vector<double> r;

	Eigen::MatrixXd A; /*!< Motion model */
	Eigen::MatrixXd C; /*!< Observation model */

	Eigen::MatrixXd Q; /*!< Process Noise */
	Eigen::MatrixXd R; /*!< Measurement Noise */

	Eigen::MatrixXd P0; /*!< Initial Covariance */

	Eigen::MatrixXd x_hat;    /*!< [x, y, z, xdot, ydot] 3D position and 2D velocity */
	Eigen::MatrixXd P_hat;    /*!< Covariance of the state */
	Eigen::MatrixXd x_check;
	Eigen::MatrixXd P_check;
	//Eigen::MatrixXd y_prev;


};
