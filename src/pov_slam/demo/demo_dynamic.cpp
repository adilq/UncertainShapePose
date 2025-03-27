#include <cmath>
#include <iostream>
#include <numeric>
#include <random>

#include "g2o/core/block_solver.h" //ok
#include "g2o/core/optimization_algorithm_levenberg.h" //ok
#include "g2o/solvers/linear_solver_eigen.h" // g2o/solvers/linear_solver_eigen.h

#include "g2o/types/types_slam.h" //moved
#include "g2o/types/edge_se3.h" //moved
#include "g2o/types/vertex_se3.h" //moved
#include "g2o/types/edge_xyz_prior.h" // moved
#include "g2o/types/edge_se3_pointxyz.h" // moved

#include "g2o/types/edge_se3_prior.h" // moved

#include "g2o/types/edge_se3_mixture.h" // moved
#include "g2o/types/edge_se3_pointxyz_noninformative_mixture_vem.h"


#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/uniform.hpp>

#include <opencv2/opencv.hpp>

#include "types/Pose2.hpp"
#include "types/kalman_filter.hpp"


using namespace std;
using namespace g2o;
using namespace g2o::internal;

typedef BlockSolver<BlockSolverTraits<-1, -1> > SlamBlockSolver;
typedef LinearSolverEigen<SlamBlockSolver::PoseMatrixType> SlamLinearSolver;

std::random_device rd;
std::mt19937 generator(5303612481028);
// 73856093275028
// 270923102481028

// 4 is good
// 10 has ambiguity at the start
// Used 12 -- good
// 14- symmetry, all objects shifted right
// 13 a dynamic goes up 
// 666 dynamic goes up
// 0 moved goes up

// std::mt19937 generator(rd());

Isometry3 Isometry3ByAngle(Eigen::Vector3d trans, double angle) {
	Eigen::Vector3d rotAxisAngle(0, 0, 1);
	Eigen::AngleAxisd rotation(angle, rotAxisAngle.normalized());
	Eigen::Isometry3d result = (Eigen::Isometry3d) rotation.toRotationMatrix();
	result.translation() = trans;
	return result;
}

double randnormal(double stddev) {
	std::normal_distribution<double> dist(0, stddev);
	return dist(generator);
}

Eigen::Vector3d CorrputedVector3d(Eigen::Vector3d trans,
		double trans_std = 0.1) {
	std::normal_distribution<double> trans_dist(0, trans_std);
	double r1 = trans_dist(generator);
	double r2 = trans_dist(generator);
	double r3 = 0; //trans_dist(generator);
	Eigen::Vector3d corrupted_trans = trans + Eigen::Vector3d(r1, r2, r3);
	return corrupted_trans;
}

Isometry3 CorrputedIsometry3ByAngle(Eigen::Vector3d trans, double angle,
		double trans_std = 0.1, double rot_std = 0.1) {
	std::normal_distribution<double> rot_dist(0, rot_std);
	Eigen::Vector3d rotAxisAngle(0, 0, 1);
	double corrupted_angle = angle + rot_dist(generator);
	Eigen::AngleAxisd rotation(corrupted_angle, rotAxisAngle.normalized());
	Eigen::Isometry3d result = (Eigen::Isometry3d) rotation.toRotationMatrix();
	result.translation() = CorrputedVector3d(trans, trans_std);
	return result;
}

inline Isometry2 Isometry2ByAngle(Eigen::Vector2d trans, double angle) {
	Eigen::Rotation2Dd rot(angle);
	Eigen::Isometry2d result = (Eigen::Isometry2d) rot;
	result.translation() = trans;
	return result;
}

struct RectObs {
	Eigen::Vector3d center;
	Eigen::Vector3d fl;
	Eigen::Vector3d fr;
	Eigen::Vector3d bl;
	Eigen::Vector3d br;
	Eigen::Vector3d vel;
};

struct RectVert {
	Pose2 center;
	Pose2 fl;
	Pose2 fr;
	Pose2 bl;
	Pose2 br;
	Pose2 vel;
};

//TODO: tracker should be yaw aware
//TODO: vertex id should not be hard-coded..

class Object {
public:
	Object(double pos_true_x, double pos_true_y, double pos_true_z,
			double pos_est_x, double pos_est_y, double pos_est_z,
			int staticness = 1) {

		confidence = a / (a + b);

		not_moving = staticness;
		// if the object is dynamic, generate a velocity
		// assume dynamic objs have good init estimate

		if (not_moving == 0) {
			std::uniform_real_distribution<> dis_vel_lin(0.1, 0.3);
			std::uniform_real_distribution<> dis_vel_ang(-0.2, 0.2);

			vel_lin = dis_vel_lin(generator);
			vel_ang = dis_vel_ang(generator);
		}

		std::uniform_real_distribution<> dis_heading(0, 2.0*M_PI);
		heading = dis_heading(generator);

		pose_true = Pose2(pos_true_x, pos_true_y, pos_true_z, heading);
		pose_est = Pose2(pos_est_x, pos_est_y, pos_est_z, heading + randnormal(0.1));
		pose_track = pose_est;

		tracker.init_state(pos_est_x, pos_est_y, pos_est_z);

		// bbox
		fr_true = pose_true + fr_local;
		fl_true = pose_true + fl_local;
		br_true = pose_true + br_local;
		bl_true = pose_true + bl_local;

		fr_est = pose_est + fr_local;
		fl_est = pose_est + fl_local;
		br_est = pose_est + br_local;
		bl_est = pose_est + bl_local;

		fr_track = fr_est;
		fl_track = fl_est;
		br_track = br_est;
		bl_track = bl_est;
	}

	void set_true_pose(Pose2 p) {
		pose_true = p;
		fr_true = pose_true + fr_local;
		fl_true = pose_true + fl_local;
		br_true = pose_true + br_local;
		bl_true = pose_true + bl_local;
	}

	void set_est_pose(Pose2 p) {
		pose_est = p;
		fr_est = pose_est + fr_local;
		fl_est = pose_est + fl_local;
		br_est = pose_est + br_local;
		bl_est = pose_est + bl_local;
	}

	void step() {
		if (not_moving == 0) {
			Pose2 curr_pose_true = Pose2(pose_true.x(), pose_true.y(), pose_true.z(), heading);
			Pose2 motion(vel_lin,0,vel_ang);
			Pose2 next_pose_true = curr_pose_true + motion;

			heading = next_pose_true.yaw();
			pose_true = Pose2(next_pose_true.x(), next_pose_true.y(), next_pose_true.z(), heading);

			fr_true = pose_true + fr_local;
			fl_true = pose_true + fl_local;
			br_true = pose_true + br_local;
			bl_true = pose_true + bl_local;
		}
	}

	void filter(Pose2 obs) {
		if (not_moving == 0) {
			tracker.filter(obs.x(), obs.y(), obs.z());
			auto res = tracker.position();
			auto res_vel = tracker.velocity();

			double dx = res(0) - pose_track.x();
			double dy = res(1) - pose_track.y();
			double dz = res(2) - pose_track.z();

			fr_track.x() += dx; fr_track.y() += dy; fr_track.z() += dz;
			fl_track.x() += dx; fl_track.y() += dy; fl_track.z() += dz;
			bl_track.x() += dx; bl_track.y() += dy; bl_track.z() += dz;
			br_track.x() += dx; br_track.y() += dy; br_track.z() += dz;

			pose_track = Pose2(res(0),res(1),res(2),0);
			// std::cout <<  "  vel true " << vel_lin << std::endl;
			// std::cout <<  "  vel est " << res_vel.norm() << std::endl;
		}
	}

	void propagate() {
		if (not_moving == 0) {
			tracker.propagate();
			auto res = tracker.pred_position();
			auto res_vel = tracker.pred_velocity();

			double dx = res(0) - pose_track.x();
			double dy = res(1) - pose_track.y();
			double dz = res(2) - pose_track.z();

			fr_track.x() += dx; fr_track.y() += dy; fr_track.z() += dz;
			fl_track.x() += dx; fl_track.y() += dy; fl_track.z() += dz;
			bl_track.x() += dx; bl_track.y() += dy; bl_track.z() += dz;
			br_track.x() += dx; br_track.y() += dy; br_track.z() += dz;

			pose_track = Pose2(res(0),res(1),res(2),0);
		}
	}

	Pose2 getTrackerPosition() {
		Eigen::Vector3d pos = tracker.position();
		Eigen::Vector2d vel = tracker.velocity();
		double dir = atan2(vel(1),vel(0));
		return Pose2(pos(0),pos(1),pos(2),dir);
	}

	Pose2 getTrackerVelocity() {
		Eigen::Vector2d vel = tracker.velocity();
		double dir = atan2(vel(1),vel(0));
		return Pose2(vel(0),vel(1),0,dir);
	}

	double getConfidence() const {
		return confidence;
	}

	double getConfidenceBK() const {
		return confidence_bk;
	}

	RectObs getMeasurement(Pose2 robot) {
		RectObs obs;
		obs.center = (pose_true - robot).position_eigen();
		obs.fl = (fl_true - robot).position_eigen();
		obs.fr = (fr_true - robot).position_eigen();
		obs.bl = (bl_true - robot).position_eigen();
		obs.br = (br_true - robot).position_eigen();

		return obs;
	}

	RectObs getVertRel() {
		RectObs obs;
		obs.center = Vector3(0,0,0);
		obs.fl = fl_local.position_eigen();
		obs.fr = fr_local.position_eigen();
		obs.bl = bl_local.position_eigen();
		obs.br = br_local.position_eigen();

		return obs;
	}

	RectObs getExpectedMeasurement(Pose2 robot) {
		RectObs obs;
		obs.center = (pose_est - robot).position_eigen();
		obs.fl = (fl_est - robot).position_eigen();
		obs.fr = (fr_est - robot).position_eigen();
		obs.bl = (bl_est - robot).position_eigen();
		obs.br = (br_est - robot).position_eigen();

		return obs;
	}

	void backup_values() {
		a_bk = a;
		b_bk = b;
		mu_bk = mu;
		sig_bk = sig;
		confidence_bk = confidence;

		reset_values();
	}

	void restore_values() {
		a = a_bk;
		b = b_bk;
		mu = mu_bk;
		sig = sig_bk;
		confidence = confidence_bk;
	}

	void reset_values() {
		a = 2.0;
		b = 2.0;
		mu = 0.0;
		sig = 0.1;
		confidence = a / (a + b);
	}

	void updateProbability(double change, double std_change) {

		double s_weight = 1;

		if (type == 0 && inlier == false)
			s_weight = 3;  // dynamic outlier, drop fast
		if (type == 0 && inlier == true)
			s_weight = 0;  // dynamic inlier, rise slow
		if (type == 1 && inlier == false)
			s_weight = 0;  // static outlier, drop slow
		if (type == 1 && inlier == true)
			s_weight = 3;  // static inlier, rise fast

		//std::cout << "Amount of change: " << change << " std " << std_change <<std::endl;
		double tolerance = 20.0 * std_change;

		double s_sq = 1.0 / (1.0 / (pow(sig, 2)) + 1.0 / (pow(std_change, 2)));
		double m = s_sq * (mu / (pow(sig, 2)) + change / (pow(std_change, 2)));

		auto Kvals = K(s_weight);
		double K1 = Kvals.first;
		double K2 = Kvals.second;

		boost::math::normal_distribution<double> norm_dist(mu, sig);
		boost::math::uniform_distribution<double> uniform_dist(0.0, tolerance);

		double C1 = K1 * std::max(boost::math::pdf(norm_dist, change), eps);
		double C2 =
				fabs(change) >= tolerance - eps ?
						K2 * boost::math::pdf(uniform_dist, tolerance) :
						K2 * boost::math::pdf(uniform_dist, fabs(change));
		C1 = max(eps, C1);
		C2 = max(eps, C2);

		double C_norm = C1 + C2;
		C1 /= C_norm;
		C2 /= C_norm;

		inlier = C1 >= C2 ? true : false;

		double mu_prime = C1 * m + C2 * mu;
		sig = sqrt(
				C1 * (s_sq + pow(m, 2)) + C2 * (pow(sig, 2) + pow(mu, 2))
						- pow(mu_prime, 2));

		double gamma = (a + type * s_weight + 1) / (a + b + s_weight + 1);
		double eta = (a + type * s_weight) / (a + b + s_weight + 1);
		double theta = C1 * gamma + C2 * eta;
		double alpha = ((a + type * s_weight + 2) * (a + type * s_weight + 1)
				/ ((a + b + s_weight + 1) * (a + b + s_weight + 2)));
		double beta = ((a + type * s_weight + 1) * (a + type * s_weight)
				/ ((a + b + s_weight + 1) * (a + b + s_weight + 2)));

		mu = mu_prime;
		a = (C1 * theta * alpha + beta * C2 * theta - pow(theta, 2))
				/ (pow(theta, 2) - C1 * alpha - C2 * beta);
		b = ((C1 * theta * alpha + beta * C2 * theta - pow(theta, 2))
				* (1 - C1 * gamma - C2 * eta)
				/ ((pow(theta, 2) - C1 * alpha - C2 * beta)
						* (C1 * gamma + C2 * eta)));

		double cap = 25;
		if (a > cap || b > cap) {
			double ratio = max(a, b) / cap;
			a /= ratio;
			b /= ratio;
		}

		confidence = a / (a + b);
		life++;
	}

	std::vector<double> getBetaParams() {
		return std::vector<double> {a, b};
	}

	bool isDynamic() const {
		return not_moving == 0;
	}

	double r_in = 0.5;

	Pose2 pose_true;
	Pose2 pose_est;
	Pose2 pose_track;


	Pose2 fr_true, fl_true, br_true, bl_true;
	Pose2 fr_est, fl_est, br_est, bl_est;
	Pose2 fr_track, fl_track, br_track, bl_track;

	Pose2 fr_local = Pose2(0.25,-0.1,0);
	Pose2 fl_local = Pose2(0.25,0.1,0);
	Pose2 br_local = Pose2(-0.25,-0.1,0);
	Pose2 bl_local = Pose2(-0.25,0.1,0);

	Kalman tracker;

private:
	std::pair<double, double> K(double k) const {
		double lk1 = (lgamma(a + b) + lgamma(a + k * type + 1)
				+ lgamma(b + k - k * type))
				- (lgamma(a) + lgamma(b) + lgamma(a + b + k + 1));
		double lk2 = (lgamma(a + b) + lgamma(a + k * type)
				+ lgamma(b + k - k * type + 1))
				- (lgamma(a) + lgamma(b) + lgamma(a + b + k + 1));

		double k1 = exp(lk1);
		double k2 = exp(lk2);

		double ks = k1 + k2;

		k1 /= ks;
		k2 /= ks;

		k1 = max(eps, k1);
		k2 = max(eps, k2);

		return std::make_pair(k1, k2);
	}

	int type = 0;
	double a = 2.0;
	double b = 2.0;
	double mu = 0.0;
	double sig = 0.1;
	double eps = 1e-5;
	bool inlier = true;
	float confidence = 0.0;
	int life = 0;

	double a_bk = 2.0;
	double b_bk = 2.0;
	double mu_bk = 0.0;
	double sig_bk = 0.1;
	float confidence_bk = 0.0;

	int not_moving = 1;
	double vel_lin = 0;
	double vel_ang = 0;
	double heading = 0;
};


class Environment {
public:
	Environment(size_t max_iter, double meas_stddev, double max_meas_range, int window_size) :
			max_iter_(max_iter),
			meas_stddev_(meas_stddev),
			max_meas_range_(max_meas_range),
			window_size_(window_size),
			num_static_objs_(0),
			num_dynamic_objs_(0) {
		init_solver();
		static_objects_ = {};
		static_object_status_ = {};
		static_object_reconsted_ = {};

		dynamic_objects_ = {};
		dynamic_object_status_ = {};

		canvas_max_x_ = std::numeric_limits<double>::min();
		canvas_min_x_ = std::numeric_limits<double>::max();
		canvas_max_y_ = std::numeric_limits<double>::min();
		canvas_min_y_ = std::numeric_limits<double>::max();
	
		file_iter.open(path1, ios::in | ios::app); //VC
		// file_iter << "Step" << "," << "Iteration" << "," << "Object ID" << "," << "Stationarity" << "," << "True Stationarity" << std::endl; //VC
	}

	std::vector<Object> static_objects_;
	std::vector<Object> dynamic_objects_;

	void set_robot_pose(double pos_true_x, double pos_true_y, double pos_true_z, double pos_true_yaw,
			 double pos_est_x, double pos_est_y, double pos_est_z, double pos_est_yaw) {
		robot_pose_true_ = Pose2(pos_true_x, pos_true_y, pos_true_z, pos_true_yaw);
		robot_pose_est_ = Pose2(pos_est_x, pos_est_y, pos_est_z, pos_est_yaw);
		robot_pose_odom_ = robot_pose_est_;
		robot_pose_est_init_ = robot_pose_est_;

		history_pose_est_.push_back(robot_pose_est_);

		update_canvas_size(robot_pose_true_);
		update_canvas_size(robot_pose_est_);
		update_canvas_size(robot_pose_odom_);
	}

	void add_object(double pos_true_x, double pos_true_y, double pos_true_z,
			 double pos_est_x, double pos_est_y, double pos_est_z, int staticness = 1) {
		Object new_obj(pos_true_x, pos_true_y, pos_true_z,
				 pos_est_x, pos_est_y, pos_est_z, staticness);

		if (staticness == 1) {
			static_objects_.push_back(new_obj);
			static_object_status_.push_back(0);
			static_object_reconsted_.push_back(0);
			num_static_objs_++;
		} else {
			dynamic_objects_.push_back(new_obj);
			dynamic_object_status_.push_back(0);
			num_dynamic_objs_++;
		}

		update_canvas_size(new_obj.pose_true);
		update_canvas_size(new_obj.pose_est);
	}

	void reconst_static_object(size_t j) {
		auto corrupted_meas_center = history_meas_static_.at(step_).at(j).center;
		auto corrupted_meas_fl = history_meas_static_.at(step_).at(j).fl;
		auto corrupted_meas_fr = history_meas_static_.at(step_).at(j).fr;
		auto corrupted_meas_bl = history_meas_static_.at(step_).at(j).bl;
		auto corrupted_meas_br = history_meas_static_.at(step_).at(j).br;

		Pose2 meas_pose_center(corrupted_meas_center(0),corrupted_meas_center(1),corrupted_meas_center(2),0);
		Pose2 new_est_center = robot_pose_est_ + meas_pose_center;

		Pose2 meas_pose_fl(corrupted_meas_fl(0),corrupted_meas_fl(1),corrupted_meas_fl(2),0);
		Pose2 new_est_fl = robot_pose_est_ + meas_pose_fl;

		Pose2 meas_pose_fr(corrupted_meas_fr(0),corrupted_meas_fr(1),corrupted_meas_fr(2),0);
		Pose2 new_est_fr = robot_pose_est_ + meas_pose_fr;

		Pose2 meas_pose_bl(corrupted_meas_bl(0),corrupted_meas_bl(1),corrupted_meas_bl(2),0);
		Pose2 new_est_bl = robot_pose_est_ + meas_pose_bl;

		Pose2 meas_pose_br(corrupted_meas_br(0),corrupted_meas_br(1),corrupted_meas_br(2),0);
		Pose2 new_est_br = robot_pose_est_ + meas_pose_br;

		static_objects_.at(j).pose_est = new_est_center;
		static_objects_.at(j).fl_est = new_est_fl;
		static_objects_.at(j).fr_est = new_est_fr;
		static_objects_.at(j).bl_est = new_est_bl;
		static_objects_.at(j).br_est = new_est_br;

		static_objects_.at(j).reset_values();
		static_object_reconsted_.at(j) = 1;
		// remove past measurements after reconstruction
		for (int i(0); i<=step_; i++) {
			history_meas_static_.at(i).erase(j);
			auto pos = std::find(history_obs_static_.at(i).begin(), history_obs_static_.at(i).end(), j);
			if (pos != history_obs_static_.at(i).end()) history_obs_static_.at(i).erase(pos);
		}
	}

	void reconst_dynamic_object(size_t j) {
		auto corrupted_meas_center = history_meas_dynamic_.at(step_).at(j).center;
		auto corrupted_meas_fl = history_meas_dynamic_.at(step_).at(j).fl;
		auto corrupted_meas_fr = history_meas_dynamic_.at(step_).at(j).fr;
		auto corrupted_meas_bl = history_meas_dynamic_.at(step_).at(j).bl;
		auto corrupted_meas_br = history_meas_dynamic_.at(step_).at(j).br;

		Pose2 meas_pose_center(corrupted_meas_center(0),corrupted_meas_center(1),corrupted_meas_center(2),0);
		Pose2 new_est_center = robot_pose_est_ + meas_pose_center;

		Pose2 meas_pose_fl(corrupted_meas_fl(0),corrupted_meas_fl(1),corrupted_meas_fl(2),0);
		Pose2 new_est_fl = robot_pose_est_ + meas_pose_fl;

		Pose2 meas_pose_fr(corrupted_meas_fr(0),corrupted_meas_fr(1),corrupted_meas_fr(2),0);
		Pose2 new_est_fr = robot_pose_est_ + meas_pose_fr;

		Pose2 meas_pose_bl(corrupted_meas_bl(0),corrupted_meas_bl(1),corrupted_meas_bl(2),0);
		Pose2 new_est_bl = robot_pose_est_ + meas_pose_bl;

		Pose2 meas_pose_br(corrupted_meas_br(0),corrupted_meas_br(1),corrupted_meas_br(2),0);
		Pose2 new_est_br = robot_pose_est_ + meas_pose_br;


		dynamic_objects_.at(j).pose_est = new_est_center;
		dynamic_objects_.at(j).fl_est = new_est_fl;
		dynamic_objects_.at(j).fr_est = new_est_fr;
		dynamic_objects_.at(j).bl_est = new_est_bl;
		dynamic_objects_.at(j).br_est = new_est_br;

		dynamic_objects_.at(j).reset_values();
	}

	void filter_dynamic_object(size_t j) {
		auto corrupted_meas = history_meas_dynamic_.at(step_).at(j).center;
		Pose2 meas_pose(corrupted_meas(0),corrupted_meas(1),corrupted_meas(2),0);
		Pose2 new_est_pose = robot_pose_est_ + meas_pose;
		new_est_pose.setYaw(0);
		dynamic_objects_.at(j).filter(new_est_pose);
		history_meas_dynamic_.at(step_).at(j).vel = dynamic_objects_.at(j).getTrackerVelocity().position_eigen();
	}

	void step() {
		// Propagate robot
		double n_x = randnormal(motion_stddev_);
		double n_y = randnormal(motion_stddev_);
		double n_yaw = randnormal(motion_stddev_);

		Pose2 motion(0.2,0,0.1);
		Pose2 noisy_motion(0.2+n_x, n_y, 0.1+n_yaw);

		robot_pose_true_ = robot_pose_true_ + motion;
		robot_pose_est_ = robot_pose_est_ + noisy_motion;
		robot_pose_odom_ = robot_pose_odom_ + noisy_motion;

		history_odom_.push_back(noisy_motion);

		// update static object measurements
		std::vector<size_t> obs_static_objs = get_observed_static_objects();
		std::map<size_t, RectObs> meas_static;

		for (auto j : obs_static_objs) {
			static_objects_.at(j).backup_values();

			auto gt_meas = static_objects_.at(j).getMeasurement(robot_pose_true_);
			meas_static[j].center = CorrputedVector3d(gt_meas.center, meas_stddev_);
			meas_static[j].fl = CorrputedVector3d(gt_meas.fl, meas_stddev_);
			meas_static[j].fr = CorrputedVector3d(gt_meas.fr, meas_stddev_);
			meas_static[j].bl = CorrputedVector3d(gt_meas.bl, meas_stddev_);
			meas_static[j].br = CorrputedVector3d(gt_meas.br, meas_stddev_);

			update_canvas_size(static_objects_.at(j).pose_true);
			update_canvas_size(static_objects_.at(j).pose_est);
		}

		history_meas_static_.push_back(meas_static);
		history_obs_static_.push_back(obs_static_objs);

		// step dynamic objects and get measurements
		std::vector<size_t> obs_dynamic_objs = get_observed_dynamic_objects();
		std::map<size_t, RectObs> meas_dynamic;

		for (auto j : obs_dynamic_objs) {
			dynamic_objects_.at(j).backup_values();
			dynamic_objects_.at(j).step(); // move actual pose
			dynamic_objects_.at(j).propagate(); // update belief

			auto gt_meas = dynamic_objects_.at(j).getMeasurement(robot_pose_true_);
			meas_dynamic[j].center = CorrputedVector3d(gt_meas.center, meas_stddev_);
			meas_dynamic[j].fl = CorrputedVector3d(gt_meas.fl, meas_stddev_);
			meas_dynamic[j].fr = CorrputedVector3d(gt_meas.fr, meas_stddev_);
			meas_dynamic[j].bl = CorrputedVector3d(gt_meas.bl, meas_stddev_);
			meas_dynamic[j].br = CorrputedVector3d(gt_meas.br, meas_stddev_);
			meas_dynamic[j].vel = dynamic_objects_.at(j).getTrackerVelocity().position_eigen();

			update_canvas_size(dynamic_objects_.at(j).pose_true);
			update_canvas_size(dynamic_objects_.at(j).pose_est);
		}

		history_meas_dynamic_.push_back(meas_dynamic);
		history_obs_dynamic_.push_back(obs_dynamic_objs);

		update_canvas_size(robot_pose_true_);
		update_canvas_size(robot_pose_est_);
		update_canvas_size(robot_pose_odom_);
	}

	void optimize() {
		std::cout << "Step " << step_ << endl;

		pose_error_file.open(path3, ios::in | ios::app); //VC
		pose_error_file << pose_error() << "," << step_ << "," << 0 << std::endl; 

		for (size_t i = 0; i < max_iter_; i++) {
			optimize_one_step(i);
			pose_error_file << pose_error() << "," <<  step_ << "," <<  i+1 << std::endl;
			
			for (int j=0; j < static_objects_.size(); j++){
				file_iter << step_ << "," << i << "," << j << "," << static_objects_.at(j).getConfidence() << "," << static_objects_.at(j).r_in << "," << 1 << std::endl;
			}
			
			for (int j=0; j < dynamic_objects_.size(); j++){
				file_iter << step_ << "," << i << "," << j + static_objects_.size()<< "," << dynamic_objects_.at(j).getConfidence() << "," << dynamic_objects_.at(j).r_in << "," << 0 << std::endl;
			}
		}

		// update conf for static objects
		for (auto j : history_obs_static_.at(step_)) {
			static_objects_.at(j).restore_values();
			auto expected_meas = static_objects_.at(j).getExpectedMeasurement(robot_pose_est_);
			auto meas = history_meas_static_.at(step_).at(j);

			static_objects_.at(j).updateProbability(calc_error(meas, expected_meas), obj_est_stddev_);
		}

		// update conf for dynamic objects
		for (auto j : history_obs_dynamic_.at(step_)) {
			dynamic_objects_.at(j).restore_values();
			auto expected_meas = dynamic_objects_.at(j).getExpectedMeasurement(robot_pose_est_);
			auto meas = history_meas_dynamic_.at(step_).at(j);

			dynamic_objects_.at(j).updateProbability(calc_error(meas, expected_meas), obj_est_stddev_);
		}

		// check if stationary objects require reconstruction
		std::cout << "    Static" << endl;
		for (size_t j(0); j< num_static_objs_; ++j) {
			std::cout << "        Object " << j << " new score: " << static_objects_.at(j).getConfidence() << std::endl;
			if (static_objects_.at(j).getConfidence() < reconst_lb) {
				reconst_static_object(j);
			}
		}

		// filter dynamic objects (and reconstruct if off by a lot)
		std::cout << "    Dynamic" << endl;
		for (size_t j(0); j< num_dynamic_objs_; ++j) {
			std::cout << "        Object " << num_static_objs_+j << " new score: " << dynamic_objects_.at(j).getConfidence() << std::endl;
			if (dynamic_objects_.at(j).getConfidence() < reconst_lb) {
				reconst_dynamic_object(j);
			}

			filter_dynamic_object(j);
		}

		history_pose_est_.push_back(robot_pose_est_);

		std::cout << "    Number of inactive objects: " << num_inactive_objects() << endl;
		std::cout << "    Robot pose error: " << pose_error() << endl;

		for (size_t j(0); j<num_dynamic_objs_; j++) {
			std::cout << "      Dynamic object " << num_static_objs_+j <<
					" pose error: " << (dynamic_objects_.at(j).pose_true - dynamic_objects_.at(j).pose_est).norm() << endl;
		}

		step_++;
		pose_error_file.close();
	}

	std::vector<size_t> get_inactive_object_ids() {
		std::vector<size_t> list;
		for (size_t id(0); id < static_object_status_.size(); ++id) {
			if (static_object_status_.at(id) > 0) list.push_back(id);
		}
		return list;
	}

	int num_inactive_objects() {
		return std::accumulate(static_object_status_.begin(), static_object_status_.end(), 0);
	}

	double pose_error() {
		return (robot_pose_true_ - robot_pose_est_).norm();
	}

private:

	double calc_error(RectObs obs, RectObs exp_obs) {
		auto error_c = pow((obs.center - exp_obs.center).norm(),2);
		auto error_fl = pow((obs.fl - exp_obs.fl).norm(),2);
		auto error_fr = pow((obs.fr - exp_obs.fr).norm(),2);
		auto error_bl = pow((obs.bl - exp_obs.bl).norm(),2);
		auto error_br = pow((obs.br - exp_obs.br).norm(),2);
		return sqrt(0.2 * (error_c + error_fl + error_fr + error_bl + error_br));
	}

	void init_solver() {
		// Init solver
		SlamLinearSolver* linearSolver = new SlamLinearSolver;
		SlamBlockSolver* blockSolver = new SlamBlockSolver(linearSolver);
		linearSolver->setBlockOrdering(false);
		OptimizationAlgorithmLevenberg *solver =
				new OptimizationAlgorithmLevenberg(blockSolver);

		optimizer_.setAlgorithm(solver);

		ParameterSE3Offset *paramOffset = new ParameterSE3Offset();
		paramOffset->setId(999999);
		optimizer_.addParameter(paramOffset);
	}

	std::vector<size_t> get_observed_static_objects() {
		std::vector<size_t> observed_objs;
		for (size_t j(0); j<num_static_objs_; ++j) {
			observed_objs.push_back(j);
		}
		return observed_objs;
	}

	std::vector<size_t> get_observed_dynamic_objects() {
		std::vector<size_t> observed_objs;
		for (size_t j(0); j<num_dynamic_objs_; ++j) {
			observed_objs.push_back(j);
		}
		return observed_objs;
	}

	void optimize_one_step(size_t iter) {

		int window_size = min(window_size_, step_+2);
		int window_start = max(0, step_-window_size+1);

		int g_id = 0;

		if (debug_) std::cout << "Iteration " << iter << std::endl;
		// All pointers
		std::vector<VertexSE3*> robots = { };
		std::vector<VertexSE3*> dynamic_objs = { };
		// std::vector<EdgeSE3NonInformativeMixture*> lm_obj_mmniedges = { }; // only for current step
		std::vector<std::pair<int, EdgeSE3PointXYZNonInformativeMixtureELBO*> > lm_edges = { }; // only for current step

		// Update POCD scores based current robot pose and object position estimates
		if (debug_) std::cout << "    Step 1: Update stationarity scores" << std::endl;
		for (size_t j : history_obs_static_.at(step_)) {
			if (debug_) std::cout << "        Obj " << j << " prev score: " << static_objects_.at(j).getConfidence() << std::endl;

			auto expected_meas = static_objects_.at(j).getExpectedMeasurement(robot_pose_est_);
			auto meas = history_meas_static_.at(step_).at(j);

			static_objects_.at(j).updateProbability(calc_error(meas, expected_meas), obj_est_stddev_);

			if (debug_) std::cout << "              new score: " << static_objects_.at(j).getConfidence() << std::endl;
		}

		for (size_t j : history_obs_dynamic_.at(step_)) {
			if (debug_) std::cout << "        Obj " << num_static_objs_ + j << " prev score: " << dynamic_objects_.at(j).getConfidence() << std::endl;

			auto expected_meas = dynamic_objects_.at(j).getExpectedMeasurement(robot_pose_est_);
			auto meas = history_meas_dynamic_.at(step_).at(j);

			dynamic_objects_.at(j).updateProbability(calc_error(meas, expected_meas), obj_est_stddev_);

			if (debug_) std::cout << "              new score: " << dynamic_objects_.at(j).getConfidence() << std::endl;
		}

		// With the new POCD score, estimate robot pose
		if (debug_) std::cout << "    Step 2: Update robot pose" << std::endl;
		//std::cout << "Iteration " << iter << std::endl;
		//std::cout << "  ### robot ### " << iter << std::endl;
		for (int i(step_-window_size+2); i<=step_+1; ++i) {
			VertexSE3* robot = new VertexSE3;
			robot->setId(i);
			g_id = max(g_id, i);
			//std::cout << "robot id: " << i << std::endl;
			if (i == step_+1) {
				robot->setEstimate(robot_pose_est_.isometry3());
			} else {
				robot->setEstimate(history_pose_est_.at(i).isometry3());
			}
			robots.push_back(robot);
			optimizer_.addVertex(robot);
		}

		// Add the odoms
		for (int i(step_-window_size+2); i<=step_; ++i) {
			EdgeSE3Orig *odom_i = new EdgeSE3Orig;
			odom_i->vertices()[0] = optimizer_.vertex(i); //from
			odom_i->vertices()[1] = optimizer_.vertex(i + 1); //to
			odom_i->setMeasurement(history_odom_.at(i).isometry3());
			odom_i->setInformation(motion_stddev_ * g2o::EdgeSE3Orig::InformationType::Identity());
			odom_i->setParameterId(0, 999999);
		}

		// Add the landmarks
		//std::cout << "  ### SO ### " << iter << std::endl;
		std::vector<std::vector<int> > static_obj_vertex_ids(static_objects_.size(), std::vector<int>(4,-1));
		for (auto j : history_obs_static_.at(step_)) {
			auto vert_local = static_objects_.at(j).getVertRel();

			VertexSE3 *landmark_c = new VertexSE3;
			landmark_c->setId(step_+5*j+2);
			optimizer_.addVertex(landmark_c);

			//static_objs.push_back(landmark_c);

			VertexPointXYZ *landmark_fl = new VertexPointXYZ;
			landmark_fl->setId(step_+5*j+3);
			landmark_fl->setEstimate(static_objects_.at(j).fl_est.position_eigen());
			optimizer_.addVertex(landmark_fl);
			static_obj_vertex_ids.at(j).at(0) = step_+5*j+3;

			EdgeXYZPrior *landmark_fl_prior = new EdgeXYZPrior();
			landmark_fl_prior->setVertex(0, optimizer_.vertex(step_+5*j+3));
			landmark_fl_prior->setMeasurement(static_objects_.at(j).fl_est.position_eigen());
			landmark_fl_prior->setInformation(0.5 * g2o::EdgeXYZPrior::InformationType::Identity());
			landmark_fl_prior->setParameterId(0, 999999);
			optimizer_.addEdge(landmark_fl_prior);

			EdgeSE3PointXYZ *landmark_fl_local = new EdgeSE3PointXYZ;
			landmark_fl_local->setVertex(0, optimizer_.vertex(step_+5*j+2));
			landmark_fl_local->setVertex(1, optimizer_.vertex(step_+5*j+3));
			landmark_fl_local->setMeasurement(vert_local.fl);
			landmark_fl_local->setInformation(20 * g2o::EdgeXYZPrior::InformationType::Identity());
			landmark_fl_local->setParameterId(0, 999999);
			optimizer_.addEdge(landmark_fl_local);

			VertexPointXYZ *landmark_fr = new VertexPointXYZ;
			landmark_fr->setId(step_+5*j+4);
			landmark_fr->setEstimate(static_objects_.at(j).fr_est.position_eigen());
			optimizer_.addVertex(landmark_fr);
			static_obj_vertex_ids.at(j).at(1) = step_+5*j+4;

			EdgeXYZPrior *landmark_fr_prior = new EdgeXYZPrior();
			landmark_fr_prior->setVertex(0, optimizer_.vertex(step_+5*j+4));
			landmark_fr_prior->setMeasurement(static_objects_.at(j).fr_est.position_eigen());
			landmark_fr_prior->setInformation(0.5 * g2o::EdgeXYZPrior::InformationType::Identity());
			landmark_fr_prior->setParameterId(0, 999999);
			optimizer_.addEdge(landmark_fr_prior);

			EdgeSE3PointXYZ *landmark_fr_local = new EdgeSE3PointXYZ;
			landmark_fr_local->setVertex(0, optimizer_.vertex(step_+5*j+2));
			landmark_fr_local->setVertex(1, optimizer_.vertex(step_+5*j+4));
			landmark_fr_local->setMeasurement(vert_local.fr);
			landmark_fr_local->setInformation(20 * g2o::EdgeXYZPrior::InformationType::Identity());
			landmark_fr_local->setParameterId(0, 999999);
			optimizer_.addEdge(landmark_fr_local);

			VertexPointXYZ *landmark_bl = new VertexPointXYZ;
			landmark_bl->setId(step_+5*j+5);
			landmark_bl->setEstimate(static_objects_.at(j).bl_est.position_eigen());
			optimizer_.addVertex(landmark_bl);
			static_obj_vertex_ids.at(j).at(2) = step_+5*j+5;

			EdgeXYZPrior *landmark_bl_prior = new EdgeXYZPrior();
			landmark_bl_prior->setVertex(0, optimizer_.vertex(step_+5*j+5));
			landmark_bl_prior->setMeasurement(static_objects_.at(j).bl_est.position_eigen());
			landmark_bl_prior->setInformation(0.5 * g2o::EdgeXYZPrior::InformationType::Identity());
			landmark_bl_prior->setParameterId(0, 999999);
			optimizer_.addEdge(landmark_bl_prior);

			EdgeSE3PointXYZ *landmark_bl_local = new EdgeSE3PointXYZ;
			landmark_bl_local->setVertex(0, optimizer_.vertex(step_+5*j+2));
			landmark_bl_local->setVertex(1, optimizer_.vertex(step_+5*j+5));
			landmark_bl_local->setMeasurement(vert_local.bl);
			landmark_bl_local->setInformation(20 * g2o::EdgeXYZPrior::InformationType::Identity());
			landmark_bl_local->setParameterId(0, 999999);
			optimizer_.addEdge(landmark_bl_local);

			VertexPointXYZ *landmark_br = new VertexPointXYZ;
			landmark_br->setId(step_+5*j+6);
			g_id = max(g_id, step_+5*(int)j+6);
			landmark_br->setEstimate(static_objects_.at(j).br_est.position_eigen());
			optimizer_.addVertex(landmark_br);
			static_obj_vertex_ids.at(j).at(3) = step_+5*j+6;

			EdgeXYZPrior *landmark_br_prior = new EdgeXYZPrior();
			landmark_br_prior->setVertex(0, optimizer_.vertex(step_+5*j+6));
			landmark_br_prior->setMeasurement(static_objects_.at(j).br_est.position_eigen());
			landmark_br_prior->setInformation(0.5 * g2o::EdgeXYZPrior::InformationType::Identity());
			landmark_br_prior->setParameterId(0, 999999);
			optimizer_.addEdge(landmark_br_prior);

			EdgeSE3PointXYZ *landmark_br_local = new EdgeSE3PointXYZ;
			landmark_br_local->setVertex(0, optimizer_.vertex(step_+5*j+2));
			landmark_br_local->setVertex(1, optimizer_.vertex(step_+5*j+6));
			landmark_br_local->setMeasurement(vert_local.br);
			landmark_br_local->setInformation(20 * g2o::EdgeXYZPrior::InformationType::Identity());
			landmark_br_local->setParameterId(0, 999999);
			optimizer_.addEdge(landmark_br_local);
		}
		//std::cout << "  ### DO ### " << iter << std::endl;
		auto dyn_obj_vertex_ids = get_dyn_obj_temporal_obs(window_start, window_size);

		for (int i(window_start); i<=step_; ++i) {
			for(size_t j : history_obs_dynamic_.at(i)) {
				auto vert_local = dynamic_objects_.at(j).getVertRel();

				VertexSE3 *object_c = new VertexSE3;
				object_c->setId(dyn_obj_vertex_ids.at(i).at(j));
				optimizer_.addVertex(object_c);

				dynamic_objs.push_back(object_c);

				VertexPointXYZ *object_fl = new VertexPointXYZ;
				object_fl->setId(dyn_obj_vertex_ids.at(i).at(j)+1);
				object_fl->setEstimate(dynamic_objects_.at(j).fl_est.position_eigen());
				optimizer_.addVertex(object_fl);

				EdgeSE3PointXYZ *object_fl_local = new EdgeSE3PointXYZ;
				object_fl_local->setVertex(0, optimizer_.vertex(dyn_obj_vertex_ids.at(i).at(j)));
				object_fl_local->setVertex(1, optimizer_.vertex(dyn_obj_vertex_ids.at(i).at(j)+1));
				object_fl_local->setMeasurement(vert_local.fl);
				object_fl_local->setInformation(20 * g2o::EdgeXYZPrior::InformationType::Identity());
				object_fl_local->setParameterId(0, 999999);
				optimizer_.addEdge(object_fl_local);

				VertexPointXYZ *object_fr = new VertexPointXYZ;
				object_fr->setId(dyn_obj_vertex_ids.at(i).at(j)+2);
				object_fr->setEstimate(dynamic_objects_.at(j).fr_est.position_eigen());
				optimizer_.addVertex(object_fr);

				EdgeSE3PointXYZ *object_fr_local = new EdgeSE3PointXYZ;
				object_fr_local->setVertex(0, optimizer_.vertex(dyn_obj_vertex_ids.at(i).at(j)));
				object_fr_local->setVertex(1, optimizer_.vertex(dyn_obj_vertex_ids.at(i).at(j)+2));
				object_fr_local->setMeasurement(vert_local.fr);
				object_fr_local->setInformation(20 * g2o::EdgeXYZPrior::InformationType::Identity());
				object_fr_local->setParameterId(0, 999999);
				optimizer_.addEdge(object_fr_local);

				VertexPointXYZ *object_bl = new VertexPointXYZ;
				object_bl->setId(dyn_obj_vertex_ids.at(i).at(j)+3);
				object_bl->setEstimate(dynamic_objects_.at(j).bl_est.position_eigen());
				optimizer_.addVertex(object_bl);

				EdgeSE3PointXYZ *object_bl_local = new EdgeSE3PointXYZ;
				object_bl_local->setVertex(0, optimizer_.vertex(dyn_obj_vertex_ids.at(i).at(j)));
				object_bl_local->setVertex(1, optimizer_.vertex(dyn_obj_vertex_ids.at(i).at(j)+3));
				object_bl_local->setMeasurement(vert_local.bl);
				object_bl_local->setInformation(20 * g2o::EdgeXYZPrior::InformationType::Identity());
				object_bl_local->setParameterId(0, 999999);
				optimizer_.addEdge(object_bl_local);

				VertexPointXYZ *object_br = new VertexPointXYZ;
				object_br->setId(dyn_obj_vertex_ids.at(i).at(j)+4);
				object_br->setEstimate(dynamic_objects_.at(j).br_est.position_eigen());
				optimizer_.addVertex(object_br);

				EdgeSE3PointXYZ *object_br_local = new EdgeSE3PointXYZ;
				object_br_local->setVertex(0, optimizer_.vertex(dyn_obj_vertex_ids.at(i).at(j)));
				object_br_local->setVertex(1, optimizer_.vertex(dyn_obj_vertex_ids.at(i).at(j)+4));
				object_br_local->setMeasurement(vert_local.br);
				object_br_local->setInformation(20 * g2o::EdgeXYZPrior::InformationType::Identity());
				object_br_local->setParameterId(0, 999999);
				optimizer_.addEdge(object_br_local);

			}
		}

		// add vertex observations
		for (int i(window_start); i<=step_; ++i) {
			for (size_t j : history_obs_static_.at(i)) {
				std::vector<double> weights_lm { static_objects_.at(j).getConfidence(),
										1 - static_objects_.at(j).getConfidence() };
				std::vector<double> beta_params = static_objects_.at(j).getBetaParams();

				// fl
				EdgeSE3PointXYZNonInformativeMixtureELBO *landmarkObservation_fl =
						new EdgeSE3PointXYZNonInformativeMixtureELBO;

				EdgeSE3PointXYZ *landmarkObservation_fl_p1 = new EdgeSE3PointXYZ;
				landmarkObservation_fl_p1->vertices()[0] = optimizer_.vertex(i+1); //from
				landmarkObservation_fl_p1->vertices()[1] = optimizer_.vertex(step_+5*j+3); //to
				landmarkObservation_fl_p1->setMeasurement(history_meas_static_.at(i).at(j).fl);
				landmarkObservation_fl_p1->setInformation(
						(1.0/obj_est_stddev_) * g2o::EdgeSE3PointXYZ::InformationType::Identity()); // Gaussian
				landmarkObservation_fl_p1->setParameterId(0, 999999);

				EdgeSE3PointXYZNonInformative *landmarkObservation_fl_p2 =
						new EdgeSE3PointXYZNonInformative;
				landmarkObservation_fl_p2->vertices()[0] = optimizer_.vertex(i+1); //from
				landmarkObservation_fl_p2->vertices()[1] = optimizer_.vertex(step_+5*j+3); //to
				landmarkObservation_fl_p2->setMeasurement(history_meas_static_.at(i).at(j).fl);
				landmarkObservation_fl_p2->setInformation(
						10 * g2o::EdgeSE3PointXYZ::InformationType::Identity()); // Uniform
				landmarkObservation_fl_p2->setParameterId(0, 999999);

				auto components_fl_lm = std::make_pair(landmarkObservation_fl_p1,
						landmarkObservation_fl_p2);
				landmarkObservation_fl->setParameterId(0, 999999);
				landmarkObservation_fl->initializeComponents(components_fl_lm,
						weights_lm, beta_params);
				optimizer_.addEdge(landmarkObservation_fl);

				// fr
				EdgeSE3PointXYZNonInformativeMixtureELBO *landmarkObservation_fr =
						new EdgeSE3PointXYZNonInformativeMixtureELBO;

				EdgeSE3PointXYZ *landmarkObservation_fr_p1 = new EdgeSE3PointXYZ;
				landmarkObservation_fr_p1->vertices()[0] = optimizer_.vertex(
						i + 1); //from
				landmarkObservation_fr_p1->vertices()[1] = optimizer_.vertex(
						step_ + 5 * j + 4); //to
				landmarkObservation_fr_p1->setMeasurement(
						history_meas_static_.at(i).at(j).fr);
				landmarkObservation_fr_p1->setInformation(
						(1.0/obj_est_stddev_) * g2o::EdgeSE3PointXYZ::InformationType::Identity()); // Gaussian
				landmarkObservation_fr_p1->setParameterId(0, 999999);

				EdgeSE3PointXYZNonInformative *landmarkObservation_fr_p2 =
						new EdgeSE3PointXYZNonInformative;
				landmarkObservation_fr_p2->vertices()[0] = optimizer_.vertex(
						i + 1); //from
				landmarkObservation_fr_p2->vertices()[1] = optimizer_.vertex(
						step_ + 5 * j + 4); //to
				landmarkObservation_fr_p2->setMeasurement(
						history_meas_static_.at(i).at(j).fr);
				landmarkObservation_fr_p2->setInformation(
						10 * g2o::EdgeSE3PointXYZ::InformationType::Identity()); // Uniform
				landmarkObservation_fr_p2->setParameterId(0, 999999);

				auto components_fr_lm = std::make_pair(
						landmarkObservation_fr_p1, landmarkObservation_fr_p2);
				landmarkObservation_fr->setParameterId(0, 999999);
				landmarkObservation_fr->initializeComponents(components_fr_lm,
						weights_lm, beta_params);
				optimizer_.addEdge(landmarkObservation_fr);

				// bl
				EdgeSE3PointXYZNonInformativeMixtureELBO *landmarkObservation_bl =
						new EdgeSE3PointXYZNonInformativeMixtureELBO;

				EdgeSE3PointXYZ *landmarkObservation_bl_p1 = new EdgeSE3PointXYZ;
				landmarkObservation_bl_p1->vertices()[0] = optimizer_.vertex(
						i + 1); //from
				landmarkObservation_bl_p1->vertices()[1] = optimizer_.vertex(
						step_ + 5 * j + 5); //to
				landmarkObservation_bl_p1->setMeasurement(
						history_meas_static_.at(i).at(j).bl);
				landmarkObservation_bl_p1->setInformation(
						(1.0/obj_est_stddev_) * g2o::EdgeSE3PointXYZ::InformationType::Identity()); // Gaussian
				landmarkObservation_bl_p1->setParameterId(0, 999999);

				EdgeSE3PointXYZNonInformative *landmarkObservation_bl_p2 =
						new EdgeSE3PointXYZNonInformative;
				landmarkObservation_bl_p2->vertices()[0] = optimizer_.vertex(
						i + 1); //from
				landmarkObservation_bl_p2->vertices()[1] = optimizer_.vertex(
						step_ + 5 * j + 5); //to
				landmarkObservation_bl_p2->setMeasurement(
						history_meas_static_.at(i).at(j).bl);
				landmarkObservation_bl_p2->setInformation(
						10 * g2o::EdgeSE3PointXYZ::InformationType::Identity()); // Uniform
				landmarkObservation_bl_p2->setParameterId(0, 999999);

				auto components_bl_lm = std::make_pair(
						landmarkObservation_bl_p1, landmarkObservation_bl_p2);
				landmarkObservation_bl->setParameterId(0, 999999);
				landmarkObservation_bl->initializeComponents(components_bl_lm,
						weights_lm, beta_params);
				optimizer_.addEdge(landmarkObservation_bl);

				// br
				EdgeSE3PointXYZNonInformativeMixtureELBO *landmarkObservation_br =
						new EdgeSE3PointXYZNonInformativeMixtureELBO;

				EdgeSE3PointXYZ *landmarkObservation_br_p1 = new EdgeSE3PointXYZ;
				landmarkObservation_br_p1->vertices()[0] = optimizer_.vertex(
						i + 1); //from
				landmarkObservation_br_p1->vertices()[1] = optimizer_.vertex(
						step_ + 5 * j + 6); //to
				landmarkObservation_br_p1->setMeasurement(
						history_meas_static_.at(i).at(j).br);
				landmarkObservation_br_p1->setInformation(
						(1.0/obj_est_stddev_) * g2o::EdgeSE3PointXYZ::InformationType::Identity()); // Gaussian
				landmarkObservation_br_p1->setParameterId(0, 999999);

				EdgeSE3PointXYZNonInformative *landmarkObservation_br_p2 =
						new EdgeSE3PointXYZNonInformative;
				landmarkObservation_br_p2->vertices()[0] = optimizer_.vertex(
						i + 1); //from
				landmarkObservation_br_p2->vertices()[1] = optimizer_.vertex(
						step_ + 5 * j + 6); //to
				landmarkObservation_br_p2->setMeasurement(
						history_meas_static_.at(i).at(j).br);
				landmarkObservation_br_p2->setInformation(
						10 * g2o::EdgeSE3PointXYZ::InformationType::Identity()); // Uniform
				landmarkObservation_br_p2->setParameterId(0, 999999);

				auto components_br_lm = std::make_pair(
						landmarkObservation_br_p1, landmarkObservation_br_p2);

				landmarkObservation_br->setParameterId(0, 999999);
				landmarkObservation_br->initializeComponents(components_br_lm,
						weights_lm, beta_params);
				optimizer_.addEdge(landmarkObservation_br);

				if (i == step_) {
					lm_edges.push_back(std::make_pair(j,landmarkObservation_br));
				}
			}

			for (size_t j : history_obs_dynamic_.at(i)) {

				std::vector<double> weights_lm { dynamic_objects_.at(j).getConfidence(),
										1 - dynamic_objects_.at(j).getConfidence() };
				std::vector<double> beta_params = dynamic_objects_.at(j).getBetaParams();

				// fl
				EdgeSE3PointXYZNonInformativeMixtureELBO *landmarkObservation_fl =
						new EdgeSE3PointXYZNonInformativeMixtureELBO;

				EdgeSE3PointXYZ *landmarkObservation_fl_p1 = new EdgeSE3PointXYZ;
				landmarkObservation_fl_p1->vertices()[0] = optimizer_.vertex(i+1); //from
				landmarkObservation_fl_p1->vertices()[1] = optimizer_.vertex(dyn_obj_vertex_ids.at(i).at(j)+1); //to
				landmarkObservation_fl_p1->setMeasurement(history_meas_dynamic_.at(i).at(j).fl);
				landmarkObservation_fl_p1->setInformation(
						(1.0/obj_est_stddev_) * g2o::EdgeSE3PointXYZ::InformationType::Identity()); // Gaussian
				landmarkObservation_fl_p1->setParameterId(0, 999999);

				EdgeSE3PointXYZNonInformative *landmarkObservation_fl_p2 =
						new EdgeSE3PointXYZNonInformative;
				landmarkObservation_fl_p2->vertices()[0] = optimizer_.vertex(i+1); //from
				landmarkObservation_fl_p2->vertices()[1] = optimizer_.vertex(dyn_obj_vertex_ids.at(i).at(j)+1); //to
				landmarkObservation_fl_p2->setMeasurement(history_meas_dynamic_.at(i).at(j).fl);
				landmarkObservation_fl_p2->setInformation(
						10 * g2o::EdgeSE3PointXYZ::InformationType::Identity()); // Uniform
				landmarkObservation_fl_p2->setParameterId(0, 999999);

				auto components_fl_lm = std::make_pair(landmarkObservation_fl_p1,
						landmarkObservation_fl_p2);
				landmarkObservation_fl->setParameterId(0, 999999);
				landmarkObservation_fl->initializeComponents(components_fl_lm,
						weights_lm, beta_params);
				optimizer_.addEdge(landmarkObservation_fl);

				// fr
				EdgeSE3PointXYZNonInformativeMixtureELBO *landmarkObservation_fr =
						new EdgeSE3PointXYZNonInformativeMixtureELBO;

				EdgeSE3PointXYZ *landmarkObservation_fr_p1 = new EdgeSE3PointXYZ;
				landmarkObservation_fr_p1->vertices()[0] = optimizer_.vertex(
						i + 1); //from
				landmarkObservation_fr_p1->vertices()[1] = optimizer_.vertex(
						dyn_obj_vertex_ids.at(i).at(j)+2); //to
				landmarkObservation_fr_p1->setMeasurement(
						history_meas_dynamic_.at(i).at(j).fr);
				landmarkObservation_fr_p1->setInformation(
						(1.0/obj_est_stddev_) * g2o::EdgeSE3PointXYZ::InformationType::Identity()); // Gaussian
				landmarkObservation_fr_p1->setParameterId(0, 999999);

				EdgeSE3PointXYZNonInformative *landmarkObservation_fr_p2 =
						new EdgeSE3PointXYZNonInformative;
				landmarkObservation_fr_p2->vertices()[0] = optimizer_.vertex(
						i + 1); //from
				landmarkObservation_fr_p2->vertices()[1] = optimizer_.vertex(
						dyn_obj_vertex_ids.at(i).at(j)+2); //to
				landmarkObservation_fr_p2->setMeasurement(
						history_meas_dynamic_.at(i).at(j).fr);
				landmarkObservation_fr_p2->setInformation(
						10 * g2o::EdgeSE3PointXYZ::InformationType::Identity()); // Uniform
				landmarkObservation_fr_p2->setParameterId(0, 999999);

				auto components_fr_lm = std::make_pair(
						landmarkObservation_fr_p1, landmarkObservation_fr_p2);
				landmarkObservation_fr->setParameterId(0, 999999);
				landmarkObservation_fr->initializeComponents(components_fr_lm,
						weights_lm, beta_params);
				optimizer_.addEdge(landmarkObservation_fr);

				// bl
				EdgeSE3PointXYZNonInformativeMixtureELBO *landmarkObservation_bl =
						new EdgeSE3PointXYZNonInformativeMixtureELBO;

				EdgeSE3PointXYZ *landmarkObservation_bl_p1 = new EdgeSE3PointXYZ;
				landmarkObservation_bl_p1->vertices()[0] = optimizer_.vertex(
						i + 1); //from
				landmarkObservation_bl_p1->vertices()[1] = optimizer_.vertex(
						dyn_obj_vertex_ids.at(i).at(j)+3); //to
				landmarkObservation_bl_p1->setMeasurement(
						history_meas_dynamic_.at(i).at(j).bl);
				landmarkObservation_bl_p1->setInformation(
						(1.0/obj_est_stddev_) * g2o::EdgeSE3PointXYZ::InformationType::Identity()); // Gaussian
				landmarkObservation_bl_p1->setParameterId(0, 999999);

				EdgeSE3PointXYZNonInformative *landmarkObservation_bl_p2 =
						new EdgeSE3PointXYZNonInformative;
				landmarkObservation_bl_p2->vertices()[0] = optimizer_.vertex(
						i + 1); //from
				landmarkObservation_bl_p2->vertices()[1] = optimizer_.vertex(
						dyn_obj_vertex_ids.at(i).at(j)+3); //to
				landmarkObservation_bl_p2->setMeasurement(
						history_meas_dynamic_.at(i).at(j).bl);
				landmarkObservation_bl_p2->setInformation(
						10 * g2o::EdgeSE3PointXYZ::InformationType::Identity()); // Uniform
				landmarkObservation_bl_p2->setParameterId(0, 999999);

				auto components_bl_lm = std::make_pair(
						landmarkObservation_bl_p1, landmarkObservation_bl_p2);
				landmarkObservation_bl->setParameterId(0, 999999);
				landmarkObservation_bl->initializeComponents(components_bl_lm,
						weights_lm, beta_params);
				optimizer_.addEdge(landmarkObservation_bl);

				// br
				EdgeSE3PointXYZNonInformativeMixtureELBO *landmarkObservation_br =
						new EdgeSE3PointXYZNonInformativeMixtureELBO;

				EdgeSE3PointXYZ *landmarkObservation_br_p1 = new EdgeSE3PointXYZ;
				landmarkObservation_br_p1->vertices()[0] = optimizer_.vertex(
						i + 1); //from
				landmarkObservation_br_p1->vertices()[1] = optimizer_.vertex(
						dyn_obj_vertex_ids.at(i).at(j)+4); //to
				landmarkObservation_br_p1->setMeasurement(
						history_meas_dynamic_.at(i).at(j).br);
				landmarkObservation_br_p1->setInformation(
						(1.0/obj_est_stddev_) * g2o::EdgeSE3PointXYZ::InformationType::Identity()); // Gaussian
				landmarkObservation_br_p1->setParameterId(0, 999999);

				EdgeSE3PointXYZNonInformative *landmarkObservation_br_p2 =
						new EdgeSE3PointXYZNonInformative;
				landmarkObservation_br_p2->vertices()[0] = optimizer_.vertex(
						i + 1); //from
				landmarkObservation_br_p2->vertices()[1] = optimizer_.vertex(
						dyn_obj_vertex_ids.at(i).at(j)+4); //to
				landmarkObservation_br_p2->setMeasurement(
						history_meas_dynamic_.at(i).at(j).br);
				landmarkObservation_br_p2->setInformation(
						10 * g2o::EdgeSE3PointXYZ::InformationType::Identity()); // Uniform
				landmarkObservation_br_p2->setParameterId(0, 999999);

				auto components_br_lm = std::make_pair(
						landmarkObservation_br_p1, landmarkObservation_br_p2);
				landmarkObservation_br->setParameterId(0, 999999);
				landmarkObservation_br->initializeComponents(components_br_lm,
						weights_lm, beta_params);
				optimizer_.addEdge(landmarkObservation_br);

				if (i == step_) {
					lm_edges.push_back(std::make_pair(num_static_objs_+j,landmarkObservation_br));
				}

				//if (i == step_) lm_mmniedges.push_back(landmarkObservation);

				/////////////////////

				// add motion edge
				if (i > window_start) {
					EdgeSE3Mixture *motion_edge = new EdgeSE3Mixture;

					EdgeSE3Orig *motion_edge_p1 = new EdgeSE3Orig;
					motion_edge_p1->vertices()[0] = optimizer_.vertex(dyn_obj_vertex_ids.at(i-1).at(j)); //from
					motion_edge_p1->vertices()[1] = optimizer_.vertex(dyn_obj_vertex_ids.at(i).at(j)); //to
					motion_edge_p1->setMeasurement(Isometry3ByAngle(history_meas_dynamic_.at(i).at(j).vel, 0));
					motion_edge_p1->setInformation((1.0/obj_est_stddev_) * g2o::EdgeSE3Orig::InformationType::Identity()); // Gaussian

					EdgeSE3Orig *motion_edge_p2 = new EdgeSE3Orig;
					motion_edge_p2->vertices()[0] = optimizer_.vertex(dyn_obj_vertex_ids.at(i-1).at(j)); //from
					motion_edge_p2->vertices()[1] = optimizer_.vertex(dyn_obj_vertex_ids.at(i).at(j)); //to
					motion_edge_p2->setMeasurement(Isometry3ByAngle(Eigen::Vector3d(0,0,0), 0));
					//motion_edge_p2->setMeasurement(Isometry3ByAngle(history_meas_dynamic_.at(i).at(j).vel, 0));
					motion_edge_p2->setInformation((1.0/obj_est_stddev_) * g2o::EdgeSE3Orig::InformationType::Identity()); // Gaussian

					auto components_motion = std::vector<EdgeSE3Orig*>{motion_edge_p1, motion_edge_p2};
					std::vector<double> weights_motion { 0.5, 0.5 };
					motion_edge->initializeComponents(components_motion, weights_motion);
					optimizer_.addEdge(motion_edge);
				}
			}
		}

		/*
		for (int i(window_start); i<=step_; ++i) {
			for (size_t j : history_obs_dynamic_.at(i)) {

			}
		}
		*/


		// Add dynamic object motion edges

		if (save_) {
			optimizer_.save(
					std::string("g2o/pg_before_step" + std::to_string(step_) +
					"_iter" + std::to_string(iter) + ".g2o").data());
		}

		VertexSE3* robot = robots.back();

		optimizer_.setVerbose(false);

		auto pose_before = robot->estimate().matrix();
		//std::cout << "Robot pose before opt: "  << std::endl << pose_before << std::endl;

		if (debug_) std::cout << "    Optimizing..." << endl;
		optimizer_.initializeOptimization();
		optimizer_.optimize(20);
		if (debug_) std::cout << "    done" << endl;

		auto pose_after = robot->estimate().matrix();
		//std::cout << "Robot pose after opt: "  << std::endl << pose_after << std::endl;


		for (auto& e : lm_edges) {
			int id = e.first;
			EdgeSE3PointXYZNonInformativeMixtureELBO* edge = e.second;
			if (id < num_static_objs_) {
				static_objects_.at(id).r_in = edge->getELBOWeights().at(0);
			} else {
				id = id - num_static_objs_;
				dynamic_objects_.at(id).r_in = edge->getELBOWeights().at(0);
			}
		}


		if (iter == max_iter_-1) {
			// update static object estimated pose
			for (int j(0); j<static_obj_vertex_ids.size(); j++) {
				auto& obj = static_obj_vertex_ids.at(j);
				if (obj.at(0) == -1) continue;

				VertexPointXYZ* lm_fl = static_cast<VertexPointXYZ*>(optimizer_.vertex(obj.at(0)));
				VertexPointXYZ* lm_fr = static_cast<VertexPointXYZ*>(optimizer_.vertex(obj.at(1)));
				VertexPointXYZ* lm_bl = static_cast<VertexPointXYZ*>(optimizer_.vertex(obj.at(2)));
				VertexPointXYZ* lm_br = static_cast<VertexPointXYZ*>(optimizer_.vertex(obj.at(3)));

				auto lm_fl_new_est = lm_fl->estimate();
				static_objects_.at(j).fl_est = Pose2(lm_fl_new_est(0), lm_fl_new_est(1), 0);

				auto lm_fr_new_est = lm_fr->estimate();
				static_objects_.at(j).fr_est = Pose2(lm_fr_new_est(0), lm_fr_new_est(1), 0);

				auto lm_bl_new_est = lm_bl->estimate();
				static_objects_.at(j).bl_est = Pose2(lm_bl_new_est(0), lm_bl_new_est(1), 0);

				auto lm_br_new_est = lm_br->estimate();
				static_objects_.at(j).br_est = Pose2(lm_br_new_est(0), lm_br_new_est(1), 0);
			}
		}

		// update dynamic object estimated pose
		for (size_t j : history_obs_dynamic_.at(step_)) {
			VertexSE3* obj = static_cast<VertexSE3*>(optimizer_.vertex(dyn_obj_vertex_ids.at(step_).at(j)));

			auto obj_new_est = obj->estimate().matrix();
			auto new_rot = obj_new_est.block<3, 3>(0, 0);
			auto new_trans = obj_new_est.block<3, 1>(0, 3);
			auto new_euler = new_rot.eulerAngles(0, 1, 2);

			Pose2 new_est = Pose2(new_trans(0), new_trans(1), new_euler(2));
			dynamic_objects_.at(j).set_est_pose(new_est);
		}



		if (save_) {
			optimizer_.save(
					std::string("g2o/pg_after_step" + std::to_string(step_) +
					"_iter" + std::to_string(iter) + ".g2o").data());
		}

		auto new_rot = pose_after.block<3, 3>(0, 0);
		auto new_trans = pose_after.block<3, 1>(0, 3);
		auto new_euler = new_rot.eulerAngles(0, 1, 2);

		robot_pose_est_ = Pose2(new_trans(0), new_trans(1), new_euler(2));
		if (debug_) robot_pose_est_.print("    New robot pose: ");

		for (auto &p : lm_edges) {
			int obj_id = (p.second->vertex(1)->id() - step_ - 2) / 5;
			if (obj_id >= num_static_objs_) break;
			static_object_status_.at(obj_id) = p.second->getBestComponent();
		}
		if (debug_)
			std::cout << "    Number of inactive observations: " << num_inactive_objects() << endl;

		/*
		for (size_t j : observed_objs) {
			VertexPointXYZ* lmk =  static_cast<VertexPointXYZ*>(optimizer_.vertex(j+1));
			HyperGraph::EdgeSet edges = lmk->edges();
			for (auto* e : edges) {
				std::cout << e->elementType() << std::endl;
			}
		}
		*/
		if (iter == max_iter_ - 1) {
			draw();
		}

		optimizer_.clear();
	}

	std::vector<size_t> get_dyn_obj_vertex_intervals(int window_start, int window_size) {
		std::vector<size_t> dyn_obj_vertex_ids(num_dynamic_objs_+1, 0);
		dyn_obj_vertex_ids.at(0) = window_start + window_size + num_static_objs_ + 1;
		if (num_dynamic_objs_ > 0) {
			// calculate indices for dynamic object temporal vertices
			for (int i(window_start); i<=step_; ++i) {
				for (size_t j : history_obs_dynamic_.at(i)) {
					dyn_obj_vertex_ids.at(j+1) += 1;
				}
			}
			for (size_t j(1); j<=num_dynamic_objs_; j++) {
				dyn_obj_vertex_ids.at(j) += dyn_obj_vertex_ids.at(j-1);
			}
		}
		return dyn_obj_vertex_ids;
	}

	std::map<int,std::map<size_t,size_t>> get_dyn_obj_temporal_obs(int window_start, int window_size) {
		std::map<int,std::map<size_t,size_t>> dyn_obj_vertex_ids = {};
		int idx = window_start + window_size + 5*(int)num_static_objs_ + 1;

		for (int i(window_start); i<=step_; ++i) {
			// calculate indices for dynamic object temporal vertices
			for (size_t j : history_obs_dynamic_.at(i)) {
				dyn_obj_vertex_ids[i][j] = idx;
				idx += 5;
			}
		}
		return dyn_obj_vertex_ids;
	}

	//std::vector<size_t> interval_to_vector(const std::vector<size_t>& intervals) {

	//}

	void update_canvas_size(Pose2 pose) {
		if (pose.x() > canvas_max_x_) canvas_max_x_ = pose.x();
		if (pose.x() < canvas_min_x_) canvas_min_x_ = pose.x();
		if (pose.y() > canvas_max_y_) canvas_max_y_ = pose.y();
		if (pose.y() < canvas_min_y_) canvas_min_y_ = pose.y();
		// std::cout << canvas_max_x_ << "  "<< canvas_min_x_ << "  "<< canvas_max_y_ << "  "<< canvas_min_y_ << std::endl;
	}

	// void draw(const std::vector<EdgeSE3PointXYZNonInformativeMixture*>& mmedges) {
	void draw() {
		double scale = 100; // cm
		double buffer = 100; // 1m buffer
		int width = (canvas_max_x_ - canvas_min_x_) * scale + 2 * buffer;
		int height = (canvas_max_y_ - canvas_min_y_) * scale + 2 * buffer;
		int cx = buffer - canvas_min_x_ * scale;
		int cy = buffer + canvas_max_y_ * scale;

		cv::Mat image(height, width, CV_8UC3, cv::Scalar(255,255,255));

		draw_triangle(image, robot_pose_true_, 0.3, 0.5, scale, cx, cy, cv::Scalar(0,150,0));
		//draw_triangle(image, robot_pose_odom_, 0.3, 0.5, scale, cx, cy, cv::Scalar(200,150,0));
		draw_triangle(image, robot_pose_est_, 0.3, 0.5, scale, cx, cy, cv::Scalar(0,150,200));

		cv::Point2f robot_pos_est_cv(robot_pose_est_.x()*scale+cx, -robot_pose_est_.y()*scale+cy);

		// draw static objects
		for (size_t j(0); j<num_static_objs_; ++j) {
			const Object& obj = static_objects_.at(j);
			int status = static_object_status_.at(j);
			int reconsted = static_object_reconsted_.at(j);
			double conf = obj.getConfidenceBK();
			cv::Scalar lmk_true_color;
			cv::Scalar lmk_est_color;

			if (!reconsted) {
				lmk_true_color = cv::Scalar(0,150,0);
				lmk_est_color = cv::Scalar(0,150,200);
			} else if (reconsted) {
				lmk_true_color = cv::Scalar(0,150,0);
				lmk_est_color = cv::Scalar(200,150,0);
			}

			cv::Point2f pos_true(obj.pose_true.x()*scale+cx, -obj.pose_true.y()*scale+cy);
			cv::Point2f pos_est(obj.pose_est.x()*scale+cx, -obj.pose_est.y()*scale+cy);

			draw_rectangle(image, obj, true, false, scale, cx, cy, lmk_true_color);
			draw_rectangle(image, obj, false, false, scale, cx, cy, lmk_est_color);

			cv::circle(image,pos_true,4,lmk_true_color,cv::FILLED,cv::LINE_AA);
			cv::circle(image,pos_est,4,lmk_est_color,cv::FILLED,cv::LINE_AA);
			cv::line(image,pos_true,pos_est,cv::LINE_AA);

			std::string text_id = "Obj ID: " + std::to_string(j);
			cv::Point2f pos_id = pos_true + cv::Point2f(-90, 0);
			std::string text_conf = "Obj Score: " + std::to_string(conf * 100).substr(0,5) + "%";
			cv::Point2f pos_conf = pos_true + cv::Point2f(-90, -14);
			cv::putText(image, text_id, pos_id, cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0,0,0), 1, cv::LINE_AA);
			cv::putText(image, text_conf, pos_conf, cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0,0,0), 1, cv::LINE_AA);

			cv::line(image,robot_pos_est_cv,pos_est,lmk_est_color,1,cv::LINE_AA);
		}


		for (size_t j(0); j<num_dynamic_objs_; ++j) {
			const Object& obj = dynamic_objects_.at(j);
			double conf = obj.getConfidenceBK();

			cv::Scalar lmk_true_color = cv::Scalar(0,150,0);
			cv::Scalar lmk_est_color = cv::Scalar(0,150,200);
			cv::Scalar lmk_track_color = cv::Scalar(0,0,0);

			cv::Point2f pos_true(obj.pose_true.x()*scale+cx, -obj.pose_true.y()*scale+cy);
			cv::Point2f pos_est(obj.pose_est.x()*scale+cx, -obj.pose_est.y()*scale+cy);
			cv::Point2f pos_track(obj.pose_track.x()*scale+cx, -obj.pose_track.y()*scale+cy);

			draw_rectangle(image, obj, true, false, scale, cx, cy, lmk_true_color);
			//draw_rectangle(image, obj, false, true, scale, cx, cy, lmk_track_color);
			draw_rectangle(image, obj, false, false, scale, cx, cy, lmk_est_color);

			cv::circle(image,pos_true,4,lmk_true_color,cv::FILLED,cv::LINE_AA);
			//cv::circle(image,pos_track,4,lmk_track_color,cv::FILLED,cv::LINE_AA);
			cv::circle(image,pos_est,4,lmk_est_color,cv::FILLED,cv::LINE_AA);
			cv::line(image,pos_true,pos_est,cv::LINE_AA);

			std::string text_id = "Obj ID: " + std::to_string(num_static_objs_ + j);
			cv::Point2f pos_id = pos_true + cv::Point2f(-90, 0);
			std::string text_conf = "Obj Score: " + std::to_string(conf * 100).substr(0,5) + "%";
			cv::Point2f pos_conf = pos_true + cv::Point2f(-90, -14);
			std::string text_dyn = "Dynamic";
			cv::Point2f pos_dyn = pos_true + cv::Point2f(-90, 14);
			cv::putText(image, text_id, pos_id, cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0,0,0), 1, cv::LINE_AA);
			cv::putText(image, text_conf, pos_conf, cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0,0,0), 1, cv::LINE_AA);
			cv::putText(image, text_dyn, pos_dyn, cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0,0,0), 1, cv::LINE_AA);
		}


		cv::Point2f robot_pos_est(robot_pose_est_.x()*scale+cx, -robot_pose_est_.y()*scale+cy);

		cv::imshow("Env", image);
		// cv::imwrite("/home/jqian/Downloads/steps/Fig_" + std::to_string(step_) + ".png", image);
		cv::waitKey(100);
	}

	void draw_triangle(cv::Mat& img, Pose2 p, double base, double height, double scale, int cx, int cy, cv::Scalar color) {
		Pose2 l = p + Pose2(0,base/2,0);
		Pose2 r = p + Pose2(0,-base/2,0);
		Pose2 f = p + Pose2(height,0,0);
		cv::Point2i lp(l.x()*scale+cx, -l.y()*scale+cy);
		cv::Point2i rp(r.x()*scale+cx, -r.y()*scale+cy);
		cv::Point2i fp(f.x()*scale+cx, -f.y()*scale+cy);

		cv::Point2i triangle[1][3];
		triangle[0][0]  = lp;
		triangle[0][1]  = rp;
		triangle[0][2]  = fp;
		const cv::Point2i* ppt[1] = { triangle[0] };

		int npt[] = { 3 };

		cv::fillPoly(img,ppt,npt,1,color,cv::LINE_AA);
	}

	void draw_rectangle(cv::Mat& img, Object obj, bool is_gt, bool is_track, double scale, int cx, int cy, cv::Scalar color) {
		Pose2 fl = is_gt ? obj.fl_true : obj.fl_est;
		Pose2 fr = is_gt ? obj.fr_true : obj.fr_est;
		Pose2 bl = is_gt ? obj.bl_true : obj.bl_est;
		Pose2 br = is_gt ? obj.br_true : obj.br_est;
		if (is_track) {
			fl = obj.fl_track;
			fr = obj.fr_track;
			bl = obj.bl_track;
			br = obj.br_track;
		}

		cv::Point2i flp(fl.x()*scale+cx, -fl.y()*scale+cy);
		cv::Point2i frp(fr.x()*scale+cx, -fr.y()*scale+cy);
		cv::Point2i blp(bl.x()*scale+cx, -bl.y()*scale+cy);
		cv::Point2i brp(br.x()*scale+cx, -br.y()*scale+cy);

		cv::Point2i rect[1][4];
		rect[0][0]  = flp;
		rect[0][1]  = frp;
		rect[0][2]  = brp;
		rect[0][3]  = blp;
		const cv::Point* ppt[1] = { rect[0] };

		int npt[] = { 4 };

		//cv::fillPoly(img,ppt,npt,1,color,cv::LINE_AA);
		cv::polylines(img,ppt,npt,1,true,color,1,cv::LINE_AA);
	}

	SparseOptimizer optimizer_;

	size_t max_iter_;
	double meas_stddev_;
	double max_meas_range_;

	double obj_est_stddev_ = 0.2;

	double motion_stddev_ = 0.03;

	double reconst_lb = 0.2;

	// 0 if active, 1 if inactive
	std::vector<int> static_object_status_;
	std::vector<int> static_object_reconsted_;

	std::vector<int> dynamic_object_status_;

	size_t num_static_objs_;
	size_t num_dynamic_objs_;

	Pose2 robot_pose_true_;
	Pose2 robot_pose_est_;
	Pose2 robot_pose_est_init_;
	Pose2 robot_pose_odom_;

	std::vector<Pose2> history_odom_ = {};
	std::vector<Pose2> history_pose_est_ = {};

	std::vector<std::map<size_t, RectObs> > history_meas_static_ = {};
	std::vector<std::vector<size_t> > history_obs_static_ = {};

	std::vector<std::map<size_t, RectObs> > history_meas_dynamic_ = {};
	std::vector<std::vector<size_t> > history_obs_dynamic_ = {};

	double canvas_max_x_;
	double canvas_min_x_;
	double canvas_max_y_;
	double canvas_min_y_;

	int step_ = 0;
	int window_size_ = 10;

	bool debug_ = false;
	bool save_ = false;

	const char* path1 = "Iterations_dynamic_example.csv"; //VC
	ofstream file_iter;

	const char* path3 = "pose_error.csv"; //VC
	ofstream pose_error_file;

};


int main() {

	int max_window_size = 6;
	int increment = 5;
	int num_obj = 10;
	double stationary_prob = 0.3; // has not moved
	double static_prob = 0.75; // not moving

	const char* path2 = "Steps_dynamic_example.csv"; //VC
	ofstream file_steps;

	file_steps.open(path2, ios::in | ios::app); //VC
	// file_steps << "Step" << "," << "Object ID" << "," << "Stationarity" << "," << "True Stationarity" << std::endl; //VC

	for (int window_size(5); window_size<max_window_size; window_size+=increment) {

		//std::srand(static_cast<unsigned int>(std::time(NULL)));

		Environment env(20, 0.03, 100, window_size); // iterations std_dev max_meas_range window_size

		env.set_robot_pose(1.0, 0.0, 0.0, 3.1, 1.15, 0.05, 0.0, -3.0);

		std::uniform_real_distribution<> dis_pos(-5.0, 5.0);
		std::uniform_real_distribution<> dis_stationary(0.0, 1.0);
		std::uniform_real_distribution<> dis_dist(0.4, 1.0);
		std::uniform_real_distribution<> dis_angle(0.0, 2.0 * M_PI);

		int arraySign[2] = { -1, 1 };

		for (int m(0); m<num_obj; ++m) {

			double obj_x = dis_pos(generator);
			double obj_y = dis_pos(generator);
			double stationary = dis_stationary(generator);

			if (stationary >= stationary_prob && stationary < static_prob) {
				double dist_error = dis_dist(generator);
				double angle_error = dis_angle(generator);
				double dx = dist_error * cos(angle_error);
				double dy = dist_error * sin(angle_error);
				env.add_object(obj_x, obj_y, 0, obj_x + dx, obj_y + dy, 0);
			} else if (stationary >=  static_prob) {
				env.add_object(obj_x, obj_y, 0, obj_x, obj_y, 0, 0);
			} else {
				env.add_object(obj_x, obj_y, 0, obj_x, obj_y, 0);
			}
		}

		for (int i(0); i < 50; i++) {
			env.step();
			env.optimize();
			for (int j=0; j < env.static_objects_.size(); j++){
				file_steps << i << "," << j << "," << env.static_objects_.at(j).getConfidence() << "," << 1 << std::endl;
			}
			
			for (int j=0; j < env.dynamic_objects_.size(); j++){
				file_steps << i << "," << j + env.static_objects_.size()<< "," << env.dynamic_objects_.at(j).getConfidence() << "," << 0 << std::endl;
			}
		}
		
		file_steps.close();

	}
	return 0;
}
