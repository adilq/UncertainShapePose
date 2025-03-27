//
// A pose representing a point and direction in 2D space
//
#pragma once

#include <string>
#include <vector>
#include <Eigen/Dense>

class Pose2 {
 public:
    // Constructors
    Pose2();
    Pose2(double x, double y, double yaw);
    Pose2(double x, double y, double z, double yaw);
    Pose2(Eigen::Matrix4f m);

    // Operators
    Pose2 operator+(const Pose2& innovation);
    Pose2 operator-(const Pose2& innovation) const;
    Pose2 operator* (const double scale);
    Pose2 operator=(const Pose2& other) const;
    bool operator==(const Pose2 other);
    void operator=(const Pose2 &position);
    friend std::ostream& operator<<(std::ostream &os, const Pose2 &pos);
        // Print function
    void print(std::string str = "") const;
    void setX(double x);
    void setY(double y);
    void setZ(double z);
    void setYaw(double yaw);
    void setVelocity(double vel);
    // Get funtions
    const double& x() const;
    const double& y() const;
    const double& z() const;
    const double& yaw() const;
    double& x();
    double& y();
    double& z();
    double& yaw();
    const double& velocity() const;
    std::vector<double> position() const;
    Eigen::Vector3d position_eigen() const;
    Eigen::Isometry3d isometry3() const;
    void rotate(double ang_rad);
    Pose2 difference(Pose2 innovation) const;
    double norm() const;
    double distance(const Pose2 &other) const;
    double squared_norm() const;
    double squared_distance(const Pose2 &other) const;
    double dot(const Pose2 &other) const;

 private:
    double wrapAngle(double angle);

    double _x;
    double _y;
    double _z;
    double _yaw;

    double _velocity;
};

class Pose2HashFunction {
 public:
    size_t operator()(const Pose2& pose2) const {
        size_t x_hash = std::hash<double>{}(pose2.x());
        size_t y_hash = std::hash<double>{}(pose2.y());
        size_t z_hash = std::hash<double>{}(pose2.y());
        size_t yaw_hash = std::hash<double>{}(pose2.yaw());
        return (x_hash ^ (y_hash << 1) ^ (z_hash << 2) >> 1) ^ (yaw_hash << 1);
    }
};
