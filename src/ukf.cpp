#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.67; // 2*sigma = 3.33 Assume 95% cars decrease speed from 60 km/h to 0 in 5 sec (and vice versa)

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.08; //Assume 95% cars change angle ~10 grad/s^2

  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
   * End DO NOT MODIFY section for measurement noise values
   */

  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
      use_laser_ = false;
      use_radar_ = true;

      // Get radar measurements.
      double rho = meas_package.raw_measurements_[0]; // Radius measurement (m).
      double phi = meas_package.raw_measurements_[1]; // Angle measurement (rad).
      double rho_d = meas_package.raw_measurements_[2]; // Velocity (m/s).

      // Transform radar measurements to state vector.
      // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
      x_ << rho*cos(phi), rho*sin(phi), rho_d, phi, 0;

      // Set state covariance matrix.
      P_.fill(0.0);
      P_(0, 0) = std_radr_*std_radr_*std_radphi_*std_radphi_/(std_radr_*std_radr_+std_radphi_*std_radphi_); // Variance of product of two normal distributions.
      P_(1, 1) = P_(0, 0);
      P_(2, 2) = std_radrd_*std_radrd_;
      P_(3, 3) = std_radphi_ * std_radphi_;
      P_(4, 4) = P_(3, 3);
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];
      P_.fill(0.0);
      P_(0, 0) = std_laspx_*std_laspx_;
      P_(1, 1) = std_laspy_*std_laspy_;
    } else {
      throw std::runtime_error("Invalid sensor type");
    }
    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
  } else {
    double dt = static_cast<double>(time_us_ - meas_package.timestamp_)*1e-6;
    Prediction(dt); // Predict state based on process model.
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
      UpdateRadar(meas_package);
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
      UpdateLidar(meas_package);
    else
      throw std::runtime_error("Invalid sensor type");
  }
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location.
   * Modify the state vector, x_. Predict sigma points, the state,
   * and the state covariance matrix.
   */
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}