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
    n_x_ = x_.size();
    n_aug_ = n_x_ + 2; // Take into account process noise.
    lambda_ = 3 - n_aug_; // Define spreading parameter.
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

  // Augmented state vector.
  VectorXd x_aug = VectorXd(7);

  // Augmented state covariance.
  MatrixXd P_aug = MatrixXd(7, 7);

  // Sigma points matrix.
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // Copy sate vector and leave tail as zero.
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_+1) = 0;

  // Copy state covariance and augment it with process noise.
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_*std_a_;
  P_aug(n_x_+1, n_x_+1) = std_yawdd_*std_yawdd_;

  // Square root matrix of covariance matrix.
  MatrixXd P_aug_sqrt = P_aug.llt().matrixL();

  // Create augmented sigma points.
  Xsig_aug.col(0) = x_aug;
  for (int c = 1; c <= n_aug_; c++){
    Xsig_aug.col(c) = x_aug + sqrt(lambda_+n_aug_)*P_aug_sqrt.col(c-1);
    Xsig_aug.col(c + n_aug_) = x_aug - sqrt(lambda_+n_aug_)*P_aug_sqrt.col(c-1);
  }

  // Set weights for each sigma point.
  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = lambda_/(lambda_+n_aug_);
  weights_.tail(2*n_aug_).fill(0.5/(lambda_ + n_aug_));

  // Matrix with predicted sigma points as columns.
  MatrixXd Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  // Predict sigma points (map them on target Gaussian after passing through non-linear function).
  for (int i = 0; i < 2*n_aug_+1; i++){ // Go through columns
      double px = Xsig_aug(0, i);
      double py = Xsig_aug(1, i);
      double v =  Xsig_aug(2, i);
      double phi =  Xsig_aug(3, i);
      double phid = Xsig_aug(4, i);
      double nu_a = Xsig_aug(5, i);
      double nu_psid = Xsig_aug(6, i);
      if (phid > 1e-3){ // Check if phi velocity is zero.
        Xsig_pred_(0, i) = px + (v/phid)*(sin(phi+phid*delta_t)-sin(phi));
        Xsig_pred_(1, i) = py + (v/phid)*(-cos(phi+phid*delta_t)+cos(phi));
      } else {
        Xsig_pred_(0, i) = px + v*cos(phi)*delta_t;
        Xsig_pred_(1, i) = py + v*sin(phi)*delta_t;
      }
      Xsig_pred_(2, i) = v;
      Xsig_pred_(3, i) = phi + phid*delta_t;
      Xsig_pred_(4, i) = phid;

      Xsig_pred_(0, i) += 0.5*delta_t*delta_t*cos(phi)*nu_a;
      Xsig_pred_(1, i) += 0.5*delta_t*delta_t*sin(phi)*nu_a;
      Xsig_pred_(2, i) += delta_t * nu_a;
      Xsig_pred_(3, i) += 0.5*delta_t*delta_t*nu_psid;
      Xsig_pred_(4, i) += delta_t * nu_psid;
  }

  // Predict state.
  x_ = (Xsig_pred_.array().rowwise() * weights_.transpose().array()).rowwise().sum();

  // Predict state covariance matrix.
  P_.fill(0);
  for (int i = 0; i < 2*n_aug_+1; i++)
  {
    VectorXd diff = Xsig_pred_.col(i) - x_;
    while(diff(3)>M_PI) diff(3) -= 2*M_PI; // Angle should be [-pi, pi].
    while(diff(3)<-M_PI) diff(3) += 2*M_PI;  // Angle should be [-pi, pi].
    P_ += weights_(i)*diff*diff.transpose();
  }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  int n_z = 2; // Dimension of lidar measurement space (x, y).

  // Matrix for sigma points in measurement space.
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  VectorXd z_pred = VectorXd(n_z); // Predicted measurement.

  MatrixXd S = MatrixXd(n_z,n_z); // Measurement covariance.

  // Transform sigma points into measurement space.
  for (int i = 0; i < 2*n_aug_+1; i++){
      Zsig(0, i) = Xsig_pred_(0, i);
      Zsig(1, i) = Xsig_pred_(1, i);
  }

  // Calculate predicted measurement.
  z_pred = (Zsig.array().rowwise() * weights_.transpose().array()).rowwise().sum();

  // Calculate covariance of predicted measurements.
  S.fill(0);
  for (int i = 0; i < 2*n_aug_+1; i++){
      VectorXd diff = Zsig.col(i) - z_pred;
      S += weights_(i)*diff*diff.transpose();
  }

  MatrixXd R = MatrixXd(n_z, n_z);
  R  << std_laspx_*std_laspx_, 0,
  std_laspy_*std_laspy_, 0,
  S = S + R;

  // Cross correlation matrix Tc.
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  Tc.fill(0);
  for (int i = 0; i < 2*n_aug_+1; i++){
      VectorXd diff_x = Xsig_pred_.col(i) - x_;
      VectorXd diff_z = Zsig.col(i) - z_pred;
      Tc += weights_(i) * diff_x * diff_z.transpose();
  }

  // Kalman gain K.
  MatrixXd K(n_x_, n_z);
  K = Tc * S.inverse();

  // Update state mean and covariance matrix.
  VectorXd z(2);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];
  VectorXd z_diff = z-z_pred;

  x_ = x_ + K*z_diff;

  P_ = P_ - K*S*K.transpose();
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  int n_z = 3; // Dimension of radar measurement space (rho, phi, rho_d).

  // Matrix for sigma points in measurement space.
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  VectorXd z_pred = VectorXd(n_z); // Predicted measurement.

  MatrixXd S = MatrixXd(n_z,n_z); // Measurement covariance.

  // Transform sigma points into measurement space.
  for (int i = 0; i < 2*n_aug_+1; i++){
      double p_x = Xsig_pred_(0, i);
      double p_y = Xsig_pred_(1, i);
      double v =  Xsig_pred_(2, i);
      double psi = Xsig_pred_(3, i);
      double psi_d = Xsig_pred_(4, i);

      double rho = sqrtf(p_x*p_x + p_y*p_y);
      double phi  = atan2(p_y, p_x);
      double rho_d = (p_x*cos(psi)*v + p_y*sin(psi)*v)/rho;

      Zsig(0, i) = rho;
      Zsig(1, i) = phi;
      Zsig(2, i) = rho_d;
  }

  // Calculate predicted measurement.
  z_pred = (Zsig.array().rowwise() * weights_.transpose().array()).rowwise().sum();

  // Calculate covariance of predicted measurements.
  S.fill(0);
  for (int i = 0; i < 2*n_aug_+1; i++){
      VectorXd diff = Zsig.col(i) - z_pred;
      while (diff(1)>M_PI) diff(1) -= 2*M_PI;  // Angle should be [-pi, pi].
      while (diff(1)<-M_PI) diff(1) += 2*M_PI;  // Angle should be [-pi, pi].
      S += weights_(i)*diff*diff.transpose();
  }

  MatrixXd R = MatrixXd(n_z, n_z);
  R  << std_radr_*std_radr_, 0, 0,
  0, std_radphi_*std_radphi_, 0,
  0, 0, std_radrd_*std_radrd_;
  S = S + R;

  // Cross correlation matrix Tc.
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  Tc.fill(0);
  for (int i = 0; i < 2*n_aug_+1; i++){
      VectorXd diff_x = Xsig_pred_.col(i) - x_;
      VectorXd diff_z = Zsig.col(i) - z_pred;
      while (diff_x(1) < -M_PI) diff_x(1)=+2*M_PI;
      while (diff_x(1) > M_PI) diff_x(1)=-2*M_PI;
      while (diff_z(1) < -M_PI) diff_z(1)=+2*M_PI;
      while (diff_z(1) > M_PI) diff_z(1)=-2*M_PI;
      Tc += weights_(i) * diff_x * diff_z.transpose();
  }

  // Kalman gain K.
  MatrixXd K(n_x_, n_z);
  K = Tc * S.inverse();

  // Update state mean and covariance matrix.
  VectorXd z(3);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];
  VectorXd z_diff = z-z_pred;
  while (z_diff(1) < -M_PI) z_diff(1)=+2*M_PI;
  while (z_diff(1) > M_PI) z_diff(1)=-2*M_PI;

  x_ = x_ + K*z_diff;

  P_ = P_ - K*S*K.transpose();
}