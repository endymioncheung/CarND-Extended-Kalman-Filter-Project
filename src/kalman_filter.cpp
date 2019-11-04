#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in; // object state
  P_ = P_in; // object covariance matrix
  F_ = F_in; // state transition matrix
  H_ = H_in; // measurement matrix
  R_ = R_in; // measurement covariance matrix
  Q_ = Q_in; // process covariance matrix
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  // KF Prediction Step
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  VectorXd z_pred = H_ * x_;
  VectorXd y      = z - z_pred;
  
  MatrixXd Ht = H_.transpose();
  MatrixXd S  = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K  = P_ * Ht * Si;

  // new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

double normalize_angle(double x){
  // Normalize angles to [-PI,PI]
  x = fmod(x+M_PI, 2*M_PI);
  if (x < 0)
      x += 2*M_PI;
  return x - M_PI;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
   
  double px = x_(0);
  double py = x_(1);
  double vx = x_(2);
  double vy = x_(3);

  // Convert cartesian coordinates to polar coordinates
  double rho     = sqrt(px*px + py*py);
  double phi     = atan2(py, px);
  double rho_dot = 0;

  // check division by zero (or close to zero)
  if (fabs(rho) < 0.0001) {
    rho_dot = 0;
    // cout << "UpdateEKF() - Error - Division by Zero" << endl;
  } else {
    rho_dot = (px*vx + py*vy) / rho;
  }
  
  VectorXd z_pred(3);
  z_pred << rho, phi, rho_dot;
  VectorXd y = z - z_pred;
  
  /* Normalizing the angles because the Kalman filter
   is expecting small angle values
   between the range -pi and pi.
  */
  y(1) = normalize_angle(y(1));

  MatrixXd Ht = H_.transpose();
  MatrixXd S  = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K  = P_ * Ht * Si;
  
  // new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
