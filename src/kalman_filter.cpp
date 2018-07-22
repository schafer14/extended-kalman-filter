#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  MatrixXd I = MatrixXd::Identity(4, 4);
	MatrixXd Ht = H_.transpose();

  MatrixXd y = z - H_ * x_;
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();

  x_ = x_ + (K * y);
  P_ = (I - K * H_) * P_;

  
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];

  MatrixXd Hj = tools.CalculateJacobian(x_);

  MatrixXd I = MatrixXd::Identity(4, 4);

  MatrixXd R = MatrixXd(3, 3);
  R << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  if (z[1] < -3.14159 || z[1] > 3.14159) {
    cout << "Bad radar measurement: " << z[1] << endl;
    return ;
  }

  VectorXd hX = VectorXd(3);
  hX << sqrt(px * px + py * py),
        atan2(py, px),
        (px * vx + py * vy) / sqrt(px * px + py * py);
  


  MatrixXd y = z - hX;

	while (y(1) > M_PI) {
		y(1) -= 2 * M_PI;
	}

	while (y(1) < -M_PI) {
		y(1) += 2 * M_PI;
	}

	MatrixXd Ht = Hj.transpose();
  MatrixXd S = Hj * P_ * Ht + R;
  MatrixXd K = P_ * Ht * S.inverse();


  x_ = x_ + (K * y);
  P_ = (I - K * Hj) * P_;
}
