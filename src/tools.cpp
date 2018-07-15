#include <iostream>
#include <stdexcept>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */

  VectorXd result = VectorXd(5);
  result << 0, 0, 0, 0, 0;
  return result;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

  if (px == 0 && py == 0){
	  throw std::overflow_error("Divide by zero exception");
	}
	
	float dpdpx = px / sqrt(px * px + py * py);
	float dpdpy = py / sqrt(px * px + py * py);
	float drdpx = - py / (px * px + py * py);
	float drdpy = px / (px * px + py * py);
	float dvdpx = py * (vx * py - vx * py) / pow(px * px + py * py, 3 / 2);
	float dvdpy = px * (vy * px - vx * py) / pow(px * px + py * py, 3 / 2);
	float dvdvx = px / sqrt(px * px + py * py);
	float dvdvy = py / sqrt(px * px + py * py);
	
	Hj << dpdpx, dpdpy, 0, 0, 
	      drdpx, drdpy, 0, 0,
	      dvdpx, dvdpy, dvdvx, dvdvy;

	return Hj;
}
