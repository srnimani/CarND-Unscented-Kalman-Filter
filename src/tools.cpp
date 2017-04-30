#include <iostream>
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
    VectorXd rmse(4);
    rmse << 0,0,0,0;
    int n = estimations.size();
    for (int i=0; i<n; i++) {
        VectorXd residual_error = estimations[i] - ground_truth[i];
        residual_error = (residual_error.array()*residual_error.array());
        rmse += residual_error;
    }
    rmse = rmse / n;
    rmse = rmse.array().sqrt();
    return rmse;
}

double Tools::NormalizeAngle(double angle)
{
    while (angle > M_PI) angle -= 2. * M_PI;
    while (angle <-M_PI) angle += 2. * M_PI;
    return angle ;
}
