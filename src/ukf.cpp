#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF()
{
    
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;
    
    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 0.50; // was 30 and we are talking about cycle

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.50; // was 30

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
     TODO: Done

     Complete the initialization. See ukf.h for other member properties.

     Hint: one or more values initialized above might be wildly off...
     */
    
    is_initialized_ = false; //

    // Set state dimension
    n_x_ = 5;
    
    // Define spreading parameter
    lambda_ = 3 - n_x_;
    
    // Augmented state dimension
    n_aug_ = n_x_ + 2;
    
    // Predicted Sigma point matrix
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    
    // State vector
    x_ = VectorXd(n_x_);
    
    // Initial covariance matrix.. start with an identity matrix..
    P_ = MatrixXd(n_x_, n_x_);
    P_ <<   1,  0,  0,  0,  0,
            0,  1,  0,  0,  0,
            0,  0,  1,  0,  0,
            0,  0,  0,  1,  0,
            0,  0,  0,  0,  1;
    
    // Set weights for sigma points
    weights_ = VectorXd(2 * n_aug_ + 1);
    weights_(0) = lambda_/(lambda_+n_aug_);
    
    for (int i=1; i< 2 * n_aug_+ 1; i++)
    {
        weights_(i) = 0.5/(lambda_ + n_aug_);
    }
    
}

UKF::~UKF() {}


/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
  /**
  TODO: Done

  Complete this function! Make sure you switch between laser and radar
  measurements.
  */
    if( !is_initialized_ )
    {
        time_us_ = meas_package.timestamp_;
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
        {
            float px    = meas_package.raw_measurements_(0) * cos(meas_package.raw_measurements_(1));
            float py    = meas_package.raw_measurements_(0) * sin(meas_package.raw_measurements_(1));
            float v     = meas_package.raw_measurements_(2);
            float yaw   = meas_package.raw_measurements_(1);
            x_ << px, py, 0, 0, 0;
            is_initialized_ = true;
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
        {
            float px    = meas_package.raw_measurements_(0);
            float py    = meas_package.raw_measurements_(1);
            x_ << px, py, 0, 0, 0;
            is_initialized_ = true;
        }
        return;
    }
    
    // Compute the time elapsed between the current and previous measurements
    double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
    
    // Prediction step ...
    Prediction(dt);
    
    // Update step ...
    if ((meas_package.sensor_type_ == MeasurementPackage::RADAR) && use_radar_)
    {
        UpdateRadar(meas_package);
    }
    else if ((meas_package.sensor_type_ == MeasurementPackage::LASER) && use_laser_ )
    {
        UpdateLidar(meas_package);
    }
    time_us_ = meas_package.timestamp_;
    
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t)
{
  /**
  TODO: Done

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  /* Prediction involves the following steps:
   
   - Creating the augmented mean state (x_aug_)
   - Building the covariance matrix (P_)
   - Generating the augmented sigma points (XSig_aug_)
   - Predicting the sigma points (Xsig_pred)
   - Predicting the state mean (X_)
   - Predicting the covariance matrix (P_)
   */
    
    //create augmented mean state
    VectorXd x_aug_     = VectorXd(n_aug_);
    x_aug_.head(n_x_)   = x_;
    x_aug_(5) = 0;
    x_aug_(6) = 0;
    
    /* Augmented Covariance Matrix */
    
    MatrixXd P_aug_     = MatrixXd(n_aug_, n_aug_);
    P_aug_.fill(0.0);
    P_aug_.topLeftCorner(n_x_, n_x_) = P_;
    
    // process noise
    
    P_aug_(5,5) = std_a_ * std_a_;
    P_aug_(6,6) = std_yawdd_ * std_yawdd_; // Process Noise
    
    //create square root matrix
    MatrixXd L = P_aug_.llt().matrixL();
    
    //create augmented sigma points
    MatrixXd Xsig_aug_  = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    Xsig_aug_.col(0)    = x_aug_;
    double sqrt_lambda_n_aug = sqrt(lambda_ + n_aug_);
    for (int i = 0; i< n_aug_; i++)
    {
        Xsig_aug_.col(i+1)          = x_aug_ + sqrt_lambda_n_aug * L.col(i);
        Xsig_aug_.col(i+1+n_aug_)   = x_aug_ - sqrt_lambda_n_aug * L.col(i);
    }
    
    //predict sigma points
    for (int i = 0; i< 2 * n_aug_ + 1; i++)
    {
        //extract values for better readability
        const double px           = Xsig_aug_(0,i);
        const double py           = Xsig_aug_(1,i);
        const double v            = Xsig_aug_(2,i);
        const double yaw          = Xsig_aug_(3,i);
        const double yawd         = Xsig_aug_(4,i);
        const double nu_a         = Xsig_aug_(5,i);
        const double nu_yawdd     = Xsig_aug_(6,i);
        
        //predicted state values
        double px_pred, py_pred;
        double v_pred       = v;
        double yaw_pred     = yaw + yawd * delta_t;
        double yawd_pred    = yawd;
        
        
        //avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_pred = px + (v/yawd) * (sin(yaw_pred) - sin(yaw));
            py_pred = py + (v/yawd) * (cos(yaw) - cos(yaw_pred));
        }
        else {
            px_pred = px + v * delta_t * cos(yaw);
            py_pred = py + v * delta_t * sin(yaw);
        }
        
        
        //add noise elements
        px_pred     += 0.5 * nu_a * delta_t * delta_t * cos(yaw);
        py_pred     += 0.5 * nu_a * delta_t * delta_t * sin(yaw);
        v_pred      += nu_a * delta_t;
        yaw_pred    += 0.5 * nu_yawdd * delta_t * delta_t;
        yawd_pred   += nu_yawdd * delta_t;
        
        //write predicted sigma point into right column
        Xsig_pred_(0,i) = px_pred;
        Xsig_pred_(1,i) = py_pred;
        Xsig_pred_(2,i) = v_pred;
        Xsig_pred_(3,i) = yaw_pred;
        Xsig_pred_(4,i) = yawd_pred;
    }
    
    
    //predicted state mean
    x_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++)
    {  //iterate over sigma points
        x_ = x_ + weights_(i) * Xsig_pred_.col(i);
    }
    
    //predicted state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) //iterate over sigma points
    {
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        
        //angle normalization
        x_diff(3) = tools_.NormalizeAngle(x_diff(3));
        P_ += weights_(i) * x_diff * x_diff.transpose() ;
    }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package)
{
  /**
  TODO: Done

  Complete this function! Use laser data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the laser NIS.
  */
    int n_z = 2; // Number of Lidar measurements - only x and y
    
    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    
    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++)
    {
        // extract values for better readibility
        double px_pred = Xsig_pred_(0,i);
        double py_pred = Xsig_pred_(1,i);
        
        // measurement model
        Zsig(0,i) = px_pred ;
        Zsig(1,i) = py_pred ;
    }
    
    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    
    for (int i=0; i < 2 * n_aug_ + 1; i++)
    {
        z_pred += weights_(i) * Zsig.col(i);
    }
   
    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++)
    {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        S += weights_(i) * z_diff * z_diff.transpose();
    }
    
    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z, n_z);
    R <<    std_laspx_ * std_laspx_, 0,
            0,                       std_laspy_ * std_laspy_;
    
    S  += R;
    
    // Calclulate cross correlation matrix
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);
    
    for (int i = 0; i < 2 * n_aug_ + 1; i++)
    {
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        VectorXd z_diff = Zsig.col(i) - z_pred;
        
        // Normalize the angle to be within -PI to +PI
        x_diff(3) = tools_.NormalizeAngle(x_diff(3));
        Tc += weights_(i) * x_diff * z_diff.transpose();
    }
    
    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    
    // Update state and covariance matrices
    VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
    
    x_  = x_ + K * z_diff;
    P_ -= K * S * K.transpose();
    
    // laser NIS ...
    NIS_laser_  = z_diff.transpose() * S.inverse() * z_diff;

}


/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package)
{
  /**
  TODO: Done

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
    
  int n_z = 3;
  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);
    
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {  //2n+1 simga points
        
    // extract values for better readibility
    double px_pred      = Xsig_pred_(0,i);
    double py_pred      = Xsig_pred_(1,i);
    double v_pred       = Xsig_pred_(2,i);
    double yaw_pred     = Xsig_pred_(3,i);
    double yaw_dot_pred = Xsig_pred_(4,i);
    double p_pred_meas  = sqrt(px_pred * px_pred + py_pred * py_pred);
    
    if (p_pred_meas < 0.00001) p_pred_meas = 0.00001 ; // Ensure no divide by 0 error
    
    // measurement model
      
      Zsig(0,i) = p_pred_meas ;                                                                             //r

    // Check for divide by 0 for phi..
      
    if (px_pred != 0)
    {
        Zsig(1,i) = atan2(py_pred, px_pred);                                                                //phi
        Zsig(2,i) = (px_pred * cos(yaw_pred) * v_pred + py_pred * sin(yaw_pred) * v_pred) / p_pred_meas ;   //r_dot
    }
      
    else
    {
        Zsig(1,i) = 0;
        Zsig(2,i) = 0;
    }
      
  }
    
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
    
  for (int i=0; i < 2 * n_aug_ + 1; i++)
  {
    z_pred += weights_(i) * Zsig.col(i);
  }
    
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
    
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {  //2n+1 simga points
    VectorXd z_diff = Zsig.col(i) - z_pred;
        
    //angle normalization
    z_diff(1) = tools_.NormalizeAngle(z_diff(1));
    S   += weights_(i) * z_diff * z_diff.transpose();
  }
    
  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_ * std_radr_,    0,                          0,
            0,                      std_radphi_ * std_radphi_,  0,
            0,                      0,                          std_radrd_ * std_radrd_;
    
  S += R;
    
  //calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
    
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {  //2n+1 simga points
        
    VectorXd z_diff = Zsig.col(i) - z_pred;
        
    //angle normalization
    z_diff(1) = tools_.NormalizeAngle(z_diff(1));
      
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
        
    //angle normalization
    x_diff(3) = tools_.NormalizeAngle(x_diff(3));
      
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
    
  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();
    
  // incoming radar measurement
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
    
  //angle normalization
  z_diff(1) = tools_.NormalizeAngle(z_diff(1));
    
  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();
    
  // Radar NIS ...
  NIS_radar_  = z_diff.transpose() * S.inverse() * z_diff;
    
 }
