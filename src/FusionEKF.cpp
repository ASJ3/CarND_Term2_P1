#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>//A: adding to be able to do cos and sin operations

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;//A:

  F_ = MatrixXd(4,4);//A: Designing a generic F matrix for the predict function
  F_ << 1,0,1,0,
        0,1,0,1,
        0,0,1,0,
        0,0,0,1;//A:
  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  KalmanFilter ekf_; //A:

  
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    previous_timestamp_ = measurement_pack.timestamp_;//A:

    
    ekf_.P_ = MatrixXd(4,4);//A:
    

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      ekf_.x_[0] = measurement_pack.raw_measurements_[0]*cos(measurement_pack.raw_measurements_[1]);//A: x
      ekf_.x_[1] = measurement_pack.raw_measurements_[0]*sin(measurement_pack.raw_measurements_[1]);//A: y
      ekf_.x_[2] = measurement_pack.raw_measurements_[2]*cos(measurement_pack.raw_measurements_[1]);//A: vx
      ekf_.x_[3] = measurement_pack.raw_measurements_[2]*sin(measurement_pack.raw_measurements_[1]);//A: vy

      //A: Initial covariance matrix with relatively high uncertainty with regards to speed
      ekf_.P_ << 1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 10, 0,
                0, 0, 0, 10;//A:

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;//A:

      //A: Initial covariance matrix with high uncertainty with regards to speed
      ekf_.P_ << 1, 0, 0, 0,
                 0, 1, 0, 0,
                 0, 0, 1000, 0,
                 0, 0, 0, 1000;//A:

    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  float timeInterval = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;;//A:
  float timeInterval_sq = timeInterval * timeInterval;//A:
  float timeInterval_cube = (timeInterval_sq * timeInterval)/2.0;//A:
  float timeInterval_quad = (timeInterval_cube * timeInterval)/4.0;//A:

  float noise_ax = 9.0;
  float noise_ay = 9.0;

  previous_timestamp_ = measurement_pack.timestamp_;//A:

  //A: Defining Q, the covariance matrix of the process noise nu
  MatrixXd Q = MatrixXd(4,4);//A:
  Q << timeInterval_quad*noise_ax, 0, timeInterval_cube*noise_ax, 0,
       0, timeInterval_quad*noise_ay, 0, timeInterval_cube*noise_ay,
       timeInterval_cube*noise_ax, 0, timeInterval_sq*noise_ax, 0,
       0, timeInterval_cube*noise_ay, 0, timeInterval_sq*noise_ay;//A:

  //A: Now that we have the time interval, updating F_
  F_(0,2) = timeInterval;//A:
  F_(1,3) = timeInterval;//A: 

  ekf_.x_ = F_*ekf_.x_;//A:
  ekf_.P_ = F_*ekf_.P_*F_.transpose() + Q;//A: result in a 4x4 matrix
  //ekf_.Predict();//:A Inactivating Predict() for now

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
      VectorXd z;//A:
      z = VectorXd(3);//A:
      z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], measurement_pack.raw_measurements_[2];//A:
      float sqr_add = ekf_.x_[0]*ekf_.x_[0] + ekf_.x_[1]*ekf_.x_[1];//A:
      float h_first = sqrt(sqr_add);//A:
      float denominator_3_2 = sqr_add*h_first;//A:
      float h_sec = atan2(ekf_.x_[1], ekf_.x_[0]);//A:
      float h_third = (ekf_.x_[0]*ekf_.x_[2] + ekf_.x_[1]*ekf_.x_[3])/h_first;//A:

      Hj_ << ekf_.x_[0]/h_first, ekf_.x_[1]/h_first, 0, 0,
             -ekf_.x_[1]/sqr_add, ekf_.x_[0]/sqr_add,0, 0,
             ekf_.x_[1]*(ekf_.x_[2]*ekf_.x_[1] - ekf_.x_[3]*ekf_.x_[0])/denominator_3_2,   ekf_.x_[0]*(ekf_.x_[3]*ekf_.x_[0] - ekf_.x_[2]*ekf_.x_[1])/denominator_3_2, ekf_.x_[0]/h_first, ekf_.x_[1]/h_first;//A:


      //A:needs to make sure phi is between -pi and pi by adding or subtracting 2pi
      float Pi = 3.14159265358979323846;
      while (h_sec < -Pi) {
        h_sec = h_sec + 2*Pi;
      }

      while (h_sec > Pi) {
        h_sec = h_sec - 2*Pi;
      }

      VectorXd hx;//A:
      hx = VectorXd(3);//A:
      hx << h_first, h_sec, h_third;//A:

      VectorXd y_polar = z - hx;//A:
      //A:needs to convert y_polar back to cartesian coordinates.
      /*VectorXd y_cartesian;
      y_cartesian = VectorXd(4);
      y_cartesian[0] = y_polar[0]*cos(y_polar[1]);//A: x
      y_cartesian[1] = y_polar[0]*sin(y_polar[1]);//A: y
      y_cartesian[2] = y_polar[2]*cos(y_polar[1]);//A: vx
      y_cartesian[3] = y_polar[2]*sin(my_polar[1]);//A: vy
      */
      

      MatrixXd S = Hj_ * ekf_.P_ * Hj_.transpose() + R_radar_;//A: result: 3x3 matrix
      MatrixXd K = ekf_.P_ * Hj_.transpose()*S.inverse();//A: result: 4x3 matrix

      ekf_.x_ = ekf_.x_ + K * y_polar;//A: result: vector of 4
      MatrixXd I = MatrixXd::Identity(4, 4);//A:
      ekf_.P_ = (I - K*Hj_) * ekf_.P_;//A: result: 4x4 matrix

  } else {
    // Laser updates
      VectorXd z;//A:
      z = VectorXd(2);//A:
      z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1];//A:

      VectorXd y = z - H_laser_*ekf_.x_;//A:
      MatrixXd Ht = H_laser_.transpose();//A: result: 4x2 matrix
      MatrixXd S = H_laser_ * ekf_.P_ * Ht + R_laser_;//A: result: 2x2 matrix
      MatrixXd K = ekf_.P_ * Ht * S.inverse();//A: result: 4x2 matrix

      ekf_.x_ = ekf_.x_ + K * y;//A: K is a matrix with 4 rows, 2 cols and y is a vector with two rows -> we get a 4X1 matrix

      MatrixXd I = MatrixXd::Identity(4, 4);//A:
      ekf_.P_ = (I - K*H_laser_) * ekf_.P_;//A: result: 4x4 matrix
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
