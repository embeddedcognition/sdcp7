/*
#######################################################
## AUTHOR: James Beasley                             ##
## DATE: May 8, 2017                                 ##
## UDACITY SDC: Project 7 (Unscented Kalman Filters) ##
#######################################################
*/

#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

//constructor
UKF::UKF()
{
    //bool to invoke initialization on first "ProcessMeasurement" call
    is_initialized_ = false;

    //if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    //if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    //initial state vector
    x_ = VectorXd(5);

    // initial covariance matrix
    P_ = MatrixXd(5, 5);

    //process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 30;

    //process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 30;

    //laser measurement noise standard deviation position1 in m
    std_laspx_ = 0.15;

    //laser measurement noise standard deviation position2 in m
    std_laspy_ = 0.15;

    //radar measurement noise standard deviation radius in m
    std_radr_ = 0.3;

    //radar measurement noise standard deviation angle in rad
    std_radphi_ = 0.03;

    //radar measurement noise standard deviation radius change in m/s
    std_radrd_ = 0.3;
}

//destructor
UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(const MeasurementPackage& measurement_pack)
{
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t)
{
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const MeasurementPackage& measurement_pack)
{
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const MeasurementPackage& measurement_pack)
{
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}

//normalize the supplied angle to be within -pi to pi
double UKF::NormalizeAngle(const double angle)
{
    //local vars
    double normalized_angle = angle;

    //adjust phi to be between -pi and pi
    //http://stackoverflow.com/questions/11980292/how-to-wrap-around-a-range
    if (fabs(angle) > PI)
    {
        double two_pi = 2 * PI;
        normalized_angle -= round(normalized_angle / two_pi) * two_pi;
    }

    return normalized_angle;
}
