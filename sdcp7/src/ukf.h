/*
#######################################################
## AUTHOR: James Beasley                             ##
## DATE: May 8, 2017                                 ##
## UDACITY SDC: Project 7 (Unscented Kalman Filters) ##
#######################################################
*/

#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF
{
    public:
        //constructor
        UKF();
        //destructor
        virtual ~UKF();

        //process the latest measurement received from the sensor
        void ProcessMeasurement(const MeasurementPackage& measurement_pack);

        //return the current state vector x --> (px, py, v, yaw, yaw_dot)
        VectorXd GetState();

        //return NIS for radar
        double GetRadarNIS();

        //return NIS for laser
        double GetLaserNIS();

    private:
        ///* initially set to false, set to true in first call of ProcessMeasurement
        bool is_initialized_;

        //the current NIS for radar
        double NIS_radar_;

        //the current NIS for laser
        double NIS_laser_;

        //state vector: px (position x), py (position y), v (speed or magnitude of velocity), yaw (orientation), yaw_dot (rate of change of orientation, if this is zero the vehicle is traveling in a straight line)
        //in the CTRV model the velocity and yaw rate are constant
        VectorXd x_;

        //state covariance matrix
        MatrixXd P_;

        //predicted sigma points matrix
        MatrixXd Xsig_pred_;

        //weights of sigma points
        VectorXd weights_;

        //previous timestamp
        long long previous_timestamp_;

        //process noise standard deviation longitudinal acceleration in m/s^2
        const double std_a_;

        //process noise standard deviation yaw acceleration in rad/s^2
        const double std_yawdd_;

        //laser measurement noise standard deviation position1 in m
        const double std_laspx_;

        //laser measurement noise standard deviation position2 in m
        const double std_laspy_;

        //radar measurement noise standard deviation radius in m
        const double std_radr_;

        //radar measurement noise standard deviation angle in rad
        const double std_radphi_;

        //radar measurement noise standard deviation radius change in m/s
        const double std_radrd_ ;

        //state dimension
        const int n_x_;

        //augmented state dimension
        const int n_aug_;

        //sigma point spreading parameter
        const double lambda_;

        //PI
        const double PI;

        //function that handles first time init
        void FirstTimeInit(const MeasurementPackage& measurement_pack);

        //perform unscented kalman prediciton steps
        void PerformPrediction(const double delta_t);

        //update state and covariance using radar measurement
        void UpdateRadar(const MeasurementPackage& measurement_pack);

        //update state and covariance using lidar measurement
        void UpdateLidar(const MeasurementPackage& measurement_pack);

        //generate the augmented sigma points
        void ComputeAugmentedSigmaPoints(MatrixXd& Xsig_aug);

        //predict sigma points by pushing the augmented state through the CTRV model
        void PredictSigmaPoints(const MatrixXd& Xsig_aug, const double delta_t);

        //last step in delivering the new predicted state and covariance, we need to compute the mean and covariance of the predicted state
        void ComputeMeanAndCovarianceofPredictedSigmaPoints();

        //normalize the supplied angle to be within -pi to pi
        double NormalizeAngle(const double angle);
};

#endif /* UKF_H */
