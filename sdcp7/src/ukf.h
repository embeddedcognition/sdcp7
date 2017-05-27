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

        //state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
        VectorXd x_;

        //the current NIS for radar
        double NIS_radar_;

        //the current NIS for laser
        double NIS_laser_;

        //process the latest measurement received from the sensor
        void ProcessMeasurement(const MeasurementPackage& measurement_pack);

    private:
        ///* initially set to false, set to true in first call of ProcessMeasurement
        bool is_initialized_;

        ///* if this is false, laser measurements will be ignored (except for init)
        bool use_laser_;

        ///* if this is false, radar measurements will be ignored (except for init)
        bool use_radar_;

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

        //perform unscented kalman prediciton steps
        void PerformPrediction(double delta_t);

        /**
         * Updates the state and the state covariance matrix using a laser measurement
         * @param meas_package The measurement at k+1
         */
        void UpdateLidar(const MeasurementPackage& measurement_pack);

        /**
         * Updates the state and the state covariance matrix using a radar measurement
         * @param meas_package The measurement at k+1
         */
        void UpdateRadar(const MeasurementPackage& measurement_pack);

        //normalize the supplied angle to be within -pi to pi
        double NormalizeAngle(const double angle);

        //function that handles first time init
        void FirstTimeInit(const MeasurementPackage& measurement_pack);

        //generate the augmented sigma points
        void ComputeAugmentedSigmaPoints(MatrixXd* Xsig_out);

        //predict sigma points by pushing the augmented state through the CTRV model
        void PredictSigmaPoints(MatrixXd& Xsig_aug_in, const double delta_t_in);

        //last step in delivering the new predicted state and covariance, we need to compute the mean and covariance of the predicted state
        void ComputeMeanAndCovarianceofPredictedSigmaPoints();

        //transform the predicted state into the radar measurement space and then predict measurement mean and covariance
        void PredictRadarMeasurement(const int& n_z_in, VectorXd* z_pred_out, MatrixXd* S_out, MatrixXd* Zsig_out);

        //transform the predicted state into the lidar measurement space and then predict measurement mean and covariance
        void PredictLidarMeasurement(const int& n_z_in, VectorXd* z_pred_out, MatrixXd* S_out, MatrixXd* Zsig_out);

        //final step for this iteration, update the predicted state and covariance using the radar/lidar measurement
        void UpdateState(const VectorXd& z_in, const VectorXd& z_pred_in, const MatrixXd& Zsig_in, const MatrixXd& S_in, const int& n_z_in);
};

#endif /* UKF_H */
