/*
#######################################################
## AUTHOR: James Beasley                             ##
## DATE: May 8, 2017                                 ##
## UDACITY SDC: Project 7 (Unscented Kalman Filters) ##
#######################################################
*/

#include "ukf.h"
#include <iostream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

//constructor
UKF::UKF()
//initialize constants through member initialization list
: std_a_(3), std_yawdd_(2), std_laspx_(0.15), std_laspy_(0.15), std_radr_(0.3),
  std_radphi_(0.03), std_radrd_(0.3), n_x_(5), n_aug_(7), lambda_(3 - n_aug_), PI(3.14159265358979)
{
    //bool to invoke initialization on first "ProcessMeasurement" call
    is_initialized_ = false;

    //initial state vector
    x_ = VectorXd(n_x_);

    //initial covariance matrix
    P_ = MatrixXd(n_x_, n_x_);

    //predicted sigma points
    Xsig_pred_ = MatrixXd(n_x_, (2 * n_aug_) + 1);

    //weights for sigma points
    weights_ = VectorXd((2 * n_aug_) + 1);
    weights_.fill(1 / (2 * (lambda_ + n_aug_)));  //set all elements
    weights_(0) = lambda_ / (lambda_ + n_aug_);   //set first element differently

    //previous measurement timestamp
    previous_timestamp_ = 0;

    //NIS
    NIS_radar_ = 0;
    NIS_laser_ = 0;
}

//destructor
UKF::~UKF() {}

//process the latest measurement received from the sensor
void UKF::ProcessMeasurement(const MeasurementPackage& measurement_pack)
{
    //local vars
    double delta_t;    //elapsed time between the current and previous measurements (in seconds)

    //if this is the first time, we need to initialize x with the current measurement,
    //then exit (filtering will start with the next iteration)
    if (!is_initialized_)
    {
        FirstTimeInit(measurement_pack);
        //no need to predict or update so we return
        return;
    }

    //compute the time elapsed between the current and previous measurements (in seconds)
    delta_t = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;

    //capture the timestamp of the measurement for the next iteration
    previous_timestamp_ = measurement_pack.timestamp_;

    //perform kalman prediction step
    PerformPrediction(delta_t);

    //perform kalman update step
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
        //RADAR
        UpdateRadar(measurement_pack);
    }
    else
    {
        //LIDAR
        UpdateLidar(measurement_pack);
    }

    // print the output
    cout << "x_ = " << x_ << endl;
    cout << "P_ = " << P_ << endl;
}

//this is executed the first time only
void UKF::FirstTimeInit(const MeasurementPackage& measurement_pack)
{
    //local vars
    //for laser: col_1 will equal px and col_2 will equal py
    //for radar: col_1 will equal r (e.g., rho) and col_2 will equal θ (e.g., phi)
    double col_1 = measurement_pack.raw_measurements_(0);
    double col_2 = measurement_pack.raw_measurements_(1);
    double px; //x position
    double py; //y position

    //if this is radar data
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
        //we need to convert from polar coordinates (r, θ) to cartesian coordinates (x, y), via x = r * cos(θ) and y = r * sin(θ)
        px = col_1 * cos(col_2);
        py = col_1 * sin(col_2);
    }
    //if this is lidar data
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
    {
        //we're already in cartesian coordinates, so just assign to px and py
        px = col_1;
        py = col_2;
    }

    //init state vector x, which is in cartesian coordinates (px, py, v, psi, psi_dot)
    //first measurement, so position is the only valuable information
    x_ << px, py, 0, 0, 0;

    //init the matrix P to the identity matrix
    P_ = MatrixXd::Identity(n_x_, n_x_);

    //capture the timestamp of the measurement for future use
    previous_timestamp_ = measurement_pack.timestamp_;

    //done initializing
    is_initialized_ = true;
}

//perform kalman prediction step
//at this point we have the posterior mean state and covariance matrix from the last iteration (representing distribution of current state)
//the model we're using (CTRV) has a state vector of size 5
void UKF::PerformPrediction(const double delta_t)
{
    //local vars
    MatrixXd Xsig_aug = MatrixXd(n_aug_, (2 * n_aug_) + 1);     //augmented sigma point matrix (state + process noise)

    //represent the uncertainty of the posterior state estimation with sigma points
    //the augmented state includes the noise vector (last two elements)
    //first sigma point is always the mean state estimate
    ComputeAugmentedSigmaPoints(Xsig_aug);

    //predict sigma points (by inserting each augmented sigma point into the CTRV process model)
    PredictSigmaPoints(Xsig_aug, delta_t);

    //compute the mean and covariance of the predicted state
    ComputeMeanAndCovarianceofPredictedSigmaPoints();
}

//update state and covariance using radar measurement
void UKF::UpdateRadar(const MeasurementPackage& measurement_pack)
{
    //local vars
    const int n_z = 3;                                    //set measurement dimension, radar can measure r, phi, and r_dot
    MatrixXd Zsig = MatrixXd(n_z, (2 * n_aug_) + 1);      //predicted sigma points mapped into the radar measurement space
    VectorXd z_pred = VectorXd(n_z);                      //mean predicted measurement
    MatrixXd S = MatrixXd(n_z, n_z);                      //measurement noise covariance
    MatrixXd R = MatrixXd(n_z, n_z);                      //measurement noise matrix
    MatrixXd Tc = MatrixXd(n_x_, n_z);                    //create cross correlation matrix
    double px, py, v, yaw, rho, phi, rho_dot;             //used to extract particular values for easier handling
    double normalized_angle; //used to normalize angles

    //init since we're summing against them
    S.fill(0);
    R.fill(0);
    Tc.fill(0);

    //transform sigma points into measurement space
    //loop through the predicted sigma points and convert them to the measurement space
    for (int i = 0; i < (2 * n_aug_) + 1; i++)
    {
        //extract values from cur predicted sigma point column
        px = Xsig_pred_(0, i);
        py = Xsig_pred_(1, i);
        v = Xsig_pred_(2, i);
        yaw = Xsig_pred_(3, i);

        //compute the rho transformation
        rho = sqrt((px * px) + (py * py));

        //check for divide by zero for rho dot computation
        if (rho < 0.0001)
        {
            rho = 0.0001;
        }

        //check for divide by zero for phi
        if (fabs(px) < 0.0001)
        {
            //compute the phi transformation
            phi = atan2(py, 0.0001);
        }
        else
        {
            //compute the phi transformation
            phi = atan2(py, px);
        }

        //compute the rho_dot transformation
        rho_dot = ((px * cos(yaw) * v) + (py * sin(yaw) * v)) / rho;

        //assign the computed factors to the current column of the measurement space matrix
        Zsig.col(i) << rho, phi, rho_dot;
    }

    //compute mean predicted measurement
    //multiply the weights with each element in the measurement sigma points matrix
    //then do a row-wise sum of the resulting matrix so that a summed vector is returned
    z_pred << (Zsig * weights_).rowwise().sum();

    //calculate measurement covariance matrix S
    //save time and space by computing the mean and transpose ahead of time
    MatrixXd mean_subtracted = Zsig.colwise() - z_pred;
    MatrixXd mean_subtracted_transposed = mean_subtracted.transpose();
    //sum the individual matrices
    for (int i = 0; i < (2 * n_aug_) + 1; i++)
    {
        //ensure phi is between -pi and pi
        normalized_angle = NormalizeAngle(mean_subtracted(1, i));
        mean_subtracted(1, i) = normalized_angle;
        mean_subtracted_transposed(i, 1) = normalized_angle;
        S += (weights_(i) * mean_subtracted.col(i) * mean_subtracted_transposed.row(i));
    }

    //compute the measurement noise matrix
    R(0, 0) = std_radr_ * std_radr_;
    R(1, 1) = std_radphi_ * std_radphi_;
    R(2, 2) = std_radrd_ * std_radrd_;

    //add the noise matrix to the covariance matrix
    S += R;

    //compute cross correlation matrix Tc
    //compute means early, saving time and space
    MatrixXd Xmean_subtracted = Xsig_pred_.colwise() - x_;
    MatrixXd Zmean_subtracted_transposed = (Zsig.colwise() - z_pred).transpose();
    //sum the individual matrices
    for (int i = 0; i < (2 * n_aug_) + 1; i++)
    {
        //ensure yaw is between -pi and pi
        Xmean_subtracted(3, i) = NormalizeAngle(Xmean_subtracted(3, i));
        //ensure phi is between -pi and pi
        Zmean_subtracted_transposed(i, 1) = NormalizeAngle(Zmean_subtracted_transposed(i, 1));
        Tc += (weights_(i) * Xmean_subtracted.col(i) * Zmean_subtracted_transposed.row(i));
    }

    //compute Kalman gain K;
    MatrixXd K = Tc * S.inverse();

    //compute error
    VectorXd z = measurement_pack.raw_measurements_;
    VectorXd z_diff = z - z_pred;
    //ensure phi is between -pi and pi
    z_diff(1) = NormalizeAngle(z_diff(1));

    //compute NIS
    NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

    //update state mean and covariance matrix
    x_ += K * z_diff;
    P_ -= K * S * K.transpose();
}

//update state and covariance using lidar measurement
void UKF::UpdateLidar(const MeasurementPackage& measurement_pack)
{
    //local vars
    const int n_z = 2;                                    //set measurement dimension, radar can measure px and py
    MatrixXd Zsig = MatrixXd(n_z, (2 * n_aug_) + 1);      //predicted sigma points mapped into the radar measurement space
    VectorXd z_pred = VectorXd(n_z);                      //mean predicted measurement
    MatrixXd S = MatrixXd(n_z, n_z);                      //measurement noise covariance
    MatrixXd R = MatrixXd(n_z, n_z);                      //measurement noise matrix
    MatrixXd Tc = MatrixXd(n_x_, n_z);                    //create cross correlation matrix
    double px, py;                                        //used to extract particular values for easier handling

    //init since we're summing against them
    S.fill(0);
    R.fill(0);
    Tc.fill(0);

    //transform sigma points into measurement space
    //loop through the predicted sigma points and convert them to the measurement space
    for (int i = 0; i < (2 * n_aug_) + 1; i++)
    {
        //extract values from cur predicted sigma point column
        px = Xsig_pred_(0, i);
        py = Xsig_pred_(1, i);

        //assign the computed factors to the current column of the measurement space matrix
        Zsig.col(i) << px, py;
    }

    //compute mean predicted measurement
    //multiply the weights with each element in the measurement sigma points matrix
    //then do a row-wise sum of the resulting matrix so that a summed vector is returned
    z_pred << (Zsig * weights_).rowwise().sum();

    //calculate measurement covariance matrix S
    //save time and space by computing the mean and transpose ahead of time
    MatrixXd mean_subtracted = Zsig.colwise() - z_pred;
    MatrixXd mean_subtracted_transposed = mean_subtracted.transpose();
    //sum the individual matrices
    for (int i = 0; i < (2 * n_aug_) + 1; i++)
    {
        S += (weights_(i) * mean_subtracted.col(i) * mean_subtracted_transposed.row(i));
    }

    //compute the measurement noise matrix
    R(0, 0) = std_laspx_ * std_laspx_;
    R(1, 1) = std_laspy_ * std_laspy_;

    //add the noise matrix to the covariance matrix
    S += R;

    //compute cross correlation matrix Tc
    //compute means early, saving time and space
    MatrixXd Xmean_subtracted = Xsig_pred_.colwise() - x_;
    MatrixXd Zmean_subtracted_transposed = (Zsig.colwise() - z_pred).transpose();
    //sum the individual matrices
    for (int i = 0; i < (2 * n_aug_) + 1; i++)
    {
        //ensure yaw is between -pi and pi
        Xmean_subtracted(3, i) = NormalizeAngle(Xmean_subtracted(3, i));
        Tc += (weights_(i) * Xmean_subtracted.col(i) * Zmean_subtracted_transposed.row(i));
    }

    //compute Kalman gain K;
    MatrixXd K = Tc * S.inverse();

    //compute error
    VectorXd z = measurement_pack.raw_measurements_;
    VectorXd z_diff = z - z_pred;

    //compute NIS
    NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;

    //update state mean and covariance matrix
    x_ += K * z_diff;
    P_ -= K * S * K.transpose();
}

//create augmented sigma points
//allows us to represent the uncertainty of the covariance matrix Q with sigma points
void UKF::ComputeAugmentedSigmaPoints(MatrixXd& Xsig_aug)
{
    //local vars
    MatrixXd Q = MatrixXd(2, 2);                                //process noise covariance matrix
    VectorXd x_aug = VectorXd(n_aug_);                          //augmented mean vector
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);                  //augmented state covariance

    //compute augmented mean state
    x_aug.head(n_x_) = x_; //copy x into first 5 elements of x_aug
    x_aug(5) = 0; //copy mean value of acceleration noise
    x_aug(6) = 0; //copy mean value of acceleration noise

    //init process noise covariance matrix
    //we're supplied the std dev for longitudinal acceleration and yaw acceleration, but we must square to get variance
    Q << (std_a_ * std_a_), 0,
         0, (std_yawdd_ * std_yawdd_);

    //compute augmented state covariance
    //init
    P_aug.fill(0); //if this isn't here, P_aug will be created with garbage in it that will drive x_, P_, and RMSE to inf/nan
    //set P to the top left corner
    P_aug.topLeftCorner(P_.rows(), P_.cols()) = P_;
    //set Q to the bottom right corner
    P_aug.bottomRightCorner(Q.rows(), Q.cols()) = Q;

    //create square root matrix
    //calculate square root of P
    MatrixXd A = P_aug.llt().matrixL();
    //get sqrt of lambda + n_aug
    double sqrt_lambda_plus_naug = sqrt(lambda_ + n_aug_);
    //compute total sqrt
    MatrixXd sqrt_total = sqrt_lambda_plus_naug * A;

    //populate the augmented mean
    Xsig_aug.col(0) = x_aug;

    //populate the sigma points
    for (int i = 0; i < n_aug_; i++)
    {
        //sigma points to the right of the mean state (directionally)
        Xsig_aug.col(i + 1) = x_aug + sqrt_total.col(i);
        //sigma points to the left of the mean (directionally)
        Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt_total.col(i);
    }
}

//perform k + 1 prediction using augmented sigma points
//process noise has a non-linear effect on the state in the process model, therefore augmented sigma points are required
void UKF::PredictSigmaPoints(const MatrixXd& Xsig_aug, const double delta_t)
{
    //local vars
    VectorXd state_pred = VectorXd(n_x_);                   //predicted state (calculated via the CTRV model)
    VectorXd noise_pred = VectorXd(n_x_);                   //predicted noise (calculated via the CTRV model)
    double delta_t_squared = delta_t * delta_t;             //delta_t squared is a multi-use operation (so we compute it once)
    double px, py, v, yaw, yaw_dot, nu_a, nu_yaw_dot_dot;   //used to hold current sigma point column values
    double predicted_px, predicted_py;                      //used to hold predicted px and py values

    //predict sigma points
    //process each sigma point prediction
    for (int i = 0; i < (2 * n_aug_) + 1; i++)
    {
        //extract augmented state column values for easier handling/tracking
        px = Xsig_aug(0, i);
        py = Xsig_aug(1, i);
        v = Xsig_aug(2, i);
        yaw = Xsig_aug(3, i);
        yaw_dot = Xsig_aug(4, i);
        nu_a = Xsig_aug(5, i);
        nu_yaw_dot_dot = Xsig_aug(6, i);

        //deal with potential of yaw_dot being zero in computing state prediction
        if (fabs(yaw_dot) > 0.0001)
        {
            //compute predicted px and py via CTRV when yaw_dot is greater than zero
            predicted_px = (v / yaw_dot) * (sin(yaw + (yaw_dot * delta_t)) - sin(yaw));
            predicted_py = (v / yaw_dot) * (-cos(yaw + (yaw_dot * delta_t)) + cos(yaw));
        }
        else
        {
            //compute predicted px and py via CTRV when yaw_dot is zero
            predicted_px = (v * cos(yaw) * delta_t);
            predicted_py = (v * sin(yaw) * delta_t);
        }

        //compute the process model transformation for each vector member (state)
        state_pred << predicted_px,
                      predicted_py,
                      0,
                      yaw_dot * delta_t,
                      0;
        //compute the process model transformation for each vector member (noise)
        noise_pred << (delta_t_squared * cos(yaw) * nu_a) / 2,
                      (delta_t_squared * sin(yaw) * nu_a) / 2,
                      delta_t * nu_a,
                      (delta_t_squared * nu_yaw_dot_dot) / 2,
                      delta_t * nu_yaw_dot_dot;

        //add state and noise prediction to the current state values to make up the new state (stored in a global var Xsig_pred_)
        Xsig_pred_.col(i) << px + state_pred(0) + noise_pred(0),
                             py + state_pred(1) + noise_pred(1),
                             v + state_pred(2) + noise_pred(2),
                             yaw + state_pred(3) + noise_pred(3),
                             yaw_dot + state_pred(4) + noise_pred(4);
    }
}

//last step in delivering the new predicted state and covariance, we need to compute the mean and covariance of the predicted state
void UKF::ComputeMeanAndCovarianceofPredictedSigmaPoints()
{
    //local vars
    double normalized_angle; //used to normalize angles

    //predict state mean
    //multiply the weights with each element in the predicted sigma points matrix
    //then do a row-wise sum of the resulting matrix so that a summed vector is returned (stored in a global var x_)
    x_ << (Xsig_pred_ * weights_).rowwise().sum();

    //predict state covariance
    //subtract the x vector from each column of Xsig_pred
    //save time and space by computing the mean and transpose ahead of time
    MatrixXd mean_subtracted = Xsig_pred_.colwise() - x_;
    MatrixXd mean_subtracted_transposed = mean_subtracted.transpose();
    //sum the individual 5x5's, result is stored in a global var P_)
    P_.fill(0); //init as we're summing
    for (int i = 0; i < (2 * n_aug_) + 1; i++)
    {
        //ensure yaw is between -pi and pi
        normalized_angle = NormalizeAngle(mean_subtracted(3, i));
        mean_subtracted(3, i) = normalized_angle;
        mean_subtracted_transposed(i, 3) = normalized_angle;
        P_ += weights_(i) * mean_subtracted.col(i) * mean_subtracted_transposed.row(i);
    }
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
