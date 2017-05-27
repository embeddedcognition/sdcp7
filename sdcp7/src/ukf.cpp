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
: std_a_(0.8), std_yawdd_(1), std_laspx_(0.15), std_laspy_(0.15), std_radr_(0.3),
  std_radphi_(0.03), std_radrd_(0.3), n_x_(5), n_aug_(7), lambda_(3 - n_aug_), PI(3.14159265358979)
{
    //bool to invoke initialization on first "ProcessMeasurement" call
    is_initialized_ = false;

    //if this is false, laser measurements will be ignored (except during init)
    use_laser_ = false;

    //if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    //initial state vector
    x_ = VectorXd(5);

    //initial covariance matrix
    P_ = MatrixXd(5, 5);

    //predicted sigma points
    Xsig_pred_ = MatrixXd(n_x_, (2 * n_aug_) + 1);

    //weights for sigma points
    weights_ = VectorXd((2 * n_aug_) + 1);
    //init weights
    weights_.fill(1 / (2 * (lambda_ + n_aug_)));  //set all elements
    weights_(0) = lambda_ / (lambda_ + n_aug_);    //set first element differently

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

//perform kalman prediction step
//at this point we have the posterior mean state and covariance matrix from the last iteration (representing distribution of current state)
//the model we're using (CTRV) has a state vector of size 5
void UKF::PerformPrediction(const double delta_t)
{
    //local vars
    MatrixXd Xsig_aug;      //augmented sigma points (augmented state)

    //represent the uncertainty of the posterior state estimation with sigma points
    //the augmented state includes the noise vector (last two elements)
    //first sigma point is always the mean state estimate
    ComputeAugmentedSigmaPoints(&Xsig_aug);

    //predict sigma points (by inserting each augmented sigma point into the CTRV process model)
    PredictSigmaPoints(Xsig_aug, delta_t);

    //compute the mean and covariance of the predicted state
    ComputeMeanAndCovarianceofPredictedSigmaPoints();
}

//update state and covariance using lidar measurement
void UKF::UpdateLidar(const MeasurementPackage& measurement_pack)
{
    //local vars
    MatrixXd S;         //measurement noise covariance
    VectorXd z_pred;    //mean predicted measurement
    MatrixXd Zsig;      //predicted sigma points mappined into the lidar measurement space
    int n_z = 2;        //set measurement dimension, lidar can measure px and py

    //transform the predicted state into the lidar measurement space and then predict measurement mean and covariance
    PredictLidarMeasurement(n_z, &z_pred, &S, &Zsig);

    //now that we have a predicted mean and covariance for the specific measurement space (radar or lidar), we can perform the kalman update step
    //computing the resulting mean and covariance of the predicted measurement
    UpdateState(measurement_pack.raw_measurements_, z_pred, Zsig, S, n_z);

    //compute radar NIS
}

//update state and covariance using radar measurement
void UKF::UpdateRadar(const MeasurementPackage& measurement_pack)
{
    //local vars
    MatrixXd S;         //measurement noise covariance
    VectorXd z_pred;    //mean predicted measurement
    MatrixXd Zsig;      //predicted sigma points mappined into the radar measurement space
    int n_z = 3;        //set measurement dimension, radar can measure r, phi, and r_dot

    //transform the predicted state into the radar measurement space and then predict measurement mean and covariance
    PredictRadarMeasurement(n_z, &z_pred, &S, &Zsig);

    //now that we have a predicted mean and covariance for the specific measurement space (radar or lidar), we can perform the kalman update step
    //computing the resulting mean and covariance of the predicted measurement
    UpdateState(measurement_pack.raw_measurements_, z_pred, Zsig, S, n_z);

    //compute radar NIS
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
    //since we'll not have enough information to initialize the velocity portion of the state vector (vx and vy),
    //i.e., we only have the current position to go on, we'll set it to zero
    x_ << px, py, 0, 0, 0;

    //init the matrix P to the identity matrix
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    //capture the timestamp of the measurement for future use
    previous_timestamp_ = measurement_pack.timestamp_;

    //done initializing
    is_initialized_ = true;
}

//create augmented sigma points
//allows us to represent the uncertainty of the covariance matrix Q with sigma points
void UKF::ComputeAugmentedSigmaPoints(MatrixXd* Xsig_aug_out)
{
    //local vars
    MatrixXd Q = MatrixXd(2, 2);                                //process noise covariance matrix
    VectorXd x_aug = VectorXd(7);                               //augmented mean vector
    MatrixXd P_aug = MatrixXd(7, 7);                            //augmented state covariance
    MatrixXd Xsig_aug = MatrixXd(n_aug_, (2 * n_aug_) + 1);     //augmented sigma point matrix (state + process noise)

    //compute augmented mean state
    x_aug.head(x_.size()) = x_; //copy x into first 5 elements of x_aug
    x_aug(5) = 0; //copy mean value of acceleration noise
    x_aug(6) = 0; //copy mean value of acceleration noise

    //init process noise covariance matrix
    //we're supplied the std dev for longitudinal acceleration and yaw acceleration, but we must square to get variance
    Q << (std_a_ * std_a_), 0,
         0, (std_yawdd_ * std_yawdd_);

    //compute augmented state covariance
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

    //create augmented sigma points
    int Xsig_aug_index = 1; //index starts at second column

    //populate the augmented mean
    Xsig_aug.col(0) = x_aug;

    //populate the sigma points
    for (int i = 0; i < sqrt_total.cols(); i++)
    {
        //sigma points to the right of the mean state (directionally)
        Xsig_aug.col(Xsig_aug_index) = x_aug + sqrt_total.col(i);
        //sigma points to the left of the mean (directionally)
        Xsig_aug.col(Xsig_aug_index + sqrt_total.cols()) = x_aug - sqrt_total.col(i);
        Xsig_aug_index++;
    }

    //output augmented sigma points
    *Xsig_aug_out = Xsig_aug;
}

//perform k + 1 prediction using augmented sigma points
//process noise has a non-linear effect on the state in the process model, therefore augmented sigma points are required
void UKF::PredictSigmaPoints(MatrixXd& Xsig_aug_in, const double delta_t_in)
{
    //local vars
    VectorXd state_pred = VectorXd(n_x_);                   //predicted state (calculated via the CTRV model)
    VectorXd noise_pred = VectorXd(n_x_);                   //predicted noise (calculated via the CTRV model)
    double delta_t_squared = delta_t_in * delta_t_in;       //delta_t squared is a multi-use operation (so we compute it once)
    double px, py, v, yaw, yaw_dot, nu_a, nu_yaw_dot_dot;   //used to hold current sigma point column values
    double predicted_px, predicted_py;                      //used to hold predicted px and py values

    //predict sigma points
    //process each sigma point prediction
    for (int i = 0; i < Xsig_aug_in.cols(); i++)
    {
        //extract augmented state column values for easier handling/tracking
        px = Xsig_aug_in(0, i);
        py = Xsig_aug_in(1, i);
        v = Xsig_aug_in(2, i);
        yaw = Xsig_aug_in(3, i);
        yaw_dot = Xsig_aug_in(4, i);
        nu_a = Xsig_aug_in(5, i);
        nu_yaw_dot_dot = Xsig_aug_in(6, i);

        //deal with potential of yaw_dot being zero in computing state prediction
        if (fabs(yaw_dot) > 0.0001)
        {
            //compute predicted px and py via CTRV when yaw_dot is greater than zero
            predicted_px = (v / yaw_dot) * (sin(yaw + (yaw_dot * delta_t_in)) - sin(yaw));
            predicted_py = (v / yaw_dot) * (-cos(yaw + (yaw_dot * delta_t_in)) + cos(yaw));
        }
        else
        {
            //compute predicted px and py via CTRV when yaw_dot is zero
            predicted_px = v * cos(yaw) * delta_t_in;
            predicted_py = v * sin(yaw) * delta_t_in;
        }

        //compute the process model transformation for each vector member (state)
        state_pred << predicted_px,
                      predicted_py,
                      0,
                      yaw_dot * delta_t_in,
                      0;
        //compute the process model transformation for each vector member (noise)
        noise_pred << (delta_t_squared * cos(yaw) * nu_a) / 2,
                      (delta_t_squared * sin(yaw) * nu_a) / 2,
                      delta_t_in * nu_a,
                      (delta_t_squared * nu_yaw_dot_dot) / 2,
                      delta_t_in * nu_yaw_dot_dot;

        //populate process model state elements from current sigma point column (i.e., extract first 5 elements of current sigma point column as they align with the process model specification)
        Xsig_pred_.col(i) << px, py, v, yaw, yaw_dot;

        //add sum of state and noise prediction to the current state values to make up the new state (stored in a global var Xsig_pred_)
        Xsig_pred_.col(i) += (state_pred + noise_pred);
    }
}

//last step in delivering the new predicted state and covariance, we need to compute the mean and covariance of the predicted state
void UKF::ComputeMeanAndCovarianceofPredictedSigmaPoints()
{
    //predict state mean
    //multiply the weights with each element in the predicted sigma points matrix
    //then do a row-wise sum of the resulting matrix so that a summed vector is returned (stored in a global var x_)
    x_ = (Xsig_pred_ * weights_).rowwise().sum();

    //predict state covariance matrix
    //set to zeros
    P_.fill(0);
    //subtract the x vector from each column of Xsig_pred
    MatrixXd mean_subtracted = Xsig_pred_.colwise() - x_;

    //sum the individual 5x5's, result is stored in a global var P_)
    for (int i = 0; i < (2 * n_aug_) + 1; i++)
    {
        P_ += weights_(i) * mean_subtracted.col(i) * (mean_subtracted.col(i)).transpose();
    }
}

//transform the predicted state into the radar measurement space and then predict measurement mean and covariance
void UKF::PredictRadarMeasurement(const int& n_z_in, VectorXd* z_pred_out, MatrixXd* S_out, MatrixXd* Zsig_out)
{
    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z_in, (2 * n_aug_) + 1);
    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z_in);
    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z_in, n_z_in);
    //used to extract particular values for easier handling
    double px, py, v, psi, rho, phi, rho_dot;
    //used for angle normalization
    double* cur_phi;

    //transform sigma points into measurement space
    //loop through the predicted sigma points and convert them to the measurement space
    for (int i = 0; i < (2 * n_aug_) + 1; i++)
    {
        //extract values from cur predicted sigma point column
        px = (Xsig_pred_.col(i))(0);
        py = (Xsig_pred_.col(i))(1);
        v = (Xsig_pred_.col(i))(2);
        psi = (Xsig_pred_.col(i))(3);

        //compute the rho transformation
        rho = sqrt((px * px) + (py * py));
        //compute the phi transformation
        phi = atan2(py, px);
        //compute the rho_dot transformation
        rho_dot = ((px * cos(psi) * v) + (py * sin(psi) * v)) / rho;

        //assign the computed factors to the current column of the measurement space matrix
        Zsig.col(i) << rho, phi, rho_dot;
    }

    //compute mean predicted measurement
    //multiply the weights with each element in the measurement sigma points matrix
    //then do a row-wise sum of the resulting matrix so that a summed vector is returned
    z_pred = (Zsig * weights_).rowwise().sum();

    //calculate measurement covariance matrix S
    //save time and space by computing the mean and transpose ahead of time
    MatrixXd mean_subtracted = Zsig.colwise() - z_pred;
    MatrixXd mean_subtracted_transposed = mean_subtracted.transpose();
    S.fill(0); //init with zeros since we're adding to it each time
    //sum the individual matrices
    for (int i = 0; i < (2 * n_aug_) + 1; i++)
    {
        //ensure phi is between -pi and pi
        cur_phi = &((mean_subtracted.col(i))(1));   //get a handle to the value of phi in the current column so it can be normalized
        *cur_phi = NormalizeAngle(*cur_phi);

        //ensure phi is between -pi and pi
        cur_phi = &((mean_subtracted_transposed.row(i))(1));   //get a handle to the value of phi in the current row (as this version is transposed) so it can be normalized
        *cur_phi = NormalizeAngle(*cur_phi);

        S += (weights_(i) * mean_subtracted.col(i) * mean_subtracted_transposed.row(i));
    }

    //create the measurment noise matrix
    MatrixXd R = MatrixXd(3, 3);
    R.fill(0); //fill with zeros
    R(0, 0) = std_radr_ * std_radr_;
    R(1, 1) = std_radphi_ * std_radphi_;
    R(2, 2) = std_radrd_ * std_radrd_;

    //add the noise matrix to the covariance matrix
    S += R;

    //output predicted measurement mean and covariance, as well as the sigma points that are now mapped into the radar measurement space
    *z_pred_out = z_pred;
    *S_out = S;
    *Zsig_out = Zsig;
}

//transform the predicted state into the lidar measurement space and then predict measurement mean and covariance
void UKF::PredictLidarMeasurement(const int& n_z_in, VectorXd* z_pred_out, MatrixXd* S_out, MatrixXd* Zsig_out)
{
    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z_in, (2 * n_aug_) + 1);
    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z_in);
    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z_in, n_z_in);
    //used to extract particular values for easier handling
    double px, py;

    //transform sigma points into measurement space
    //loop through the predicted sigma points and convert them to the measurement space
    for (int i = 0; i < (2 * n_aug_) + 1; i++)
    {
        //extract values from cur predicted sigma point column
        px = (Xsig_pred_.col(i))(0);
        py = (Xsig_pred_.col(i))(1);

        //assign the computed factors to the current column of the measurement space matrix
        Zsig.col(i) << px, py;
    }

    //compute mean predicted measurement
    //multiply the weights with each element in the measurement sigma points matrix
    //then do a row-wise sum of the resulting matrix so that a summed vector is returned
    z_pred = (Zsig * weights_).rowwise().sum();

    //calculate measurement covariance matrix S
    //save time and space by computing the mean and transpose ahead of time
    MatrixXd mean_subtracted = Zsig.colwise() - z_pred;
    MatrixXd mean_subtracted_transposed = mean_subtracted.transpose();
    S.fill(0); //init with zeros since we're adding to it each time
    //sum the individual matrices
    for (int i = 0; i < (2 * n_aug_) + 1; i++)
    {
        S += (weights_(i) * mean_subtracted.col(i) * mean_subtracted_transposed.row(i));
    }

    //create the measurment noise matrix
    MatrixXd R = MatrixXd(2, 2);
    R.fill(0); //fill with zeros
    R(0, 0) = std_laspx_ * std_laspx_;
    R(1, 1) = std_laspy_ * std_laspy_;

    //add the noise matrix to the covariance matrix
    S += R;

    //output predicted measurement mean and covariance, as well as the sigma points that are now mapped into the radar measurement space
    *z_pred_out = z_pred;
    *S_out = S;
    *Zsig_out = Zsig;
}

//final step for this iteration, update the predicted state and covariance using the radar/lidar measurement
void UKF::UpdateState(const VectorXd& z_in, const VectorXd& z_pred_in, const MatrixXd& Zsig_in, const MatrixXd& S_in, const int& n_z_in)
{
    //compute means
    MatrixXd Xmean_subtracted = Xsig_pred_.colwise() - x_;
    MatrixXd Zmean_subtracted = (Zsig_in.colwise() - z_pred_in).transpose();

    //compute cross correlation matrix
    MatrixXd Tc = MatrixXd(n_x_, n_z_in);
    Tc.fill(0); //init to zeros as we're adding to existing values
    //sum the individual matrices
    for (int i = 0; i < (2 * n_aug_) + 1; i++)
    {
        //ensure psi is between -pi and pi
        double* cur_psi = &((Xmean_subtracted.col(i))(3));   //get a handle to the value of psi in the current column so it can be normalized
        *cur_psi = NormalizeAngle(*cur_psi);

        //ensure phi is between -pi and pi
        double* cur_phi = &((Zmean_subtracted.row(i))(1));   //get a handle to the value of phi in the current row (as it is transposed) so it can be normalized
        *cur_phi = NormalizeAngle(*cur_phi);

        Tc += weights_(i) * Xmean_subtracted.col(i) * Zmean_subtracted.row(i);
    }

    //compute Kalman gain K;
    MatrixXd K = Tc * S_in.inverse();

    //update state mean and covariance matrix
    x_ += K * (z_in - z_pred_in);
    P_ -= K * S_in * K.transpose();
}
