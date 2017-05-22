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
: std_a_(30), std_yawdd_(30), std_laspx_(0.15), std_laspy_(0.15), std_radr_(0.3),
  std_radphi_(0.03), std_radrd_(0.3), n_x_(5), n_aug_(7), lambda_(3 - n_aug_), PI(3.14159265358979)
{
    //bool to invoke initialization on first "ProcessMeasurement" call
    is_initialized_ = false;

    //if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

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
    double delta_t_;    //elapsed time between the current and previous measurements (in seconds)

    //if this is the first time, we need to initialize x with the current measurement,
    //then exit (filtering will start with the next iteration)
    if (!is_initialized_)
    {
        FirstTimeInit(measurement_pack);
        //no need to predict or update so we return
        return;
    }

    //compute the time elapsed between the current and previous measurements (in seconds)
    delta_t_ = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;

    //capture the timestamp of the measurement for future use
    previous_timestamp_ = measurement_pack.timestamp_;

    //update the state transition matrix (F) according to the new elapsed time
    kf_.UpdateStateTransitionMatrix(dt);

    //update the process (noise) covariance matrix according to the new elapsed time
    kf_.UpdateProcessCovarianceMatrix(dt);

    //perform kalman prediction step
    kf_.PerformPredict();

    //perform kalman update step
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
        //compute the jacobian
        VectorXd x_state = kf_.GetState();
        Hj_ = tools_.ComputeJacobian(x_state);
        //set the appropriate H & R matrices based upon the sensor type
        kf_.SetMeasurementMatrices(Hj_, R_radar_);
        //for radar, use extended kalman filter equations
        kf_.PerformUpdateEKF(measurement_pack.raw_measurements_);
    }
    else
    {
        //set the appropriate H & R matrices based upon the sensor type
        kf_.SetMeasurementMatrices(H_laser_, R_laser_);
        //for lidar, use normal kalman filter equations
        kf_.PerformUpdate(measurement_pack.raw_measurements_);
    }

    // print the output
    //cout << "x_ = " << kf_.GetState() << endl;
    //cout << "P_ = " << kf_.GetStateCovariance() << endl;
}

//perform kalman prediction step
void UKF::Prediction(const double delta_t)
{
    //generate sigma points

    //predict sigma points

    //predict the state and state covariance matrix
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

//this is executed the first time only
void UKF::FirstTimeInit(const MeasurementPackage& measurement_pack)
{
    //local vars
    //for laser: col_1 will equal px and col_2 will equal py
    //for radar: col_1 will equal r (e.g., rho) and col_2 will equal θ (e.g., phi)
    float col_1 = measurement_pack.raw_measurements_(0);
    float col_2 = measurement_pack.raw_measurements_(1);
    float px; //x position
    float py; //y position

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
void UKF::ComputeAugmentedSigmaPoints(MatrixXd* Xsig_out)
{
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

    //populate the "right" (positive) and "left" (negative) sigma points at the same time
    for (int i = 0; i < sqrt_total.cols(); i++)
    {
        Xsig_aug.col(Xsig_aug_index) = x_aug + sqrt_total.col(i);
        Xsig_aug.col(Xsig_aug_index + sqrt_total.cols()) = x_aug - sqrt_total.col(i);
        Xsig_aug_index++;
    }

    //print result
    std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

    //write result
    *Xsig_out = Xsig_aug;
}

void UKF::SigmaPointPrediction(MatrixXd& Xsig_aug_in, const double delta_t_in, MatrixXd* Xsig_out)
{
    //create matrix with predicted sigma points as columns
     MatrixXd Xsig_pred = MatrixXd(n_x_, (2 * n_aug_) + 1);
  VectorXd state_pred = VectorXd(n_x_);
  VectorXd noise_pred = VectorXd(n_x_);
  VectorXd cur_sig;
  double delta_t_squared = delta_t_in * delta_t_in;

  //predict sigma points
  //process each sigma point prediction
  for (int i = 0; i < Xsig_aug_in.cols(); i++)
  {
      //get a handle to the first 5 elelments of the current sigma point
      cur_sig = (Xsig_aug_in.col(i)).head(5);
      //compute the process model transformation for each vector member (state)
      state_pred << ((cur_sig(2) / cur_sig(4)) * (sin(cur_sig(3) + (cur_sig(4) * delta_t_in)) - sin(cur_sig(3)))),
                    ((cur_sig(2) / cur_sig(4)) * (-cos(cur_sig(3) + (cur_sig(4) * delta_t_in)) + cos(cur_sig(3)))),
                    0,
                    cur_sig(4) * delta_t_in,
                    0;
      //compute the process model transformation for each vector member (noise)
      noise_pred << ((delta_t_squared * cos(cur_sig(3)) * (Xsig_aug_in.col(i))(5)) / 2),
                    ((delta_t_squared * sin(cur_sig(3)) * (Xsig_aug_in.col(i))(5)) / 2),
                    delta_t_in * (Xsig_aug_in.col(i))(5),
                    ((delta_t_squared * (Xsig_aug_in.col(i))(6)) / 2),
                    delta_t_in * (Xsig_aug_in.col(i))(6);
      //compute the new state for the current column
      Xsig_pred.col(i) = cur_sig + state_pred + noise_pred;
  }
  //avoid division by zero
  //write predicted sigma points into right column

  //######### need to deal with the case when psi dot is zero ##########



/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;

  //write result
  *Xsig_out = Xsig_pred;

}


void UKF::PredictMeanAndCovariance(MatrixXd& Xsig_pred_in, VectorXd* x_out, MatrixXd* P_out)
{
  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  //predict state mean
  //multiply the weights with each element in the predicted sigma points matrix
  //then do a row-wise sum of the resulting matrix so that a summed vector is returned
  x_ = (Xsig_pred_in * weights_).rowwise().sum();

  //predict state covariance matrix
  //set to zeros
  P.fill(0);
  //subtract the x vector from each column of Xsig_pred
  MatrixXd mean_subtracted = Xsig_pred_in.colwise() - x_;

  //sum the individual 5x5's
  for (int i = 0; i < (2 * n_aug_) + 1; i++)
  {
       P += weights_(i) * mean_subtracted.col(i) * (mean_subtracted.col(i)).transpose();
  }

  //print result
  std::cout << "Predicted state" << std::endl;
  std::cout << x << std::endl;
  std::cout << "Predicted covariance matrix" << std::endl;
  std::cout << P << std::endl;

  //write result
  *x_out = x;
  *P_out = P;
}

//predict measurement mean and covariance for radar
void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out)
{
    //set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;

    //init weights
  weights.fill(1 / (2 * (lambda + n_aug)));  //set all elements
  weights(0) = lambda / (lambda + n_aug);    //set first element differently

  //radar measurement noise standard deviation radius in m
  double std_radr = 0.3;

  //radar measurement noise standard deviation angle in rad
  double std_radphi = 0.0175;

  //radar measurement noise standard deviation radius change in m/s
  double std_radrd = 0.1;

  //create example matrix with predicted sigma points
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  Xsig_pred <<
         5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
           1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
          2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
         0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
          0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  double px;
  double py;
  double v;
  double psi;
  double rho;
  double phi;
  double rho_dot;

  //transform sigma points into measurement space
  //loop through the predicted sigma points and convert them to the measurement space
  for (int i = 0; i < 2 * n_aug + 1; i++)
  {
       //extract values from cur predicted sigma point column
       px = (Xsig_pred.col(i))(0);
       py = (Xsig_pred.col(i))(1);
       v = (Xsig_pred.col(i))(2);
       psi = (Xsig_pred.col(i))(3);
       //compute the rho transformation
       rho = sqrt((px * px) + (py * py));
       //compute the phi transformation
       phi = atan2(py, px);
       //compute the rho_dot transformation
       rho_dot = ((px * cos(psi) * v) + (py * sin(psi) * v)) / rho;
       //assign the computed factors to the current column of the measurement space matrix
       Zsig.col(i) << rho, phi, rho_dot;
  }

  //calculate mean predicted measurement
  //multiply the weights with each element in the measurement sigma points matrix
  //then do a row-wise sum of the resulting matrix so that a summed vector is returned
  z_pred = (Zsig * weights).rowwise().sum();

  //calculate measurement covariance matrix S
  MatrixXd mean_subtracted = Zsig.colwise() - z_pred;
  //sum the individual matrices
  for (int i = 0; i < 2 * n_aug + 1; i++)
  {
       S += (weights(i) * mean_subtracted.col(i) * (mean_subtracted.col(i)).transpose());
  }

  //create the measurment noise matrix
  MatrixXd R = MatrixXd(3, 3);
  R.fill(0); //fill with zeros
  R(0, 0) = std_radr * std_radr;
  R(1, 1) = std_radphi * std_radphi;
  R(2, 2) = std_radrd * std_radrd;

  //add the noise matrix to the covariance matrix
  S += R;



/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  std::cout << "S: " << std::endl << S << std::endl;

  //write result
  *z_out = z_pred;
  *S_out = S;
}


void UKF::UpdateState(VectorXd* x_out, MatrixXd* P_out) {

  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //define spreading parameter
  double lambda = 3 - n_aug;

  //set vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);
  //init weights
  weights.fill(1 / (2 * (lambda + n_aug)));  //set all elements
  weights(0) = lambda / (lambda + n_aug);    //set first element differently

  //create example matrix with predicted sigma points in state space
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  Xsig_pred <<
         5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
           1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
          2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
         0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
          0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

  //create example vector for predicted state mean
  VectorXd x = VectorXd(n_x);
  x <<
     5.93637,
     1.49035,
     2.20528,
    0.536853,
    0.353577;

  //create example matrix for predicted state covariance
  MatrixXd P = MatrixXd(n_x,n_x);
  P <<
  0.0054342,  -0.002405,  0.0034157, -0.0034819, -0.00299378,
  -0.002405,    0.01084,   0.001492,  0.0098018,  0.00791091,
  0.0034157,   0.001492,  0.0058012, 0.00077863, 0.000792973,
 -0.0034819,  0.0098018, 0.00077863,   0.011923,   0.0112491,
 -0.0029937,  0.0079109, 0.00079297,   0.011249,   0.0126972;

  //create example matrix with sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);
  Zsig <<
      6.1190,  6.2334,  6.1531,  6.1283,  6.1143,  6.1190,  6.1221,  6.1190,  6.0079,  6.0883,  6.1125,  6.1248,  6.1190,  6.1188,  6.12057,
     0.24428,  0.2337, 0.27316, 0.24616, 0.24846, 0.24428, 0.24530, 0.24428, 0.25700, 0.21692, 0.24433, 0.24193, 0.24428, 0.24515, 0.245239,
      2.1104,  2.2188,  2.0639,   2.187,  2.0341,  2.1061,  2.1450,  2.1092,  2.0016,   2.129,  2.0346,  2.1651,  2.1145,  2.0786,  2.11295;

  //create example vector for mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred <<
      6.12155,
     0.245993,
      2.10313;

  //create example matrix for predicted measurement covariance
  MatrixXd S = MatrixXd(n_z,n_z);
  S <<
      0.0946171, -0.000139448,   0.00407016,
   -0.000139448,  0.000617548, -0.000770652,
     0.00407016, -0.000770652,    0.0180917;

  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z <<
      5.9214,   //rho in m
      0.2187,   //phi in rad
      2.0062;   //rho_dot in m/s

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //calculate cross correlation matrix
  MatrixXd Xmean_subtracted = Xsig_pred.colwise() - x;
  MatrixXd Zmean_subtracted = (Zsig.colwise() - z_pred).transpose();
  Tc.fill(0); //init to zeros as we're adding to existing values
  //sum the individual matrices
  for (int i = 0; i < 2 * n_aug + 1; i++)
  {
       Tc += weights(i) * Xmean_subtracted.col(i) * Zmean_subtracted.row(i);
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //update state mean and covariance matrix
  x += K * (z - z_pred);
  P -= K * S * K.transpose();

/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  std::cout << "Updated state x: " << std::endl << x << std::endl;
  std::cout << "Updated state covariance P: " << std::endl << P << std::endl;

  //write result
  *x_out = x;
  *P_out = P;
}
