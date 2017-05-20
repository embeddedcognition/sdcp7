/*
#######################################################
## AUTHOR: James Beasley                             ##
## DATE: May 8, 2017                                 ##
## UDACITY SDC: Project 7 (Unscented Kalman Filters) ##
#######################################################
*/

#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;

//constructor
Tools::Tools() {}

//destructor
Tools::~Tools() {}

//compute the rmse of our estimations vs. ground truth
VectorXd Tools::ComputeRMSE(const vector<VectorXd>& estimations, const vector<VectorXd>& ground_truth)
{
    //local vars
    VectorXd rmse(4);

    //init rmse (returned if we have invalid input)
    rmse << 0, 0, 0, 0;

    //check the validity of the following inputs:
    //* the estimation vector size should not be zero
    //* the estimation vector size should equal ground truth vector size
    if ((estimations.size() != ground_truth.size()) || (estimations.size() == 0))
    {
        cout << "Invalid estimation or ground_truth data." << endl;
        return rmse;
    }

    //accumulate squared residuals
    for (unsigned int i = 0; i < estimations.size(); ++i)
    {
        //compute residual
        VectorXd residual = estimations[i] - ground_truth[i];

        //coefficient-wise multiplication
        residual = residual.array() * residual.array();
        //add residual
        rmse += residual;
    }

    //calculate the mean
    rmse = rmse / estimations.size();

    //calculate the squared root
    rmse = rmse.array().sqrt();

    //return the resulting rmse
    return rmse;
}
