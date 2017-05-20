/*
#######################################################
## AUTHOR: James Beasley                             ##
## DATE: May 8, 2017                                 ##
## UDACITY SDC: Project 7 (Unscented Kalman Filters) ##
#######################################################
*/

#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "Eigen/Dense"

class Tools
{
    public:
        //constructor
        Tools();
        //destructor
        virtual ~Tools();

        //compute RMSE
        Eigen::VectorXd ComputeRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);
};

#endif /* TOOLS_H_ */
