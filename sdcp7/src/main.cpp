/*
#######################################################
## AUTHOR: James Beasley                             ##
## DATE: May 8, 2017                                 ##
## UDACITY SDC: Project 7 (Unscented Kalman Filters) ##
#######################################################
*/

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include "Eigen/Dense"
#include "ukf.h"
#include "ground_truth_package.h"
#include "measurement_package.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

//function declarations
void check_arguments(int argc, char* argv[]);
void check_files(ifstream& in_file, string& in_name, ofstream& out_file, string& out_name);

//main thread
int main(int argc, char* argv[])
{
    //argument check
    check_arguments(argc, argv);

    string in_file_name_ = argv[1];
    ifstream in_file_(in_file_name_.c_str(), ifstream::in);
    string out_file_name_ = argv[2];
    ofstream out_file_(out_file_name_.c_str(), ofstream::out);

    //file check
    check_files(in_file_, in_file_name_, out_file_, out_file_name_);

    vector<MeasurementPackage> measurement_pack_list;
    vector<GroundTruthPackage> gt_pack_list;

    string line;

    //prep the measurement packages (each line represents a measurement at a timestamp)
    while (getline(in_file_, line))
    {
        string sensor_type;
        MeasurementPackage meas_package;
        GroundTruthPackage gt_package;
        istringstream iss(line);
        long long timestamp;

        //reads first element from the current line
        iss >> sensor_type;
        if (sensor_type.compare("L") == 0)
        {
            //LASER MEASUREMENT
            //read measurements at this timestamp
            meas_package.sensor_type_ = MeasurementPackage::LASER;
            meas_package.raw_measurements_ = VectorXd(2);
            float x;
            float y;
            iss >> x;
            iss >> y;
            meas_package.raw_measurements_ << x, y;
            iss >> timestamp;
            meas_package.timestamp_ = timestamp;
            measurement_pack_list.push_back(meas_package);
        }
        else if (sensor_type.compare("R") == 0)
        {
            //RADAR MEASUREMENT
            //read measurements at this timestamp
            meas_package.sensor_type_ = MeasurementPackage::RADAR;
            meas_package.raw_measurements_ = VectorXd(3);
            float ro;
            float phi;
            float ro_dot;
            iss >> ro;
            iss >> phi;
            iss >> ro_dot;
            meas_package.raw_measurements_ << ro, phi, ro_dot;
            iss >> timestamp;
            meas_package.timestamp_ = timestamp;
            measurement_pack_list.push_back(meas_package);
        }

        //read ground truth data to compare later
        float x_gt;
        float y_gt;
        float vx_gt;
        float vy_gt;
        iss >> x_gt;
        iss >> y_gt;
        iss >> vx_gt;
        iss >> vy_gt;
        gt_package.gt_values_ = VectorXd(4);
        gt_package.gt_values_ << x_gt, y_gt, vx_gt, vy_gt;
        gt_pack_list.push_back(gt_package);
    }

    //create a UKF instance
    UKF ukf;

    //used to compute the RMSE later
    vector<VectorXd> estimations;
    vector<VectorXd> ground_truth;

    //get measurement count
    size_t N = measurement_pack_list.size();

    // column names for output file
    out_file_ << "time_stamp" << "\t";
    out_file_ << "px_state" << "\t";
    out_file_ << "py_state" << "\t";
    out_file_ << "v_state" << "\t";
    out_file_ << "yaw_angle_state" << "\t";
    out_file_ << "yaw_rate_state" << "\t";
    out_file_ << "sensor_type" << "\t";
    out_file_ << "NIS" << "\t";
    out_file_ << "px_measured" << "\t";
    out_file_ << "py_measured" << "\t";
    out_file_ << "px_ground_truth" << "\t";
    out_file_ << "py_ground_truth" << "\t";
    out_file_ << "vx_ground_truth" << "\t";
    out_file_ << "vy_ground_truth" << "\n";

    //start the UKF-based fusion
    for (size_t k = 0; k < N; ++k)
    {
        //start filtering from the second frame (the speed is unknown in the first frame)
        ukf.ProcessMeasurement(measurement_pack_list[k]);

        //get x state vector
        VectorXd x_state = ukf.GetState();

        // timestamp
        out_file_ << measurement_pack_list[k].timestamp_ << "\t"; // pos1 - est

        // output the state vector
        out_file_ << x_state(0) << "\t"; // pos1 - est
        out_file_ << x_state(1) << "\t"; // pos2 - est
        out_file_ << x_state(2) << "\t"; // vel_abs -est
        out_file_ << x_state(3) << "\t"; // yaw_angle -est
        out_file_ << x_state(4) << "\t"; // yaw_rate -est

        //output the measurements
        if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::LASER)
        {
            //sensor type
            out_file_ << "lidar" << "\t";

            //get NIS value
            double nis_laser = ukf.GetLaserNIS();

            //NIS value
            out_file_ << nis_laser << "\t";

            //output the lidar measurement
            out_file_ << measurement_pack_list[k].raw_measurements_(0) << "\t";
            out_file_ << measurement_pack_list[k].raw_measurements_(1) << "\t";
        }
        else if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::RADAR)
        {
            //sensor type
            out_file_ << "radar" << "\t";

            //get NIS value
            double nis_radar = ukf.GetRadarNIS();

            //NIS value
            out_file_ << nis_radar << "\t";

            //output the radar measurement in the cartesian coordinates
            float ro = measurement_pack_list[k].raw_measurements_(0);
            float phi = measurement_pack_list[k].raw_measurements_(1);
            out_file_ << ro * cos(phi) << "\t"; // p1_meas
            out_file_ << ro * sin(phi) << "\t"; // ps_meas
        }

        //output the ground truth packages
        out_file_ << gt_pack_list[k].gt_values_(0) << "\t";
        out_file_ << gt_pack_list[k].gt_values_(1) << "\t";
        out_file_ << gt_pack_list[k].gt_values_(2) << "\t";
        out_file_ << gt_pack_list[k].gt_values_(3) << "\n";

        //convert ukf x vector to cartesian to compare to ground truth
        VectorXd ukf_x_cartesian_ = VectorXd(4);

        float x_estimate_ = x_state(0);
        float y_estimate_ = x_state(1);
        float vx_estimate_ = x_state(2) * cos(x_state(3));
        float vy_estimate_ = x_state(2) * sin(x_state(3));

        ukf_x_cartesian_ << x_estimate_, y_estimate_, vx_estimate_, vy_estimate_;

        estimations.push_back(ukf_x_cartesian_);

        ground_truth.push_back(gt_pack_list[k].gt_values_);
    }

    //compute the accuracy (RMSE)
    Tools tools;
    cout << "Accuracy - RMSE:" << endl << tools.ComputeRMSE(estimations, ground_truth) << endl;

    //close files
    if (out_file_.is_open())
    {
        out_file_.close();
    }

    if (in_file_.is_open())
    {
        in_file_.close();
    }

    cout << "Done!" << endl;
    return 0;
}

//check input arguments
void check_arguments(int argc, char* argv[])
{
    //local vars
    string usage_instructions = "Usage instructions: ";
    usage_instructions += argv[0];
    usage_instructions += " path/to/input.txt output.txt";
    bool has_valid_args = false;

    // make sure the user has provided input and output files
    if (argc == 1)
    {
        cerr << usage_instructions << endl;
    }
    else if (argc == 2)
    {
        cerr << "Please include an output file.\n" << usage_instructions << endl;
    }
    else if (argc == 3)
    {
        has_valid_args = true;
    }
    else if (argc > 3)
    {
        cerr << "Too many arguments.\n" << usage_instructions << endl;
    }

    if (!has_valid_args)
    {
        exit(EXIT_FAILURE);
    }
}

//check input and output files
void check_files(ifstream& in_file, string& in_name, ofstream& out_file, string& out_name)
{
    if (!in_file.is_open())
    {
        cerr << "Cannot open input file: " << in_name << endl;
        exit(EXIT_FAILURE);
    }

    if (!out_file.is_open())
    {
        cerr << "Cannot open output file: " << out_name << endl;
        exit(EXIT_FAILURE);
    }
}
