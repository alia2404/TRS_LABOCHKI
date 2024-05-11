#include <iostream>
#include "eigen-3.4.0\Eigen\Core"
#include "my_data.h"
#include "task1.h"
#include "task2.h"
#include "task3.h"
#include "task4.h"

using namespace std;
using namespace Eigen;

/*
* 
* 
*/

int main()
{
    std::cout << "Hello World!\n";
    std::cout << "Eigen version : " << EIGEN_MAJOR_VERSION << "."
        << EIGEN_MINOR_VERSION << endl;

    int iter_num = 5;
    double eps = 1e-9;
    vector<double> x_arr = { 0.1, 0.05, 0.02, 0.01 };
    
    for (int i = 0; i < x_arr.size(); i++) {
        double x_step = x_arr[i];
        double y_step = x_step * LY / LX;
        std::string filename = "1mur_" + std::to_string(i);
        t1(x_step, y_step, eps, iter_num, ERR_CONST, filename);
    }
    for (int i = 0; i < x_arr.size(); i++) {
        double x_step = x_arr[i];
        double y_step = x_step * LY / LX;
        std::string filename = "2mur_" + std::to_string(i);
        t2(x_step, y_step, eps, iter_num, ERR_CONST, filename);
    }
    /*
    for (int i = 0; i < x_arr.size() - 2; i++) {
        double x_step = x_arr[i];
        double y_step = x_step * LY / LX;
        std::string filename = "3_" + std::to_string(i);
        t3(x_step, y_step, err, iter_num, ERR_CONST, filename);
    }
    for (int i = 2; i < x_arr.size(); i++) {
        double x_step = x_arr[i];
        double y_step = x_step * LY / LX;
        std::string filename = "4_" + std::to_string(i);
        t4(x_step, y_step, err, iter_num, ERR_CONST, filename);
    }
    */
    /*
    for (int i = 2; i < 3; i++) {
        double x_step = x_arr[i];
        double y_step = x_step * LY / LX;
        std::string filename = "2_" + std::to_string(i);
        t2(x_step, y_step, err, iter_num, ERR_CONST, filename);
    }
    for (int i = 2; i < 3; i++) {
        double x_step = x_arr[i];
        double y_step = x_step * LY / LX;
        std::string filename = "3_" + std::to_string(i);
        t3(x_step, y_step, err, iter_num, ERR_CONST, filename);
    }
    */
    /*
    for (int i = 1; i < 4; i++) {
        double x_step = x_arr[i];
        double y_step = x_step * LY / LX;
        std::string filename = "2nom_" + std::to_string(i);
        err = 0.000001;
        t2(x_step, y_step, err, iter_num, ERR_CONST, filename);
    }
    */

    //t1(x_step, y_step, err, iter_num, ITER_CONST, "1_check");
    return 0;
}