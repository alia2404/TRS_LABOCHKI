// trs_lab1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "task1.h"
#include "task2.h"
#include "task3.h"
#include "task4.h"
#include "task5.h"

int main()
{
    
    std::vector<double> x_step_arr = { 1e-1, 1e-2, 1e-3 };
    std::vector<double> t_step_arr = { 1e-1, 1e-2, 1e-3/*, 1e-4*/};
    std::vector<double> alpha_arr = { 0.8, 0.9, 1.0, 1.1, 1.2};
    /*
    // Task 1
    for (int i = 0; i < x_step_arr.size()-1; i++) {
        for (int j = 0; j < alpha_arr.size(); j++) {
            double x_step = x_step_arr[i];
            double alpha = alpha_arr[j];
            double t_step = alpha * x_step * x_step / 2;
            double stoptime1 = 1;
            std::string filename = "1_" + std::to_string(i) + "_" + std::to_string(j);
            t1(x_step, t_step, stoptime1, filename);
            std::cout << filename << " saved" << std::endl;
        }
    }
    */

    /*
    // Task 2 (implicit method)
    for (int i = 0; i < x_step_arr.size(); i++) {
        for (int j = 0; j < t_step_arr.size(); j++) {

            //double stoptime2 = 1;
            //double t_step2 = 1e-2;
            double x_step = x_step_arr[i];
            double t_step2 = t_step_arr[j];
            double stoptime2 = 1;
            std::string filename = "2_" + std::to_string(i) + "_" + std::to_string(j);
            t2(x_step, t_step2, stoptime2, filename);
            std::cout << filename << " saved" << std::endl;
        }
    }
    */

    /*
    // Task 2 (Crank-Nicolson method)
    for (int i = 0; i < x_step_arr.size(); i++) {
        for (int j = 0; j < t_step_arr.size(); j++) {
            //double stoptime2 = 1;
            //double t_step2 = 1e-2;
            double x_step = x_step_arr[i];
            double t_step3 = t_step_arr[j];
            double stoptime3 = 1;
            std::string filename = "3_" + std::to_string(i) + "_" + std::to_string(j);
            t3(x_step, t_step3, stoptime3, filename);
            std::cout << filename << " saved" << std::endl;
        }
    }
    */

    /*
    // Task 3. Check linear
    std::vector<double> sigma_arr = { 0., 0.5, 1. };
    for (int i = 0; i < sigma_arr.size(); i++) {
        double sigma = sigma_arr[i];
        double x_step = 1e-2;
        double t_step1 = 0.9 * x_step * x_step / 2;
        double stoptime = 1.;
        std::string filename = "4_check_" + std::to_string(i);
        t4(x_step, t_step1, stoptime, sigma, filename);
    }
    */

    /*
    // Task 3. Nonlinear
    for (int i = 0; i < x_step_arr.size()-1; i++) {
        for (int j = 0; j < alpha_arr.size()-2; j++) {
            double x_step = x_step_arr[i];
            double alpha = alpha_arr[j];
            double t_step = alpha * x_step * x_step / 2;
            double stoptime1 = 1;
            std::string filename = "4_" + std::to_string(i) + "_" + std::to_string(j);
            t4(x_step, t_step, stoptime1, 0.5, filename);
            std::cout << filename << " saved" << std::endl;
        }
    }
    */
    
    /*
    // Task 4    
    // check on linear
    double delta_abs = 0.0001;
    for (int i = 0; i < t_step_arr.size()-1; i++) {
        double x_step = 0.01;
        double t_step = t_step_arr[i+1];
        double stoptime = 1.;
        std::string filename = "5_check_" + std::to_string(i);
        t5(x_step, t_step, stoptime, delta_abs, filename);
        std::cout << filename << " saved" << std::endl;
    }
    */

    /*
    // Task 4
    // search steps on linear
    for (int i = 0; i < x_step_arr.size()-1; i++) {
        for (int j = 0; j < t_step_arr.size()-1; j++) {
            double x_step = x_step_arr[i+1];
            double t_step = t_step_arr[j+1];
            double stoptime1 = 1.;
            std::string filename = "4_search_" + std::to_string(i) + "_" + std::to_string(j);
            t4(x_step, t_step, stoptime1, 0.5, filename);
            std::cout << filename << " saved" << std::endl;
        }
    }
    double delta_abs = 0.001;
    for (int i = 0; i < x_step_arr.size() - 1; i++) {
        for (int j = 0; j < t_step_arr.size() - 1; j++) {
            double x_step = x_step_arr[i + 1];
            double t_step = t_step_arr[j + 1];
            double stoptime1 = 1;
            std::string filename = "5_search_" + std::to_string(i) + "_" + std::to_string(j);
            t5(x_step, t_step, stoptime1, delta_abs, filename);
            std::cout << filename << " saved" << std::endl;
        }
    }
    */


    // Task 4
    // compare methods
    double delta_rel = 0.00001;
    t4(0.001, 0.001, 1., 0.5, "4_speed");
    std::cout << "4 saved" << std::endl;
    t5(0.001, 0.001, 1., delta_rel*10, "5_speed");
    std::cout << "5 saved" << std::endl;
    t5(0.001, 0.001, 1., delta_rel, "5_speed2");
    std::cout << "5_s2 saved" << std::endl;

}