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
    // ������� 1�� ���
    // �������� ���
    double step1 = 1e-10;
    // ������������ ��� (��� ��. �����)
    double step2 = step1 * 0.5;
    double stoptime = 4e-5;

    t1(step1, stoptime, "../../../../my_girls_2.0/LAB_1/1_step1");
    // t1(step2, stoptime, "../../../../my_girls_2.0/LAB_1/1_step2");
  
    // t2(step1, stoptime, "../../../../my_girls_2.0/LAB_1/2");
    
    // t3(step1, stoptime, "../../../../my_girls_2.0/LAB_1/3");
    
    int n = 1e2;
    double x_step;
    // t4(n, -10, 10, "../../../../my_girls_2.0/LAB_1/4.csv", x_step);

    n = 1e4;
    t5(n, -10, 10, "../../../../my_girls_2.0/LAB_1/5.csv", x_step);

}