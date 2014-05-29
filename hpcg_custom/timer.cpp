/*
 * File:   main.cpp
 * Author: Grady W. Ellison
 *
 * Created on February 7, 2014, 1:38 PM
 */

#include "timer.h"

using namespace std;

void tick(timeval &t)
{
	gettimeofday(&t, NULL);
}

double tock(timeval &t)
{
	timeval t2;
	gettimeofday(&t2, NULL);
	double time = (t2.tv_sec - t.tv_sec) * 1000.0f;
	time += (t2.tv_usec - t.tv_usec) / 1000.0f;
	return time;
}
