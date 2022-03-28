#pragma once
#include<string>
#include<iostream>

#ifndef Object_parameter_H_
#define Object_parameter_H_



class Object_parameter
{
private:
public:
	// members

	int obj_size1, obj_size2;
	double mig_F, mig_E, emit_F1, emit_E1, emit_F2, emit_E2, rotate_F, rotate_E, TM_F, TM_E, radius;


	//functions
	Object_parameter() { std::cout << "The defaulted construct function is invoked for the Object_parameter.\n"; }
	Object_parameter(int, int, double, double, double, double, double, double, double, double, double, double, double);
	~Object_parameter() {}
};

#endif // !Object_parameter

