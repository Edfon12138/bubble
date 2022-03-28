#pragma once

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include<ctime>
#include"object.h"
#include"setting.h"

void initialize_output_file(std::string _filename);

void initialize_react_file();

// 在固定的时刻，将所有object的信息，以dump的形式输出
void output_dump(const std::vector<Object *> &_obj_ptr_list, const Setting &_set, bool _initialize);

//在固定时刻，在vac.txt, sia.txt, fia.txt, vac_fia.txt, sia_fia.txt中输出相应Object的尺寸信息
void output_txt(const std::vector<Object *> &_obj_ptr_list, const Setting &_set, bool _initialize);

//输出某一帧的cascade
void output_cascade(const std::vector<Object *> &_obj_ptr_list, const Setting &_set);

void write_content(std::string _filename, std::ostringstream &_oss, bool);

//cascade insert
void output_react_0_0(const Setting &_set, std::string str);

//obj reach bound in absorb boundary condition
void output_react_1_0(const Object &_obj1, const Setting &_set);

//obj被晶界吸收
void output_react_1_0_GB(const Object &_obj1, const Setting &_set);

//obj被位错吸收
void output_react_1_0_dislocation(const Object &_obj1, const Setting &_set);

//两个Object完全湮灭，产生0个新object
//输入两个object，将反应信息写入"react.txt"
void output_react_2_0(const Object &_obj1, const Object &_obj2, const Setting &_set);

// obj1 + obj2 = obj3
void output_react_2_1(const Object &_obj1, const Object &_obj2, const Object &_obj3, const Setting &_set);

//obj1 -> obj2 + obj3
void output_react_1_2(const Object &_obj1, const Object &_obj2, const Object &_obj3, const Setting &_set);

//TM
void output_react_TM(const Object& _obj1, const Object& _obj2, const Object& _obj3, const Setting& _set);

inline double cal_duration();


// 计算程序的cpu时间
double cal_duration()
{
	return (double)clock() / CLOCKS_PER_SEC;
}

#endif // !OUTPUT_H_
