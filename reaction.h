#pragma once

#ifndef REACTION_H_
#define REACTION_H_

#include"object.h"
#include"output.h"	//发生reaction之后，需要调用Output中的函数，更新"react.txt"


/***************************************函数声明区***************************************/

//适用于检查从pos1之后的obj（含pos1指向的obj）是否和obj_list中其他所有obj的结合情况
//输入vector和两个object在vector中的位置pos1和pos2
//对整个vector，外循环pos1，内循环pos2（通过check_dist_inner_loop函数实现）
//比较两个Object之间的距离，如果所有object都不挨着，返回false
//如果监测到任意两个object挨着，立即终止循环，返回true，且返回两个object在vector中的最新位置pos1和pos2（通过引用实现）
//只有当pos1越界时才会返回false，从而确保所有的Object都两两比较
bool check_dist(std::vector<Object *> &_obj_ptr_list, int &pos1, int &pos2);


//pos1是不动的，只动Pos2，所以Pos2要传引用，pos1不需要
//pos2从初始位置开始，依此递增，如果Pos2不等于pos1，则比较它们指向的Object之间的距离
//如果距离大于半径之和，则无事发生；否则，返回true
bool check_dist_inner_loop(const std::vector<Object *> &_obj_ptr_list, const int &_pos1, int &_pos2);

bool check_dist_inner_loop_abc(const std::vector<Object *> &_obj_ptr_list, const int &_pos1, int &_pos2, const Setting &_set);



//check_dist函数的升级版
//判断跨周期结合
bool check_dist_cross_bound(std::vector<Object *> &_obj_ptr_list, int &pos1, int &pos2, const Setting &_set);


//check_dist_inner_loop函数的升级版
//判断跨周期结合
bool check_dist_cross_bound_inner_loop(const std::vector<Object *> &_obj_ptr_list, const int &_pos1, int &_pos2, 
	const Setting &_set);



//输入pos1和pos2指向的两个object
//如果不是完全湮灭，则动态生成一个新的oject,并将指针压入obj_ptr_list
//删除两个反应object，释放它们的内存，并将指向它们的指针从vector中去除
//该函数不会改变pos1,pos2的值
void carry_out_reaction(std::vector<Object *> &_obj_ptr_list, int &_pos1, const int &_pos2, Object &_obj1, Object &_obj2,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set);


// 输入两个object的size和type1，返回新产生的object的类型
// -1代表生成了没考虑到的类型；0代表完全湮灭，不产生新object
int judge_new_obj_type(int _type1, int _obj1_s1, int _type2, int _obj2_s1);


//根据type, x, y, z, size1, size2, dir，创建一个新的object
void create_object(int _type, Object &_obj1, Object &_obj2, std::vector<Object *> &_obj_ptr_list,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set);


//计算新object的位置
//应该用size加权平均还是用radius加权平均？？目前用total_size 加权平均
void generate_new_pos(Object &_obj1, Object &_obj2, std::array<double, 3> &_new_pos);


#endif // !REACTION_H_
