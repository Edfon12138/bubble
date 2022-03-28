#pragma once


#ifndef OBJECT_H_
#define OBJECT_H_

#include<vector>
#include <array>
#include<iostream>
#include<algorithm>
#include<cmath>
#include<set>

#include"database.h"
#include"setting.h"
#include"random.h"


#define PI 3.1415926

class Object
{
private:
public:
	// members
	static int ID;
	std::array<double, 3> pos;	//实时坐标
	std::array<double, 3> displacement;	//位移

	//记录object在哪个网格
	std::array<int, 3> grid_pos;

	//Record how many times the object cross the bound
	std::array<int, 3> cross_bound;	

	//记录各向异性迁移object穿边界历史信息的数组
	std::vector<std::array<double, 3>> cross_bound_information;


	//no direction~0;	111~1;	-111~2;	1-11~3;	11-1~4;
	int dir;	

	//type: vac~1;	sia~2;	fia~3;	vac_fia~4;	sia_fia~5
	int type;

	//包含空位，即为1；包含自间隙，即为2；都不包含则为0
	int type_plus;

	size_t id;
	int size1, size2;
	//指前因子和对应的能垒
	double mig_F, mig_E, emit_F1, emit_E1, emit_F2, emit_E2, rotate_F, rotate_E, TM_F, TM_E;

	//根据指前因子、能垒、温度、kB算出来的频率
	double mig_Freq, emit_Freq1, emit_Freq2, rotate_Freq, TM_Freq;
	double radius;

	//记录位于tis的object，哪一个坐标是整数，哪一个坐标是半整数，哪一个坐标是1/4整数
	std::array<int, 3> tis_quarter_half_int_pos;
	std::map<int, double> frequency;	//mig~0, emit1~1, emit2~2, rotate~3 and Trap mutation~4

	//被虚拟位错吸收的概率
	double p_meet_dislocation;

	//functions
	Object() { std::cout << "The object is contructed in dafault style.\n"; }

	Object(int _type, double _x, double _y, double _z, int _size1, int _size2, int _dir,
		const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
		const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set);

	~Object() {};

	//根据Object的size1和size2，对Object进行参数化
	void parameterization(const Object_parameter & _parameter);

	//先将Object平移到原点，寻找离他最近的bcc格点，再将他平移回原来的位置附近
	void on_bcc_lattice(const Setting &_set);
	//先将Object平移到原点，寻找离他最近的tis格点，再将他平移回原来的位置附近
	void on_tis_site(const Setting &_set);

	//根据object的类型、半径，位错的密度、半径，计算object被位错吸收的概率
	void cal_p_meet_dislocation(const Setting &_set);

	//计算object具体在哪个网格
	inline void refresh_grid_pos(const Setting &_set);

	// 初始生成Object之后执行；或者pbc条件下，每次object迁移之后执行
	// 如果越过边界，则将object平移到盒子内
	inline void check_pbc(const Setting &_set);

	//检查object是否越界，从-x越界返回1；-y返回2；-z返回3；+z返回4；+y返回5；+x返回6；否则返回0
	inline int check_cross_bound_xyz(const Setting &_set);

	//如果z越界，即返回true，否则返回false
	inline bool check_cross_bound_z(const Setting &_set);

	//检查obj是否到达预设的晶粒边界
	inline bool check_grain_radius(const Setting &_set);

	//自适应搜寻可能发生的事件，记录在frequency map中，后续进入event list
	void cal_frequency(const Setting &_set);
	
	//该Obj执行一步扩散
	//fia暂时用vac的三维扩散方式
	//vac_fia用vac的三维扩散方式；sia_fia用sia的一维扩散方式
	inline void mig(const Setting &_set);

	// 自间隙迁移
	inline void mig_sia(const Setting &_set);

	// 空位迁移
	inline void mig_vac(const Setting &_set);

	//fia迁移，从tis到tis
	inline void mig_fia_tis_to_tis(const Setting &_set);

	// 自间隙转向，从剩余的3个方向中随机任选一个方向，作为sia的新方向
	inline void rotate_sia();

	inline void swap(int &a, int &b);



};


/****************************普通函数声明****************************/

//输入两个原子的三维坐标，返回它们之间的距离
double cal_dist(std::array<double, 3> a, std::array<double, 3> b);

// 读cascade.txt，每一行生成一个object（动态分配内存），将指向它的指针存在obj_ptr_list中
void read_cascade(std::string _filename, std::vector<Object *> &_obj_ptr_list,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set);

//旋转矩阵左乘空间坐标
void multiply(std::array<std::array<double, 3>, 3> & Rx, double & x, double & y, double & z);

/****************************内联函数定义****************************/

//生成object，或者object迁移并检查完边界条件之后执行，更新object所在的网格
void Object::refresh_grid_pos(const Setting &_set)
{
	grid_pos.at(0) = int(pos[0]) / _set.box_grid;
	grid_pos.at(1) = int(pos[1]) / _set.box_grid;
	grid_pos.at(2) = int(pos[2]) / _set.box_grid;
}

//检查object是否越界，从-x越界返回1；-y返回2；-z返回3；+z返回4；+y返回5；+x返回6；否则返回0
int Object::check_cross_bound_xyz(const Setting &_set)
{
	if (pos[0] < _set.box_min[0])	return 1;
	if (pos[0] > _set.box_max[0])	return 2;
	if (pos[1] < _set.box_min[1])	return 3;
	if (pos[1] > _set.box_max[1])	return 4;
	if (pos[2] < _set.box_min[2])	return 5;
	if (pos[2] > _set.box_max[2])	return 6;
	return 0;
}

//如果x,y,z有一个越界，即返回true，否则返回false
/*
int Object::check_cross_bound_xyz(const Setting &_set)
{
	for (int i = 0; i < 3; i++)
	{
		if (pos.at(i) < _set.box_min.at(i) || pos.at(i) > _set.box_max.at(i))
		{
			return true;
		}
	}
	return false;
}
*/

// 初始生成Object之后执行；或者pbc条件下，每次object迁移之后执行
// 如果越过边界，则将object平移到盒子内
void Object::check_pbc(const Setting &_set)
{
	for (int i = 0; i < 3; i++)
	{
		while (pos.at(i) < _set.box_min.at(i))
		{
			pos.at(i) += _set.box_length.at(i);
			cross_bound.at(i)--;
		}
		while (pos.at(i) > _set.box_max.at(i))
		{
			pos.at(i) -= _set.box_length.at(i);
			cross_bound.at(i)++;
		}
	}
}

bool Object::check_cross_bound_z(const Setting &_set)
{
	if (pos.at(2) < _set.box_min.at(2) || pos.at(2) > _set.box_max.at(2))
	{
		return true;
	}
	else return false;
}

bool Object::check_grain_radius(const Setting &_set)
{
	double r2 = pow(displacement.at(0), 2) + pow(displacement.at(1), 2) + pow(displacement.at(2), 2);
	if (r2 > _set.grain_radius_square) return true;
	else return false;
}

//vac_fia用vac的三维扩散方式；sia_fia用sia的一维扩散方式
void Object::mig(const Setting &_set)
{
	switch (type)
	{
	case 1:
		mig_vac(_set);
		break;
	case 2:
		mig_sia(_set);
		break;
	case 3:
		if (_set.fia_site == 1)
			mig_fia_tis_to_tis(_set);
		else if (_set.fia_site == 0)
			mig_vac(_set);
		break;
	case 4:
		mig_vac(_set);
		break;
	case 5:
		mig_sia(_set);
		break;
	}
}

// 自间隙迁移
void Object::mig_sia(const Setting &_set)
{
	int rand_num = 0;
	for (int i = 0; i < _set.mig_degeneration; i++)
	{
		rand_num += (uni() < 0.5) ? 1 : -1;
	}
	switch (dir)
	{
	case 1:
		pos.at(0) += 0.5 * rand_num;
		pos.at(1) += 0.5 * rand_num;
		pos.at(2) += 0.5 * rand_num;
		displacement.at(0) += 0.5 * rand_num;
		displacement.at(1) += 0.5 * rand_num;
		displacement.at(2) += 0.5 * rand_num;
		break;
	case 2:
		pos.at(0) -= 0.5 * rand_num;
		pos.at(1) += 0.5 * rand_num;
		pos.at(2) += 0.5 * rand_num;
		displacement.at(0) -= 0.5 * rand_num;
		displacement.at(1) += 0.5 * rand_num;
		displacement.at(2) += 0.5 * rand_num;
		break;
	case 3:
		pos.at(0) += 0.5 * rand_num;
		pos.at(1) -= 0.5 * rand_num;
		pos.at(2) += 0.5 * rand_num;
		displacement.at(0) += 0.5 * rand_num;
		displacement.at(1) -= 0.5 * rand_num;
		displacement.at(2) += 0.5 * rand_num;
		break;
	case 4:
		pos.at(0) += 0.5 * rand_num;
		pos.at(1) += 0.5 * rand_num;
		pos.at(2) -= 0.5 * rand_num;
		displacement.at(0) += 0.5 * rand_num;
		displacement.at(1) += 0.5 * rand_num;
		displacement.at(2) -= 0.5 * rand_num;
		break;
	default:
		std::cout << "Fail to choose a direction for sia (id is " << id << ") migration.\n";
	}
}

// 空位迁移
void Object::mig_vac(const Setting &_set)
{
	std::array<int, 3> rand_num = { 0,0,0 };
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < _set.mig_degeneration; j++)
		{
			rand_num.at(i) += (uni() < 0.5) ? 1 : -1;
		}
	}
	for (int i = 0; i < 3; i++)
	{
		pos.at(i) += rand_num.at(i) * 0.5;
		displacement.at(i) += rand_num.at(i) * 0.5;
	}
}

void Object::mig_fia_tis_to_tis(const Setting &_set)
{
	//迁移mig_degeneration次
	for (int j = 0; j < _set.mig_degeneration; j++)
	{
		double d1, d2;
		d1 = (uni() < 0.5) ? 0.25 : -0.25;
		d2 = (uni() < 0.5) ? 0.25 : -0.25;
		double changed_quarter = pos.at(tis_quarter_half_int_pos.at(0)) + d1;
		if ((changed_quarter - floor(changed_quarter)) < 0.1)	//如果quarter变成了int
		{
			pos.at(tis_quarter_half_int_pos.at(0)) += d1;
			pos.at(tis_quarter_half_int_pos.at(2)) += d2;	//则int变quarter
			swap(tis_quarter_half_int_pos.at(0), tis_quarter_half_int_pos.at(2));	//交换int和quarter
		}
		else
			//如果quarter变成了half
		{
			pos.at(tis_quarter_half_int_pos.at(0)) += d1;
			pos.at(tis_quarter_half_int_pos.at(1)) += d2;	//则half变quarter
			swap(tis_quarter_half_int_pos.at(0), tis_quarter_half_int_pos.at(1));	//交换half和quarter
		}
	}
}

void Object::swap(int &a, int &b)
{
	int temp;
	temp = a;
	a = b;
	b = temp;
}

// 自间隙转向，从剩余的3个方向中随机任选一个方向，作为sia的新方向
void Object::rotate_sia()
{
	double rand_num = uni();
	std::array<int, 3> possible_dir = { 0,0,0 };
	switch (dir)
	{
	case 1:
		possible_dir.at(0) = 2;
		possible_dir.at(1) = 3;
		possible_dir.at(2) = 4;
		break;
	case 2:
		possible_dir.at(0) = 3;
		possible_dir.at(1) = 4;
		possible_dir.at(2) = 1;
		break;
	case 3:
		possible_dir.at(0) = 4;
		possible_dir.at(1) = 1;
		possible_dir.at(2) = 2;
		break;
	case 4:
		possible_dir.at(0) = 1;
		possible_dir.at(1) = 2;
		possible_dir.at(2) = 3;
		break;
	}
	if (rand_num < double(1) / 3) { dir = possible_dir.at(0); }
	else if (rand_num < double(2) / 3) { dir = possible_dir.at(1); }
	else { dir = possible_dir.at(2); }
}

#endif // !OBJECT_H_
