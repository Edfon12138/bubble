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
	std::array<double, 3> pos;	//ʵʱ����
	std::array<double, 3> displacement;	//λ��

	//��¼object���ĸ�����
	std::array<int, 3> grid_pos;

	//Record how many times the object cross the bound
	std::array<int, 3> cross_bound;	

	//��¼��������Ǩ��object���߽���ʷ��Ϣ������
	std::vector<std::array<double, 3>> cross_bound_information;


	//no direction~0;	111~1;	-111~2;	1-11~3;	11-1~4;
	int dir;	

	//type: vac~1;	sia~2;	fia~3;	vac_fia~4;	sia_fia~5
	int type;

	//������λ����Ϊ1�������Լ�϶����Ϊ2������������Ϊ0
	int type_plus;

	size_t id;
	int size1, size2;
	//ָǰ���ӺͶ�Ӧ������
	double mig_F, mig_E, emit_F1, emit_E1, emit_F2, emit_E2, rotate_F, rotate_E, TM_F, TM_E;

	//����ָǰ���ӡ����ݡ��¶ȡ�kB�������Ƶ��
	double mig_Freq, emit_Freq1, emit_Freq2, rotate_Freq, TM_Freq;
	double radius;

	//��¼λ��tis��object����һ����������������һ�������ǰ���������һ��������1/4����
	std::array<int, 3> tis_quarter_half_int_pos;
	std::map<int, double> frequency;	//mig~0, emit1~1, emit2~2, rotate~3 and Trap mutation~4

	//������λ�����յĸ���
	double p_meet_dislocation;

	//functions
	Object() { std::cout << "The object is contructed in dafault style.\n"; }

	Object(int _type, double _x, double _y, double _z, int _size1, int _size2, int _dir,
		const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
		const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set);

	~Object() {};

	//����Object��size1��size2����Object���в�����
	void parameterization(const Object_parameter & _parameter);

	//�Ƚ�Objectƽ�Ƶ�ԭ�㣬Ѱ�����������bcc��㣬�ٽ���ƽ�ƻ�ԭ����λ�ø���
	void on_bcc_lattice(const Setting &_set);
	//�Ƚ�Objectƽ�Ƶ�ԭ�㣬Ѱ�����������tis��㣬�ٽ���ƽ�ƻ�ԭ����λ�ø���
	void on_tis_site(const Setting &_set);

	//����object�����͡��뾶��λ����ܶȡ��뾶������object��λ�����յĸ���
	void cal_p_meet_dislocation(const Setting &_set);

	//����object�������ĸ�����
	inline void refresh_grid_pos(const Setting &_set);

	// ��ʼ����Object֮��ִ�У�����pbc�����£�ÿ��objectǨ��֮��ִ��
	// ���Խ���߽磬��objectƽ�Ƶ�������
	inline void check_pbc(const Setting &_set);

	//���object�Ƿ�Խ�磬��-xԽ�緵��1��-y����2��-z����3��+z����4��+y����5��+x����6�����򷵻�0
	inline int check_cross_bound_xyz(const Setting &_set);

	//���zԽ�磬������true�����򷵻�false
	inline bool check_cross_bound_z(const Setting &_set);

	//���obj�Ƿ񵽴�Ԥ��ľ����߽�
	inline bool check_grain_radius(const Setting &_set);

	//����Ӧ��Ѱ���ܷ������¼�����¼��frequency map�У���������event list
	void cal_frequency(const Setting &_set);
	
	//��Objִ��һ����ɢ
	//fia��ʱ��vac����ά��ɢ��ʽ
	//vac_fia��vac����ά��ɢ��ʽ��sia_fia��sia��һά��ɢ��ʽ
	inline void mig(const Setting &_set);

	// �Լ�϶Ǩ��
	inline void mig_sia(const Setting &_set);

	// ��λǨ��
	inline void mig_vac(const Setting &_set);

	//fiaǨ�ƣ���tis��tis
	inline void mig_fia_tis_to_tis(const Setting &_set);

	// �Լ�϶ת�򣬴�ʣ���3�������������ѡһ��������Ϊsia���·���
	inline void rotate_sia();

	inline void swap(int &a, int &b);



};


/****************************��ͨ��������****************************/

//��������ԭ�ӵ���ά���꣬��������֮��ľ���
double cal_dist(std::array<double, 3> a, std::array<double, 3> b);

// ��cascade.txt��ÿһ������һ��object����̬�����ڴ棩����ָ������ָ�����obj_ptr_list��
void read_cascade(std::string _filename, std::vector<Object *> &_obj_ptr_list,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set);

//��ת������˿ռ�����
void multiply(std::array<std::array<double, 3>, 3> & Rx, double & x, double & y, double & z);

/****************************������������****************************/

//����object������objectǨ�Ʋ������߽�����֮��ִ�У�����object���ڵ�����
void Object::refresh_grid_pos(const Setting &_set)
{
	grid_pos.at(0) = int(pos[0]) / _set.box_grid;
	grid_pos.at(1) = int(pos[1]) / _set.box_grid;
	grid_pos.at(2) = int(pos[2]) / _set.box_grid;
}

//���object�Ƿ�Խ�磬��-xԽ�緵��1��-y����2��-z����3��+z����4��+y����5��+x����6�����򷵻�0
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

//���x,y,z��һ��Խ�磬������true�����򷵻�false
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

// ��ʼ����Object֮��ִ�У�����pbc�����£�ÿ��objectǨ��֮��ִ��
// ���Խ���߽磬��objectƽ�Ƶ�������
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

//vac_fia��vac����ά��ɢ��ʽ��sia_fia��sia��һά��ɢ��ʽ
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

// �Լ�϶Ǩ��
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

// ��λǨ��
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
	//Ǩ��mig_degeneration��
	for (int j = 0; j < _set.mig_degeneration; j++)
	{
		double d1, d2;
		d1 = (uni() < 0.5) ? 0.25 : -0.25;
		d2 = (uni() < 0.5) ? 0.25 : -0.25;
		double changed_quarter = pos.at(tis_quarter_half_int_pos.at(0)) + d1;
		if ((changed_quarter - floor(changed_quarter)) < 0.1)	//���quarter�����int
		{
			pos.at(tis_quarter_half_int_pos.at(0)) += d1;
			pos.at(tis_quarter_half_int_pos.at(2)) += d2;	//��int��quarter
			swap(tis_quarter_half_int_pos.at(0), tis_quarter_half_int_pos.at(2));	//����int��quarter
		}
		else
			//���quarter�����half
		{
			pos.at(tis_quarter_half_int_pos.at(0)) += d1;
			pos.at(tis_quarter_half_int_pos.at(1)) += d2;	//��half��quarter
			swap(tis_quarter_half_int_pos.at(0), tis_quarter_half_int_pos.at(1));	//����half��quarter
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

// �Լ�϶ת�򣬴�ʣ���3�������������ѡһ��������Ϊsia���·���
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
