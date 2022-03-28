#pragma once


#ifndef SETTING_H_
#define SETTING_H_

#include<string>
#include<fstream>
#include<sstream>
#include<iostream>
#include<map>
#include<cmath>
#include<array>
#include<vector>
#include"random.h"

/*
主要功能：
读取setting.txt里的所有参数，并且记录体系当前的time,step
判断主进程是否终止
时间和步数实现递增
*/
// 存储setting.txt中所有参数，以及当前时间和步数的类
class Setting
{
private:
public:
	/***************************************members***************************************/
	//进度条参数
	int progress_percentage;

	//随机数种子
	long seed;

	//物理参数：温度，晶格常数，玻尔兹曼常数，普朗克常量
	double temperature;
	double a0, kb, hp;

	//盒子相关
	std::array<double, 3> box_min;	//单位为晶格常数
	std::array<double, 3> box_max;	//单位为晶格常数
	std::array<double, 3> box_length;	//单位为晶格常数

	int boundary_condition;	//1~absorb;	2~pbc; 3~z free
	int box_grid;	//每隔多少个晶格常数划分一个网格
	std::array<int, 3> box_grid_num;	//box在每一边被划分的网格数减1

	//最大步数，最大时间，终止判据
	unsigned long long max_step;
	double max_time;
	int stop_criteria;	//1~time	2~loop

	//当前时间、步数
	unsigned long long step;
	double time;
	
	//object相关

	int fia_site;	//0~bcc, 1~tis, 2~ois
	// 9 sites in a bcc unit cell.
	std::array<std::array<double, 3>, 9> bcc_lattice_site;
	std::array<std::array<double, 3>, 24> tis_site;

	double He_flux;

	//cascade入射
	double rate_cascade;
	double cascade_interval;
	double cascade_insert_time;
	std::vector<double> cascade_energy, cascade_possibility;
	std::vector<int> cascade_number;
	int rand_translate_cascade;	//control whether to randomly translate the inserted cascade (0~false, 1~true)
	int rand_rotate_cascade;	//(0~false1, ~true)


	//特殊参数
	double max_energy;
	int mig_degeneration;
	int grain_radius;
	int grain_radius_square;
	double default_dt;

	//sink strength相关参数
	double ss_grain_radius;
	double dislocation_density;
	double dislocation_radius;
	double p_sia_meet_GB;
	double p_vac_meet_GB;




	//输出
	int total_output;	//总共需要输出的次数
	int n_output;		//当前已经输出了多少次
	//int output_c_a_e_react;	//是否输出结合、湮灭和发射反应，1是输出，-1是不输出
	int output_cascade_injection;
	int output_combine_reaction;
	int output_annihilation_reaction;
	int output_emit_reaction;
	int output_out_bound;
	int output_GB_absorption;
	int output_disloaction_absorption;
	int output_trap_mutation;

	//临时存储setting.txt中的值
	std::map<std::string, double> setting_map;

	/***************************************functions***************************************/

	Setting() { std::cout<<"The defaulted construct function is invoked.\n"; }

	// 输入文件名，调用read_setting函数，生成setting_map，并用setting_map初始化得到一个setting对象
	Setting(std::string filename);

	~Setting(){}

	// 输入文件名，读取文件中所有非空，非#开头的行，以第一个关键字为键，第二个关键字为值，构造setting map
	void read_setting(std::string filename);

	void cal_sink_strength();

	void cal_cascade_insert_rate();

	//读取cascade_possibility.txt
	void read_cascade_possibility();

	// 如果step或者time小于截止值，返回true，继续循环
	inline bool check_end();

	//检查是否需要输出
	inline bool check_output();

	// 执行一步后，根据total_rate，增加系统时间，step++
	inline void increase_time_step(double _totalRate);
};


//definitions

// 如果step或者time小于截止值，返回true，继续循环
bool Setting::check_end()
{
	//根据时间停止
	if (stop_criteria == 1)
	{
		if (time > max_time / 100 * progress_percentage) {
			std::cout << progress_percentage << "%completed\n";
			progress_percentage++;
		}
		if (time < max_time) { return true; }
	}
	//根据步数停止
	else if (stop_criteria == 2)
	{
		if (step > max_step / 100.0 * progress_percentage) {
			std::cout << progress_percentage << "%completed\n";
			progress_percentage++;
		}

		if (step < max_step) { return true; }
	}
	else { 
		std::cout << "stop_criteria has a problem.\n"; }
	return false;
}

// 如果step或者time大于输出节点，返回true，准备输出
bool Setting::check_output()
{
	//根据时间停止
	if (stop_criteria == 1)
	{
		if (time > max_time / total_output * n_output){ 
			n_output++;
			return true; 
		}
		else {
			return false;
		}
	}
	//根据步数停止
	else if (stop_criteria == 2)
	{
		// 对于两个整数，'/'表示整除
		if (step > (max_step / total_output * n_output - 1) )
		{
			n_output++;
			return true; 
		}
		else return false;
	}
	else {
		std::cout << "stop_criteria has a problem.\n";
		return false;
	}
}

// 执行一步后，根据total_rate，增加系统时间，并且step++
void Setting::increase_time_step(double _total_rate)
{
	double rand_num = uni();
	while (rand_num < 1e-20)
	{
		rand_num = uni();
	}
	double delta_time = -log(rand_num) / _total_rate;
	time += delta_time;
	step += 1;
}

#endif // !SETTING_H_
