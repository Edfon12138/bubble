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
��Ҫ���ܣ�
��ȡsetting.txt������в��������Ҽ�¼��ϵ��ǰ��time,step
�ж��������Ƿ���ֹ
ʱ��Ͳ���ʵ�ֵ���
*/
// �洢setting.txt�����в������Լ���ǰʱ��Ͳ�������
class Setting
{
private:
public:
	/***************************************members***************************************/
	//����������
	int progress_percentage;

	//���������
	long seed;

	//����������¶ȣ��������������������������ʿ˳���
	double temperature;
	double a0, kb, hp;

	//�������
	std::array<double, 3> box_min;	//��λΪ������
	std::array<double, 3> box_max;	//��λΪ������
	std::array<double, 3> box_length;	//��λΪ������

	int boundary_condition;	//1~absorb;	2~pbc; 3~z free
	int box_grid;	//ÿ�����ٸ�����������һ������
	std::array<int, 3> box_grid_num;	//box��ÿһ�߱����ֵ���������1

	//����������ʱ�䣬��ֹ�о�
	unsigned long long max_step;
	double max_time;
	int stop_criteria;	//1~time	2~loop

	//��ǰʱ�䡢����
	unsigned long long step;
	double time;
	
	//object���

	int fia_site;	//0~bcc, 1~tis, 2~ois
	// 9 sites in a bcc unit cell.
	std::array<std::array<double, 3>, 9> bcc_lattice_site;
	std::array<std::array<double, 3>, 24> tis_site;

	double He_flux;

	//cascade����
	double rate_cascade;
	double cascade_interval;
	double cascade_insert_time;
	std::vector<double> cascade_energy, cascade_possibility;
	std::vector<int> cascade_number;
	int rand_translate_cascade;	//control whether to randomly translate the inserted cascade (0~false, 1~true)
	int rand_rotate_cascade;	//(0~false1, ~true)


	//�������
	double max_energy;
	int mig_degeneration;
	int grain_radius;
	int grain_radius_square;
	double default_dt;

	//sink strength��ز���
	double ss_grain_radius;
	double dislocation_density;
	double dislocation_radius;
	double p_sia_meet_GB;
	double p_vac_meet_GB;




	//���
	int total_output;	//�ܹ���Ҫ����Ĵ���
	int n_output;		//��ǰ�Ѿ�����˶��ٴ�
	//int output_c_a_e_react;	//�Ƿ������ϡ�����ͷ��䷴Ӧ��1�������-1�ǲ����
	int output_cascade_injection;
	int output_combine_reaction;
	int output_annihilation_reaction;
	int output_emit_reaction;
	int output_out_bound;
	int output_GB_absorption;
	int output_disloaction_absorption;
	int output_trap_mutation;

	//��ʱ�洢setting.txt�е�ֵ
	std::map<std::string, double> setting_map;

	/***************************************functions***************************************/

	Setting() { std::cout<<"The defaulted construct function is invoked.\n"; }

	// �����ļ���������read_setting����������setting_map������setting_map��ʼ���õ�һ��setting����
	Setting(std::string filename);

	~Setting(){}

	// �����ļ�������ȡ�ļ������зǿգ���#��ͷ���У��Ե�һ���ؼ���Ϊ�����ڶ����ؼ���Ϊֵ������setting map
	void read_setting(std::string filename);

	void cal_sink_strength();

	void cal_cascade_insert_rate();

	//��ȡcascade_possibility.txt
	void read_cascade_possibility();

	// ���step����timeС�ڽ�ֵֹ������true������ѭ��
	inline bool check_end();

	//����Ƿ���Ҫ���
	inline bool check_output();

	// ִ��һ���󣬸���total_rate������ϵͳʱ�䣬step++
	inline void increase_time_step(double _totalRate);
};


//definitions

// ���step����timeС�ڽ�ֵֹ������true������ѭ��
bool Setting::check_end()
{
	//����ʱ��ֹͣ
	if (stop_criteria == 1)
	{
		if (time > max_time / 100 * progress_percentage) {
			std::cout << progress_percentage << "%completed\n";
			progress_percentage++;
		}
		if (time < max_time) { return true; }
	}
	//���ݲ���ֹͣ
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

// ���step����time��������ڵ㣬����true��׼�����
bool Setting::check_output()
{
	//����ʱ��ֹͣ
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
	//���ݲ���ֹͣ
	else if (stop_criteria == 2)
	{
		// ��������������'/'��ʾ����
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

// ִ��һ���󣬸���total_rate������ϵͳʱ�䣬����step++
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
