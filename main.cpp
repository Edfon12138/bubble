#include<cstdlib>
#include<ctime>
#include<iostream>
#include<fstream>
#include<iomanip>

#include"setting.h"
#include"database.h"
#include"object_parameter.h"
#include"object.h"
#include"output.h"
#include"event.h"
#include"random.h"
#include"reaction.h"

/*
�β�ȫ����_��ͷ
�������δ������ã������������Ϳ��Դ��ݿ�����
*/

//��ʼ��ѭ��֮ǰ��ȫ����������Ƿ�����ڽ��
//����cascade֮��ȫ����������Ƿ�����ڽ��
//������������ͨ��ϣ�ֻ��鶯���Ǹ�obj��

int main()
{	
	//��ʼ������ļ�
	initialize_react_file();

	//��ȡ���ã����������á����󣬳�ʼ�������
	Setting settings("input_setting.txt");

	//���input_setting��û��ָ����������ӣ�����ϵͳ��ǰʱ��Ϊ���������
	if (settings.seed < 0)
	{
		settings.seed = (long)time(0);
	}
	std::cout << "The seed is " << settings.seed << std::endl;
	suni(settings.seed);


	//��ȡobject����������object�������ݿ�
	Database *ptr_vac_database = new Database("database_vac.txt");
	Database *ptr_sia_database = new Database("database_sia.txt");
	Database *ptr_fia_database = new Database("database_fia.txt");
	Database *ptr_vac_fia_database = new Database("database_vac_fia.txt");
	Database *ptr_sia_fia_database = new Database("database_sia_fia.txt");

	//��ȡobjectλ�ã���obj_ptr_list�洢ָ��objectָ��
	//������������֮��õ�ָ�룬�ٴν����ò��ܵõ�Object����
	std::vector<Object *> obj_ptr_list;	 	
	read_cascade("input_cascade.txt", obj_ptr_list, *ptr_vac_database, *ptr_sia_database, *ptr_fia_database,
		*ptr_vac_fia_database, *ptr_sia_fia_database, settings);


	//����ʼ������
	int pos1 = 0, pos2 = 0;
	while (check_dist_cross_bound(obj_ptr_list, pos1, pos2, settings))
	{
		carry_out_reaction(obj_ptr_list, pos1, pos2, *obj_ptr_list.at(pos1), *obj_ptr_list.at(pos2),
			*ptr_vac_database, *ptr_sia_database, *ptr_fia_database, *ptr_vac_fia_database, *ptr_sia_fia_database, settings);
		pos2 = 0;
	}


	//��ʼ��ѭ��֮ǰ���һ��
	output_dump(obj_ptr_list, settings, true);
	output_txt(obj_ptr_list, settings, true);
	output_cascade(obj_ptr_list, settings);


	//��ѭ��
	std::cout << "Start main loop!\n" << std::endl;

	bool flag_refresh_event_list = true;
	double rate_sum;
	std::vector<Event> event_list;
	while (settings.check_end())
	{
		//���obj_ptr_list�����仯�����¹����¼��б����Ҹ���rate_sum
		if (flag_refresh_event_list) {
			rate_sum = build_event_list(obj_ptr_list, event_list,settings);
			flag_refresh_event_list = false;
		}

		//����¼��б�Ϊ�գ�rate_sum��Ϊ0�����ȡ�¼���ִ���¼����ۼ�ʱ��Ͳ���
		if (rate_sum > 0.0)
		{
			//�����ȡ�¼�
			//���е�event��event_list�е�λ��
			int event_pos = find_event(event_list, rate_sum);

			//����鵽����cascade�����¼�
			if (event_list.at(event_pos).which_event == -1)
			{
				//����cascade
				std::string cascade_filename = choose_cascade_to_insert(settings);
				if (settings.output_cascade_injection > 0) {
					output_react_0_0(settings, cascade_filename);
				}
				int pos1 = 0, pos2 = 0;
				read_cascade(cascade_filename, obj_ptr_list, *ptr_vac_database, *ptr_sia_database,
					*ptr_fia_database, *ptr_vac_fia_database, *ptr_sia_fia_database, settings);
				while (check_dist_cross_bound(obj_ptr_list, pos1, pos2, settings))
				{
					carry_out_reaction(obj_ptr_list, pos1, pos2, *obj_ptr_list.at(pos1), *obj_ptr_list.at(pos2),
						*ptr_vac_database, *ptr_sia_database, *ptr_fia_database, *ptr_vac_fia_database, *ptr_sia_fia_database,
						settings);
					pos2 = 0;
				}
				flag_refresh_event_list = true;
			}
			//��������鵽����object���¼�
			else
			{
				//ִ���¼�
				//���ܵ��¼���Ǩ�ơ���ת�����䡢cascade����
				//�������Ǩ�ƻ��߷��䣬�����Ƿ�ײ�ϱ��obj
				int obj_pos = event_list.at(event_pos).pos;	//��obj��obj_ptr_list�е�λ��
				int which_event = event_list.at(event_pos).which_event;
				carry_out_event(which_event, obj_pos, obj_ptr_list, flag_refresh_event_list,
					*ptr_vac_database, *ptr_sia_database, *ptr_fia_database, *ptr_vac_fia_database, *ptr_sia_fia_database, settings);

			}
			settings.increase_time_step(rate_sum);
		}
		//����¼��б�Ϊ�գ��������ܣ�ÿ��ʱ���ۼ�default_dt�������ۼ�1
		else {
			settings.time += settings.default_dt;
			settings.step += 1;
		}

		//output_dump(obj_ptr_list, settings, false);

		//����Ƿ���Ҫ���
		if (settings.check_output())
		{
			output_dump(obj_ptr_list, settings, false);
			output_txt(obj_ptr_list, settings, false);
			output_cascade(obj_ptr_list, settings);
		}
	}

	//�����ʱ��Ϊ������׼�������һ��ִ����֮��ʱ�䳬��max_time�����ʱ������max_time
	if (settings.stop_criteria == 1 && settings.time > settings.max_time) {
		settings.time = settings.max_time;
	}

	//��ֹ��ѭ��֮�������һ��
	output_dump(obj_ptr_list, settings, false);
	output_txt(obj_ptr_list, settings, false);
	output_cascade(obj_ptr_list, settings);


	//���Դ���

	//�ͷ�����obj���ڴ�
	for (auto iter : obj_ptr_list) {
		delete iter;
	}
	//�ͷ�����database���ڴ�
	delete ptr_vac_database;
	delete ptr_sia_database;
	delete ptr_fia_database;
	delete ptr_vac_fia_database;
	delete ptr_sia_fia_database;

	return 0;
}