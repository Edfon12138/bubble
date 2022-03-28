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
形参全部用_打头
函数传参传递引用（基础数据类型可以传递拷贝）
*/

//开始主循环之前，全局两两检查是否跨周期结合
//插入cascade之后，全局两两检查是否跨周期结合
//其他情况检查普通结合（只检查动的那个obj）

int main()
{	
	//初始化输出文件
	initialize_react_file();

	//读取设置，建立“设置”对象，初始化随机数
	Setting settings("input_setting.txt");

	//如果input_setting中没有指定随机数种子，则以系统当前时间为随机数种子
	if (settings.seed < 0)
	{
		settings.seed = (long)time(0);
	}
	std::cout << "The seed is " << settings.seed << std::endl;
	suni(settings.seed);


	//读取object参数，建立object参数数据库
	Database *ptr_vac_database = new Database("database_vac.txt");
	Database *ptr_sia_database = new Database("database_sia.txt");
	Database *ptr_fia_database = new Database("database_fia.txt");
	Database *ptr_vac_fia_database = new Database("database_vac_fia.txt");
	Database *ptr_sia_fia_database = new Database("database_sia_fia.txt");

	//读取object位置，用obj_ptr_list存储指向object指针
	//迭代器解引用之后得到指针，再次解引用才能得到Object对象
	std::vector<Object *> obj_ptr_list;	 	
	read_cascade("input_cascade.txt", obj_ptr_list, *ptr_vac_database, *ptr_sia_database, *ptr_fia_database,
		*ptr_vac_fia_database, *ptr_sia_fia_database, settings);


	//检查初始结合情况
	int pos1 = 0, pos2 = 0;
	while (check_dist_cross_bound(obj_ptr_list, pos1, pos2, settings))
	{
		carry_out_reaction(obj_ptr_list, pos1, pos2, *obj_ptr_list.at(pos1), *obj_ptr_list.at(pos2),
			*ptr_vac_database, *ptr_sia_database, *ptr_fia_database, *ptr_vac_fia_database, *ptr_sia_fia_database, settings);
		pos2 = 0;
	}


	//开始主循环之前输出一次
	output_dump(obj_ptr_list, settings, true);
	output_txt(obj_ptr_list, settings, true);
	output_cascade(obj_ptr_list, settings);


	//主循环
	std::cout << "Start main loop!\n" << std::endl;

	bool flag_refresh_event_list = true;
	double rate_sum;
	std::vector<Event> event_list;
	while (settings.check_end())
	{
		//如果obj_ptr_list发生变化，重新构建事件列表，并且更新rate_sum
		if (flag_refresh_event_list) {
			rate_sum = build_event_list(obj_ptr_list, event_list,settings);
			flag_refresh_event_list = false;
		}

		//如果事件列表不为空，rate_sum不为0，则抽取事件、执行事件、累加时间和步数
		if (rate_sum > 0.0)
		{
			//随机抽取事件
			//抽中的event在event_list中的位置
			int event_pos = find_event(event_list, rate_sum);

			//如果抽到的是cascade入射事件
			if (event_list.at(event_pos).which_event == -1)
			{
				//插入cascade
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
			//否则如果抽到的是object的事件
			else
			{
				//执行事件
				//可能的事件：迁移、旋转、发射、cascade入射
				//如果抽中迁移或者发射，则检查是否撞上别的obj
				int obj_pos = event_list.at(event_pos).pos;	//该obj在obj_ptr_list中的位置
				int which_event = event_list.at(event_pos).which_event;
				carry_out_event(which_event, obj_pos, obj_ptr_list, flag_refresh_event_list,
					*ptr_vac_database, *ptr_sia_database, *ptr_fia_database, *ptr_vac_fia_database, *ptr_sia_fia_database, settings);

			}
			settings.increase_time_step(rate_sum);
		}
		//如果事件列表为空，则程序空跑，每步时间累加default_dt；步数累加1
		else {
			settings.time += settings.default_dt;
			settings.step += 1;
		}

		//output_dump(obj_ptr_list, settings, false);

		//检查是否需要输出
		if (settings.check_output())
		{
			output_dump(obj_ptr_list, settings, false);
			output_txt(obj_ptr_list, settings, false);
			output_cascade(obj_ptr_list, settings);
		}
	}

	//如果是时间为截至标准，且最后一步执行完之后，时间超过max_time，则把时间拉回max_time
	if (settings.stop_criteria == 1 && settings.time > settings.max_time) {
		settings.time = settings.max_time;
	}

	//终止主循环之后再输出一次
	output_dump(obj_ptr_list, settings, false);
	output_txt(obj_ptr_list, settings, false);
	output_cascade(obj_ptr_list, settings);


	//测试代码

	//释放所有obj的内存
	for (auto iter : obj_ptr_list) {
		delete iter;
	}
	//释放所有database的内存
	delete ptr_vac_database;
	delete ptr_sia_database;
	delete ptr_fia_database;
	delete ptr_vac_fia_database;
	delete ptr_sia_fia_database;

	return 0;
}