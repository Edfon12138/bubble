#include"output.h"


//初始化output files
void initialize_output_file(std::string _filename)
{
	std::ostringstream oss;
	std::ofstream fout;
	fout.open(_filename);
	fout << "step\ttime\tid\ttype\tx\ty\tz\tsize1\tsize2\tdir\n";
	fout.close();
}


//初始化"output_react.txt"
void initialize_react_file()
{
	std::ostringstream oss;
	std::ofstream fout;
	fout.open("output_react.txt");
	fout << "#step\ttime\tcpu_time\t" << "id\ttype\tsize1\tsize2\t" << "id\ttype\tsize1\tsize2\t" << "react\t"
		<< "id\ttype\tsize1\tsize2\t" << "id\ttype\tsize1\tsize2\n";
	fout.close();
}


//两个Object完全湮灭，产生0个新object
//输入两个object，将反应信息写入"react.txt"
void output_react_2_0(const Object &_obj1, const Object &_obj2, const Setting &_set)
{
	double duration = cal_duration();
	std::ostringstream oss;
	oss << _set.step << '\t' << _set.time << '\t' << duration << '\t'
		<< _obj1.id << '\t' << _obj1.type << '\t' << _obj1.size1 << '\t' << _obj1.size2 << '\t'
		<< _obj2.id << '\t' << _obj2.type << '\t' << _obj2.size1 << '\t' << _obj2.size2 << '\t'
		<< "anneal\t-\t-\t-\t-\t-\t-\t-\t-\n";
	write_content("output_react.txt", oss, false);
}

// 输出合并反应：obj1 + obj2 = obj3
void output_react_2_1(const Object &_obj1, const Object &_obj2, const Object &_obj3, const Setting &_set)
{
	double duration = cal_duration();
	std::ostringstream oss;
	oss << _set.step << '\t'<<  _set.time << '\t' << duration << '\t'
		<< _obj1.id << '\t' << _obj1.type << '\t' << _obj1.size1 << '\t' << _obj1.size2 << '\t'
		<< _obj2.id << '\t' << _obj2.type << '\t' << _obj2.size1 << '\t' << _obj2.size2 << '\t'
		<< "react\t"
		<< _obj3.id << '\t' << _obj3.type << '\t' << _obj3.size1 << '\t' << _obj3.size2 << '\t'
		<< "-\t-\t-\t-\n";
	write_content("output_react.txt", oss, false);
}

//输出emit反应
void output_react_1_2(const Object &_obj1, const Object &_obj2, const Object &_obj3, const Setting &_set)
{
	double duration = cal_duration();
	std::ostringstream oss;
	oss << _set.step << '\t' << _set.time << '\t' << duration << '\t'
		<< _obj1.id << '\t' << _obj1.type << '\t' << _obj1.size1 << '\t' << _obj1.size2 << '\t'
		<< "-\t-\t-\t-\temit1\t"
		<< _obj2.id << '\t' << _obj2.type << '\t' << _obj2.size1 << '\t' << _obj2.size2 << '\t'
		<< _obj3.id << '\t' << _obj3.type << '\t' << _obj3.size1 << '\t' << _obj3.size2 << '\n';
	write_content("output_react.txt", oss, false);
}

//obj reach bound in absorb boundary condition
void output_react_1_0(const Object &_obj1, const Setting &_set)
{
	double duration = cal_duration();
	std::ostringstream oss;
	oss << _set.step << '\t' << _set.time << '\t' << duration << '\t'
		<< _obj1.id << '\t' << _obj1.type << '\t' << _obj1.size1 << '\t' << _obj1.size2 << '\t'
		<< "-\t-\t-\t-\t"
		<< "out_bound\t"
		<< "-\t-\t-\t-\t-\t-\t-\t-\n";
	write_content("output_react.txt", oss, false);
}

//obj被晶界吸收
void output_react_1_0_GB(const Object &_obj1, const Setting &_set)
{
	double duration = cal_duration();
	std::ostringstream oss;
	oss << _set.step << '\t' << _set.time << '\t' << duration << '\t'
		<< _obj1.id << '\t' << _obj1.type << '\t' << _obj1.size1 << '\t' << _obj1.size2 << '\t'
		<< "-\t-\t-\t-\t"
		<< "GB\t"
		<< "-\t-\t-\t-\t-\t-\t-\t-\n";
	write_content("output_react.txt", oss, false);
}

//obj被位错吸收
void output_react_1_0_dislocation(const Object &_obj1, const Setting &_set)
{
	double duration = cal_duration();
	std::ostringstream oss;
	oss << _set.step << '\t' << _set.time << '\t' << duration << '\t'
		<< _obj1.id << '\t' << _obj1.type << '\t' << _obj1.size1 << '\t' << _obj1.size2 << '\t'
		<< "-\t-\t-\t-\t"
		<< "dislocation\t"
		<< "-\t-\t-\t-\t-\t-\t-\t-\n";
	write_content("output_react.txt", oss, false);
}


//cascade insert
void output_react_0_0(const Setting &_set, std::string str)
{
	double duration = cal_duration();
	std::ostringstream oss;
	oss << _set.step << '\t' << _set.time << '\t' << duration << '\t'
		<< "-\t-\t-\t-\t-\t-\t-\t-\t"
		<< str 
		<< "\t-\t-\t-\t-\t-\t-\t-\t-\n";
	write_content("output_react.txt", oss, false);
}

// 在固定的时刻，将所有object的信息，以dump的形式输出
void output_dump(const std::vector<Object *> &_obj_ptr_list, const Setting &_set, bool _initialize = false)
{
	std::ostringstream oss;
	oss << "ITEM: TIMESTEP\n\t" << _set.step << "\nITEM: NUMBER OF ATOMS\n\t" << _obj_ptr_list.size() << "\nITEM: BOX BOUNDS\n\t"
		<< _set.box_min.at(0) << '\t' << _set.box_max.at(0) << "\n\t"
		<< _set.box_min.at(1) << '\t' << _set.box_max.at(1) << "\n\t"
		<< _set.box_min.at(2) << '\t' << _set.box_max.at(2) << '\n'
		<< "ITEM: ATOMS id type x y z size1 size2 radius time\n";
	for (auto iter = _obj_ptr_list.cbegin(); iter != _obj_ptr_list.cend(); iter++)
	{
		oss << '\t' << (*iter)->id << '\t' << (*iter)->type << '\t' <<
			(*iter)->pos.at(0) << '\t' << (*iter)->pos.at(1) << '\t' << (*iter)->pos.at(2) << '\t' <<
			(*iter)->size1 << '\t' << (*iter)->size2 << '\t' << (*iter)->radius << '\t' <<
			_set.time << "\n";
	}
	write_content("output_dump.dump", oss, _initialize);
}

//在固定时刻，输出每一个Object的尺寸信息，true代表初始化输出文件
void output_txt(const std::vector<Object *> &_obj_ptr_list, const Setting &_set, bool _initialize = false)
{
	std::ostringstream oss_vac, oss_sia, oss_fia, oss_vac_fia, oss_sia_fia;
	double duration = cal_duration();
	oss_vac << "step\t" << _set.step << "\ttime\t" << _set.time << "\ttemperature\t" << _set.temperature 
		<< "\tcpu_time\t" << duration << std::endl;
	oss_sia << "step\t" << _set.step << "\ttime\t" << _set.time << "\ttemperature\t" << _set.temperature 
		<< "\tcpu_time\t" << duration << std::endl;
	oss_fia << "step\t" << _set.step << "\ttime\t" << _set.time << "\ttemperature\t" << _set.temperature 
		<< "\tcpu_time\t" << duration << std::endl;
	oss_vac_fia << "step\t" << _set.step << "\ttime\t" << _set.time << "\ttemperature\t" << _set.temperature 
		<< "\tcpu_time\t" << duration << std::endl;
	oss_sia_fia << "step\t" << _set.step << "\ttime\t" << _set.time << "\ttemperature\t" << _set.temperature 
		<< "\tcpu_time\t" << duration << std::endl;
	for (auto iter : _obj_ptr_list)
	{
		switch (iter->type)
		{
		case 1:
			oss_vac << iter->size1 << std::endl;
			break;
		case 2:
			oss_sia << iter->size1 << std::endl;
			break;
		case 3:
			oss_fia << iter->size2 << std::endl;
			break;
		case 4:
			oss_vac_fia << iter->size1 << '\t' << iter->size2 << std::endl;
			break;
		case 5:
			oss_sia_fia << iter->size1 << '\t' << iter->size2 << std::endl;
			break;
		}
	}
	write_content("time_vac.txt", oss_vac, _initialize);
	write_content("time_sia.txt", oss_sia, _initialize);
	write_content("time_fia.txt", oss_fia, _initialize);
	write_content("time_vac_fia.txt", oss_vac_fia, _initialize);
	write_content("time_sia_fia.txt", oss_sia_fia, _initialize);
}

void output_cascade(const std::vector<Object *> &_obj_ptr_list, const Setting &_set)
{
	std::ostringstream oss;
	oss << "#steps\t" << _set.step << "\ttime\t" << _set.time << "\tseed\t" << _set.seed << '\n' <<
		"#type 1~vac, 2~sia, 3~fia, 4~vac_fia, 5~sia_fia\n"
		"#no direction~0;	111~1;	-111~2;	1 - 11~3;	11 - 1~4\n"<<
		"#id	type	x	y	z	size1	size2	dir\n";
	for (auto iter : _obj_ptr_list)
	{
		Object &obj = *iter;
		oss << obj.id << '\t' << obj.type << '\t' << 
			obj.pos.at(0) << '\t' << obj.pos.at(1) << '\t' << obj.pos.at(2) << '\t' << 
			obj.size1 << '\t' << obj.size2 << '\t' << obj.dir << std::endl;
	}
	write_content("output_cascade.txt", oss, true);
}

//输出trap_mutation反应
void output_react_TM(const Object& _obj1, const Object& _obj2, const Object& _obj3, const Setting& _set)
{

	double duration = cal_duration();
	std::ostringstream oss;
	oss << _set.step << '\t' << _set.time << '\t' << duration << '\t'
		<< _obj1.id << '\t' << _obj1.type << '\t' << _obj1.size1 << '\t' << _obj1.size2 << '\t'
		<< "-\t-\t-\t-\ttrap_mutation\t"
		<< _obj2.id << '\t' << _obj2.type << '\t' << _obj2.size1 << '\t' << _obj2.size2 << '\t'
		<< _obj3.id << '\t' << _obj3.type << '\t' << _obj3.size1 << '\t' << _obj3.size2 << '\n';
	write_content("output_react.txt", oss, false);
}


//输入文件名，字符串流，即可将内容添加到相应文件中
//第三个参数为true，对应'w'模式；第三个参数为false，对应'a'模式
void write_content(std::string _filename, std::ostringstream &_oss, bool _initialize)
{
	std::ofstream out_file;
	// 覆盖模式写入，第一次调用时使用
	if (_initialize) {
		out_file.open(_filename);
	}
	// 追加模式写入
	else {
		out_file.open(_filename, std::ios::out | std::ios::app);
	}
	out_file << _oss.str();
	out_file.close();
}

