#include"setting.h"


// 输入文件名，调用read_setting函数，生成setting_map，并用setting_map初始化得到一个setting对象
Setting::Setting(std::string filename)
{
	read_setting(filename);
	temperature = setting_map["temperature"];
	box_min.at(0) = int(setting_map["x_min"]);
	box_min.at(1) = int(setting_map["y_min"]);
	box_min.at(2) = int(setting_map["z_min"]);
	box_max.at(0) = int(setting_map["x_max"]);
	box_max.at(1) = int(setting_map["y_max"]);
	box_max.at(2) = int(setting_map["z_max"]);

	box_grid = int(setting_map["box_grid"]);
	for (int i = 0; i < 3; i++)
	{
		box_grid_num[i] = (int)ceil((box_max[i] - box_min[i]) / (box_grid)) - 1;
	}

	for (int i = 0; i < 3; i++) { box_length[i] = box_max.at(i) - box_min.at(i); }
	max_step = (unsigned long long)setting_map["max_step"];
	max_time = setting_map["max_time"];
	a0 = setting_map["a0"];
	kb = setting_map["kb"];
	hp = setting_map["hp"];
	stop_criteria = int(setting_map["stop_criteria"]);
	boundary_condition = int(setting_map["boundary_condition"]);
	output_cascade_injection = setting_map["output_cascade_injection"];
	output_combine_reaction = setting_map["output_combine_reaction"];
	output_annihilation_reaction = setting_map["output_annihilation_reaction"];
	output_emit_reaction = setting_map["output_emit_reaction"];
	output_out_bound = setting_map["output_out_bound"];
	output_GB_absorption = setting_map["output_GB_absorption"];
	output_disloaction_absorption = setting_map["output_disloaction_absorption"];
	output_trap_mutation = setting_map["output_trap_mutation"];

	total_output = int(setting_map["total_output"]);
	n_output = 1;
	fia_site = int(setting_map["fia_site"]);
	cascade_interval = setting_map["cascade_interval"];
	rand_translate_cascade = (int)setting_map["rand_translate_cascade"];
	rand_rotate_cascade = (int)setting_map["rand_rotate_cascade"];
	mig_degeneration = (int)setting_map["mig_degeneration"];
	max_energy = setting_map["max_energy"];
	grain_radius = setting_map["grain_radius"];
	grain_radius_square = pow(grain_radius, 2);
	default_dt = setting_map["default_dt"];
	ss_grain_radius = setting_map["ss_grain_radius"];
	dislocation_density = setting_map["dislocation_density"];
	dislocation_radius = setting_map["dislocation_radius"];
	seed = setting_map["seed"];
	cal_sink_strength();

	He_flux = setting_map["He_flux"];

	cal_cascade_insert_rate();


	//下一次插入cascade的时间
	cascade_insert_time = cascade_interval;
	//程序运行的百分比进度
	progress_percentage = 1;
	
	step = 0;
	time = 0.0;
	
	bcc_lattice_site = {
		0.0, 0.0, 0.0,
		0.5, 0.5, 0.5,
		1.0, 0.0, 0.0,
		0.0, 0.0, 1.0,
		0.0, 0.0, 1.0,
		1.0, 1.0, 0.0,
		1.0, 0.0, 1.0,
		0.0, 1.0, 1.0,
		1.0, 1.0, 1.0 };

	tis_site = {
		0.25, 0.50, 0.0,   0.25, 0.50, 1.0,
		0.50, 0.25, 0.0,   0.50, 0.25, 1.0,
		0.75, 0.50, 0.0,   0.75, 0.50, 1.0,
		0.50, 0.75, 0.0,   0.50, 0.75, 1.0,

		0.25, 0.0, 0.50,   0.25, 1.0, 0.50,
		0.50, 0.0, 0.25,   0.50, 1.0, 0.25,
		0.75, 0.0, 0.50,   0.75, 1.0, 0.50,
		0.50, 0.0, 0.75,   0.50, 1.0, 0.75,

		0.0, 0.50, 0.25,   1.0, 0.50, 0.25,
		0.0, 0.25, 0.50,   1.0, 0.25, 0.50,
		0.0, 0.50, 0.75,   1.0, 0.50, 0.75,
		0.0, 0.75, 0.50,   1.0, 0.75, 0.50, };


	read_cascade_possibility();
}

// 输入文件名，读取文件中所有非空，非#开头的行，以第一个关键字为键，第二个关键字为值，构造setting map
void Setting::read_setting(std::string _filename)
{
	std::fstream inFile(_filename);
	if (inFile.is_open())
	{
		std::string str;
		while (std::getline(inFile, str))
		{
			if (!str.empty() && str.at(0) != '#')	//skip the empty line and the line started with '#'
			{
				std::istringstream ss(str);
				std::string temp;
				ss >> temp;
				ss >> setting_map[temp];
			}
		}
		std::cout << _filename << " is succssefully read.\n";
		inFile.close();
	}
	else {
		std::cout << "Fail to open " << _filename << std::endl;
	}
}

//读取cascade_possibility文件，记录待插入cascade的信息
void Setting::read_cascade_possibility()
{
	std::string str;
	std::ifstream inFile("input_cascade_possibility.txt");
	if (inFile.is_open()) {
		while (std::getline(inFile, str))
		{
			if (str.empty() || str.at(0) == '#') { continue; }	//skip the empty line and the line started with '#'
			std::istringstream ss(str);
			double e, p;
			int n;
			ss >> e >> p >> n;
			cascade_energy.push_back(e);
			cascade_possibility.push_back(p);
			cascade_number.push_back(n);
		}
		std::cout << "input_cascade_possibility.txt is successfully read\n";
		inFile.close();
	}
	else {
		std::cout << "Fail to open input_cascade_possibility.txt\n";
	}
}

void Setting::cal_sink_strength()
{
	if (ss_grain_radius > 0)
	{
		p_sia_meet_GB = 5.625 * a0 * a0 / ss_grain_radius / ss_grain_radius;
		p_vac_meet_GB = 1.8 * a0 * a0 / ss_grain_radius / ss_grain_radius;
		std::cout << "P_sia_meet_GB is " << p_sia_meet_GB << std::endl;
		std::cout << "P_vac_meet_GB is " << p_vac_meet_GB << std::endl;
	}
}

void Setting::cal_cascade_insert_rate()
{
	if (cascade_interval > 0)
	{
		rate_cascade = 1.0 / cascade_interval;
	}
	else
		rate_cascade = -1.0;
	std::cout << "The cascade insert rate is: " << rate_cascade << " per second." << std::endl;
}