#include"event.h"



//输入obj_ptr_list和空的event_list，输出装满所有可能发生的事件的event list，并且返回总rate
//event_list中的每个元素包含三个信息：
//1. 截止到该event的rate总和
//2. 执行该event的object在obj_ptr_list的位置
//3. 该object应该执行哪个事件（插入cascade~-1, mig~0, emit1~1, emit2~2 and rotate~3）
double build_event_list(std::vector<Object *> &_obj_ptr_list, std::vector<Event> &_event_list, Setting _set)
{
	//先清空之前的event_list，再重新往里面添加event
	_event_list.clear();
	double rate_sum = 0.0;
	for (auto iter = _obj_ptr_list.cbegin(); iter != _obj_ptr_list.cend(); iter++)
	{
		Event temp_event;
		//对整个frequency map循环
		for (auto &map_event_freq : (*iter)->frequency)
		{
			rate_sum += map_event_freq.second;
			temp_event.rate = rate_sum;
			temp_event.pos = iter - _obj_ptr_list.cbegin();
			temp_event.which_event = map_event_freq.first;
			_event_list.push_back(temp_event);
		}
	}

	if (_set.rate_cascade > 0)
	{
		Event temp_event;
		rate_sum += _set.rate_cascade;
		temp_event.rate = rate_sum;
		//插入cascade这个事件，Pos = -1, which_event = -1
		temp_event.pos = -1;
		temp_event.which_event = -1;
		_event_list.push_back(temp_event);
	}

	//在事件列表最后加上He入射的事件
	if (_set.He_flux > 0)
	{
		Event He_Insert;
		rate_sum += _set.He_flux * _set.box_length.at(0) * _set.box_length.at(1) * _set.a0 * _set.a0 * 1E-20;
		He_Insert.rate = rate_sum;
		He_Insert.pos = _event_list.size();
		He_Insert.which_event = 5;
		_event_list.push_back(He_Insert);
	}

	if (_event_list.size() == 0) {
		std::cout << "The event list is empty!\n";
	}
	return rate_sum;
}


//输入event list，产生一个随机数rand_num，返回一个event，它在event list中的位置为pos_i，该event满足如下条件
//event_list.at(pos_i-1).rate < rand_num < event_list.at(pos_i).rate
int find_event(std::vector<Event> &_event_list, double _rate_sum)
{
	double rand_num = uni();
	double rand_rate = rand_num * _rate_sum;
	//binary search 二分查找
	int start_pos = 0, end_pos = _event_list.size() - 1;
	while (start_pos < end_pos)
	{
		int mid_pos = (start_pos + end_pos) / 2;
		double product = (rand_rate - _event_list.at(mid_pos).rate) * (rand_rate - _event_list.at(mid_pos+1).rate);
		if (product < 0) {
			return mid_pos+1;
		}
		else if (rand_rate > _event_list.at(mid_pos+1).rate)
		{
			start_pos = mid_pos + 1;
		}
		else if (rand_rate <= _event_list.at(mid_pos).rate)
		{
			end_pos = mid_pos;
		}
	}
	return start_pos;
}


//输入object再obj_ptr_list中的位置，以及它要执行的反应，执行该反应
//如果是迁移事件，迁移之后，应检查碰撞，如果碰撞，_flag_refresh_event_list重置为true，准备重建event_list
//两种发射事件，均需要重建event_list
void carry_out_event(int _which_event, int _obj_pos, std::vector<Object *> &_obj_ptr_list, bool &_flag_refresh_event_list,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set)
{
	switch (_which_event)
	{
	//迁移
	case 0:
	{
		//执行迁移操作
		carry_out_mig(_obj_pos, _obj_ptr_list, _flag_refresh_event_list,
			_vac_database, _sia_database, _fia_database, _vac_fia_database, _sia_fia_database, _set);
		break;
	}
	case 1:		
		carry_out_emit1(_obj_ptr_list, *_obj_ptr_list.at(_obj_pos), _obj_pos,
			_vac_database, _sia_database, _fia_database, _vac_fia_database, _sia_fia_database, _set);
		_flag_refresh_event_list = true;
		break;
	case 2:
		carry_out_emit2(_obj_ptr_list, *_obj_ptr_list.at(_obj_pos), _obj_pos,
			_vac_database, _sia_database, _fia_database, _vac_fia_database, _sia_fia_database, _set);
		_flag_refresh_event_list = true;
		break;
	case 3:
		_obj_ptr_list.at(_obj_pos)->rotate_sia();
		break;
	case 4:
		carry_out_TM(_obj_ptr_list, *_obj_ptr_list.at(_obj_pos), _obj_pos,
			_vac_database, _sia_database, _fia_database, _vac_fia_database, _sia_fia_database, _set);
		_flag_refresh_event_list = true;
		//std::cout << "tm happens" << std::endl;
		break;
	case 5:
		carry_out_He_insert(_obj_ptr_list,
			_vac_database, _sia_database, _fia_database, _vac_fia_database, _sia_fia_database, _set);
		_flag_refresh_event_list = true;
		//std::cout << "He_insert happens" << std::endl;
		//		_flag_trap_mutation = true;		

		break;
	default:
		break;
	}
}


//执行迁移事件，迁移完之后
//1.检查边界条件，pbc下越界，则移回盒子内；abc越界，则kill
//2.如果是pbc，迁移完之后，根据sink strength检查是否撞到晶界
void carry_out_mig(int _obj_pos, std::vector<Object *> &_obj_ptr_list, bool &_flag_refresh_event_list,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set)
{

	_obj_ptr_list.at(_obj_pos)->mig(_set);

	//执行完迁移事件后，是否需要check_reaction
	bool flag_check_reaction = true;
	//如果是周期边界条件
	if (_set.boundary_condition == 2) {
		_obj_ptr_list.at(_obj_pos)->check_pbc(_set);
		_obj_ptr_list.at(_obj_pos)->refresh_grid_pos(_set);
		//设置一个开关，如果它被任一sink吸收之后，将其关闭
		bool check_sink_strength=true;
		//通过检查位移，判断Obj是否穿过晶界
		if (_set.grain_radius > 0)
		{
			//如果位移超过晶粒半径
			if (_obj_ptr_list.at(_obj_pos)->check_grain_radius(_set)) {
				if (_set.output_GB_absorption > 0) {
					output_react_1_0_GB(*(_obj_ptr_list.at(_obj_pos)), _set);
				}
				delete _obj_ptr_list.at(_obj_pos);
				_obj_ptr_list.erase(std::begin(_obj_ptr_list) + _obj_pos);
				//迁移的obj被杀掉了，所以要更新event list
				_flag_refresh_event_list = true;
				//迁移的那个Obj被杀掉了，所以不用check_reaction
				flag_check_reaction = false;
			}
		}
		//通过sink strength，判断Obj是否穿过晶界
		if (_set.ss_grain_radius > 0 && check_sink_strength)
		{
			int type = _obj_ptr_list.at(_obj_pos)->type;
			//如果Obj是空位，那么抽一个随机数检查是否小于p_vac_meet_GB，如果小于，kill该obj
			if (type == 1)
			{
				if (uni() < _set.p_vac_meet_GB)
				{
					if (_set.output_GB_absorption > 0) {
						output_react_1_0_GB(*(_obj_ptr_list.at(_obj_pos)), _set);
					}
					delete _obj_ptr_list.at(_obj_pos);
					_obj_ptr_list.erase(std::begin(_obj_ptr_list) + _obj_pos);
					//迁移的obj被杀掉了，所以要更新event list
					_flag_refresh_event_list = true;
					//迁移的那个Obj被杀掉了，所以不用check_reaction
					flag_check_reaction = false;
					//迁移的那个Obj被杀掉了，所以不用再判断是否被其他sink吸收
					check_sink_strength = false;
				}
			}
			//如果Obj是sia，那么抽一个随机数检查是否小于p_sia_meet_GB，如果小于，kill该obj
			else if (type == 2)
			{
				if (uni() < _set.p_sia_meet_GB)
				{
					if (_set.output_GB_absorption > 0) {
						output_react_1_0_GB(*(_obj_ptr_list.at(_obj_pos)), _set);
					}
					delete _obj_ptr_list.at(_obj_pos);
					_obj_ptr_list.erase(std::begin(_obj_ptr_list) + _obj_pos);
					//迁移的obj被杀掉了，所以要更新event list
					_flag_refresh_event_list = true;
					//迁移的那个Obj被杀掉了，所以不用check_reaction
					flag_check_reaction = false;
					//迁移的那个Obj被杀掉了，所以不用再判断是否被其他sink吸收
					check_sink_strength = false;
				}
			}
		}
		if (_set.dislocation_density > 0 && check_sink_strength)
		{
			if (uni() < _obj_ptr_list.at(_obj_pos)->p_meet_dislocation)
			{
				if (_set.output_disloaction_absorption > 0) {
					output_react_1_0_dislocation(*(_obj_ptr_list.at(_obj_pos)), _set);
				}
				delete _obj_ptr_list.at(_obj_pos);
				_obj_ptr_list.erase(std::begin(_obj_ptr_list) + _obj_pos);
				//迁移的obj被杀掉了，所以要更新event list
				_flag_refresh_event_list = true;
				//迁移的那个Obj被杀掉了，所以不用check_reaction
				flag_check_reaction = false;
				//迁移的那个Obj被杀掉了，所以不用再判断是否被其他sink吸收
				check_sink_strength = false;

			}

		}
		//迁移完之后，看看和其他Obj是否撞上
		if (flag_check_reaction)
		{
			int pos2 = 0;
			if (check_dist_cross_bound_inner_loop(_obj_ptr_list, _obj_pos, pos2, _set))
			{
				carry_out_reaction(_obj_ptr_list, _obj_pos, pos2, *_obj_ptr_list.at(_obj_pos), *_obj_ptr_list.at(pos2),
					_vac_database, _sia_database, _fia_database, _vac_fia_database, _sia_fia_database, _set);
				//如果执行了反应，代表有obj发生变化，因此需要更新event_list
				_flag_refresh_event_list = true;
			}
		}
	}
	//如果是吸收边界条件
	else if (_set.boundary_condition == 1) {
		//如果撞到边界了，就把obj去掉，并且不再check_reaction
		if (_obj_ptr_list.at(_obj_pos)->check_cross_bound_xyz(_set)>0) {
			if (_set.output_out_bound > 0) {
				output_react_1_0(*(_obj_ptr_list.at(_obj_pos)), _set);
			}
			delete _obj_ptr_list.at(_obj_pos);
			_obj_ptr_list.erase(std::begin(_obj_ptr_list) + _obj_pos);
			//迁移的obj被杀掉了，所以要更新event list
			_flag_refresh_event_list = true;
			//迁移的那个Obj被杀掉了，所以不用check_reaction
			flag_check_reaction = false;
		}
		//如果没有撞到边界，就应该更新以下Object所在的网格
		else {
			_obj_ptr_list.at(_obj_pos)->refresh_grid_pos(_set);
		}
		//迁移完之后，看看和其他Obj是否撞上
		if (flag_check_reaction)
		{
			int pos2 = 0;
			if (check_dist_inner_loop_abc(_obj_ptr_list, _obj_pos, pos2, _set))
			{
				carry_out_reaction(_obj_ptr_list, _obj_pos, pos2, *_obj_ptr_list.at(_obj_pos), *_obj_ptr_list.at(pos2),
					_vac_database, _sia_database, _fia_database, _vac_fia_database, _sia_fia_database, _set);
				//如果执行了反应，代表有obj发生变化，因此需要更新event_list
				_flag_refresh_event_list = true;
			}
		}
	}
	//如果是XY周期，Z开放
	else if (_set.boundary_condition == 3)
	{
		if (_obj_ptr_list.at(_obj_pos)->check_cross_bound_z(_set))
		{
			if (_set.output_out_bound > 0) {
				output_react_1_0(*(_obj_ptr_list.at(_obj_pos)), _set);
			}
			delete _obj_ptr_list.at(_obj_pos);
			_obj_ptr_list.erase(std::begin(_obj_ptr_list) + _obj_pos);
			//迁移的obj被杀掉了，所以要更新event list
			_flag_refresh_event_list = true;
			//迁移的那个Obj被杀掉了，所以不用check_reaction
			flag_check_reaction = false;
		}
		else
		{
			_obj_ptr_list.at(_obj_pos)->check_pbc(_set);
			_obj_ptr_list.at(_obj_pos)->refresh_grid_pos(_set);
		}
		if (flag_check_reaction)
		{
			int pos2 = 0;
			if (check_dist_cross_bound_inner_loop(_obj_ptr_list, _obj_pos, pos2, _set))
			{
				carry_out_reaction(_obj_ptr_list, _obj_pos, pos2, *_obj_ptr_list.at(_obj_pos), *_obj_ptr_list.at(pos2),
					_vac_database, _sia_database, _fia_database, _vac_fia_database, _sia_fia_database, _set);
				//如果执行了反应，代表有obj发生变化，因此需要更新event_list
				_flag_refresh_event_list = true;
			}
		}
	}

}


//发射主元素（空位或者自间隙）
//obj1是被发射的元素，它的size为1
//大元素Obj2位置不变，obj1随机发射到球面上（球半径为母obj的半径加3）
void carry_out_emit1(std::vector<Object *> &_obj_ptr_list, Object &_obj, int _obj_pos,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set)
{
	//创建两个新的obj
	//确定两个待创建的obj的构造参数
	int o1_type, o1_size1, o1_size2, o1_dir;
	double x1, y1, z1;
	int o2_type, o2_size1, o2_size2, o2_dir;
	double x2 = _obj.pos.at(0), y2 = _obj.pos.at(1), z2 = _obj.pos.at(2);
	//产生一个0-2pi之间的随机弧度值，确定被发射元素的坐标
	double angle1 = uni() * 2 * 3.1415926;
	double angle2 = uni() * 2 * 3.1415926;
	//球半径为母obj的半径加3
	double r = _obj.radius + 3;
	x1 = x2 + r * sin(angle1) * cos(angle2);
	y1 = y2 + r * sin(angle1) * sin(angle2);
	z1 = z2 + r * cos(angle1);

	if (_obj.type == 1) {
		o1_type = 1; o1_size1 = 1;				o1_size2 = 0; o1_dir = 0;
		o2_type = 1; o2_size1 = _obj.size1 - 1; o2_size2 = 0; o2_dir = 0;
	}
	else if (_obj.type == 2) {
		o1_type = 2; o1_size1 = 1;				o1_size2 = 0; o1_dir = (int)(uni() * 4) + 1;
		o2_type = 2; o2_size1 = _obj.size1 - 1; o2_size2 = 0; o2_dir = _obj.dir;
	}
	else if (_obj.type == 4) {
		//vac + fia团簇
		if (_obj.size1 == 1) {
			o1_type = 1; o1_size1 = 1; o1_size2 = 0;          o1_dir = 0;
			o2_type = 3; o2_size1 = 0; o2_size2 = _obj.size2; o2_dir = 0;
		}
		// vac + vac_fia
		else {
			o1_type = 1; o1_size1 = 1;             o1_size2 = 0;          o1_dir = 0;
			o2_type = 4; o2_size1 = _obj.size1 - 1; o2_size2 = _obj.size2; o2_dir = 0;
		}
	}
	else if (_obj.type == 5) {
		//sia + fia团簇
		if (_obj.size1 == 1) {
			o1_type = 2; o1_size1 = 1; o1_size2 = 0;          o1_dir = (int)(uni() * 4) + 1;
			o2_type = 3; o2_size1 = 0; o2_size2 = _obj.size2; o2_dir = _obj.dir;
		}
		// sia + sia_fia
		else {
			o1_type = 2; o1_size1 = 1;              o1_size2 = 0;          o1_dir = (int)(uni() * 4) + 1;
			o2_type = 5; o2_size1 = _obj.size1 - 1; o2_size2 = _obj.size2; o2_dir = _obj.dir;
		}

	}

	//创建两个新的Obj，obj1是被吐出来的
	Object *ptr_obj1 = new Object(o1_type, x1, y1, z1, o1_size1, o1_size2, o1_dir, _vac_database, _sia_database,
		_fia_database, _vac_fia_database, _sia_fia_database, _set);
	Object *ptr_obj2 = new Object(o2_type, x2, y2, z2, o2_size1, o2_size2, o2_dir, _vac_database, _sia_database,
		_fia_database, _vac_fia_database, _sia_fia_database, _set);
	_obj_ptr_list.push_back(ptr_obj1);
	_obj_ptr_list.push_back(ptr_obj2);

	//输出
	if (_set.output_emit_reaction > 0) {
		output_react_1_2(*_obj_ptr_list.at(_obj_pos), *ptr_obj1, *ptr_obj2, _set);
	}

	//删除老的Obj
	delete _obj_ptr_list.at(_obj_pos);
	_obj_ptr_list.erase(std::begin(_obj_ptr_list) + _obj_pos);

	//check_reaction，对吐出来的obj循环一轮即可，pos1是被吐出来的obj的索引
	int pos1 = _obj_ptr_list.size() - 2, pos2 = 0;
	if (check_dist_cross_bound_inner_loop(_obj_ptr_list, pos1, pos2,_set))
	{
		carry_out_reaction(_obj_ptr_list, pos1, pos2, *_obj_ptr_list.at(pos1), *_obj_ptr_list.at(pos2),
			_vac_database, _sia_database, _fia_database, _vac_fia_database, _sia_fia_database, _set);
	}
}


//发射副元素（fia）
//obj1是被发射的元素，它的size为1
//大元素Obj2位置不变，obj1随机发射到球面上（球半径为母obj的半径加2）
void carry_out_emit2(std::vector<Object *> &_obj_ptr_list, Object &_obj, int _obj_pos,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set)
{
	//创建两个新的obj
	//确定两个待创建的obj的构造参数
	int o1_type, o1_size1, o1_size2, o1_dir;
	double x1, y1, z1;
	int o2_type, o2_size1, o2_size2, o2_dir;
	double x2 = _obj.pos.at(0), y2 = _obj.pos.at(1), z2 = _obj.pos.at(2);
	//产生一个0-2pi之间的随机弧度值，确定被发射元素的坐标
	double angle1 = uni() * 2 * 3.1415926;
	double angle2 = uni() * 2 * 3.1415926;
	//球半径为母obj的半径加3
	double r = _obj.radius + 3;
	x1 = x2 + r * sin(angle1) * cos(angle2);
	y1 = y2 + r * sin(angle1) * sin(angle2);
	z1 = z2 + r * cos(angle1);

	// fia + fia
	if (_obj.type == 3) {
		o1_type = 3; o1_size1 = 0; o1_size2 = 1;              o1_dir = 0;
		o2_type = 3; o2_size1 = 0; o2_size2 = _obj.size2 - 1; o2_dir = 0;
	}
	else if (_obj.type == 4) {
		//fia + vac团簇
		if (_obj.size2 == 1) {
			o1_type = 3; o1_size1 = 0;          o1_size2 = 1; o1_dir = 0;
			o2_type = 1; o2_size1 = _obj.size1; o2_size2 = 0; o2_dir = 0;
		}
		// fia + vac_fia
		else {
			o1_type = 3; o1_size1 = 0;          o1_size2 = 1;            o1_dir = 0;
			o2_type = 4; o2_size1 = _obj.size1; o2_size2 = _obj.size2-1; o2_dir = 0;
		}
	}
	else if (_obj.type == 5) {
		//fia + sia团簇
		if (_obj.size2 == 1) {
			o1_type = 3; o1_size1 = 0;          o1_size2 = 1; o1_dir = 0;
			o2_type = 2; o2_size1 = _obj.size1; o2_size2 = 0; o2_dir = _obj.dir;
		}
		// fia + sia_fia
		else {
			o1_type = 3; o1_size1 = 0;              o1_size2 = 1;        o1_dir = 0;
			o2_type = 5; o2_size1 = _obj.size1; o2_size2 = _obj.size2-1; o2_dir = _obj.dir;
		}
	}

	//创建两个新的Obj，obj1是被吐出来的元素
	Object *ptr_obj1 = new Object(o1_type, x1, y1, z1, o1_size1, o1_size2, o1_dir, _vac_database, _sia_database,
		_fia_database, _vac_fia_database, _sia_fia_database, _set);
	Object *ptr_obj2 = new Object(o2_type, x2, y2, z2, o2_size1, o2_size2, o2_dir, _vac_database, _sia_database,
		_fia_database, _vac_fia_database, _sia_fia_database, _set);
	_obj_ptr_list.push_back(ptr_obj1);
	_obj_ptr_list.push_back(ptr_obj2);

	//输出
	if (_set.output_emit_reaction > 0) {
		output_react_1_2(*_obj_ptr_list.at(_obj_pos), *ptr_obj1, *ptr_obj2, _set);
	}

	//删除老的Obj
	delete _obj_ptr_list.at(_obj_pos);
	_obj_ptr_list.erase(std::begin(_obj_ptr_list) + _obj_pos);

	//check_reaction，对吐出来的obj循环一轮即可
	int pos1 = _obj_ptr_list.size() - 2, pos2 = 0;
	if (check_dist_cross_bound_inner_loop(_obj_ptr_list, pos1, pos2,_set))
	{
		carry_out_reaction(_obj_ptr_list, pos1, pos2, *_obj_ptr_list.at(pos1), *_obj_ptr_list.at(pos2),
			_vac_database, _sia_database, _fia_database, _vac_fia_database, _sia_fia_database, _set);
	}
}


//随机选择一个cascade
std::string choose_cascade_to_insert(Setting &_set)
{
	double rand_num = uni();
	for (auto iter = _set.cascade_possibility.cbegin(); iter != _set.cascade_possibility.cend(); iter++)
	{
		if (*iter > rand_num)
		{
			int num = iter - _set.cascade_possibility.cbegin();
			std::ostringstream os;
			os << "cascade." << _set.cascade_energy.at(num);
			double x = uni();
			int n = (int)(x * _set.cascade_number.at(num)) + 1;
			os << '.' << n << ".txt";
			return os.str();
		}
	}
	return "fail to choose cascade_to_insert\n";
}


//trap mutation事件
//He-V就是当可以发生trap mutation事件的时候，来一个一个发生SIA原子。
//He-SIA就是当可以发生trap mutation事件的时候，将整个SIA团簇一下子都发射出去。

void carry_out_TM(std::vector<Object*>& _obj_ptr_list, Object& _obj, int _obj_pos,
	const Database& _vac_database, const Database& _sia_database, const Database& _fia_database,
	const Database& _vac_fia_database, const Database& _sia_fia_database, const Setting& _set)
{
	//创建两个新的obj
	//确定两个待创建的obj的构造参数
	//SIA
	int ol_size1(_obj.size1), ol_size2(_obj.size2), tm_sia(0);
	int o1_type, o1_size1, o1_size2, o1_dir;
	o1_type = 2;
	//o1_size1 = 1;//o1_size1 = _obj.size2/4 - _obj.size1;
	if (_obj.type == 5)
	{
		o1_size1 = _obj.size1 + 1;
	}
	else if (_obj.type == 4)
	{
		while (ol_size2 > 5 * pow(ol_size1, 0.86))
		{
			ol_size1++;
			tm_sia++;
		}
		o1_size1 = tm_sia;
	}
	else
		o1_size1 = 1;
	o1_size2 = 0; o1_dir = (int)floor(uni() * 4) + 1;
	//He-V
	double x1, y1, z1;
	int o2_type, o2_size1, o2_size2, o2_dir;
	o2_type = 4;
	if (_obj.type == 5)
	{
		o2_size1 = 1;
	}
	else if (_obj.type == 4)
	{
		while (ol_size2 > 5 * pow(ol_size1, 0.86))
		{
			ol_size1++;
			tm_sia++;
		}
		o2_size1 = tm_sia + _obj.size1;
	}
	else
		o2_size1 = _obj.size1 + 1;
	o2_size2 = _obj.size2, o2_dir = 0;
	double x2 = _obj.pos.at(0), y2 = _obj.pos.at(1), z2 = _obj.pos.at(2);
	//产生一个0-2pi之间的随机弧度值，确定被发射元素的坐标
	double angle1 = uni() * 2 * 3.1415926;
	double angle2 = uni() * 2 * 3.1415926;
	//球半径为母obj的半径加3
	double o1_r = pow(o1_size1 / pow(3, 0.5) / PI, 0.5) + 0.4813; //SIA 半径
	double o2_r = pow((3 * o2_size1 / 8.0 / PI), 1 / 3.0) + 0.34; //V-He 半径
	double r = o1_r + o2_r + 3;//double r = _obj.radius + 3;
	x1 = x2 + r * sin(angle1) * cos(angle2);
	y1 = y2 + r * sin(angle1) * sin(angle2);
	z1 = z2 + r * cos(angle1);


	//创建两个新的Obj，obj1是被吐出来的
	Object* ptr_obj1 = new Object(o1_type, x1, y1, z1, o1_size1, o1_size2, o1_dir, _vac_database, _sia_database,
		_fia_database, _vac_fia_database, _sia_fia_database, _set);
	Object* ptr_obj2 = new Object(o2_type, x2, y2, z2, o2_size1, o2_size2, o2_dir, _vac_database, _sia_database,
		_fia_database, _vac_fia_database, _sia_fia_database, _set);

	_obj_ptr_list.push_back(ptr_obj1);
	_obj_ptr_list.push_back(ptr_obj2);

	//输出
	
	if (_set.output_trap_mutation > 0) {
		output_react_TM(*_obj_ptr_list.at(_obj_pos), *ptr_obj1, *ptr_obj2, _set);
	}
	

	//删除老的Obj
	delete _obj_ptr_list.at(_obj_pos);
	_obj_ptr_list.erase(std::begin(_obj_ptr_list) + _obj_pos);

	//check_reaction，对新生成的obj都需要循环一轮
	int pos1 = _obj_ptr_list.size() - 2, pos2 = 0;
	if (check_dist_cross_bound_inner_loop(_obj_ptr_list, pos1, pos2, _set))
	{
		carry_out_reaction(_obj_ptr_list, pos1, pos2, *_obj_ptr_list.at(pos1), *_obj_ptr_list.at(pos2),
			_vac_database, _sia_database, _fia_database, _vac_fia_database, _sia_fia_database, _set);
	}

	pos1 = _obj_ptr_list.size() - 1; pos2 = 0;
	if (check_dist_cross_bound_inner_loop(_obj_ptr_list, pos1, pos2, _set))
	{
		carry_out_reaction(_obj_ptr_list, pos1, pos2, *_obj_ptr_list.at(pos1), *_obj_ptr_list.at(pos2),
			_vac_database, _sia_database, _fia_database, _vac_fia_database, _sia_fia_database, _set);
	}
}


//He离子注入
void carry_out_He_insert(std::vector<Object*>& _obj_ptr_list,
	const Database& _vac_database, const Database& _sia_database, const Database& _fia_database,
	const Database& _vac_fia_database, const Database& _sia_fia_database, const Setting& _set)
{
	int o_type, o_size1, o_size2, o_dir;
	double x, y, z;

	//产生一个He
	o_type = 3; o_size1 = 0; o_size2 = 1, o_dir = 0;
	//位置 x,y 方向均匀分布  z方向满足一个分布函数
	x = uni() * _set.box_length.at(0);
	y = uni() * _set.box_length.at(1);
	//z = (0.16/3.16 * 10) - (1.65/3.16 * 10) * log(1.0-uni());
	// z = (0.16 / _set.a0 * 10) - (1.65 / _set.a0 * 10) * log(1.0 - uni());
	z = uni() * _set.box_length.at(2);

	Object* ptr_obj = new Object(o_type, x, y, z, o_size1, o_size2, o_dir, _vac_database, _sia_database,
		_fia_database, _vac_fia_database, _sia_fia_database, _set);
	_obj_ptr_list.push_back(ptr_obj);

	////输出He植入事件
	//if (_set.output_c_a_e_react > 0) {
	//	output_react_He_Insert(*ptr_obj, _set);
	//}

	//check_reaction，对新生成的He执行一轮即可
	int pos1 = _obj_ptr_list.size() - 1, pos2 = 0;
	if (check_dist_cross_bound_inner_loop(_obj_ptr_list, pos1, pos2, _set))
	{
		carry_out_reaction(_obj_ptr_list, pos1, pos2, *_obj_ptr_list.at(pos1), *_obj_ptr_list.at(pos2),
			_vac_database, _sia_database, _fia_database, _vac_fia_database, _sia_fia_database, _set);
	}
}
