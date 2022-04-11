#include"event.h"



//����obj_ptr_list�Ϳյ�event_list�����װ�����п��ܷ������¼���event list�����ҷ�����rate
//event_list�е�ÿ��Ԫ�ذ���������Ϣ��
//1. ��ֹ����event��rate�ܺ�
//2. ִ�и�event��object��obj_ptr_list��λ��
//3. ��objectӦ��ִ���ĸ��¼�������cascade~-1, mig~0, emit1~1, emit2~2 and rotate~3��
double build_event_list(std::vector<Object *> &_obj_ptr_list, std::vector<Event> &_event_list, Setting _set)
{
	//�����֮ǰ��event_list�����������������event
	_event_list.clear();
	double rate_sum = 0.0;
	for (auto iter = _obj_ptr_list.cbegin(); iter != _obj_ptr_list.cend(); iter++)
	{
		Event temp_event;
		//������frequency mapѭ��
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
		//����cascade����¼���Pos = -1, which_event = -1
		temp_event.pos = -1;
		temp_event.which_event = -1;
		_event_list.push_back(temp_event);
	}

	//���¼��б�������He������¼�
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


//����event list������һ�������rand_num������һ��event������event list�е�λ��Ϊpos_i����event������������
//event_list.at(pos_i-1).rate < rand_num < event_list.at(pos_i).rate
int find_event(std::vector<Event> &_event_list, double _rate_sum)
{
	double rand_num = uni();
	double rand_rate = rand_num * _rate_sum;
	//binary search ���ֲ���
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


//����object��obj_ptr_list�е�λ�ã��Լ���Ҫִ�еķ�Ӧ��ִ�и÷�Ӧ
//�����Ǩ���¼���Ǩ��֮��Ӧ�����ײ�������ײ��_flag_refresh_event_list����Ϊtrue��׼���ؽ�event_list
//���ַ����¼�������Ҫ�ؽ�event_list
void carry_out_event(int _which_event, int _obj_pos, std::vector<Object *> &_obj_ptr_list, bool &_flag_refresh_event_list,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set)
{
	switch (_which_event)
	{
	//Ǩ��
	case 0:
	{
		//ִ��Ǩ�Ʋ���
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


//ִ��Ǩ���¼���Ǩ����֮��
//1.���߽�������pbc��Խ�磬���ƻغ����ڣ�abcԽ�磬��kill
//2.�����pbc��Ǩ����֮�󣬸���sink strength����Ƿ�ײ������
void carry_out_mig(int _obj_pos, std::vector<Object *> &_obj_ptr_list, bool &_flag_refresh_event_list,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set)
{

	_obj_ptr_list.at(_obj_pos)->mig(_set);

	//ִ����Ǩ���¼����Ƿ���Ҫcheck_reaction
	bool flag_check_reaction = true;
	//��������ڱ߽�����
	if (_set.boundary_condition == 2) {
		_obj_ptr_list.at(_obj_pos)->check_pbc(_set);
		_obj_ptr_list.at(_obj_pos)->refresh_grid_pos(_set);
		//����һ�����أ����������һsink����֮�󣬽���ر�
		bool check_sink_strength=true;
		//ͨ�����λ�ƣ��ж�Obj�Ƿ񴩹�����
		if (_set.grain_radius > 0)
		{
			//���λ�Ƴ��������뾶
			if (_obj_ptr_list.at(_obj_pos)->check_grain_radius(_set)) {
				if (_set.output_GB_absorption > 0) {
					output_react_1_0_GB(*(_obj_ptr_list.at(_obj_pos)), _set);
				}
				delete _obj_ptr_list.at(_obj_pos);
				_obj_ptr_list.erase(std::begin(_obj_ptr_list) + _obj_pos);
				//Ǩ�Ƶ�obj��ɱ���ˣ�����Ҫ����event list
				_flag_refresh_event_list = true;
				//Ǩ�Ƶ��Ǹ�Obj��ɱ���ˣ����Բ���check_reaction
				flag_check_reaction = false;
			}
		}
		//ͨ��sink strength���ж�Obj�Ƿ񴩹�����
		if (_set.ss_grain_radius > 0 && check_sink_strength)
		{
			int type = _obj_ptr_list.at(_obj_pos)->type;
			//���Obj�ǿ�λ����ô��һ�����������Ƿ�С��p_vac_meet_GB�����С�ڣ�kill��obj
			if (type == 1)
			{
				if (uni() < _set.p_vac_meet_GB)
				{
					if (_set.output_GB_absorption > 0) {
						output_react_1_0_GB(*(_obj_ptr_list.at(_obj_pos)), _set);
					}
					delete _obj_ptr_list.at(_obj_pos);
					_obj_ptr_list.erase(std::begin(_obj_ptr_list) + _obj_pos);
					//Ǩ�Ƶ�obj��ɱ���ˣ�����Ҫ����event list
					_flag_refresh_event_list = true;
					//Ǩ�Ƶ��Ǹ�Obj��ɱ���ˣ����Բ���check_reaction
					flag_check_reaction = false;
					//Ǩ�Ƶ��Ǹ�Obj��ɱ���ˣ����Բ������ж��Ƿ�����sink����
					check_sink_strength = false;
				}
			}
			//���Obj��sia����ô��һ�����������Ƿ�С��p_sia_meet_GB�����С�ڣ�kill��obj
			else if (type == 2)
			{
				if (uni() < _set.p_sia_meet_GB)
				{
					if (_set.output_GB_absorption > 0) {
						output_react_1_0_GB(*(_obj_ptr_list.at(_obj_pos)), _set);
					}
					delete _obj_ptr_list.at(_obj_pos);
					_obj_ptr_list.erase(std::begin(_obj_ptr_list) + _obj_pos);
					//Ǩ�Ƶ�obj��ɱ���ˣ�����Ҫ����event list
					_flag_refresh_event_list = true;
					//Ǩ�Ƶ��Ǹ�Obj��ɱ���ˣ����Բ���check_reaction
					flag_check_reaction = false;
					//Ǩ�Ƶ��Ǹ�Obj��ɱ���ˣ����Բ������ж��Ƿ�����sink����
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
				//Ǩ�Ƶ�obj��ɱ���ˣ�����Ҫ����event list
				_flag_refresh_event_list = true;
				//Ǩ�Ƶ��Ǹ�Obj��ɱ���ˣ����Բ���check_reaction
				flag_check_reaction = false;
				//Ǩ�Ƶ��Ǹ�Obj��ɱ���ˣ����Բ������ж��Ƿ�����sink����
				check_sink_strength = false;

			}

		}
		//Ǩ����֮�󣬿���������Obj�Ƿ�ײ��
		if (flag_check_reaction)
		{
			int pos2 = 0;
			if (check_dist_cross_bound_inner_loop(_obj_ptr_list, _obj_pos, pos2, _set))
			{
				carry_out_reaction(_obj_ptr_list, _obj_pos, pos2, *_obj_ptr_list.at(_obj_pos), *_obj_ptr_list.at(pos2),
					_vac_database, _sia_database, _fia_database, _vac_fia_database, _sia_fia_database, _set);
				//���ִ���˷�Ӧ��������obj�����仯�������Ҫ����event_list
				_flag_refresh_event_list = true;
			}
		}
	}
	//��������ձ߽�����
	else if (_set.boundary_condition == 1) {
		//���ײ���߽��ˣ��Ͱ�objȥ�������Ҳ���check_reaction
		if (_obj_ptr_list.at(_obj_pos)->check_cross_bound_xyz(_set)>0) {
			if (_set.output_out_bound > 0) {
				output_react_1_0(*(_obj_ptr_list.at(_obj_pos)), _set);
			}
			delete _obj_ptr_list.at(_obj_pos);
			_obj_ptr_list.erase(std::begin(_obj_ptr_list) + _obj_pos);
			//Ǩ�Ƶ�obj��ɱ���ˣ�����Ҫ����event list
			_flag_refresh_event_list = true;
			//Ǩ�Ƶ��Ǹ�Obj��ɱ���ˣ����Բ���check_reaction
			flag_check_reaction = false;
		}
		//���û��ײ���߽磬��Ӧ�ø�������Object���ڵ�����
		else {
			_obj_ptr_list.at(_obj_pos)->refresh_grid_pos(_set);
		}
		//Ǩ����֮�󣬿���������Obj�Ƿ�ײ��
		if (flag_check_reaction)
		{
			int pos2 = 0;
			if (check_dist_inner_loop_abc(_obj_ptr_list, _obj_pos, pos2, _set))
			{
				carry_out_reaction(_obj_ptr_list, _obj_pos, pos2, *_obj_ptr_list.at(_obj_pos), *_obj_ptr_list.at(pos2),
					_vac_database, _sia_database, _fia_database, _vac_fia_database, _sia_fia_database, _set);
				//���ִ���˷�Ӧ��������obj�����仯�������Ҫ����event_list
				_flag_refresh_event_list = true;
			}
		}
	}
	//�����XY���ڣ�Z����
	else if (_set.boundary_condition == 3)
	{
		if (_obj_ptr_list.at(_obj_pos)->check_cross_bound_z(_set))
		{
			if (_set.output_out_bound > 0) {
				output_react_1_0(*(_obj_ptr_list.at(_obj_pos)), _set);
			}
			delete _obj_ptr_list.at(_obj_pos);
			_obj_ptr_list.erase(std::begin(_obj_ptr_list) + _obj_pos);
			//Ǩ�Ƶ�obj��ɱ���ˣ�����Ҫ����event list
			_flag_refresh_event_list = true;
			//Ǩ�Ƶ��Ǹ�Obj��ɱ���ˣ����Բ���check_reaction
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
				//���ִ���˷�Ӧ��������obj�����仯�������Ҫ����event_list
				_flag_refresh_event_list = true;
			}
		}
	}

}


//������Ԫ�أ���λ�����Լ�϶��
//obj1�Ǳ������Ԫ�أ�����sizeΪ1
//��Ԫ��Obj2λ�ò��䣬obj1������䵽�����ϣ���뾶Ϊĸobj�İ뾶��3��
void carry_out_emit1(std::vector<Object *> &_obj_ptr_list, Object &_obj, int _obj_pos,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set)
{
	//���������µ�obj
	//ȷ��������������obj�Ĺ������
	int o1_type, o1_size1, o1_size2, o1_dir;
	double x1, y1, z1;
	int o2_type, o2_size1, o2_size2, o2_dir;
	double x2 = _obj.pos.at(0), y2 = _obj.pos.at(1), z2 = _obj.pos.at(2);
	//����һ��0-2pi֮����������ֵ��ȷ��������Ԫ�ص�����
	double angle1 = uni() * 2 * 3.1415926;
	double angle2 = uni() * 2 * 3.1415926;
	//��뾶Ϊĸobj�İ뾶��3
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
		//vac + fia�Ŵ�
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
		//sia + fia�Ŵ�
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

	//���������µ�Obj��obj1�Ǳ��³�����
	Object *ptr_obj1 = new Object(o1_type, x1, y1, z1, o1_size1, o1_size2, o1_dir, _vac_database, _sia_database,
		_fia_database, _vac_fia_database, _sia_fia_database, _set);
	Object *ptr_obj2 = new Object(o2_type, x2, y2, z2, o2_size1, o2_size2, o2_dir, _vac_database, _sia_database,
		_fia_database, _vac_fia_database, _sia_fia_database, _set);
	_obj_ptr_list.push_back(ptr_obj1);
	_obj_ptr_list.push_back(ptr_obj2);

	//���
	if (_set.output_emit_reaction > 0) {
		output_react_1_2(*_obj_ptr_list.at(_obj_pos), *ptr_obj1, *ptr_obj2, _set);
	}

	//ɾ���ϵ�Obj
	delete _obj_ptr_list.at(_obj_pos);
	_obj_ptr_list.erase(std::begin(_obj_ptr_list) + _obj_pos);

	//check_reaction�����³�����objѭ��һ�ּ��ɣ�pos1�Ǳ��³�����obj������
	int pos1 = _obj_ptr_list.size() - 2, pos2 = 0;
	if (check_dist_cross_bound_inner_loop(_obj_ptr_list, pos1, pos2,_set))
	{
		carry_out_reaction(_obj_ptr_list, pos1, pos2, *_obj_ptr_list.at(pos1), *_obj_ptr_list.at(pos2),
			_vac_database, _sia_database, _fia_database, _vac_fia_database, _sia_fia_database, _set);
	}
}


//���丱Ԫ�أ�fia��
//obj1�Ǳ������Ԫ�أ�����sizeΪ1
//��Ԫ��Obj2λ�ò��䣬obj1������䵽�����ϣ���뾶Ϊĸobj�İ뾶��2��
void carry_out_emit2(std::vector<Object *> &_obj_ptr_list, Object &_obj, int _obj_pos,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set)
{
	//���������µ�obj
	//ȷ��������������obj�Ĺ������
	int o1_type, o1_size1, o1_size2, o1_dir;
	double x1, y1, z1;
	int o2_type, o2_size1, o2_size2, o2_dir;
	double x2 = _obj.pos.at(0), y2 = _obj.pos.at(1), z2 = _obj.pos.at(2);
	//����һ��0-2pi֮����������ֵ��ȷ��������Ԫ�ص�����
	double angle1 = uni() * 2 * 3.1415926;
	double angle2 = uni() * 2 * 3.1415926;
	//��뾶Ϊĸobj�İ뾶��3
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
		//fia + vac�Ŵ�
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
		//fia + sia�Ŵ�
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

	//���������µ�Obj��obj1�Ǳ��³�����Ԫ��
	Object *ptr_obj1 = new Object(o1_type, x1, y1, z1, o1_size1, o1_size2, o1_dir, _vac_database, _sia_database,
		_fia_database, _vac_fia_database, _sia_fia_database, _set);
	Object *ptr_obj2 = new Object(o2_type, x2, y2, z2, o2_size1, o2_size2, o2_dir, _vac_database, _sia_database,
		_fia_database, _vac_fia_database, _sia_fia_database, _set);
	_obj_ptr_list.push_back(ptr_obj1);
	_obj_ptr_list.push_back(ptr_obj2);

	//���
	if (_set.output_emit_reaction > 0) {
		output_react_1_2(*_obj_ptr_list.at(_obj_pos), *ptr_obj1, *ptr_obj2, _set);
	}

	//ɾ���ϵ�Obj
	delete _obj_ptr_list.at(_obj_pos);
	_obj_ptr_list.erase(std::begin(_obj_ptr_list) + _obj_pos);

	//check_reaction�����³�����objѭ��һ�ּ���
	int pos1 = _obj_ptr_list.size() - 2, pos2 = 0;
	if (check_dist_cross_bound_inner_loop(_obj_ptr_list, pos1, pos2,_set))
	{
		carry_out_reaction(_obj_ptr_list, pos1, pos2, *_obj_ptr_list.at(pos1), *_obj_ptr_list.at(pos2),
			_vac_database, _sia_database, _fia_database, _vac_fia_database, _sia_fia_database, _set);
	}
}


//���ѡ��һ��cascade
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


//trap mutation�¼�
//He-V���ǵ����Է���trap mutation�¼���ʱ����һ��һ������SIAԭ�ӡ�
//He-SIA���ǵ����Է���trap mutation�¼���ʱ�򣬽�����SIA�Ŵ�һ���Ӷ������ȥ��

void carry_out_TM(std::vector<Object*>& _obj_ptr_list, Object& _obj, int _obj_pos,
	const Database& _vac_database, const Database& _sia_database, const Database& _fia_database,
	const Database& _vac_fia_database, const Database& _sia_fia_database, const Setting& _set)
{
	//���������µ�obj
	//ȷ��������������obj�Ĺ������
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
	//����һ��0-2pi֮����������ֵ��ȷ��������Ԫ�ص�����
	double angle1 = uni() * 2 * 3.1415926;
	double angle2 = uni() * 2 * 3.1415926;
	//��뾶Ϊĸobj�İ뾶��3
	double o1_r = pow(o1_size1 / pow(3, 0.5) / PI, 0.5) + 0.4813; //SIA �뾶
	double o2_r = pow((3 * o2_size1 / 8.0 / PI), 1 / 3.0) + 0.34; //V-He �뾶
	double r = o1_r + o2_r + 3;//double r = _obj.radius + 3;
	x1 = x2 + r * sin(angle1) * cos(angle2);
	y1 = y2 + r * sin(angle1) * sin(angle2);
	z1 = z2 + r * cos(angle1);


	//���������µ�Obj��obj1�Ǳ��³�����
	Object* ptr_obj1 = new Object(o1_type, x1, y1, z1, o1_size1, o1_size2, o1_dir, _vac_database, _sia_database,
		_fia_database, _vac_fia_database, _sia_fia_database, _set);
	Object* ptr_obj2 = new Object(o2_type, x2, y2, z2, o2_size1, o2_size2, o2_dir, _vac_database, _sia_database,
		_fia_database, _vac_fia_database, _sia_fia_database, _set);

	_obj_ptr_list.push_back(ptr_obj1);
	_obj_ptr_list.push_back(ptr_obj2);

	//���
	
	if (_set.output_trap_mutation > 0) {
		output_react_TM(*_obj_ptr_list.at(_obj_pos), *ptr_obj1, *ptr_obj2, _set);
	}
	

	//ɾ���ϵ�Obj
	delete _obj_ptr_list.at(_obj_pos);
	_obj_ptr_list.erase(std::begin(_obj_ptr_list) + _obj_pos);

	//check_reaction���������ɵ�obj����Ҫѭ��һ��
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


//He����ע��
void carry_out_He_insert(std::vector<Object*>& _obj_ptr_list,
	const Database& _vac_database, const Database& _sia_database, const Database& _fia_database,
	const Database& _vac_fia_database, const Database& _sia_fia_database, const Setting& _set)
{
	int o_type, o_size1, o_size2, o_dir;
	double x, y, z;

	//����һ��He
	o_type = 3; o_size1 = 0; o_size2 = 1, o_dir = 0;
	//λ�� x,y ������ȷֲ�  z��������һ���ֲ�����
	x = uni() * _set.box_length.at(0);
	y = uni() * _set.box_length.at(1);
	//z = (0.16/3.16 * 10) - (1.65/3.16 * 10) * log(1.0-uni());
	// z = (0.16 / _set.a0 * 10) - (1.65 / _set.a0 * 10) * log(1.0 - uni());
	z = uni() * _set.box_length.at(2);

	Object* ptr_obj = new Object(o_type, x, y, z, o_size1, o_size2, o_dir, _vac_database, _sia_database,
		_fia_database, _vac_fia_database, _sia_fia_database, _set);
	_obj_ptr_list.push_back(ptr_obj);

	////���Heֲ���¼�
	//if (_set.output_c_a_e_react > 0) {
	//	output_react_He_Insert(*ptr_obj, _set);
	//}

	//check_reaction���������ɵ�Heִ��һ�ּ���
	int pos1 = _obj_ptr_list.size() - 1, pos2 = 0;
	if (check_dist_cross_bound_inner_loop(_obj_ptr_list, pos1, pos2, _set))
	{
		carry_out_reaction(_obj_ptr_list, pos1, pos2, *_obj_ptr_list.at(pos1), *_obj_ptr_list.at(pos2),
			_vac_database, _sia_database, _fia_database, _vac_fia_database, _sia_fia_database, _set);
	}
}
