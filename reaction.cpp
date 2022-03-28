#include"reaction.h"


//�����ڼ���pos1֮���obj����pos1ָ���obj���Ƿ��obj_list����������obj�Ľ�����
//����vector������object��vector�е�λ��pos1��pos2
//������vector����ѭ��pos1����ѭ��pos2��ͨ��check_dist_inner_loop����ʵ�֣�
//�Ƚ�����Object֮��ľ��룬�������object�������ţ�����false
//�����⵽��������object���ţ�������ֹѭ��������true���ҷ�������object��vector�е�����λ��pos1��pos2��ͨ������ʵ�֣�
//ֻ�е�pos1Խ��ʱ�Ż᷵��false���Ӷ�ȷ�����е�Object�������Ƚ�
bool check_dist(std::vector<Object *> &_obj_ptr_list, int &pos1, int &pos2)
{
	//ȷ�� _obj_ptr_list.at(_pos1)��Խ��
	for (; (size_t)pos1 < _obj_ptr_list.size(); pos1++)
	{
		if (check_dist_inner_loop(_obj_ptr_list, pos1, pos2)) {
			return true;
		};
		pos2 = 0;
	}
	return false;
}

//��������object��vector�е�λ��pos1��pos2����pos2ѭ��
//������ǵ�grid_pos�����ڵģ��������ǵľ���
//�������С�ڵ��ڰ뾶�ͣ���������true
bool check_dist_inner_loop(const std::vector<Object *> &_obj_ptr_list, const int &_pos1, int &_pos2)
{
	auto iter1 = _obj_ptr_list.cbegin() + _pos1;
	for (auto iter2 = _obj_ptr_list.cbegin() + _pos2; iter2 != _obj_ptr_list.cend(); iter2++, _pos2++)
	{
		if (iter2 != iter1)
		{
			//�����ж��������ڵ������Ƿ����ڣ�������ڣ��űȽϾ����Ƿ�С�ڰ뾶֮��
			if (abs((*iter1)->grid_pos[0] - (*iter2)->grid_pos[0]) < 2 && abs((*iter1)->grid_pos[1] - (*iter2)->grid_pos[1]) < 2
				&& abs((*iter1)->grid_pos[2] - (*iter2)->grid_pos[2]) < 2)
			{
				double dist = cal_dist((*iter1)->pos, (*iter2)->pos);
				double total_radius = (*iter1)->radius + (*iter2)->radius;
				if (dist <= total_radius) { return true; }
			}
		}
	}
	return false;
}

//���ձ߽������£��������񻮷�˼�룬�ж�����obj�Ƿ�Ӵ�
bool check_dist_inner_loop_abc(const std::vector<Object *> &_obj_ptr_list, const int &_pos1, int &_pos2, const Setting &_set)
{
	auto iter1 = _obj_ptr_list.cbegin() + _pos1;
	for (auto iter2 = _obj_ptr_list.cbegin() + _pos2; iter2 != _obj_ptr_list.cend(); iter2++, _pos2++)
	{
		if (iter2 != iter1)
		{
			//�����ж��Ƿ���Ҫ������
			std::array<int, 3> delta_grid;
			delta_grid[0] = (*iter1)->grid_pos[0] - (*iter2)->grid_pos[0];
			delta_grid[1] = (*iter1)->grid_pos[1] - (*iter2)->grid_pos[1];
			delta_grid[2] = (*iter1)->grid_pos[2] - (*iter2)->grid_pos[2];
			//������ߵ�grid֮��С��2������Ҫ������
			if (abs(delta_grid[0]) < 2 && abs(delta_grid[1]) < 2 && abs(delta_grid[2]) < 2)
			{
				std::array<double, 3> dr;
				dr[0] = (*iter1)->pos[0] - (*iter2)->pos[0];
				dr[1] = (*iter1)->pos[1] - (*iter2)->pos[1];
				dr[2] = (*iter1)->pos[2] - (*iter2)->pos[2];
				double dist = pow(pow(dr[0], 2) + pow(dr[1], 2) + pow(dr[2], 2), 0.5);
				double total_radius = (*iter1)->radius + (*iter2)->radius;
				if (dist < total_radius) {
					return true;
				}
			}
		}
	}
	return false;
}


bool check_dist_cross_bound(std::vector<Object *> &_obj_ptr_list, int &pos1, int &pos2, const Setting &_set)
{
	//ȷ�� _obj_ptr_list.at(_pos1)��Խ��
	for (; (size_t)pos1 < _obj_ptr_list.size(); pos1++)
	{
		if (check_dist_cross_bound_inner_loop(_obj_ptr_list, pos1, pos2,_set)) {
			return true;
		};
		pos2 = 0;
	}
	return false;

}

//�жϿ����ڽ�ϣ��ϰ汾
/*
bool check_dist_cross_bound_inner_loop(const std::vector<Object *> &_obj_ptr_list, const int &_pos1, int &_pos2, 
	const Setting &_set)
{
	auto iter1 = _obj_ptr_list.cbegin() + _pos1;
	for (auto iter2 = _obj_ptr_list.cbegin() + _pos2; iter2 != _obj_ptr_list.cend(); iter2++, _pos2++)
	{
		if (iter2 != iter1)
		{
			//��һ�����������ڵĽ�ϣ���֤drʼ����-L/2~L/2֮��
			//dr�������Obj��ֱ��λ�ƣ�dr_new�������obj�����ھ��������λ��
			std::array<double, 3> dr,dr_new;
			dr[0] = (*iter1)->pos[0] - (*iter2)->pos[0];
			dr[1] = (*iter1)->pos[1] - (*iter2)->pos[1];
			dr[2] = (*iter1)->pos[2] - (*iter2)->pos[2];
			for (int i = 0; i < 3; i++)
			{
				if (dr[i] > _set.box_length[i] / 2)
				{
					dr_new[i] = dr[i] - _set.box_length[i];
				}
				else if (dr[i] < _set.box_length[i] / -2)
				{
					dr_new[i] = dr[i] + _set.box_length[i];
				}
				else {
					dr_new[i] = dr[i];
				}
			}
			double dist = pow(pow(dr_new[0], 2) + pow(dr_new[1], 2) + pow(dr_new[2], 2), 0.5);
			double total_radius = (*iter1)->radius + (*iter2)->radius;

			//����ǿ����ڵĽ�ϣ����Object1ƽ��һ�£�������true
			if (dist <= total_radius)
			{
				for (int i = 0; i < 3; i++)
				{
					if (dr[i] > _set.box_length[i] / 2)
					{
						(*iter1)->pos[i] -= _set.box_length[i];
					}
					else if (dr[i] < _set.box_length[i] / -2)
					{
						(*iter1)->pos[i] += _set.box_length[i];
					}
				}
				return true;
			}
		}
	}
	return false;
}
*/

//�жϿ����ڽ�ϣ��������񻮷ֵ��°汾
bool check_dist_cross_bound_inner_loop(const std::vector<Object *> &_obj_ptr_list, const int &_pos1, int &_pos2,
	const Setting &_set)
{
	auto iter1 = _obj_ptr_list.cbegin() + _pos1;
	for (auto iter2 = _obj_ptr_list.cbegin() + _pos2; iter2 != _obj_ptr_list.cend(); iter2++, _pos2++)
	{
		if (iter2 != iter1)
		{
			//�����ж��Ƿ���Ҫ������
			std::array<int, 3> delta_grid;
			delta_grid[0] = (*iter1)->grid_pos[0] - (*iter2)->grid_pos[0];
			delta_grid[1] = (*iter1)->grid_pos[1] - (*iter2)->grid_pos[1];
			delta_grid[2] = (*iter1)->grid_pos[2] - (*iter2)->grid_pos[2];
			//������ߵ�grid֮��С��2��˵��ֱ�����ڣ������ߵ���box_grid_num��˵�������п��ܿ����ڽ�ϣ�������Ҫ������
			if ((abs(delta_grid[0]) < 2 || abs(delta_grid[0]) == _set.box_grid_num[0])
				&& (abs(delta_grid[1]) < 2 || abs(delta_grid[1]) == _set.box_grid_num[1])
				&& (abs(delta_grid[2]) < 2 || abs(delta_grid[2]) == _set.box_grid_num[2]))
			{
				std::array<double, 3> dr, dr_new;
				dr[0] = (*iter1)->pos[0] - (*iter2)->pos[0];
				dr[1] = (*iter1)->pos[1] - (*iter2)->pos[1];
				dr[2] = (*iter1)->pos[2] - (*iter2)->pos[2];
				for (int i = 0; i < 3; i++)
				{
					if (dr[i] > _set.box_length[i] / 2)
					{
						dr_new[i] = dr[i] - _set.box_length[i];
					}
					else if (dr[i] < _set.box_length[i] / -2)
					{
						dr_new[i] = dr[i] + _set.box_length[i];
					}
					else {
						dr_new[i] = dr[i];
					}
				}
				double dist = pow(pow(dr_new[0], 2) + pow(dr_new[1], 2) + pow(dr_new[2], 2), 0.5);
				double total_radius = (*iter1)->radius + (*iter2)->radius;

				//����ǿ����ڵĽ�ϣ����Object1ƽ��һ�£�������true
				if (dist <= total_radius)
				{
					for (int i = 0; i < 3; i++)
					{
						if (dr[i] > _set.box_length[i] / 2)
						{
							(*iter1)->pos[i] -= _set.box_length[i];
						}
						else if (dr[i] < _set.box_length[i] / -2)
						{
							(*iter1)->pos[i] += _set.box_length[i];
						}
					}
					return true;
				}
			}
		}
	}
	return false;
}

//����pos1��pos2ָ�������object
//���������ȫ������̬����һ���µ�oject,����ָ��ѹ��obj_ptr_list
//ɾ��������Ӧobject���ͷ����ǵ��ڴ棬����ָ�����ǵ�ָ���vector��ȥ��
//�ú�������ı�pos1,pos2��ֵ
void carry_out_reaction(std::vector<Object *> &_obj_ptr_list, int &_pos1, const int &_pos2, Object &_obj1, Object &_obj2,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set)
{
	//�ж��²�����object������
	int type = judge_new_obj_type(_obj1.type, _obj1.size1, _obj2.type, _obj2.size1);
	//����������µ�Object
	if (type > 0)
	{
		create_object(type, _obj1, _obj2, _obj_ptr_list, _vac_database, _sia_database, _fia_database,
			_vac_fia_database, _sia_fia_database, _set);
		if (_set.output_combine_reaction > 0)
		{
			output_react_2_1(_obj1, _obj2, **(_obj_ptr_list.end() - 1), _set);
		}
	}
	//�����Ӧ������ȫ�����ˣ�û�в����µ�Object
	else if (type == 0) 
	{ 
		if (_set.output_annihilation_reaction > 0)
		{
			output_react_2_0(_obj1, _obj2, _set);
		}
	}
	else if (type == -1) { std::cout << "Unexpected reaction is detected!\n"; }

	//��ɱ��Object���ͷ��ڴ�
	delete _obj_ptr_list.at(_pos1);
	delete _obj_ptr_list.at(_pos2);
	//��erase vector����Ӧλ�õ�ָ�룬��ɾ�����Ԫ�أ���ɾǰ���Ԫ��
	if (_pos2 > _pos1) {
		_obj_ptr_list.erase(std::begin(_obj_ptr_list) + _pos2);
		_obj_ptr_list.erase(std::begin(_obj_ptr_list) + _pos1);
	}
	else if (_pos2 < _pos1){
		_obj_ptr_list.erase(std::begin(_obj_ptr_list) + _pos1);
		_obj_ptr_list.erase(std::begin(_obj_ptr_list) + _pos2);
		//���pos2��pos1ǰ�棬��ɾ��������object����һ����ѭ����object��������pos1-1
		_pos1--;
	}
}

// ��������object��size��type1�������²�����object������
// -1����������û���ǵ������ͣ�0������ȫ���𣬲�������object
int judge_new_obj_type(int _type1, int _obj1_s1, int _type2, int _obj2_s1)
{
	if (
		(_type1 == 1 && _type2 == 2 && _obj1_s1 == _obj2_s1) ||		//vac + sia����size1���
		(_type1 == 2 && _type2 == 1 && _obj1_s1 == _obj2_s1))		//sia + vac����size1���
		return 0;
	else if (
		(_type1 == 1 && _type2 == 1) ||								//vac + vac
		(_type1 == 1 && _type2 == 2 && _obj1_s1 > _obj2_s1) ||		//vac + sia����obj1��size1��
		(_type1 == 2 && _type2 == 1 && _obj1_s1 < _obj2_s1))		//sia + vac����obj2��size1��
		return 1;
	else if (
		(_type1 == 2 && _type2 == 2) ||								//sia + sia
		(_type1 == 1 && _type2 == 2 && _obj1_s1 < _obj2_s1) ||		//vac + sia����obj2��size1��
		(_type1 == 2 && _type2 == 1 && _obj1_s1 > _obj2_s1))		//sia + vac����obj1��size1��
		return 2;
	else if (
		(_type1 == 3 && _type2 == 3) ||								//fia + fia
		(_type1 == 1 && _type2 == 5 && _obj1_s1 == _obj2_s1) ||		//vac + sia_fia����size1���
		(_type1 == 5 && _type2 == 1 && _obj1_s1 == _obj2_s1) ||		//sia_fia + vac����size1���
		(_type1 == 2 && _type2 == 4 && _obj1_s1 == _obj2_s1) ||		//sia + vac_fia����size1���
		(_type1 == 4 && _type2 == 2 && _obj1_s1 == _obj2_s1) ||		//vac_fia + fia����size1���
		(_type1 == 4 && _type2 == 5 && _obj1_s1 == _obj2_s1) ||		//vac_fia + sia_fia����size1���
		(_type1 == 5 && _type2 == 4 && _obj1_s1 == _obj2_s1))		//sia_fia + vac_fia����size1���
		return 3;
	else if (
		(_type1 == 1 && _type2 == 3) ||								//vac + fia
		(_type1 == 3 && _type2 == 1) ||								//fia + vac
		(_type1 == 1 && _type2 == 4) ||								//vac + vac_fia
		(_type1 == 4 && _type2 == 1) ||								//vac_fia + vac
		(_type1 == 3 && _type2 == 4) ||								//fia + vac_fia
		(_type1 == 4 && _type2 == 3) ||								//vac_fia + fia
		(_type1 == 4 && _type2 == 4) ||								//vac_fia + vac_fia
		(_type1 == 1 && _type2 == 5 && _obj1_s1 > _obj2_s1) ||		//vac + sia_fia����obj1��size1��
		(_type1 == 5 && _type2 == 1 && _obj1_s1 < _obj2_s1) ||		//sia_fia + vac����obj2��size1��
		(_type1 == 2 && _type2 == 4 && _obj1_s1 < _obj2_s1) ||		//sia + vac_fia����obj2��size1��
		(_type1 == 4 && _type2 == 2 && _obj1_s1 > _obj2_s1) ||		//vac_fia + sia����obj1��size1��
		(_type1 == 4 && _type2 == 5 && _obj1_s1 > _obj2_s1) ||		//vac_fia + sia_fia����obj1��size1��
		(_type1 == 5 && _type2 == 4 && _obj1_s1 < _obj2_s1))		//sia_fia + vac_fia����obj2��size1��
		return 4;
	else if (
		(_type1 == 2 && _type2 == 3) ||								//sia + fia
		(_type1 == 3 && _type2 == 2) ||								//fia + sia
		(_type1 == 2 && _type2 == 5) ||								//sia + sia_fia
		(_type1 == 5 && _type2 == 2) ||								//sia_fia + sia
		(_type1 == 3 && _type2 == 5) ||								//fia + sia_fia
		(_type1 == 5 && _type2 == 3) ||								//sia_fia + fia
		(_type1 == 5 && _type2 == 5) ||								//sia_fia + sia_fia
		(_type1 == 1 && _type2 == 5 && _obj1_s1 < _obj2_s1) ||		//vac + sia_fia����obj2��size1��
		(_type1 == 5 && _type2 == 1 && _obj1_s1 > _obj2_s1) ||		//sia_fia + vac����obj1��size1��
		(_type1 == 2 && _type2 == 4 && _obj1_s1 > _obj2_s1) ||		//sia + vac_fia����obj1��size1��
		(_type1 == 4 && _type2 == 2 && _obj1_s1 < _obj2_s1) ||		//vac_fia + sia����obj2��size1��
		(_type1 == 4 && _type2 == 5 && _obj1_s1 < _obj2_s1) ||		//vac_fia + sia_fia����obj2��size1��
		(_type1 == 5 && _type2 == 4 && _obj1_s1 > _obj2_s1))		//sia_fia + vac_fia����obj1��size1��
		return 5;
	else
		return -1;
}

//��������ײ�ϵ�object������һ���µ�object
//����type, x, y, z, size1, size2, dir������һ���µ�object
void create_object(int _type, Object &_obj1, Object &_obj2, std::vector<Object *> &_obj_ptr_list,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set)
{
	std::array<double, 3> new_pos = { 0,0,0 };
	generate_new_pos(_obj1, _obj2, new_pos);
	int size1;
	//�����Ԫ��������ͬ����size1���
	if (_obj1.type_plus == _obj2.type_plus) {
		size1 = _obj1.size1 + _obj2.size1;
	}
	//�����Ԫ�����Ͳ�ͬ����size1�������ȡ����ֵ
	else {
		size1 = abs(_obj1.size1 - _obj2.size1);
	}
	int size2 = _obj1.size2 + _obj2.size2;

	int dir = 0;
	//���������sia����sia_fia
	if (_type == 2 || _type == 5)
	{
		//�����sia+sia������sia�ķ�����ϴ���Ǹ���ͬ
		if (_obj1.type == 2 && _obj2.type == 2)
		{
			if (_obj1.size1 > _obj2.size1)
			{
				dir = _obj1.dir;
			}
			else dir = _obj2.dir;
		}
		//������dir�͵�һ������sia�ĸ�obj��dir��ͬ
		else
		{
			if (_obj1.type_plus == 2)
			{
				dir = _obj1.dir;
			}
			else
				dir = _obj2.dir;
		}
	}
	Object *ptr_obj = new Object(_type, new_pos.at(0), new_pos.at(1), new_pos.at(2), size1, size2, dir,
		_vac_database, _sia_database, _fia_database, _vac_fia_database, _sia_fia_database, _set);
	_obj_ptr_list.push_back(ptr_obj);
}

//������object��λ��
//���������λ��obj�Ͱ����Լ�϶��obj����������Obj�ڴ�obj�Ļ����ϣ���Զ��Сobj
//������obj�ڴ�obj��Сobj֮�䣬Ŀǰ��total_size ��Ȩƽ��
void generate_new_pos(Object &_obj1, Object &_obj2, std::array<double, 3> &_new_pos)
{
	double dist = pow(pow((_obj1.pos.at(0) - _obj2.pos.at(0)), 2) + pow((_obj1.pos.at(1) - _obj2.pos.at(1)), 2) 
		+ pow((_obj1.pos.at(2) - _obj2.pos.at(2)), 2), 0.5);
	//�������obj��λ�������ص�������obj��λ������obj��ͬ
	if (dist < 1e-6)
	{
		for (int i = 0; i < 3; i++) {
			_new_pos.at(i) = _obj1.pos.at(i);
		}
	}
	else 
	{
		//��λ��Obj���Լ�϶��obj����
		if ((_obj1.type_plus + _obj2.type_plus) == 3)
		{
			//��obj2��size=0ʱ��obj1��λ�ò��䣻��obj2��size=obj1��sizeʱ��Obj1��λ��Ӧ��ƽ��obj1��radius
			if (_obj1.size1 > _obj2.size1) {
				for (int i = 0; i < 3; i++) {
					_new_pos.at(i) = _obj1.pos.at(i) + (_obj1.pos.at(i) - _obj2.pos.at(i)) / dist *
						(double)_obj2.size1  / (double)_obj1.size1   * _obj1.radius;
					/*
					_new_pos.at(i) = _obj1.pos.at(i) +
						(_obj1.pos.at(i) - _obj2.pos.at(i)) / (_obj1.size1 + _obj2.size1) * _obj2.size1;
					*/
				}
			}
			else {
				for (int i = 0; i < 3; i++) {
					_new_pos.at(i) = _obj2.pos.at(i) + (_obj2.pos.at(i) - _obj1.pos.at(i)) / dist *
						(double)_obj1.size1  / (double)_obj2.size1  * _obj2.radius;
					/*
					_new_pos.at(i) = _obj2.pos.at(i) +
						(_obj2.pos.at(i) - _obj1.pos.at(i)) / (_obj1.size1 + _obj2.size1)*_obj1.size1;
						*/
				}
			}
		}
		//ͬ����obj����
		else {
			for (int i = 0; i < 3; i++)
			{
				_new_pos.at(i) = (_obj1.pos.at(i) * (_obj1.size1 + _obj1.size2) + _obj2.pos.at(i) * (_obj2.size1 + _obj2.size2)) /
					(_obj1.size1 + _obj1.size2 + _obj2.size1 + _obj2.size2);
			}
		}
	}
}
