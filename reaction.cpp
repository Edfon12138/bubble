#include"reaction.h"


//适用于检查从pos1之后的obj（含pos1指向的obj）是否和obj_list中其他所有obj的结合情况
//输入vector和两个object在vector中的位置pos1和pos2
//对整个vector，外循环pos1，内循环pos2（通过check_dist_inner_loop函数实现）
//比较两个Object之间的距离，如果所有object都不挨着，返回false
//如果监测到任意两个object挨着，立即终止循环，返回true，且返回两个object在vector中的最新位置pos1和pos2（通过引用实现）
//只有当pos1越界时才会返回false，从而确保所有的Object都两两比较
bool check_dist(std::vector<Object *> &_obj_ptr_list, int &pos1, int &pos2)
{
	//确保 _obj_ptr_list.at(_pos1)不越界
	for (; (size_t)pos1 < _obj_ptr_list.size(); pos1++)
	{
		if (check_dist_inner_loop(_obj_ptr_list, pos1, pos2)) {
			return true;
		};
		pos2 = 0;
	}
	return false;
}

//输入两个object在vector中的位置pos1和pos2，对pos2循环
//如果它们的grid_pos是相邻的，则检查它们的距离
//如果距离小于等于半径和，立即返回true
bool check_dist_inner_loop(const std::vector<Object *> &_obj_ptr_list, const int &_pos1, int &_pos2)
{
	auto iter1 = _obj_ptr_list.cbegin() + _pos1;
	for (auto iter2 = _obj_ptr_list.cbegin() + _pos2; iter2 != _obj_ptr_list.cend(); iter2++, _pos2++)
	{
		if (iter2 != iter1)
		{
			//首先判断它们所在的网格是否相邻，如果相邻，才比较距离是否小于半径之和
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

//吸收边界条件下，基于网格划分思想，判断两个obj是否接触
bool check_dist_inner_loop_abc(const std::vector<Object *> &_obj_ptr_list, const int &_pos1, int &_pos2, const Setting &_set)
{
	auto iter1 = _obj_ptr_list.cbegin() + _pos1;
	for (auto iter2 = _obj_ptr_list.cbegin() + _pos2; iter2 != _obj_ptr_list.cend(); iter2++, _pos2++)
	{
		if (iter2 != iter1)
		{
			//首先判断是否需要检查距离
			std::array<int, 3> delta_grid;
			delta_grid[0] = (*iter1)->grid_pos[0] - (*iter2)->grid_pos[0];
			delta_grid[1] = (*iter1)->grid_pos[1] - (*iter2)->grid_pos[1];
			delta_grid[2] = (*iter1)->grid_pos[2] - (*iter2)->grid_pos[2];
			//如果三边的grid之差小于2，则需要检查距离
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
	//确保 _obj_ptr_list.at(_pos1)不越界
	for (; (size_t)pos1 < _obj_ptr_list.size(); pos1++)
	{
		if (check_dist_cross_bound_inner_loop(_obj_ptr_list, pos1, pos2,_set)) {
			return true;
		};
		pos2 = 0;
	}
	return false;

}

//判断跨周期结合，老版本
/*
bool check_dist_cross_bound_inner_loop(const std::vector<Object *> &_obj_ptr_list, const int &_pos1, int &_pos2, 
	const Setting &_set)
{
	auto iter1 = _obj_ptr_list.cbegin() + _pos1;
	for (auto iter2 = _obj_ptr_list.cbegin() + _pos2; iter2 != _obj_ptr_list.cend(); iter2++, _pos2++)
	{
		if (iter2 != iter1)
		{
			//第一步，检查跨周期的结合，保证dr始终在-L/2~L/2之间
			//dr存放两个Obj的直接位移；dr_new存放两个obj在周期镜像间的最近位移
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

			//如果是跨周期的结合，则把Object1平移一下，并返回true
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

//判断跨周期结合，基于网格划分的新版本
bool check_dist_cross_bound_inner_loop(const std::vector<Object *> &_obj_ptr_list, const int &_pos1, int &_pos2,
	const Setting &_set)
{
	auto iter1 = _obj_ptr_list.cbegin() + _pos1;
	for (auto iter2 = _obj_ptr_list.cbegin() + _pos2; iter2 != _obj_ptr_list.cend(); iter2++, _pos2++)
	{
		if (iter2 != iter1)
		{
			//首先判断是否需要检查距离
			std::array<int, 3> delta_grid;
			delta_grid[0] = (*iter1)->grid_pos[0] - (*iter2)->grid_pos[0];
			delta_grid[1] = (*iter1)->grid_pos[1] - (*iter2)->grid_pos[1];
			delta_grid[2] = (*iter1)->grid_pos[2] - (*iter2)->grid_pos[2];
			//如果三边的grid之差小于2（说明直接相邻），或者等于box_grid_num（说明两者有可能跨周期结合），则需要检查距离
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

				//如果是跨周期的结合，则把Object1平移一下，并返回true
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

//输入pos1和pos2指向的两个object
//如果不是完全湮灭，则动态生成一个新的oject,并将指针压入obj_ptr_list
//删除两个反应object，释放它们的内存，并将指向它们的指针从vector中去除
//该函数不会改变pos1,pos2的值
void carry_out_reaction(std::vector<Object *> &_obj_ptr_list, int &_pos1, const int &_pos2, Object &_obj1, Object &_obj2,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set)
{
	//判断新产生的object的类型
	int type = judge_new_obj_type(_obj1.type, _obj1.size1, _obj2.type, _obj2.size1);
	//如果产生了新的Object
	if (type > 0)
	{
		create_object(type, _obj1, _obj2, _obj_ptr_list, _vac_database, _sia_database, _fia_database,
			_vac_fia_database, _sia_fia_database, _set);
		if (_set.output_combine_reaction > 0)
		{
			output_react_2_1(_obj1, _obj2, **(_obj_ptr_list.end() - 1), _set);
		}
	}
	//如果反应产物完全湮灭了，没有产生新的Object
	else if (type == 0) 
	{ 
		if (_set.output_annihilation_reaction > 0)
		{
			output_react_2_0(_obj1, _obj2, _set);
		}
	}
	else if (type == -1) { std::cout << "Unexpected reaction is detected!\n"; }

	//先杀掉Object，释放内存
	delete _obj_ptr_list.at(_pos1);
	delete _obj_ptr_list.at(_pos2);
	//再erase vector中相应位置的指针，先删后面的元素，再删前面的元素
	if (_pos2 > _pos1) {
		_obj_ptr_list.erase(std::begin(_obj_ptr_list) + _pos2);
		_obj_ptr_list.erase(std::begin(_obj_ptr_list) + _pos1);
	}
	else if (_pos2 < _pos1){
		_obj_ptr_list.erase(std::begin(_obj_ptr_list) + _pos1);
		_obj_ptr_list.erase(std::begin(_obj_ptr_list) + _pos2);
		//如果pos2在pos1前面，则删除了两个object后，下一个待循环的object的索引是pos1-1
		_pos1--;
	}
}

// 输入两个object的size和type1，返回新产生的object的类型
// -1代表生成了没考虑到的类型；0代表完全湮灭，不产生新object
int judge_new_obj_type(int _type1, int _obj1_s1, int _type2, int _obj2_s1)
{
	if (
		(_type1 == 1 && _type2 == 2 && _obj1_s1 == _obj2_s1) ||		//vac + sia，且size1相等
		(_type1 == 2 && _type2 == 1 && _obj1_s1 == _obj2_s1))		//sia + vac，且size1相等
		return 0;
	else if (
		(_type1 == 1 && _type2 == 1) ||								//vac + vac
		(_type1 == 1 && _type2 == 2 && _obj1_s1 > _obj2_s1) ||		//vac + sia，且obj1的size1大
		(_type1 == 2 && _type2 == 1 && _obj1_s1 < _obj2_s1))		//sia + vac，且obj2的size1大
		return 1;
	else if (
		(_type1 == 2 && _type2 == 2) ||								//sia + sia
		(_type1 == 1 && _type2 == 2 && _obj1_s1 < _obj2_s1) ||		//vac + sia，且obj2的size1大
		(_type1 == 2 && _type2 == 1 && _obj1_s1 > _obj2_s1))		//sia + vac，且obj1的size1大
		return 2;
	else if (
		(_type1 == 3 && _type2 == 3) ||								//fia + fia
		(_type1 == 1 && _type2 == 5 && _obj1_s1 == _obj2_s1) ||		//vac + sia_fia，且size1相等
		(_type1 == 5 && _type2 == 1 && _obj1_s1 == _obj2_s1) ||		//sia_fia + vac，且size1相等
		(_type1 == 2 && _type2 == 4 && _obj1_s1 == _obj2_s1) ||		//sia + vac_fia，且size1相等
		(_type1 == 4 && _type2 == 2 && _obj1_s1 == _obj2_s1) ||		//vac_fia + fia，且size1相等
		(_type1 == 4 && _type2 == 5 && _obj1_s1 == _obj2_s1) ||		//vac_fia + sia_fia，且size1相等
		(_type1 == 5 && _type2 == 4 && _obj1_s1 == _obj2_s1))		//sia_fia + vac_fia，且size1相等
		return 3;
	else if (
		(_type1 == 1 && _type2 == 3) ||								//vac + fia
		(_type1 == 3 && _type2 == 1) ||								//fia + vac
		(_type1 == 1 && _type2 == 4) ||								//vac + vac_fia
		(_type1 == 4 && _type2 == 1) ||								//vac_fia + vac
		(_type1 == 3 && _type2 == 4) ||								//fia + vac_fia
		(_type1 == 4 && _type2 == 3) ||								//vac_fia + fia
		(_type1 == 4 && _type2 == 4) ||								//vac_fia + vac_fia
		(_type1 == 1 && _type2 == 5 && _obj1_s1 > _obj2_s1) ||		//vac + sia_fia，且obj1的size1大
		(_type1 == 5 && _type2 == 1 && _obj1_s1 < _obj2_s1) ||		//sia_fia + vac，且obj2的size1大
		(_type1 == 2 && _type2 == 4 && _obj1_s1 < _obj2_s1) ||		//sia + vac_fia，且obj2的size1大
		(_type1 == 4 && _type2 == 2 && _obj1_s1 > _obj2_s1) ||		//vac_fia + sia，且obj1的size1大
		(_type1 == 4 && _type2 == 5 && _obj1_s1 > _obj2_s1) ||		//vac_fia + sia_fia，且obj1的size1大
		(_type1 == 5 && _type2 == 4 && _obj1_s1 < _obj2_s1))		//sia_fia + vac_fia，且obj2的size1大
		return 4;
	else if (
		(_type1 == 2 && _type2 == 3) ||								//sia + fia
		(_type1 == 3 && _type2 == 2) ||								//fia + sia
		(_type1 == 2 && _type2 == 5) ||								//sia + sia_fia
		(_type1 == 5 && _type2 == 2) ||								//sia_fia + sia
		(_type1 == 3 && _type2 == 5) ||								//fia + sia_fia
		(_type1 == 5 && _type2 == 3) ||								//sia_fia + fia
		(_type1 == 5 && _type2 == 5) ||								//sia_fia + sia_fia
		(_type1 == 1 && _type2 == 5 && _obj1_s1 < _obj2_s1) ||		//vac + sia_fia，且obj2的size1大
		(_type1 == 5 && _type2 == 1 && _obj1_s1 > _obj2_s1) ||		//sia_fia + vac，且obj1的size1大
		(_type1 == 2 && _type2 == 4 && _obj1_s1 > _obj2_s1) ||		//sia + vac_fia，且obj1的size1大
		(_type1 == 4 && _type2 == 2 && _obj1_s1 < _obj2_s1) ||		//vac_fia + sia，且obj2的size1大
		(_type1 == 4 && _type2 == 5 && _obj1_s1 < _obj2_s1) ||		//vac_fia + sia_fia，且obj2的size1大
		(_type1 == 5 && _type2 == 4 && _obj1_s1 > _obj2_s1))		//sia_fia + vac_fia，且obj1的size1大
		return 5;
	else
		return -1;
}

//输入两个撞上的object，创建一个新的object
//根据type, x, y, z, size1, size2, dir，创建一个新的object
void create_object(int _type, Object &_obj1, Object &_obj2, std::vector<Object *> &_obj_ptr_list,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set)
{
	std::array<double, 3> new_pos = { 0,0,0 };
	generate_new_pos(_obj1, _obj2, new_pos);
	int size1;
	//如果主元素类型相同，则size1相加
	if (_obj1.type_plus == _obj2.type_plus) {
		size1 = _obj1.size1 + _obj2.size1;
	}
	//如果主元素类型不同，则size1相减，并取绝对值
	else {
		size1 = abs(_obj1.size1 - _obj2.size1);
	}
	int size2 = _obj1.size2 + _obj2.size2;

	int dir = 0;
	//如果产生了sia或者sia_fia
	if (_type == 2 || _type == 5)
	{
		//如果是sia+sia，则新sia的方向与较大的那个相同
		if (_obj1.type == 2 && _obj2.type == 2)
		{
			if (_obj1.size1 > _obj2.size1)
			{
				dir = _obj1.dir;
			}
			else dir = _obj2.dir;
		}
		//否则，新dir和第一个包含sia的父obj的dir相同
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

//计算新object的位置
//如果包含空位的obj和包含自间隙的obj相遇，则新Obj在大obj的基础上，更远离小obj
//否则，新obj在大obj和小obj之间，目前用total_size 加权平均
void generate_new_pos(Object &_obj1, Object &_obj2, std::array<double, 3> &_new_pos)
{
	double dist = pow(pow((_obj1.pos.at(0) - _obj2.pos.at(0)), 2) + pow((_obj1.pos.at(1) - _obj2.pos.at(1)), 2) 
		+ pow((_obj1.pos.at(2) - _obj2.pos.at(2)), 2), 0.5);
	//如果两个obj的位置正好重叠，则新obj的位置与老obj相同
	if (dist < 1e-6)
	{
		for (int i = 0; i < 3; i++) {
			_new_pos.at(i) = _obj1.pos.at(i);
		}
	}
	else 
	{
		//空位型Obj与自间隙型obj相遇
		if ((_obj1.type_plus + _obj2.type_plus) == 3)
		{
			//当obj2的size=0时，obj1的位置不变；当obj2的size=obj1的size时，Obj1的位置应该平移obj1的radius
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
		//同类型obj相遇
		else {
			for (int i = 0; i < 3; i++)
			{
				_new_pos.at(i) = (_obj1.pos.at(i) * (_obj1.size1 + _obj1.size2) + _obj2.pos.at(i) * (_obj2.size1 + _obj2.size2)) /
					(_obj1.size1 + _obj1.size2 + _obj2.size1 + _obj2.size2);
			}
		}
	}
}
