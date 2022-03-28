#pragma once

#ifndef EVENT_H_
#define EVENT_H_

#include<cmath>

#include"object.h"
#include"reaction.h"
#include"random.h"

class Event
{
public:
	/***************************************members***************************************/
	double rate;
	int pos;
	short which_event;
};


/***************************************函数声明区***************************************/

//输入obj_ptr_list和空的event_list，输出装满所有可能发生的事件的event list，并且返回总rate
//event_list中的每个元素包含三个信息：
//1. 截止到该event的rate总和
//2. 执行该event的object在obj_ptr_list的位置
//3. 该object应该执行哪个事件（mig~0, emit1~1, emit2~2 and rotate~3）
double build_event_list(std::vector<Object *> &_obj_ptr_list, std::vector<Event> &_event_list, Setting _set);


//输入event list，产生一个随机数rand_num，返回一个event，它在event list中的位置为pos_i，该event满足如下条件
//event_list.at(pos_i-1).rate < rand_num < event_list.at(pos_i).rate
int find_event(std::vector<Event> &_event_list, double _rate_sum);


//输入object在obj_ptr_list中的位置，以及它要执行的反应，执行该反应
//如果是迁移事件，迁移之后，应检查碰撞，如果碰撞，_flag_refresh_event_list重置为true，准备重建event_list
//两种发射事件，均需要重建event_list
void carry_out_event(int _which_event, int _obj_pos, std::vector<Object *> &_obj_ptr_list, bool &_flag_refresh_event_list,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set);

//执行迁移事件，检查边界
//如果有obj发生变化（产生或者消失），则将_flag_refresh_event_list重置为true
void carry_out_mig(int _obj_pos, std::vector<Object *> &_obj_ptr_list, bool &_flag_refresh_event_list,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set);

//obj1是被发射的元素，它的size为1
//大元素Obj2位置不变，obj1随机发射到球面上（球半径为母obj的半径加2）
void carry_out_emit1(std::vector<Object *> &_obj_ptr_list, Object &_obj, int _obj_pos,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set);

//obj1是被发射的元素，它的size为1
//大元素Obj2位置不变，obj1随机发射到球面上（球半径为母obj的半径加2）
void carry_out_emit2(std::vector<Object *> &_obj_ptr_list, Object &_obj, int _obj_pos,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set);


//随机选择一个cascade
std::string choose_cascade_to_insert(Setting &_set);


//trap_mutation
void carry_out_TM(std::vector<Object*>& _obj_ptr_list, Object& _obj, int _obj_pos,
	const Database& _vac_database, const Database& _sia_database, const Database& _fia_database,
	const Database& _vac_fia_database, const Database& _sia_fia_database, const Setting& _set);


void carry_out_He_insert(std::vector<Object*>& _obj_ptr_list,
	const Database& _vac_database, const Database& _sia_database, const Database& _fia_database,
	const Database& _vac_fia_database, const Database& _sia_fia_database, const Setting& _set);

#endif // !EVENT_H_
