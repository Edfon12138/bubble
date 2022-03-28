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


/***************************************����������***************************************/

//����obj_ptr_list�Ϳյ�event_list�����װ�����п��ܷ������¼���event list�����ҷ�����rate
//event_list�е�ÿ��Ԫ�ذ���������Ϣ��
//1. ��ֹ����event��rate�ܺ�
//2. ִ�и�event��object��obj_ptr_list��λ��
//3. ��objectӦ��ִ���ĸ��¼���mig~0, emit1~1, emit2~2 and rotate~3��
double build_event_list(std::vector<Object *> &_obj_ptr_list, std::vector<Event> &_event_list, Setting _set);


//����event list������һ�������rand_num������һ��event������event list�е�λ��Ϊpos_i����event������������
//event_list.at(pos_i-1).rate < rand_num < event_list.at(pos_i).rate
int find_event(std::vector<Event> &_event_list, double _rate_sum);


//����object��obj_ptr_list�е�λ�ã��Լ���Ҫִ�еķ�Ӧ��ִ�и÷�Ӧ
//�����Ǩ���¼���Ǩ��֮��Ӧ�����ײ�������ײ��_flag_refresh_event_list����Ϊtrue��׼���ؽ�event_list
//���ַ����¼�������Ҫ�ؽ�event_list
void carry_out_event(int _which_event, int _obj_pos, std::vector<Object *> &_obj_ptr_list, bool &_flag_refresh_event_list,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set);

//ִ��Ǩ���¼������߽�
//�����obj�����仯������������ʧ������_flag_refresh_event_list����Ϊtrue
void carry_out_mig(int _obj_pos, std::vector<Object *> &_obj_ptr_list, bool &_flag_refresh_event_list,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set);

//obj1�Ǳ������Ԫ�أ�����sizeΪ1
//��Ԫ��Obj2λ�ò��䣬obj1������䵽�����ϣ���뾶Ϊĸobj�İ뾶��2��
void carry_out_emit1(std::vector<Object *> &_obj_ptr_list, Object &_obj, int _obj_pos,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set);

//obj1�Ǳ������Ԫ�أ�����sizeΪ1
//��Ԫ��Obj2λ�ò��䣬obj1������䵽�����ϣ���뾶Ϊĸobj�İ뾶��2��
void carry_out_emit2(std::vector<Object *> &_obj_ptr_list, Object &_obj, int _obj_pos,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set);


//���ѡ��һ��cascade
std::string choose_cascade_to_insert(Setting &_set);


//trap_mutation
void carry_out_TM(std::vector<Object*>& _obj_ptr_list, Object& _obj, int _obj_pos,
	const Database& _vac_database, const Database& _sia_database, const Database& _fia_database,
	const Database& _vac_fia_database, const Database& _sia_fia_database, const Setting& _set);


void carry_out_He_insert(std::vector<Object*>& _obj_ptr_list,
	const Database& _vac_database, const Database& _sia_database, const Database& _fia_database,
	const Database& _vac_fia_database, const Database& _sia_fia_database, const Setting& _set);

#endif // !EVENT_H_
