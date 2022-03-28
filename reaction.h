#pragma once

#ifndef REACTION_H_
#define REACTION_H_

#include"object.h"
#include"output.h"	//����reaction֮����Ҫ����Output�еĺ���������"react.txt"


/***************************************����������***************************************/

//�����ڼ���pos1֮���obj����pos1ָ���obj���Ƿ��obj_list����������obj�Ľ�����
//����vector������object��vector�е�λ��pos1��pos2
//������vector����ѭ��pos1����ѭ��pos2��ͨ��check_dist_inner_loop����ʵ�֣�
//�Ƚ�����Object֮��ľ��룬�������object�������ţ�����false
//�����⵽��������object���ţ�������ֹѭ��������true���ҷ�������object��vector�е�����λ��pos1��pos2��ͨ������ʵ�֣�
//ֻ�е�pos1Խ��ʱ�Ż᷵��false���Ӷ�ȷ�����е�Object�������Ƚ�
bool check_dist(std::vector<Object *> &_obj_ptr_list, int &pos1, int &pos2);


//pos1�ǲ����ģ�ֻ��Pos2������Pos2Ҫ�����ã�pos1����Ҫ
//pos2�ӳ�ʼλ�ÿ�ʼ�����˵��������Pos2������pos1����Ƚ�����ָ���Object֮��ľ���
//���������ڰ뾶֮�ͣ������·��������򣬷���true
bool check_dist_inner_loop(const std::vector<Object *> &_obj_ptr_list, const int &_pos1, int &_pos2);

bool check_dist_inner_loop_abc(const std::vector<Object *> &_obj_ptr_list, const int &_pos1, int &_pos2, const Setting &_set);



//check_dist������������
//�жϿ����ڽ��
bool check_dist_cross_bound(std::vector<Object *> &_obj_ptr_list, int &pos1, int &pos2, const Setting &_set);


//check_dist_inner_loop������������
//�жϿ����ڽ��
bool check_dist_cross_bound_inner_loop(const std::vector<Object *> &_obj_ptr_list, const int &_pos1, int &_pos2, 
	const Setting &_set);



//����pos1��pos2ָ�������object
//���������ȫ������̬����һ���µ�oject,����ָ��ѹ��obj_ptr_list
//ɾ��������Ӧobject���ͷ����ǵ��ڴ棬����ָ�����ǵ�ָ���vector��ȥ��
//�ú�������ı�pos1,pos2��ֵ
void carry_out_reaction(std::vector<Object *> &_obj_ptr_list, int &_pos1, const int &_pos2, Object &_obj1, Object &_obj2,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set);


// ��������object��size��type1�������²�����object������
// -1����������û���ǵ������ͣ�0������ȫ���𣬲�������object
int judge_new_obj_type(int _type1, int _obj1_s1, int _type2, int _obj2_s1);


//����type, x, y, z, size1, size2, dir������һ���µ�object
void create_object(int _type, Object &_obj1, Object &_obj2, std::vector<Object *> &_obj_ptr_list,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set);


//������object��λ��
//Ӧ����size��Ȩƽ��������radius��Ȩƽ������Ŀǰ��total_size ��Ȩƽ��
void generate_new_pos(Object &_obj1, Object &_obj2, std::array<double, 3> &_new_pos);


#endif // !REACTION_H_
