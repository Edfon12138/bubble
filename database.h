#pragma once

#include<string>
#include<vector>
#include<vector>
#include<sstream>
#include<fstream>

#include"object_parameter.h"

#ifndef DATABASE_H_
#define DATABASE_H_

// object�������ݿ⣬ÿһ���������һ��vector��vector�д洢�˸���object�����в���
class Database
{
public:
	// members
	size_t row, column;
	std::vector<Object_parameter> database;

	// functions
	Database() { std::cout << "The defaulted construct function is invoked for the Database.\n"; }
	Database(std::string _filename);
	~Database() {}

	// ��txt�ļ��ж�ȡ������ֵ
	void read_parameter(std::string filename);

	// ��pure object���õ�����find�������������ݿ��ж�Ӧ��λ��
	const Object_parameter &find(int) const;

	// ��complex object����˫����find�������������ݿ��ж�Ӧ��λ��
	const Object_parameter &find(int, int) const;
};




#endif // !DATABASE_H_
