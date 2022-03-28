#pragma once

#include<string>
#include<vector>
#include<vector>
#include<sstream>
#include<fstream>

#include"object_parameter.h"

#ifndef DATABASE_H_
#define DATABASE_H_

// object参数数据库，每一个对象包含一个vector，vector中存储了该类object的所有参数
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

	// 从txt文件中读取参数的值
	void read_parameter(std::string filename);

	// 对pure object，用单参数find函数找它在数据库中对应的位置
	const Object_parameter &find(int) const;

	// 对complex object，用双参数find函数找它在数据库中对应的位置
	const Object_parameter &find(int, int) const;
};




#endif // !DATABASE_H_
