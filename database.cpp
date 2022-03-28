#include"database.h"

Database::Database(std::string _filename)
{
	std::ifstream inFile(_filename);
	if (inFile.is_open())
	{
		std::string str;
		//读取数据库的row和column
		while (std::getline(inFile, str))
		{
			if (!str.empty() && str.at(0) == '!') //skip the empty line and the line started with '#'
			{
				std::istringstream ss(str);
				std::string temp;
				ss >> temp >> row >> column;
				break;
			}
		}
	}
	else {
		std::cout << "Fail to open " << _filename;
	}
	inFile.close();
	//读取数据库的具体参数
	read_parameter(_filename);
}


// 从txt文件中读取参数的值
void Database::read_parameter(std::string _filename)
{
	std::ifstream inFile(_filename);
	if (inFile.is_open())
	{
		std::string str;
		while (std::getline(inFile, str))
		{
			if (!str.empty() && str.at(0) != '#' && str.at(0) != '!') //skip the empty line and the line started with '#'/'!'
			{
				std::istringstream ss(str);
				int _size1, _size2;
				double _mig_F, _mig_E, _emit_F1, _emig_E1, _bind_E1, _emit_F2, _emig_E2, _bind_E2,
					_rotate_E, _rotate_F, _TM_F, _TM_E, _radius;

				ss >> _size1 >> _size2 >> _mig_F >> _mig_E >> _emit_F1 >> _emig_E1 >> _bind_E1
					>> _emit_F2 >> _emig_E2 >> _bind_E2 >> _rotate_F >> _rotate_E >> _TM_F >> _TM_E >> _radius;

				database.emplace_back(_size1, _size2, _mig_F, _mig_E, _emit_F1, _emig_E1 + _bind_E1,
					_emit_F2, _emig_E2 + _bind_E2, _rotate_F, _rotate_E, _TM_F, _TM_E, _radius);
			}
		}
		std::cout << "The parameter file " << _filename << " is succssefully read.\n";
		inFile.close();
	}
	else {
		std::cout << "Fail to open" << _filename << std::endl;
	}

}


// 对pure object，用单参数find函数找它在数据库中对应的位置
const Object_parameter &Database::find(int _size) const
{
	size_t pos = _size - 1;	//m行n列的数据库，则AxBy是第(x-1)*n+y个
	return database.at(pos);
}


// 对complex object，用双参数find函数找它在数据库中对应的位置
const Object_parameter &Database::find(int _size1, int _size2) const
{
	size_t pos = (_size1-1) * column + _size2 - 1;	// m行n列的数据库，则AxBy是第(x-1)*n+y个
	return database.at(pos);
}
