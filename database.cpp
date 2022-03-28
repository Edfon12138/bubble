#include"database.h"

Database::Database(std::string _filename)
{
	std::ifstream inFile(_filename);
	if (inFile.is_open())
	{
		std::string str;
		//��ȡ���ݿ��row��column
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
	//��ȡ���ݿ�ľ������
	read_parameter(_filename);
}


// ��txt�ļ��ж�ȡ������ֵ
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


// ��pure object���õ�����find�������������ݿ��ж�Ӧ��λ��
const Object_parameter &Database::find(int _size) const
{
	size_t pos = _size - 1;	//m��n�е����ݿ⣬��AxBy�ǵ�(x-1)*n+y��
	return database.at(pos);
}


// ��complex object����˫����find�������������ݿ��ж�Ӧ��λ��
const Object_parameter &Database::find(int _size1, int _size2) const
{
	size_t pos = (_size1-1) * column + _size2 - 1;	// m��n�е����ݿ⣬��AxBy�ǵ�(x-1)*n+y��
	return database.at(pos);
}
