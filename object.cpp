#include"object.h"
#include"database.h"

/****************************��̬��Ա������ʼ��****************************/
//��̬��Ա������Դ�ļ��ж��壬��ͷ�ļ��ж���ᱨ��
int Object::ID = 1;

Object::Object(int _type, double _x, double _y, double _z, int _size1, int _size2, int _dir,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set)
{
	id = ID; type = _type; size1 = _size1; size2 = _size2; dir = _dir;
	pos.at(0) = _x; pos.at(1) = _y; pos.at(2) = _z;
	
	p_meet_dislocation = 0.0;

	//���sia����fia_sia��dir����1-4֮�䣬���������Ϣ�����������һ��1-4֮���dir
	if ((type == 2 || type == 5) && (dir < 1 || dir > 4))
	{
		std::cout << "The dir of object " << id << " is wrong!\n";
		dir = (int)(uni() * 4) + 1;
	}

	if (type == 1 || type == 4) { type_plus = 1; }
	else if (type == 2 || type == 5) { type_plus = 2; }
	else { type_plus = 0; }

	displacement = { 0,0,0 };
	cross_bound = { 0,0,0 };

	//һ�㳬��database�ķ�Χ������Ϊ��obj��Ǩ�ơ������䡢��ת��
	mig_F = 6e12;
	mig_E = 10;
	emit_F1 = 6e12;
	emit_E1 = 10;
	emit_F2 = 6e12;
	emit_E2 = 10;
	rotate_F = 6e12;
	rotate_E = 10;
	TM_F = 6e12;
	TM_E = 10;

	//����Object�����ͣ�ѡ����Ӧ�����ݿ����Object������
	//����ߴ糬�����ݿⷶΧ����ʹ��Ĭ�ϲ���
	switch(type) 
	{
	case 1:
		//������������ݿ�ķ�Χ�������ݿ��ȡǨ�ơ����䡢��ת���뾶����
		if ((size_t)size1 <= _vac_database.row) {
			parameterization(_vac_database.find(size1));
		}
		//������������ݿ�ķ�Χ�������ṩĬ��ֵ
		else {
			/*
			// Rn = 3^0.5/4+(3*n/8/PI)^(1/3)-(3/8/PI)^(1/3)+0.01
			radius = pow((3 * size1 / 8.0 / PI), 1 / 3.0) - 0.0493;
			//sizeԽ�󣬷���ָǰ����Խ��������n^2/3
			emit_F1 = 6e12 * pow(size1, 0.666667);
			*/
			//��һ�׷����������˿�λ�İ뾶��V1-V1��4NN���Խ�ϣ����Ҵ��Ŵ�Ҳ�п��ܽ���
			radius = 1.3739 + pow(size1 * 0.11936, 1 / 3.0) - pow(0.11936, 1 / 3.0);
			if (size1 % 16 == 1) {
				emit_E1 = 1.18 + 1.66;
			}
			else {
				emit_E1 = 2.63 + 1.66;
			}
			//std::cout << id << "��object��Ĭ�ϲ������ɣ�����������" << type << std::endl;
		}
		break;
	case 2:
		//������������ݿ�ķ�Χ�������ݿ��ȡǨ�ơ����䡢��ת���뾶����
		if ((size_t)size1 <= _sia_database.row) {
			parameterization(_sia_database.find(size1));
		}
		//������������ݿ�ķ�Χ�������ṩĬ��ֵ
		else {
			mig_F = 1419300000000 * pow(size1, -0.365);
			mig_E = 10;
			
			radius = 1.15 * (1.3739 + pow(size1 * 0.11936, 1 / 3.0) - pow(0.11936, 1 / 3.0));
			/*
			// Rn =1/2+0.2+(B4/3^0.5/PI)^0.5-(1/3^0.5/PI)^0.5
			radius = pow(size1 / pow(3, 0.5) / PI, 0.5) + 0.2713;
			*/
			//std::cout << id << "��object��Ĭ�ϲ������ɣ�����������" << type << std::endl;
		}
		break;
	case 3:
		//������������ݿ�ķ�Χ�������ݿ��ȡǨ�ơ����䡢��ת���뾶����
		if ((size_t)size2 <= _fia_database.column) {
			parameterization(_fia_database.find(size2));
		}
		//������������ݿ�ķ�Χ�������ṩĬ��ֵ
		else {
			//�ʹ���λ�İ뾶��ͬ
			radius = 3 / 3.17 + pow(size1 * 0.11936, 1 / 3.0) - pow(0.11936, 1 / 3.0);
			//std::cout << id << "��object��Ĭ�ϲ������ɣ�����������" << type << std::endl;
			emit_F2 = 6e10; 
			emit_E2 = 1.99;
			TM_E = 10.0;
		}
		break;
	case 4:
		//������������ݿ�ķ�Χ�������ݿ��ȡǨ�ơ����䡢��ת���뾶����
		if ((size_t)size1 <= _vac_fia_database.row && (size_t)size2 <= _vac_fia_database.column){
			parameterization(_vac_fia_database.find(size1, size2));
		}
		//������������ݿ�ķ�Χ�������ṩĬ��ֵ
		else {
			//�ʹ���λ�İ뾶��ͬ
			radius =1.3739 + pow(size1 * 0.11936, 1 / 3.0) - pow(0.11936, 1 / 3.0);
			//std::cout << id << "��object��Ĭ�ϲ������ɣ�����������" << type << std::endl;
			/*
			if (double(size2) / size1 <= 3.0) {
				emit_E2 = 3.5;
			}
			else if (double(size2) / size1 <= 4.0) {
				emit_E2 = 2.1;
			}
			else {
				emit_E2 = 1.5;
			}
			*/
			mig_E = 10;
			double x = (double)size2 / (double)size1;
			if (size1 < 200)
			{
				if (x < 8) {
					double temp_E1_1 = -0.1497 * pow(x, 2) + 2.443 * x + 2.3739; //2.3739 = 1.66 + 0.7139
					double temp_E1_2 = -5.66904 * (pow((double)size1, 2 / 3) - pow((double)size1 - 1, 2 / 3)) + 4.89; // 4.89 = 3.23 + 1.66
					emit_E1 = temp_E1_1 > temp_E1_2 ? temp_E1_1 : temp_E1_2;
					emit_E2 = 0.0703 * pow(x, 2) - 0.9317 * x + 6.1644;
				}
				else {
					emit_E1 = 10;
					emit_E2 = 2;
				}
			}
			else {
				if (0.2 < x < 8)
				{
					emit_E1 = 4.29;
					emit_E2 = 0.0703 * pow(x, 2) - 0.9317 * x + 6.1644;
				}
				else if (x > 8) {
					emit_E1 = 10;
					emit_E2 = 2;
				}
				else {
					if (size1 % 16 == 1) {
						emit_E1 = 1.18 + 1.66;
					}
					else {
						emit_E1 = 2.63 + 1.66;
					}
					emit_E2 = 0.0703 * pow(x, 2) - 0.9317 * x + 6.1644;
				}
			}
			//����trap mutation�¼�
			TM_E = 10.0;
		}
		break;

	case 5:
		//������������ݿ�ķ�Χ�������ݿ��ȡǨ�ơ����䡢��ת���뾶����
		if ((size_t)size1 <= _sia_fia_database.row && (size_t)size2 <= _sia_fia_database.column)
		{
			parameterization(_sia_fia_database.find(size1, size2));
		}
		//������������ݿ�ķ�Χ�������ṩĬ��ֵ
		else {
			//R������N^1/2
			//R������N^1/3
			radius = 1.15 * (1.3739 + pow(size1 * 0.11936, 1 / 3.0) - pow(0.11936, 1 / 3.0));
			//radius = pow(0.119366 * size1, 1 / 3.0) + 0.38;
			//std::cout << id << "��object��Ĭ�ϲ������ɣ�����������" << type << std::endl;
			//����
			if (size1 > 35)
			{
				mig_E = 10;
				mig_F = 1419300000000 * pow(size1, -0.365) * pow(10, -size2);
			}
			//mig_E = 10;
			emit_E1 = 10;
			emit_E2 = 3;

			//trap mutation ���He/SIA����3����He����������5���Ϳ��Է���trap mutation�����ֵ�ǿ��Ա��
				TM_E = 10.0;

		}
		break;
	}

	/**/
	//��Obj��������վ��λ��
	//�����fia�������fia_siteѡ������Ӧ��ռλģʽ
	if (type == 3)
	{
		switch (_set.fia_site)
		{
		case 0:	//��bcc lattice site
			on_bcc_lattice(_set);
			break;
		case 1:	//��tis site
			on_tis_site(_set);
			break;
		default:
			std::cout << "Unspecified fia site\n";
			break;
		}
	}

	/*
	else {
		on_bcc_lattice(_set);
	}
	*/

	//�����fia������TIS-TISǨ�Ʒ�ʽ����ȷ����xyz�����ĸ����������ĸ��ǰ��������ĸ���1/4����
	if (type == 3 && _set.fia_site == 1)
	{
		for (int i = 0; i < 3; i++)
		{
			double fractional_part = pos.at(i) - floor(pos.at(i));
			if (fractional_part < 0.1)	//�����һλ��int���������λ�ô���tis_quarter_half_int_pos[2]
			{
				tis_quarter_half_int_pos.at(2) = i;
			}
			else if (fabs(fractional_part - 0.5) < 0.1)	//�����һλ��half���������λ�ô���tis_quarter_half_int_pos[1]
			{
				tis_quarter_half_int_pos.at(1) = i;
			}
			else {	//�����һλ��quarter���������λ�ô���tis_quarter_half_int_pos[0]
				tis_quarter_half_int_pos.at(0) = i;
			}
		}
	}
	else {
		tis_quarter_half_int_pos = { 0,0,0 };
	}

	check_pbc(_set);
	refresh_grid_pos(_set);

	cal_frequency(_set);
	cal_p_meet_dislocation(_set);
	ID++;
}

void Object::cal_p_meet_dislocation(const Setting &_set)
{
	double k_square;
	double dislocation_density = _set.dislocation_density;
	double dislocation_radius = _set.dislocation_radius;
	double a0 = _set.a0 * 1e-10;
	double a0_square = a0 * a0;

	//��3D��obj
	if (type_plus == 1)
	{
		double rou, rou_square;
		rou = (_set.dislocation_radius + radius)* a0 * sqrt(PI*dislocation_density);
		rou_square = rou * rou;
		k_square = 2 * PI * dislocation_density * (1 - rou_square) / (-log(rou) - 0.75 + 0.25*rou_square*(4 - rou_square));
		p_meet_dislocation = k_square * 0.125 * a0_square;
	}
	//��1D��obj
	else if (type_plus == 2) 
	{
		k_square = 2 * pow(PI*(dislocation_radius + radius) * a0 * dislocation_density, 2);
		p_meet_dislocation = k_square * 0.375 * a0_square;
	}
	else{
		p_meet_dislocation = 0.0;
	}
	/*
	if (type == 1) {
		std::cout << "��λ��λ�����յĸ��ʣ�" << p_meet_dislocation << std::endl;
	}
	else if (type == 2) {
		std::cout << "�Լ�϶��λ�����յĸ��ʣ�" << p_meet_dislocation << std::endl;
	}
	*/
}

//����Object��size1��size2����Object���в�����
void Object::parameterization(const Object_parameter & _parameter)
{
	mig_F = _parameter.mig_F;
	mig_E = _parameter.mig_E;
	emit_F1 = _parameter.emit_F1;
	emit_E1 = _parameter.emit_E1;
	emit_F2 = _parameter.emit_F2;
	emit_E2 = _parameter.emit_E2;
	rotate_F = _parameter.rotate_F;
	rotate_E = _parameter.rotate_E;
	TM_F = _parameter.TM_F;
	TM_E = _parameter.TM_E;
	radius = _parameter.radius;
}

//�Ƚ�Objectƽ�Ƶ�ԭ�㣬Ѱ�����������bcc��㣬�ٽ���ƽ�ƻ�ԭ����λ�ø���
void Object::on_bcc_lattice(const Setting &_set)
//3 steps
{
	std::array<double, 3> pos_ = { 0.0, 0.0, 0.0, };
	std::array<int, 3> times = { 0, 0, 0, };
	//step1 transform to the neighborhood of origin
	for (int i = 0; i < 3; i++)
	{
		times.at(i) = int(floor(pos.at(i)));
		pos_.at(i) = pos.at(i) - times.at(i);
	}
	//step2 find the 1NN lattice site
	std::array<double, 9> dist;
	for (int i = 0; i < 9; i++) 
	{
		dist.at(i) = cal_dist(pos_, _set.bcc_lattice_site.at(i));
	}
	//�ҳ�dist�����У���Сֵ������
	int min_pos = 0;
	for (int i = 1; i < 9; i++)
	{
		if (dist.at(i) < dist.at(min_pos))
		{
			min_pos = i;
		}
	}
	//step3 transform back to original location
	for (int i = 0; i < 3; i++)
	{
		pos.at(i) = times.at(i) + _set.bcc_lattice_site.at(min_pos).at(i);
	}
}

void Object::on_tis_site(const Setting &_set)
{
	std::array<double, 3> pos_ = { 0.0, 0.0, 0.0, };
	std::array<int, 3> times = { 0, 0, 0, };
	//step1 transform to the neighborhood of origin
	for (int i = 0; i < 3; i++)
	{
		times.at(i) = int(floor(pos.at(i)));
		pos_.at(i) = pos.at(i) - times.at(i);
	}
	//step2 find the 1NN tis site
	std::array<double, 24> dist;
	for (int i = 0; i < 24; i++)
	{
		dist.at(i) = cal_dist(pos_, _set.tis_site.at(i));
	}
	//�ҳ�dist�����У���Сֵ������
	int min_pos = 0;
	for (int i = 1; i < 24; i++)
	{
		if (dist.at(i) < dist.at(min_pos))
		{
			min_pos = i;
		}
	}
	//step3 transform back to original location
	for (int i = 0; i < 3; i++)
	{
		pos.at(i) = times.at(i) + _set.tis_site.at(min_pos).at(i);
	}
}

//��������ԭ�ӵ���ά���꣬��������֮��ľ���
double cal_dist(std::array<double, 3> a, std::array<double, 3> b)
{
	return pow(pow((a.at(0) - b.at(0)), 2) + pow((a.at(1) - b.at(1)), 2) + pow((a.at(2) - b.at(2)), 2), 0.5);
}

//����С��max_energy���¼����ż�¼��frequency map�У���������event list
void Object::cal_frequency(const Setting &_set)
{
	if (mig_E < _set.max_energy)
	{
		mig_Freq = mig_F * exp(-mig_E / _set.kb / _set.temperature);
		frequency[0] = mig_Freq / _set.mig_degeneration;
	}
	if (emit_E1 < _set.max_energy)
	{
		emit_Freq1 = emit_F1 * exp(-emit_E1 / _set.kb / _set.temperature);
		frequency[1] = emit_Freq1;
	}
	if (emit_E2 < _set.max_energy)
	{
		emit_Freq2 = emit_F2 * exp(-emit_E2 / _set.kb / _set.temperature);
		frequency[2] = emit_Freq2;
	}
	if (rotate_E < _set.max_energy)
	{
		rotate_Freq = rotate_F * exp(-rotate_E / _set.kb / _set.temperature);
		frequency[3] = rotate_Freq;
	}
	if (TM_E < _set.max_energy)
	{
		TM_Freq = TM_F * exp(-TM_E / _set.kb / _set.temperature);
		frequency[4] = TM_Freq;
	}
}

// ��cascade.txt��ÿһ������һ��object����̬�����ڴ棩����ָ������ָ�����obj_ptr_list��
void read_cascade(std::string _filename, std::vector<Object *> &_obj_ptr_list,
	const Database &_vac_database, const Database &_sia_database, const Database &_fia_database,
	const Database &_vac_fia_database, const Database &_sia_fia_database, const Setting &_set)
{
	std::ifstream inFile(_filename);
	if (inFile.is_open())
	{
		//Ԥ�Ȳ���һ��ƽ��λ�ƣ������������cascade��ʱ�����ƽ��
		double dx = uni()*_set.box_length.at(0);
		double dy = uni()*_set.box_length.at(1);
		double dz = uni()*_set.box_length.at(2);
		//Ԥ�Ȳ��������Ƕȣ������������cascade��ʱ�������ת
		double alpha = uni()*(2 * 3.14159265359);
		double beta = uni()*(2 * 3.14159265359);
		double gamma = uni()*(2 * 3.14159265359);
		int num = int(uni() * 3);
		
		std::array<std::array<double, 3>, 3> rx = { 1.0, 0.0, 0.0, 0.0};

		//Ԥ�Ȳ���������ת����
		//��x����תalpha
		std::array<std::array<double, 3>, 3> Rx = { 1.0, 0.0, 0.0, 0.0, cos(alpha), -sin(alpha), 0.0, sin(alpha), cos(alpha) };
		//��y����תbeta
		std::array<std::array<double, 3>, 3> Ry = { cos(beta), 0.0, sin(beta), 0.0, 1.0, 0.0, -sin(beta), 0.0, cos(beta) };
		//��z����תgamma
		std::array<std::array<double, 3>, 3> Rz = { cos(gamma), -sin(gamma), 0.0, sin(gamma), cos(gamma), 0.0, 0.0, 0.0, 1.0 };
		std::string str;
		while (std::getline(inFile, str))
		{
			if (!str.empty() && str.at(0) != '#') //skip the empty line and the line started with '#'
			{
				std::istringstream ss(str);
				int temp, type, size1, size2, dir;
				double x, y, z;
				ss >> temp >> type >> x >> y >> z >> size1 >> size2 >> dir;
				if (_set.step > 0 && _set.rand_translate_cascade)
				{
					x += dx;
					y += dy;
					z += dz;
				}
				if (_set.step > 0 && _set.rand_rotate_cascade > 0)
				{
					switch (num)
					{
					case 0:
						multiply(Rx, x, y, z);
						multiply(Ry, x, y, z);
						multiply(Rz, x, y, z);
						break;
					case 1:
						multiply(Ry, x, y, z);
						multiply(Rz, x, y, z);
						multiply(Rx, x, y, z);
						break;
					case 2:
						multiply(Rz, x, y, z);
						multiply(Rx, x, y, z);
						multiply(Ry, x, y, z);
						break;
					default:
						std::cout << "There is a problem in read_cascade function\n";
						break;
					}
					if (type == 2){ dir = (int)(uni() * 4) + 1; 
					}

				}
				Object *ptr_obj = new Object(type, x, y, z, size1, size2, dir, _vac_database, _sia_database,
					_fia_database, _vac_fia_database, _sia_fia_database, _set);
				_obj_ptr_list.push_back(ptr_obj);
			}
		}
		inFile.close();
		//std::cout << "The " << _filename << " is successfully read.\n";
	}
	else {
		std::cout << "Fail to open " << _filename << std::endl;
	}
}

//��ת������˿ռ�����
void multiply(std::array<std::array<double, 3>, 3> & R, double & x, double & y, double & z)
{
	double x_new = R[0][0] * x + R[0][1] * y + R[0][2] * z;
	double y_new = R[1][0] * x + R[1][1] * y + R[1][2] * z;
	double z_new = R[2][0] * x + R[2][1] * y + R[2][2] * z;

	x = x_new;
	y = y_new;
	z = z_new;
}


