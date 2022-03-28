#include"object_parameter.h"



Object_parameter::Object_parameter(int _size1,int _size2,
	double _mig_F, double _mig_E, double _emit_F1, double _emit_E1, double _emit_F2, double _emit_E2,
	double _rotate_F, double _rotate_E, double _TM_F, double _TM_E, double _radius)
{
	obj_size1 = _size1;		obj_size2 = _size2;
	mig_F = _mig_F;			mig_E = _mig_E;
	emit_F1 = _emit_F1;		emit_E1 = _emit_E1;
	emit_F2 = _emit_F2;		emit_E2 = _emit_E2;
	rotate_F = _rotate_F;	rotate_E = _rotate_E;
	TM_F = _TM_F;			TM_E = _TM_E;
	radius = _radius;
}

