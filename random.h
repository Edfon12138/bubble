#pragma once

#ifndef RANDOM_H_
#define RANDOM_H_

#define BIGCOUNT 10000000    /* how many to do without re-initializing */

unsigned long poww2(int j);

static unsigned long m1 = 32767;
static long int mb[607];
/*=
	{
	30788, 23052, 2053, 19346, 10646, 19427, 23975,
	19049, 10949, 19693, 29746, 26748, 2796, 23890,
	29168, 31924, 16499
	};
*/
static int mdig = 32;
static unsigned long m2 = 256;
static int i = 272;
static int j = 606;


double suni(unsigned long jseed);

//generate a random number in [0,1)
double uni();

/*
		poww2
		*               calculate 2**j
		*/



#endif // !RANDOM_H_
