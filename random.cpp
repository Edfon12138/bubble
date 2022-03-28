#include"random.h"

double suni(unsigned long jseed)
{
	long int j0, j1, k0, k1;
	double uni(void);

	/*      printf(" suni %ld\n", jseed);*/
	m1 = poww2(mdig - 2) - 1;     /* avoid overflow if m1 is full size */
	m1 += m1;
	m1++;
	/* printf(" m1 %lu, m2 %lu, mdig %d, jseed %u\n", m1, m2, mdig, jseed); */
	m2 = poww2((int)(mdig / 2));
	/* printf(" m1 %lu, m2 %lu, mdig %d, jseed %u\n", m1, m2, mdig, jseed); */
	jseed %= m1;                    /* jseed should less than m1 */
	if ((jseed & 1) == 0)           /* jseed should be odd */
		jseed--;
	k0 = 9069 % m2;                 /* simple congruential generator */
	k1 = 9069 / m2;                 /* the fanciness avoids overflow */
	j0 = jseed % m2;
	j1 = jseed / m2;
	for (i = 0; i < 607; i++)
	{
		jseed = j0 * k0;
		j1 = (jseed / m2 + j0 * k1 + j1 * k0) % (m2 / 2);
		j0 = jseed % m2;
		mb[i] = j0 + m2 * j1;
		/* printf("%2d %10u\n", i, mb.at(i)); */
	}
	i = 272;
	j = 606;
	return uni();
}

double uni()
{
	long int k;
	static unsigned long count = 0;
	k = mb[i] - mb[j];
	if (k < 0)
		k += m1;
	/* printf(" In UNI -- k = %ld\n",k); */
	if (++count >= BIGCOUNT)
	{
		count = 0;
		suni(k);
	}
	else
	{
		mb[j] = k;
		if (--i < 0)
			i = 606;
		if (--j < 0)
			j = 606;
	}
	/* printf("%lf\n", (double)k/(double)m1); */
	/* putchar(k%2 ?'+':'-');*/
	double rand_num = (double)k / (double)m1;
	return rand_num;
}

unsigned long poww2(int j)
{
	unsigned long x = 1;
	/* printf("poww2--j= %d\n",j); */
	while (j--)
		x *= 2;
	/* printf("poww2--x= %lu\n",x); */
	return (unsigned long)x;
}
