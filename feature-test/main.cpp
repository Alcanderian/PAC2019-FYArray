#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "DataStruct_Array.h"
#define F 2.2E3
#define Time 1E6
using namespace std;
using namespace FYSPACE;

const int ONE_D   = 1;
const int TWO_D   = 2;
const int THREE_D = 3;
const int ni      = 500;
const int nj      = 400;
const int nk      = 300;

typedef double RDouble;
typedef FYArray<RDouble ,3> RDouble3D;
typedef FYArray<RDouble ,4> RDouble4D;

inline unsigned long long rdtsc(void)
{
	unsigned long hi = 0, lo = 0;

	__asm__ __volatile__ ("lfence;rdtsc" : "=a"(lo), "=d"(hi));

	return (((unsigned long long)lo))|(((unsigned long long)hi)<<32);
}

int main()
{
	double start,end,elapsed;
	const int nDim = THREE_D;
	const double fourth = 0.25;
	int mst = 0;
	int med = 3;


	Range I(-1,ni+1);
	Range J(-1,nj+1);
	Range K(-1,nk+1);
	// RDouble3D x(I, J, K, fortranArray);
	// RDouble3D y(I, J, K, fortranArray);
	// RDouble3D z(I, J, K, fortranArray);

	printf("[fortranArray]\n");
	RDouble3D ftest(3, 4, 1, fortranArray);
	int fnnn = 0;
	for(int i = 0; i < 3; ++i) {
		for(int j = 0; j < 4; ++j) {
			ftest(i, j, 0) = (double)(++fnnn);
		}
	}

	double *Pftest = &ftest[0]; // to get the correct pointer
	printf("shape: 3, 4, 1\n");
	printf("stride: %d, %d, %d\n", ftest.stride(0), ftest.stride(1), ftest.stride(2));
	printf("j * 3 + i: ");
	for(int i = 0; i < 3; ++i) {
		for(int j = 0; j < 4; ++j) {
			printf("%02.0lf ", Pftest[j * 3 + i]);
		}
		printf("| ");
	}
	printf("\n");
	printf("i * 4 + j: ");
	for(int i = 0; i < 3; ++i) {
		for(int j = 0; j < 4; ++j) {
			printf("%02.0lf ", Pftest[i * 4 + j]);
		}
		printf("| ");
	}
	printf("\n");

	printf("[cArray]\n");
	RDouble3D ctest(3, 4, 1);
	int cnnn = 0;
	for(int i = 0; i < 3; ++i) {
		for(int j = 0; j < 4; ++j) {
			ctest(i, j, 0) = (double)(++cnnn);
		}
	}

	double *Pctest = &ctest[0]; // to get the correct pointer
	printf("shape: 3, 4, 1\n");
	printf("stride: %d, %d, %d\n", ctest.stride(0), ctest.stride(1), ctest.stride(2));
	printf("j * 3 + i: ");
	for(int i = 0; i < 3; ++i) {
		for(int j = 0; j < 4; ++j) {
			printf("%02.0lf ", Pctest[j * 3 + i]);
		}
		printf("| ");
	}
	printf("\n");
	printf("i * 4 + j: ");
	for(int i = 0; i < 3; ++i) {
		for(int j = 0; j < 4; ++j) {
			printf("%02.0lf ", Pctest[i * 4 + j]);
		}
		printf("| ");
	}
	printf("\n");

	exit(0);
	return 0;
}
