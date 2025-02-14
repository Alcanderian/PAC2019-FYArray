#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstring>
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

int preccheck(RDouble4D dqdx_4d,RDouble4D dqdy_4d,RDouble4D dqdz_4d);

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
	RDouble3D x(I, J, K, fortranArray);
	RDouble3D y(I, J, K, fortranArray);
	RDouble3D z(I, J, K, fortranArray);
	for ( int k = -1; k <= nk+1; ++ k )
	{
		for ( int j = -1; j <= nj+1; ++ j )
		{
			for ( int i = -1; i <= ni+1; ++ i )
			{
				x(i,j,k) = i*0.1;
				y(i,j,k) = j*0.2;
				z(i,j,k) = k*0.3;
			}
		}
	}
	

	// 申请变量空间
	I = Range(-1,ni+1);
	J = Range(-1,nj+1);
        K = Range(-1,nk+1);
        Range D(1,3);
	RDouble4D xfn (I,J,K,D,fortranArray);  // 网格单元↙左下面法向，D为方向
	RDouble4D yfn (I,J,K,D,fortranArray);
	RDouble4D zfn (I,J,K,D,fortranArray);
	RDouble4D area(I,J,K,D,fortranArray);  // 网格单元↙左下面面积
	RDouble3D vol (I,J,K,  fortranArray);  // 网格单元体积

        Range M(0,3); // 4个变量：速度u、v、w，温度T
        RDouble4D q_4d(I,J,K,M,fortranArray); // 存储流场量，位置在单元中心
	RDouble4D dqdx_4d(I,J,K,M,fortranArray); // 存储流场量计算得到的梯度偏x
	RDouble4D dqdy_4d(I,J,K,M,fortranArray); // 存储流场量计算得到的梯度偏y
	RDouble4D dqdz_4d(I,J,K,M,fortranArray); // 存储流场量计算得到的梯度偏z

	// 计算网格单元几何数据 xfn、fn、zfn、area、vol
	// 速度u、v、w，温度T 流场变量赋值，存储在q_4d中，便于后面速度、温度界面梯度计算
	// 程序每执行一个迭代步，流场变量被更新。此处给初场值u=1.0，v=0.0，w=0.0，T=1.0
	for ( int k = -1; k <= nk+1; ++ k )
	{
		for ( int j = -1; j <= nj+1; ++ j )
		{
			for ( int i = -1; i <= ni+1; ++ i )
			{
				xfn(i,j,k,1) = 1.0;
				xfn(i,j,k,2) = 0.0;
				xfn(i,j,k,3) = 0.0;
				yfn(i,j,k,1) = 0.0;
				yfn(i,j,k,2) = 1.0;
				yfn(i,j,k,3) = 0.0;
				zfn(i,j,k,1) = 0.0;
				zfn(i,j,k,2) = 0.0;
				zfn(i,j,k,3) = 1.0;
				area(i,j,k,1) = 0.06;
				area(i,j,k,2) = 0.03;
				area(i,j,k,3) = 0.02;
				vol(i,j,k) = 0.006;
			}
		}
	}
	for ( int k = -1; k <= nk+1; ++ k )
	{
		for ( int j = -1; j <= nj+1; ++ j )
		{
			for ( int i = -1; i <= ni+1; ++ i )
			{
				q_4d(i,j,k,0) = (x(i,j,k) * x(i,j,k) + y(i,j,k)*y(i,j,k)- 1.3164) / 2.1547; // u = a*x*x+b*y*y
				q_4d(i,j,k,1) = (z(i,j,k)*z(i,j,k) - 0.2157 ) * 0.137; // v=c*z*z
				q_4d(i,j,k,2) = (2.0*x(i,j,k) +  1.737) / 3.14; // w=d*x
				q_4d(i,j,k,3) = x(i,j,k) + y(i,j,k) + 1.3765; // T = x + y
			}
		}
	}
	start=rdtsc();
	//以上为数据初始化部分，不可修改！
	// --------------------------------------------------------------------
	// 求解速度、温度在“单元界面”上的梯度，i、j、k三个方向依次求解
	// 在程序中是“耗时部分”，每一个迭代步都会求解，以下为未优化代码
	// 希望参赛队伍在理解该算法的基础上，实现更高效的界面梯度求解，提升程序执行效率
	// --------------------------------------------------------------------
	// 此处开始统计计算部分代码运行时间
	RDouble3D worksx(I,J,K,fortranArray);
	RDouble3D worksy(I,J,K,fortranArray);
	RDouble3D worksz(I,J,K,fortranArray);
	RDouble3D workqm(I,J,K,fortranArray);

	for ( int nsurf = 1; nsurf <= THREE_D; ++ nsurf )
	{
#define EXPANDED
#ifdef EXPANDED
#pragma omp parallel
{
#endif
#define LOC4D(i0, i1, i2, i3)	((i0) * s0 + (i1) * s1 + (i2) * s2 + (i3) * s3)
#define LOC3D(i0, i1, i2)		((i0) * s0 + (i1) * s1 + (i2) * s2)
		Range I(1,ni+1);
		Range J(1,nj+1);
		Range K(1,nk+1);

		int index[] = {1,2,3,1,2};

		int ns1 = nsurf;
		int ns2 = index[nsurf  ];
		int ns3 = index[nsurf+1];

		int il1 = 0;
		int il2 = 0;
		int il3 = 0;
		int jl1 = 0;
		int jl2 = 0;
		int jl3 = 0;
		int kl1 = 0;
		int kl2 = 0;
		int kl3 = 0;

		if ( nsurf == 1 )
		{
			il1 = 1;
			jl2 = 1;
			kl3 = 1;
		}
		else if ( nsurf == 2 )
		{
			jl1 = 1;
			kl2 = 1;
			il3 = 1;
		}
		else if ( nsurf == 3 )
		{
			kl1 = 1;
			il2 = 1;
			jl3 = 1;
		}

		Range M(mst,med);

#ifdef EXPANDED
		double* Pworksx = &worksx[0];
		double* Pworksy = &worksy[0];
		double* Pworksz = &worksz[0];
		double* Pworkqm = &workqm[0];
		double* Pdqdx_4d = &dqdx_4d[0];
		double* Pdqdy_4d = &dqdy_4d[0];
		double* Pdqdz_4d = &dqdz_4d[0];
		const double* Pxfn = &xfn[0];
		const double* Pyfn = &yfn[0];
		const double* Pzfn = &zfn[0];
		const double* Parea = &area[0];
		const double* Pq_4d = &q_4d[0];
		const double* Pvol = &vol[0];

		const int s0 = 1;
		const int s1 = s0 * (ni + 3);
		const int s2 = s1 * (nj + 3);
		const int s3 = s2 * (nk + 3);
#endif

#ifdef EXPANDED
		for ( int m = mst; m <= med; ++ m ) {
#pragma omp for
			for(int k = 1; k <= nk+1; ++k) {
				for(int j = 1; j <= nj+1; ++j) {
#pragma ivdep
					for(int i = 1; i <= ni+1; ++i) {
						Pdqdx_4d[LOC4D(i,j,k,m)] = 0.0;
						Pdqdx_4d[LOC4D(i,j,k,m)] = 0.0;
						Pdqdx_4d[LOC4D(i,j,k,m)] = 0.0;
					}
				}
			}
		}
#else
		dqdx_4d(I,J,K,M) = 0.0;
		dqdy_4d(I,J,K,M) = 0.0;
		dqdz_4d(I,J,K,M) = 0.0;
#endif

#ifdef EXPANDED
#pragma omp for
		for(int k = 1; k <= nk+1; ++k) {
			for(int j = 1; j <= nj+1; ++j) {
#pragma ivdep
				for(int i = 1; i <= ni+1; ++i) {
					Pworksx[LOC3D(i,j,k)] = \
						Pxfn[LOC4D(i,j,k,ns1)] * Parea[LOC4D(i,j,k,ns1)] + \
						Pxfn[LOC4D(i-il1,j-jl1,k-kl1,ns1)] * Parea[LOC4D(i-il1,j-jl1,k-kl1,ns1)];
					Pworksy[LOC3D(i,j,k)] = \
						Pyfn[LOC4D(i,j,k,ns1)] * Parea[LOC4D(i,j,k,ns1)] + \
						Pyfn[LOC4D(i-il1,j-jl1,k-kl1,ns1)] * Parea[LOC4D(i-il1,j-jl1,k-kl1,ns1)];
					Pworksz[LOC3D(i,j,k)] = \
						Pzfn[LOC4D(i,j,k,ns1)] * Parea[LOC4D(i,j,k,ns1)] + \
						Pzfn[LOC4D(i-il1,j-jl1,k-kl1,ns1)] * Parea[LOC4D(i-il1,j-jl1,k-kl1,ns1)];
				}
			}
		}
#else
		worksx(I,J,K) = xfn(I,J,K,ns1) * area(I,J,K,ns1) + xfn(I-il1,J-jl1,K-kl1,ns1) * area(I-il1,J-jl1,K-kl1,ns1);
		worksy(I,J,K) = yfn(I,J,K,ns1) * area(I,J,K,ns1) + yfn(I-il1,J-jl1,K-kl1,ns1) * area(I-il1,J-jl1,K-kl1,ns1);
		worksz(I,J,K) = zfn(I,J,K,ns1) * area(I,J,K,ns1) + zfn(I-il1,J-jl1,K-kl1,ns1) * area(I-il1,J-jl1,K-kl1,ns1);
#endif
		for ( int m = mst; m <= med; ++ m )
		{
#ifdef EXPANDED
#pragma omp for
			for(int k = 1; k <= nk+1; ++k) {
				for(int j = 1; j <= nj+1; ++j) {
#pragma ivdep
					for(int i = 1; i <= ni+1; ++i) {
						Pdqdx_4d[LOC4D(i,j,k,m)] = \
							- Pworksx[LOC3D(i,j,k)] * Pq_4d[LOC4D(i-il1,j-jl1,k-kl1,m)];
						Pdqdy_4d[LOC4D(i,j,k,m)] = \
							- Pworksy[LOC3D(i,j,k)] * Pq_4d[LOC4D(i-il1,j-jl1,k-kl1,m)];
						Pdqdz_4d[LOC4D(i,j,k,m)] = \
							- Pworksz[LOC3D(i,j,k)] * Pq_4d[LOC4D(i-il1,j-jl1,k-kl1,m)];
					}
				}
			}
#else
			dqdx_4d(I,J,K,m) = - worksx(I,J,K) * q_4d(I-il1,J-jl1,K-kl1,m);
			dqdy_4d(I,J,K,m) = - worksy(I,J,K) * q_4d(I-il1,J-jl1,K-kl1,m);
			dqdz_4d(I,J,K,m) = - worksz(I,J,K) * q_4d(I-il1,J-jl1,K-kl1,m);
#endif
		}

		for ( int m = mst; m <= med; ++ m )
		{
#ifdef EXPANDED
#pragma omp for
			for(int k = 1; k <= nk+1; ++k) {
				for(int j = 1; j <= nj+1; ++j) {
#pragma ivdep
					for(int i = 1; i <= ni+1; ++i) {
						Pdqdx_4d[LOC4D(i-il1,j-jl1,k-kl1,m)] += \
							Pworksx[LOC3D(i,j,k)] * Pq_4d[LOC4D(i-il1,j-jl1,k-kl1,m)];
						Pdqdy_4d[LOC4D(i-il1,j-jl1,k-kl1,m)] += \
							Pworksy[LOC3D(i,j,k)] * Pq_4d[LOC4D(i-il1,j-jl1,k-kl1,m)];
						Pdqdz_4d[LOC4D(i-il1,j-jl1,k-kl1,m)] += \
							Pworksz[LOC3D(i,j,k)] * Pq_4d[LOC4D(i-il1,j-jl1,k-kl1,m)];
					}
				}
			}
#else
			dqdx_4d(I-il1,J-jl1,K-kl1,m) += worksx(I,J,K) * q_4d(I-il1,J-jl1,K-kl1,m);
			dqdy_4d(I-il1,J-jl1,K-kl1,m) += worksy(I,J,K) * q_4d(I-il1,J-jl1,K-kl1,m);
			dqdz_4d(I-il1,J-jl1,K-kl1,m) += worksz(I,J,K) * q_4d(I-il1,J-jl1,K-kl1,m);
#endif
		}

		if ( ( nsurf != 2 ) || ( nDim != TWO_D ) )
		{
#ifdef EXPANDED
#pragma omp for
			for(int k = 1; k <= nk+1; ++k) {
				for(int j = 1; j <= nj+1; ++j) {
#pragma ivdep
					for(int i = 1; i <= ni+1; ++i) {
						Pworksx[LOC3D(i,j,k)] = \
							Pxfn[LOC4D(i,j,k,ns2)] * Parea[LOC4D(i,j,k,ns2)] + \
							Pxfn[LOC4D(i-il1,j-jl1,k-kl1,ns2)] * Parea[LOC4D(i-il1,j-jl1,k-kl1,ns2)];
						Pworksy[LOC3D(i,j,k)] = \
							Pyfn[LOC4D(i,j,k,ns2)] * Parea[LOC4D(i,j,k,ns2)] + \
							Pyfn[LOC4D(i-il1,j-jl1,k-kl1,ns2)] * Parea[LOC4D(i-il1,j-jl1,k-kl1,ns2)];
						Pworksz[LOC3D(i,j,k)] = \
							Pzfn[LOC4D(i,j,k,ns2)] * Parea[LOC4D(i,j,k,ns2)] + \
							Pzfn[LOC4D(i-il1,j-jl1,k-kl1,ns2)] * Parea[LOC4D(i-il1,j-jl1,k-kl1,ns2)];
					}
				}
			}
#else
			worksx(I,J,K) = xfn(I,J,K,ns2) * area(I,J,K,ns2) + xfn(I-il1,J-jl1,K-kl1,ns2) * area(I-il1,J-jl1,K-kl1,ns2);
			worksy(I,J,K) = yfn(I,J,K,ns2) * area(I,J,K,ns2) + yfn(I-il1,J-jl1,K-kl1,ns2) * area(I-il1,J-jl1,K-kl1,ns2);
			worksz(I,J,K) = zfn(I,J,K,ns2) * area(I,J,K,ns2) + zfn(I-il1,J-jl1,K-kl1,ns2) * area(I-il1,J-jl1,K-kl1,ns2);
#endif

			for ( int m = mst; m <= med; ++ m )
			{
#ifdef EXPANDED
#pragma omp for
				for(int k = 1; k <= nk+1; ++k) {
					for(int j = 1; j <= nj+1; ++j) {
#pragma ivdep
						for(int i = 1; i <= ni+1; ++i) {
							double temp = fourth * ( \
								Pq_4d[LOC4D(i,j,k,m)] + Pq_4d[LOC4D(i-il1,j-jl1,k-kl1,m)] + \
								Pq_4d[LOC4D(i-il2,j-jl2,k-kl2,m)] + Pq_4d[LOC4D(i-il1-il2,j-jl1-jl2,k-kl1-kl2,m)] );

							Pdqdx_4d[LOC4D(i,j,k,m)] -= \
								Pworksx[LOC3D(i,j,k)] * temp;
							Pdqdy_4d[LOC4D(i,j,k,m)] -= \
								Pworksy[LOC3D(i,j,k)] * temp;
							Pdqdz_4d[LOC4D(i,j,k,m)] -= \
								Pworksz[LOC3D(i,j,k)] * temp;

							Pdqdx_4d[LOC4D(i-il2,j-jl2,k-kl2,m)] += \
								Pworksx[LOC3D(i,j,k)] * temp;
							Pdqdy_4d[LOC4D(i-il2,j-jl2,k-kl2,m)] += \
								Pworksy[LOC3D(i,j,k)] * temp;
							Pdqdz_4d[LOC4D(i-il2,j-jl2,k-kl2,m)] += \
								Pworksz[LOC3D(i,j,k)] * temp;
						}
					}
				}
#else
				workqm(I,J,K) = fourth * ( q_4d(I,J,K,m) + q_4d(I-il1,J-jl1,K-kl1,m) + q_4d(I-il2,J-jl2,K-kl2,m) + q_4d(I-il1-il2,J-jl1-jl2,K-kl1-kl2,m) );

				dqdx_4d(I,J,K,m) -= worksx(I,J,K) * workqm(I,J,K);
				dqdy_4d(I,J,K,m) -= worksy(I,J,K) * workqm(I,J,K);
				dqdz_4d(I,J,K,m) -= worksz(I,J,K) * workqm(I,J,K);

				dqdx_4d(I-il2,J-jl2,K-kl2,m) += worksx(I,J,K) * workqm(I,J,K);
				dqdy_4d(I-il2,J-jl2,K-kl2,m) += worksy(I,J,K) * workqm(I,J,K);
				dqdz_4d(I-il2,J-jl2,K-kl2,m) += worksz(I,J,K) * workqm(I,J,K);
#endif
			}
		}

		if ( ( nsurf != 1 ) || ( nDim != TWO_D ) )
		{
#ifdef EXPANDED
#pragma omp for
			for(int k = 1; k <= nk+1; ++k) {
				for(int j = 1; j <= nj+1; ++j) {
#pragma ivdep
					for(int i = 1; i <= ni+1; ++i) {
						Pworksx[LOC3D(i,j,k)] = \
							Pxfn[LOC4D(i,j,k,ns3)] * Parea[LOC4D(i,j,k,ns3)] + \
							Pxfn[LOC4D(i-il1,j-jl1,k-kl1,ns3)] * Parea[LOC4D(i-il1,j-jl1,k-kl1,ns3)];
						Pworksy[LOC3D(i,j,k)] = \
							Pyfn[LOC4D(i,j,k,ns3)] * Parea[LOC4D(i,j,k,ns3)] + \
							Pyfn[LOC4D(i-il1,j-jl1,k-kl1,ns3)] * Parea[LOC4D(i-il1,j-jl1,k-kl1,ns3)];
						Pworksz[LOC3D(i,j,k)] = \
							Pzfn[LOC4D(i,j,k,ns3)] * Parea[LOC4D(i,j,k,ns3)] + \
							Pzfn[LOC4D(i-il1,j-jl1,k-kl1,ns3)] * Parea[LOC4D(i-il1,j-jl1,k-kl1,ns3)];
					}
				}
			}
#else
			worksx(I,J,K) = xfn(I,J,K,ns3) * area(I,J,K,ns3) + xfn(I-il1,J-jl1,K-kl1,ns3) * area(I-il1,J-jl1,K-kl1,ns3);
			worksy(I,J,K) = yfn(I,J,K,ns3) * area(I,J,K,ns3) + yfn(I-il1,J-jl1,K-kl1,ns3) * area(I-il1,J-jl1,K-kl1,ns3);
			worksz(I,J,K) = zfn(I,J,K,ns3) * area(I,J,K,ns3) + zfn(I-il1,J-jl1,K-kl1,ns3) * area(I-il1,J-jl1,K-kl1,ns3);
#endif

			for ( int m = mst; m <= med; ++ m )
			{
#ifdef EXPANDED
#pragma omp for
				for(int k = 1; k <= nk+1; ++k) {
					for(int j = 1; j <= nj+1; ++j) {
#pragma ivdep
						for(int i = 1; i <= ni+1; ++i) {
							double temp = fourth * ( \
								Pq_4d[LOC4D(i,j,k,m)] + Pq_4d[LOC4D(i-il1,j-jl1,k-kl1,m)] + \
								Pq_4d[LOC4D(i-il3,j-jl3,k-kl3,m)] + Pq_4d[LOC4D(i-il1-il3,j-jl1-jl3,k-kl1-kl3,m)] );

							Pdqdx_4d[LOC4D(i,j,k,m)] -= \
								Pworksx[LOC3D(i,j,k)] * temp;
							Pdqdy_4d[LOC4D(i,j,k,m)] -= \
								Pworksy[LOC3D(i,j,k)] * temp;
							Pdqdz_4d[LOC4D(i,j,k,m)] -= \
								Pworksz[LOC3D(i,j,k)] * temp;

							Pdqdx_4d[LOC4D(i-il3,j-jl3,k-kl3,m)] += \	
								Pworksx[LOC3D(i,j,k)] * temp;
							Pdqdy_4d[LOC4D(i-il3,j-jl3,k-kl3,m)] += \
								Pworksy[LOC3D(i,j,k)] * temp;
							Pdqdz_4d[LOC4D(i-il3,j-jl3,k-kl3,m)] += \
								Pworksz[LOC3D(i,j,k)] * temp;
						}
					}
				}
#else
				workqm(I,J,K) = fourth * ( q_4d(I,J,K,m) + q_4d(I-il1,J-jl1,K-kl1,m) + q_4d(I-il3,J-jl3,K-kl3,m) + q_4d(I-il1-il3,J-jl1-jl3,K-kl1-kl3,m) );

				dqdx_4d(I,J,K,m) -= worksx(I,J,K) * workqm(I,J,K);
				dqdy_4d(I,J,K,m) -= worksy(I,J,K) * workqm(I,J,K);
				dqdz_4d(I,J,K,m) -= worksz(I,J,K) * workqm(I,J,K);

				dqdx_4d(I-il3,J-jl3,K-kl3,m) += worksx(I,J,K) * workqm(I,J,K);
				dqdy_4d(I-il3,J-jl3,K-kl3,m) += worksy(I,J,K) * workqm(I,J,K);
				dqdz_4d(I-il3,J-jl3,K-kl3,m) += worksz(I,J,K) * workqm(I,J,K);
#endif
			}
		}

		Range I0(1,ni);
		Range J0(1,nj);
		Range K0(1,nk);

#ifdef EXPANDED
#else
		workqm(I0,J0,K0) = 1.0 / (  vol(I0, J0, K0) + vol(I0-il1, J0-jl1, K0-kl1) );
#endif

		for ( int m = mst; m <= med; ++ m )
		{
#ifdef EXPANDED
#pragma omp for
			for(int k = 1; k <= nk; ++k) {
				for(int j = 1; j <= nj; ++j) {
#pragma ivdep
					for(int i = 1; i <= ni; ++i) {
						double temp = 1.0 / (  Pvol[LOC3D(i,j,k)] + Pvol[LOC3D(i-il1, j-jl1, k-kl1)] );
						Pdqdx_4d[LOC4D(i,j,k,m)] *= temp;
						Pdqdy_4d[LOC4D(i,j,k,m)] *= temp;
						Pdqdz_4d[LOC4D(i,j,k,m)] *= temp;
					}
				}
			}
#else
			dqdx_4d(I0,J0,K0,m) *= workqm(I0,J0,K0);
			dqdy_4d(I0,J0,K0,m) *= workqm(I0,J0,K0);
			dqdz_4d(I0,J0,K0,m) *= workqm(I0,J0,K0);
#endif
		}
#ifdef EXPANDED
}
#endif
	// 该方向界面梯度值被计算出来后，会用于粘性通量计算，该值使用后下一方向会重新赋0计算
#undef LOC3D
#undef LOC4D
	}
	//----------------------------------------------------
	//以下为正确性对比部分，不可修改！
	//----------------------------------------------------
	end=rdtsc();
	elapsed= (end - start)/(F*Time);
	cout<<"The programe elapsed "<<elapsed<<setprecision(8)<<" s"<<endl;
	if(!preccheck(dqdx_4d,dqdy_4d,dqdz_4d))
		cout<<"Result check passed!"<<endl;
	return 0;
}

int preccheck(RDouble4D dqdx_4d,RDouble4D dqdy_4d,RDouble4D dqdz_4d)
{
	double tmp,real;
	ifstream file("check.txt",std::ofstream::binary);
	if ( !file )
	{
		cout << "Error opening check file! ";
		exit(1);
	}
    	//for ( int i = 0; i < ni; ++ i )
		for ( int i = 0; i < 2; ++ i )
	{
    		for ( int j = 0; j < nj; ++ j )
		{
			for ( int k = 0; k < nk; ++ k )
    			{
				for (int m = 0; m < 3; ++ m)
    				{
					file.read(reinterpret_cast<char*>(&tmp), sizeof(double));
					if(fabs(dqdx_4d(i,j,k,m) - tmp) > 1e-6)
					{
						real = dqdx_4d(i,j,k,m);
						cout<<"Precision check failed !"<<endl;
						cout<<"Your result is "<<setprecision(15)<<real<<endl;
						cout<<"The Standard result is "<<setprecision(15)<<tmp<<endl;
						cout<<"The wrong position is "<<endl;
						cout<<"i="<<i<<",j="<<j<<",k="<<k<<",m="<<m<<endl;
						exit(1);
					}

					file.read(reinterpret_cast<char*>(&tmp), sizeof(double));
					if(fabs(dqdy_4d(i,j,k,m) - tmp) > 1e-6)
					{
						real = dqdy_4d(i,j,k,m);
						cout<<"Precision check failed !"<<endl;
						cout<<"Your result is "<<setprecision(15)<<real<<endl;
						cout<<"The Standard result is "<<setprecision(15)<<tmp<<endl;
						cout<<"The wrong position is "<<endl;
						cout<<"i="<<i<<",j="<<j<<",k="<<k<<",m="<<m<<endl;
						exit(1);
					}

					file.read(reinterpret_cast<char*>(&tmp), sizeof(double));
					if(fabs(dqdz_4d(i,j,k,m) - tmp) >1e-6)
					{
						real = dqdz_4d(i,j,k,m);
						cout<<"Precision check failed !"<<endl;
						cout<<"Your result is "<<setprecision(15)<<real<<endl;
						cout<<"The Standard result is "<<setprecision(15)<<tmp<<endl;
						cout<<"The wrong position is "<<endl;
						cout<<"i="<<i<<",j="<<j<<",k="<<k<<",m="<<m<<endl;
						exit(1);
					}
				}
			}
		}
	}
	file.close();
	return 0;
}
