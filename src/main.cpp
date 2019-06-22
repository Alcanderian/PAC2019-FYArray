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

	I = Range(-1,ni+1);
	J = Range(-1,nj+1);
	K = Range(-1,nk+1);
	Range D(1,3);
	RDouble4D xfn (I,J,K,D,fortranArray);
	RDouble4D yfn (I,J,K,D,fortranArray);
	RDouble4D zfn (I,J,K,D,fortranArray);
	RDouble4D area(I,J,K,D,fortranArray);
	RDouble3D vol (I,J,K,  fortranArray);

	Range M(0,3);
	RDouble4D q_4d(I,J,K,M,fortranArray);
	RDouble4D dqdx_4d(I,J,K,M,fortranArray);
	RDouble4D dqdy_4d(I,J,K,M,fortranArray);
	RDouble4D dqdz_4d(I,J,K,M,fortranArray);

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
	printf("init finished\n");
	start=rdtsc();


	{
	Range I0(1,ni);
	Range J0(1,nj);
	Range K0(1,nk);

	Range I_(1,ni+1);
	Range J_(1,nj+1);
	Range K_(1,nk+1);

	const RDouble4D& xfn_const_ref = xfn;
	const RDouble4D& yfn_const_ref = yfn;
	const RDouble4D& zfn_const_ref = zfn;
	const RDouble4D& area_const_ref = area;
	const RDouble3D& vol_const_ref = vol;

	const RDouble4D& q_4d_const_ref = q_4d;

	RDouble4D xfn_area_mul(I,J,K,D,fortranArray);
	RDouble4D yfn_area_mul(I,J,K,D,fortranArray);
	RDouble4D zfn_area_mul(I,J,K,D,fortranArray);
	xfn_area_mul(I,J,K,D) = xfn_const_ref(I,J,K,D) * area_const_ref(I,J,K,D);
	yfn_area_mul(I,J,K,D) = yfn_const_ref(I,J,K,D) * area_const_ref(I,J,K,D);
	zfn_area_mul(I,J,K,D) = zfn_const_ref(I,J,K,D) * area_const_ref(I,J,K,D);

	RDouble4D* x_dim_sums[2][2][2];
	RDouble4D* y_dim_sums[2][2][2];
	RDouble4D* z_dim_sums[2][2][2];
	RDouble4D xmul_dim1_sum(I_,J_,K_,D,fortranArray);
	RDouble4D ymul_dim1_sum(I_,J_,K_,D,fortranArray);
	RDouble4D zmul_dim1_sum(I_,J_,K_,D,fortranArray);
	RDouble4D xmul_dim2_sum(I_,J_,K_,D,fortranArray);
	RDouble4D ymul_dim2_sum(I_,J_,K_,D,fortranArray);
	RDouble4D zmul_dim2_sum(I_,J_,K_,D,fortranArray);
	RDouble4D xmul_dim3_sum(I_,J_,K_,D,fortranArray);
	RDouble4D ymul_dim3_sum(I_,J_,K_,D,fortranArray);
	RDouble4D zmul_dim3_sum(I_,J_,K_,D,fortranArray);
	xmul_dim1_sum(I_,J_,K_,D) = xfn_area_mul(I_,J_,K_,D) + xfn_area_mul(I_-1,J_,K_,D);
	ymul_dim1_sum(I_,J_,K_,D) = yfn_area_mul(I_,J_,K_,D) + yfn_area_mul(I_-1,J_,K_,D);
	zmul_dim1_sum(I_,J_,K_,D) = zfn_area_mul(I_,J_,K_,D) + zfn_area_mul(I_-1,J_,K_,D);
	xmul_dim2_sum(I_,J_,K_,D) = xfn_area_mul(I_,J_,K_,D) + xfn_area_mul(I_,J_-1,K_,D);
	ymul_dim2_sum(I_,J_,K_,D) = yfn_area_mul(I_,J_,K_,D) + yfn_area_mul(I_,J_-1,K_,D);
	zmul_dim2_sum(I_,J_,K_,D) = zfn_area_mul(I_,J_,K_,D) + zfn_area_mul(I_,J_-1,K_,D);
	xmul_dim3_sum(I_,J_,K_,D) = xfn_area_mul(I_,J_,K_,D) + xfn_area_mul(I_,J_,K_-1,D);
	ymul_dim3_sum(I_,J_,K_,D) = yfn_area_mul(I_,J_,K_,D) + yfn_area_mul(I_,J_,K_-1,D);
	zmul_dim3_sum(I_,J_,K_,D) = zfn_area_mul(I_,J_,K_,D) + zfn_area_mul(I_,J_,K_-1,D);
	x_dim_sums[1][0][0] = &xmul_dim1_sum;
	x_dim_sums[0][1][0] = &xmul_dim2_sum;
	x_dim_sums[0][0][1] = &xmul_dim3_sum;
	y_dim_sums[1][0][0] = &ymul_dim1_sum;
	y_dim_sums[0][1][0] = &ymul_dim2_sum;
	y_dim_sums[0][0][1] = &ymul_dim3_sum;
	z_dim_sums[1][0][0] = &zmul_dim1_sum;
	z_dim_sums[0][1][0] = &zmul_dim2_sum;
	z_dim_sums[0][0][1] = &zmul_dim3_sum;

	RDouble4D q_4d_xy_conv(I_, J_, K_, M, fortranArray);
	RDouble4D q_4d_yz_conv(I_, J_, K_, M, fortranArray);
	RDouble4D q_4d_xz_conv(I_, J_, K_, M, fortranArray);
	// this can also be optimized
	q_4d_xy_conv(I_, J_, K_, M) = fourth * (q_4d_const_ref(I_, J_, K_, M) + q_4d_const_ref(I_ - 1, J_, K_, M) + q_4d_const_ref(I_, J_ - 1, K_, M) + q_4d_const_ref(I_ - 1, J_ - 1, K_, M));
	q_4d_yz_conv(I_, J_, K_, M) = fourth * (q_4d_const_ref(I_, J_, K_, M) + q_4d_const_ref(I_, J_ - 1, K_, M) + q_4d_const_ref(I_, J_, K_ - 1, M) + q_4d_const_ref(I_, J_ - 1, K_ - 1, M));
	q_4d_xz_conv(I_, J_, K_, M) = fourth * (q_4d_const_ref(I_, J_, K_, M) + q_4d_const_ref(I_ - 1, J_, K_, M) + q_4d_const_ref(I_, J_, K_ - 1, M) + q_4d_const_ref(I_ - 1, J_, K_ - 1, M));
	RDouble4D* q_4d_convs[2][2][2];
	q_4d_convs[1][1][0] = &q_4d_xy_conv;
	q_4d_convs[0][1][1] = &q_4d_yz_conv;
	q_4d_convs[1][0][1] = &q_4d_xz_conv;

	RDouble3D rev_vol_sum_dim1(I0, J0, K0, fortranArray);
	RDouble3D rev_vol_sum_dim2(I0, J0, K0, fortranArray);
	RDouble3D rev_vol_sum_dim3(I0, J0, K0, fortranArray);
	rev_vol_sum_dim1(I0, J0, K0) = 1.0 / (vol_const_ref(I0, J0, K0) + vol_const_ref(I0 - 1, J0, K0));
	rev_vol_sum_dim2(I0, J0, K0) = 1.0 / (vol_const_ref(I0, J0, K0) + vol_const_ref(I0, J0 - 1, K0));
	rev_vol_sum_dim3(I0, J0, K0) = 1.0 / (vol_const_ref(I0, J0, K0) + vol_const_ref(I0, J0, K0 - 1));
	RDouble3D* rev_vol_sums[2][2][2];
	rev_vol_sums[1][0][0] = &rev_vol_sum_dim1;
	rev_vol_sums[0][1][0] = &rev_vol_sum_dim2;
	rev_vol_sums[0][0][1] = &rev_vol_sum_dim3;

	RDouble3D worksx(I_,J_,K_,fortranArray);
	RDouble3D worksy(I_,J_,K_,fortranArray);
	RDouble3D worksz(I_,J_,K_,fortranArray);
	RDouble3D workqm(I_,J_,K_,fortranArray);

	for ( int nsurf = 1; nsurf <= THREE_D; ++ nsurf )
	{
		int ns1 = nsurf;
		int ns2 = (nsurf + 1 - 1) % THREE_D + 1;
		int ns3 = (nsurf + 2 - 1) % THREE_D + 1;

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

		// part 1, step 1
		worksx(I_,J_,K_) = (*x_dim_sums[il1][jl1][kl1])(I_, J_, K_, ns1);
		worksy(I_,J_,K_) = (*y_dim_sums[il1][jl1][kl1])(I_, J_, K_, ns1);
		worksz(I_,J_,K_) = (*z_dim_sums[il1][jl1][kl1])(I_, J_, K_, ns1);

		// part 1, step 2
		for ( int m = mst; m <= med; ++ m )
		{
			dqdx_4d(I_,J_,K_,m) = - worksx(I_,J_,K_) * q_4d_const_ref(I_-il1,J_-jl1,K_-kl1,m);  // identical
			dqdy_4d(I_,J_,K_,m) = - worksy(I_,J_,K_) * q_4d_const_ref(I_-il1,J_-jl1,K_-kl1,m);
			dqdz_4d(I_,J_,K_,m) = - worksz(I_,J_,K_) * q_4d_const_ref(I_-il1,J_-jl1,K_-kl1,m);
		}

		// part 1, step 3
		for ( int m = mst; m <= med; ++ m )
		{
			dqdx_4d(I_-il1,J_-jl1,K_-kl1,m) += worksx(I_,J_,K_) * q_4d_const_ref(I_-il1,J_-jl1,K_-kl1,m);  // identical
			dqdy_4d(I_-il1,J_-jl1,K_-kl1,m) += worksy(I_,J_,K_) * q_4d_const_ref(I_-il1,J_-jl1,K_-kl1,m);
			dqdz_4d(I_-il1,J_-jl1,K_-kl1,m) += worksz(I_,J_,K_) * q_4d_const_ref(I_-il1,J_-jl1,K_-kl1,m);
		}

		// part 2, step 1
		worksx(I_,J_,K_) = (*x_dim_sums[il1][jl1][kl1])(I_, J_, K_, ns2);
		worksy(I_,J_,K_) = (*y_dim_sums[il1][jl1][kl1])(I_, J_, K_, ns2);
		worksz(I_,J_,K_) = (*z_dim_sums[il1][jl1][kl1])(I_, J_, K_, ns2);

		for ( int m = mst; m <= med; ++ m )
		{
			workqm(I_,J_,K_) = (*q_4d_convs[il1+il2][jl1+jl2][kl1+kl2])(I_,J_,K_,m);

			// part 2, step 2
			dqdx_4d(I_,J_,K_,m) -= worksx(I_,J_,K_) * workqm(I_,J_,K_);
			dqdy_4d(I_,J_,K_,m) -= worksy(I_,J_,K_) * workqm(I_,J_,K_);
			dqdz_4d(I_,J_,K_,m) -= worksz(I_,J_,K_) * workqm(I_,J_,K_);

			// part 2, step 3
			dqdx_4d(I_-il2,J_-jl2,K_-kl2,m) += worksx(I_,J_,K_) * workqm(I_,J_,K_);
			dqdy_4d(I_-il2,J_-jl2,K_-kl2,m) += worksy(I_,J_,K_) * workqm(I_,J_,K_);
			dqdz_4d(I_-il2,J_-jl2,K_-kl2,m) += worksz(I_,J_,K_) * workqm(I_,J_,K_);
		}

		// part 3, step 1
		worksx(I_,J_,K_) = (*x_dim_sums[il1][jl1][kl1])(I_, J_, K_, ns3);
		worksy(I_,J_,K_) = (*y_dim_sums[il1][jl1][kl1])(I_, J_, K_, ns3);
		worksz(I_,J_,K_) = (*z_dim_sums[il1][jl1][kl1])(I_, J_, K_, ns3);

		for ( int m = mst; m <= med; ++ m )
		{
			workqm(I_,J_,K_) = (*q_4d_convs[il1+il3][jl1+jl3][kl1+kl3])(I_,J_,K_,m);

			// part 3, step 2
			dqdx_4d(I_,J_,K_,m) -= worksx(I_,J_,K_) * workqm(I_,J_,K_);
			dqdy_4d(I_,J_,K_,m) -= worksy(I_,J_,K_) * workqm(I_,J_,K_);
			dqdz_4d(I_,J_,K_,m) -= worksz(I_,J_,K_) * workqm(I_,J_,K_);

			// part 3, step 3
			dqdx_4d(I_-il3,J_-jl3,K_-kl3,m) += worksx(I_,J_,K_) * workqm(I_,J_,K_);
			dqdy_4d(I_-il3,J_-jl3,K_-kl3,m) += worksy(I_,J_,K_) * workqm(I_,J_,K_);
			dqdz_4d(I_-il3,J_-jl3,K_-kl3,m) += worksz(I_,J_,K_) * workqm(I_,J_,K_);
		}

		// part 4, step 2
		for ( int m = mst; m <= med; ++ m )
		{
			dqdx_4d(I0,J0,K0,m) *= (*rev_vol_sums[il1][jl1][kl1])(I0,J0,K0);
			dqdy_4d(I0,J0,K0,m) *= (*rev_vol_sums[il1][jl1][kl1])(I0,J0,K0);
			dqdz_4d(I0,J0,K0,m) *= (*rev_vol_sums[il1][jl1][kl1])(I0,J0,K0);
		}

	}
	}

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
		// #pragma omp parallel for
		// for ( int i = 0; i < ni; ++ i )
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
