#pragma once
#define GSL_DLL
#include "tanhFit.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include<assert.h>
#include <iostream>
#include <windows.h>
#include <omp.h>
#include<random>
#include "opmat.h"

//#include "nlinfit.h"

using namespace std;

class cut_pro_ransac
{
public:
	double fx(double x, double* a)
	{
		double a1 = a[0];
		double a2 = a[1];
		double a3 = a[2];
		double a4 = a[3];
		double a5 = a[4];
		double a6 = a[5];
		//double a6 = 0;
		double m = (a4 - x) / a5;
		double fz = (1 + a3 * m) * exp(m) - (1 + a6 * m) * exp(-m);
		double fm = exp(m) + exp(-m);
		double Yi = a1 * fz / fm + a2;
		return Yi;
	}
	void Norm(double* x, double* y, double*& xnorm, double*& ynorm,double &kx,double &ky)
	{
		double xMin = 100;
		double xMax = 0;
		double yMin = 100;
		double yMax = 0;
		const int nx = _msize(x) / sizeof(*x);
		for (int i = 0; i < nx; i++)
		{
			if (xMin > x[i])
			{
				xMin = x[i];
			}
			if (xMax < x[i])
			{
				xMax = x[i];
			}
			if (yMin > y[i])
			{
				yMin = y[i];
			}
			if (yMax < y[i])
			{
				yMax = y[i];
			}
		}
		/*
		if (x[0]>x[nx-1])
		{
			xMin = x[nx - 1];
			xMax = x[0];
			yMin = y[nx - 1];
			yMax = y[0];
		}
		else
		{
			xMin = x[0];
			xMax = x[nx - 1];
			yMin = y[0];
			yMax = y[nx - 1];
		}*/
		kx = 1 / (xMax - xMin);
		ky = 1/ (yMax - yMin);
		for (int i = 0; i < nx; i++)
		{
			xnorm[i] = (x[i] - xMin) / (xMax - xMin);
			ynorm[i] = (y[i] - yMin) / (yMax - yMin);
		}

	}
	double caucldis(double* fit_x, double* fit_y, double point_x, double point_y)
	{
		const int nx = _msize(fit_x) / sizeof(*fit_x);
		double distance = 0;
		double Mindistance = 100;
		for (int i = 0; i < nx; i++)
		{
			distance = pow(fit_x[i] - point_x, 2) + pow(fit_y[i] - point_y, 2);
			distance = pow(distance, 0.5);
			if (distance < Mindistance)
			{
				Mindistance = distance;
			}
		}
		//delete fit_x;
		//fit_x = NULL;
		//delete fit_y;
		//fit_y = NULL;
		return Mindistance;
	}
	double distance_p(double* fit_x, double* fit_y, double point_x, double point_y, double* als)
	{
		const int nx = _msize(fit_x) / sizeof(*fit_x);
		double distance = 0;
		double Mindistance = 100;
		int MindisIdex = 0;
		int idx1 = 0;
		int idx2 = 0;
		const int PointNums = 1500;

		for (int i = 0; i < nx; i++)
		{
			distance = pow(fit_x[i] - point_x, 2) + pow(fit_y[i] - point_y, 2);
			distance = pow(distance, 0.5);
			if (distance < Mindistance)
			{
				Mindistance = distance;
				MindisIdex = i;
			}
		}
		//cout << Mindistance << endl;
		//return Mindistance;

		if (Mindistance <= inside_dist)
		{
			return Mindistance;
		}
		else
		{
			double* x = new double[PointNums];
			double* y = new double[PointNums];
			for (int i = 0; i < PointNums; i++)
			{
				x[i] = 0;
				y[i] = 0;
			}
			idx2 = MindisIdex + 3;
			idx1 = MindisIdex - 3;
			if (idx2 > nx - 1)
			{
				idx2 = nx - 1;
			}
			if (idx1 < 0)
			{
				idx1 = 0;
			}
			const double interdis = (fit_x[idx2] - fit_x[idx1]) / (PointNums - 1);
			//cout << fit_x[MindisIdex] << '\t' << fit_x[idx1] << '\t' << fit_x[idx2] <<'\t'<< point_x<<endl;
			for (int i = 0; i < PointNums; i++)
			{
				x[i] = fit_x[idx1] + interdis * double(i);
				y[i] = fx(x[i], als);
			}
			Mindistance = caucldis(x, y, point_x, point_y);
			//cout << Mindistance << '\n';
			delete x;
			x = NULL;
			delete y;
			y = NULL;
			return Mindistance;
		}


		/*
		int idx = -1;
		double fitx, fity;
		for (int i = 0; idx ==-1; i++)
		{
			//cout << fit_x[i] << endl;
			if (fit_x[i] >= point_x)
			{
				idx = i;
			}
		}
		if (idx == -1)
		{
			fitx = fit_x[nx-1];
			fity = fit_y[nx-1];
		}
		else
		{
			fitx = fit_x[idx];
			fity = fit_y[idx];
		}
		double distance = pow(fitx - point_x, 2) + pow(fity - point_y, 2);
		distance = pow(distance, 0.5);
		return distance;*/

	}
	double reC(double* distance)
	{
		const int nx = _msize(distance) / sizeof(*distance);
		//for (int i = 0; i < nx; i++)
		//{
			//cout << distance[i] << '\n';
		//}
		//system("pause");
		double cout_dist = 0;
		double insidenums = 0;
		for (int i = 0; i < nx; i++)
		{
			cout_dist = cout_dist + distance[i] * distance[i];
			if (distance[i] <= inside_dist)
			{
				insidenums = insidenums + 1;
			}
		}
		double C = (insidenums / nx) / cout_dist;
		//system("cls");
		return C;
	}
	template <typename Te>
	void ransac(double* R0, double* ne0, int num,double* Randx0,double* Randy0,int isa, double& C,double* &apar,double** Randidxs,Te ainitd)
	{
		//int leng = sizeof(ainitd) / sizeof(ainitd[0]);
		//int leng = _msize(ainitd) / sizeof(*ainitd);
		//std::cout << leng << std::endl;
		//system("pause");
		const int nr0 = _msize(R0) / sizeof(*R0);
		const int nx = 100;

		//double* rand_x = Randx0;
		//double* rand_y = Randy0;
		double* rand_x = new double[num];
		double* rand_y = new double[num];
		random_point(R0, ne0, num, rand_x, rand_y, isa ,Randidxs);
		const int PointNums = nx;
		double* fit_x = new double[nx];
		double* fit_y = new double[nx];
	
		/*
		for (int i = 0; i < 6; i++)
		{
			//rand_x1[i] = 0;
			//rand_y1[i] = 1;
			cout << rand_x[i] << '\t' << rand_y[i] << endl;
		}
		cout << endl;*/
		
		//system("pause");
		tanhFit nft;
		nft.refl_fitting(rand_x, rand_y, fit_x, fit_y, apar, ainitd,nx);
		//cout << 555 << endl;
		//apars = als;
		//apars[isa] = apar;
		double* distance_list = new double[nr0];
		for (int i = 0; i < nr0; i++)
		{
			distance_list[i] = distance_p(fit_x, fit_y, R0[i], ne0[i], apar);
			//cout << distance_list[i] << '\t' << endl;
		}
		C = reC(distance_list);



		delete[] rand_x;
		rand_x = NULL;
		delete[] rand_y;
		rand_y = NULL;
		delete []fit_x;
		fit_x = NULL;
		delete []fit_y;
		fit_y = NULL;
		delete []distance_list;
		distance_list = NULL;
	}
	void random_point(double* R0, double* ne0, int num, double*& data_x0, double*& data_y0, int srandi,double** Randidxs)
	{
		const int nx = _msize(R0) / sizeof(*R0);
		//double* checkrandidx = new double[num];
		int idx = 0;
		//srand(srandi);
		
		for (int i = 0; i < num; i++)
		{
			//idx = rand() % nx;
			//idx = rands() % nx;
			idx = Randidxs[i][srandi];
			data_x0[i] = R0[idx];
			data_y0[i] = ne0[idx];


		}
	}
	void ransac_main(double* mean_R, double* mean_n,int idx1,double* & pedsInfo,int ompthreads)
	{
		
		const int nr = _msize(mean_R) / sizeof(*mean_R);
		//double* pedsInfo = new double[3];
		double* mean_R1 = new double[nr];
		double* mean_n1 = new double[nr];
		for (int i = 0; i < nr; i++)
		{
			mean_R1[i] = mean_R[i];
			mean_n1[i] = mean_n[i];
		}
		double* xnorm = new double[nr];
		double* ynorm = new double[nr];
		double kx = 0;
		double ky = 0;
		Norm(mean_R1, mean_n1,xnorm,ynorm,kx,ky);

		const int num = 7;
		//inside_dist = 0.015;
		const int iter_num = 2000;

		double** Randidxs = new double* [num];
		for (int i = 0; i < num; i++)
		{
			Randidxs[i] = new double[iter_num];

		}
		for (int i = 0; i < num; i++)
		{
			for (int j = 0; j < iter_num; j++)
			{
				Randidxs[i][j] = rands() % nr;
			}
		}
		//double* C = new double[iter_num];
		//double** apars = new double* [iter_num];

		int judge = 0;

		//char* file = "F:/小黑云同步/学术/硕士文件/ransac算法/c++/ransac/SavedataRandixyM.mat";
		//char* pname = "SavedataRandixy";
		//double **Rands=opt.readmat2d(file, pname);
		const int randN = 7;
		double* Randx = new double[randN];
		double* Randy = new double[randN];

		double x_init[6] = { 0.5, 0.5, 0.05 ,0.8,0.09,-0.07 };
		double C[iter_num];
		//double* C = new double[iter_num];
		double** apars = new double* [iter_num];
		for (int i = 0; i < iter_num; i++)
			apars[i] = new double[6];
		//omp_lock_t lock;
		//omp_init_lock(&lock); //初始化互斥锁
		omp_set_num_threads(ompthreads);
		#pragma omp parallel for
		for (int i = 0; i < iter_num; i++)
		{
			int judge = 0;
			double C0 = 0;
			double* apar = new double[6];
			//double apar[6];

			while (judge == 0)
			{
				
				try
				{
					//std::cout << i << std::endl;
					ransac(xnorm, ynorm, num,Randx, Randy,i,C0,apar, Randidxs, x_init);
					judge = 1;
					//cout << C0 << '\t' << Ci << endl;
					//omp_set_lock(&lock); //获得互斥器
					//cout << C0 << '\t' << Ci << endl;
					//omp_unset_lock(&lock); //释放互斥器
					C[i] = C0;
					apars[i][0] = apar[0];
					apars[i][1] = apar[1];
					apars[i][2] = apar[2];
					apars[i][3] = apar[3];
					apars[i][4] = apar[4];
					apars[i][5] = apar[5];
				}
				catch (exception e)
				{
					cout << e.what() << endl;   //捕获异常，然后程序结束
				}
			}
			delete []apar;
			apar = NULL;
		}
		//std::cout << "......" << std::endl;
		//double* aparc = x_init;

		double max_C = 0;
		int idx = 0;
		for (int i = 0; i < iter_num; i++)
		{
			if (C[i] >= max_C)
			{
				max_C = C[i];
				idx = i;
			}
		}
		//double* a1 = apars[idx];
		double a1[6];
		for (int i = 0; i < 6; i++)
		{
			a1[i] = apars[idx][i];
		}
		double* x1 = new double[1000];
		double* y1 = new double[1000];
		const int Nx = _msize(x1) / sizeof(*x1);
		double x00 = 0.5;
		const double interdis = (1 - x00) / (Nx - 1);
		double gradit0 = 0;
		double graditmax = 0;
		for (int i = 0; i < Nx; i++)
		{
			x1[i] = x00 + interdis * double(i);
			y1[i] = fx(x1[i], a1);
			if (i > 0)
			{
				gradit0 = abs((y1[i] - y1[i - 1]) / (x1[i] - x1[i - 1]));
				if (gradit0 > graditmax)
				{
					graditmax = gradit0;
				}
			}
			//cout << graditmax << endl;
			//cout << x[i] << '\t' << y[i] << endl;
		}
		pedsInfo[0] = (a1[0] + a1[1])/ky;//height
		pedsInfo[1] = (2 * a1[4])/kx;
		pedsInfo[2] = kx*graditmax/ky;
		//cout << pedsInfo[0] << '\t' << pedsInfo[1] << '\t' << pedsInfo[2] << endl;


		//delete[] a1;
		//a1 = NULL;
		delete []mean_R1;
		mean_R1 = NULL;
		delete []mean_n1;
		mean_n1 = NULL;
		delete []xnorm;
		xnorm = NULL;
		delete []ynorm;
		ynorm = NULL;
		for (int i = 0; i < num; i++)
		{
			delete[] Randidxs[i];
			Randidxs[i] = NULL;
		}
		//Randidxs = nullptr;
		delete[] Randidxs;
		Randidxs = NULL;


		for (int i = 0; i < iter_num; i++)
		{
			delete[] apars[i];
			apars[i] = NULL;
		}
		//Randidxs = nullptr;
		delete[] apars;
		apars = NULL;


		delete []Randx;
		Randx = NULL;
		delete []Randy;
		Randy = NULL;
		delete []x1;
		x1 = NULL;
		delete []y1;
		y1 = NULL;

		
		//char* s1 = "SavedataResultC";
		//opt.writemat2d(SavedataResult,s1);
		//char* s = "SavedataRandidx";
		//opt.writemat2d(SavedataRandidx,s);
		//system("pause");





	}
	private:
		double ds = 0;
		double inside_dist = 0.015;
		//tanhFit nft;
		//const int PointNums = 700;
		opmat opt;
		random_device rands;
};

