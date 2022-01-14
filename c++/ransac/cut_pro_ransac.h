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
	void Norm(double* x, double* y, double*& xnorm, double*& ynorm)
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
	double reC(double* distance, double inside_dist)
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
	void ransac(double* R0, double* ne0, int num, double inside_dist,double* Randx0,double* Randy0,int isa, double& C,double* &apars,double* &checkrandidx)
	{
		const int nr0 = _msize(R0) / sizeof(*R0);
		const int nx = 100;

		//double* rand_x = Randx0;
		//double* rand_y = Randy0;
		double* rand_x = new double[6];
		double* rand_y = new double[6];
		random_point(R0, ne0, 6, rand_x, rand_y, isa,checkrandidx);
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
		double* als;
		tanhFit nft;
		nft.refl_fitting(rand_x, rand_y, fit_x, fit_y,als,nx);
		apars = als;
		double* distance_list = new double[nr0];
		for (int i = 0; i < nr0; i++)
		{
			distance_list[i] = distance_p(fit_x, fit_y, R0[i], ne0[i], als);
			//cout << distance_list[i] << '\t' << endl;
		}
		C = reC(distance_list, inside_dist);

		delete rand_x;
		rand_x = NULL;
		delete rand_y;
		rand_y = NULL;
		delete fit_x;
		fit_x = NULL;
		delete fit_y;
		fit_y = NULL;
		delete distance_list;
		distance_list = NULL;
	}
	void random_point(double* R0, double* ne0, int num, double*& data_x0, double*& data_y0, int srandi,double* &checkrandidx)
	{
		const int nx = _msize(R0) / sizeof(*R0);
		//double* checkrandidx = new double[num];
		int idx = 0;
		srand(srandi);
		random_device e;
		for (int i = 0; i < num; i++)
		{
			//idx = rand() % nx;
			idx = e() % nx;
			data_x0[i] = R0[idx];
			data_y0[i] = ne0[idx];
			checkrandidx[i] = idx;

		}
	}
	void ransac_main(double* mean_R, double* mean_n,double* &x1,double* &y1)
	{
		
		const int nr = _msize(mean_R) / sizeof(*mean_R);
		double* xnorm = new double[nr];
		double* ynorm = new double[nr];
		Norm(mean_R, mean_n,xnorm,ynorm);

		const int num = 6;
		//inside_dist = 0.015;
		int iter_num = 2000;
		double* C = new double[iter_num];
		double** apars = new double* [iter_num];
		double** SavedataRandidx = new double* [num];
		for (int i = 0; i < num; i++)
		{
			SavedataRandidx[i] = new double[iter_num];

		}
		int judge = 0;

		char* file = "F:/小黑云同步/学术/硕士文件/ransac算法/c++/ransac/SavedataRandixyM.mat";
		char* pname = "SavedataRandixy";
		double **Rands=opt.readmat2d(file, pname);
		const int randN = 6;
		double* Randx = new double[randN];
		double* Randy = new double[randN];
		#pragma omp parallel for
		for (int i = 0; i < iter_num; i++)
		{
			int judge = 0;
			double C0 = 0;
			double* apar;
			double* checkrandidx = new double[num];
			for (int ro = 0; ro < randN; ro++)
			{
				Randx[ro] = Rands[ro][i];
				Randy[ro] = Rands[ro+ randN][i];
			}

			while (judge == 0)
			{
				
				try
				{
					//std::cout << i << std::endl;
					ransac(xnorm, ynorm, num, inside_dist,Randx, Randy,i,C0,apar, checkrandidx);
					judge = 1;
					C[i] = C0;
					apars[i] = apar;
					for (int ix = 0; ix < num; ix++)
					{
						SavedataRandidx[ix][i] = checkrandidx[ix];
					}
				}
				catch (exception e)
				{
					cout << e.what() << endl;   //捕获异常，然后程序结束
				}
			}
		}

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
		double* a1 = apars[idx];

		const int Nx = 200;
		double x00 = 0.5;
		const double interdis = (1 - x00) / (Nx - 1);
		for (int i = 0; i < Nx; i++)
		{
			x1[i] = x00 + interdis * double(i);
			y1[i] = fx(x1[i], a1);
			//cout << x[i] << '\t' << y[i] << endl;
		}

		//char* s = "SavedataResultC";
		//opt.writemat2d(SavedataResult,s);
		char* s = "SavedataRandidx";
		opt.writemat2d(SavedataRandidx,s);
		//system("pause");





	}
	private:
		double ds = 0;
		double inside_dist = 0.015;
		//tanhFit nft;
		//const int PointNums = 700;
		opmat opt;
};

