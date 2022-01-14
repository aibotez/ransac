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

	void ransac(double* R0, double* ne0, int num, double inside_dist,double* Randx0,double* Randy0)
	{
		//const int nx = _msize(R0) / sizeof(*R0);
		const int nx = 200;

		//double* rand_x = Randx0;
		//double* rand_y = Randy0;
		const int PointNums = nx;
		double* fit_x = new double[nx];
		double* fit_y = new double[nx];
		double nped = 0;
		double nwidth = 0;
		double maxgrad = 0;
	
				
		for (int i = 0; i < 6; i++)
		{
			//rand_x1[i] = 0;
			//rand_y1[i] = 1;
			cout << Randx0[i] << '\t' << Randx0[i] << endl;
		}
		cout << endl;
		//system("pause");
		double* als;
		tanhFit nft;
		nft.refl_fitting(Randx0, Randy0, fit_x, fit_y, nped, nwidth,als,nx);
	}
	void random_point(double* R0, double* ne0, int num, double*& data_x0, double*& data_y0, int srandi)
	{
		const int nx = _msize(R0) / sizeof(*R0);
		//double* checkrandidx = new double[num];
		int idx = 0;
		srand(srandi);
		//random_device e;
		for (int i = 0; i < num; i++)
		{
			idx = rand() % nx;
			data_x0[i] = R0[i];
			data_y0[i] = ne0[i];

		}
	}
	void ransac_main(double* mean_R, double* mean_n,double* &x1,double* &y1)
	{
		const int nr = _msize(mean_R) / sizeof(*mean_R);
		const int num = 6;
		//inside_dist = 0.015;
		int iter_num = 2000;
		//int judge = 0;

		char* file = "F:/小黑云同步/学术/硕士文件/ransac算法/c++/ransac/SavedataRandixyM.mat";
		char* pname = "SavedataRandixy";
		double **Rands=opt.readmat2d(file, pname);
		const int randN = 6;
		double* Randx = new double[randN];
		double* Randy = new double[randN];
		//#pragma omp parallel for
		for (int i = 0; i < iter_num; i++)
		{
			int judge = 0;
			for (int ro = 0; ro < randN; ro++)
			{
				Randx[ro] = Rands[ro][i];
				Randy[ro] = Rands[ro+ randN][i];
			}
				//double* checkrandidx = new double[num];
			while (judge == 0)
			{
				
				try
				{
					//std::cout << i << std::endl;
					ransac(mean_R, mean_n, num, inside_dist,Randx, Randy);
					judge = 1;
				}
				catch (exception e)
				{
					cout << e.what() << endl;   //捕获异常，然后程序结束
				}
			}
		}
	}
	private:
		double ds = 0;
		double inside_dist = 0.015;
		//tanhFit nft;
		//const int PointNums = 700;
		opmat opt;
};

