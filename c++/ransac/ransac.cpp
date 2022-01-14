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
#include "cut_pro_ransac.h"
#include "opmat.h"
//#include "testl.h"

using namespace std;


void read_data(char* fie, double*& x0, double*& y0)
{
	ifstream infile;
	infile.open(fie);   //将文件流对象与文件连接起来 
	assert(infile.is_open());   //若失败,则输出错误消息,并终止程序运行 
	string s, sz, sne;
	double b = 0;
	double bz = 0;
	double bne = 0;
	int st = 0;
	int col = 1;
	int row = 0;
	vector<vector<double>> data;
	vector<double> vec10;
	while (getline(infile, s))
	{
		if (col == 1)
		{
			for (int idx = 0; idx < s.length(); idx++)
			{
				//cout << s[idx];
				if (s[idx] == '\t')
				{
					col = col + 1;
				}
			}
		}
		col = 2;
		stringstream stream;
		stream << s;
		for (int j = 0; j < col; j++)
		{
			stream >> b;
			//cout << b << endl;
			vec10.push_back(b);
		}
		data.push_back(vec10);
		vec10.clear();
		row = row + 1;
	}
	infile.close();
	int rows = row;
	int cols = col;
	//cout << data[0][1]<<endl;
	double* R1 = new double[rows];
	double* ne1 = new double[rows];
	double* ne2 = new double[rows];
	for (int i = 0; i < row; i++)
	{
		//cout << i<<'\t'<<data[i][0] << endl;
		R1[i] = data[i][0];
		//ne1[i] = data[i][1] * pow(10, 19);
		ne2[i] = data[i][1];
	}
	x0 = R1;
	y0 = ne2;
}


int main()
{
	cut_pro_ransac rans;
	opmat opt;
	//testl rans;
	char fie[] = "F:/小黑云同步/学术/硕士文件/ransac算法/pro_fit_ransac/86787.txt";
	double* xx;
	double* yy;
	double a[6];
	read_data(fie, xx, yy);
	//rans.test();
	DWORD start, end;
	
	//tanhFit nfi;
	double* x = new double[8];
	double* y = new double[8];
	x[0] = 2.26258;
	x[1] = 2.29984;
	x[2] = 1.91715;
	x[3] = 2.31139;
	x[4] = 2.30675;
	x[5] = 2.35;
	x[6] = 2.16339998567287;
	x[7] = 2.25807009854049;

	y[7] = 2.70213356133346;
	y[0] = 2.24956;
	y[1] = 0.813882;
	y[2] = 4.46345;
	y[3] = 0.259411;
	y[4] = 0.367118;
	y[5] = 0;
	x[6] = 3.57429916352006;//10.3056 -5.49436        -57.6974        1.9125  0.130059        0.385792

//	nfi.tannlinfit(xx, yy);
	double* fitx;
	double* fity;
	double nped0;
	double nwidth0;
	double maxgrad0;
	start = GetTickCount();
	/*
#pragma omp parallel for
	for (int i = 0; i < 3000; i++)
	{
		nfi.refl_fitting(xx, yy, fitx, fity, nped0, nwidth0, maxgrad0);
	}*/
	//nfi.refl_fitting(xx, yy, fitx, fity, nped0, nwidth0, maxgrad0);
	//const int N = _msize(fitx) / sizeof(*fitx);
	//for (int i = 0; i < N; i++)
	//{
		//cout << fitx[i] << '\t' << fity[i] << endl;
	//}

	//system("pause");
	const int Nx = 200;
	double* x1 = new double[Nx];
	double* y1 = new double[Nx];
	const int datanums = 10;
	double** data = new double*[Nx];
	for (int i = 0; i < Nx; i++)
	{
		data[i] = new double[datanums + 1];
	}
//#pragma omp parallel for
	for (int i = 0; i < datanums; i++)
	{
		cout << i << endl;
		rans.ransac_main(xx, yy, x1, y1);
		for (int j = 0; j < Nx; j++)
		{
			data[j][0] = x1[j];
			data[j][i+1] = y1[j];
		}
	}
	end = GetTickCount() - start;
	cout << end/1000 << endl;
	char* s = "dataC";
	opt.writemat2d(data, s);
	system("pause");
	return 0;
}

