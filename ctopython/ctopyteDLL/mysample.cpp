#define BOOST_PYTHON_STATIC_LIB
#define BOOST_LIB_NAME "boost_numpy"
#include "cut_pro_ransac.h"
//#include <boost/config/auto_link.hpp>
// 当使用预编译的头时，需要使用此源文件，编译才能成功。
#include <iostream>
#include <string>
#define BOOST_PYTHON_STATIC_LIB
#include <boost/python.hpp>   // 必须引入这个头文件
//#include <boost/python/numpy.hpp>
#include "boost/multi_array.hpp"
using namespace std;
using namespace boost::python;

static omp_lock_t lock;


int CurEcxuteNums = 0;
class CResult
{
public:

	int GetCurProcess()
	{
		//CurEcxuteNums = CurEcxuteNums + 1;
		cout << "Cur: " << CurEcxuteNums0 << endl;
		return CurEcxuteNums0;
	}
	void AddProcess(int icurs)
	{
		CurEcxuteNums0 = icurs;
		//cout << "Cur: " << CurEcxuteNums << endl;
		//return CurEcxuteNums;
	}
private:
	int CurEcxuteNums0 = 0;
};


class Test
{
public:
	CResult CR;
	//int CurEcxuteNums = 0;
	int Add(const int x, const int y)
	{

		std::cout << x + y << std::endl;
		return x + y;
	}
	void tst()
	{
		data0 = new double* [Row];
		for (int i = 0; i < Row; i++)
		{
			data0[i] = new double[Col];
		}
	}
	int Del(const int x, const int y)
	{
		return x - y;
	}
	int GetCurProcess()
	{
		int CurEcxuteNums1 = 0;
		CurEcxuteNums1 = CR.GetCurProcess();
		cout << "Cur1: "<< CR.GetCurProcess() << endl;
		return CurEcxuteNums1;
	}

	void Ranact()
	{ 
		omp_init_lock(&lock);
		const int datanums = Col;
		//double** pedsInfo = new double* [3];
		//double** pedsInfo = new double* [3];
		for (int i = 0; i < 3; i++)
		{
			pedsInfo[i] = new double[datanums];
		}
		//#pragma omp parallel for
		for (int i = 0; i < datanums; i++)
		{
			omp_set_lock(&lock);
			//cout << i << endl;
			CurEcxuteNums += 1;
			CR.AddProcess(CurEcxuteNums);
			CR.GetCurProcess();
			cout << CurEcxuteNums << endl;
			omp_unset_lock(&lock);
			double* peds = new double[3];
			double* xx1 = new double[Row];
			double* yy1 = new double[Row];
			for (int ro = 0; ro < Row; ro++)
			{
				xx1[ro] = R0[ro][i];
				yy1[ro] = Z0[ro][i];
			}
			cut_pro_ransac rans;
			rans.ransac_main(xx1, yy1, i, peds);
			for (int m = 0; m < 3; m++)
			{
				pedsInfo[m][i] = peds[m];
			}
		}
		omp_destroy_lock(&lock);

	}
	boost::python::list ExcuteRan(boost::python::list& R, boost::python::list& Z)
	{
		boost::python::list peds;
		int row = len(R);
		int col = len(R[0]);
		Row = row;
		Col = col;
		R0 = new double* [Row];
		Z0 = new double* [Row];
		for (int i = 0; i < Row; i++)
		{
			R0[i] = new double[Col];
			Z0[i] = new double[Col];
		}
		
		for (int i = 0; i < row; ++i)
		{
			for (int j = 0; j < col; j++)
			{
				R0[i][j] = boost::python::extract<double>(R[i][j]);
				Z0[i][j] = boost::python::extract<double>(Z[i][j]);
			}
		}
		Ranact();
		
		for (int i = 0; i < 3; ++i)
		{
			boost::python::list pedstemp;
			for (int j = 0; j < Col; j++)
			{
				pedstemp.append(pedsInfo[i][j]);
			}
			peds.append(pedstemp);
		}
		return peds;
	}
	boost::python::list Square(boost::python::list& data, boost::python::list& lst)
	{
		boost::python::list ret;
		boost::python::list ret1;
		std::cout << len(lst) << std::endl;
		std::cout << len(lst[0]) << std::endl;
		int row = len(lst);
		int col = len(lst[0]);
		Row = row;
		Col = col;
		tst();
		for (int i = 0; i < row; ++i)
		{
			for (int j = 0; j < col; j++)
			{
				data0[i][j] = boost::python::extract<double>(lst[i][j]);
				ret1.append(lst[i][j]);
			}
		}
		for (int i = 0; i < len(data); ++i)
		{
			ret.append(data[i] * data[i]);
		}

		return ret1;
	}
private:
	double** data0;
	double** R0;
	double** Z0;
	double** pedsInfo = new double* [3];
	int Row = 0;
	int Col = 0;
	
};

BOOST_PYTHON_MODULE(ctopytest)
{
	using namespace boost::python;
	//Py_Initialize();
	//np::initialize();
	class_<Test>("Test")
		.def("Add", &Test::Add)
		.def("Square", &Test::Square)
		.def("ExcuteRan", &Test::ExcuteRan)
		.def("GetCurProcess", &Test::GetCurProcess)
		.def("Del", &Test::Del);
}
