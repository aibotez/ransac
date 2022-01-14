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




class Test
{
public:
	//int CurEcxuteNums = 0;
	int Add(const int x, const int y)
	{

		std::cout << x + y << std::endl;
		return x + y;
	}

	int Del(const int x, const int y)
	{
		return x - y;
	}
	int GetCurProcess()
	{
		int CurEcxuteNums1 = 0;
		return CurEcxuteNums1;
	}

	void Ranact(const int ompthreads,double* R0,double* Z0)
	{

		const int datanums = Col;
		//double** pedsInfo = new double* [3];
		//double** pedsInfo = new double* [3];
		for (int i = 0; i < 3; i++)
		{
			pedsInfo[i] = 0;
		}
		//#pragma omp parallel for
		for (int i = 0; i < datanums; i++)
		{
			
			cut_pro_ransac rans;
			rans.ransac_main(R0, Z0, i, pedsInfo, ompthreads);

		}

	}
	boost::python::list ExcuteRan(boost::python::list& R, boost::python::list& Z,const int ompthreads)
	{
		int row = len(R);
		int col = 1;
		Row = row;
		Col = col;
		double* R0 = new double[row];
		double* Z0 = new double[row];
		for (int i = 0; i < Row; i++)
		{
			R0[i] = 0;
			Z0[i] = 0;
		}

		for (int i = 0; i < row; ++i)
		{
			R0[i] = boost::python::extract<double>(R[i]);
			Z0[i] = boost::python::extract<double>(Z[i]);
		}
		Ranact(ompthreads,R0,Z0);
		boost::python::list peds;
		for (int i = 0; i < 3; ++i)
		{
			peds.append(pedsInfo[i]);
		}
		delete []R0;
		R0 = NULL;
		delete []Z0;
		Z0 = NULL;
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


		return ret1;
	}
private:
	//double* data0;
	//double* R0;
	//double* Z0;
	double* pedsInfo = new double [3];
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
		.def("ExcuteRan", &Test::ExcuteRan)
		.def("GetCurProcess", &Test::GetCurProcess)
		.def("Del", &Test::Del);
}
