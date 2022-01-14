#pragma once
#include <stdlib.h>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <gsl/gsl_linalg.h>
#include <vector>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_rng.h>

#pragma comment(lib, "gsl.lib")
#pragma comment(lib, "gslcblas.lib")


//int n_nums = 20;//对e指数展开，n_nums为最大次幂
struct data {
	size_t n;
	double* t;
	double* y;
};

double factorial(int n)
{
	double rst = 1;
	for (int i = 1; i <= n; i++)
	{
		rst = rst * i;
	}
	return rst;
}
double fx(double x, double a[6])
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
int expb_f(const gsl_vector* x, void* data, gsl_vector* f)
{
	size_t n = ((struct data*)data)->n;
	double* t = ((struct data*)data)->t;
	double* y = ((struct data*)data)->y;

	//double A = gsl_vector_get(x, 0);
	//double lambda = gsl_vector_get(x, 1);
	//double b = gsl_vector_get(x, 0);
	double a1 = gsl_vector_get(x, 0);
	double a2 = gsl_vector_get(x, 1);
	double a3 = gsl_vector_get(x, 2);
	double a4 = gsl_vector_get(x, 3);
	double a5 = gsl_vector_get(x, 4);
	double a6 = gsl_vector_get(x, 5);
	//double a6 = 0;

	size_t i;
	for (i = 0; i < n; i++)
	{
		double m = (a4 - t[i]) / a5;
		double fz = (1 + a3 * m) * exp(m) - (1 + a6 * m) * exp(-m);
		double fm = exp(m) + exp(-m);
		double Yi = a1 * fz / fm + a2;

		double dif = abs(Yi - y[i]);
		/* Model Yi = A * exp(-lambda * t_i) + b */
		//double Yi = A * exp(-lambda * t[i]) + b+a4*t[i];
		//cout << dif << endl;
		gsl_vector_set(f, i, Yi - y[i]);
		//gsl_vector_set(f, i, dif);
	}

	return GSL_SUCCESS;
}

int expb_df(const gsl_vector* x, void* data, gsl_matrix* J)
{
	size_t n = ((struct data*)data)->n;
	double* t = ((struct data*)data)->t;

	//double A = gsl_vector_get(x, 0);
	//double lambda = gsl_vector_get(x, 1);
	double a1 = gsl_vector_get(x, 0);
	double a2 = gsl_vector_get(x, 1);
	double a3 = gsl_vector_get(x, 2);
	double a4 = gsl_vector_get(x, 3);
	double a5 = gsl_vector_get(x, 4);
	double a6 = gsl_vector_get(x, 5);
	//double a6 = 0;

	size_t i;
	//int jn = int(n_nums / 2);
	double m;
	for (i = 0; i < n; i++)
	{
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi = A * exp(-lambda * t_i) + b  */
		/* and the xj are the parameters (A,lambda,b) */
		m = (a4 - t[i]) / a5;
		double fma = exp(m) + exp(-m);
		double dya1 = ((1 + a3 * m) * exp(m) - (1 + a6 * m) * exp(-m)) / fma;
		double dya2 = 1;
		double dya3 = a1 * m * exp(m) / fma;
		double dya6 = -a1 * m * exp(-m) / fma;



		double fz = ((1 + a3 * m) * exp(m) - (1 + a6 * m) * exp(-m)) / 1;
		double fm = fma;
		double dfz = (a3 * m + a3 + 1) * exp(m) + (a6 * m - a6 + 1) * exp(-m);
		double dfm = exp(m) - exp(-m);
		double ds = (dfz / fm) - fz * dfm / fm / fm;
		ds = a1 * ds;

		/*
		double ym = 2 + m * m;
		double yz = (a3 - a6 + 2) * m + (a3 + a6) * m * m;
		double dyz = (a3 - a6 + 2) + (a3 + a6) * 2 * m;
		double dym = 2 * m;
		for (int n = 2; n <= jn; n++)
		{
			double n1 = n;
			yz = yz + ((a3 - a6) / factorial(2 * n - 2) + 2 / factorial(2 * n - 1)) * pow(m, 2 * n - 1) + ((a3 + a6) / factorial(2 * n - 1)) * pow(m, 2 * n);
			dyz = dyz + ((a3 - a6) / factorial(2 * n - 2) + 2 / factorial(2 * n - 1)) * (2 * n1 - 1) * pow(m, 2 * n - 2) + ((a3 + a6) / factorial(2 * n - 1)) * (2 * n1) * pow(m, 2 * n - 1);
			ym = ym + pow(m, 2 * n) / factorial(2 * n - 1);
			dym = dym + 2 * n1 * pow(m, 2 * n - 1) / factorial(2 * n - 1);
		}
		double ds = dyz / ym - yz * dym / pow(ym, 2);
		*/


		double dya4 = ds / a5;
		double dya5 = ds * (t[i] - a4) / pow(a5, 2);
		//cout << dya1 << '\t' << dya2 << '\t' << dya3 << '\t' << dya4 << '\t' << dya4 << '\t' << dya5 << '\t' << dya6 << endl;
		//double e = exp(-lambda * t[i]);
		gsl_matrix_set(J, i, 0, dya1);
		gsl_matrix_set(J, i, 1, dya2);
		gsl_matrix_set(J, i, 2, dya3);
		gsl_matrix_set(J, i, 3, dya4);
		gsl_matrix_set(J, i, 4, dya5);
		gsl_matrix_set(J, i, 5, dya6);
	}

	return GSL_SUCCESS;
}


class tanhFit
{
public:
	void callback(const size_t iter, void* params, const gsl_multifit_nlinear_workspace* w)
	{
		gsl_vector* f = gsl_multifit_nlinear_residual(w);
		gsl_vector* x = gsl_multifit_nlinear_position(w);
		double rcond;

		/* compute reciprocal condition number of J(x) */
		gsl_multifit_nlinear_rcond(&rcond, w);

		fprintf(stderr, "iter %2zu: A = %.4f, lambda = %.4f, b = %.4f, a4 = %.4f, cond(J) = %8.4f, |f(x)| = %.4f\n",
			iter,
			gsl_vector_get(x, 0),
			gsl_vector_get(x, 1),
			gsl_vector_get(x, 2),
			gsl_vector_get(x, 3),
			1.0 / rcond,
			gsl_blas_dnrm2(f));
	}
	void tannlinfit(double* x0, double* y0)
	{
		const int N = _msize(x0) / sizeof(*x0);

		const gsl_multifit_nlinear_type* T = gsl_multifit_nlinear_trust;
		gsl_multifit_nlinear_workspace* w;
		gsl_multifit_nlinear_fdf fdf;
		gsl_multifit_nlinear_trs* rr;

		//gsl_multifit_nlinear_parameters fdf_params;
		gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
		const size_t n = N;
		const size_t p = 6;

		gsl_vector* f;
		gsl_matrix* J;
		gsl_matrix* covar = gsl_matrix_alloc(p, p);
		gsl_matrix_set_zero(covar);
		//double t[N], y[N], weights[N];
		double* t = x0;
		double* y = y0;
		double* weights = new double[N];
		struct data d = { n, t, y };
		double h = 2;
		//double w = 0.8;
		double slope1 = 0.03;
		double ped_pos = 2.25;
		double w1 = 0.8;
		double slope2 = 7;
		double x_init[6] = { h / 2, h / 2, slope1 ,ped_pos,w1,slope2 }; /* starting values */
		//double x_init[5] = { h / 2, h / 2, slope1 ,ped_pos,w1};
		gsl_vector_view x = gsl_vector_view_array(x_init, p);
		//gsl_vector_view wts = gsl_vector_view_array(weights, n);
		gsl_rng* r;
		double chisq, chisq0;
		int status, info;
		size_t i;
		const double xtol = 1e-7;
		const double gtol = 1e-7;
		const double ftol = 0.0;
		gsl_rng_env_setup();
		r = gsl_rng_alloc(gsl_rng_default);
		/* define the function to be minimized */
		//fdf.f=&nlinfit.expb_f;
		fdf.f = expb_f;
		fdf.df = expb_df;   /* set to NULL for finite-difference Jacobian */
		fdf.fvv = NULL;     /* not using geodesic acceleration */
		fdf.n = n;
		fdf.p = p;
		fdf.params = &d;
		/* allocate workspace with default parameters */
		w = gsl_multifit_nlinear_alloc(T, &fdf_params, n, p);
		/* initialize solver with starting point and weights */
		gsl_multifit_nlinear_winit(&x.vector, NULL, &fdf, w);
		//gsl_multifit_nlinear_winit(&x.vector, &wts.vector, &fdf, w);
		/* compute initial cost function */
		f = gsl_multifit_nlinear_residual(w);
		gsl_blas_ddot(f, f, &chisq0);
		/* solve the system with a maximum of 100 iterations */
		status = gsl_multifit_nlinear_driver(20000, xtol, gtol, ftol, NULL
			, NULL, &info, w);
		//status = gsl_multifit_nlinear_driver(10000, xtol, gtol, ftol, NULL
			//, NULL, &info, w);
		/* compute covariance of best fit parameters */
		J = gsl_multifit_nlinear_jac(w);
		gsl_multifit_nlinear_covar(J, 0.0, covar);
		/* compute final cost */
		gsl_blas_ddot(f, f, &chisq);

		for (int i = 0; i < 6; i++)
		{
			//a[i] = gsl_vector_get(w->x, i);
			as[i] = gsl_vector_get(w->x, i);
			//std::cout << gsl_vector_get(w->x, i) << '\t';
		}
		//std::cout << std::endl;


		gsl_multifit_nlinear_free(w);
		gsl_matrix_free(covar);
		gsl_rng_free(r);

	}
	void refl_fitting(double* x0, double* y0, double*& fitx, double*& fity, double& nped0, double& nwidth0,double* &a,int nx)
	{
		//const int N = _msize(x0) / sizeof(*x0);
		//std::cout << N;
		tannlinfit(x0, y0);
		a = as;

		const int PointNums = nx;
		double x00 = 2.12;
		//double* x = new double[PointNums];
		//double* y = new double[PointNums];
		
		const double interdis = (2.35 - x00) / (PointNums-1);
		//double* grad = new double[PointNums - 1];
		double maxgrad = 0;
		for (int i = 0; i < PointNums; i++)
		{
			fitx[i] = x00 + interdis * double(i);
			fity[i] = fx(fitx[i], as);
		}
		nped0 = as[0] + as[1];
		nwidth0 = 2 * as[4];
		//std::cout << interdis<<'\t'<<x[0] << std::endl;
	}
private:
	double *as = new double[6];
	/*
	int n_nums = 20;//对e指数展开，n_nums为最大次幂
	struct data {
		size_t n;
		double* t;
		double* y;
	};*/
};



