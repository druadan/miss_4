#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;
#define FUNC_2

double f(double y, double t)
{

	#ifdef FUNC_1
    return pow(y,4) + 2*pow(y,2)*pow(t,2) +  pow(t,4);
	#endif

	#ifdef FUNC_2
    return pow(y,3) + pow(y,2)*t - pow(t,2) - 4*t;
	#endif

	#ifdef FUNC_3
    return sin(y*t);
	#endif

	#ifdef FUNC_4
    return sin(y*t);
	#endif
}

double a_derivative(double y, double t)
{
	#ifdef FUNC_1
    return 4*pow(y,3) + 4*y*pow(t,2);
	#endif

	#ifdef FUNC_2
    return 3*pow(y,2) + 2*y*t;
	#endif

	#ifdef FUNC_3
    return cos(y*t) * t;
	#endif

	#ifdef FUNC_4
    return sin(y*t);
	#endif
}

double b_derivative(double y, double t)
{
	#ifdef FUNC_1
    return 4*pow(y,2)*t + 4*pow(t,3);
	#endif

	#ifdef FUNC_2
    return pow(y,2) - 2*t - 4;
	#endif

	#ifdef FUNC_3
    return cos(y*t) * y;
	#endif

	#ifdef FUNC_4
    return sin(y*t);
	#endif
}

double a_hessian(double y, double t)
{
	#ifdef FUNC_1
    return 12*pow(y,2) + 4*pow(t,2);
	#endif

	#ifdef FUNC_2
    return 6*y + 2*t;
	#endif

	#ifdef FUNC_3
    return -1 * pow(t,2) * sin(y*t);
	#endif

	#ifdef FUNC_4
    return sin(y*t);
	#endif
}

double b_hessian(double y, double t)
{
	#ifdef FUNC_1
    return 8*y*t;
	#endif

	#ifdef FUNC_2
    return 2*y;
	#endif

	#ifdef FUNC_3
    return cos(y*t) - y*t*sin(y*t);
	#endif

	#ifdef FUNC_4
    return sin(y*t);
	#endif
}

double c_hessian(double y, double t)
{
	#ifdef FUNC_1
    return 8*y*t;
	#endif

	#ifdef FUNC_2
    return 2*y;
	#endif

	#ifdef FUNC_3
    return cos(y*t) - y*t*sin(y*t);
	#endif

	#ifdef FUNC_4
    return sin(y*t);
	#endif
}

double d_hessian(double y, double t)
{
	#ifdef FUNC_1
    return 4*pow(y,2) + 12*pow(t,2);
	#endif

	#ifdef FUNC_2
    return -2;
	#endif

	#ifdef FUNC_3
    return -1 * pow(y,2) * sin(y*t);
	#endif

	#ifdef FUNC_4
    return sin(y*t);
	#endif
}


void draw_plot(vector<double> *x, vector<double> *y)
{
	ofstream file;
	file.open("data.txt");

	for(unsigned int i = 0; i < x->size(); ++i)
	{
		file << x->at(i) << " " << y->at(i) << endl;
	}
	file.close();
	system("gnuplot < plot_script");
	
}

int main() {
	const double tmin = 0;
	const double tmax = 1;

	//wartosc poczatkowa
	double x0;
	//dokladnosc
	double e;

    double h = 0.01;
    
	//cin >> x0 >> e;
    
	x0 = 0.1;
	e = 0.01;

	vector<double> t;
	vector<double> x;

	t.push_back(tmin);
	x.push_back(x0);
    
    for(double tv = tmin +h ; tv <= tmax; tv += h)
	{
        double der_a = a_derivative(x.back(),tv);
        double der_b = b_derivative(x.back(),tv);
        
        double a = a_hessian(x.back(),tv);
        double b = b_hessian(x.back(),tv);
        double c = c_hessian(x.back(),tv);
        double d = d_hessian(x.back(),tv);
        
        double inverse_prefix = 1/(a*d - b*c);
        
        double inv_a = inverse_prefix * d;
        double inv_b = inverse_prefix * -1 * c;
        double inv_c = inverse_prefix * -1 * b;
        double inv_d = inverse_prefix * a;
        
		t.push_back(tv);
		x.push_back(x.back() - (inv_c*der_a + inv_d*der_b));
	}

	draw_plot(&t, &x);

	return 0;
}
