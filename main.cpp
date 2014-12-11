#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;
typedef double (*fp)(double, double);


double t_start = 0;
double t_end = 1;


double f(double y, double t)
{
    return pow(y,4) + 2*pow(y,2)*pow(t,2) +  pow(t,4);
	//return y + t;
}



void plot(vector<double> *x, vector<double> *y)
{
	ofstream file;
	file.open("data.txt");

	for(unsigned int i = 0; i < x->size(); i++)
	{
		file << x->at(i) << " " << y->at(i) << endl;
	}
	file.close();
	system("gnuplot < plot_script");
	
}

double improved_euler(fp f, double y, double t, double h)
{

	double k1 = f(y, t);
	double k2 = f(y + h*k1, t + h );
	
	double result = y + h*(k1 + k2) / 2;

	return result;

}

double rk3(fp f, double y, double t, double h)
{
	/*
	double k1 = h * f(y, t);
	double k2 = h * f(y + k1 / 2, t + h / 2);
	double k3 = h * f(y - k1 + 2*k2 , t + h);
	double result = y + (k1 + 4*k2 + k3 ) / 6;
	*/
	
	double k1 =  f(y, t);
	double k2 =  f(y + h*k1, t + h );
	double k3 =  f(y + h*(k1+k2)/4 , t + h/2);
	double result = y + h*(k1 + k2 + 4*k3 ) / 6;
	
	return result;

}

double get_h(double h, double e, double euler, double rk){
	double chicken_factor = 0.9;
	double new_h = chicken_factor * h * pow(e/(abs(euler - rk)),1.0/3);
	return new_h;
}

int main() {

	double x0;
	// epsilon
	double e;
	// inital step
    double h = 0.1;
    
	cin >> x0 ;
    cin >> e;

	//x0 = 0.1;
	//e = 0.001;

	vector<double> t_values;
	vector<double> x_values;

	t_values.push_back(t_start);
	x_values.push_back(x0);
    
    for(double t = t_start;  e <= abs( t- t_end); )
	{
		double rk = rk3(f, x_values.back(), t, h);
		double euler = improved_euler(f, x_values.back(), t, h);

		// if error is acceptable
		if(abs(euler - rk) < e){
			t += h;
			t_values.push_back(t);
			x_values.push_back(rk);	
		}
	
		double new_h = get_h(h,e,euler,rk);

		// check if we got outside t_end range, if yes then reduce the step
		if(t_values.back() + new_h > t_end){
			h = t_end - t_values.back();
		}else{
			h = new_h;
		}

	}

	plot(&t_values, &x_values);

	return 0;
}
