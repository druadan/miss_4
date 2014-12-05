#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

typedef double (*fpointer)(double, double);

/*
 * Funkcja reprezentujaca funkcje wymuszajaca
 */
double f(double y, double t)
{
	//return sin(t) - y * y - 2 + t * t * t;
	return -2 * t * t * t + 12 * t * t - 20 * t + 8.5 + y;
}

double runge_kutta(fpointer f, double y, double t, double T)
{
	double F1 = T * f(y, t);
	double F2 = T * f(y + F1 / 2, t + T / 2);
	double F3 = T * f(y + F2 / 2, t + T / 2);
	double F4 = T * f(y + F3, t + T);
	double result = y + (F1 + 2 * F2 + 2 * F3 + F4) / 6;
	return result;
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
	system("gnome-open wykres.jpg");
}

int main() {
	const double tmin = 0;
	const double tmax = 1;

	//wartosc poczatkowa
	double x0;
	//dokladnosc
	double e;
	//aktualna dlugosc kroku
	double h;
	//rzad metody
	double p = 4.0;

	cin >> x0 >> e;

	h = 0.01;

	vector<double> t;
	vector<double> x;

	t.push_back(tmin);
	x.push_back(x0);

	/*h = e;
	for(double tv = tmin + h; tv <= tmax; tv += h)
	{
		t.push_back(tv);
		x.push_back(runge_kutta(f, x.back(), tv, h));
	}*/

	/*for(double tv = tmin + h; tv <= (tmax + (h / 2)); tv += h)
	{
		double y1 = x.back();
		double y2 = x.back();

		y1 = runge_kutta(f, y1, tv, 2 * h);
		y2 = runge_kutta(f, y2, tv, h);
		y2 = runge_kutta(f, y2, tv + h, h);

		double delta = fabs(y2 - y1);

		if(delta <= e)
		{
			t.push_back(tv);
			x.push_back(y1);
			h *= 2.0;
		}
		else
		{
			tv -= h;
			h /= 2.0;
		}
	}*/

	for(double tv = tmin; tv < tmax;)
	{
		if(tv + 2 * h > tmax)
		{
			h = (tmax - tv) / 2;
		}

		double xk1 = runge_kutta(f, x.back(), tv, h);
		double xk2 = runge_kutta(f, xk1, tv + h, h);
		double wk2 = runge_kutta(f, x.back(), tv, 2 * h);

		double gamma = ((e * h) / (tmax - tmin)) * ((pow(2.0,p) - 1) / fabs(xk2 - wk2));
		gamma = pow(gamma, 1.0 / p);

		double h1 = 0.8 * gamma * h;

		if(h1 > (5 * h))
		{
			h1 = 5 * h;
		}
		else if(h1 < (h / 5))
		{
			h1 = h / 5;
		}

		if(h1 < h)
		{
			if((fabs(xk2 - wk2) / (pow(2.0, p) - 1)) <= (e * h / (tmax - tmin)))
			{
				t.push_back(tv + h);
				x.push_back(xk1);

				t.push_back(tv + 2 * h);
				x.push_back(xk2);

				tv += (2 * h);
			}
			else
			{
				h = h1;
			}
		}
		else
		{
			t.push_back(tv + h);
			x.push_back(xk1);

			t.push_back(tv + 2 * h);
			x.push_back(xk2);

			tv += (2 * h);

			h = h1;
		}

	}

	draw_plot(&t, &x);

	return EXIT_SUCCESS;
}
