#include<iostream>
#include<string>
#include<cmath>
#include<utility>
#include<vector>
#include<ctime>
#include<random>

using namespace std;

double f(double x)
{
	return pow(2, x) - 2 * cos(x);
}

double f_d(double x)
{
	return pow(2, x) * log(2) + 2 * sin(x);
}

double f_dd(double x)
{
	return pow(2, x) * log(2) * log(2) + 2 * cos(x);
}

void root_separation(double A, double B, int N, vector<pair<double, double>>& segments)
{
	double h = (B - A) / N;
	int counter = 0;
	double X1 = A;
	double X2 = X1 + h;
	double Y1 = f(X1);
	double Y2 = 0;
	pair<double, double> segment;
	while (X2 <= B)
	{
		Y2 = f(X2);
		if (Y1 * Y2 <= 0)
		{
			segment.first = X1;
			segment.second = X2;
			segments.push_back(segment);
		}
		X1 = X2;
		X2 = X1 + h;
		Y1 = Y2;
	}
}


void print(vector<pair<double, double>> segments)
{
	for (int i = 0; i < segments.size(); ++i)
	{
		cout << "[" << segments[i].first << ", " << segments[i].second << "] ";
	}
	cout << endl;
}


template <typename T>
void print(vector<T> X)
{
	for (int i = 0; i < X.size(); ++i)
	{
		cout << X[i] << " ";
	}
	cout << endl;
}

void bisection(double epsilon, vector<pair<double, double>> segments, vector<double>& X, vector<double>& delta, vector<int>& M)
{
	double x = 0;
	double d = 0;
	double a = 0;
	double b = 0;
	double c = 0;
	int m = 0;
	for (int i = 0; i < segments.size(); ++i)
	{
		m = 0;
		a = segments[i].first;
		b = segments[i].second;
		while (b - a > 2 * epsilon)
		{
			c = (a + b) / 2;
			if (f(a) * f(c) <= 0)
			{
				b = c;
			}
			else
			{
				a = c;
			}
			++m;
		}
		x = (a + b) / 2;
		d = (b - a) / 2;
		X.push_back(x);
		delta.push_back(d);
		M.push_back(m);
	}
}

void Newtons_method(double epsilon, vector<pair<double, double>> segments, vector<double>& X, vector<double>& delta, vector<int>& M)
{
	double x = 0;
	double d = 0;
	int m = 0;
	int k = 0;
	for (int i = 0; i < segments.size(); ++i)
	{
		m = 0;
		x = (segments[i].first + segments[i].second) / 2;
		if (f(x) * f_dd(x) > 0)
		{
			while (abs(f(x) / f_d(x)) >= epsilon)
			{
				x = x - f(x) / f_d(x);
				++m;
			}

			d = abs(f(x) / f_d(x));
			x = x - f(x) / f_d(x);
			X.push_back(x);
			delta.push_back(d);
			M.push_back(m);
		}
		else
		{
			k = 0;
			while (f(x) * f_dd(x) <= 0 && k <= 100)
			{
				x = (segments[i].first + (rand() % 100) / 100 * (segments[i].second - segments[i].first + 1));
				++k;
			}
			if (k == 100)
			{
				cout << "не выполненеы условия теоремы о сходимости для отрезка №" << i << "!" << endl;
			}
			else
			{
				while (abs(f(x) / f_d(x)) >= epsilon)
				{
					x = x - f(x) / f_d(x);
					++m;
				}

				d = abs(f(x) / f_d(x));
				x = x - f(x) / f_d(x);
				X.push_back(x);
				delta.push_back(d);
				M.push_back(m);
			}
		}
	}
}



void Newtons_method_modified(double epsilon, vector<pair<double, double>> segments, vector<double>& X, vector<double>& delta, vector<int>& M)
{
	double x = 0;
	double x0 = 0;
	double d = 0;
	int m = 0;
	int k = 0;

	for (int i = 0; i < segments.size(); ++i)
	{
		m = 1;
		x0 = (segments[i].first + segments[i].second) / 2;
		if (f(x0) * f_dd(x0) > 0)
		{
			x = x0 - f(x0) / f_d(0);
			while (abs(f(x) / f_d(x0)) >= epsilon)
			{
				x = x - f(x) / f_d(x0);
				++m;
			}

			d = abs(f(x) / f_d(x));
			x = x - f(x) / f_d(x);
			X.push_back(x);
			delta.push_back(d);
			M.push_back(m);
		}
		else
		{
			k = 0;
			while (f(x0) * f_dd(x0) <= 0 && k <= 100)
			{
				x0 = (segments[i].first + (rand() % 100) / 100 * (segments[i].second - segments[i].first + 1));
				++k;
			}
			if (k == 100)
			{
				cout << "не выполненеы условия теоремы о сходимости для отрезка №" << i << "!" << endl;
			}
			else
			{
				x = x0 - f(x0) / f_d(0);
				while (abs(f(x) / f_d(x0)) >= epsilon)
				{
					x = x - f(x) / f_d(x0);
					++m;
				}

				d = abs(f(x) / f_d(x));
				x = x - f(x) / f_d(x);
				X.push_back(x);
				delta.push_back(d);
				M.push_back(m);
			}
		}

	}
}

void secant_method(double epsilon, vector<pair<double, double>> segments, vector<double>& X, vector<double>& delta, vector<int>& M)
{
	double x1 = 0;
	double x2 = 0;
	double x = 0;
	double d = 0;
	int m = 0;
	for (int i = 0; i < segments.size(); ++i)
	{
		m = 1;
		x1 = segments[i].first;
		x2 = segments[i].second;
		while (abs(f(x2) - f(x1)) >= epsilon)
		{
			x = x2;
			x2 = x2 - (f(x2) / (f(x2) - f(x1)) * (x2 - x1));
			x1 = x;
			++m;
		}

		d = abs(f(x) / f_d(x));
		x = x - f(x) / f_d(x);
		X.push_back(x);
		delta.push_back(d);
		M.push_back(m);
	}
}

void print_npk(vector<pair<double, double>> segments)
{
	for (int i = 0; i < segments.size(); ++i)
	{
		cout << (segments[i].first + segments[i].second) / 2 << " ";
	}
	cout << endl;
}

void print_discrepancy(vector<double> X)
{
	for (int i = 0; i < X.size(); ++i)
	{
		cout << abs(f(X[i])) << " ";
	}
	cout << endl;
}

int main()
{
	setlocale(LC_ALL, "Russian");
	srand(time(NULL));
	cout << "ЧИСЛЕННЫЕ МЕТОДЫ РЕШЕНИЯ НЕЛИНЕЙНЫХ УРАВНЕНИЙ" << endl;
	double A = 0;
	double B = 0;
	double epsilon = 0;
	cout << "Введите А:" << endl;
	cin >> A;
	cout << "Введите B:" << endl;
	cin >> B;
	cout << "Введите epsilon:" << endl;
	cin >> epsilon;

	string function = "2^x-2cos(x)";

	cout << "[" << A << ", " << B << "], epsilon = " << epsilon << ", f(x) = " << function << endl;

	int N = 0;
	cout << "Введите N:" << endl;
	cin >> N;

	vector<pair<double, double>> segments;

	root_separation(A, B, N, segments);
	print(segments);
	cout << "counts = " << segments.size() << endl;
	vector<double> X;
	vector<double> delta;
	vector<int> m;
	bisection(epsilon, segments, X, delta, m);
	cout << "-------------------------------------" << endl;
	cout << "Метод бисекции" << endl;
	cout << "Начальные приближения к корню: ";
	print_npk(segments);
	cout << "Количество шагов: ";
	print(m);
	cout << "Приближенное решение с точность " << epsilon << ": ";
	print(X);
	cout << "длины отрезков: ";
	print(delta);
	cout << "абсолютная величина невязки: ";
	print_discrepancy(X);

	X.clear();
	delta.clear();
	m.clear();

	cout << "-------------------------------------" << endl;
	cout << "Метод Ньютона" << endl;
	Newtons_method(epsilon, segments, X, delta, m);
	cout << "Начальные приближения к корню: ";
	print_npk(segments);
	cout << "Количество шагов: ";
	print(m);
	cout << "Приближенное решение с точность " << epsilon << ": ";
	print(X);
	cout << "длины отрезков: ";
	print(delta);
	cout << "абсолютная величина невязки: ";
	print_discrepancy(X);


	X.clear();
	delta.clear();
	m.clear();

	cout << "-------------------------------------" << endl;
	cout << "Модифицированный метод Ньютона" << endl;
	Newtons_method_modified(epsilon, segments, X, delta, m);
	cout << "Начальные приближения к корню: ";
	print_npk(segments);
	cout << "Количсетво шагов: ";
	print(m);
	cout << "Приближенное решение с точность " << epsilon << ": ";
	print(X);
	cout << "длины отрезков: ";
	print(delta);
	cout << "абсолютная величина невязки: ";
	print_discrepancy(X);

	X.clear();
	delta.clear();
	m.clear();

	secant_method(epsilon, segments, X, delta, m);
	cout << "-------------------------------------" << endl;
	cout << "Метод секущих" << endl;
	cout << "Начальные приближения к корню: ";
	print_npk(segments);
	cout << "Количество шагов: ";
	print(m);
	cout << "Приближенное решение с точность " << epsilon << ": ";
	print(X);
	cout << "длины отрезков: ";
	print(delta);
	cout << "абсолютная величина невязки: ";
	print_discrepancy(X);

	return EXIT_SUCCESS;
}


