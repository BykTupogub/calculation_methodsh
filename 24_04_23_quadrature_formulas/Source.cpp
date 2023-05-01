#include<iostream>
#include<cmath>
#include <iomanip> 
#include<vector>
#include<utility>

using namespace std;


double f1(double x);
double f2(double x);
double f3(double x);
double f4(double x);
double f5(double x);
double if1(double a, double b);
double if2(double a, double b);
double if3(double a, double b);
double if4(double a, double b);
double if5(double a, double b);
vector <double> df1(double a, double b);
vector <double> df2(double a, double b);
vector <double> df3(double a, double b);
vector <double> df4(double a, double b);
vector <double> df5(double a, double b);


double left_rectangle(double(*f)(double), double a, double b);
double right_rectangle(double(*f)(double), double a, double b);
double medium_rectangle(double(*f)(double), double a, double b);
double trapezoid(double(*f)(double), double a, double b);
double simpson(double(*f)(double), double a, double b);
double three_eighths(double(*f)(double), double a, double b);

void exercise1(double(*f)(double), double(*fi)(double, double), double a, double b);
double kf(double(*f)(double), double(*kf)(double(*)(double), double, double), double a, double b, int m);
void exercise2(double(*f)(double), double(*fi)(double, double), vector <double>(*df)(double, double), double a, double b, int m);
void exercise3(double(*f)(double), double(*fi)(double, double), vector <double>(*df)(double, double), double a, double b, int m, int l);

int main()
{
	setlocale(LC_ALL, "Russian");
	cout << setprecision(16);

	cout << "Задание 4. Приближённое вычисление интеграла по квадратурным формулам" << endl;

	int menu = 6;
	double a = 0;
	double b = 0;
	int fnum = 0;
	int m = 10;
	double l = 1;

		cout << "Введите пределы интегрирования" << endl;
		cout << "a: " << endl;
		cin >> a;
		cout << "b: " << endl;
		cin >> b;

		cout << "Выберите функцию" << endl;
		cout << "1. f(x) = cos(x)" << endl << "2. f(x) = 2" << endl << "3. f(x) = 2 + x" << endl << "4. f(x) = 2 + x + x^2" << endl << "5. f(x) = 2 + x + x^2 + x^3" << endl;
		cin >> fnum;
		while (fnum < 1 || fnum > 5)
		{
			cin >> fnum;
			if (fnum < 1 || fnum > 5)
			{
				cout << "Выберите значение от 1 до 5" << endl;
			}
		}
		cout << endl;

		while (menu != 0)
		{
			cout << "Введите:" << endl << "0 для выхода из программы" << endl << "1 для выбора других пределов интегрирования" << endl
				<< "2 для выбора другой функции" << endl << "3 для приближенного вычисления по несоставным КФ"
				<< endl << "4 для приближенного вычисленя по составным КФ" << endl;
			cin >> menu;
			while (fnum < 0 || fnum > 5)
			{
				cin >> fnum;
				if (fnum < 0 || fnum > 5)
				{
					cout << "Выберите значение от 0 до 5" << endl;
				}
			}
			if (menu == 0)
			{
				return EXIT_SUCCESS;
			}
			if (menu == 1)
			{
				cout << "Введите пределы интегрирования" << endl;
				cout << "a: " << endl;
				cin >> a;
				cout << "b: " << endl;
				cin >> b;
				cout << endl;
			}
			if (menu == 2)
			{
				cout << "Выберите функцию" << endl;
				cout << "1. f(x) = cos(x)" << endl << "2. f(x) = 2" << endl << "3. f(x) = 2 + x" << endl << "4. f(x) = 2 + x + x^2" << endl << "5. f(x) = 2 + x + x^2 + x^3" << endl;
				cin >> fnum;
				while (fnum < 1 || fnum > 5)
				{
					cin >> fnum;
					if (fnum < 1 || fnum > 5)
					{
						cout << "Выберите значение от 1 до 5" << endl;
					}
				}
				cout << endl;
			}
			if (menu == 3)
			{
				switch (fnum)
				{
				case 1:
				{
					cout << "1. f(x) = cos(x)" << endl;
					exercise1(*f1, *if1, a, b);
					break;
				}

				case 2:
				{
					cout << "2. f(x) = 2" << endl;
					exercise1(*f2, *if2, a, b);
					break;
				}

				case 3: 
				{
					cout << "2. f(x) = 2 + x" << endl;
					exercise1(*f3, *if3, a, b);
					break;
				}

				case 4:
				{
					cout << "2. f(x) = 2 + x + x^2" << endl;
					exercise1(*f4, *if4, a, b);
					break;
				}

				case 5:
				{
					cout << "2. f(x) = 2 + x + x^2 + x^3" << endl;
					exercise1(*f5, *if5, a, b);
					break;
				}

				default: break;
				}
			}

			if (menu == 4)
			{
				cout << "Введите число проежутков деления" << endl;
				cin >> m;
				while (m < 1)
				{
					cout << "Введите значение >= 1" << endl;
					cin >> m;
				}
				cout << "a = " << a << ", b = " << b << ", m = " << m << ", h =" << (b - a) / m << endl;
				switch (fnum)
				{
				case 1:
				{
					cout << "1. f(x) = cos(x)" << endl;
					exercise2(*f1, *if1, *df1, a, b, m);
					break;
				}

				case 2:
				{
					cout << "2. f(x) = 2" << endl;
					exercise2(*f2, *if2, *df2, a, b, m);
					break;
				}

				case 3:
				{
					cout << "2. f(x) = 2 + x" << endl;
					exercise2(*f3, *if3, *df3, a, b, m);
					break;
				}

				case 4:
				{
					cout << "2. f(x) = 2 + x + x^2" << endl;
					exercise2(*f4, *if4, *df4, a, b, m);
					break;
				}

				case 5:
				{
					cout << "2. f(x) = 2 + x + x^2 + x^3" << endl;
					exercise2(*f5, *if5, *df5, a, b, m);
					break;
				}

				default: break;
				}
				while (menu != 0)
				{
					cout << "Введите 1 для того, чтобы увеличить m в l раз" << endl
						<< "Введите 0 для выхода в гланое меню" << endl << endl;
					cin >> menu;
					if(menu == 1)
					{
						cout << "Введите l:" << endl;
						cin >> l;
						cout << "m = " << m << ", l = " << l << endl;
						cout << "a = " << a << ", b = " << b << ", m = " << m << ", h =" << (b - a) / m / l<< endl;
						switch (fnum)
						{
						case 1:
						{
							cout << "1. f(x) = cos(x)" << endl;
							exercise3(*f1, *if1, *df1, a, b, m, l);
							break;
						}

						case 2:
						{
							cout << "2. f(x) = 2" << endl;
							exercise3(*f2, *if2, *df2, a, b, m, l);
							break;
						}

						case 3:
						{
							cout << "2. f(x) = 2 + x" << endl;
							exercise3(*f3, *if3, *df3, a, b, m, l);
							break;
						}

						case 4:
						{
							cout << "2. f(x) = 2 + x + x^2" << endl;
							exercise3(*f4, *if4, *df4, a, b, m, l);
							break;
						}

						case 5:
						{
							cout << "2. f(x) = 2 + x + x^2 + x^3" << endl;
							exercise3(*f5, *if5, *df5, a, b, m, l);
							break;
						}

						default: break;
						}
					}
				}
				menu = 6;
			}
			
			
			menu = 6;
		}
	
	return EXIT_SUCCESS;
}

double f1(double x)
{
	return cos(x);
	//return sin(x);
}
double f2(double x)
{
	return 2;
}
double f3(double x)
{
	return 2 + x;
}
double f4(double x)
{
	return 2 + x + x*x;
}
double f5(double x)
{
	return 2 + x + x*x + x*x*x;
}
double if1(double a, double b)
{
	return sin(b) - sin(a);
	//return -cos(b) + cos(a);
}
double if2(double a, double b)
{
	return 2*b - 2*a;
}
double if3(double a, double b)
{
	return 2 * b + b * b / 2 - (2*a + a*a/2);
}
double if4(double a, double b)
{
	return 2 * b + b * b / 2 + b * b * b / 3 - (2 * a + a * a / 2 + a * a * a / 3);
}
double if5(double a, double b)
{
	return 2 * b + b * b / 2 + b * b * b / 3 + b * b * b * b / 4 - (2 * a + a * a / 2 + a * a * a / 3 + a * a * a * a / 4);
}

vector <double> df1(double a, double b)
{
	vector <double> result;
	result.push_back(0);
	for (int i = 0; i < 1000; ++i)
	{
		if (result[0] < abs(sin(a + i * (b - a) / 1000))) 
		{
			result[0] = abs(sin(a + i * (b - a) / 1000));
		}
	}
	result.push_back(0);
	for (int i = 0; i < 1000; ++i)
	{
		if (result[1] < abs(cos(a + i * (b - a) / 1000)))
		{
			result[1] = abs(cos(a + i * (b - a) / 1000));
		}
	}
	result.push_back(result[0]);
	result.push_back(result[1]);
	return result;
}
vector <double> df2(double a, double b)
{
	vector <double> result;
	result.push_back(0);
	result.push_back(0);
	result.push_back(0);
	result.push_back(0);
	return result;
}
vector <double> df3(double a, double b)
{
	vector <double> result;
	result.push_back(1);
	result.push_back(0);
	result.push_back(0);
	result.push_back(0);
	return result;
}
vector <double> df4(double a, double b)
{
	vector <double> result;
	result.push_back(abs(1 + 2 * b));
	if (result[0] < abs(1 + 2 * a))
	{
		result[0] = abs(1 + 2 * a);
	}
	result.push_back(2);
	result.push_back(0);
	result.push_back(0);

	return result;
}
vector <double> df5(double a, double b)
{
	vector <double> result;
	result.push_back(1 + 2 * b + 3 * b * b);
	if (result[0] < (1 + 2 * a + 3 * a * a))
	{
		result[0] = 1 + 2 * b + 3 * b * b;
	}
	result.push_back(abs(2 + 6 * b));
	if (result[0] < abs(2 + 6 * a))
	{
		result[0] = abs(2 + 6 * a);
	}
	result.push_back(6);
	result.push_back(0);
	return result;
}

double left_rectangle(double(*f)(double) , double a, double b)
{
	return (b-a)*f(a);
}

double right_rectangle(double(*f)(double), double a, double b)
{
	return (b - a) * f(b);
}

double medium_rectangle(double(*f)(double), double a, double b)
{
	return (b - a) * f((a + b) / 2);
}

double trapezoid(double(*f)(double), double a, double b)
{
	return ((b - a) / 2)*(f(a) + f(b));
}

double simpson(double(*f)(double), double a, double b)
{
	return ((b - a) / 6) * (f(a) + 4 * f((a + b) / 2) + f(b));
}

double three_eighths(double(*f)(double), double a, double b)
{
	double h = (b - a) / 3;
	return (b - a) * (f(a) + 3 * f(a + h) + 3 * f(a + 2 * h) + f(b)) / 8;
}

void exercise1(double(*f)(double), double(*fi)(double, double), double a, double b)
{
	cout << "Приближенные значения" << endl;
	cout << "КФ левого прямоугольника: " << left_rectangle(*f, a, b) << endl;
	cout << "КФ правого прямоугольника: " << right_rectangle(*f, a, b) << endl;
	cout << "КФ среднего прямоугольника: " << medium_rectangle(*f, a, b) << endl;
	cout << "КФ трапеции: " << trapezoid(*f, a, b) << endl;
	cout << "КФ Симпсона (или парабол): " << simpson(*f, a, b) << endl;
	cout << "КФ 3/8: " << three_eighths(*f, a, b) << endl;
	cout << endl;
	cout << "Абсолютные фактические погрешности" << endl;
	cout << "КФ левого прямоугольника: " << abs(left_rectangle(*f, a, b) - fi(a, b)) << endl;
	cout << "КФ правого прямоугольника: " << abs(right_rectangle(*f, a, b) - fi(a, b)) << endl;
	cout << "КФ среднего прямоугольника: " << abs(medium_rectangle(*f, a, b) - fi(a, b)) << endl;
	cout << "КФ трапеции: " << abs(trapezoid(*f, a, b) - fi(a, b)) << endl;
	cout << "КФ Симпсона (или парабол): " << abs(simpson(*f, a, b) - fi(a, b)) << endl;
	cout << "КФ 3/8: " << abs(three_eighths(*f, a, b) - fi(a, b)) << endl;
	cout << endl;
}

double kf(double(*f)(double), double(*kf)(double(*)(double), double, double), double a, double b, int m)
{
	double result = 0;
	double h = (b - a) / m;
	for (int i = 0; i < m; ++i)
	{
		result += kf(*f, a + h * i, a + h * (i + 1));
		//cout << result << endl;
	}
	return result;
}

void exercise2(double(*f)(double), double(*fi)(double, double), vector <double>(*df)(double, double), double a, double b, int m)
{
	double h = (b - a) / m;
	cout << "Приближенные значения" << endl;
	double J1 = kf(*f, *left_rectangle, a, b, m);
	cout << "КФ левого прямоугольника: " << J1 << endl;
	double J2 = kf(*f, *right_rectangle, a, b, m);
	cout << "КФ правого прямоугольника: " << J2 << endl;
	double J3 = kf(*f, *medium_rectangle, a, b, m);
	cout << "КФ среднего прямоугольника: " << J3 << endl;
	double J4 = kf(*f, *trapezoid, a, b, m);
	cout << "КФ трапеции: " << J4 << endl;
	double J5 = kf(*f, *simpson, a, b, m);
	cout << "КФ Симпсона (или парабол): " << J5 << endl;
	cout << endl;
	cout << "Абсолютные фактические погрешности" << endl;
	cout << "КФ левого прямоугольника: " << abs(J1 - fi(a, b)) << endl;
	cout << "КФ правого прямоугольника: " << abs(J2 - fi(a, b)) << endl;
	cout << "КФ среднего прямоугольника: " << abs(J3 - fi(a, b)) << endl;
	cout << "КФ трапеции: " << abs(J4 - fi(a, b)) << endl;
	cout << "КФ Симпсона (или парабол): " << abs(J5 - fi(a, b)) << endl;
	cout << endl;
	cout << "Теоретическая погрешность" << endl;
	cout << "КФ левого прямоугольника: " << (b-a)*(df(a,b)[0])*h/2 << endl;
	cout << "КФ правого прямоугольника: " << (b - a) * (df(a, b)[0]) * h / 2 << endl;
	cout << "КФ среднего прямоугольника: " << (b - a) * (df(a, b)[1]) * h * h / 24 << endl;
	cout << "КФ трапеции: " << (b - a) * (df(a, b)[1]) * h * h / 12 << endl;
	cout << "КФ Симпсона (или парабол): " << (b - a) * (df(a, b)[3]) * h * h * h * h / 12 << endl;
	cout << endl;
}

void exercise3(double(*f)(double), double(*fi)(double, double), vector<double>(*df)(double, double), double a, double b, int m, int l)
{
	double h = (b - a) / m / l;
	cout << "Приближенные значения" << endl;
	double J1 = kf(*f, *left_rectangle, a, b, m*l);
	cout << "КФ левого прямоугольника: " << J1 << endl;
	double J2 = kf(*f, *right_rectangle, a, b, m*l);
	cout << "КФ правого прямоугольника: " << J2 << endl;
	double J3 = kf(*f, *medium_rectangle, a, b, m*l);
	cout << "КФ среднего прямоугольника: " << J3 << endl;
	double J4 = kf(*f, *trapezoid, a, b, m*l);
	cout << "КФ трапеции: " << J4 << endl;
	double J5 = kf(*f, *simpson, a, b, m*l);
	cout << "КФ Симпсона (или парабол): " << J5 << endl;
	cout << endl;
	cout << "Абсолютные фактические погрешности" << endl;
	cout << "КФ левого прямоугольника: " << abs(J1 - fi(a, b)) << endl;
	cout << "КФ правого прямоугольника: " << abs(J2 - fi(a, b)) << endl;
	cout << "КФ среднего прямоугольника: " << abs(J3 - fi(a, b)) << endl;
	cout << "КФ трапеции: " << abs(J4 - fi(a, b)) << endl;
	cout << "КФ Симпсона (или парабол): " << abs(J5 - fi(a, b)) << endl;
	cout << endl;
	cout << "Теоретическая погрешность" << endl;
	cout << "КФ левого прямоугольника: " << (b - a) * (df(a, b)[0]) * h / 2 << endl;
	cout << "КФ правого прямоугольника: " << (b - a) * (df(a, b)[0]) * h / 2 << endl;
	cout << "КФ среднего прямоугольника: " << (b - a) * (df(a, b)[1]) * h * h / 24 << endl;
	cout << "КФ трапеции: " << (b - a) * (df(a, b)[1]) * h * h / 12 << endl;
	cout << "КФ Симпсона (или парабол): " << (b - a) * (df(a, b)[3]) * h * h * h * h / 12 << endl;
	cout << endl;

	double JJ1 = kf(*f, *left_rectangle, a, b, m);
	double JJ2 = kf(*f, *right_rectangle, a, b, m);
	double JJ3 = kf(*f, *medium_rectangle, a, b, m);
	double JJ4 = kf(*f, *trapezoid, a, b, m);
	double JJ5 = kf(*f, *simpson, a, b, m);
	double U1 = (l * J1 - JJ1) / (l - 1);
	double U2 = (l * J2 - JJ2) / (l - 1);
	double U3 = (l * l * J3 - JJ3) / (l * l - 1);
	double U4 = (l * l * J4 - JJ4) / (l * l - 1);
	double U5 = (pow(l, 4) * J5 - JJ5) / (pow(l, 4) - 1);

	cout << "Уточнённые значения по принципу Рунге" << endl;
	cout << "КФ левого прямоугольника: " << U1 << endl;
	cout << "КФ правого прямоугольника: " << U2 << endl;
	cout << "КФ среднего прямоугольника: " << U3 << endl;
	cout << "КФ трапеции: " << U4 << endl;
	cout << "КФ Симпсона (или парабол): " << U5 << endl;
	cout << endl;
	cout << "Абсолютные фактические погрешности уточнённых значений" << endl;
	cout << "КФ левого прямоугольника: " << abs(U1 - fi(a, b)) << endl;
	cout << "КФ правого прямоугольника: " << abs(U2 - fi(a, b)) << endl;
	cout << "КФ среднего прямоугольника: " << abs(U3 - fi(a, b)) << endl;
	cout << "КФ трапеции: " << abs(U4 - fi(a, b)) << endl;
	cout << "КФ Симпсона (или парабол): " << abs(U5 - fi(a, b)) << endl;
	cout << endl;
}
