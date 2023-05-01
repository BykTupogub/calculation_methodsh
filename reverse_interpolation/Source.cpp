#include<iostream>
#include<cmath>
#include <iomanip> 
#include<vector>
#include<utility>

using namespace std;


double f(double x);

void table_generation(vector<pair<double, double>>& nodes, double a, double b, int m);

void print_table(vector <pair<double, double>> nodes);

void sort_table(vector <pair<double, double>>& nodes, double x);

void swap_table(vector <pair<double, double>>& nodes);

double Lagrange(vector <pair<double, double>> nodes, double x, int n);

vector<double> Newton_coefficient(vector <pair<double, double>> nodes, int n);

double Newton(vector <pair<double, double>> nodes, vector<double> coefficient, double x, int n);

double dist(double a, double b);

vector<double> areas_of_monotonicity(vector <pair<double, double>> nodes, double a, double b);

vector<vector<pair<double, double>>> cut(vector <pair<double, double>> nodes, vector<double> areas_of_monotonicity);

void print_sliced_table(vector<vector<pair<double, double>>> sliced_​​table);

void print_segments(vector<pair<double, double>> segments);

void print_npk(vector<pair<double, double>> segments);

template <typename T>
void print(vector<T> X);

void root_separation(double a, double b, vector<pair<double, double>>& segments, vector<double> areas_of_monotonicity);

void bisection(double epsilon, vector<pair<double, double>> segments, vector<double>& X, vector<double>& delta, vector<int>& M, vector<vector<pair<double, double>>> sliced_table, double F, int n);

void print_discrepancy(vector<double> X, double F);



int main()
{
	setlocale(LC_ALL, "Russian");
	cout << setprecision(16);

	cout << "ЗАДАЧА ОБРАТНОГО ИНТЕРПОЛИРОВАНИЯ" << endl;
	cout << "ВАРИАНТ 3" << endl;
	cout << "f(x) = exp(x) - x" << endl;
	int m = 0;
	cout << "Введите число значений в таблице:" << endl;
	cin >> m;
	double a = 0;
	double b = 0;
	cout << "Введите границы интервала:" << endl;
	cin >> a >> b;
	vector<pair<double, double>> nodes;
	table_generation(nodes, a, b, m);
	print_table(nodes);
	int menu = 3;
	double F = 0;
	int n = m + 1;
	double L = 0;
	double N = 0;
	
	vector<double> areas = areas_of_monotonicity(nodes, a, b);
	vector<vector<pair<double, double>>> sliced_​​table = cut(nodes, areas);
	print_sliced_table(sliced_​​table);

	double epsilon = 0;
	vector<pair<double, double>> segments;
	vector<double> X;
	vector<double> delta;
	vector<int> M;


	while (menu != 0)
	{
		if (menu == 1)
		{
			cout << "Введите значени функции:" << endl;
			cin >> F;
			for (int i = 0; i < sliced_​​table.size(); ++i)
			{
				if (sliced_​​table[i][0].second < sliced_​​table[i][sliced_​​table[i].size() - 1].second && (F < sliced_​​table[i][0].second || F > sliced_​​table[i][sliced_​​table[i].size() - 1].second)
					|| sliced_​​table[i][0].second > sliced_​​table[i][sliced_​​table[i].size() - 1].second && (F > sliced_​​table[i][0].second || F < sliced_​​table[i][sliced_​​table[i].size() - 1].second))
				{
					sliced_​​table.erase(sliced_​​table.begin() + i);
					--i;
					continue;
				}
				swap_table(sliced_​​table[i]);
				sort_table(sliced_​​table[i], F);
			}
			if (sliced_​​table.size() == 0)
			{
				cout << "ФУНКЦИЯ НЕ ПРИНИМАЕТ ДАННОГО ЗАЧЕНИЯ НА ВСЁМ ПРОМЕЖУТКЕ!!!" << endl;
			}

			cout << "Введите степень интерполяционного многчлена:" << endl;
			while (n + 1 > m || n <= 1)
			{
				cin >> n;
				if (n + 1 > m)
				{
					cout << "степень интерполяционнго многочлена должна быть меньше " << m << endl;
				}
				if (n <= 1)
				{
					cout << "степень интерполяционнго многочлена должна быть больше нуля" << endl;
				}
			}

			cout << "Отсортированная таблица:" << endl;
			print_sliced_table(sliced_​​table);
			for (int i = 0; i < sliced_​​table.size(); ++i)
			{
				cout << "Значение интерполяционного многочлена в точке " << F << ", найденное при помощи представления в форме Лагранжа на учасатке №" << i + 1 << ":" << endl;
				if (n > sliced_​​table[i].size())
				{
					L = Lagrange(sliced_​​table[i], F, sliced_​​table[i].size());
				}
				else
				{
					L = Lagrange(sliced_​​table[i], F, n);
				}
				cout << L << endl;
				cout << "значение абсолютной фактической погрешности для формы Лагранжа на учасатке №" << i + 1 << ":" << endl;
				cout << abs(f(L) - F) << endl << endl;
			}
			for (int i = 0; i < sliced_​​table.size(); ++i)
			{
				sliced_​​table.erase(sliced_​​table.begin());
			}
			sliced_​​table = cut(nodes, areas);
		}
		
		if (menu == 2)
		{
			cout << "Введите значени функции:" << endl;
			cin >> F;
			root_separation(a, b, segments, areas);
			
			for (int i = 0; i < sliced_​​table.size(); ++i)
			{
				if (sliced_​​table[i][0].second < sliced_​​table[i][sliced_​​table[i].size() - 1].second && (F < sliced_​​table[i][0].second || F > sliced_​​table[i][sliced_​​table[i].size() - 1].second)
					|| sliced_​​table[i][0].second > sliced_​​table[i][sliced_​​table[i].size() - 1].second && (F > sliced_​​table[i][0].second || F < sliced_​​table[i][sliced_​​table[i].size() - 1].second))
				{
					sliced_​​table.erase(sliced_​​table.begin() + i);
					segments.erase(segments.begin() + i);
					--i;
					continue;
				}
				sort_table(sliced_​​table[i], F);
			}
			if (sliced_​​table.size() == 0)
			{
				cout << "ФУНКЦИЯ НЕ ПРИНИМАЕТ ДАННОГО ЗАЧЕНИЯ НА ВСЁМ ПРОМЕЖУТКЕ!!!" << endl;
			}
			print_segments(segments);
			cout << "counts = " << segments.size() << endl;

			cout << "Введите epsilon:" << endl;
			cin >> epsilon;
			cout << "Введите степень интерполяционного многчлена:" << endl;
			while (n + 1 > m || n <= 1)
			{
				cin >> n;
				if (n + 1 > m)
				{
					cout << "степень интерполяционнго многочлена должна быть меньше " << m << endl;
				}
				if (n <= 1)
				{
					cout << "степень интерполяционнго многочлена должна быть больше нуля" << endl;
				}
			}

			bisection(epsilon, segments, X, delta, M, sliced_​​table, F, n);
			cout << "-------------------------------------" << endl;
			cout << "Метод бисекции" << endl;
			cout << "Количество шагов: ";
			print(M);
			cout << "Приближенное решение с точность " << epsilon << ": ";
			print(X);
			cout << "длины отрезков: ";
			print(delta);
			cout << "абсолютная величина невязки: ";
			print_discrepancy(X, F);

			for (int i = 0; i < sliced_​​table.size(); ++i)
			{
				sliced_​​table.erase(sliced_​​table.begin());
				segments.erase(segments.begin());
			}
			sliced_​​table = cut(nodes, areas);
		}

		cout << "Введите 0 для выхода из программы" << endl << "1 для вычисления аргумента для значения функции первым способом" << endl << "2 для решения вторым способом" << endl;
		cin >> menu;
		while (menu != 0 && menu != 1 && menu != 2)
		{
			cout << "Введите 0 для выхода из программы" << endl <<"1 для вычисления аргумента для значения функции первым способом" << endl << "2 для решения вторым способом" << endl;
			cin >> menu;
		}
		n = m + 1;
	}

	return EXIT_SUCCESS;
}


double f(double x)
{
	return exp(x) - x;
	//return sin(x);
	//return pow(x, 5) + pow(x, 4) - 10 * pow(x, 2);
}

void table_generation(vector<pair<double, double>>& nodes, double a, double b, int m)
{
	pair<double, double> node;
	--m;
	double h = (b - a) / m;
	for (int i = 0; i < m; ++i)
	{
		node.first = a;
		node.second = f(a);
		nodes.push_back(node);
		a += h;
	}
	node.first = b;
	node.second = f(b);
	nodes.push_back(node);
}

void print_table(vector<pair<double, double>> nodes)
{
	cout << '\t' << "x" << '\t' << '\t' << '\t' << "f(x)" << endl;
	for (int i = 0; i < nodes.size(); ++i)
	{
		cout << i + 1 << ':' << '\t' << nodes[i].first << '\t' << '\t' << '\t' << nodes[i].second << endl;
	}
}

void sort_table(vector<pair<double, double>>& nodes, double x)
{
	pair<double, double> temp;
	for (int i = 0; i < nodes.size() - 1; ++i)
	{
		for (int j = 0; j < nodes.size() - i - 1; ++j)
		{
			if (abs(x - nodes[j].first) > abs(x - nodes[j + 1].first))
			{
				temp = nodes[j];
				nodes[j] = nodes[j + 1];
				nodes[j + 1] = temp;
			}
		}
	}
}

void swap_table(vector<pair<double, double>>& nodes)
{
	for (int i = 0; i < nodes.size(); ++i)
	{
		swap(nodes[i].first, nodes[i].second);
	}
}

double Lagrange(vector<pair<double, double>> nodes, double x, int n)
{
	if (n > nodes.size())
	{
		n = nodes.size() - 1;
	}
	double result = 0;
	for (int i = 0; i < n; ++i)
	{
		double a = nodes[i].second;
		for (int j = 0; j < n; ++j)
		{
			if (i != j)
			{
				a *= (x - nodes[j].first) / (nodes[i].first - nodes[j].first);
			}
		}
		result += a;
	}
	return result;
}

vector<double> Newton_coefficient(vector<pair<double, double>> nodes, int n)
{
	if (n > nodes.size())
	{
		n = nodes.size() - 1;
	}
	vector<double> coefficient;
	coefficient.push_back(nodes[0].second);
	double A = 0;
	double w = 1;
	for (int i = 1; i <= n; ++i)
	{
		for (int k = 0; k <= i; ++k)
		{
			for (int j = 0; j <= i; ++j)
			{
				if (j != k)
				{
					w *= (nodes[k].first - nodes[j].first);
				}
			}
			A += nodes[k].second / w;
			w = 1;
		}
		coefficient.push_back(A);
		A = 0;
	}
	return coefficient;
}

double Newton(vector <pair<double, double>> nodes, vector<double> coefficient, double x, int n)
{
	if (n > nodes.size())
	{
		n = nodes.size() - 1;
	}
	double result = coefficient[0];
	double a = 0;
	for (int i = 1; i < n; ++i)
	{

		a = coefficient[i];
		for (int j = 0; j < i; ++j)
		{
			a *= (x - nodes[j].first);

		}
		result += a;
	}
	return result;
}

double dist(double a, double b)
{
	return abs(a - b);
}

vector<double> areas_of_monotonicity(vector<pair<double, double>> nodes, double a, double b)
{
	int m = nodes.size()*10;
	int n = 5;
	double h = (b - a) / m;
	sort_table(nodes, a);
	vector<double> coefficient = Newton_coefficient(nodes, 5);
	bool mon = Newton(nodes, coefficient, a, 5) < Newton(nodes, coefficient, a + h, 5);
	vector<double> areas_of_monotonicity;
	int k = 0;
	vector<pair<double, double>> tab;
	for (int j = 0; j <= n; ++j)
	{
		tab.push_back(nodes[j]);
	}
	int i = 1;
	while (a + h * (i + 1) <= b)
	{
		
		if (n + 1 + k != nodes.size())
		{
			if (dist(a + h * (i + 1), nodes[n + 1 + k].first) < dist(a + h * (i + 1), nodes[k].first))
			{
				tab.erase(tab.begin());
				tab.push_back(nodes[n + 1 + k]);
				coefficient = Newton_coefficient(tab, n);
				++k;
			}
		}
		
		if (mon != Newton(tab, coefficient, a + h * i, 5) < Newton(tab, coefficient, a + h * (i + 1), 5))
		{
			areas_of_monotonicity.push_back(a + h * i + h / 2);
			cout << "areas_of_monotonicity		" << a + h * i + h / 2 << "		-----------------------------" << endl;
			mon = Newton(tab, coefficient, a + h * i, 5) < Newton(tab, coefficient, a + h * (i + 1), 5);
		}
		++i;
	}

	return areas_of_monotonicity;
}

vector<vector<pair<double, double>>> cut(vector<pair<double, double>> nodes, vector<double> areas_of_monotonicity)
{
	vector<vector<pair<double, double>>> sliced_​​table;
	int k = 0;
	vector<pair<double, double>> plus;
	sliced_​​table.push_back(plus);
	for (int i = 0; i < nodes.size(); ++i)
	{
		if (k < areas_of_monotonicity.size())
		{
			if (nodes[i].first <= areas_of_monotonicity[k])
			{
				sliced_​​table[k].push_back(nodes[i]);
				
			}
			else
			{
				sliced_​​table.push_back(plus);
				++k;
				sliced_​​table[k].push_back(nodes[i]);
			}
		}
		else
		{
			sliced_​​table[k].push_back(nodes[i]);
		}
		
	}

	return sliced_​​table;
}

void print_sliced_table(vector<vector<pair<double, double>>> sliced_​​table)
{
	for (int i = 0; i < sliced_​​table.size(); ++i)
	{
		print_table(sliced_​​table[i]);
		cout << "-------------------------------------" << endl;
	}
}

void print_segments(vector<pair<double, double>> segments)
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

void print_npk(vector<pair<double, double>> segments)
{
	for (int i = 0; i < segments.size(); ++i)
	{
		cout << (segments[i].first + segments[i].second) / 2 << " ";
	}
	cout << endl;
}

void root_separation(double a, double b, vector<pair<double, double>>& segments, vector<double> areas_of_monotonicity)
{
	pair<double, double> seg;
	seg.first = a;
	for (int i = 0; i < areas_of_monotonicity.size(); ++i)
	{
		seg.second = areas_of_monotonicity[i];
		segments.push_back(seg);
		seg.first = areas_of_monotonicity[i];
	}
	seg.second = b;
	segments.push_back(seg);
}

void bisection(double epsilon, vector<pair<double, double>> segments, vector<double>& X, vector<double>& delta, vector<int>& M, vector<vector<pair<double, double>>> sliced_table, double F, int n)
{
	double x = 0;
	double d = 0;
	double a = 0;
	double fa = 0;
	double b = 0;
	double c = 0;
	double fc = 0;
	int m = 0;
	int k = 0;
	vector<double> coefficient;

	for (int i = 0; i < segments.size(); ++i)
	{
		m = 0;
		a = segments[i].first;
		b = segments[i].second;
		while (b - a > 2 * epsilon)
		{
			c = (a + b) / 2;
			//тут можно оптимизировать
			sort_table(sliced_table[i], a);
			coefficient = Newton_coefficient(sliced_table[i], n);
			fa = Newton(sliced_table[i], coefficient, a, n) - F;
			sort_table(sliced_table[i], c);
			coefficient = Newton_coefficient(sliced_table[i], n);
			fc = Newton(sliced_table[i], coefficient, c, n) - F;
			//
			if (fa * fc <= 0)
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

void print_discrepancy(vector<double> X, double F)
{
	for (int i = 0; i < X.size(); ++i)
	{
		cout << abs(f(X[i]) - F)  << " ";
	}
	cout << endl;
}
