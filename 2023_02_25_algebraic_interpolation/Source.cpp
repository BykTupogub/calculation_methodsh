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

double Lagrange(vector <pair<double, double>> nodes, double x, int n);

double divided_difference(vector <pair<double, double>> nodes, int a, int b);

double Newton(vector <pair<double, double>> nodes, double x, int n);

int main()
{
	setlocale(LC_ALL, "Russian");
	cout << setprecision(16);

	cout << "«јƒј„ј јЋ√≈Ѕ–ј»„≈— ќ√ќ »Ќ“≈–ѕќЋ»–ќ¬јЌ»я" << endl;
	cout << "¬ј–»јЌ“ 3" << endl;
	int m = 0;
	cout << "¬ведите число значений в таблице:" << endl;
	cin >> m;
	double a = 0;
	double b = 0;
	cout << "¬ведите границы интервала:" << endl;
	cin >> a >> b;
	vector<pair<double, double>> nodes;
	table_generation(nodes, a, b, m);
	print_table(nodes);
	int menu = 1;
	double x = 0;
	int n = m + 1;
	double L = 0;
	double N = 0;
	while (menu != 0)
	{
		cout << "¬ведите точку интерполировани€:" << endl;
		cin >> x;
		cout << "¬ведите степень интерпол€ционного многчлена:" << endl;
		while (n > m || n < 1)
		{
			cin >> n;
			if (n > m)
			{
				cout << "степень интерпол€ционнго многочлена должна быть меньше" << m << endl;
			}
			if (n < 1)
			{
				cout << "степень интерпол€ционнго многочлена должна быть больше нул€" << endl;
			}
		}
		sort_table(nodes, x);
		cout << "ќтсортированна€ таблица:" << endl;
		print_table(nodes);
		cout << "«начение интерпол€ционного многочлена в точке " << x << ", найденное при помощи представлени€ в форме Ћагранжа:" << endl;
		L = Lagrange(nodes, x, n);
		cout << L << endl;
		cout << "значение абсолютной фактической погрешности дл€ формы Ћагранжа:" << endl;
		cout << abs(f(x) - L) << endl;
		cout << "«начение интерпол€ционного многочлена в точке " << x << ", найденное при помощи представлени€ в форме Ќьютона:" << endl;
		N = Newton(nodes, x, n);
		cout << N << endl;
		cout << "значение абсолютной фактической погрешности дл€ формы Ќьютона:" << endl;
		cout << abs(f(x) - N) << endl;
		cout << "¬ведите 0 дл€ выхода из программы или 1 дл€ вычислени€ многочленов в новой точке" << endl;
		cin >> menu;
		while (menu != 0 && menu != 1)
		{
			cout << "¬ведите 0 дл€ выхода из программы или 1 дл€ вычислени€ многочленов в новой точке" << endl;
			cin >> menu;
		}
	}
	
	

	return EXIT_FAILURE;
}


double f(double x)
{
	return exp(x) - x;
	//return pow(x, 5) + pow(x, 4) - 10 * pow(x, 2);
}

void table_generation(vector<pair<double, double>>& nodes, double a, double b, int m)
{
	pair<double, double> node;
	double h = (b - a) / m;
	a += h/2;
	for (int i = 0; i < m; ++i)
	{
		node.first = a;
		node.second = f(a);
		nodes.push_back(node);
		a += h;
	}
}

void print_table(vector<pair<double, double>> nodes)
{
	cout << '\t' << "x" << '\t' << '\t' << '\t' << "f(x)" << endl;
	for (int i = 0; i < nodes.size(); ++i)
	{
		cout << i+1 << ':' << '\t' <<nodes[i].first << '\t' << '\t' << '\t' << nodes[i].second << endl;
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

double Lagrange(vector<pair<double, double>> nodes, double x, int n)
{
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

double divided_difference(vector<pair<double, double>> nodes, int a, int b)
{
	if (b - a == 1)
	{
		return (nodes[b].second - nodes[a].second) / (nodes[b].first - nodes[a].first);
	}
	else
	{
		double fa = divided_difference(nodes, a, b - 1);
		double fb = divided_difference(nodes, a + 1, b);
		return (fb - fa) / (nodes[b].first - nodes[a].first);
	}
}

double Newton(vector<pair<double, double>> nodes, double x, int n)
{
	double result = nodes[0].second;
	double a = 0;
	for (int i = 1; i < n; ++i)
	{
		a = divided_difference(nodes, 0, i);
		for (int j = 0; j < i; ++j)
		{
			a *= (x - nodes[j].first);
			
		}
		result += a;
	}
	return result;
}
