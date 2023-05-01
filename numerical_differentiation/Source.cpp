#include<iostream>
#include<cmath>
#include <iomanip> 
#include<vector>
#include<utility>

using namespace std;


double f(double x);

double f_1(double x);

double f_2(double x);

void table_generation(vector<pair<double, double>>& nodes, double a, double h, int m);

void print_table(vector <pair<double, double>> nodes);

vector<double> f_d1(vector<pair<double, double>>& node, double h);

vector<double> f_d2(vector<pair<double, double>>& node, double h);

int main()
{
	setlocale(LC_ALL, "Russian");
	cout << setprecision(16);

	cout << "Задача нахождения производных таблично-заданной функции по формулам численного дифференцирования" << endl;
	cout << "ВАРИАНТ 3" << endl;
	cout << "f(x) = exp(1,5*4*x)" << endl;
	int m = 0;
	double a = 0;
	double h = 0;
	double b = 0;
	vector<pair<double, double>> nodes;
	vector<double> d1;
	vector<double> d2;

	int menu = 3;
	while (menu != 0)
	{
		cout << "Введите число значений в таблице:" << endl;
		cin >> m;
		cout << "Введите a:" << endl;
		cin >> a;
		cout << "Введите h:" << endl;
		cin >> h;
		b = a + m * h;
		table_generation(nodes, a, h, m);
		//print_table(nodes);

		d1 = f_d1(nodes, h); 
		d2 = f_d2(nodes, h);

		for (int i = 0; i < nodes.size(); ++i)
		{
			cout << nodes[i].first << "	" << nodes[i].second << "	" << d1[i] << "	" << d2[i] << "	" << endl;
		}
		cout << endl;
		for (int i = 0; i < nodes.size(); ++i)
		{
			cout << scientific <<  abs(f_1(nodes[i].first) - d1[i]) << "	" << abs(f_2(nodes[i].first) - d2[i]) << endl;
		}


		cout << "Введите 0 для выхода из программы" << endl << "1 для ввода новых параметров" << endl;
		cin >> menu;
		while (menu != 0 && menu != 1)
		{
			cout << "Введите 0 для выхода из программы" << endl << "1 для вычисления аргумента для значения функции первым способом" << endl;
			cin >> menu;
		}
		nodes.clear();
		d1.clear();
		d2.clear();
		/*
		for (int i = 0; i < nodes.size(); ++i)
		{
			//nodes[nodes.begin()].first.erase();
			nodes.erase(nodes.begin());
			d1.erase(d1.begin());
			d2.erase(d2.begin());
		}
		*/
	}

	return EXIT_SUCCESS;
}

double f(double x)
{
	return exp(1.5*4*x);
	//return x * x / 2;
}

double f_1(double x)
{
	return 6*f(x);
	//return x;
}

double f_2(double x)
{
	return 36*f(x);
	//return 1;
}

vector<double> f_d1(vector<pair<double, double>>& node, double h)
{
	vector<double> result;
	result.push_back((-3*node[0].second + 4*node[1].second - node[2].second) / (2 * h));
	for (int i = 1; i < node.size() - 1; ++i)
	{
		result.push_back((node[i + 1].second - node[i - 1].second) / (2 * h));
	}
	result.push_back((3 * node[node.size() - 1].second - 4 * node[node.size() - 2].second + node[node.size() - 3].second) / (2 * h));

	return result;
}

vector<double> f_d2(vector<pair<double, double>>& node, double h)
{
	vector<double> result;
	result.push_back((2 * node[0].second - 5 * node[1].second + 4 * node[2].second - node[3].second) / (h * h));
	for (int i = 1; i < node.size() - 1; ++i)
	{
		result.push_back((node[i + 1].second - 2 * node[i].second + node[i - 1].second) / (h * h));
	}
	result.push_back((2 * node[node.size() - 1].second - 5 * node[node.size() - 2].second + 4 * node[node.size() - 3].second - node[node.size() - 4].second) / (h * h));

	return result;
}

void table_generation(vector<pair<double, double>>& nodes, double a, double h, int m)
{
	pair<double, double> node;
	--m;
	for (int i = 0; i < m + 1; ++i)
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
		cout << i + 1 << ':' << '\t' << nodes[i].first << '\t' << '\t' << '\t' << nodes[i].second << endl;
	}
}
