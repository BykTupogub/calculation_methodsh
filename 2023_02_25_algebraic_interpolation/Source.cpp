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

	cout << "������ ��������������� ����������������" << endl;
	cout << "������� 3" << endl;
	int m = 0;
	cout << "������� ����� �������� � �������:" << endl;
	cin >> m;
	double a = 0;
	double b = 0;
	cout << "������� ������� ���������:" << endl;
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
		cout << "������� ����� ����������������:" << endl;
		cin >> x;
		cout << "������� ������� ����������������� ���������:" << endl;
		while (n > m || n < 1)
		{
			cin >> n;
			if (n > m)
			{
				cout << "������� ���������������� ���������� ������ ���� ������" << m << endl;
			}
			if (n < 1)
			{
				cout << "������� ���������������� ���������� ������ ���� ������ ����" << endl;
			}
		}
		sort_table(nodes, x);
		cout << "��������������� �������:" << endl;
		print_table(nodes);
		cout << "�������� ����������������� ���������� � ����� " << x << ", ��������� ��� ������ ������������� � ����� ��������:" << endl;
		L = Lagrange(nodes, x, n);
		cout << L << endl;
		cout << "�������� ���������� ����������� ����������� ��� ����� ��������:" << endl;
		cout << abs(f(x) - L) << endl;
		cout << "�������� ����������������� ���������� � ����� " << x << ", ��������� ��� ������ ������������� � ����� �������:" << endl;
		N = Newton(nodes, x, n);
		cout << N << endl;
		cout << "�������� ���������� ����������� ����������� ��� ����� �������:" << endl;
		cout << abs(f(x) - N) << endl;
		cout << "������� 0 ��� ������ �� ��������� ��� 1 ��� ���������� ����������� � ����� �����" << endl;
		cin >> menu;
		while (menu != 0 && menu != 1)
		{
			cout << "������� 0 ��� ������ �� ��������� ��� 1 ��� ���������� ����������� � ����� �����" << endl;
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
