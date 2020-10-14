#include <iostream>
#include <iomanip>

using namespace std;

const int n = 2;

//void Yravnenie(double* F, double* x);
void Jacobi(double** Jac, double* F, double* x, int n, double apr);
void Metod_Nutona(double** Jac, double* F, double* x, int n, double apr, double* dx, double eps, int iter);
void Metod_Gaussa(double** Jac, double* F, int n, double* dx, double* x);
// ====================================================================================================
int main()
{
	int iter;
	double apr = 1e-9;
	double eps = 1e-9;

	double* F;
	F = new double[n];

	double* x;
	x = new double[n];

	double** Jac;
	Jac = new double* [n];
	for (int i = 0; i < n; i++)
		Jac[i] = new double[n + 1];

	double* dx;
	dx = new double[n];

	cout << "This program can you help to calculate the Newton method" << endl;
	cout << "==================================================" << endl;
	cout << "Vvedite chislo iteracii" << endl;
	cin >> iter;

	cout << endl << "Enter " << n << " value " << endl;
	for (int j = 0; j < n; j++)
		cin >> x[j];
	cout << "==================================================" << endl;

	Metod_Nutona(Jac, F, x, n, apr, dx, eps, iter);
	cout << "==================================================" << endl;

	system("pause");
	return 0;
}

void MyYravnenie(double* F, double* x)
{
	F[0] = x[0] * x[0] - x[1] * x[1] - 1;
	F[1] = x[0] * x[1] * x[1] * x[1] - x[1] - 3;
}

void MyJacobi(double** Jac, double* F,double* x)				// work clear
{
	Jac[0][0] = 2 * x[0];
	Jac[0][1] = -2 * x[1];
	Jac[0][2] = -F[0];
	Jac[1][0] = x[1] * x[1] * x[1];
	Jac[1][1] = 3 * x[0] * x[1] * x[1] - 1;
	Jac[1][2] = -F[1];
}

void Jacobi(double** Jac, double* F, double* x, int n, double apr)
{
	//конечно-разностный способ
	double f1, f2;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			x[j] += apr;
			MyYravnenie(F, x);
			f1 = F[i];
			x[j] -= apr;
			MyYravnenie(F, x);
			f2 = F[i];
			Jac[i][j] = (f1 - f2) / apr;
			Jac[i][n] = -F[i];
		}
	}
}

void Metod_Gaussa(double** Jac, double* F, int n, double* dx, double* x)
{
	//Pryamou xod Metoda Gaysa
	for (int i = 0; i < n; i++)										//Poisk maximalnogo elementa(v pervom stolbike)
	{	      													//glavny element i stolbca
		double max = abs(Jac[i][i]);							//perestanovka i stroki c naidenoi(pereprisvaivaem novera strok)
		int my = i;
		for (int t = i; t < n; t++)
			if (abs(Jac[t][i]) > max)
			{
				max = abs(Jac[t][i]);
				my = t;
			}

		//Peremeschenie strok
		if (my != i)
		{
			double* per = Jac[i];
			Jac[i] = Jac[my];
			Jac[my] = per;
		}

		//delenie g stroki 
		double amain = Jac[i][i];
		for (int z = 0; z < n + 1; z++)
		{
			Jac[i][z] = Jac[i][z] / amain;
		}

		//vychitanie iz i=1 i dalee strok i-y * na i-y koefficient sootvetstvuy stroki
		for (int j = i + 1; j < n; j++)
		{
			double b = Jac[j][i];
			for (int z = i; z < n + 1; z++)
				Jac[j][z] = Jac[j][z] - Jac[i][z] * b;
		}
	}

	//obratny hod gausa
	for (int i = n - 1; i > 0; i--)
	{
		for (int j = i - 1; j >= 0; j--)
			Jac[j][n] = Jac[j][n] - Jac[j][i] * Jac[i][n];
	}

	for (int i = 0; i < n; ++i)
		dx[i] = Jac[i][n];

	for (int i = 0; i < n; i++)
		x[i] += Jac[i][n];
}

void Metod_Nutona(double** Jac, double* F, double* x, int n, double apr, double* dx, double eps, int iter)
{
	double D1, D2;
	double max;
	int k = 0;
	cout << "k_iter" << setw(12) << "del1" << setw(16) << "del2" << endl;
	cout << "==================================================" << endl;

	while (true)
	{
		MyYravnenie(F, x);
		/*MyJacobi(Jac,F, x);*/
		Jacobi(Jac, F, x, n, apr);
		/*матрица якоби считается верно*/
		//for (int i = 0; i < n; i++) {
		//	for (int j = 0; j < n + 1; j++) {
		//		cout << setw(5) << Jac[i][j];
		//	}
		//	cout << endl;
		//}

		for (int j = 0; j < n; j++)
			F[j] *= -1;			                        // < ----- вектор невязки 

		Metod_Gaussa(Jac, F, n, dx, x);
		max = 0;

		MyYravnenie(F, x);
		for (int q = 0; q < n; q++)
		{
			if (abs(F[q]) > max)
				max = abs(F[q]);
		}
		D1 = max;
		max = 0;

		for (int z = 0; z < n; z++)
		{
			if (abs(x[z]) < 1 && abs(dx[z]) > max)
				max = abs(dx[z]);

			if (abs(x[z]) >= 1 && abs(dx[z] / x[z]) > max)
				max = abs(dx[z] / x[z]);
		}

		D2 = max;
		cout << setw(4) << k + 1 << "\t" << setw(10) << D1 << "\t" << setw(10) << D2 << endl;
		k++;

		if (D1 <= eps && D2 <= eps || k >= iter) {
			cout << "IER = 2" << endl;
			break;
		}
	}
	cout << "==================================================" << endl;
	cout << "X1 = " << x[0] << endl;
	cout << "X2 = " << x[1] << endl;
}