#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	// Zmienne wejœciowe:
	// ff - wskaŸnik do funkcji celu
	// N - liczba zmiennych funkcji celu
	// lb, ub - dolne i górne ograniczenie
	// epslion - zak³¹dana dok³adnoœæ rozwi¹zania
	// Nmax - maksymalna liczba wywo³añ funkcji celu
	// ud1, ud2 - user data
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);									// losujemy macierz Nx1 stosuj¹c rozk³ad jednostajny na przedziale [0,1]
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);// przeskalowywujemy rozwi¹zanie do przedzia³u [lb, ub]
			Xopt.fit_fun(ff, ud1, ud2);							// obliczmy wartoœæ funkcji celu
			if (Xopt.y < epsilon)								// sprawdzmy 1. kryterium stopu
			{
				Xopt.flag = 1;									// flaga = 1 ozancza znalezienie rozwi¹zanie z zadan¹ dok³adnoœci¹
				break;
			}
			if (solution::f_calls > Nmax)						// sprawdzmy 2. kryterium stopu
			{
				Xopt.flag = 0;									// flaga = 0 ozancza przekroczenie maksymalne liczby wywo³añ funkcji celu
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}


double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2] { 0, 0 };
		int i = 0;
		double x = x0 + d;
		solution X0, X1;

		X0.x = x0;
		double f0 = m2d(X0.fit_fun(ff, ud1, ud2));
		X1.x = x;
		double f1 = m2d(X1.fit_fun(ff, ud1, ud2));
		if (f0 == f1) {
			p[0] = x0;
			p[1] = x;
			return p;
		}

		if (f1 > f0) {
			d = -d;
			x = x0 + d;
			X1.x = x;
			f1 = m2d(X1.fit_fun(ff, ud1, ud2));
			if (f1 >= f0) {
				p[0] = x;
				p[1] = x0 - d;
				return p;
			}
		}

		double xPrev;
		do {
			if (solution::f_calls > Nmax) {
				throw string("Przekroczono maksymalna liczbe wywolan\n");
			}

			i++;
			xPrev = x;
			x = x0 + pow(alpha, i) * d;
			f0 = f1;
			X1.x = x;
			f1 = m2d(X1.fit_fun(ff, ud1, ud2));
		} while (f0 > f1);

		if (d > 0) {
			p[0] = xPrev;
			p[1] = x;
			return p;
		}

		p[0] = x;
		p[1] = xPrev;
		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2) {

	try
		{

			std::vector<double> F;
			F.push_back(1.0);
			F.push_back(1.0);

			int k = 1;
			double ratio = (b - a) / epsilon;

			while (F[k] <= ratio)
			{
				F.push_back(F[k] + F[k-1]);
				k++;
			}


			int n = k;
			std::vector<double> a_vec(n), b_vec(n), c_vec(n), d_vec(n);

			a_vec[0] = a;
			b_vec[0] = b;

			c_vec[0] = b_vec[0] - (F[k-1] / F[k]) * (b_vec[0] - a_vec[0]);
			d_vec[0] = a_vec[0] + b_vec[0] - c_vec[0];


			for (int i = 0; i <= k - 3; i++)
			{

				matrix fc = ff(matrix(c_vec[i]), ud1, ud2);
				matrix fd = ff(matrix(d_vec[i]), ud1, ud2);

				if (fc(0, 0) < fd(0, 0))
				{

					a_vec[i+1] = a_vec[i];
					b_vec[i+1] = d_vec[i];
				}
				else
				{

					b_vec[i+1] = b_vec[i];
					a_vec[i+1] = c_vec[i];
				}

				c_vec[i+1] = b_vec[i+1] - (F[k-i-2] / F[k-i-1]) * (b_vec[i+1] - a_vec[i+1]);
				d_vec[i+1] = a_vec[i+1] + b_vec[i+1] - c_vec[i+1];
			}


			solution Xopt;
		double x_final = (a_vec[k-2] + b_vec[k-2]) / 2.0;

			Xopt= c_vec[k-2];

		Xopt.x = x_final;
		matrix y_final = ff(matrix(x_final), ud1, ud2);
		Xopt.y = y_final(0, 0);
		return Xopt;
		}
		catch (string ex_info) {
			throw ("solution fib(...):\n" + ex_info);
		}
	}


solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji
		solution::clear_calls();

		double c = (a + b) / 2.0;
		double d = c, dPrev = c;

		do {
			Xopt.x = matrix(a);
			double fa = m2d(Xopt.fit_fun(ff, ud1, ud2));

			Xopt.x = matrix(b);
			double fb = m2d(Xopt.fit_fun(ff, ud1, ud2));

			Xopt.x = matrix(c);
			double fc = m2d(Xopt.fit_fun(ff, ud1, ud2));

			double l = fa * (b * b - c * c) + fb * (c * c - a * a) + fc * (a * a - b * b);
			double m = fa * (b - c) + fb * (c - a) + fc * (a - b);

			if (m <= 0) {
				throw string("m <= 0");
			}

			dPrev = d;
			d = 0.5 * (l / m);
			Xopt.x = matrix(d);
			double fd = m2d(Xopt.fit_fun(ff, ud1, ud2));

			if (a < d && d < c) {
				if (fd < fc) {
					b = c;
					c = d;
				}
				else {
					a = d;
				}
			}
			else if (c < d && d < b) {
				if (fd < fc) {
					a = c;
					c = d;
				}
				else {
					b = d;
				}
			}

			if (solution::f_calls > Nmax)
				throw string("Przekroczono maksymalną liczbę wywołań funkcji celu");
			
		} while (((b - a) > epsilon) && (fabs(d - dPrev) > gamma));

		Xopt.x = matrix(d);
		Xopt.y = ff(Xopt.x, ud1, ud2);
		Xopt.flag = 0;

		return Xopt;
	}
	catch (string ex_info) {
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}

