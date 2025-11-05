#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	// Zmienne wejściowe:
	// ff - wskaźnik do funkcji celu
	// N - liczba zmiennych funkcji celu
	// lb, ub - dolne i górne ograniczenie
	// epslion - zakłądana dokładność rozwiązania
	// Nmax - maksymalna liczba wywołań funkcji celu
	// ud1, ud2 - user data
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);									// losujemy macierz Nx1 stosując rozkład jednostajny na przedziale [0,1]
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);// przeskalowywujemy rozwiązanie do przedziału [lb, ub]
			Xopt.fit_fun(ff, ud1, ud2);							// obliczmy wartość funkcji celu
			if (Xopt.y < epsilon)								// sprawdzmy 1. kryterium stopu
			{
				Xopt.flag = 1;									// flaga = 1 ozancza znalezienie rozwiązanie z zadaną dokładnością
				break;
			}
			if (solution::f_calls > Nmax)						// sprawdzmy 2. kryterium stopu
			{
				Xopt.flag = 0;									// flaga = 0 ozancza przekroczenie maksymalne liczby wywołań funkcji celu
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
		solution::clear_calls();
		double* p = new double[3] { 0, 0, 0 };
		int i = 0;
		double x = x0 + d;
		solution Xopt;

		Xopt.x = x0;
		double f0 = m2d(Xopt.fit_fun(ff, ud1, ud2));
		Xopt.x = x;
		double f1 = m2d(Xopt.fit_fun(ff, ud1, ud2));
		if (f0 == f1) {
			p[0] = x0;
			p[1] = x;
			p[2] = Xopt.f_calls;
			return p;
		}

		if (f1 > f0) {
			d = -d;
			x = x0 + d;
			Xopt.x = x;
			f1 = m2d(Xopt.fit_fun(ff, ud1, ud2));
			if (f1 >= f0) {
				p[0] = x;
				p[1] = x0 - d;
				p[2] = Xopt.f_calls;
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
			Xopt.x = x;
			f1 = m2d(Xopt.fit_fun(ff, ud1, ud2));
		} while (f0 > f1);

		if (d > 0) {
			p[0] = xPrev;
			p[1] = x;
			p[2] = Xopt.f_calls;
			return p;
		}

		p[0] = x;
		p[1] = xPrev;
		p[2] = Xopt.f_calls;
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
		solution::clear_calls();
		std::vector<double> F;
		F.push_back(1.0);
		F.push_back(1.0);

		int k = 1;
		double ratio = (b - a) / epsilon;

		while (F[k] <= ratio)
		{
			F.push_back(F[k] + F[k - 1]);
			k++;
		}


		int n = k;
		std::vector<double> a_vec(n), b_vec(n), c_vec(n), d_vec(n);

		a_vec[0] = a;
		b_vec[0] = b;

		c_vec[0] = b_vec[0] - (F[k - 1] / F[k]) * (b_vec[0] - a_vec[0]);
		d_vec[0] = a_vec[0] + b_vec[0] - c_vec[0];

		solution Xopt;
		for (int i = 0; i <= k - 3; i++)
		{
			Xopt.x = c_vec[i];
			matrix fc = Xopt.fit_fun(ff, ud1, ud2);
			Xopt.x = d_vec[i];
			matrix fd = Xopt.fit_fun(ff, ud1, ud2);

			if (fc(0, 0) < fd(0, 0))
			{

				a_vec[i + 1] = a_vec[i];
				b_vec[i + 1] = d_vec[i];
			}
			else
			{

				b_vec[i + 1] = b_vec[i];
				a_vec[i + 1] = c_vec[i];
			}

			c_vec[i + 1] = b_vec[i + 1] - (F[k - i - 2] / F[k - i - 1]) * (b_vec[i + 1] - a_vec[i + 1]);
			d_vec[i + 1] = a_vec[i + 1] + b_vec[i + 1] - c_vec[i + 1];
		}

		double x_final = (a_vec[k - 2] + b_vec[k - 2]) / 2.0;
		Xopt = c_vec[k - 2];

		Xopt.x = x_final;
		matrix y_final = Xopt.fit_fun(ff, ud1, ud2);
		Xopt.y = y_final(0, 0);
		Xopt.flag = 1;
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

			if (m <= 0.0) {
				cout << "m <= 0\n";
				Xopt.x = NAN;
				Xopt.y = NAN;
				Xopt.flag = -1;
				return Xopt;
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
				throw string("Przekroczono maksymalna liczbe wywolan");
			
		} while (((b - a) > epsilon) && (fabs(d - dPrev) > gamma));

		Xopt.x = matrix(d);
		Xopt.y = ff(Xopt.x, ud1, ud2);
		Xopt.flag = 1;

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
		solution Xopt, XB;
		XB.x = x0;

		int n = get_size(XB.x)[0];
		matrix dirs = ident_mat(n);

		while (s > epsilon) {
			Xopt = HJ_trial(ff, XB, s, dirs);

			matrix f1 = Xopt.fit_fun(ff);
			matrix f2 = XB.fit_fun(ff);

			if (m2d(f1) < m2d(f2)) {
				while (m2d(f1) < m2d(f2)) {
					matrix xTmp = XB.x;
					XB.x = Xopt.x;
					Xopt.x = 2 * XB.x - xTmp;

					Xopt = HJ_trial(ff, Xopt, s, dirs);

					f1 = Xopt.fit_fun(ff);
					f2 = XB.fit_fun(ff);

					if (solution::f_calls > Nmax)
						throw string("Error f_calls > Nmax");
				}
			}
			else {
				s *= alpha;
			}

			if (solution::f_calls > Nmax)
				throw string("Error f_calls > Nmax");
		}

		XB.flag = 1;
		return XB;
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
		int n = get_size(XB.x)[0];
		matrix f0 = XB.fit_fun(ff);

		for (int i = 0; i < n; ++i) {
			matrix ei = get_col(ud1, i); 

			solution Xopt = XB;
			Xopt.x = XB.x + s * ei;
			matrix f1 = Xopt.fit_fun(ff);
			if (m2d(f1) < m2d(f0)) {
				XB.x = Xopt.x;
				f0 = f1;
			}
			else {
				Xopt.x = XB.x - s * ei;
				matrix f3 = Xopt.fit_fun(ff);
				if (m2d(f3) < m2d(f0)) {
					XB.x = Xopt.x;
					f0 = f3;
				}
			}
		}
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
		int n = get_len(x0);

		int i = 0;
		matrix d = ident_mat(n);      
		matrix lambda(n, 1);          
		matrix p(n, 1);               
		matrix xB = x0;               
		matrix s = s0;                
		matrix x_curr = x0;           

		solution XB_sol(xB);
		XB_sol.fit_fun(ff, ud1, ud2);
		double fB = m2d(XB_sol.y);

		while (true) {

			for (int j = 0; j < n; j++) {
				matrix dj = d[j];  

				matrix x_trial = xB + s(j) * dj;
				solution X_trial(x_trial);
				X_trial.fit_fun(ff, ud1, ud2);
				double f_trial = m2d(X_trial.y);

				if (f_trial < fB) {
					xB = x_trial;           
					fB = f_trial;
					lambda(j) = lambda(j) + s(j);  
					s(j) = alpha * s(j);    
				}
				else {
					s(j) = -beta * s(j);    
					p(j) = p(j) + 1;       
				}
			}


			i = i + 1;

			x_curr = xB;

			bool all_lambda_nonzero = true;
			bool all_p_nonzero = true;
			for (int j = 0; j < n; j++) {
				if (lambda(j) == 0) all_lambda_nonzero = false;
				if (p(j) == 0) all_p_nonzero = false;
			}

			if (all_lambda_nonzero && all_p_nonzero) {
				matrix d_new = ident_mat(n);

				matrix v1(n, 1);
				for (int j = 0; j < n; j++) {
					v1 = v1 + lambda(j) * d[j];
				}

				double norm_v1 = norm(v1);
				if (norm_v1 > epsilon) {
					d_new.set_col(v1 * (1.0 / norm_v1), 0);

					for (int j = 1; j < n; j++) {
						matrix vj = d[j];  

						for (int k = 0; k < j; k++) {
							matrix dk = d_new[k];
							double proj = m2d(trans(vj) * dk);
							vj = vj - proj * dk;
						}

						double norm_vj = norm(vj);
						if (norm_vj > epsilon) {
							d_new.set_col(vj * (1.0 / norm_vj), j);
						}
						else {
							matrix ej(n, 1);
							ej(j, 0) = 1.0;
							d_new.set_col(ej, j);
						}
					}

					d = d_new; 
				}

				lambda = matrix(n, 1);

				p = matrix(n, 1);

				s = s0;
			}

			if (solution::f_calls >= Nmax) {
				Xopt.x = x_curr;
				Xopt.y = fB;
				Xopt.flag = 0;  
				return Xopt;
			}

			double max_step = 0;
			for (int j = 0; j < n; j++) {
				double abs_s = s(j) >= 0 ? s(j) : -s(j); 
				if (abs_s > max_step) {
					max_step = abs_s;
				}
			}

			if (max_step < epsilon) {
				Xopt.x = x_curr;
				Xopt.fit_fun(ff, ud1, ud2);
				Xopt.flag = 1;  
				return Xopt;
			}
		}

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



