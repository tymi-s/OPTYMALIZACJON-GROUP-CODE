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

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2, ofstream& path_file)
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

			// Zapisz sprawdzany punkt
			path_file << Xopt.x(0) << ";" << Xopt.x(1) << ";" << f1(0) << "\n";

			if (m2d(f1) < m2d(f0)) {
				XB.x = Xopt.x;
				f0 = f1;
			}
			else {
				Xopt.x = XB.x - s * ei;
				matrix f3 = Xopt.fit_fun(ff);

				// Zapisz sprawdzany punkt
				path_file << Xopt.x(0) << ";" << Xopt.x(1) << ";" << f3(0) << "\n";

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

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		ofstream path_file("sciezka_HJ.csv", ios::app);

		solution Xopt, XB;
		XB.x = x0;

		// Zapisz punkt startowy
		XB.fit_fun(ff);
		path_file << XB.x(0) << ";" << XB.x(1) << ";" << XB.y(0) << "\n";

		int n = get_size(XB.x)[0];
		matrix dirs = ident_mat(n);

		while (s > epsilon) {
			Xopt = HJ_trial(ff, XB, s, dirs, ud2, path_file);  // path_file na końcu

			matrix f1 = Xopt.fit_fun(ff);
			matrix f2 = XB.fit_fun(ff);

			if (m2d(f1) < m2d(f2)) {
				while (m2d(f1) < m2d(f2)) {
					matrix xTmp = XB.x;
					XB.x = Xopt.x;

					// Zapisz nowy punkt bazowy
					path_file << XB.x(0) << ";" << XB.x(1) << ";" << f1(0) << "\n";

					Xopt.x = 2 * XB.x - xTmp;

					Xopt = HJ_trial(ff, Xopt, s, dirs, ud2, path_file);  // path_file na końcu

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

		path_file.close();
		XB.flag = 1;
		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		ofstream path_file("sciezka_Rosen.csv", ios::app);

		solution::clear_calls();
		solution Xopt;
		int n = get_len(x0);

		int i = 0;
		matrix d = ident_mat(n);
		matrix lambda(n, 1);
		matrix p(n, 1);
		matrix s = s0;
		matrix xB = x0;

		// Zapisz punkt startowy
		solution XB_start(xB);
		XB_start.fit_fun(ff, ud1, ud2);
		path_file << xB(0) << ";" << xB(1) << ";" << XB_start.y(0) << "\n";

		while (true) {
			for (int j = 0; j < n; j++) {
				matrix dj = d[j];

				solution X_trial(xB + s(j) * dj);
				X_trial.fit_fun(ff, ud1, ud2);
				double f_trial = m2d(X_trial.y);

				// Zapisz sprawdzany punkt
				path_file << X_trial.x(0) << ";" << X_trial.x(1) << ";" << X_trial.y(0) << "\n";

				solution XB_sol(xB);
				XB_sol.fit_fun(ff, ud1, ud2);
				double fB = m2d(XB_sol.y);

				if (f_trial < fB) {
					xB = xB + s(j) * dj;
					lambda(j) = lambda(j) + s(j);
					s(j) = alpha * s(j);
				}
				else {
					s(j) = -beta * s(j);
					p(j) = p(j) + 1;
				}
			}

			i = i + 1;

			bool all_lambda_nonzero = true;
			bool all_p_nonzero = true;
			for (int j = 0; j < n; j++) {
				if (lambda(j) == 0) all_lambda_nonzero = false;
				if (p(j) == 0) all_p_nonzero = false;
			}

			if (all_lambda_nonzero && all_p_nonzero) {
				matrix Q(n, n);

				for (int row = 0; row < n; row++) {
					for (int col = 0; col < n; col++) {
						if (col == 0) {
							Q(row, col) = lambda(row);
						}
						else {
							Q(row, col) = lambda(row);
						}
					}
				}

				matrix Lambda_matrix(n, n);
				for (int row = 0; row < n; row++) {
					for (int col = 0; col < n; col++) {
						if (col <= row) {
							Lambda_matrix(row, col) = lambda(row);
						}
						else {
							Lambda_matrix(row, col) = 0.0;
						}
					}
				}

				Q = d * Lambda_matrix;

				matrix d_new = ident_mat(n);

				matrix v1 = get_col(Q, 0);
				double norm_v1 = norm(v1);
				if (norm_v1 > epsilon) {
					d_new.set_col(v1 * (1.0 / norm_v1), 0);
				}
				else {
					matrix e1(n, 1);
					e1(0) = 1.0;
					d_new.set_col(e1, 0);
				}

				for (int j = 1; j < n; j++) {
					matrix vj = get_col(Q, j);

					for (int k = 0; k < j; k++) {
						matrix dk = get_col(d_new, k);
						double proj = m2d(trans(vj) * dk);
						vj = vj - proj * dk;
					}

					double norm_vj = norm(vj);
					if (norm_vj > epsilon) {
						d_new.set_col(vj * (1.0 / norm_vj), j);
					}
					else {
						matrix ej(n, 1);
						ej(j) = 1.0;
						d_new.set_col(ej, j);
					}
				}

				d = d_new;

				lambda = matrix(n, 1);
				p = matrix(n, 1);
				s = s0;
			}

			if (solution::f_calls > Nmax) {
				Xopt.x = xB;
				Xopt.fit_fun(ff, ud1, ud2);
				path_file.close();
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
				Xopt.x = xB;
				Xopt.fit_fun(ff, ud1, ud2);
				path_file.close();
				Xopt.flag = 1;
				return Xopt;
			}
		}

		path_file.close();
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
        int i =0;
        solution x_prev, x_curr;
        double c_current = c;

        while(true){
            i++;

            matrix params = ud1;

            if( get_size(params)[0] == 0 || get_size(params)[1]==0){
                params = matrix(1,1);
            }
            params(0) = c_current;
            double s = 0.5;           // długość boku sympleksu początkowego
            double alpha = 1.0;       // współczynnik odbicia
            double beta = 0.5;        // współczynnik zawężenia
            double gamma = 2.0;       // współczynnik ekspansji
            double delta = 0.5;       // współczynnik redukcji

            x_curr = sym_NM(ff, x_prev.x, s, alpha, beta, gamma, delta, epsilon, Nmax, params, ud2);
            Xopt = x_curr;

            if(solution::f_calls > Nmax){
                Xopt.flag = 0;
                return Xopt;

            }

            // warunek stopu:
            double distance = norm(x_curr.x - x_prev.x);
            if(distance < epsilon){
                Xopt.flag =1;
                return Xopt;
            }

            // aktualizacja współczynnika kary
            c_current = dc * c_current;

            //przejscie do nastepnej interacji
            x_prev= x_curr;

        }
        return Xopt;


	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution pen_rzeczywisty(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
        solution Xopt;
        solution x_prev, x_curr;
        x_prev.x = x0;
        double c_current = c;

        while(true){
            // Rozszerz ud1 o współczynnik kary
            matrix params(get_size(ud1)[0] + 1, 1);
            for (int i = 0; i < get_size(ud1)[0]; i++) {
                params(i) = ud1(i);
            }
            params(get_size(ud1)[0]) = c_current;  // Dodaj c_current na końcu

            double s = 0.1;
            double alpha = 1.0;
            double beta = 0.5;
            double gamma = 2.0;
            double delta = 0.5;

            x_curr = sym_NM(ff, x_prev.x, s, alpha, beta, gamma, delta, epsilon, Nmax, params, ud2);
            Xopt = x_curr;

            if(solution::f_calls > Nmax){
                Xopt.flag = 0;
                return Xopt;
            }

            // warunek stopu:
            double distance = norm(x_curr.x - x_prev.x);
            if(distance < epsilon){
                Xopt.flag = 1;
                return Xopt;
            }

            // aktualizacja współczynnika kary
            c_current = dc * c_current;

            // przejście do następnej iteracji
            x_prev = x_curr;
        }

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

		int n = get_len(x0); 
		int N_v = n + 1;     

		vector<solution> P(N_v);

		P[0].x = x0;
		P[0].fit_fun(ff, ud1, ud2); 

		for (int i = 1; i <= n; ++i)
		{
			matrix ei(n, 1, 0.0);
			ei(i - 1, 0) = 1.0;

			P[i].x = P[0].x + matrix(s) * ei;

			P[i].fit_fun(ff, ud1, ud2);
		}

		do
		{
			int i_min = 0;
			int i_max = 0;

			for (int i = 1; i < N_v; ++i)
			{
				if (m2d(P[i].y) < m2d(P[i_min].y))
					i_min = i;
				if (m2d(P[i].y) > m2d(P[i_max].y))
					i_max = i;
			}

			solution p_min = P[i_min];
			solution p_max = P[i_max];

			matrix sum_p(n, 1, 0.0);
			for (int i = 0; i < N_v; ++i)
			{
				if (i != i_max)
					sum_p = sum_p + P[i].x;
			}
			matrix p_bar_x = matrix(1.0 / n) * sum_p;

			solution p_odb;
			p_odb.x = p_bar_x + matrix(alpha) * (p_bar_x - p_max.x);
			p_odb.fit_fun(ff, ud1, ud2);

			double f_odb = m2d(p_odb.y);
			double f_min = m2d(p_min.y);
			double f_max = m2d(p_max.y);

			if (f_odb < f_min)
			{
				solution p_e;
				p_e.x = p_bar_x + matrix(gamma) * (p_odb.x - p_bar_x);
				p_e.fit_fun(ff, ud1, ud2);

				if (m2d(p_e.y) < f_odb)
				{
					P[i_max] = p_e; 
				}
				else
				{
					P[i_max] = p_odb; 
				}
			}
			else 
			{
				if (f_min <= f_odb && f_odb < f_max)
				{
					P[i_max] = p_odb; 
				}
				else 
				{
					solution p_z;
					p_z.x = p_bar_x + matrix(beta) * (p_max.x - p_bar_x);
					p_z.fit_fun(ff, ud1, ud2);

					if (m2d(p_z.y) >= f_max)
					{
						for (int i = 0; i < N_v; ++i)
						{
							if (i != i_min)
							{
								P[i].x = matrix(delta) * (P[i].x + P[i_min].x);
								P[i].fit_fun(ff, ud1, ud2);
							}
						}
					}
					else 
					{
						P[i_max] = p_z; 
					}
				}
			}

			if (solution::f_calls > Nmax)
			{
				Xopt = P[i_min]; 
				Xopt.flag = 1;   
				return Xopt;
			}

			i_min = 0;
			for (int i = 1; i < N_v; ++i)
			{
				if (m2d(P[i].y) < m2d(P[i_min].y))
					i_min = i;
			}
			p_min = P[i_min];

			double max_dist = 0.0;
			for (int i = 0; i < N_v; ++i)
			{
				double d = norm(P[i].x - p_min.x);
				if (d > max_dist)
					max_dist = d;
			}

			if (max_dist < epsilon)
			{
				Xopt = p_min;
				Xopt.flag = 0; 
				return Xopt; 
			}

		} while (true); 

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
		matrix x = x0;
		matrix x_prev;
		matrix  d;
		double norm;
		double h = h0;
       do {
        	x_prev = x;
            d = -gf(x,ud1,ud2);
        	x = x + h * d;

        	double diff0 = x(0) - x_prev(0);
        	double diff1 = x(1) - x_prev(1);
        	 norm = sqrt(diff0*diff0 + diff1*diff1);

        	if (solution::f_calls > Nmax) {
        		Xopt.flag = 0;
        		return Xopt;
        	}
        }  while ( norm >= epsilon);
        Xopt = x;
		Xopt.flag = 1;
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
		matrix x = x0;
		matrix x_prev;

		matrix g = gf(x0,ud1,ud2);
		double B_B = pow(g(0), 2) + pow(g(1), 2);
		matrix  d = -g;
		double norm;
		double h = h0;
		double beta;
		do {
			x_prev = x;
			x = x + h * d;

			matrix g_new = gf(x,ud1,ud2);
			double B_T = pow(g_new(0), 2) + pow(g_new(1), 2);
			beta = B_T/ B_B;
			d = -g_new + beta * d;
			B_B=B_T;
			double diff0 = x(0) - x_prev(0);
			double diff1 = x(1) - x_prev(1);
			norm = sqrt(diff0*diff0 + diff1*diff1);

			if (solution::f_calls > Nmax) {
				Xopt.flag = 0;
				return Xopt;
			}
		}  while ( norm >= epsilon);
		Xopt = x;
		Xopt.flag = 1;
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
		matrix x = x0;
		matrix x_prev;
		matrix  d,g,H;
		double norm;
		double h = h0;
		do {
			x_prev = x;
			g = gf(x,ud1,ud2);
			H = Hf(x,ud1,ud2);
			double det = H(0,0)*H(1,1) - H(0,1)*H(1,0);
			matrix H_inv(2, 2);
			H_inv(0,0) = H(1,1) / det;
			H_inv(0,1) = -H(0,1) / det;
			H_inv(1,0) = -H(1,0) / det;
			H_inv(1,1) = H(0,0) / det;
			d = -(H_inv * g);

			x = x + h * d;

			double diff0 = x(0) - x_prev(0);
			double diff1 = x(1) - x_prev(1);
			norm = sqrt(diff0*diff0 + diff1*diff1);

			if (solution::f_calls > Nmax) {
				Xopt.flag = 0;
				return Xopt;
			}
		}  while ( norm >= epsilon);
		Xopt = x;
		Xopt.flag = 1;
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
		double alpha = (pow(5,0.5)-1)/2;
		double c = b -alpha*(b-a);
		double d = a + alpha*(b -a);
		matrix X_c(2,1);
		X_c(1) =c;
		X_c(2) =0.0;
		matrix X_b(2,1);
		X_b(1) =b;
		X_b(2) =0.0;
		while (b-a<epsilon) {
			if (ff4T(X_c,ud1,ud2)<ff4T(X_b,ud1,ud2)) {
				b = d;
				d = c;
				c = b - alpha*(b - a);
			}
			else {
				a = c;
				c = d;
				d = a + alpha*(b-a);
			}

			if (solution::f_calls > Nmax) {
				Xopt.flag = 0;
				return Xopt;
			}
		}
		Xopt = (a+b)/2;
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






