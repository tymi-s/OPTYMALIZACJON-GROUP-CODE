#include"user_funs.h"
#include <cmath>

matrix ff0T(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla przypadku testowego
{
	matrix y;												// y zawiera warto funkcji celu
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);		// ud1 zawiera wsp馧rz璠ne szukanego optimum
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla problemu rzeczywistego
{
	matrix y;												// y zawiera warto funkcji celu
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz靖kowe
		MT = matrix(2, new double[2] { m2d(x), 0.5 });		// MT zawiera moment si造 dzia豉j鉍y na wahad這 oraz czas dzia豉nia
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);	// rozwi頊ujemy r闚nanie r騜niczkowe
	int n = get_len(Y[0]);									// d逝go rozwi頊ania
	double teta_max = Y[1](0, 0);							// szukamy maksymalnego wychylenia wahad豉
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));							// warto funkcji celu (ud1 to za這穎ne maksymalne wychylenie)
	Y[0].~matrix();											// usuwamy z pami璚i rozwi頊anie RR
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);										// definiujemy wektor pochodnych szukanych funkcji
	double m = 1, l = 0.5, b = 0.5, g = 9.81;				// definiujemy parametry modelu
	double I = m * pow(l, 2);
	dY(0) = Y(1);																// pochodna z po這瞠nia to pr璠ko
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;	// pochodna z pr璠koi to przyspieszenie
	return dY;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2) {
	return -cos(0.1 * m2d(x)) * exp(-pow(0.1 * m2d(x) - 2 * M_PI, 2)) + 0.002 * (0.1 * m2d(x) * 0.1 * m2d(x));

}

matrix ff2T(matrix X, matrix ud1, matrix ud2)
{
	try {
		double x1 = m2d(X(0, 0));
		double x2 = m2d(X(1, 0));

		double y = pow(x1, 2) + pow(x2, 2) - cos(2.5 * M_PI * x1) - cos(2.5 * M_PI * x2) + 2;

		matrix F(1, 1, y);
		return F;
	}
	catch (string ex_info) {
		throw ("matrix ff2T(...):\n" + ex_info);
	}
}

matrix df(double t, matrix Y, matrix ud1, matrix ud2) {
    double alpha = Y(0);
    double omega = Y(1);

    double k1 = ud1(0);
    double k2 = ud1(1);

    double l = 2.0;
    double mr = 1.0;
    double mc = 5.0;
    double b = 0.25;
    double alpha_ref = 3.14159265358979323846;
    double omega_ref = 0.0;

    double I = (1.0 / 3.0) * mr * l * l + mc * l * l;
    double M = k1 * (alpha_ref - alpha) + k2 * (omega_ref - omega);

    matrix dY(2, 1);
    dY(0) = omega;
    dY(1) = (M - b * omega) / I;

    return dY;
}


matrix ff2R(matrix x, matrix ud1, matrix ud2) {
    double k1 = x(0);
    double k2 = x(1);
    double t0 = 0.0;
    double tend = 100.0;
    double dt = 0.1;
    double l = 2.0;
    double mr = 1.0;
    double mc = 5.0;
    double b = 0.25;
    double alpha_ref = 3.14159265358979323846;
    double omega_ref = 0.0;
    double I = (1.0 / 3.0) * mr * l * l + mc * l * l;

    matrix Y0(2, 1);
    Y0(0) = 0.0;
    Y0(1) = 0.0;

    matrix k_params(2, 1);
    k_params(0) = k1;
    k_params(1) = k2;

    matrix* S = solve_ode(df, t0, dt, tend, Y0, k_params, ud2);

    double Q = 0.0;
    int N = get_size(S[0])[0];

    for (int i = 0; i < N; i++) {
        double alpha = S[1](i, 0);
        double omega = S[1](i, 1);
        double M = k1 * (alpha_ref - alpha) + k2 * (omega_ref - omega);

        double integrand = 10.0 * (alpha_ref - alpha) * (alpha_ref - alpha) +
            (omega_ref - omega) * (omega_ref - omega) +
            M * M;

        Q += integrand * dt;
    }

    delete[] S;
    return Q;
}

matrix ff3T(matrix x, matrix ud1, matrix ud2) {

	double x1 = m2d(x(0, 0));
	double x2 = m2d(x(1, 0));

	return (sin(M_PI*sqrt(pow(x1/M_PI,2)+pow(x2/M_PI,2))))/(M_PI*(sqrt(pow(x1/M_PI,2)+pow(x2/M_PI,2))));

}


/// ta funckcja opisuje układ równań różniczkowych ruhu piłki z efektem Magnusa
matrix df3R(double t, matrix Y, matrix ud1, matrix ud2)
{
    // Y(0) = x (pozycja pozioma)
    // Y(1) = y (pozycja pionowa)
    // Y(2) = vx (prÄ™dkoÅ›Ä‡ pozioma)
    // Y(3) = vy (prÄ™dkoÅ›Ä‡ pionowa)

    // ud1 zawiera parametry fizyczne: [v0x, omega, m, r, C, rho, g]
    double v0x = ud1(0);
    double omega = ud1(1);
    double m = ud1(2);
    double r = ud1(3);
    double C = ud1(4);
    double rho = ud1(5);
    double g = ud1(6);

    double vx = Y(2);
    double vy = Y(3);

    double S = M_PI * r * r;  // pole powierzchni

    // PrÄ™dkoÅ›Ä‡ caÅ‚kowita
    double v = sqrt(vx * vx + vy * vy);

    // SiÅ‚a oporu (w kierunku przeciwnym do ruchu)
    double Dx = 0.0;
    double Dy = 0.0;

    if (v > 1e-10) {  // Unikaj dzielenia przez zero
        Dx = 0.5 * C * rho * S * vx * abs(vx);  // Proporcjonalna do vx
        Dy = 0.5 * C * rho * S * vy * abs(vy);  // Proporcjonalna do vy
    }

    // SiÅ‚a Magnusa (prostopadÅ‚a do prÄ™dkoÅ›ci i osi rotacji)
    double FMx = rho * vy * omega * M_PI * pow(r, 3);
    double FMy = rho * vx * omega * M_PI * pow(r, 3);  // â† UWAGA: zmiana znaku!

    matrix dY(4, 1);
    dY(0) = vx;  // dx/dt = vx
    dY(1) = vy;  // dy/dt = vy
    dY(2) = -(Dx + FMx) / m;  // dvx/dt
    dY(3) = -(Dy + FMy) / m - g;  // dvy/dt

    return dY;
}

/// funkcja celu do problemu rzeczywistego lab3
// Funkcja pomocnicza do obliczania tylko x_end i x_at_y50 BEZ kary
matrix ff3R_base(matrix x, matrix ud1, matrix ud2)
{
    try {
        double v0x = x(0);
        double omega = x(1);

        // Parametry fizyczne z ud1
        double m = ud1(0);
        double r = ud1(1);
        double C = ud1(2);
        double rho = ud1(3);
        double g = ud1(4);

        // Warunki początkowe: [x, y, vx, vy]
        matrix Y0(4, 1);
        Y0(0) = 0.0;
        Y0(1) = 100.0;
        Y0(2) = v0x;
        Y0(3) = 0.0;

        // Parametry dla df3R
        matrix params(7, 1);
        params(0) = v0x;
        params(1) = omega;
        params(2) = m;
        params(3) = r;
        params(4) = C;
        params(5) = rho;
        params(6) = g;

        // Symulacja
        matrix* Y = solve_ode(df3R, 0.0, 0.01, 7.0, Y0, params, ud2);

        int n = get_size(Y[0])[0];
        double x_end = 0.0;
        double x_at_y50 = -1.0;
        bool found_y50 = false;
        bool found_ground = false;

        for (int i = 0; i < n - 1; i++) {
            double y_curr = Y[1](i, 1);
            double y_next = Y[1](i + 1, 1);

            if (!found_y50 && y_curr >= 50.0 && y_next < 50.0) {
                x_at_y50 = Y[1](i, 0);
                found_y50 = true;
            }

            if (y_curr > 0.0 && y_next <= 0.0) {
                double t_interp = y_curr / (y_curr - y_next);
                x_end = Y[1](i, 0) + t_interp * (Y[1](i + 1, 0) - Y[1](i, 0));
                found_ground = true;
                break;
            }
        }

        if (!found_ground && n > 0) {
            x_end = Y[1](n - 1, 0);
        }

        delete[] Y;

        // Zwróć [x_end, x_at_y50, found_y50]
        matrix result(3, 1);
        result(0) = x_end;
        result(1) = x_at_y50;
        result(2) = found_y50 ? 1.0 : 0.0;

        return result;
    }
    catch (string ex_info) {
        throw("matrix ff3R_base(...):\n" + ex_info);
    }
}

// Funkcja celu Z KARĄ
matrix ff3R(matrix x, matrix ud1, matrix ud2)
{
    try {
        // Oblicz podstawowe wartości
        matrix base_result = ff3R_base(x, ud1, ud2);

        double x_end = base_result(0);
        double x_at_y50 = base_result(1);
        bool found_y50 = (base_result(2) > 0.5);

        // Współczynnik kary
        double c_penalty = ud1(5);

        // Funkcja kary
        double penalty = 0.0;

        if (found_y50) {
            // Kara za wyjście poza przedział [3, 7]
            if (x_at_y50 < 3.0) {
                penalty += 1000.0 * pow(3.0 - x_at_y50, 2);
            } else if (x_at_y50 > 7.0) {
                penalty += 100.0 * pow(x_at_y50 - 5.0, 2);
            }

            // Dodatkowa kara za odległość od centrum (x=5)
            // Im dalej od 5, tym większa kara
            penalty += 1000.0 * pow(x_at_y50 - 5.0, 2);

        } else {
            penalty += 1000000.0;
        }

        // Minimalizujemy: -x_end + kara
        return matrix(-x_end + c_penalty * penalty);
    }
    catch (string ex_info) {
        throw("matrix ff3R(...):\n" + ex_info);
    }


matrix ff4T(matrix x, matrix ud1, matrix ud2) {

    double x1= x(0);
    double x2 = x(1);
    return 1.0/6.0 * pow(x1,6) - 1.05*pow(x1,4) + 2*x1*x1 + x2*x2 + x1*x2;
}

matrix gradient(matrix x,matrix ud1,matrix ud2) {
    double x1 = x(0);
    double x2 = x(1);

    // 1. Obliczamy pochodne cząstkowe (sam gradient)
    double g1 = pow(x1, 5) - 4.2 * pow(x1, 3) + 4 * x1 + x2;
    double g2 = x1 + 2 * x2;

    // 2. Tworzymy wektor wynikowy d
    matrix d(2, 1);


    d(0) = g1;
    d(1) = g2;

    return d;
}

matrix hessian(matrix x, matrix ud1, matrix ud2) {
    double x1 = x(0);
    double x2 = x(1);

    // Tworzymy macierz 2x2
    matrix H(2, 2);

    // Wypełniamy wartościami wg wzorów wyżej
    H(0, 0) = 5 * pow(x1, 4) - 12.6 * pow(x1, 2) + 4; // d2f / dx1^2
    H(0, 1) = 1.0;                                    // d2f / dx1 dx2
    H(1, 0) = 1.0;                                    // d2f / dx2 dx1 (symetryczna)
    H(1, 1) = 2.0;                                    // d2f / dx2^2

    return H;
}

// matrix gradeint_spr(matrix x) {
//     double x1 = x(0);
//     double x2 = x(1);
//
//     double g1 =
//     matrix d(2, 1);
//
//
// }













