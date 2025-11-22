#include"user_funs.h"

#define M_PI 3.141592

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

	return (sin(M_PI*sqrt(pow(x1/M_PI,2)+pow(x2/M_PI,2))))/(M_PI(sqrt(pow(x1/M_PI,2)+pow(x2/M_PI,2))));

}

matrix ff3T(matrix x, matrix ud1, matrix ud2) {
    double x1 = m2d(x(0, 0));
    double x2 = m2d(x(1, 0));

    double y = sin(M_PI * sqrt(pow(x1 / M_PI, 2) + pow(x2 / M_PI, 2)) / M_PI * sqrt(pow(x1 / M_PI, 2) + pow(x2 / M_PI, 2)));

    return matrix(1, 1, y);
}
