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
matrix ff3R(matrix x, matrix ud1, matrix ud2) {
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
        throw ("matrix ff3R(...):\n" + ex_info);
    }

}
matrix ff4T(matrix x, matrix ud1, matrix ud2) {

    double x1= x(0);
    double x2 = x(1);
    return 1.0/6.0 * pow(x1,6) - 1.05*pow(x1,4) + 2*x1*x1 + x2*x2 + x1*x2;
}

matrix gradient(matrix x,matrix ud1,matrix ud2) {
    double x1 = x(0);
    double x2 = x(1);


    double g1 = pow(x1, 5) - 4.2 * pow(x1, 3) + 4 * x1 + x2;
    double g2 = x1 + 2 * x2;


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

matrix f_line(matrix h, matrix xk, matrix dk) {
    return ff4T(xk + m2d(h) * dk, NAN, NAN);
}

matrix ff4R(matrix theta, matrix ud1, matrix ud2) {
    // ud1 = macierz X (3x100) - dane wejściowe
    // ud2 = macierz Y (1x100) - etykiety

    int m = get_size(ud2)[1]; // liczba próbek (100)
    double cost = 0.0;

    for (int i = 0; i < m; i++) {
        matrix xi = ud1[i]; // kolumna i z macierzy X
        double yi = ud2(0, i);

        // h_theta(xi) = 1 / (1 + e^(-theta^T * xi))
        double z = m2d(trans(theta) * xi);
        double h = 1.0 / (1.0 + exp(-z));

        // Dodaj małą wartość aby uniknąć log(0)
        h = max(h, 1e-15);
        h = min(h, 1.0 - 1e-15);

        cost += yi * log(h) + (1 - yi) * log(1 - h);
    }

    return matrix(-cost / m);
}

matrix gradient4R(matrix theta, matrix ud1, matrix ud2) {
    // ud1 = macierz X (3x100)
    // ud2 = macierz Y (1x100)

    int m = get_size(ud2)[1]; // 100
    int n = get_len(theta);    // 3

    matrix grad(n, 1, 0.0);

    for (int i = 0; i < m; i++) {
        matrix xi = ud1[i];
        double yi = ud2(0, i);

        double z = m2d(trans(theta) * xi);
        double h = 1.0 / (1.0 + exp(-z));

        grad = grad + (h - yi) * xi;
    }

    return grad * (1.0 / m);
}
void load_data(matrix& X, matrix& Y) {
    ifstream xfile("XData.txt");
    ifstream yfile("YData.txt");

    if (!xfile.good() || !yfile.good()) {
        throw string("Nie można otworzyć plików danych");
    }

    X = matrix(3, 100);
    Y = matrix(1, 100);

    string line, value;




    getline(xfile, line);
    for (int j = 0; j < 100; j++) {
        X(0, j) = 1.0;
    }

    // Wiersz 2: przedmiot 1
    getline(xfile, line);
    stringstream ss1(line);
    for (int j = 0; j < 100; j++) {
        getline(ss1, value, ';');
        // Usuń białe znaki
        value.erase(remove(value.begin(), value.end(), ' '), value.end());
        if (!value.empty()) {
            X(1, j) = stod(value);
        }
    }


    getline(xfile, line);
    stringstream ss2(line);
    for (int j = 0; j < 100; j++) {
        getline(ss2, value, ';');
        value.erase(remove(value.begin(), value.end(), ' '), value.end());
        if (!value.empty()) {
            X(2, j) = stod(value);
        }
    }


    getline(yfile, line);
    stringstream ss3(line);
    for (int j = 0; j < 100; j++) {
        getline(ss3, value, ';');
        value.erase(remove(value.begin(), value.end(), ' '), value.end());
        if (!value.empty()) {
            Y(0, j) = stod(value);
        }
    }

    xfile.close();
    yfile.close();

    //  sprawdzenie czy sa poprawnie pliki wczytae
    cout << "Pierwsze 5 probek:" << endl;
    for (int i = 0; i < 5; i++) {
        cout << "X[" << i << "] = [" << X(0,i) << ", " << X(1,i)
             << ", " << X(2,i) << "], Y=" << Y(0,i) << endl;
    }

    cout << "\nOstatnie 3 probki:" << endl;
    for (int i = 97; i < 100; i++) {
        cout << "X[" << i << "] = [" << X(0,i) << ", " << X(1,i)
             << ", " << X(2,i) << "], Y=" << Y(0,i) << endl;
    }
}

double calculate_accuracy(matrix theta, matrix X, matrix Y) {
    int m = get_size(Y)[1]; // 100
    int correct = 0;

    for (int i = 0; i < m; i++) {
        matrix xi = X[i];
        double yi = Y(0, i);

        double z = m2d(trans(theta) * xi);
        double h = 1.0 / (1.0 + exp(-z));

        int prediction = (h >= 0.5) ? 1 : 0;

        if (prediction == (int)yi) {
            correct++;
        }
    }

    return (100.0 * correct) / m;
}


//lab5
static std::vector<double> w(101);
static int wi = -1;

void setW(int i) {
    wi = i;
}

void makeW() {
    for (int i = 0; i < 101; i++) {
        w[i] = static_cast<double>(i) * 0.01;
    }
}

double getW() {
    return w[wi];
}

matrix ff5T1_1(matrix x, matrix ud1, matrix ud2) {
    if (isnan(ud2(0)))
        return (pow(x(0) - 3.0, 2) + pow(x(1) - 3.0, 2));
    else
        return (pow(ud1(0) + x(0) * ud2(0) - 3.0, 2) + pow(ud1(1) + x(0) * ud2(1) - 3.0, 2));
}

matrix ff5T2_1(matrix x, matrix ud1, matrix ud2) {
    if (isnan(ud2(0)))
        return (pow(x(0) + 3.0, 2) + pow(x(1) + 3.0, 2));
    else
        return (pow(ud1(0) + x(0) * ud2(0) + 3.0, 2) + pow(ud1(1) + x(0) * ud2(1) + 3.0, 2));
}

matrix ff5T3_1(matrix x, matrix ud1, matrix ud2) {
    return w[wi] * ff5T1_1(x, ud1, ud2) + (1.0 - w[wi]) * ff5T2_1(x, ud1, ud2);
}

matrix ff5T1_10(matrix x, matrix ud1, matrix ud2) {
    if (isnan(ud2(0)))
        return 10.0 * (pow(x(0) - 3.0, 2) + pow(x(1) - 3.0, 2));
    else
        return 10.0 * (pow(ud1(0) + x(0) * ud2(0) - 3.0, 2) + pow(ud1(1) + x(0) * ud2(1) - 3.0, 2));
}

matrix ff5T2_10(matrix x, matrix ud1, matrix ud2) {
    if (isnan(ud2(0)))
        return 0.1 * (pow(x(0) + 3.0, 2) + pow(x(1) + 3.0, 2));
    else
        return 0.1 * (pow(ud1(0) + x(0) * ud2(0) + 3.0, 2) + pow(ud1(1) + x(0) * ud2(1) + 3.0, 2));
}
matrix ff5T3_10(matrix x, matrix ud1, matrix ud2) {
    return w[wi] * ff5T1_10(x, ud1, ud2) + (1.0 - w[wi]) * ff5T2_10(x, ud1, ud2);
}

matrix ff5T1_100(matrix x, matrix ud1, matrix ud2) {
    if (isnan(ud2(0)))
        return 100.0 * (pow(x(0) - 3.0, 2) + pow(x(1) - 3.0, 2));
    else
        return 100.0 * (pow(ud1(0) + x(0) * ud2(0) - 3.0, 2) + pow(ud1(1) + x(0) * ud2(1) - 3.0, 2));
}

matrix ff5T2_100(matrix x, matrix ud1, matrix ud2) {
    if (isnan(ud2(0)))
        return 0.01 * (pow(x(0) + 3.0, 2) + pow(x(1) + 3.0, 2));
    else
        return 0.01 * (pow(ud1(0) + x(0) * ud2(0) + 3.0, 2) + pow(ud1(1) + x(0) * ud2(1) + 3.0, 2));
}

matrix ff5T3_100(matrix x, matrix ud1, matrix ud2) {
    return w[wi] * ff5T1_100(x, ud1, ud2) + (1.0 - w[wi]) * ff5T2_100(x, ud1, ud2);
}





// ==================== LAB 5 - PROBLEM RZECZYWISTY (BELKA) ====================

// masę belki
matrix ff5R_masa(matrix x, matrix ud1, matrix ud2) {
    double l = x(0) / 1000.0;  // konwersja mm -> m
    double d = x(1) / 1000.0;  // konwersja mm -> m
    double rho = 8920.0;       // gęstość w kg/m³

    double V = M_PI * pow(d / 2.0, 2.0) * l;
    double masa = rho * V;

    return matrix(masa);
}

// ugięcie belki
matrix ff5R_ugiecie(matrix x, matrix ud1, matrix ud2) {
    double l = x(0) / 1000.0;  // konwersja mm -> m
    double d = x(1) / 1000.0;  // konwersja mm -> m
    double P = 2000.0;         // siła w N (2 kN)
    double E = 120e9;          // moduł Younga w Pa (120 GPa)

    double u = (64.0 * P * pow(l, 3)) / (3.0 * E * M_PI * pow(d, 4));
    u = u * 1000.0;  // Konwersja na mm

    return matrix(u);
}

// naprężenie w belce
matrix ff5R_naprezenie(matrix x, matrix ud1, matrix ud2) {
    double l = x(0) / 1000.0;  // konwersja mm -> m
    double d = x(1) / 1000.0;  // konwersja mm -> m
    double P = 2000.0;         // siła w N (2 kN)

    double sigma = (32.0 * P * l) / (M_PI * pow(d, 3));
    sigma = sigma / 1e6;  // Konwersja na MPa

    return matrix(sigma);
}

// Funkcja celu (do testowania)
matrix ff5R_base(matrix x, matrix ud1, matrix ud2) {
    double masa = m2d(ff5R_masa(x, ud1, ud2));
    double ugiecie = m2d(ff5R_ugiecie(x, ud1, ud2));

    matrix result(3, 1);
    result(0) = masa;
    result(1) = ugiecie;
    result(2) = m2d(ff5R_naprezenie(x, ud1, ud2));

    return result;
}

// Główna funkcja celu z karą
matrix ff5R(matrix x, matrix ud1, matrix ud2) {
    try {
        // Sprawdź zakres zmiennych
        double l = x(0);  // mm
        double d = x(1);  // mm


        double penalty_box = 0.0;
        if (l < 200.0) penalty_box += 1e10 * pow(200.0 - l, 2);
        if (l > 1000.0) penalty_box += 1e10 * pow(l - 1000.0, 2);
        if (d < 10.0) penalty_box += 1e10 * pow(10.0 - d, 2);
        if (d > 50.0) penalty_box += 1e10 * pow(d - 50.0, 2);


        double masa = m2d(ff5R_masa(x, ud1, ud2));
        double ugiecie = m2d(ff5R_ugiecie(x, ud1, ud2));
        double naprezenie = m2d(ff5R_naprezenie(x, ud1, ud2));


        double u_max = 2.5;      // mm
        double sigma_max = 300.0; // MPa

        double penalty_constraints = 0.0;


        if (ugiecie > u_max) {
            penalty_constraints += 1e8 * pow(ugiecie - u_max, 2);
        }


        if (naprezenie > sigma_max) {
            penalty_constraints += 1e8 * pow(naprezenie - sigma_max, 2);
        }


        double f1_norm = masa / 10.0;
        double f2_norm = ugiecie / 100.0;


        double w_val;
        if (wi >= 0 && wi < 101) {
            w_val = w[wi];
        } else {
            w_val = 0.5;
        }

        // Funkcja celu: kryterium ważone + kara
        double f = w_val * f1_norm + (1.0 - w_val) * f2_norm + penalty_box + penalty_constraints;

        return matrix(f);
    }
    catch (string ex_info) {
        throw ("matrix ff5R(...):\n" + ex_info);
    }
}

matrix f_line_powell(matrix h, matrix xk, matrix dk) {
    // h - krok wzdłuż kierunku
    // xk - punkt bazowy
    // dk - kierunek poszukiwań
    return ff5R(xk + m2d(h) * dk, matrix(), matrix());
}


matrix ff6T(matrix x, matrix ud1, matrix ud2) {
    double x1 = x(0);
    double x2 = x(1);

    return x1 * x1 + x2 * x2 - cos(2.5 * M_PI * x1) - cos(2.5 * M_PI * x2) + 2;
}


// Układ równań różniczkowych
matrix df6(double t, matrix Y, matrix ud1, matrix ud2) {
    double m1 = 1.0, m2 = 2.0;
    double k1 = 4.0, k2 = 6.0;
    double F = 5.0;
    double b1 = ud1(0); // Optymalizowane b1
    double b2 = ud1(1); // Optymalizowane b2

    matrix dY(4, 1);
    dY(0) = Y(1); // dx1/dt
    dY(1) = (-b1 * Y(1) - b2 * (Y(1) - Y(3)) - k1 * Y(0) - k2 * (Y(0) - Y(2))) / m1; // d^2x1/dt^2
    dY(2) = Y(3); // dx2/dt
    dY(3) = (F + b2 * (Y(1) - Y(3)) + k2 * (Y(0) - Y(2))) / m2; // d^2x2/dt^2
    return dY;
}

// Funkcja celu dla problemu rzeczywistego
matrix ff6R(matrix x, matrix ud1, matrix ud2) {
    // ud1 będzie przechowywać dane z polozenia.txt wczytane raz
    // ud2 - pusta

    // x(0) = b1, x(1) = b2
    double t0 = 0, dt = 0.1, tend = 100;
    matrix Y0(4, 1, 0.0); // startujemy z nieruchomych ciężarków

    // Symulacja
    matrix* Y_sim = solve_ode(df6, t0, dt, tend, Y0, x, matrix());

    // Obliczanie błędu (Suma kwadratów różnic między symulacją a danymi z ud1)
    double error = 0;
    int n = get_size(Y_sim[0])[0]; // liczba kroków czasowych

    for (int i = 0; i < n; ++i) {
        // Dane z pliku (wczytane do ud1 w lab6()):
        // ud1(i, 0) - czas, ud1(i, 1) - x1_exp, ud1(i, 2) - x2_exp
        double x1_sim = Y_sim[1](i, 0);
        double x2_sim = Y_sim[1](i, 2);

        double x1_exp = ud1(i, 0);
        double x2_exp = ud1(i, 1);

        error += pow(x1_sim - x1_exp, 2) + pow(x2_sim - x2_exp, 2);
    }

    delete[] Y_sim;
    return matrix(sqrt(error / n)); // Błąd średniokwadratowy
}


