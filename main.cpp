/*********************************************
Kod stanowi uzupełnienie materiałów do ćwiczeń
w ramach przedmiotu metody optymalizacji.
Kod udostępniony na licencji CC BY-SA 3.0
Autor: dr inż. Łukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 30.09.2025
*********************************************/
#include"opt_alg.h"
#include <math.h>
void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

// Stałe
const double g = 9.81;
const double a_coef = 0.98;
const double b_coef = 0.63;
matrix ff(matrix x, matrix ud1, matrix ud2) {
    double D_A = x(0, 0) * 1e-4;  // konwersja cm² na m²

    // Parametry symulacji
    double dt = 1.0;           // krok czasu [s]
    double t_end = 2000.0;     // czas symulacji [s]

    // Parametry zbiorników
    double P_A = 2.0;          // pole podstawy A [m²]
    double P_B = 1.0;          // pole podstawy B [m²]
    double T_A_in = 95.0;      // temperatura wlotu A [°C]
    double T_B_in = 20.0;      // temperatura wlotu B [°C]
    double F_B_in = 0.01;      // natężenie wlotu B [m³/s] = 10 l/s
    double D_B = 36.5665e-4;   // pole otworu B [m²]

    // Stany początkowe
    double V_A_curr = 5.0;     // V_A(0) = 5 m³
    double T_A_curr = 95.0;    // T_A(0) = 95°C
    double V_B_curr = 1.0;     // V_B(0) = 1 m³
    double T_B_curr = 20.0;    // T_B(0) = 20°C

    double T_B_max = T_B_curr; // Śledzenie maksymalnej temperatury

    // Symulacja
    for (double t = 0; t < t_end; t += dt) {
        // Wysokość w zbiornikach: h = V/P
        double h_A = V_A_curr / P_A;
        double h_B = V_B_curr / P_B;

        if (h_A < 0.01) h_A = 0.01;  // zabezpieczenie przed zerem
        if (h_B < 0.01) h_B = 0.01;

        // Natężenia przepływu: F = a*b*D*sqrt(2*g*h)
        double F_A_out = a_coef * b_coef * D_A * sqrt(2.0 * g * h_A);
        double F_B_out = a_coef * b_coef * D_B * sqrt(2.0 * g * h_B);

        // Równania różniczkowe dla objętości
        double dV_A_dt = -F_A_out;
        double dV_B_dt = F_A_out + F_B_in - F_B_out;

        // Równania różniczkowe dla temperatury
        // Dla zbiornika A: wpływa woda o temp T_A_in = 95°C
        double dT_A_dt = 0.0;  // Brak zewnętrznego dopływu do A, więc temperatura się nie zmienia
        // (lub można założyć że wpływa woda o tej samej temp co już jest)

        // Dla zbiornika B: wpływa z A (F_A_out, T_A_curr) i z zewnątrz (F_B_in, T_B_in)
        double dT_B_dt = (F_A_out / V_B_curr) * (T_A_curr - T_B_curr) +
                         (F_B_in / V_B_curr) * (T_B_in - T_B_curr);

        // Aktualizacja stanów (metoda Eulera)
        V_A_curr += dV_A_dt * dt;
        V_B_curr += dV_B_dt * dt;
        T_A_curr += dT_A_dt * dt;  // W tym przypadku = 0, więc T_A nie zmienia się
        T_B_curr += dT_B_dt * dt;

        // Zabezpieczenia
        if (V_A_curr < 0.01) V_A_curr = 0.01;
        if (V_B_curr < 0.01) V_B_curr = 0.01;

        // Śledź maksymalną temperaturę w zbiorniku B
        if (T_B_curr > T_B_max) {
            T_B_max = T_B_curr;
        }
    }

    // Funkcja celu: (T_B_max - 50)²
    matrix y(1, 1);
    y(0, 0) = (T_B_max - 50.0) * (T_B_max - 50.0);

    return y;
}

int main()
{
	try
	{

		lab1();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	return 0;
}

void lab0()
{


    //expansion

    double x0 = 0.0;
    double d = 1.0;
    double alpha = 2.0;
    int Nmaxx = 1000;

    double* result = expansion(ff1T, x0, d, alpha, Nmaxx   );

    cout << "Przedzial [a, b]:" << endl;
    cout << "a = " << result[0] << endl;
    cout << "b = " << result[1]<< endl;
    cout << "Liczba wywolan: " << solution::f_calls << endl;

    solution::clear_calls();

	//Funkcja testowa
	double epsilon = 1e-2;									// dokładność
	int Nmax = 10000;										// maksymalna liczba wywołań funkcji celu
	matrix lb(2, 1, -5), ub(2, 1, 5),						// dolne oraz górne ograniczenie
		a(2, 1);											// dokładne rozwiązanie optymalne
	solution opt;											// rozwiązanie optymalne znalezione przez algorytm
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);			// wywołanie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie liczników

	//Wahadlo
	Nmax = 1000;											// dokładność
	epsilon = 1e-2;											// maksymalna liczba wywołań funkcji celu
	lb = 0, ub = 5;											// dolne oraz górne ograniczenie
	double teta_opt = 1;									// maksymalne wychylenie wahadła
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);		// wywołanie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie liczników

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki początkowe
		MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment siły działający na wahadło oraz czas działania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);	// rozwiązujemy równanie różniczkowe
	ofstream Sout("symulacja_lab0.csv");					// definiujemy strumień do pliku .csv
	Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
	Sout.close();											// zamykamy strumień
	Y[0].~matrix();											// usuwamy z pamięci rozwiązanie RR
	Y[1].~matrix();
}


void lab1()
{
	//Funkcja testowa
	srand(time(NULL));
//	double epsilon = 1e-8;
//	double gamma = 1e-10;
//	double alpha = 1.5;
	double x0 = (rand() % 100) - 50, d = 2.0;
//	int Nmax = 1000;
	int a = 0, b = 100;
	ofstream Sout("symulacja_lab1.csv");

//	for (int i = 0; i < 100; ++i) {
//		cout << "Punkt startowy: " << x0 << endl;
//		cout << "Expansion" << endl;
//		solution opt1;
//		double* p = expansion(ff1T, x0, d, alpha, Nmax);
//		opt1.x = p[0]; opt1.y = p[1]; opt1.f_calls = p[2];
//		cout << "[" << p[0] << "," << p[1] << "]" << endl << "f_calls = " << p[2] << endl << endl;
//		Sout << x0 << ";" << opt1.x << ";" << opt1.y << ";" << opt1.f_calls << ";";
//
//		cout << "Fibonacci" << endl;
//		solution opt3 = fib(ff1T, p[0], p[1], epsilon, gamma, Nmax);
//		cout << opt3 << endl << endl;
//		Sout << opt3.x << ";" << opt3.y << ";" << opt3.f_calls << ";";
//
//		cout << "Lagrange" << endl;
//		solution opt = lag(ff1T, p[0], p[1], epsilon, gamma, Nmax);
//		cout << opt << endl << endl;
//		Sout << opt.x << ";" << opt.y << ";" << opt.f_calls << ";\n";
//
//		++x0;
//		delete[] p;
//	}
//	Sout.close();





   // === TEST SYMULACJI dla D_A = 50 cm² ===
    cout << "=== WERYFIKACJA MODELU ===" << endl;
    matrix x_test(1, 1);
    x_test(0, 0) = 50.0;  // D_A = 50 cm²

    matrix result = ff(x_test, NAN, NAN);
    double error = result(0, 0);  // (T_B - 50)²

    // Obliczamy T_B z funkcji celu
    // (T_B - 50)² = error
    // T_B = 50 ± sqrt(error)
    // Ponieważ oczekujemy T_B = 62.5°C > 50°C, to:
    double T_B_final = 50.0 + sqrt(error);

    cout << "Test dla D_A = 50 cm2:" << endl;
    cout << "T_B koncowe = " << T_B_final << " C" << endl;
    cout << "Oczekiwane: ~62.5 C" << endl;
    cout << "Blad: " << abs(T_B_final - 62.5) << " C" << endl;
    cout << "Funkcja celu  = " << error << endl << endl;

    // === OPTYMALIZACJA ===
    cout << "=== OPTYMALIZACJA ===" << endl;
    double epsilon = 1e-6;
    double gamma = 1e-10;
    int Nmax = 1000;

    // Optymalizacja D_A w zakresie [1, 100] cm²
    solution sol_fib = fib(ff, 1.0, 100.0, epsilon, gamma, Nmax);
    solution sol_lag = lag(ff, 1.0, 100.0, epsilon, gamma, Nmax);

    cout << "Fibonacci - Optymalne D_A = " << sol_fib.x(0, 0) << " cm" << endl;
    cout << "Fibonacci - Funkcja celu = " << sol_fib.y(0, 0) << endl;
    cout << "Fibonacci - T_B = " << 50.0 + sqrt(sol_fib.y(0, 0)) << " C" << endl << endl;

    cout << "Lagrange - Optymalne D_A = " << sol_lag.x(0, 0) << " cm2" << endl;
    cout << "Lagrange - Funkcja celu = " << sol_lag.y(0, 0) << endl;
    cout << "Lagrange - T_B = " << 50.0 + sqrt(sol_lag.y(0, 0)) << " C" << endl;
}

void lab2()
{

}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}

