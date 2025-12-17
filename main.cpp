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
#include "user_funs.h"

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

string toExcel(double value) {
	string s = to_string(value);
	replace(s.begin(), s.end(), '.', ',');
	return s;
}
                                                                                    
int main()
{


	try
	{

		lab4();


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

    //solution result = expansion(ff1T, x0, d, alpha, Nmaxx   );

    cout << "Przedzial [a, b]:" << endl;
    //cout << "a = " << result.x(0) << endl;
    //cout << "b = " << result.x(1) << endl;
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


void lab1(){
	//Funkcja testowa
	srand(time(NULL));
	double epsilon = 1e-8;
	double gamma = 1e-10;
	double alpha = 1.5;
	double x0 = (rand() % 100) - 50, d = 2.0;
	int Nmax = 1000;
	int a = 0, b = 100;
	ofstream Sout("symulacja_lab1.csv");

	for (int i = 0; i < 100; ++i) {
		cout << "Punkt startowy: " << x0 << endl;
		cout << "Expansion" << endl;
		solution opt1;
		double* p = expansion(ff1T, x0, d, alpha, Nmax);
		opt1.x = p[0]; opt1.y = p[1]; opt1.f_calls = p[2];
		cout << "[" << p[0] << "," << p[1] << "]" << endl << "f_calls = " << p[2] << endl << endl;
		Sout << x0 << ";" << opt1.x << ";" << opt1.y << ";" << opt1.f_calls << ";";

		cout << "Fibonacci" << endl;
		solution opt3 = fib(ff1T, p[0], p[1], epsilon, gamma, Nmax);
		cout << opt3 << endl << endl;
		Sout << opt3.x << ";" << opt3.y << ";" << opt3.f_calls << ";";

		cout << "Lagrange" << endl;
		solution opt = lag(ff1T, p[0], p[1], epsilon, gamma, Nmax);
		cout << opt << endl << endl;
		Sout << opt.x << ";" << opt.y << ";" << opt.f_calls << ";\n";

		++x0;
		delete[] p;
	}
	Sout.close();
}

void lab2()
{
	srand(time(NULL));

	double s = 0.08, alpha = 0.5, epsilon = 1e-6;
	double beta = 0.5;
	matrix s0(2, 1);
	s0(0, 0) = 0.08;
	s0(1, 0) = 0.08;

	// Stały punkt startowy dla wizualizacji ścieżki
	double x0data[2];
	x0data[0] = 0.5;  // Możesz zmienić na inne wartości
	x0data[1] = 0.5;
	matrix x0(2, x0data);

	cout << "Punkt startowy: (" << x0data[0] << ", " << x0data[1] << ")" << endl << endl;

	// Hook-Jeeves
	cout << "=== Hook-Jeeves ===" << endl;
	ofstream file_HJ("sciezka_HJ.csv");
	file_HJ << "X1;X2;Y\n";
	solution::clear_calls();
	solution opt1 = HJ(ff2T, x0, s, alpha, epsilon, 10000);
	file_HJ.close();

	cout << "Wynik koncowy:" << endl;
	cout << "  X = (" << opt1.x(0) << ", " << opt1.x(1) << ")" << endl;
	cout << "  Y = " << opt1.y << endl;
	cout << "  Wywolania funkcji: " << opt1.f_calls << endl << endl;

	// Rosenbrock
	cout << "=== Rosenbrock ===" << endl;
	ofstream file_Rosen("sciezka_Rosen.csv");
	file_Rosen << "X1;X2;Y\n";
	solution::clear_calls();
	solution opt2 = Rosen(ff2T, x0, s0, alpha, beta, epsilon, 10000);
	file_Rosen.close();

	cout << "Wynik koncowy:" << endl;
	cout << "  X = (" << opt2.x(0) << ", " << opt2.x(1) << ")" << endl;
	cout << "  Y = " << opt2.y << endl;
	cout << "  Wywolania funkcji: " << opt2.f_calls << endl << endl;

	cout << "Sciezki zapisane do plikow: sciezka_HJ.csv i sciezka_Rosen.csv" << endl;
}

void lab3()
{
    ////////////////////////////////////////////// PROBLEM RZECZYWISTY ///////////////////////////////////////

    /// parametry fizyczne:
    double masa  = 0.6;
    double promien = 0.12;
    double g = 9.81;
    double wspolczynnik_oporu = 0.47;
    double gestosc_powietrza = 1.2;

    /// parametry optymalizacji:
    double epsilon = 1e-4;
    int Nmax = 50000;

    double c = 2000.0;
    double dc = 10.0;


    matrix x0(2, 1);
    x0(0) = 5.0;
    x0(1) = 10.0;
 //   x0(0) = 0.0;
 //   x0(1) = 45.0;
    /// Przygotowanie parametrdla funkcji celu
    matrix ud1(5, 1);  // parametry fizyczne (bez wspÃ³Å‚czynnika kary)
    ud1(0) = masa;
    ud1(1) = promien;
    ud1(2) = wspolczynnik_oporu;
    ud1(3) = gestosc_powietrza;
    ud1(4) = g;

    matrix ud2(1, 1, 0.0);  // Pusta macierz


        // TEST WERYFIKACYJNY Z INSTRUKCJI
    cout << "=== TEST WERYFIKACYJNY MODELU ===" << endl;

    matrix x_test(2, 1);
    x_test(0) = 5.0;   // v0x
    x_test(1) = 10.0;  // omega

    matrix ud_test(5, 1);
    ud_test(0) = 0.6;    // masa
    ud_test(1) = 0.12;   // promien
    ud_test(2) = 0.47;   // C
    ud_test(3) = 1.2;    // rho
    ud_test(4) = 9.81;   // g

    matrix result_test = ff3R_base(x_test, ud_test, matrix(1,1,0.0));

    cout << "Dla v0x=5, omega=10:" << endl;
    cout << "  x_end = " << result_test(0) << " m (powinno ~41.41)" << endl;
    cout << "  x przy y=50 = " << result_test(1) << " m (powinno ~21.61)" << endl;
    cout << endl;

    cout << "=== OPTYMALIZACJA PROBLEMU RZECZYWISTEGO ===" << endl;
    cout << "Punkt startowy: v0x = " << x0(0) << " m/s, omega = " << x0(1) << " rad/s" << endl << endl;

    solution::clear_calls();

    /// OPTYMALIZACJA
    solution opt = pen_rzeczywisty(ff3R, x0, c, dc, epsilon, Nmax, ud1, ud2);

    cout << "=== WYNIKI OPTYMALIZACJI ===" << endl;
    cout << opt << endl;
    cout << "v0x_opt = " << opt.x(0) << " m/s" << endl;
    cout << "omega_opt = " << opt.x(1) << " rad/s" << endl;
    cout << "Wartosc funkcji celu (z kara): " << opt.y << endl << endl;
    cout << "Liczba wywolan funkcji celu: " << solution::f_calls << endl;  // ← DODAJ TO
    cout << "Liczba wywolan z opt.f_calls: " << opt.f_calls << endl << endl;  // ← I TO

    /// SYMULACJA DLA OPTYMALNYCH PARAMETRÃ“W
    cout << "=== SYMULACJA TRAJEKTORII ===" << endl;

    matrix Y0_sim(4, 1);
    Y0_sim(0) = 0.0;
    Y0_sim(1) = 100.0;
    Y0_sim(2) = opt.x(0);  // v0x optymalne
    Y0_sim(3) = 0.0;

    matrix params_sim(7, 1);
    params_sim(0) = opt.x(0);
    params_sim(1) = opt.x(1);
    params_sim(2) = masa;
    params_sim(3) = promien;
    params_sim(4) = wspolczynnik_oporu;
    params_sim(5) = gestosc_powietrza;
    params_sim(6) = g;

    matrix* Y = solve_ode(df3R, 0.0, 0.01, 7.0, Y0_sim, params_sim, ud2);

    /// Zapis do pliku CSV
    ofstream Sout("zadanie3_SYMULACJA_backspin.csv");
    Sout << "t;x;y;vx;vy\n";

    int n = get_size(Y[0])[0];
    double x_end_actual = 0.0;
    double x_at_y50_actual = 0.0;

    for (int i = 0; i < n; i++) {
        double t = Y[0](i, 0);
        double x = Y[1](i, 0);
        double y = Y[1](i, 1);
        double vx = Y[1](i, 2);
        double vy = Y[1](i, 3);

        Sout << t << ";" << x << ";" << y << ";" << vx << ";" << vy << "\n";


        if (i > 0 && Y[1](i-1, 1) >= 50.0 && y < 50.0) {
            x_at_y50_actual = x;
        }


        if (y <= 0.0 && i > 0) {
            x_end_actual = x;
            break;
        }
    }

    Sout.close();

    cout << "\n=== SZCZEGOLOWE WYNIKI SYMULACJI ===" << endl;
    cout << "Rzeczywisty x_end = " << x_end_actual << " m" << endl;
    cout << "Pozycja x przy y=50m: " << x_at_y50_actual << " m" << endl;

    // Sprawdzenie ograniczenia
    bool trafiono_kosz = (x_at_y50_actual >= 3.0 && x_at_y50_actual <= 7.0);

    cout << "\n--- SPRAWDZENIE OGRANICZENIA ---" << endl;
    cout << "Wymagany przedzial przy y=50m: [3.0, 7.0] m" << endl;
    cout << "Aktualna pozycja przy y=50m: " << x_at_y50_actual << " m" << endl;

    if (trafiono_kosz) {
        cout << " TRAFIONO DO KOSZA! Ograniczenie spelnione." << endl;
    } else {
        cout << " NIE TRAFIONO DO KOSZA! Ograniczenie NIE spelnione." << endl;
        if (x_at_y50_actual < 3.0) {
            cout << "  Pilka przeleciala " << (3.0 - x_at_y50_actual) << " m ZA BLISKO" << endl;
        } else {
            cout << "  Pilka przeleciala " << (x_at_y50_actual - 7.0) << " m ZA DALEKO" << endl;
        }
    }



    delete[] Y;
}

void lab4()
{
	double steps[] = { 0.05, 0.25 }; 
	int Nmax = 1000;
	double epsilon = 1e-3; 

	srand(time(NULL));
	
	ofstream lab4SD("lab4Static.csv");
	if (!lab4SD.good()) {
		cout << "Cant open file";
		return;
	}

	ofstream lab4Golden("lab4Golden.csv");
	if (!lab4Golden.good()) {
		cout << "Cant open file";
		return;
	}

	lab4SD << "SD;X1;" << "X2;" << "x1_res;" << "x2_res;" << "y;" << "f_calls;" << "g_calls;"
		<< "CG;" << "x1_res;" << "x2_res;" << "y;" << "f_calls;" << "g_calls;"
		<< "Newton;" << "x1_res;" << "x2_res;" << "y;" << "f_calls;" << "g_calls\n";
	lab4Golden << "SD;X1;" << "X2;" << "x1_res;" << "x2_res;" << "y;" << "f_calls;" << "g_calls;"
		<< "CG;" << "x1_res;" << "x2_res;" << "y;" << "f_calls;" << "g_calls;"
		<< "Newton;" << "x1_res;" << "x2_res;" << "y;" << "f_calls;" << "g_calls\n";

	for (int step = 0; step < 2; ++step) {
		//Static steps
		for (int i = 0; i < 100; ++i) {
			solution::clear_calls();
			matrix x0(2, 1);
			x0(0) = ((double)rand() / RAND_MAX) * 4.0 - 2.0;
			x0(1) = ((double)rand() / RAND_MAX) * 4.0 - 2.0;

			//SD
			solution sd = SD(ff4T, gradient, x0, steps[step], epsilon, Nmax);
			sd.fit_fun(ff4T, NAN, NAN);
			lab4SD << ';' << toExcel(x0(0)) << ";" << toExcel(x0(1)) << ";" << toExcel(sd.x(0)) << ";" << toExcel(sd.x(1)) << ";"
				<< sd.y << sd.f_calls << ";" << sd.g_calls << ";";


			//CG
			solution::clear_calls();
			solution cg = CG(ff4T, gradient, x0, steps[step], epsilon, Nmax);
			cg.fit_fun(ff4T, NAN, NAN);
			lab4SD << ";" << toExcel(cg.x(0)) << ";" << toExcel(cg.x(1)) << ";"
				<< cg.y << cg.f_calls << ";" << cg.g_calls << ";";

			//Newton
			solution::clear_calls();
			solution nt = Newton(ff4T, gradient, hessian, x0, steps[step], epsilon, Nmax);
			nt.fit_fun(ff4T, NAN, NAN);
			lab4SD << ";" << toExcel(nt.x(0)) << ";" << toExcel(nt.x(1)) << ";"
				<< nt.y << nt.f_calls << ";" << nt.g_calls << "\n";

			//SD Golden
			solution::clear_calls();
			solution sdGolden = SD(ff4T, gradient, x0, 0.0, epsilon, Nmax);
			sdGolden.fit_fun(ff4T, NAN, NAN);
			lab4Golden << ";" << toExcel(x0(0)) << ";" << toExcel(x0(1)) << ";" << toExcel(sdGolden.x(0)) << ";" << toExcel(sdGolden.x(1)) << ";"
				<< sdGolden.y << sdGolden.f_calls << ";" << sdGolden.g_calls << ";";

			//CG Golden
			solution::clear_calls();
			solution cgGolden = CG(ff4T, gradient, x0, 0.0, epsilon, Nmax);
			cgGolden.fit_fun(ff4T, NAN, NAN);
			lab4Golden << ";" << toExcel(cgGolden.x(0)) << ";" << toExcel(cgGolden.x(1)) << ";"
				<< cgGolden.y << cgGolden.f_calls << ";" << cgGolden.g_calls << ";";

			//Newtwon Godlen
			solution::clear_calls();
			solution ntGolden = Newton(ff4T, gradient, hessian, x0, 0.0, epsilon, Nmax);
			ntGolden.fit_fun(ff4T, NAN, NAN);
			lab4Golden << ";" << toExcel(ntGolden.x(0)) << ";" << toExcel(ntGolden.x(1)) << ";"
				<< ntGolden.y << ntGolden.f_calls << ";" << ntGolden.g_calls << "\n";
		}
		cout << "\nnext step\n";
		lab4SD << "\nNext Step\n";
		lab4Golden << "\nNext step\n";
	}
	lab4SD.close();
	lab4Golden.close();
}

void lab5()
{

}

void lab6()
{

}




