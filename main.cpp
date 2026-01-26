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

		lab6();


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
    ////////////////////////////////////////////////////////  PROBLEM RZECZYWISTY //////////////////////////////////////////////////
    cout << "\n=== PROBLEM RZECZYWISTY ===" << endl;

    // Wczytaj dane
    matrix X, Y;
    load_data(X, Y);

    // Punkt startowy
    matrix theta0(3, 1, 0.0);

    // Parametry
//    double steps[] = {0.01, 0.001, 0.0001};
//    int Nmax = 10000;
//    double epsilon = 1e-5;

    ofstream lab4Real("lab4Real.csv");
    lab4Real << "h;theta0;theta1;theta2;J(theta*);P(theta*);g_calls\n";

    for (int i = 0; i < 3; i++) {
        solution::clear_calls();

        solution cg = CG(ff4R, gradient4R, theta0, steps[i], epsilon, Nmax, X, Y);
        cg.fit_fun(ff4R, X, Y);

        double accuracy = calculate_accuracy(cg.x, X, Y);

        lab4Real << toExcel(steps[i]) << ";"
                 << toExcel(cg.x(0)) << ";"
                 << toExcel(cg.x(1)) << ";"
                 << toExcel(cg.x(2)) << ";"
                 << toExcel(cg.y(0)) << ";"
                 << toExcel(accuracy) << ";"
                 << cg.g_calls << "\n";

        cout << "h = " << steps[i] << ", Accuracy = " << accuracy << "%" << endl;
    }

    lab4Real.close();

        // ========== WYKRES DLA NAJLEPSZEGO PRZYPADKU ==========
    cout << "\n=== GENEROWANIE WYKRESU ===" << endl;

    // Weź najlepsze θ (z h=0.001)
    solution::clear_calls();
    solution best_cg = CG(ff4R, gradient4R, theta0, 0.001, epsilon, Nmax, X, Y);

    double theta0_val = best_cg.x(0);
    double theta1_val = best_cg.x(1);
    double theta2_val = best_cg.x(2);

    cout << "Najlepsze theta: [" << theta0_val << ", " << theta1_val << ", " << theta2_val << "]" << endl;

    ofstream wykres("lab4_wykres.csv");
    wykres << "x1;x2;typ\n";

    // Punkty przyjętych (zielone)
    for (int i = 0; i < 100; i++) {
        if (Y(0, i) == 1) {
            wykres << X(1, i) << ";" << X(2, i) << ";Przyjety\n";
        }
    }

    // Punkty odrzuconych (czerwone)
    for (int i = 0; i < 100; i++) {
        if (Y(0, i) == 0) {
            wykres << X(1, i) << ";" << X(2, i) << ";Odrzucony\n";
        }
    }

    // Granica klasyfikacji (niebieska linia)
    for (double x1 = 20; x1 <= 100; x1 += 1.0) {
        double x2 = -(theta0_val + theta1_val * x1) / theta2_val;
        wykres << x1 << ";" << x2 << ";Granica\n";
    }

    wykres.close();
    cout << "Wykres zapisany do: lab4_wykres.csv" << endl;
}

using FunctionPtr = matrix (*)(matrix, matrix, matrix);
void uruchom_optymalizacje(
	const std::vector<matrix>& start_points, // Przekazujemy punkty przez referencję (bez kopiowania)
	FunctionPtr func_celu,                   // Funkcja optymalizowana (np. ff5T3_1)
	FunctionPtr func_f1,                     // Funkcja składowa f1 (np. ff5T1_1)
	FunctionPtr func_f2,                     // Funkcja składowa f2 (np. ff5T2_1)
	double current_a,                        // Wartość parametru a (tylko do zapisu)
	std::ofstream& file                      // Otwarty plik do zapisu wyników
) {
	double epsilon = 1e-3;
	int Nmax = 1000;

	std::cout << "Przetwarzanie dla a = " << current_a << "..." << std::endl;

	// Pętla po wagach 0..100
	for (int i = 0; i < 101; i++) {
		setW(i);
		solution::clear_calls();

		// Używamy i-tego punktu z wektora start_points
		// Przekazujemy odpowiednią funkcję celu (func_celu)
		solution xOpt = Powell(func_celu, start_points[i], epsilon, Nmax, NAN, NAN);

		// Pobranie wyników
		double x1_res = xOpt.x(0);
		double x2_res = xOpt.x(1);


		double val_f1 = func_f1(xOpt.x, NAN, NAN)(0);
		double val_f2 = func_f2(xOpt.x, NAN, NAN)(0);
		int calls = solution::f_calls;

		// Zapis do wspólnego pliku
		// Format: a; w; x1_opt; x2_opt; f1; f2; calls
		file << x1_res << ";" << x2_res << ";"
			 << val_f1 << ";" << val_f2 << ";"
			 << calls << "\n";
	}
}
void lab5()
{


	//notatki
    // optymalna belka wyjdzie krótka
    // 101 optymalizacji - dla problemu rzeczywistego  i ten wykres to funckja f1 (oś x) względem funkcji f2 (oś y)
    // w przypadku testowym szukamy minimum f1 i f2
    // a w przypadku rzeczywistym f1 to jest masa a f2 to ugięcie
    // a = 100 sprawia że funkcja f1 jest o wiele ważniejsza dla tego dla w =0 mamy rozwiązanie f1 a dla w >0 mamy rozwiązanie f2
    // żeby nie dostać punktów tylko grupę punktów to trzeba f1 pomnożyć razy 1000 żeby obie funkcje przyjmowały wartości mniej więcej tego samego rzędu

	//==================================================================================== TEST ALGORYTMU  ==============================================================================================
	srand(time(NULL));
	double epsilon = 1e-3;
	int Nmax = 1000;
	double max_values[2]{ 10.0, 10.0 };
	double min_values[2]{ -10.0, -10.0 };

	makeW();

//	for (int i = 0; i < 101; i++)
//	{
//		matrix test = matrix(2, new double[2]{
//			min_values[0] + static_cast<double>(rand()) / RAND_MAX * (max_values[0] - min_values[0]),
//				min_values[1] + static_cast<double>(rand()) / RAND_MAX * (max_values[1] - min_values[1])
//			});
//		setW(i);
//		solution::clear_calls();
//		solution xOpt = Powell(ff5T3_1, test, epsilon, Nmax, NAN, NAN);
//		cout << xOpt << endl;

	//==================================================================================== OPTYMALIZACJA ==============================================================================================
//	srand(time(NULL));
//    makeW();
//    std::cout << "Generowanie stalej puli 101 punktow startowych..." << std::endl;
//    std::vector<matrix> starting_points;
//    starting_points.reserve(101);
//    ofstream file_start("x_startowe.csv");
//    file_start << "x_1;x_2\n";
//    for (int i = 0; i < 101; i++) {
//       double r1 = static_cast<double>(rand()) / RAND_MAX;
//       double r2 = static_cast<double>(rand()) / RAND_MAX;
//
//       double x1 = min_values[0] + r1 * (max_values[0] - min_values[0]);
//       double x2 = min_values[1] + r2 * (max_values[1] - min_values[1]);
//
//
//       starting_points.push_back(matrix(2, new double[2]{ x1, x2 }));
//
//
//       file_start <<  x1 << ";" << x2 << "\n";
//    }
//    file_start.close();
//
//
//    ofstream file_results("wyniki_zadanie_5a.csv");
//    if (!file_results.is_open()) {
//       std::cerr << "Nie udalo sie otworzyc pliku do zapisu!" << std::endl;
//       exit(0);
//    }
//
//    file_results << "x1_opt;x2_opt;f1_val;f2_val;calls\n";
//
//    uruchom_optymalizacje(starting_points, ff5T3_1, ff5T1_1, ff5T2_1, 1.0, file_results);
//
//    uruchom_optymalizacje(starting_points, ff5T3_10, ff5T1_10, ff5T2_10, 10.0, file_results);
//
//    uruchom_optymalizacje(starting_points, ff5T3_100, ff5T1_100, ff5T2_100, 100.0, file_results);
//
//
//    file_results.close();
//    std::cout << "Zakonczono. Wszystkie wyniki zapisane w 'wyniki_zadanie_5a.csv'." << std::endl;

	//==================================================================================== PROBLEM RZECZYWISTY  ==============================================================================================
	// ============== PROBLEM TESTOWY  ==============
	//dla a=1
	/*
	cout << "\n=== PROBLEM TESTOWY (a=1) ===" << endl;
	ofstream test1("lab5_test_a1.csv");
	test1 << "w;x1;x2;f1;f2;f_combined;f_calls\n";

	for (int i = 0; i < 101; i++)
	{
		matrix x0(2, 1);
		x0(0) = -10.0 + static_cast<double>(rand()) / RAND_MAX * 20.0;  // [-10, 10]
		x0(1) = -10.0 + static_cast<double>(rand()) / RAND_MAX * 20.0;  // [-10, 10]

		setW(i);
		solution::clear_calls();
		solution xOpt = Powell(ff5T3_1, x0, epsilon, Nmax, matrix(), matrix());


		double f1 = m2d(ff5T1_1(xOpt.x, matrix(), matrix()));
		double f2 = m2d(ff5T2_1(xOpt.x, matrix(), matrix()));
		double w_current = static_cast<double>(i) * 0.01;

		test1 << toExcel(w_current) << ";"
		      << toExcel(xOpt.x(0)) << ";" << toExcel(xOpt.x(1)) << ";"
		      << toExcel(f1) << ";" << toExcel(f2) << ";"
		      << toExcel(m2d(xOpt.y)) << ";" << xOpt.f_calls << "\n";

		cout << "w=" << w_current << " -> f1=" << f1 << ", f2=" << f2 << endl;
	}
	test1.close();
	cout << "Wyniki zapisane do: lab5_test_a1.csv\n" << endl;
	*/


	cout << "\n=== PROBLEM RZECZYWISTY - BELKA ===" << endl;

	// Test weryfikacyjny z instrukcji
	cout << "\n--- TEST WERYFIKACYJNY ---" << endl;
	matrix x_test(2, 1);
	x_test(0) = 500.0;  // l = 500 mm
	x_test(1) = 25.0;   // d = 25 mm

	matrix test_result = ff5R_base(x_test, matrix(), matrix());
	cout << "Dla l=500mm, d=25mm:" << endl;
	cout << "  Masa = " << test_result(0) << " kg (powinno ~2.19)" << endl;
	cout << "  Ugiecie = " << test_result(1) << " mm (powinno ~36.22)" << endl;
	cout << "  Naprezenie = " << test_result(2) << " MPa (powinno ~651.9)" << endl;

	// Optymalizacja dla różnych wag
	cout << "\n--- OPTYMALIZACJA ---" << endl;
	ofstream real("lab5_real.csv");
	real << "w;l0;d0;l_star;d_star;masa_kg;ugiecie_mm;naprezenie_MPa;f_combined;f_calls\n";

	for (int i = 0; i < 101; i++)
	{
		// Losowy punkt startowy w dozwolonym zakresie
		matrix x0(2, 1);
		x0(0) = 200.0 + static_cast<double>(rand()) / RAND_MAX * (1000.0 - 200.0);  // l ∈ [200, 1000]
		x0(1) = 10.0 + static_cast<double>(rand()) / RAND_MAX * (50.0 - 10.0);      // d ∈ [10, 50]

		setW(i);  // KLUCZOWE: ustaw wagę przed optymalizacją
		solution::clear_calls();

		// Wywołaj Powell z PUSTYMI macierzami (konstrukcja matrix() tworzy pustą macierz)
		solution xOpt = Powell(ff5R, x0, epsilon, Nmax, matrix(), matrix());

		// Oblicz rzeczywiste wartości dla znalezionego rozwiązania
		matrix result = ff5R_base(xOpt.x, matrix(), matrix());
		double masa = result(0);
		double ugiecie = result(1);
		double naprezenie = result(2);

		double w_current = static_cast<double>(i) * 0.01;

		real << toExcel(w_current) << ";"
			 << toExcel(x0(0)) << ";" << toExcel(x0(1)) << ";"
			 << toExcel(xOpt.x(0)) << ";" << toExcel(xOpt.x(1)) << ";"
			 << toExcel(masa) << ";" << toExcel(ugiecie) << ";" << toExcel(naprezenie) << ";"
			 << toExcel(m2d(xOpt.y)) << ";" << xOpt.f_calls << "\n";

		cout << "w=" << w_current << " -> l=" << xOpt.x(0) << "mm, d=" << xOpt.x(1)
		     << "mm, m=" << masa << "kg, u=" << ugiecie << "mm, sigma=" << naprezenie << "MPa";

		// Sprawdź ograniczenia
		bool violated = false;
		if (ugiecie > 2.5) {
			cout << " [UWAGA: Ugiecie > 2.5mm!]";
			violated = true;
		}
		if (naprezenie > 300.0) {
			cout << " [UWAGA: Naprezenie > 300MPa!]";
			violated = true;
		}
		if (!violated && ugiecie <= 2.5 && naprezenie <= 300.0) {
			cout << " [OK]";
		}
		cout << endl;
	}
	real.close();

	cout << "\n=== PODSUMOWANIE ===" << endl;
	cout << "Wyniki zapisane do: lab5_real.csv" << endl;
	cout << "\nAby zwizualizowac front Pareto:" << endl;
	cout << "1. Otworz lab5_real.csv w Excelu" << endl;
	cout << "2. Stworz wykres rozrzutu (scatter plot) z:" << endl;
	cout << "   - Os X: masa_kg" << endl;
	cout << "   - Os Y: ugiecie_mm" << endl;
	cout << "3. Wykres pokaze front Pareto - kompromis miedzy masa a ugieciem" << endl;
	cout << "\nInterpretacja wynikww:" << endl;
	cout << "- w=0.0: minimalizacja TYLKO ugiecia (masa moze byc duze)" << endl;
	cout << "- w=1.0: minimalizacja TYLKO masy (ugiecie może byc duze)" << endl;
	cout << "- w=0.5: kompromis 50/50 miedzy masa a ugieciem" << endl;
	cout << "- Rozwi zania optymalne w sensie Pareto leza na krzywej" << endl;
}






void lab6()
{
	/*
	// --- STARY KOD TESTOWY (ZAKOMENTOWANY) ---
	srand(time(NULL));
	double epsilon = 1e-3;
	int Nmax = 10000;
	int mi = 5;
	int lambd = 10;
	double sigma[] = { 0.01, 0.1, 1.0, 10.0, 100.0 };
	matrix lb(2, std::unique_ptr<double[]>(new double[2]{ -5.0, -5.0 }).get()),
		   ub(2, std::unique_ptr<double[]>(new double[2]{ 5.0, 5.0 }).get());

	for (int i = 0; i < 100; i++)
	{
		solution::clear_calls();
		solution xOpt = EA(ff6T, 2, lb, ub, mi, lambd, sigma[4], epsilon, Nmax);
		cout << xOpt << endl;
	}
	*/

	// =================================================================
	// === PROBLEM RZECZYWISTY - OPTYMALIZACJA b1, b2 ===
	// =================================================================

	 // cout << "=== LAB 6: PROBLEM RZECZYWISTY (IDENTYFIKACJA b1, b2) ===" << endl;
  //
  //   // 1. Wczytywanie danych pomiarowych (POPRAWIONE)
  //   ifstream file("polozenia.txt");
  //   if (!file.is_open()) {
  //       cerr << "BLAD: Brak pliku polozenia.txt!" << endl;
  //       return;
  //   }
  //
  //   matrix data_exp(1001, 2);
  //   string line;
  //   int row = 0;
  //   while (getline(file, line) && row < 1001) {
  //       // Zamieniamy przecinki na kropki i sredniki na spacje
  //       replace(line.begin(), line.end(), ',', '.');
  //       replace(line.begin(), line.end(), ';', ' ');
  //
  //       stringstream ss(line);
  //       double val1, val2;
  //       // Wczytujemy TYLKO dwie wartosci, bo tak wyglada Twoj plik
  //       if (ss >> val1 >> val2) {
  //           data_exp(row, 0) = val1; // x1_exp
  //           data_exp(row, 1) = val2; // x2_exp
  //           row++;
  //       }
  //   }
  //   file.close();
  //   cout << "Wczytano " << row << " wierszy danych eksperymentalnych." << endl;

    // 2. Ustawienia EA
	int N = 2;
	// Zakładam metodę wypełniania lub dostęp indeksowy
	matrix lb(N, 1);
	matrix ub(N, 1);

	lb(0, 0) = -5.0;
	lb(1, 0) = -5.0;
	ub(0, 0) = 5.0;
	ub(1, 0) = 5.0;
    // matrix lb(2, 1, 0.1);
    // matrix ub(2, 1, 3.0);
    int mi = 20;
    int lambda = 60;
    //double sigma = 0.5;
    double epsilon = 1e-6;
    int Nmax = 15000;

	double sigma_values[] = {0.01, 0.1, 1.0, 10.0, 100.0};
	int num_sigma = 5;
	int num_runs = 100;

	// Generator liczb losowych
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dis(-5.0, 5.0);

	ofstream outfile("wyniki_optymalizacji_3.txt");
	outfile << fixed << setprecision(10);
	outfile << "Sigma\tRun\tx1_start\tx2_start\tx1_opt\tx2_opt\tf_opt\tIteracje\n";

	// Pętla po wartościach sigma
	for (int s = 0; s < num_sigma; s++) {
		double sigma = sigma_values[s];
		matrix sigma0(N, 1);
		sigma0(0) = sigma;
		sigma0(1) = sigma;

		cout << "Testowanie sigma = " << sigma << endl;

		// 100 optymalizacji dla danej wartości sigma
		for (int run = 0; run < num_runs; run++) {
			// Losowy punkt startowy
			matrix x0(N, 1);
			x0(0) = dis(gen);
			x0(1) = dis(gen);

			// Zapisz punkt startowy
			double x1_start = x0(0);
			double x2_start = x0(1);

			// Wykonaj optymalizację
			solution result = EA(ff6T, N, lb, ub, mi, lambda, sigma0, epsilon, Nmax);
             int iteracje = solution::f_calls/ lambda;
			// Zapisz wyniki do pliku
			outfile << sigma << "\t"
					<< run + 1 << "\t"
					<< x1_start << "\t"
					<< x2_start << "\t"
					<< result.x(0) << "\t"
					<< result.x(1) << "\t"
					<< result.y(0) << "\t"
					<< iteracje << "\n";

			// Informacja o postępie co 10 uruchomień
			if ((run + 1) % 10 == 0) {
				cout << "  Ukończono " << run + 1 << " / " << num_runs << " uruchomień" << endl;
			}
		}

		cout << "Zakończono dla sigma = " << sigma << "\n" << endl;
	}

	outfile.close();
	cout << "Wyniki zapisane do pliku: wyniki_optymalizacji.txt" << endl;

	// Opcjonalnie: oblicz i wyświetl statystyki
	cout << "\nOptymalizacja zakończona!" << endl;

    // // 3. Start optymalizacji
    // solution::clear_calls();
    // cout << "Optymalizacja trwa (to moze potrwac ok. 30-60 sekund)..." << endl;
    // solution opt = EA(ff6R, N, lb, ub, mi, lambda, sigma, epsilon, Nmax, data_exp);
    //
    // // 4. Wyniki do Tabeli 3
    // cout << "\n--- WYNIKI DO TABELI 3 ---" << endl;
    // cout << "b1_opt = " << opt.x(0) << " Ns/m" << endl;
    // cout << "b2_opt = " << opt.x(1) << " Ns/m" << endl;
    // cout << "Blad (y*) = " << opt.y(0) << endl;
    // cout << "Liczba wywolan = " << solution::f_calls << endl;
    //
    // // 5. Zapis danych do wykresu (do Excela)
    // matrix Y0(4, 1, 0.0);
    // matrix* Y = solve_ode(df6, 0.0, 0.1, 100.0, Y0, opt.x);
    // ofstream Sout("wynik_symulacji_lab6.csv");
    // Sout << "t;x1_sim;x2_sim;x1_exp;x2_exp\n";
    // for (int i = 0; i < 1001; i++) {
    //     Sout << toExcel(Y[0](i, 0)) << ";"
    //          << toExcel(Y[1](i, 0)) << ";"
    //          << toExcel(Y[1](i, 2)) << ";"
    //          << toExcel(data_exp(i, 0)) << ";"
    //          << toExcel(data_exp(i, 1)) << "\n";
    // }
    // Sout.close();
    // cout << "\nZapisano 'wynik_symulacji_lab6.csv'. Otworz go w Excelu." << endl;
    // delete[] Y;
}


