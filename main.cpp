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

    solution result = expansion(ff1T, x0, d, alpha, Nmaxx   );

    cout << "Przedzial [a, b]:" << endl;
    cout << "a = " << result.x(0) << endl;
    cout << "b = " << result.x(1) << endl;
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

