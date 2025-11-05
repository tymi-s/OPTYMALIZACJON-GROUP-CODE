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

int main()
{
	try
	{

		lab2();


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
	/*srand(time(NULL));

	double s = 0.15, alpha = 0.5, epsilon = 1e-6;
	double beta= 0.5;
	matrix s0(2, 1);
	s0(0, 0) = 0.15;
	s0(1, 0) = 0.15;
	//ofstream file("wyniki_Rosenbrock_3.csv");
	//file << "x1;x2;X1_min;X2_min;Y;Wywolania\n";

	for (int i=0;i<100;++i) {
		double x0data[2];
		x0data[0]=((rand()%200)/100.0)-1.0;
		x0data[1]=((rand()%200)/100.0)-1.0;
		matrix x0(2,x0data);
		cout << "Punkty startowe:( " << x0data[0] <<" , " << x0data[1] <<")"<< endl;
		solution::clear_calls();
		 solution opt1 = HJ(ff2T,x0, s, alpha, epsilon, 10000);
		//cout << "Hook-Jeeves" << endl;
		// cout <<"X = ("<<  opt1.x(0) <<","<< opt1.x(1) << ")" << endl;
		// cout <<"Y = " <<  opt1.y << endl;
		// cout <<"Wywolania funkcji: " << opt1.f_calls << ";" << endl;
		solution opt2 = Rosen(ff2T,x0, s0, alpha,beta, epsilon, 10000);
		cout << "Rosenbrock" << endl;
		cout <<"X = ("<<  opt2.x(0) <<","<< opt2.x(1) << ")" << endl;
		cout <<"Y = " <<  opt2.y << endl;
		cout <<"Wywolania funkcji: " << opt2.f_calls << ";" << endl;

		// file << x0data[0] << ";"
		//  << x0data[1] << ";"
		//  << opt2.x(0) << ";"
		//  << opt2.x(1) << ";"
		//  << opt2.y << ";"
		//  << opt2.f_calls << "\n";

	}*/

	solution::clear_calls();

	try {
		cout << "=== PROBLEM RZECZYWISTY - OPTYMALIZACJA REGULATORA ===" << endl;

		double s_hj = 1.0;          
		double alpha_hj = 0.5;       
		double epsilon = 1e-4;       
		int Nmax = 10000;           

		matrix s0_rosen(2, 1);
		s0_rosen(0) = 1.0;          
		s0_rosen(1) = 1.0;          
		double alpha_rosen = 2.0;   
		double beta_rosen = 0.5;    

		matrix x0(2, 1);
		x0(0) = 10.0;  
		x0(1) = 10.0;  

		cout << "\nPunkt startowy: k1 = " << x0(0) << ", k2 = " << x0(1) << endl;

		cout << "\n--- Weryfikacja dla k1=5, k2=5 ---" << endl;
		matrix x_test(2, 1);
		x_test(0) = 5.0;
		x_test(1) = 5.0;
		matrix Q_test = ff2R(x_test, NAN, NAN);
		cout << "Q(5, 5) = " << Q_test(0) << " (oczekiwane: ~775.229)" << endl;

		cout << "\n--- Metoda Hooke'a-Jeevesa ---" << endl;
		solution Xopt_HJ = HJ(ff2R, x0, s_hj, alpha_hj, epsilon, Nmax);
		cout << "Wynik: " << Xopt_HJ << endl;

		cout << "\n--- Metoda Rosenbrocka ---" << endl;
		solution Xopt_Rosen = Rosen(ff2R, x0, s0_rosen, alpha_rosen, beta_rosen, epsilon, Nmax);
		cout << "Wynik: " << Xopt_Rosen << endl;

		cout << "\n--- Symulacja z optymalnymi parametrami (HJ) ---" << endl;
		double k1_opt = Xopt_HJ.x(0);
		double k2_opt = Xopt_HJ.x(1);

		matrix k_opt(2, 1);
		k_opt(0) = k1_opt;
		k_opt(1) = k2_opt;

		matrix Y0(2, 1);
		Y0(0) = 0.0;
		Y0(1) = 0.0;

		matrix* S = solve_ode(df, 0.0, 0.1, 100.0, Y0, k_opt, NAN);

		ofstream file("symulacja.csv");
		file << "t;alpha;omega" << endl;

		int N = get_size(S[0])[0];
		for (int i = 0; i < N; i++) {
			file << S[0](i) << ";" << S[1](i, 0) << ";" << S[1](i, 1) << endl;
		}
		file.close();

		cout << "Wyniki symulacji zapisane do pliku symulacja.csv" << endl;

		delete[] S;

	}
	catch (string ex_info) {
		cout << ex_info << endl;
	}

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



