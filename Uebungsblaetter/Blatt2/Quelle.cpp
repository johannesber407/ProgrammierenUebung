#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>
#include <sstream>

using namespace std;

void a_4a() {
	int N;
	float p = 1;
	cout << "Aufgabe 4a:" << endl;
	cout << "Geben sie eine ganze Zahl N ein:";
	cin >> N;
	for (int i = 0; i <= N; i++) {
		p *= ((2 * i) + 1);
	}
	cout << "Das Ergebnis ist: " << p;
}
void a_4b() {
	int N;
	cout << "Aufgabe 4b:" << endl;
	cout << "Geben sie eine ganze Zahl N ein:";
	cin >> N;
	double s=0.0;
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= i; j++) {
			s +=( 1.0 / (i * j));
		}
	}
	cout << "Das Ergebnis ist: " << s;
}

void aufgabe4() {
	char t;
	cout << "welche Teilaufgabe? ";
	cin >> t;
	switch (t) {
	case('a'):a_4a();
		break;
	case('b'):a_4b();
		break;
	}
	
}
void aufgabe5() {
	cout << "Aufgabe 5:" << endl;
	int N;
	cout << "Geben sie eine ganze Zahl N ein:";
	cin >> N;
	double s=1;
	for (int i = N; i >= 1; i--) {
		s = 1.0 +  (i / (2 * i + 1.0))*s;
	}
	cout << "Das Ergebnis ist: p_n =" << setprecision(10) << s <<" und 2p_n = " <<2*s;

}
void aufgabe6() {
	cout << "Aufgabe 6:" << endl;
	cout << "Betrachte float:" << endl;
	float eps=1.0;
	while (0 != (1 - (1 + eps))) {
		eps = eps / 2;
	}
	cout << "Die Groessenornung ist:" << 2 * eps <<endl<<endl;
	cout << "Betrachte double:" << endl;
	double eps1 = 1.0;
	while (0 != (1 - (1 + eps1))) {
		eps1 = eps1 / 2;
	}
	cout << "Die Groessenornung ist:" << 2 * eps1 << endl;
}

/*
 * Programmieren fuer Physiker, SS 2023, Blatt 03, Aufg. 7
 *
 * kubische Gleichung loesen. Fragment.
 */


// reelle dritte Wurzel berechnen, auch fuer neg. Argumente
double cubicroot(double x)
{
	if (x > 0) return exp(log(x) / 3.0);
	if (x < 0) return -exp(log(-x) / 3.0);
	return 0.0;
}
double pi(int N) {
	double s = 1;
	for (int i = N; i >= 1; i--) {
		s = 1.0 + (i / (2 * i + 1.0)) * s;
	}
	return(2*s);
}

void aufgabe7() {
	double a, b, c, d, p, q, D, x1, x2, x3;
	double Pi = pi(10000);
	cout << "Aufgabe 7" << endl;
	cout << "Bitte Parameter a eingeben: ";
	cin >> a;
	cout << "Bitte Parameter b eingeben: ";
	cin >> b;
	cout << "Bitte Parameter c eingeben: ";
	cin >> c;
	cout << "Bitte Parameter d eingeben: ";
	cin >> d;

	p = (3 * a * c - b * b) / ((3 * a * a));
	q = ((2 * b * b * b) / (27 * a * a * a)) - ((b * c) / (3 * a * a)) + (d / a);
	D = pow((q / 2), 2) + pow((p / 3), 3);

	if (D > 0) {
		x1 = cubicroot(((-1) * q / 2) + sqrt(D)) + cubicroot(((-1) * q / 2) - sqrt(D)) - (b / (3 * a));
		cout << "Es gibt die einfache Nullstelle x1 = " << x1;
	}
	else if (D == 0) {
		if (q == 0) {
			x1 = -b / (3 * a);
			cout << "Es gibt die dreifache Nullstelle x1 = " << x1;
		}
		else {
			x1 = cubicroot(q / 2) - (b / (3 * a));
			x2 = cubicroot(-4 * q) - (b / (3 * a));
			cout << "Es gibt die doppelte Nullstelle x1 = " << x1 << " und die einfache Nullstelle x2 = " << x2;
		}
	}
	else {
		double h = acos((-q / 2) * sqrt(-27 / (p * p * p)));
		x1 = -sqrt(-4 * p / 3) * cos((h / 3) + (Pi / 3)) - b / (3 * a);
		x2 = -sqrt(-4 * p / 3) * cos((h / 3) - (Pi / 3)) - b / (3 * a);
		x3 = sqrt((-4 * p) / 3) * cos(h / 3) - b / (3 * a);
		cout << "Es gibt die einfache Nullstelle x1 = " << x1 << ", die einfachen Nullstelle x2 = " << "und die einfache Nullstelle x3 = " << x3;
	}
}
int binomi(int n, int k) {
	if (k == 0 || n == k) {
		return 1;
	}
	else {
		return (binomi(n - 1, k - 1) + binomi(n - 1, k));
	}
}
void aufgabe8() {
	//const int l = 30;
	cout << "Aufgabe 8: " << endl;
	for (int i = 0; i <= 10; i++) {
		//int output[l];
		for (int n = 1; n <= 10 - i; n++) {
			cout << "   ";
		}
		for (int j = 0; j <= i; j++) {
			int b = binomi(i, j);
			cout <<setw(5)<< b << " ";
			//output[10 - i + j]=b;
		}
		/*for (int k = l - 1; k >= 0; i--) {
			cout << setw(5) << output[k];
		}*/
		cout << endl;
	}
}
int* rekursion(int a, int n, int max) {
	if (a == 1) {
		int e[3] = { 1,n, max };
		return e;
	}
	n = n + 1;
	if (a % 2 == 0) {
		a = a / 2;
	}
	else {
		a = 3*a + 1;
	}
	if (a > max) {
		max = a;
	}
	int* test = rekursion(a, n, max);
	int res[3] = { test[0], test[1], test[2] };// keine ahnung, warum, aber anders hat es nicht funktioniert.
	return res;
}
void aufgabe10() {
	int a, a_init, n, n_max, a_max; 
	int* res1;
	a_max = 0;
	n_max = 0;
	for (int i = 1; i <= 100; i++) {
		
		res1 = rekursion(i, 0, i);
		if (res1[1] > n_max) {
			n_max = res1[1];
			a_max = i;
		}

	}
	res1 = rekursion(a_max, 0, a_max);
	cout << "Die längste Sequenz hat a_initial= " << a_max<< " Schritte n= " << res1[1] << " a_max= " << res1[2] << endl;
	cout << "Bitte eine natuerliche Zahl eingeben: ";
	cin >> a;
	a_init = a;
	n = 0;
	int* res;
	//int x, y, z;
	 res = rekursion(a, n, a);
	 cout << "a_initial= " << a_init << " Schritte n= " << res[1] << " a_max= " << res[2] << endl;
}

int aufgabe11() {
	vector<vector<float>> res;
	string filename("a11-weitsprung.txt");
	//float number;
	
	ifstream input_file(filename);
	if (!input_file.is_open()) {
		cerr << "Could not open the file - '"
			<< filename << "'" << endl;
		return EXIT_FAILURE;
	}

	string line;
	while (getline(input_file, line)) {
		vector<float> row; // Individual row vector
		istringstream iss(line);
		float value;

		while (iss >> value) {
			row.push_back(value);
		}
		res.push_back(row);
	}
	/*while (input_file >> number) {
		vector<float> row;
		row.push_back(number);
		cout << number << "; ";
		res[0].push_back(number);
	}*/
	cout << endl;
	input_file.close();
	vector <float> best;
	int champion = 0;
	float global_max= 0;
	for (int i = 0; i < res.size(); i++) {
		float max = 0;
		for (int j = 1; j < res[i].size(); j++) {
			if (res[i][j] > max) {
				max = res[i][j];
			}
		}
		if (global_max < max) {
			global_max = max;
			champion = i;
		}
		best.push_back(max);
	}
	
	for (int i = 0; i < res.size(); i++) {
		for (int j = 0; j < res[i].size(); j++) {
			cout << setw(6) << res[i][j] << "  "; 
		}
		cout << "Bestweite: " << best[i];
		if (i == champion) {
			cout << " Champion";
		}
		cout << endl;
	}
	return EXIT_SUCCESS;
}
int aufgabe13() {
	int test;
	string filename;
	string pivot;
	//einlesen
	cout << "Gauss - Algorithmus" << endl << "Welche Testdatei? ";
	cin >> test;
	if (test == 1) {
		filename = "a13-lgs1.txt";
	}
	else {
		filename = "a13-lgs2.txt";
	}
	vector<vector<float>> data;
	
	//float number;

	ifstream input_file(filename);
	if (!input_file.is_open()) {
		cerr << "Could not open the file - '"
			<< filename << "'" << endl;
		return EXIT_FAILURE;
	}

	string line;
	while (getline(input_file, line)) {
		vector<float> row; // Individual row vector
		istringstream iss(line);
		float value;

		while (iss >> value) {
			row.push_back(value);
		}
		data.push_back(row);
	}
	int dimension = data[0][0];
	vector<vector<float>> A;
	vector<float> b;
	for (int i = 1; i <= dimension; i++) {
		A.push_back(data[i]);
	}
	b = data[dimension + 1];
	cout << "Dimension = " << dimension << endl;

	for (int i = 0; i < A.size(); i++) {
		for (int j = 0; j < A[i].size(); j++) {
			cout << setw(6) << A[i][j] << "  ";
		}
		cout << setw(10) << "|";
		cout << setw(6) << b[i];
		cout << endl;
	}
	cout << endl;
	cout << "Spaltpiovotisierung? [j/n]";
	cin >> pivot;
	if (pivot == "j") {
		cout << "Gaussalgorithmus mit Pivotisierung" << endl << endl;
		int spalte = 0;
		int zeile = 0;
		while (spalte != (dimension - 1)) {
			//Pivotisierung
			float max= 0;
			int pivot;
			for (int i = zeile; i != dimension; i++) {
				if (abs(A[i][spalte]) > max) {
					max = abs(A[i][spalte]);
					pivot = i;
				}
			}
			// pivot-Element nach oben tauschen
			vector<float> Atemp = A[pivot]; //Zeilen tauschen
			float btemp = b[pivot];
			A[pivot] = A[zeile];
			b[pivot] = b[zeile];
			A[zeile] = Atemp;
			b[zeile] = btemp;


			for (int i = zeile + 1; i != (dimension); i++) {
				float factor = (-1 * A[i][spalte]) / A[zeile][spalte]; //factor herausfinden
				for (int j = spalte; j != (dimension); j++) {
					A[i][j] += factor * A[zeile][j];                         //Elemente umformen
				}
				b[i] += factor*b[zeile];


				/*for (int i = 0; i < A.size(); i++) {            // print result
					for (int j = 0; j < A[i].size(); j++) {
						cout << setw(6) << A[i][j] << "  ";
					}
					cout << setw(10) << "|";
					cout << setw(6) << b[i];
					cout << endl;
				}
				cout << endl;*/

			}
			spalte++;
			zeile = spalte;
		}
		for (int i = 0; i < A.size(); i++) {            // print result
			for (int j = 0; j < A[i].size(); j++) {
				cout << setw(10) << A[i][j] << "  ";
			}
			cout << setw(10) << "|";
			cout << setw(10) << b[i];
			cout << endl;
		}
		cout << endl;
	}
	else {
		cout << "Gaussalgorithmus ohne Pivotisierung" << endl << endl;
		int spalte = 0;
		int zeile = 0;
		while (spalte != (dimension - 1)) {
			//falls eine 0 an führender Stelle steht
			if (A[zeile][spalte] == 0) {
				for (int i = zeile; i != (dimension - 1); i++) {
					if (A[i][spalte] == 0) {
						vector<float> Atemp = A[i]; //Zeilen tauschen
						float btemp = b[i];
						A[i] = A[zeile];
						b[i] = b[zeile];
						A[zeile] = Atemp;
						b[zeile] = btemp;
						i = dimension - 1; //um die for-Schleife zu unterbrechen
					}
				}
			}

			
			for (int i = zeile+1; i != (dimension); i++) {
				float factor = (-1 * A[i][spalte]) / A[zeile][spalte]; //factor herausfinden
				for (int j = spalte; j != (dimension); j++) {
					A[i][j] += factor * A[zeile][j];                         //Elemente umformen
				}
				b[i] += factor * b[zeile];


				/*for (int i = 0; i < A.size(); i++) {            // print result
					for (int j = 0; j < A[i].size(); j++) {
						cout << setw(6) << A[i][j] << "  ";
					}
					cout << setw(10) << "|";
					cout << setw(6) << b[i];
					cout << endl;
				}
				cout << endl;*/

			}
			spalte++; 
			zeile = spalte;
		}
		for (int i = 0; i < A.size(); i++) {            // print result
			for (int j = 0; j < A[i].size(); j++) {
				cout << setw(10) << A[i][j] << "  ";
			}
			cout << setw(10) << "|";
			cout << setw(10) << b[i];
			cout << endl;
		}
		cout << endl;
	}
	//Lösungsvektor finden
	vector<float> X(dimension);
	for (int i = dimension-1; i >= 0; i--) {
		float sum= 0;
		for (int j=0; j != dimension; j++) {
			sum += X[j] * A[i][j];
		}
		X[i] = (b[i] - sum) / A[i][i];
	}
	cout << "Loesungsvektor X: " << endl;
	for (int i = 0; i <= dimension - 1; i++) {
		cout << X[i] << " | ";
	}
	cout << endl << endl;
}
void aufgabe15() {

}
int value_of_roem(char c) {
	if(c== 'M') {
		return 1000;
	}
	else if( c== 'D') {
		return 500;
	}
	else if (c == 'C') {
		return 100;
	}
	else if (c == 'L') {
		return 50;
	}
	else if (c == 'X') {
		return 10;
	}
	else if (c == 'V') {
		return 5;
	}
	else if (c == 'I') {
		return 1;
	}
	else {
		return 0;
	}
}
void aufgabe16() {
	cout << "Aufgabe 16; Roemische zu Arabisch" << endl;
	cout << "Bitte geben sie eine roemische Zahl ein: ";
	string s;
	int len, res, temp;
	cin >> s;
	char* roem = &(s[0]);
	len = strlen(roem);
// To upper case
	for (int i = 0; i <= len - 1; i++) {
		roem[i] = toupper(roem[i]);
	}
	res = 0;
	for (int i = 0; i <= len - 2; i++) {
		if(value_of_roem(roem[i]) < value_of_roem(roem[i + 1])) {
			res -= value_of_roem(roem[i]);
		}
		else {
			res += value_of_roem(roem[i]);
		}
	}
	res += value_of_roem(roem[len - 1]);
	cout << "Die Zahl ist: " << res << endl;
	
	return;
}

int main() {
	while (true) {
		int x;
		cout << "Welche Aufgabe soll ausgefuehrt werden? ";
		cin >> x;
		switch (x) {
		case(4):
			aufgabe4();
			break;
		case(5):
			aufgabe5();
			break;
		case(6):
			aufgabe6();
			break;
		case(7):
			aufgabe7();
			break;

		case(8):
			aufgabe8();
			break;
		case(10):
			aufgabe10();
			break;
		case(11):
			aufgabe11();
			break;

		case(13):
			aufgabe13();
			break;

		case(15):
			aufgabe15();
			break;
		case(16):
			aufgabe16();
			break;
		}
	}
	cout << endl;
}