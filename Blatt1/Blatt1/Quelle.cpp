#include <cmath>
#include <iostream>
using namespace std;

int main() {
	
	for (int i = 12; i > 0; i--) {
		int k = -1 * pow(10, i);
		cout << "exp(" << k << ") = " << exp(k) << endl;
	}
	for (int i = 0; i < 12; i++) {
		int k = pow(10, i);
		cout << "exp(" << k << ") = " << exp(k) << endl;
	}
	float x, res;
	cout << "Bitte eine Zahl eingeben: ";
	cin >> x;
	res = exp(x);
	cout << "exp(" << x << ") = " << res;

}