#include "matrix_2d.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main() {
	try {
		string filename;
		cout << "Input filename with extension: ";
		getline(cin, filename);
		ifstream fin(filename);

		if (fin.good()) {
			size_t sz = 0;
			char fac_mode = 0;
			fin >> sz >> fac_mode;
			m2d::Matrix2D data(sz, sz);
			m2d::InputMatrix(fin, data);
			m2d::Matrix2D L(sz, sz), U(sz, sz);
			if (fac_mode == 'd') m2d::LUFactorizeDoolittle(data, L, U);
			else m2d::LUFactorizeCrout(data, L, U);
			cout << "L:\n";
			L.print();
			cout << "U:\n";
			U.print();
			char wait;
			cin >> wait;
			return 0;
		}
		else {
			char wait;
			cin >> wait;
			return -1;
		}
	}
	catch (exception & e) { // most exceptions thrown by m2d are of default c++ ones
		cerr << e.what() << '\n'; // and use what() for error messages.
		char wait;
		cin >> wait;
		return -2;
	}
}