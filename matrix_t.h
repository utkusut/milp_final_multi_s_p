//////////////////////////////////////////////////
//
// This is the header file for allocating
// Dynamic Memory for arrays 
//
//		matrix_t.h
//
//
//	This program allocates dynamic memory
//
//	By Salih Ramazan
//
///////////////////////////////////////////////////

#ifndef MATRIX_T_H
#define MATRIX_T_H

#include <iostream>

using namespace std;

template <class T>
class MATRIX{
private:
public:
	void newMTX(T* &matrix, unsigned n1max);
	void newMTX(T** &matrix, unsigned n1max, unsigned n2max);
	void newMTX(T*** &matrix, unsigned n1max, unsigned n2max, unsigned n3max);
	void newMTX(T**** &matrix, unsigned n1max, unsigned n2max,
		unsigned n3max, unsigned n4max);
	void newMTX(T***** &matrix, unsigned n1max, unsigned n2max,
		unsigned n3max, unsigned n4max, unsigned n5max);
	void delMTX(T* &matrix, unsigned n1max=0);
	void delMTX(T** &matrix, unsigned n1max, unsigned n2max=0);
	void delMTX(T*** &matrix, unsigned n1max, unsigned n2max, unsigned n3max=0);
	void delMTX(T**** &matrix, unsigned n1max, unsigned n2max,
		unsigned n3max, unsigned n4max=0);
	void delMTX(T***** &matrix, unsigned n1max, unsigned n2max,
		unsigned n3max, unsigned n4max, unsigned n5max=0);
};


// 1D Matrix
template <class T>
void MATRIX<T>::newMTX(T* &matrix, unsigned n1max)
{
	if (( matrix = new T [n1max]) == NULL) {
		cerr << "Error in allocating n1max array" <<endl;
	}
}

template <class T>
void MATRIX<T>::delMTX(T* &matrix, unsigned n1max)
{
	n1max=0;
	delete [] matrix;
	matrix = NULL;
}


// 2D Matrix
template <class T>
void MATRIX<T>::newMTX(T** &matrix, unsigned n1max, unsigned n2max)
{
	if (( matrix = new T *[n1max]) == NULL) {
		cerr << "Error in allocating n1max array" <<endl;
	}
	for (unsigned n1 = 0; n1 < n1max; ++n1) {
		if (( matrix[n1] = new T [n2max]) == NULL) {
			cerr << "Error in allocating n1max array" << endl;
		}
	}
}

template <class T>
void MATRIX<T>::delMTX(T** &matrix, unsigned n1max, unsigned n2max)
{
	n2max=0;
	for (unsigned n1 = n1max; n1 > 0; --n1) {
		delete [] matrix[n1-1];
		matrix[n1-1] = NULL;
	}
	delete [] matrix;
	matrix = NULL;
}


// 3D Matrix
template <class T>
void MATRIX<T>::newMTX(T*** &matrix, unsigned n1max,
						  unsigned n2max, unsigned n3max)
{
	if (( matrix = new T **[n1max]) == NULL) {
		cerr << "Error in allocating n1max array" <<endl;
	}
	for (unsigned n1 = 0; n1 < n1max; ++n1) {
		if (( matrix[n1] = new T *[n2max]) == NULL) {
			cerr << "Error in allocating n1max array" << endl;
		}
		for (unsigned n2 = 0; n2 < n2max; ++n2) {
			if (( matrix[n1][n2] = new T [n3max]) == NULL) {
				cerr << "Error in allocating n1max array" << endl;
			}
		}	// n2 loop
	}	// n1 loop
}

template <class T>
void MATRIX<T>::delMTX(T*** &matrix, unsigned n1max,
						   unsigned n2max, unsigned n3max)
{
	n3max=0;
	for (unsigned n1 = n1max; n1 > 0; --n1) {
		for (unsigned n2 = n2max; n2 > 0; --n2) {
			delete [] matrix[n1-1][n2-1];
			matrix[n1-1][n2-1] = NULL;
		}	// n2 loop
		delete [] matrix[n1 - 1];
		matrix[n1-1] = NULL;
	}	// n1 loop
	delete [] matrix;
	matrix = NULL;
}


// 4D Matrix
template <class T>
void MATRIX<T>::newMTX(T**** &matrix, unsigned n1max, unsigned n2max,
						  unsigned n3max, unsigned n4max)
{
	if (( matrix = new T ***[n1max]) == NULL) {
		cerr << "Error in allocating n1max array" <<endl;
	}
	for (unsigned n1 = 0; n1 < n1max; ++n1) {
		if (( matrix[n1] = new T **[n2max]) == NULL) {
			cerr << "Error in allocating n1max array" << endl;
		}
		for (unsigned n2 = 0; n2 < n2max; ++n2) {
			if (( matrix[n1][n2] = new T *[n3max]) == NULL) {
				cerr << "Error in allocating n1max array" << endl;
			}
			for (unsigned n3 = 0; n3 < n3max; ++n3) {
				if (( matrix[n1][n2][n3] = new T [n4max]) == NULL) {
					cerr << "Error in allocating n1max array" << endl;
				}
			}	// n3 loop
		}	// n2 loop
	}	// n1 loop
}

template <class T>
void MATRIX<T>::delMTX(T**** &matrix, unsigned n1max, unsigned n2max,
						   unsigned n3max, unsigned n4max)
{
	n4max=0;
	for (unsigned n1 = n1max; n1 > 0; --n1) {
		for (unsigned n2 = n2max; n2 > 0; --n2) {
			for (unsigned n3 = n3max; n3 > 0; --n3) {
				delete [] matrix[n1-1][n2-1][n3-1];
				matrix[n1-1][n2-1][n3-1] = NULL;
			}	// n3 loop
			delete [] matrix[n1-1][n2-1];
			matrix[n1-1][n2-1] = NULL;
		}	// n2 loop
		delete [] matrix[n1-1];
		matrix[n1-1] = NULL;
	}	// n1 loop
	delete [] matrix;
	matrix = NULL;
}


// 5D Matrix
template <class T>
void MATRIX<T>::newMTX(T***** &matrix, unsigned n1max, unsigned n2max,
						  unsigned n3max, unsigned n4max, unsigned n5max)
{
	if (( matrix = new T ****[n1max]) == NULL) {
		cerr << "Error in allocating n1max array" <<endl;
	}
	for (unsigned n1 = 0; n1 < n1max; ++n1) {
		if (( matrix[n1] = new T ***[n2max]) == NULL) {
			cerr << "Error in allocating n1max array" << endl;
		}
		for (unsigned n2 = 0; n2 < n2max; ++n2) {
			if (( matrix[n1][n2] = new T **[n3max]) == NULL) {
				cerr << "Error in allocating n1max array" << endl;
			}
			for (unsigned n3 = 0; n3 < n3max; ++n3) {
				if (( matrix[n1][n2][n3] = new T *[n4max]) == NULL) {
					cerr << "Error in allocating n1max array" << endl;
				}
				for (unsigned n4 = 0; n4 < n4max; ++n4) {
					if (( matrix[n1][n2][n3][n4] = new T [n5max]) == NULL) {
						cerr << "Error in allocating n1max array" << endl;
					}
				}	// n5 loop
			}	// n3 loop
		}	// n2 loop
	}	// n1 loop
}

template <class T>
void MATRIX<T>::delMTX(T***** &matrix, unsigned n1max, unsigned n2max,
						   unsigned n3max, unsigned n4max, unsigned n5max)
{
	for (unsigned n1 = n1max; n1 > 0; --n1) {
		for (unsigned n2 = n2max; n2 > 0; --n2) {
			for (unsigned n3 = n3max; n3 > 0; --n3) {
				for (unsigned n4 = n4max; n4 > 0; --n4) {
					delete [] matrix[n1-1][n2-1][n3-1][n4-1];
					matrix[n1-1][n2-1][n3-1][n4-1] = NULL;
				}	// n4 loop
				delete [] matrix[n1-1][n2-1][n3-1];
				matrix[n1-1][n2-1][n3-1] = NULL;
			}	// n3 loop
			delete [] matrix[n1-1][n2-1];
			matrix[n1-1][n2-1] = NULL;
		}	// n2 loop
		delete [] matrix[n1-1];
		matrix[n1-1] = NULL;
	}	// n1 loop
	delete [] matrix;
	matrix = NULL;
}

#endif