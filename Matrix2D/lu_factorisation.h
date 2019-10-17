/**
 * @file lu_factorisation.h
 * Interfaces to the LU-factorisation functions.
 */
#ifndef _MATRIX2D_NMETHODS_LUFACT_
#define _MATRIX2D_NMETHODS_LUFACT_

#include "matrix_2d.h"
#include <exception>
using namespace m2d;
 /** LU-Factoriser implementing Doolittle's method.
 * In accordance with Doolittle's method, it assumes the diagonal of the lower matrix L to be 1s (ones).
 * \param A: The source matrix, which must be square.
 * \param L: The lower triangular matrix, passed as reference to be written to. In this method, its diagonal is assumed to be 1s.
 * \param U: The upper triangular matrix, also passed as reference.
 * \exception range_error(): Throws when matrix A is not factorizable.
 */
void MATRIX2D_LIB LUFactorizeDoolittle(const Matrix2D& A, Matrix2D& L, Matrix2D& U);
/** LU-Factoriser implementing Crout's method.
* In accordance with Crout's method, it assumes the diagonal of the upper matrix U to be 1s (ones).
* \param A: The source matrix, which must be square.
* \param L: The lower triangular matrix, passed as reference to be written to.
* \param U: The upper triangular matrix, also passed as reference. In this method, its diagonal is assumed to be 1s.
* \exception range_error(): Throws when matrix A is not factorizable.
*/
void MATRIX2D_LIB LUFactorizeCrout(const Matrix2D& A, Matrix2D& L, Matrix2D& U);
#endif //_MATRIX2D_NMETHODS_LUFACT_
