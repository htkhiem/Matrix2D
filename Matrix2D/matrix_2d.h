/**
 * @file matrix_2d.h
 * Interface to general double-precision 2D matrices.
 */

#ifndef MATRIX2D_BASE
#define MATRIX2D_BASE

#ifdef MATRIX2D_EXPORTS
#define MATRIX2D_LIB __declspec(dllexport) // dll building
#else
#define MATRIX2D_LIB __declspec(dllimport) // usage
#endif

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <iomanip>

using namespace std;
namespace m2d {
	/** Main 2D Matrix class definition.
	* This matrix class is range-checked and compatible with const parameters.
	*/
	class MATRIX2D_LIB Matrix2D {
		double **elem; /** A 2D array of double-precision floating point numbers. */
		size_t size_x, size_y; /** Size of this matrix, which must always be initialised. */
	public:
		/** Simple constructor.
		* Initialises all elements to zero.
		* @param size_x: Can be understood as number of rows.
		* @param size_y: Can be understood as number of columns.
		*/
		Matrix2D(size_t size_x, size_t size_y);
		/** Copy constructor
		* @param src: The source Matrix2D object to copy from.
		*/
		Matrix2D(const Matrix2D& src);
		/** Destructor.
		* Deallocates 2D double array.
		*/
		~Matrix2D();
		/** Getter method
		* Returns the value at the specified location in the 2D double matrix.
		* @param pos_x: Row index of desired value.
		* @param pos_y: Column index of desired value.
		* @return a copy of the value at that location.
		* @exception: out_of_range() if supplied pos_x and/or pos_y are out-of-range of this matrix.
		*/
		double getAt(size_t pos_x, size_t pos_y) const;
		/** Setter method
		* Sets the cell at the specified location to the given value.
		* @param pos_x: Row index of target cell.
		* @param pos_y: Column index of target cell.
		* @param val: Value to be set to target cell.
		* @exception: out_of_range() if supplied pos_x and/or pos_y are out-of-range of this matrix. 
		*/
		void setAt(size_t pos_x, size_t pos_y, double val);
		/** Method to get row count.
		* @return The matrix's row count as std::size_t.
		*/
		size_t getSizeX() const { return size_x; }
		/** Method to get column count.
		* @return The matrix's column count as std::size_t.
		*/
		size_t getSizeY() const { return size_y; }
		/** Method to extract a submatrix from the current matrix.
		* The submatrix is defined by the position of its top-left element and its size (counting from that element).
		* @param pos_x: x-position of the top-left element of the desired submatrix.
		* @param pos_y: x-position of the top-left element of the desired submatrix.
		* @param sub_size_x: Row count of desired submatrix.
		* @param sub_size_y: Column count of desired submatrix.
		* @return A reference to the resulting submatrix.
		* @exception Throws range_error() if the submatrix exceeds the source's range, for example when its
		* top-left element is at the bottom-right of the source and its size is larger than 1 in any direction.
		*/
		Matrix2D subMatrix(size_t pos_x, size_t pos_y, size_t sub_size_x, size_t sub_size_y) const;
		/** Returns the cofactor of an element specified by its indices.
		* @param pos_x: The row-index of the specified element.
		* @param pos_y: The column-index of the specified element.
		* @return The cofactor of that element, as a matrix.
		* @exception: Throws range_error() if the supplied indices are out-of-range.
		*/
		Matrix2D cofactorOf(size_t pos_x, size_t pos_y) const;
		/** Computes the cofactor matrix of this matrix.
		* @return The cofactor matrix, through a copy constructor.
		*/
		Matrix2D cofactorMatrix();
		/** Prints the matrix.
		* Format: space-separated for now, iomanip later.
		*/
		void print() const;
		/** Checks if this is a square matrix.
		* @return True if matrix is square and false otherwise.
		*/
		bool isSquare() const { return (size_x == size_y); }
		/** Checks if this is an upper triangular matrix. Only square matrices can be triangular.
		* @return True if matrix is upper-triangular, false otherwise.
		*/
		bool isUpperTriangular() const;
		/** Checks if this is a lower triangular matrix. Only square matrices can be triangular.
		* @return True if matrix is lower-triangular, false otherwise.
		*/
		bool isLowerTriangular() const;
		/** Checks if this is a diagonal matrix. Only square matrices can be diagonal.
		* @return True if matrix is diagonal, false otherwise.
		*/
		bool isDiagonal() const;
		/**
		* Checks if this matrix is diagonally-dominant, that is, the sum of the absolute values of all elements
		* except the diagonal one in a row is less than the diagonal element.
		* @return True if it is diagonally dominant, false if it's not, or if it's non-square.
		*/
		bool isDiagonallyDominant() const;
		/** Computes and returns the determinant of this matrix, if it's square.
		* Uses a basic form of LU decomposition (assumes det(P) = 1).
		* @return The determinant of this matrix if it's square.
		* @exception Throws range_error() if this matrix isn't square.
		*/

		/// Matrix arithmetics
		/** Adds two matrices together, provided they're of the same size.
		* @return The resulting matrix.
		* @exception Throws invalid_argument() if the size of the matrices do not match.
		*/
		Matrix2D operator+(const Matrix2D& other) const;
		/** Subtract another matrix from this matrix, provided they're of the same size.
		* @return The resulting matrix.
		* @exception Throws invalid_argument() if the size of the matrices do not match.
		*/
		Matrix2D operator-(const Matrix2D& other) const;
		/** Multiply two matrices, provided this matrix's size_y is the same as the other's size_x.
		* @return The resulting matrix. Note that the resulting matrix might have a different size.
		* @exception Throws invalid_argument() if the sizes do not match the above criterion.
		*/
		Matrix2D operator*(const Matrix2D& other) const;
		/** Inverts this matrix. Not always possible.
		* @return The inverted matrix. Note that it also modifies the current matrix's value to that of the inverted one.
		*/
		void invert();


		
		double det() const;
		/** Method to transpose current matrix.
		* It writes to a new matrix, deallocates the current one then point *data to the new one.
		* @return The determinant, as a double-precision value.
		* @exception invalid_argument() if this matrix isn't square.
		*/
		void transpose() noexcept;
	};

	// Non-member functions

	/** Method to create a matrix from an std::ifstream.
	* The ifstream must contain a matrix written row-by-row and separated by spaces.
	* @param ifs: Reference to the ifstream containing input.
	* @param m: Reference to the target matrix. This function directly modifies its parameter instead of returning its own matrix.
	*/
	extern "C" void MATRIX2D_LIB InputMatrix(ifstream &ifs, Matrix2D &m);
	/** This function concatenates two given matrices horizontally.
	* The two matrices must have the same row count.
	* @param left: Reference to the first matrix.
	* @param right: Reference to the second matrix, which will be concatenated to the right of the first matrix.
	* @exception range_error() if the two matrices do not have the same row count.
	* @return A reference to the resulting matrix.
	*/
	MATRIX2D_LIB Matrix2D ConcatenateHorizontally(const Matrix2D &left, const Matrix2D &right);
	/** This function concatenates two given matrices vertically.
	* The two matrices must have the same column count.
	* @param top: Reference to the first matrix.
	* @param bottom: Reference to the second matrix, which will be concatenated below the first matrix.
	* @return A reference to the resulting matrix.
	* @exception range_error() if the two matrices do not have the same column count.
	*/
	MATRIX2D_LIB Matrix2D ConcatenateVertically(const Matrix2D &top, const Matrix2D &bottom);
	/** LU-Factoriser implementing Doolittle's method.
	* In accordance with Doolittle's method, it assumes the diagonal of the lower matrix L to be 1s (ones).
	* \param A: The source matrix, which must be square.
	* \param L: The lower triangular matrix, passed as reference to be written to. In this method, its diagonal is assumed to be 1s.
	* \param U: The upper triangular matrix, also passed as reference.
	* \exception range_error(): Throws when matrix A is not factorizable.
	*/
	extern "C" MATRIX2D_LIB void LUFactorizeDoolittle(const Matrix2D& A, Matrix2D& L, Matrix2D& U);
	/** LU-Factoriser implementing Crout's method.
	* In accordance with Crout's method, it assumes the diagonal of the upper matrix U to be 1s (ones).
	* \param A: The source matrix, which must be square.
	* \param L: The lower triangular matrix, passed as reference to be written to.
	* \param U: The upper triangular matrix, also passed as reference. In this method, its diagonal is assumed to be 1s.
	* \exception range_error(): Throws when matrix A is not factorizable.
	*/
	extern "C" MATRIX2D_LIB void LUFactorizeCrout(const Matrix2D& A, Matrix2D& L, Matrix2D& U);
}

#endif // MATRIX2D_BASE
 