Matrix2D
========

Basic 2D double-precision matrix data structure with some useful methods
packaged into a DLL.

Originally from the LU_factoriser repository, now made universal for use
with all my hobby projects.

Building instructions
---------------------

Simply clone the repository to your machine and use Visual Studio 2017 and up to
build. After building, there should be a .DLL file in the output folder.

Usage for simple projects
-------------------------

Drop the matrix_2d.h header file into your project and \#include it in your
source code. The Matrix2D class and its associated functions live in the m2d
namespace, as a note.

After compiling your code, drop matrix_2d.dll into the same folder as your
executable.

// TODO: add examples and a GitHub Pages wiki.
