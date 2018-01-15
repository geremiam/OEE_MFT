// IO.h
/* Input and output in different formats */
#ifndef IO_H
#define IO_H

void PrintMatrix(const int num_rows, const int num_cols, const double*const mat, std::ostream& out);
void PrintMatrix(const int num_rows, const int num_cols, const double*const*const mat, std::ostream& out);
void PrintMatrix(const int num_rows, const int num_cols, const std::complex<double>*const mat, std::ostream& out);
void PrintMatrix(const int num_rows, const int num_cols, const std::complex<double>*const*const mat, std::ostream& out);

#endif
