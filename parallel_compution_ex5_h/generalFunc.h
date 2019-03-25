#include<stdio.h>
#include<mpi.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include<omp.h>
#include<time.h>
#include <float.h>
#include<math.h>
#ifndef __generalFunc
#define  __generalFunc


inline void allocSuccess(unsigned int allocsNum, ...) {
	va_list allcList;
	va_start(allcList, allocsNum);
	for (unsigned int i = 0; i < allocsNum; i++)
		if (va_arg(allcList, int*) == NULL) {
			printf("Error:one of your dynamic allocation has faild");
			MPI_Abort(MPI_COMM_WORLD, MPI_ERR_NO_MEM);
		}
	va_end(allcList);
}



inline  void freeAll(unsigned int allocsNum, ...) {

	va_list allcList;
	va_start(allcList, allocsNum);
	for (unsigned int i = 0; i < allocsNum; i++)
		free(va_arg(allcList, int*));
	va_end(allcList);


}
inline void forwardLineInFile(FILE* f, unsigned int numOfLine) {
	char ch = fgetc(f);
	for (unsigned int i = 0; i < numOfLine; i++) {
		if (ch == EOF)
			break;
		while (ch != EOF && ch != '\n')
			ch = fgetc(f);
		if (i != numOfLine - 1)
			ch = fgetc(f);
	}
}


#endif // !__generalFunc


