#include "generalFunc.h"
#ifndef __Kmeans_structs
#define __Kmeans_structs

typedef struct processDeatiels {
	unsigned int myId, numOfSlaves;
	MPI_Status* status;
	unsigned initPointsIndex, pointsToCalc;
}process;
typedef struct centroid {
	double x, y, z;
	unsigned int index;
	unsigned int numOfPoints;
}centroid;
typedef struct point {
	unsigned int index;
	double xi, yi, zi;
	double vxi, vyi, vzi;
	double x, y, z;
	int closestCentroidIndex;
	double distanceFromCentroid;
}point;

///File opening. all processes read the first line of Constants, including Master
inline FILE* defineKmeansConstsFromFile(char* fileName,unsigned * N, unsigned* K, unsigned* LIMIT, double* T, double* dT, double* QM) {
	FILE* f;
	fopen_s(&f, fileName, "r");
	allocSuccess(1, f);
	fscanf_s(f, "%d%d%lf%lf%d%lf", N, K, T, dT, LIMIT, QM);
	return f;

}

inline void printOutput(char* fileName, centroid* centroids, unsigned int numOfCentroids, double q, double t) {
	FILE* f;
	fopen_s(&f, fileName, "w");
	fprintf(f, "----------t=%lf,q=%lf----------\n", t, q);
	for (unsigned int i = 0; i < numOfCentroids; i++)
		fprintf(f, "c%d (%lf,%lf%lf),numOfPoints=%d\n", i, centroids[i].x, centroids[i].y, centroids[i].z, centroids[i].numOfPoints);
	fclose(f);

}

inline centroid convertPointToCentroid(point* p,unsigned int index) {
	centroid chosenCentroid;
	chosenCentroid.x = p->x;
	chosenCentroid.y = p->y;
	chosenCentroid.z = p->z;
	chosenCentroid.numOfPoints = 0;
	chosenCentroid.index = index;
	return chosenCentroid;
}
inline double euclideanDistancePP(point * p0, point * p1) {
	double distance = sqrt(pow((p0->x - p1->x), 2) + pow((p0->y - p1->y), 2) + pow((p0->z - p1->z), 2));
	return distance;
}
inline double euclideanDistancePC(point * p, centroid * c) {
	double distance = sqrt(pow((p->x - c->x), 2) + pow((p->y - c->y), 2) + pow((p->z - c->z), 2));
	return distance;
}
inline double euclideanDistanceCC(centroid * c0, centroid * c1) {
	double distance = sqrt(pow((c0->x - c1->x), 2) + pow((c0->y - c1->y), 2) + pow((c0->z - c1->z), 2));
	return distance;
}
inline void printCentroids(centroid* centroids, unsigned int numOfCentroids) {
	printf("\n");
	for (unsigned int i = 0; i <numOfCentroids; i++) //centroid print
		printf("c%d=(%lf,%lf,%lf),numOfPoint=%d\n", centroids[i].index, centroids[i].x, centroids[i].y, centroids[i].z, centroids[i].numOfPoints);

	printf("\n");
}
inline void printPoints(point* points, unsigned int numOfPoints) {
	printf("\n");
	for (unsigned int i = 0; i < numOfPoints; i++) { //point print
		if (points[i].distanceFromCentroid != DBL_MAX)
			printf("p%d=initial:(%lf,%lf,%lf),current:(%lf,%lf,%lf),dis c%d=%lf\n", points[i].index, points[i].xi, points[i].yi, points[i].zi, points[i].x, points[i].y, points[i].z, points[i].closestCentroidIndex, points[i].distanceFromCentroid);
		else
			printf("p%d=initial(%lf,%lf,%lf),current:(%lf,%lf,%lf),dis c%d=Infinity\n", points[i].index, points[i].xi, points[i].yi, points[i].zi, points[i].x, points[i].y, points[i].z, points[i].closestCentroidIndex);
	}
	printf("\n");
}
inline centroid* initCentroids(unsigned int numOfCentroids) {
	centroid* centroids = (centroid*)calloc(numOfCentroids, sizeof(centroid));
	int i;
	#pragma omp parallel for
	for ( i = 0; i < numOfCentroids; i++) {
		centroids[i].index = i;
		centroids[i].numOfPoints = 0;
	}
	return centroids;

}

inline void updatePoints(point* points, unsigned int numOfPoints, centroid* centroids, unsigned int numOfCentroids) {//works only if points dis init to infinity
	int i;
	#pragma omp parallel for
	for (i=0 ; i < numOfPoints; i++){
		double minDistance = points[i].distanceFromCentroid;//all points distance initialize with Infinity

		for (unsigned int j = 0; j < numOfCentroids; j++){
			double newDis = euclideanDistancePC(&points[i], &centroids[j]);
			if (newDis < minDistance) {
				minDistance = newDis;
				points[i].distanceFromCentroid = minDistance;
				points[i].closestCentroidIndex = centroids[j].index;
			}
		}
	}
}
inline void updateCentroids(point* points, unsigned int numOfPoints, centroid* centroids) {//works only if centroids numOfPoints init to 0
	for (unsigned int i = 0; i < numOfPoints; i++)
		centroids[points[i].closestCentroidIndex].numOfPoints++;

	



}


inline void pointsInit(FILE* f, point* points, process* proc) {
	int i;
	
	#pragma omp parallel for
	for (i = 0; i <(int)proc->pointsToCalc; i++) {
		fscanf_s(f, "%lf%lf%lf%lf%lf%lf", &points[proc->initPointsIndex + i].xi, &points[proc->initPointsIndex + i].yi, &points[proc->initPointsIndex + i].zi, &points[proc->initPointsIndex + i].vxi, &points[proc->initPointsIndex + i].vyi, &points[proc->initPointsIndex + i].vzi);
		points[proc->initPointsIndex+i].index = proc->initPointsIndex + i;
		points[proc->initPointsIndex+i].closestCentroidIndex = -1;
		points[proc->initPointsIndex+i].distanceFromCentroid = DBL_MAX;//DBL_MAX is like Infinity
	}
}
inline point * findMaxDisPoint(point* points, unsigned int numOfPoints) {
	double max = points[0].distanceFromCentroid;
	point * maxDisPoint = &points[0];
	for (unsigned int i = 1; i < numOfPoints; i++)
		if (max < points[i].distanceFromCentroid) {
			max = points[i].distanceFromCentroid;
			maxDisPoint = &points[i];
		}

	return maxDisPoint;
}

inline double calc_q(point * points, unsigned int numOfPoints, centroid* centroids, unsigned int numOfcentroids) {
	if (numOfcentroids == 1)
		return -1;
	point ** pointsGropByCentoids = (point **)calloc(numOfcentroids, sizeof(point*));
	unsigned int* counters = (unsigned int*)calloc(numOfcentroids, sizeof(unsigned int));
	double* diameters = (double*)calloc(numOfcentroids, sizeof(double));
	double q = 0;

	allocSuccess(3, pointsGropByCentoids, counters, diameters);

	int i;
	#pragma omp parallel for
	for (i = 0; i < numOfcentroids; i++) {
		pointsGropByCentoids[i] = (point *)calloc(centroids[i].numOfPoints, sizeof(point));
		allocSuccess(1, pointsGropByCentoids[i]);
	}


	for (i = 0; i < numOfPoints; i++) {//sorted all points to matrix by centroids
		pointsGropByCentoids[points[i].closestCentroidIndex][counters[points[i].closestCentroidIndex]] = points[i];
		counters[points[i].closestCentroidIndex]++;
	}


	#pragma omp parallel for firstprivate(counters)
	for (i = 0; i < numOfcentroids; i++) {//find diameter& store in diametrs
		double diameter = 0;
		for (int j = 0; j < counters[i]; j++) {
			for (int k = j + 1; k < counters[i]; k++) {
				double newDiameter = euclideanDistancePP(&pointsGropByCentoids[i][j], &pointsGropByCentoids[i][k]);
				if (diameter < newDiameter)
					diameter = newDiameter;
			}
		}
		diameters[i] = diameter;
	}
	//for (int i = 0; i < numOfcentroids; i++)
	//	printf("c%d diameter=%lf,", i, diameters[i]);
	//printf("\n");
	#pragma omp parallel for firstprivate(diameters,centroids) shared(q)
	for (i = 0; i < numOfcentroids; i++)
		for (unsigned int j = 0; j < numOfcentroids; j++)
			if (i != j)
				q += (diameters[i] / euclideanDistanceCC(&centroids[i], &centroids[j]));

	q /= (numOfcentroids*(numOfcentroids - 1));

	#pragma omp parallel for
	for (i = 0; i < numOfcentroids; i++)
		free(pointsGropByCentoids[i]);
	freeAll(3, pointsGropByCentoids, counters, diameters);




	return q;



}



#endif // __Kmeans_structs


