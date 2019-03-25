#include <cuda_runtime_api.h>
#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <device_functions.h>

#include<stdio.h>
#include"K-means_structs.h"
#include "constants.h"
#ifndef __cudaSection
#define __cudaSection
extern "C" unsigned int cudaMaxThredsPerBlock() {
	int nDevices;
	cudaGetDeviceCount(&nDevices);

	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, 0);
	return prop.maxThreadsDim[0];

}





//__global__ void cudaFindMaxArr(double * pointsDis, unsigned int*  numOfPoints) {// logN finding max insted of N
//	__shared__ double tempNumOfPoints;
//	tempNumOfPoints = *numOfPoints;
//	double maxDis = 0;
//	unsigned int threadID = blockIdx.x*blockDim.x + threadIdx.x;
//	while (threadID<tempNumOfPoints) {
//		if (threadID % 2 == 0 && threadID != tempNumOfPoints - 1) {
//			maxDis = (pointsDis[threadID] > pointsDis[threadID + 1]) ? pointsDis[threadID] : pointsDis[threadID + 1];
//			pointsDis[threadID / 2] = maxDis;
//		}
//		__syncthreads();
//
//		if (threadID == tempNumOfPoints - 1) {
//			maxDis = (pointsDis[threadID] > pointsDis[threadID - 1]) ? pointsDis[threadID] : pointsDis[threadID - 1];
//			unsigned int index = ((int)tempNumOfPoints % 2 == 1) ? (threadID / 2) : ((threadID - 1) / 2);
//			pointsDis[index] = maxDis;
//		}
//
//
//		tempNumOfPoints = ceil(tempNumOfPoints / 2);
//
//		__syncthreads();
//		if (tempNumOfPoints == 1)
//			break;
//	}
//}
//extern "C" void cudaFindMaxArrMediator(double * arr, unsigned int  arrSize, double* max) {
//	double *d_arr;
//	unsigned  int* d_arrSize;
//
//	cudaMalloc(&d_arr, arrSize * sizeof(double));
//	cudaMalloc(&d_arrSize, sizeof(unsigned int));
//
//	cudaMemcpy(d_arr, arr, arrSize * sizeof(double), cudaMemcpyHostToDevice);
//	cudaMemcpy(d_arrSize, &arrSize, sizeof(unsigned int), cudaMemcpyHostToDevice);
//
//	unsigned int numOfBlocks = (unsigned int)ceil(arrSize*1.0 / cudaMaxThredsPerBlock());
//	unsigned int blocksPerThreds = (numOfBlocks>1) ? cudaMaxThredsPerBlock() : arrSize;
//	cudaFindMaxArr << <numOfBlocks, blocksPerThreds >> > (d_arr, d_arrSize);
//
//	cudaMemcpy(max, &d_arr[0], sizeof(double), cudaMemcpyDeviceToHost);//update point
//
//	cudaFree(d_arr);
//	cudaFree(d_arrSize);
//
//}


__global__ void cudaCalcPoints(point* points, unsigned int*  numOfPoints, double* t) {

	unsigned int threadID = blockIdx.x*blockDim.x + threadIdx.x;

	if (threadID <*numOfPoints) {
		points[threadID].x = points[threadID].xi + points[threadID].vxi*(*t);
		points[threadID].y = points[threadID].yi + points[threadID].vyi*(*t);
		points[threadID].z = points[threadID].zi + points[threadID].vzi*(*t);

	}



}
extern "C" void cudaCalcPointsMediator(point *points, unsigned  int numOfPoints, double t) {
	point *d_points;
	unsigned  int* d_numOfPoints;
	double* d_t;

	cudaMalloc(&d_points, numOfPoints * sizeof(point));
	cudaMalloc(&d_numOfPoints, sizeof(unsigned int));
	cudaMalloc(&d_t, sizeof(double));

	cudaMemcpy(d_points, points, numOfPoints * sizeof(point), cudaMemcpyHostToDevice);
	cudaMemcpy(d_numOfPoints, &numOfPoints, sizeof(unsigned int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_t, &t, sizeof(double), cudaMemcpyHostToDevice);


	unsigned int numOfBlocks = (unsigned int)ceil(numOfPoints*1.0 / cudaMaxThredsPerBlock());
	unsigned int blocksPerThreds = (numOfBlocks>1) ? cudaMaxThredsPerBlock() : numOfPoints;
	cudaCalcPoints << <numOfBlocks, blocksPerThreds >> > (d_points, d_numOfPoints, d_t);//calc the GPU 2 qurter

	cudaMemcpy(points, d_points, numOfPoints * sizeof(point), cudaMemcpyDeviceToHost);//update point

	cudaFree(d_points);
	cudaFree(d_numOfPoints);
	cudaFree(d_t);
}

//__global__ void cudaCalcPointsDis(point *points, unsigned  int* numOfPoints, double * distances) {
//	unsigned int threadID = blockIdx.x*blockDim.x + threadIdx.x;
//	if (threadID != 0 && threadID < *numOfPoints) {//threadID != 0 because i don't want to calc dis of point from itself
//		double dis=cudaCalcEuclideanDistancePP(&points[0],&points[threadID]);
//		distances[threadID - 1] = dis;
//	}
//
//}

//extern "C" void cudaCalcPointsMaxDisMediator(point *points, unsigned  int numOfPoints,double* maxDis) {
//
//	if (numOfPoints == 2) {
//		*maxDis = euclideanDistancePP(&points[0], &points[1]);
//		return;
//	}
//	
//
//	point *d_points;
//	unsigned  int* d_numOfPoints;
//	double *d_distances;
//
//	cudaMalloc(&d_points, numOfPoints * sizeof(point));
//	cudaMalloc(&d_numOfPoints, sizeof(unsigned int));
//	cudaMalloc(&d_distances, sizeof(double)*(numOfPoints - 1));
//
//
//	cudaMemcpy(d_points, points, numOfPoints * sizeof(point), cudaMemcpyHostToDevice);
//	cudaMemcpy(d_numOfPoints, &numOfPoints, sizeof(unsigned int), cudaMemcpyHostToDevice);
//
//
//	unsigned int numOfBlocks = (unsigned int)ceil(numOfPoints*1.0 / cudaMaxThredsPerBlock());
//	unsigned int blocksPerThreds = (numOfBlocks>1) ? cudaMaxThredsPerBlock() : numOfPoints;
//
//
//	cudaCalcPointsDis << <numOfBlocks, blocksPerThreds >> > (d_points, d_numOfPoints, d_distances);////calc in O(1) dis of point[i] with all the points ahead
//
//
//	
//
//	numOfPoints--;
//	cudaMemcpy(d_numOfPoints, &numOfPoints, sizeof(unsigned int), cudaMemcpyHostToDevice);
//	cudaFindMaxArr << <numOfBlocks, blocksPerThreds >> > (d_distances, d_numOfPoints);//log(n)
//
//	cudaMemcpy(maxDis, &d_distances[0], sizeof(double), cudaMemcpyDeviceToHost);
//
//
//	cudaFree(d_points);
//	cudaFree(d_numOfPoints);
//	cudaFree(d_distances);
//
//}



//__global__ void  cudaSortPointsToGroups(point * points, unsigned int* numOfPoints, centroid* centroids, point ** pointsGropByCentoids) {
//	unsigned int threadID = blockIdx.x*blockDim.x + threadIdx.x;
//
//	if (threadID < *numOfPoints) {
//
//		pointsGropByCentoids[points[threadID].closestCentroidIndex][points[threadID].index%]
//	
//	}
//}
//extern "C" point ** cudaSortPointsToGroupsMediator(point * points, unsigned int numOfPoints, centroid* centroids, unsigned int numOfcentroids) {
//
//
//	point *d_points,**d_pointsGropByCentoids;
//	centroid * d_centroids;
//	unsigned  int* d_numOfPoints;
//
//	cudaMalloc(&d_points, numOfPoints * sizeof(point));
//	cudaMalloc(&d_numOfPoints, sizeof(unsigned int));
//	cudaMalloc(&d_pointsGropByCentoids, sizeof(point*)*numOfcentroids);
//	cudaMalloc(&d_centroids, sizeof(centroid)*numOfcentroids);
//
//	for (unsigned int i = 0; i < numOfcentroids; i++){
//		cudaMalloc(&d_pointsGropByCentoids[i], sizeof(point)*centroids[i].numOfPoints);
//	}
//
//
//	cudaMemcpy(d_points, points, numOfPoints * sizeof(point), cudaMemcpyHostToDevice);
//	cudaMemcpy(d_numOfPoints, &numOfPoints, sizeof(unsigned int), cudaMemcpyHostToDevice);
//	cudaMemcpy(d_centroids, centroids, sizeof(centroid)*numOfcentroids, cudaMemcpyHostToDevice);
//
//
//	unsigned int numOfBlocks = (unsigned int)ceil(numOfPoints*1.0 / cudaMaxThredsPerBlock());
//	unsigned int blocksPerThreds = (numOfBlocks>1) ? cudaMaxThredsPerBlock() : numOfPoints;
//	cudaSortPointsToGroups << <numOfBlocks, blocksPerThreds >> > (d_points, d_numOfPoints, d_centroids, d_pointsGropByCentoids);//calc the GPU 2 qurter
//
//}
//__global__ void cudaCalcMaxPointDis(point* points, unsigned int*  numOfPoints,double * pointsDis) {
//
//	unsigned int threadID = blockIdx.x*blockDim.x + threadIdx.x;
//	double dis = 0;
//	if(threadID<*numOfPoints)
//	for (unsigned int i = 0; i < *numOfPoints; i++){
//		if (threadID > i) {
//			dis = cudaCalcEuclideanDistancePP(&points[i], &points[threadID]);
//			pointsDis[threadID] = dis;
//		}
//		
//		__syncthreads();
//		if (threadID == 0) {
//			double localMax;
//			//cudaFindMaxArr<<<>>>(pointsDis, numOfPoints);
//			if (localMax > globalMaxDis)
//				globalMaxDis = localMax;
//		}
//		__syncthreads();
//
//	}
//
//		
//}


//__device__ double cudaCalcEuclideanDistancePP(point * p1, point* p2) {
//	double distance = sqrt(pow(p1->x - p2->x, 2) + pow(p1->y - p2->y, 2) + pow(p1->z - p2->z, 2));
//	return distance;
//}











extern "C" void printCudaDitels() {
	int nDevices;
	cudaGetDeviceCount(&nDevices);
	for (int i = 0; i < nDevices; i++) {
		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, i);
		//printf("Device Number: %d\n", i);
		//printf("Device name: %s\n", prop.name);
		//printf("Memory Clock Rate (KHz): %d\n",
		//	prop.memoryClockRate);
		//printf("Memory Bus Width (bits): %d\n",
		//	prop.memoryBusWidth);
		//printf("Peak Memory Bandwidth (GB/s): %f\n",
		//	2.0*prop.memoryClockRate*(prop.memoryBusWidth / 8) / 1.0e6);
		printf("Num of multi processor:%d\n", prop.multiProcessorCount);
		printf("Max Threads Per Multi Processor:%d\n", prop.maxThreadsPerMultiProcessor);
		printf("\n");

		printf("Max size of each dim of grid(%d,%d,%d)\n", prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
		printf("Max thread Dim (%d,%d,%d)\n", prop.maxThreadsDim[0], prop.maxThreadsDim[1], prop.maxThreadsDim[2]);
		printf("Max threads per block:%d\n", prop.maxThreadsPerBlock);

		printf("Max texture 3D:(%d,%d,%d)\n", prop.maxTexture3D[0], prop.maxTexture3D[1], prop.maxTexture3D[2]);
		printf("Max texture 2D:(%d,%d)\n", prop.maxTexture2D[0], prop.maxTexture2D[1]);

	}

}





#endif // !__cudaSection



