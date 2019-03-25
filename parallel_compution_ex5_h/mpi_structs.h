#include<mpi.h>
#include <stddef.h>
#include <stdio.h>
#include "K-means_structs.h"
#ifndef __mpi_structs
#define __mpi_structs


MPI_Datatype MPI_CENTROID, MPI_POINT;


inline void createPoint(MPI_Datatype* MPI_POINT) {
	MPI_Datatype point_data_types[] = { MPI_UNSIGNED,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE
		,MPI_DOUBLE,MPI_DOUBLE,MPI_INT,MPI_DOUBLE };
	int  blocklen_point[12] = { 1,1,1,1,1,1,1,1,1,1,1,1 };
	MPI_Aint disp_point[12] = { offsetof(point,index),offsetof(point,xi),offsetof(point,yi),offsetof(point,zi),offsetof(point,vxi),offsetof(point,vyi),
		offsetof(point,vzi), offsetof(point,x), offsetof(point,y), offsetof(point,z), offsetof(point,closestCentroidIndex),
		offsetof(point,distanceFromCentroid) };
	MPI_Type_create_struct(12, blocklen_point, disp_point, point_data_types, MPI_POINT);
	MPI_Type_commit(MPI_POINT);
}
inline void createCentroid(MPI_Datatype* MPI_CENTROID) {
	MPI_Datatype centroid_data_types[] = { MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_UNSIGNED,MPI_UNSIGNED };
	int blocklen_centroid[5] = { 1,1,1,1,1 };
	MPI_Aint disp_centroid[5] = { offsetof(centroid,x) ,offsetof(centroid, y),offsetof(centroid, z),offsetof(centroid,index), offsetof(centroid, numOfPoints) };
	MPI_Type_create_struct(5, blocklen_centroid, disp_centroid, centroid_data_types, MPI_CENTROID);
	MPI_Type_commit(MPI_CENTROID);
}

inline void MPI_MyInit(int* argc, char *** argv, unsigned int *numprocs, unsigned int *myId, unsigned int *numOfSlaves, MPI_Datatype* MPI_CENTROID, MPI_Datatype* MPI_POINT) {

	///////////MPI Master-slave standart init
	MPI_Init(argc, argv);
	MPI_Comm_size(MPI_COMM_WORLD, numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, myId);
	if (*numprocs == 1) {
		printf("Erorr! number of processes have to be bigger than 1.");
		MPI_Abort(MPI_COMM_WORLD, MPI_ERR_COUNT);
	}
	*numOfSlaves = *numprocs - 1;
	///////////End of MPI Master-slave standart init



	///building MPI_CENTROID&MPI_POINT struct
	createCentroid(MPI_CENTROID);
	createPoint(MPI_POINT);
	///end of building MPI_CENTROID&MPI_POINT struct


}

inline void MPI_RandomGatherPointsToMaster(process* proc, point * points, unsigned int numOfPoints) {
	if (proc->myId == MASTER) {
		int bufferSize = 0;
		for (unsigned int i = 0; i < numOfPoints;) {
			MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, proc->status);
			int bufferSize = (proc->status->MPI_SOURCE == proc->numOfSlaves) ? (numOfPoints / proc->numOfSlaves + numOfPoints%proc->numOfSlaves) : (numOfPoints / proc->numOfSlaves);
			MPI_Recv(&points[i], bufferSize, MPI_POINT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, proc->status);//store points at Master
			i += bufferSize;
		}

	}

	else
		MPI_Send(&points[proc->initPointsIndex], proc->pointsToCalc, MPI_POINT, MASTER, 0, MPI_COMM_WORLD);



}
inline void KMeansPlusPlusInit(process* proc, point * points, unsigned int numOfPoints, centroid * centroids, unsigned int numOfCentorids) {

	centroids[0] = convertPointToCentroid(&points[0], 0);
	if (proc->myId == MASTER) {
		point * slavesRemotePoints = calloc(proc->numOfSlaves, sizeof(point));
		allocSuccess(1, slavesRemotePoints);
		for (unsigned int i = 1; i < numOfCentorids; i++) {
			for (unsigned int j = 0; j < proc->numOfSlaves; j++)
				MPI_Recv(&slavesRemotePoints[j], 1, MPI_POINT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, proc->status);

			point * finalRemotePoint = findMaxDisPoint(slavesRemotePoints, proc->numOfSlaves);

			centroids[i] = convertPointToCentroid(finalRemotePoint, i);
			MPI_Bcast(&centroids[i], 1, MPI_CENTROID, MASTER, MPI_COMM_WORLD);
		}
		free(slavesRemotePoints);


	}
	else {
		updatePoints(&points[proc->initPointsIndex], proc->pointsToCalc, centroids, 1);//updates slave's chunck points acording first centroid.

		for (unsigned int i = 1; i < numOfCentorids; i++) {
			point * slaveRemotePoint = findMaxDisPoint(&points[proc->initPointsIndex], proc->pointsToCalc);
			MPI_Send(slaveRemotePoint, 1, MPI_POINT, MASTER, 0, MPI_COMM_WORLD);
			MPI_Bcast(&centroids[i], 1, MPI_CENTROID, MASTER, MPI_COMM_WORLD);
			updatePoints(&points[proc->initPointsIndex], proc->pointsToCalc, &centroids[i], 1);//updates slave's chunck points acording first centroid.

		}
		updateCentroids(&points[proc->initPointsIndex], proc->pointsToCalc, centroids);
	}
	MPI_RandomGatherPointsToMaster(proc, points, numOfPoints);

	MPI_Bcast(points, numOfPoints, MPI_POINT, MASTER, MPI_COMM_WORLD);

	for (unsigned int i = 0; i < numOfCentorids; i++) {
		unsigned int centroidNumOfPoints = 0;
		MPI_Allreduce(&centroids[i].numOfPoints, &centroidNumOfPoints, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
		centroids[i].numOfPoints = centroidNumOfPoints;
	}


}
inline process  processInit(unsigned int myId, unsigned int numOfSlaves, MPI_Status* status, unsigned int numOfPoints) {
	unsigned int initPointIndex = (numOfPoints / numOfSlaves)*(myId - 1);
	unsigned int pointsToCalc = (myId == numOfSlaves) ? (numOfPoints / numOfSlaves + numOfPoints%numOfSlaves) : (numOfPoints / numOfSlaves);
	process proc = { myId,numOfSlaves,status,initPointIndex,pointsToCalc };
	if (myId == MASTER) {//master don't calc points
		proc.pointsToCalc = 0;
		proc.initPointsIndex = 0;
	}
	return proc;
}
inline void Kmeans(process* proc, point * points, unsigned int numOfPoints, centroid * centroids, unsigned int numOfCentorids, double LIMIT) {
	centroid * oldCentroids = (centroid*)calloc(numOfCentorids, sizeof(centroid));
	allocSuccess(1, oldCentroids);
	for (unsigned int i = 0; i < LIMIT; i++) {
		int j;

		memcpy(oldCentroids, centroids, numOfCentorids * sizeof(centroid));
		unsigned int allOptimize = 0;
		#pragma omp parallel for
		for (j = 0; j < numOfCentorids; j++) {
			centroids[j].x = 0;
			centroids[j].y = 0;
			centroids[j].z = 0;
			centroids[j].numOfPoints = 0;
		}	

		for (j = 0; j < proc->pointsToCalc; j++) {//sum all cordinate and init points
			centroids[points[proc->initPointsIndex + j].closestCentroidIndex].x += points[proc->initPointsIndex + j].x;
			centroids[points[proc->initPointsIndex + j].closestCentroidIndex].y += points[proc->initPointsIndex + j].y;
			centroids[points[proc->initPointsIndex + j].closestCentroidIndex].z += points[proc->initPointsIndex + j].z;
			points[proc->initPointsIndex + j].distanceFromCentroid = DBL_MAX;//later in update point dis to new centroid val
																			 //point won't update dis if dis dcreased,so init to Infinity
		}
		for (j = 0; j < numOfCentorids; j++) {//send to master the sum of all local centroids, and sum in master to global val before avrge
			double tempX = centroids[j].x, tempY = centroids[j].y, tempZ = centroids[j].z;
			MPI_Reduce(&tempX, &centroids[j].x, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
			MPI_Reduce(&tempY, &centroids[j].y, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
			MPI_Reduce(&tempZ, &centroids[j].z, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
		}
		if (proc->myId == MASTER) {
			#pragma omp parallel for
			for (j = 0; j < numOfCentorids; j++) {// avrage sum in master

				if (oldCentroids[j].numOfPoints != 0) {
					centroids[j].x /= oldCentroids[j].numOfPoints;
					centroids[j].y /= oldCentroids[j].numOfPoints;
					centroids[j].z /= oldCentroids[j].numOfPoints;
				}
			}
		}
		
		MPI_Bcast(centroids, numOfCentorids, MPI_CENTROID, MASTER, MPI_COMM_WORLD);//update everyone in new centroid cordinates
		updatePoints(&points[proc->initPointsIndex], proc->pointsToCalc, centroids, numOfCentorids);//update points new dis from centroids
		updateCentroids(&points[proc->initPointsIndex], proc->pointsToCalc, centroids);// update each local centroids.numOfPoints
		MPI_RandomGatherPointsToMaster(proc, points, numOfPoints);//send to master updated points dis

		for (unsigned int j = 0; j < numOfCentorids; j++) {//update centroids.numOfPoints to global scale at Master
			unsigned int temp = centroids[j].numOfPoints;
			MPI_Reduce(&temp, &centroids[j].numOfPoints, 1, MPI_UNSIGNED, MPI_SUM, MASTER, MPI_COMM_WORLD);
		}
		if (proc->myId == MASTER) {
			#pragma omp parallel for shared(allOptimize)
			for (j = 0; j < numOfCentorids; j++) {//check convertion in master only
				double move = euclideanDistanceCC(&centroids[j], &oldCentroids[j]);
				if (move == 0)
					allOptimize++;
			}
		}

		MPI_Bcast(&allOptimize, 1, MPI_UNSIGNED, MASTER, MPI_COMM_WORLD);//if convertion aplly notify all slaves to break loop
		if (allOptimize == numOfCentorids) //break condition
			break;

	}
	free(oldCentroids);

}

#endif // !__mpi_structs