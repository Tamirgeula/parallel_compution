#include "constants.h"
#include "K-means_structs.h"
#include "mpi_structs.h"
#ifndef __main
#define __main

unsigned int* openMpInit(process* proc) {
	unsigned int* procsThreads = (unsigned int*)calloc(proc->numOfSlaves+1, sizeof(unsigned int));
	allocSuccess(1, procsThreads);
	procsThreads[proc->myId] = omp_get_max_threads();//each procces get and set his max threads abilities!

	omp_set_num_threads(procsThreads[proc->myId]);
	for (unsigned int i = 0; i < proc->numOfSlaves && proc->myId == MASTER; i++) {
		MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, proc->status);
		MPI_Recv(&procsThreads[proc->status->MPI_SOURCE], 1, MPI_UNSIGNED, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, proc->status);
	}
	if (proc->myId != MASTER)
		MPI_Send(&procsThreads[proc->myId], 1, MPI_UNSIGNED, MASTER, 0, MPI_COMM_WORLD);
	MPI_Bcast(procsThreads, proc->numOfSlaves+1, MPI_UNSIGNED, MASTER, MPI_COMM_WORLD);
	return procsThreads;
}



int main(int argc, char * argv[])
{
	
	/////////////k mean constns
	unsigned int N,K,LIMIT;
	double T,dT,QM;
	FILE* f = defineKmeansConstsFromFile(INPUT_FILE_NAME, &N, &K, &LIMIT, &T, &dT, &QM);
	/////////////end of k mean constns


	///////////MPI Master-slave standart init
	unsigned int numOfProcs , myId , numOfSlaves;
	MPI_Status status;
	MPI_MyInit(&argc, &argv, &numOfProcs, &myId,&numOfSlaves,&MPI_CENTROID,&MPI_POINT);//in this func if numOfSlave=1, then program Aborts+cunstructs speacial structs!
	process p = processInit(myId,numOfSlaves, &status,N);//static divide for points resp for each proccess+convinient for function argument transfer
	///////////End of MPI Master-slave standart init


	/////OpenMp init
	unsigned int* procsThreads= openMpInit(&p);//arr of each proc active threads acording to location
	/////end of OpenMp init

	///structs for k means
	point * points = (point*)calloc(N, sizeof(point));
	centroid* centroids = initCentroids(K);
	allocSuccess(2, points, centroids);
	int reachedQM;
	///structs for k means



	if (myId == MASTER) {

		printf("\nN=%d,K=%d,dT=%lf,T=%lf,LIMIT=%d,QM=%lf\n",N,K,dT,T,LIMIT,QM);
		for (double i = 0; i <=T; i += dT) {
				MPI_RandomGatherPointsToMaster(&p, points, N);
				MPI_Bcast(points, N, MPI_POINT, MASTER, MPI_COMM_WORLD);
				KMeansPlusPlusInit(&p, points, N, centroids, K);//smart initialize
				Kmeans(&p, points, N, centroids, K, LIMIT);//each proc responsiable for chunck of point
				double q=calc_q(points, N, centroids, K);
				printf("----------dT=%lf----------\n",i);
				printf("q=%lf\n", q);
				printCentroids(centroids, K);
				printf("----------dT=%lf----------\n\n",i);
				if (q <= QM) {
					reachedQM = 1;
					MPI_Bcast(&reachedQM, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
					printOutput(OUTPUT_FILE_NAME,centroids, K, q, i);
					break;
				}
				MPI_Bcast(&reachedQM, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
		}//end of time interval

	}//end of master

	else {//slaves
		//////file points reading
		forwardLineInFile(f, p.initPointsIndex +1);//aimes each slave's pointer to it's location
	    pointsInit(f, points, &p);
		////end of file points reading 

		for (double i = 0; i <=T; i+=dT){
			cudaCalcPointsMediator(&points[p.initPointsIndex], p.pointsToCalc, i);
			MPI_RandomGatherPointsToMaster(&p, points, N);
			MPI_Bcast(points, N, MPI_POINT, MASTER, MPI_COMM_WORLD);////all slaves get global points from Master
			KMeansPlusPlusInit(&p, points, N, centroids, K);//smart initialize
			Kmeans(&p, points, N, centroids, K, LIMIT);//each proc responsiable for chunck of point
			MPI_Bcast(&reachedQM, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
			if (reachedQM == 1)
				break;
		}
	}//end of slaves section

	freeAll(3, points, centroids, procsThreads);
	fclose(f);
	MPI_Finalize();
	return 0;
}

#endif // !1
