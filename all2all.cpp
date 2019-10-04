#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

#define BUFFER_SIZE 16400

int main(int argc,char *argv[])
{

/*
    int my_rank, num_workers;
    MPI_Comm SLAVES_WORLD;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_workers);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    //createSlavesCommunicator(&SLAVES_WORLD);

    char send_msg[20*num_workers], recv_buf[20*num_workers];
    int i;
    for(i=0;i<num_workers;i++){
        sprintf(&send_msg[i*20], "test from %d to %d", my_rank,i);
    } 

    MPI_Barrier(MPI_COMM_WORLD);
    //MPI_Alltoall(send_msg, strlen(send_msg), MPI_CHAR, recv_buf, 20, MPI_CHAR, MPI_COMM_WORLD);
    MPI_Alltoall(send_msg, 20, MPI_CHAR, recv_buf, 20, MPI_CHAR, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);


    for(i=0;i<num_workers;i++){
        printf("slave %d recvd message %s\n", my_rank, &recv_buf[20*i]);
    }
*/

	int size;
	int rank;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int* intSendPack = new int[BUFFER_SIZE]();// multiplier 3 for 3 different elements of index, newModules, oldModules and +1 for elementCount
	int* intReceivePack = new int[BUFFER_SIZE]();

	for (int prId = 0; prId < size; prId++) {
		if (prId != rank) {
			printf("testing message send from sender:%d and receiver:%d\n",
					rank, prId);
			MPI_Send(intSendPack, BUFFER_SIZE, MPI_INT, prId, 0,
			MPI_COMM_WORLD);
		  }
	}

	for (int sId = 0; sId < size; sId++) {
		if (sId != rank) {
			MPI_Recv(intReceivePack, BUFFER_SIZE, MPI_INT, sId, 0,
			MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			printf("testing message receive from sender:%d and receiver:%d\n",
					sId, rank);
		}
	}


    MPI_Finalize();
    return 0;
}
