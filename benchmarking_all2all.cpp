#include <iostream>
#include <algorithm>
#include <mpi.h>

#define BUFFER_SIZE 16384

void point2point(int*, int*, int, int);

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int rank_id = 0, com_sz = 0;
    double t0 = 0.0, tf = 0.0;
    MPI_Comm_size(MPI_COMM_WORLD, &com_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);

    int* intSendPack = new int[BUFFER_SIZE]();
    int* result = new int[BUFFER_SIZE*com_sz]();
    std::fill(intSendPack, intSendPack + BUFFER_SIZE, rank_id);
    std::fill(result + BUFFER_SIZE*rank_id, result + BUFFER_SIZE*(rank_id+1), rank_id);

    // Send-Receive
    t0 = MPI_Wtime();
    point2point(intSendPack, result, rank_id, com_sz);
    MPI_Barrier(MPI_COMM_WORLD);
    tf = MPI_Wtime();
    if (!rank_id)
        std::cout << "Send-receive time: " << tf - t0 << std::endl;

    // Collective
    std::fill(result, result + BUFFER_SIZE*com_sz, 0);
    std::fill(result + BUFFER_SIZE*rank_id, result + BUFFER_SIZE*(rank_id+1), rank_id);
    t0 = MPI_Wtime();
    MPI_Allgather(intSendPack, BUFFER_SIZE, MPI_INT, result, BUFFER_SIZE, MPI_INT, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    tf = MPI_Wtime();
    if (!rank_id)
        std::cout << "Allgather time: " << tf - t0 << std::endl;

    MPI_Finalize();
    delete[] intSendPack;
    delete[] result;
    return 0;
}

// Send/receive communication
void point2point(int* send_buf, int* result, int rank_id, int com_sz)
{
    MPI_Status status;
    // Exchange and store the data
    for (int i=0; i<com_sz; i++){
        if (i != rank_id){
            MPI_Sendrecv(send_buf, BUFFER_SIZE, MPI_INT, i, 0, 
                result + i*BUFFER_SIZE, BUFFER_SIZE, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
        }
    }
}
