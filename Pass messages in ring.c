
#define K 1024            /* One Kilobyte */
#define M K*K             /* One Megabyte */
#define MAXSIZE K*M       /* One Gigabyte */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>


int main(int argc, char *argv[]) {
    const int tag = 42;
    int id, ntasks, source_id, i;
    MPI_Status status;
    int inputSize;
    double starttime, endtime;
    int process_stop_flag = 0;

    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
        printf("MPI_init failed!\n");
        exit(1);
    }

    // Number of processes
    if (MPI_Comm_size(MPI_COMM_WORLD, &ntasks) != MPI_SUCCESS) {
        printf("MPI_Comm_size failed!\n");
        exit(1);
    }

    // Get id of this process
    if (MPI_Comm_rank(MPI_COMM_WORLD, &id) != MPI_SUCCESS) {
        printf("MPI_Comm_rank failed!\n");
        exit(1);
    }

    // Check that we run on at least two processors
    if (ntasks < 2) {
        printf("You have to use at least 2 processors to run this program\n");
        MPI_Finalize();
        exit(0);
    }

    while(1){
        // Process 0 is the sender of inputSize
        if (id == 0) {
            printf("Please give an input size in bytes:\n");
            fflush(stdout);
            scanf("%d", &inputSize);

            if(inputSize <= 0){
                process_stop_flag = 1;
                MPI_Finalize();
                exit(0);
            }

            if (inputSize > MAXSIZE) {
                printf("Input size is too large, maximum value is %d\n", MAXSIZE);
                inputSize = 0;
            }

            printf("Process %d - sending %d\n", id, inputSize);

            for (i = 1; i < ntasks; i++) {
                if (MPI_Send(&inputSize, 1, MPI_INT, i, tag, MPI_COMM_WORLD) != MPI_SUCCESS) {
                    printf("Process %i: Error in MPI_Send!\n", id);
                    exit(1);
                }
            }

            char *buffer = (char*) malloc(inputSize*sizeof(char));
            if (buffer == NULL) {
                printf("Could not allocate memory, exiting\n");
                MPI_Finalize();
                exit(0);
            }
            // else{
            //     printf("Allocation is successful for process: %d\n", id);
            // }

            // Sending actual msg to 1
            if (MPI_Send(buffer, inputSize, MPI_CHAR, id+1, tag, MPI_COMM_WORLD) != MPI_SUCCESS) {
                printf("Process %d : Error in MPI_Send !\n", id);
                exit(1);
            }
            else{
                starttime = MPI_Wtime();    //starting time after sending msg to process 1
                //printf("Sending msg from process %d to process %d\n", id, id+1);
            } 

            if (MPI_Recv(buffer, inputSize, MPI_CHAR, ntasks-1, tag, MPI_COMM_WORLD, &status) != MPI_SUCCESS) {
                printf("Error in actual msg MPI_Recv!\n");
                exit(1);
            }
            else{
                endtime = MPI_Wtime();
                //source_id = status.MPI_SOURCE;
                //printf("Great !! Recived msg from process %d to process %d\n", source_id, id);
            }

            printf("That took %f seconds\n", endtime-starttime);
        }
        else {  // Processes 1 to N-1 are the recievers

            if (process_stop_flag == 1){
                MPI_Finalize();
                exit(0);
            }

            if (MPI_Recv(&inputSize, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status) != MPI_SUCCESS) {
                printf("Error in MPI_Recv!\n");
                exit(1);
            }

            // Get id of sender
            source_id = status.MPI_SOURCE;
            //printf("Process %d - got %d from %d\n", id, inputSize, source_id);

            char *buffer = (char*) malloc(inputSize*sizeof(char));
            if (buffer == NULL) {
                printf("Could not allocate memory, exiting\n");
                MPI_Finalize();
                exit(0);
            }
            // else {
            //     printf("Allocation is successful for process: %d\n", id);
            // }

            if (MPI_Recv(buffer, inputSize, MPI_CHAR, id-1, tag, MPI_COMM_WORLD, &status) != MPI_SUCCESS) {
                printf("Error in msg MPI_Recv!\n");
                exit(1);
            }
            else{
                source_id = status.MPI_SOURCE;
                //printf("Recived msg from process %d to process %d\n", source_id, id);
            }

            if (id != ntasks-1) {
                if (MPI_Send(buffer, inputSize, MPI_CHAR, id+1, tag, MPI_COMM_WORLD) != MPI_SUCCESS) {
                    printf("Process %d : Error in MPI_Send !\n", id);
                    exit(1);
                }
                // else{
                //     printf("Sending msg from process %d to process %d\n", id, id+1);
                // }
            }
            else{
                if (MPI_Send(buffer, inputSize, MPI_CHAR, 0, tag, MPI_COMM_WORLD) != MPI_SUCCESS) {
                    printf("Process %d : Error in MPI_Send !\n", id);
                    exit(1);
                }
                // else {
                //     printf("Sending msg from process %d to process 0\n", id);
                // }
            }
        }
    }
    

    // All processes do this
    if (MPI_Finalize() != MPI_SUCCESS) {
        printf("Error in MPI_Finalize!\n");
        exit(1);
    }
    exit(0);
}