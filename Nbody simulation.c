#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <unistd.h>

#define MAXPROC 384

const double G = 6.67259e-7;
const double dt = 1.0;


int write_particles(int N, double *X, double *Y, char *fn) {
  FILE *fp;

  if ((fp = fopen(fn, "w")) == NULL) {
    printf("Couldn't open file %s\n", fn);
    return 0;
  }

  for (int i = 0; i < N; i++) {
    fprintf(fp, "%3.2f %3.2f \n", X[i], Y[i]);
  }
  
  fprintf(fp, "\n");
  fclose(fp);
  return (1);
}


double dist(double px, double py, double qx, double qy) {
  return sqrt(pow(px - qx, 2) + pow(py - qy, 2));
}


void compute_force(int N, double *X, double *Y, double *mass, double *Fx, double *Fy, int start_index, int end_index) {
  const double mindist = 0.0001;

  for (int i = start_index; i < end_index; i++) {                       
    Fx[i] = Fy[i] = 0.0;
    for (int j = 0; j < N; j++) {
      if (i != j) {
        double r = dist(X[i], Y[i], X[j], Y[j]);
        if (r > mindist) {
          double r3 = pow(r, 3);
          Fx[i] += G * mass[i] * mass[j] * (X[j] - X[i]) / r3;
          Fy[i] += G * mass[i] * mass[j] * (Y[j] - Y[i]) / r3;
        }
      }
    }
  }
}


void check_process_count(int number_of_processes, int id) {
  if (number_of_processes < 2 || number_of_processes > MAXPROC) {
    if (id == 0) {
      printf("\n\nYou have to use at lest 2 and at most %d processes\n", MAXPROC);
    }

    MPI_Finalize();
    exit(0);
  }
}


void generate_mass_and_position(double *mass, double *X, double *Y, int id, int bodies_per_process) {
  const double size = 100.0;
  short int seedval[3] = {5 + id, 5 + id, 5 + id};
  seed48(seedval);

  for (int i = 0; i < bodies_per_process; i++) {
    mass[i] = 1000.0 * drand48();     // 0 <= mass < 1000
    X[i] = size * drand48();          // 0 <= X < 100
    Y[i] = size * drand48();          // 0 <= Y < 100
  }
}


int main(int argc, char **argv) {
  int number_of_processes, index, id;

  const int N = 100;               // Number of bodies
  const int timesteps = 1000;     // Number of timesteps
  const double size = 100.0; 

  MPI_Init(&argc, &argv);        
  MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  MPI_Status status;      
  MPI_Request recv_req[MAXPROC];
  
  check_process_count(number_of_processes, id);

  double *mass;          /* mass of bodies */
  double *X;             /* x-positions of bodies */
  double *Y;             /* y-positions of bodies */
  double *Vx;            /* velocities on x-axis of bodies */
  double *Vy;            /* velocities on y-axis of bodies */
  double *Fx;            /* forces on x-axis of bodies */
  double *Fy;            /* forces on y-axis of bodies */

  const int comm_tag = 12;

  int bodies_per_process = (N + number_of_processes - 1) / number_of_processes;

  if (id == 0) {
    printf("\nN-body simulation, number of bodies = %d \n", N);
    printf("\nNo of Bodies per Process = %d \n", bodies_per_process);
  }

  mass     = (double *) calloc(N, sizeof(double));      // Mass 
  X        = (double *) calloc(N, sizeof(double));      // Position (x,y) at current time step
  Y        = (double *) calloc(N, sizeof(double));
  Vx       = (double *) calloc(N, sizeof(double));      // Velocities
  Vy       = (double *) calloc(N, sizeof(double));
  Fx       = (double *) calloc(N, sizeof(double));      // Forces
  Fy       = (double *) calloc(N, sizeof(double));


  short int seedval[3] = {7, 7, 7};
  seed48(seedval);

  for (int i = 0; i < N ; i++) {
    mass[i]  = 1000.0 * drand48();    // 0 <= mass < 1000
    X[i] = size * drand48();          // 0 <= X < 100
    Y[i] = size * drand48();          // 0 <= Y < 100
  }

  if(id == 0) {
    write_particles(N, X, Y, "initial_pos.txt");
  }
  
  int start_index = id * bodies_per_process;
  int end_index = start_index + bodies_per_process;

  // computing local forces
  compute_force(N, X, Y, mass, Fx, Fy, start_index, end_index);
  
  // Set up the velocity vectors caused by initial forces for Leapfrog method
  for(int i = start_index; i < end_index; i++){
    Vx[i] = 0.5 * dt * Fx[i] / mass[i];
    Vy[i] = 0.5 * dt * Fy[i] / mass[i];
  }

  // creating structure to sned X and Y coordinates of the bodies
  struct particles_elements {
    double X_arr[bodies_per_process];
    double Y_arr[bodies_per_process];
  };
  
  // copying structure to memory
  struct particles_elements p1; 
  memcpy(p1.X_arr, X, sizeof(double[bodies_per_process]));
  memcpy(p1.Y_arr, Y, sizeof(double[bodies_per_process]));

  int t = 0;
  while (t < timesteps) {    // Loop for this many timesteps
    t++;

    if (id == 0) { 
      printf("%d ", t); 
      fflush(stdout);  // Print out the timestep
    }
    
    // Calculate new positions
    int k = 0;
    for (int i = start_index; i < end_index; i++){
      X[i] = X[i] + Vx[i] * dt;
      Y[i] = Y[i] + Vy[i] * dt;

      // copying positions back to structure
      p1.X_arr[k] = X[i];
      p1.Y_arr[k] = Y[i];
      k++;
    }

    // sending positions to all other processes
    for (int i = 0; i < number_of_processes; i++) {
      if (i != id) {
        MPI_Send(&p1, sizeof(struct particles_elements), MPI_BYTE, i, comm_tag, MPI_COMM_WORLD);
      }
    }

    // receiving positions from all other processes
    for(int i = 0; i < number_of_processes; i++) {
      if(i != id) {
        struct particles_elements p3;  
        MPI_Recv(&p3, sizeof(struct particles_elements), MPI_BYTE, i, comm_tag, MPI_COMM_WORLD, &status);

        // updating local positions array
        for(int j = 0; j < bodies_per_process; j++){
            X[i * bodies_per_process + j] = p3.X_arr[j];
            Y[i * bodies_per_process + j] = p3.Y_arr[j];
        }
      }
    }

    // computing global forces
    compute_force(N, X, Y, mass, Fx, Fy, start_index, end_index);

    // updating velocities
    for(int i = start_index; i < end_index; i++){
      Vx[i] = Vx[i] + Fx[i] * dt / mass[i];
      Vy[i] = Vy[i] + Fy[i] * dt / mass[i];
    }
  }

  if(id == 0) {
    write_particles(N, X, Y, "final_pos.txt");
  }

  printf("\n\n");
  exit(0);
}
