#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

#include <cstdint>
#include <mpi.h>
#include <vector>

// Program Constants
#define nsteps   1000
#define savefreq 10
#define density  0.0005
#define mass     0.01
#define cutoff   0.01
#define min_r    (cutoff / 100)
#define dt       0.0005
#define parts_thresh 50000
// Particle Data Structure
typedef struct particle_t {
    uint64_t id; // Particle ID
    double x;    // Position X
    double y;    // Position Y
    double vx;   // Velocity X
    double vy;   // Velocity Y
    double ax;   // Acceleration X
    double ay;   // Acceleration Y
} particle_t;

typedef struct message_container_t {
    int particle_going_out_num;
    std::vector<particle_t> particle_going_out;

    int particle_coming_in_num;
    std::vector<particle_t> particle_coming_in;

    int num_ghost_particles_coming_in;
    std::vector<particle_t> ghost_particles_coming_in;

    int num_ghost_particles_going_out;
    std::vector<particle_t> ghost_particles_going_out;
} message_container_t;

extern MPI_Datatype PARTICLE;

// Simulation routine
void init_simulation(particle_t* parts, int num_parts, double size, int rank, int num_procs);
void simulate_one_step(particle_t* parts, int num_parts, double size, int rank, int num_procs);
void gather_for_save(particle_t* parts, int num_parts, double size, int rank, int num_procs);

#endif
