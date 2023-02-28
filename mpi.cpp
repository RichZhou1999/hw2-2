/*
The idea of this mpi code:
    -- split the whole region into sub regions. The number of subregions
    is equal to the number of threads given. Each region has the width of
    "size" and the height of ("size"/ number of threads)
    -- the boundary of each subregion requires communication. Each region has
    two boundaries, upper and lower. After one step, transefer the particles on
    the boundary to nearby subregions. These particles are needed to calculate
    forces
-- region scheme:
    __________________
    |________________| <- subregion 1
    |________________| <- subregion 2
    |________________| <- subregion 3
    |________________| <- subregion 4
                                 */

#include "common.h"
#include <mpi.h>
#include <cmath>
#include <set>
#include <unordered_map>
#include <iostream>
#include <bits/stdc++.h>
#include <fstream>
#include <sstream>

std::unordered_map<int, std::unordered_set<particle_t *>> bins;
std::unordered_map<int, int> particle_to_bin;
double zone_size;
double bin_size = cutoff;
int row_lda;
int column_lda;
double upper_boundary;
double lower_boundary;
int particle_going_out_num;
std::vector<particle_t> particle_going_out;

int particle_going_out_upper_num;
std::vector<particle_t> particle_going_out_upper;

int particle_going_out_lower_num;
std::vector<particle_t> particle_going_out_lower;

int particle_possible_coming_in_upper_num;
std::vector<particle_t> particle_possible_coming_in_upper;
int particle_possible_coming_in_lower_num;
std::vector<particle_t> particle_possible_coming_in_lower;

int particle_beyond_upper_boundary_num;
std::vector<particle_t> particle_beyond_upper_boundary;
std::unordered_map<int, std::unordered_set<particle_t *>> particle_beyond_upper_boundary_bins;

int particle_beyond_lower_boundary_num;
std::vector<particle_t> particle_beyond_lower_boundary;
std::unordered_map<int, std::unordered_set<particle_t *>> particle_beyond_lower_boundary_bins;

int upper_boudary_particle_num;
std::vector<particle_t> particle_upper_boundary;
int lower_boudary_particle_num;
std::vector<particle_t> particle_lower_boundary;

std::vector<particle_t *> flatten_particles;

std::vector<particle_t > particle_send_gathering;
int number_particles_sending;
int number_particles_receiving_sum;
std::vector<particle_t > particle_receive_gathering;
std::vector<int > displacement;
std::vector<int> number_particles_receiving;
int step = 0;

// std::vector<int > displacement;
// Put any static global variables here that you will use throughout the simulation.


int calculate_bin_number(double x, double y, double size, double bin_size, int row_lda, int column_lda){
    if ((y < upper_boundary)){
        return -1;
    }
    if (( y >= lower_boundary)){
        return -2;
    }
    y = y - upper_boundary;
    double quotient;
    quotient = y/bin_size;
    int row = (int) quotient;

    quotient = x/bin_size;
    int column = (int) quotient;

    return row* column_lda + column;

}

void update_boundary_particles(){
    particle_upper_boundary.clear();
    upper_boudary_particle_num = 0;
    particle_lower_boundary.clear();
    lower_boudary_particle_num = 0;

    for(int i = 0; i < column_lda; i++){
        for (auto it = bins[i].begin(); it != bins[i].end(); ++it){
            particle_upper_boundary.push_back(**it);
            upper_boudary_particle_num += 1;
        }
        // for (auto it = bins[i + column_lda].begin(); it != bins[i+ column_lda].end(); ++it){
        //     particle_upper_boundary.push_back(**it);
        //     upper_boudary_particle_num += 1;
        // }
        for (auto it = bins[i + (row_lda-2)*column_lda].begin(); it != bins[i+ (row_lda-2)*column_lda].end(); ++it){
            particle_lower_boundary.push_back(**it);
            lower_boudary_particle_num += 1;
        }
        for (auto it = bins[i + (row_lda-1)*column_lda].begin(); it != bins[i+ (row_lda-1)*column_lda].end(); ++it){
            particle_lower_boundary.push_back(**it);
            lower_boudary_particle_num += 1;
        }
    }
}



void apply_force_one_direction(particle_t& particle, particle_t& neighbor) {
    // Calculate Distance
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = (1 - cutoff / r) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
    // return std::make_tuple(coef * dx, coef * dy);
}

void apply_force_bi_direction(particle_t& particle, particle_t& neighbor) {
    // Calculate Distance
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = (1 - cutoff / r) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
    neighbor.ax -= coef * dx;
    neighbor.ay -= coef * dy;
    // return std::make_tuple(coef * dx, coef * dy);
}


bool check_boundary(int row , int column){
    if ((row < 0) or (row>=row_lda)){
        return false;
    }
    if ((column < 0) or (column >=column_lda)){
        return false;
    }
    return true;
}


void move(particle_t& p, double size) {
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method

    int origin_bin = calculate_bin_number(p.x,p.y, size,bin_size,row_lda, column_lda);
    // int origin_bin = particle_to_bin[p.id];
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;

    // Bounce from walls
    while (p.x < 0 || p.x > size) {
        p.x = p.x < 0 ? -p.x : 2 * size - p.x;
        p.vx = -p.vx;
    }

    while (p.y < 0 || p.y > size) {
        p.y = p.y < 0 ? -p.y : 2 * size - p.y;
        p.vy = -p.vy;
    }

    int new_bin =calculate_bin_number(p.x, p.y, size, bin_size, row_lda, column_lda);
    if(origin_bin == new_bin) return;

    if (new_bin >= 0 ){
        // std::cout<< new_bin << "\n";
        bins[origin_bin].erase(&p);
        bins[new_bin].insert(&p);
    }
    // else{
    //     bins[origin_bin].erase(&p);
    //     particle_going_out.push_back(p);
    //     particle_going_out_num += 1;
    // }
    if (new_bin == -1){
        bins[origin_bin].erase(&p);
        particle_going_out_upper.push_back(p);
        particle_going_out_upper_num += 1;
    }

    if (new_bin == -2){
        bins[origin_bin].erase(&p);
        particle_going_out_lower.push_back(p);
        particle_going_out_lower_num += 1;
    }

}


void apply_force_bin(int row, int column, int row2, int column2){
    // Force every particle in each bin to interact
    if (!check_boundary(row, column) && !check_boundary(row2, column2)){
        return;
    }
    for (auto it2 = bins[column2+row2*column_lda].begin(); it2 != bins[column2+row2*column_lda].end(); ++it2){
        for (auto it = bins[column+row*column_lda].begin(); it != bins[column+row*column_lda].end(); ++it){
            // Interact particles
            apply_force_bi_direction(**it, **it2);
        }
    }
}

void apply_force_bin_self(int row, int column, int row2, int column2){
    // Force every particle in each bin to interact
    if (!check_boundary(row, column) && !check_boundary(row2, column2)){
        return;
    }
    for (auto it2 = bins[column2+row2*column_lda].begin(); it2 != bins[column2+row2*column_lda].end(); ++it2){
        for (auto it = bins[column+row*column_lda].begin(); it != bins[column+row*column_lda].end(); ++it){
            // Interact particles
            apply_force_one_direction(**it, **it2);
        }
    }
}

// void apply_force_bin_upper_boundary(int row, int column, int row2, int column2){
//     // Force every particle in each bin to interact
//     if (!check_boundary(row, column) && !check_boundary(row2, column2)){
//         return;
//     }
//     for (auto it2 = particle_beyond_upper_boundary.begin(); it2 != particle_beyond_upper_boundary.end(); ++it2){
//         for (auto it = bins[column+row*column_lda].begin(); it != bins[column+row*column_lda].end(); ++it){
//             // Interact particles
//             apply_force_bi_direction(**it,*it2);
//         }
//     }
// }

void apply_force_bin_upper_boundary(int row, int column, int row2, int column2){
    // Force every particle in each bin to interact
    if (!check_boundary(row, column) && !check_boundary(row2, column2)){
        return;
    }
    if(column != 0){
        for (auto it2 = particle_beyond_upper_boundary_bins[column-1].begin(); it2 != particle_beyond_upper_boundary_bins[column-1].end(); ++it2){
            for (auto it = bins[column+row*column_lda].begin(); it != bins[column+row*column_lda].end(); ++it){
                // Interact particles
                apply_force_bi_direction(**it,**it2);
            }
        }
    }
    for (auto it2 = particle_beyond_upper_boundary_bins[column].begin(); it2 != particle_beyond_upper_boundary_bins[column].end(); ++it2){
        for (auto it = bins[column+row*column_lda].begin(); it != bins[column+row*column_lda].end(); ++it){
            // Interact particles
            apply_force_bi_direction(**it,**it2);
        }
    }
    if(column != (column_lda -1)){
        for (auto it2 = particle_beyond_upper_boundary_bins[column+1].begin(); it2 != particle_beyond_upper_boundary_bins[column+1].end(); ++it2){
            for (auto it = bins[column+row*column_lda].begin(); it != bins[column+row*column_lda].end(); ++it){
                // Interact particles
                apply_force_bi_direction(**it,**it2);
            }
        }
    }

}

// void apply_force_bin_lower_boundary(int row, int column, int row2, int column2){
//     // Force every particle in each bin to interact
//     if (!check_boundary(row, column) && !check_boundary(row2, column2)){
//         return;
//     }
//     for (auto it2 = particle_beyond_lower_boundary.begin(); it2 != particle_beyond_lower_boundary.end(); ++it2){
//         for (auto it = bins[column+row*column_lda].begin(); it != bins[column+row*column_lda].end(); ++it){
//             // Interact particles
//             apply_force_bi_direction(**it,*it2);
//         }
//     }
// }


void apply_force_bin_lower_boundary(int row, int column, int row2, int column2){
    // Force every particle in each bin to interact
    if (!check_boundary(row, column) && !check_boundary(row2, column2)){
        return;
    }
    if (column != 0){
        for (auto it2 = particle_beyond_lower_boundary_bins[column-1].begin(); it2 != particle_beyond_lower_boundary_bins[column-1].end(); ++it2){
            for (auto it = bins[column+row*column_lda].begin(); it != bins[column+row*column_lda].end(); ++it){
                // Interact particles
                apply_force_bi_direction(**it,**it2);
            }
        }
    }
    for (auto it2 = particle_beyond_lower_boundary_bins[column].begin(); it2 != particle_beyond_lower_boundary_bins[column].end(); ++it2){
        for (auto it = bins[column+row*column_lda].begin(); it != bins[column+row*column_lda].end(); ++it){
            // Interact particles
            apply_force_bi_direction(**it,**it2);
        }
    }
    if (column != (column_lda -1)){
        for (auto it2 = particle_beyond_lower_boundary_bins[column+1].begin(); it2 != particle_beyond_lower_boundary_bins[column+1].end(); ++it2){
            for (auto it = bins[column+row*column_lda].begin(); it != bins[column+row*column_lda].end(); ++it){
                // Interact particles
                apply_force_bi_direction(**it,**it2);
            }
        }
    }

}

// void test_send_recv(int rank, int num_procs){
//     int temp = rank;
//     int recv_temp;
//     if (rank != (num_procs -1)){
//            MPI_Send(
//             &temp,
//             1,
//             MPI_INT,
//             (rank + 1)%num_procs,
//             0,
//             MPI_COMM_WORLD);
//     }
//     if (rank !=0 ){
//         MPI_Recv(&recv_temp,
//         1,
//         MPI_INT,
//         (rank - 1),
//         0,
//         MPI_COMM_WORLD,
//         MPI_STATUS_IGNORE);
//         std::cout<<"Processor"<< rank << "received" << recv_temp <<"\n";
//     }
// }

void send_recv_boundary_particles(int rank , int num_procs){
    if(rank != (num_procs -1)){
        MPI_Send(&lower_boudary_particle_num,
                 1,
                 MPI_INT,
                 rank+1,
                 0,
                 MPI_COMM_WORLD);
        MPI_Send(&particle_lower_boundary[0],
                 lower_boudary_particle_num,
                 PARTICLE,
                 rank+1,
                 0,
                 MPI_COMM_WORLD);
    }
    if(rank != 0){
        MPI_Recv(&particle_beyond_upper_boundary_num,
                 1,
                 MPI_INT,
                 rank-1,
                 0,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        //std::cout<<"zone: "<< rank << "received upper" << particle_beyond_upper_boundary_num<< " paticles" <<"\n";
        particle_beyond_upper_boundary.resize(particle_beyond_upper_boundary_num);
        MPI_Recv(&particle_beyond_upper_boundary[0],
                 particle_beyond_upper_boundary_num,
                 PARTICLE,
                 rank-1,
                 0,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

    }

    if(rank != 0){
        MPI_Send(&upper_boudary_particle_num,
                 1,
                 MPI_INT,
                 rank-1,
                 0,
                 MPI_COMM_WORLD);
        MPI_Send(&particle_upper_boundary[0],
                 upper_boudary_particle_num,
                 PARTICLE,
                 rank-1,
                 0,
                 MPI_COMM_WORLD);
    }


    if(rank != (num_procs -1)){
        MPI_Recv(&particle_beyond_lower_boundary_num,
                 1,
                 MPI_INT,
                 rank+1,
                 0,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        //std::cout<<"zone: "<< rank << "received lower" << particle_beyond_lower_boundary_num<< " paticles" <<"\n";
        particle_beyond_lower_boundary.resize(particle_beyond_lower_boundary_num);
        MPI_Recv(&particle_beyond_lower_boundary[0],
                 particle_beyond_lower_boundary_num,
                 PARTICLE,
                 rank+1,
                 0,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    }

}
void init_simulation(particle_t* parts, int num_parts, double size, int rank, int num_procs) {
    zone_size = size / num_procs;


    double quotient = zone_size / bin_size;
    row_lda = (int) ceil(quotient);

    quotient = size / bin_size;
    column_lda = (int) ceil(quotient);

    upper_boundary = zone_size * rank;
    lower_boundary = zone_size * (rank + 1);

    // std::cout<<"upper_boundary" <<upper_boundary << "rank" << rank << "\n";
    // std::cout<<"lower_boundary" <<lower_boundary << "rank" << rank << "\n";
    const int space = ceil(1.5 * bin_size * bin_size * 1. / density);
    for(int i = 0; i< row_lda*column_lda; ++i){
        bins[i].reserve(space);
    }
    // std::ofstream myfile;

    // std::ostringstream oss;
    // oss << "example" << rank << ".txt";
    // std::string name = oss.str();

    // myfile.open (name);
    for (int i = 0; i < num_parts; ++i){
        int index;
        index = calculate_bin_number(parts[i].x,parts[i].y, size, bin_size, row_lda, column_lda);
        // myfile<<" y " << parts[i].y <<  " rank " << rank << " index " << index << "\n";
        if (index >= 0){
            bins[index].insert(&parts[i]);
        }
        // bins[index].insert(&parts[i]);
        // particle_to_bin[parts[i].id] = index;
    }
    // myfile.close();

    number_particles_sending = 0;
    for( int i =0 ; i < row_lda; i++){
        for( int j =0; j < column_lda; j++){
            for (auto it = bins[j+i*column_lda].begin(); it != bins[j+i*column_lda].end(); ++it){
                number_particles_sending+=1;
            }
        }
    }
    // std::cout << "rank "<<rank<< "num" << number_particles_sending << "\n";

    update_boundary_particles();
    send_recv_boundary_particles(rank, num_procs);
    // if(rank != (num_procs -1)){
    //     MPI_Send(&lower_boudary_particle_num,
    //             1,
    //             MPI_INT,
    //             rank+1,
    //             0,
    //             MPI_COMM_WORLD);
    //     MPI_Send(&particle_lower_boundary[0],
    //     lower_boudary_particle_num,
    //     PARTICLE,
    //     rank+1,
    //     0,
    //     MPI_COMM_WORLD);
    // }
    // if(rank != 0){
    //     MPI_Recv(&particle_beyond_upper_boundary_num,
    //             1,
    //             MPI_INT,
    //             rank-1,
    //             0,
    //             MPI_COMM_WORLD,
    //             MPI_STATUS_IGNORE);
    //     std::cout<<"zone: "<< rank << "received upper" << particle_beyond_upper_boundary_num<< " paticles" <<"\n";
    //     particle_beyond_upper_boundary.resize(particle_beyond_upper_boundary_num);
    //     MPI_Recv(&particle_beyond_upper_boundary[0],
    //             particle_beyond_upper_boundary_num,
    //             PARTICLE,
    //             rank-1,
    //             0,
    //             MPI_COMM_WORLD,
    //             MPI_STATUS_IGNORE);

    // }

    // if(rank != 0){
    //     MPI_Send(&upper_boudary_particle_num,
    //             1,
    //             MPI_INT,
    //             rank-1,
    //             0,
    //             MPI_COMM_WORLD);
    //     MPI_Send(&particle_upper_boundary[0],
    //     upper_boudary_particle_num,
    //     PARTICLE,
    //     rank-1,
    //     0,
    //     MPI_COMM_WORLD);
    // }


    // if(rank != (num_procs -1)){
    //     MPI_Recv(&particle_beyond_lower_boundary_num,
    //             1,
    //             MPI_INT,
    //             rank+1,
    //             0,
    //             MPI_COMM_WORLD,
    //             MPI_STATUS_IGNORE);
    //     std::cout<<"zone: "<< rank << "received lower" << particle_beyond_lower_boundary_num<< " paticles" <<"\n";
    //     particle_beyond_lower_boundary.resize(particle_beyond_lower_boundary_num);
    //     MPI_Recv(&particle_beyond_lower_boundary[0],
    //             particle_beyond_lower_boundary_num,
    //             PARTICLE,
    //             rank+1,
    //             0,
    //             MPI_COMM_WORLD,
    //             MPI_STATUS_IGNORE);
    // }


    // if(rank != (num_procs -1)){
    //     MPI_Send(&particle_lower_boundary[0],
    //     lower_boudary_particle_num,
    //     PARTICLE,
    //     rank+1,
    //     0,
    //     MPI_COMM_WORLD);
    // }
    // if(rank != 0){
    //     particle_beyond_upper_boundary.resize(particle_beyond_upper_boundary_num);
    //     MPI_Recv(&particle_beyond_upper_boundary[0],
    //             particle_beyond_upper_boundary_num,
    //             PARTICLE,
    //             rank-1,
    //             0,
    //             MPI_COMM_WORLD,
    //             MPI_STATUS_IGNORE);
    // }
    // for(int i=0; i < particle_beyond_upper_boundary.size(); i++){
    //     std::cout <<particle_beyond_upper_boundary[i].id <<"\n";
    // }

    // if(rank != (num_procs -1)){
    //     MPI_Send(&particle_lower_boundary,
    //             lower_boudary_particle_num,
    //             PARTICLE,
    //             rank+1,
    //             0,
    //             MPI_COMM_WORLD);
    // }



    // You can use this space to initialize data objects that you may need
    // This function will be called once before the algorithm begins
    // Do not do any particle simulation here
}
void apply_force_bins(int row, int column){
    // Create a matrix of force for each bin to bin

    // Loop over 5 forward neighbor bins (3 below, 1 to left, and itself)
    apply_force_bin_self(row, column, row, column);
    // right neighbor
    apply_force_bin(row, column, row, column + 1);
    // bottom left neighbor
    apply_force_bin(row, column, row + 1, column - 1);
    // bottom neighbor
    apply_force_bin(row, column, row + 1, column);
    // bottom right neighbor
    apply_force_bin(row, column, row + 1, column + 1);
}

void apply_force_bins_upper_boundary(int row, int column){
    // Create a matrix of force for each bin to bin

    // Loop over 5 forward neighbor bins (3 below, 1 to left, and itself)
    apply_force_bin_self(row, column, row, column);
    // right neighbor
    apply_force_bin(row, column, row, column + 1);
    // bottom left neighbor
    apply_force_bin(row, column, row + 1, column - 1);
    // bottom neighbor
    apply_force_bin(row, column, row + 1, column);
    // bottom right neighbor
    apply_force_bin(row, column, row + 1, column + 1);
    apply_force_bin_upper_boundary(row, column, row, column);
}

void apply_force_bins_lower_boundary(int row, int column){
    // Create a matrix of force for each bin to bin

    // Loop over 5 forward neighbor bins (3 below, 1 to left, and itself)
    apply_force_bin_self(row, column, row, column);
    // right neighbor
    apply_force_bin(row, column, row, column + 1);
    apply_force_bin_lower_boundary(row, column, row, column);
}
void update_flatten_particles(){
    flatten_particles.clear();
    flatten_particles.reserve(ceil(1.5 * zone_size * bin_size*column_lda * 1. / density));
    for( int i =0 ; i < row_lda; i++){
        for( int j =0; j < column_lda; j++){
            for (auto it = bins[j+i*column_lda].begin(); it != bins[j+i*column_lda].end(); ++it){
                // particle_t* ptr = *it;
                flatten_particles.push_back(*it);
            }
        }
    }
}


void generate_particle_beyond_boundary_bins(){
    const int space = ceil(1.5 * bin_size * bin_size * 1. / density);
    for(int i = 0; i<column_lda; ++i){
        particle_beyond_upper_boundary_bins[i].clear();
        particle_beyond_upper_boundary_bins[i].reserve(space);
        particle_beyond_lower_boundary_bins[i].clear();
        particle_beyond_lower_boundary_bins[i].reserve(space);
    }
    // int sum = 0;
    for (auto it = particle_beyond_upper_boundary.begin(); it != particle_beyond_upper_boundary.end(); ++it){
        int index;
        double quotient = ((*it).x)/bin_size;
        index = int(quotient);
        particle_beyond_upper_boundary_bins[index].insert(&(*it));
    }
    // std::cout<<"sum"<<sum<<"\n";
    for (auto it = particle_beyond_lower_boundary.begin(); it != particle_beyond_lower_boundary.end(); ++it){
        int index;
        double quotient = ((*it).x)/bin_size;
        index = int(quotient);
        particle_beyond_lower_boundary_bins[index].insert(&(*it));
    }
}

void simulate_one_step(particle_t* parts, int num_parts, double size, int rank, int num_procs) {
    step += 1;
    generate_particle_beyond_boundary_bins();
    // Write this function
    particle_possible_coming_in_lower.clear();
    particle_possible_coming_in_upper.clear();
    // particle_going_out.clear();
    // particle_going_out_num = 0;
    particle_going_out_lower.clear();
    particle_going_out_lower_num = 0;
    particle_going_out_upper.clear();
    particle_going_out_upper_num = 0;


    for( int i =0 ; i < row_lda; i++){
        for( int j =0; j < column_lda; j++){
            for (auto it = bins[j+i*column_lda].begin(); it != bins[j+i*column_lda].end(); ++it){
                particle_t* ptr = *it;
                ptr->ax = 0;
                ptr->ay = 0;
            }
        }
    }

    for( int i =1 ; i < row_lda - 1; i++){
        for( int j =0; j < column_lda; j++){
            apply_force_bins(i, j);
        }
    }

    for( int j =0; j < column_lda; j++){
        apply_force_bins_upper_boundary(0,j);
        apply_force_bins_lower_boundary(row_lda-1,j);
    }


    for( int j =0; j < column_lda; j++){
        apply_force_bin_upper_boundary(1, j, 1, j);
        apply_force_bin_lower_boundary(row_lda-2, j, row_lda-2, j);
    }

    // for (int i = 0; i < num_parts; ++i) {
    //     move(parts[i], size);
    // }

    int going_out_space = ceil(1.5 * bin_size * zone_size * 1. / density);
    particle_going_out.reserve(going_out_space);
    update_flatten_particles();
    for (auto it = flatten_particles.begin(); it != flatten_particles.end(); ++it) {
        move(**it, size);
    }
    // int going_out_space = ceil(1.5 * bin_size * zone_size * 1. / density);
    // particle_going_out.reserve(going_out_space);

    // for (auto it = flatten_particles.begin(); it != flatten_particles.end(); ++it) {
    //     int index;
    //     index = calculate_bin_number((*it)->x,(*it)->y,size,bin_size,row_lda,column_lda);
    //     if (index == -1){
    //         particle_going_out.push_back(**it);
    //         particle_going_out_num += 1;
    //     }
    // }

    flatten_particles.clear();

    if(rank != (num_procs -1)){
        MPI_Send(&particle_going_out_lower_num,
                 1,
                 MPI_INT,
                 rank+1,
                 0,
                 MPI_COMM_WORLD);
        // if (particle_going_out_lower_num != 0){
        //     std::cout << "p"<<rank<<" send " << particle_going_out_lower_num << " to " << rank+1 << " step " << step << "\n";
        // }
        // std::cout << "p"<<rank<<" send " << particle_going_out_lower_num << "to" << rank+1<< "\n";
    }
    if (rank != 0){
        MPI_Send(&particle_going_out_upper_num,
                 1,
                 MPI_INT,
                 rank-1,
                 0,
                 MPI_COMM_WORLD);
        // if (particle_going_out_upper_num != 0){
        //     std::cout << "p"<<rank<<" send " << particle_going_out_upper_num << " to" << rank-1  << " step " << step << "\n";
        // }
        // std::cout << "p"<<rank<<" send " << particle_going_out_upper_num << "to" << rank-1<< "\n";
    }


    if(rank != (num_procs -1)){
        MPI_Recv(&particle_possible_coming_in_lower_num,
                 1,
                 MPI_INT,
                 rank+1,
                 0,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        // if (particle_possible_coming_in_lower_num !=0){
        //     std::cout << "p"<<rank<<" recv " << particle_possible_coming_in_lower_num << " from " << rank+1 << " step " << step << "\n";
        // }
        // std::cout << "p"<<rank<<" recv " << particle_possible_coming_in_lower_num << "from" << rank+1<< "\n";
    }
    if(rank != 0){
        MPI_Recv(&particle_possible_coming_in_upper_num,
                 1,
                 MPI_INT,
                 rank-1,
                 0,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        // if(particle_possible_coming_in_upper_num != 0){
        //     std::cout << "p"<<rank<<" recv " << particle_possible_coming_in_upper_num << " from " << rank-1<< " step " << step << "\n";
        // }
        // std::cout << "p"<<rank<<" recv " << particle_possible_coming_in_upper_num << "from" << rank-1<< "\n";
    }

    if(rank != (num_procs -1)){
        MPI_Send(&particle_going_out_lower[0],
                 particle_going_out_lower_num,
                 PARTICLE,
                 rank+1,
                 1,
                 MPI_COMM_WORLD);
    }
    if (rank != 0){
        MPI_Send(&particle_going_out_upper[0],
                 particle_going_out_upper_num,
                 PARTICLE,
                 rank-1,
                 1,
                 MPI_COMM_WORLD);
    }

    if(rank != (num_procs -1)){
        particle_possible_coming_in_lower.resize(particle_possible_coming_in_lower_num);
        MPI_Recv(&particle_possible_coming_in_lower[0],
                 particle_possible_coming_in_lower_num,
                 PARTICLE,
                 rank+1,
                 1,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        // if (particle_possible_coming_in_lower.size()!=0){
        //     std::cout << "p"<<rank<<" recv " << particle_possible_coming_in_lower.size() << " from " << rank+1<< " step " << step << "\n";
        //     std::cout << "y" << particle_possible_coming_in_lower[0].y << " step " << step << "\n";
        // }
        // std::cout << "p"<<rank<<" recv " << particle_possible_coming_in_lower.size() << "from" << rank+1<< "\n";
    }
    if(rank != 0){
        particle_possible_coming_in_upper.resize(particle_possible_coming_in_upper_num);
        MPI_Recv(&particle_possible_coming_in_upper[0],
                 particle_possible_coming_in_upper_num,
                 PARTICLE,
                 rank-1,
                 1,
                 MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        // if (particle_possible_coming_in_upper.size() != 0){
        //     std::cout << "p"<<rank<<" recv " << particle_possible_coming_in_upper.size() << " from " << rank-1 << " step " << step << "\n";
        //     std::cout << "y " << particle_possible_coming_in_upper[0].y << " step " << step << "\n";
        // }
        // std::cout << "p"<<rank<<" recv " << particle_possible_coming_in_upper.size() << "from" << rank-1 << "\n";
        // std::cout << " --------------------------------------" << "\n";
    }

    // bins.clear();

    // for( int i =0 ; i < row_lda; i++){
    //     for( int j =0; j < column_lda; j++){
    //         for (auto it = bins[j+i*column_lda].begin(); it != bins[j+i*column_lda].end(); ++it){
    //             int index;
    //             index = calculate_bin_number((*it)->x,(*it)->y, size, bin_size, row_lda, column_lda);
    //             if(index != 1){
    //                 bins[index].insert(*it);
    //             }
    //         }
    //     }
    // }


    for (auto it = particle_possible_coming_in_upper.begin(); it != particle_possible_coming_in_upper.end(); ++it) {
        particle_t* temp = new particle_t;
        temp -> id = (*it).id;
        temp -> x = (*it).x;
        temp -> y = (*it).y;
        temp -> vx = (*it).vx;
        temp -> vy = (*it).vy;
        temp -> ax = (*it).ax;
        temp -> ay = (*it).ay;
        int index;
        index = calculate_bin_number(temp->x,temp->y, size, bin_size, row_lda, column_lda);
        // std::cout << "rank" << rank << "index "<< index << " step " << step << "\n";
        if (index >= 0){
            bins[index].insert(temp);
        }
    }

    for (auto it = particle_possible_coming_in_lower.begin(); it != particle_possible_coming_in_lower.end(); ++it) {
        particle_t* temp = new particle_t;
        temp -> id = (*it).id;
        temp -> x = (*it).x;
        temp -> y = (*it).y;
        temp -> vx = (*it).vx;
        temp -> vy = (*it).vy;
        temp -> ax = (*it).ax;
        temp -> ay = (*it).ay;
        int index;
        index = calculate_bin_number(temp->x,temp->y, size, bin_size, row_lda, column_lda);
        // std::cout <<  "rank" <<rank <<"index "<< index << " step " << step << "\n";
        if (index >= 0){
            bins[index].insert(temp);
        }
    }
    // int temp_sum = 0;
    // for( int i =0 ; i < row_lda; i++){
    //     for( int j =0; j < column_lda; j++){
    //         temp_sum += bins[j+i*column_lda].size();
    //     }
    // }
    // std::cout <<  "rank" <<rank <<"tempsum "<< temp_sum << " step " << step << "\n";



    particle_possible_coming_in_lower.clear();
    particle_possible_coming_in_lower_num = 0;
    particle_possible_coming_in_upper.clear();
    particle_possible_coming_in_upper_num = 0;
    update_boundary_particles();
    send_recv_boundary_particles(rank, num_procs);
    // int sum = 0;
    // for( int i =0 ; i < row_lda; i++){
    //     for( int j =0; j < column_lda; j++){
    //         sum += bins[j+i*column_lda].size();
    //     }
    // }

}

void gather_for_save(particle_t* parts, int num_parts, double size, int rank, int num_procs) {
    // Write this function such that at the end of it, the master (rank == 0)
    // processor has an in-order view of all particles. That is, the array
    // parts is complete and sorted by particle id.

    // std::cout << "123" << "\n";
    number_particles_sending = 0;
    particle_send_gathering.clear();
    for( int i =0 ; i < row_lda; i++){
        for( int j =0; j < column_lda; j++){
            for (auto it = bins[j+i*column_lda].begin(); it != bins[j+i*column_lda].end(); ++it){
                // particle_t* ptr = *it;
                particle_send_gathering.push_back(**it);
                number_particles_sending+=1;
            }
        }
    }
    // std::cout <<"p "<< rank <<" " << number_particles_sending << " step "<< step << "\n";
    // std::cout << number_particles_sending << "\n";
    number_particles_receiving.clear();
    number_particles_receiving.reserve(num_procs);
    if(rank ==0){
        MPI_Gather(&number_particles_sending,
                   1,
                   MPI_INT,
                   &number_particles_receiving[0],
                   1,
                   MPI_INT,
                   0,
                   MPI_COMM_WORLD);

    }else{
        MPI_Gather(&number_particles_sending,
                   1,
                   MPI_INT,
                   nullptr,
                   0,
                   MPI_INT,
                   0,
                   MPI_COMM_WORLD);
        // std::cout << number_particles_sending << "\n";
    }

    // if(rank == 0){
    //     number_particles_receiving_sum = 0;
    //     int displacement[num_procs];
    //     for(int i = 0; i < num_procs; i++){
    //         displacement[i] = number_particles_receiving_sum;
    //         number_particles_receiving_sum = number_particles_receiving_sum + number_particles_receiving[i];
    //     }
    //     std::cout << number_particles_receiving_sum << "\n";
    // }

    if (rank == 0){
        particle_receive_gathering.clear();
        particle_receive_gathering.resize(num_parts);
        number_particles_receiving_sum = 0;
        displacement.clear();
        displacement.resize(num_procs);
        for(int i = 0; i < num_procs; i++){
            displacement[i] = number_particles_receiving_sum;
            number_particles_receiving_sum = number_particles_receiving_sum + number_particles_receiving[i];
        }
        std::cout << number_particles_receiving_sum << " step "<< step << "\n";

        MPI_Gatherv(&particle_send_gathering[0],
                    number_particles_sending,
                    PARTICLE,
                    &particle_receive_gathering[0],
                    &number_particles_receiving[0],
                    &displacement[0],
                    PARTICLE,
                    0,
                    MPI_COMM_WORLD);
        for( int i =0; i < num_parts; i++){
            parts[particle_receive_gathering[i].id-1] = particle_receive_gathering[i];
        }
    }
    else{
        MPI_Gatherv(&particle_send_gathering[0],
                    number_particles_sending,
                    PARTICLE,
                    nullptr,
                    nullptr,
                    nullptr,
                    PARTICLE,
                    0,
                    MPI_COMM_WORLD);
    }

}