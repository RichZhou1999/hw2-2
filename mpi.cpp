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
double zone_size;
double bin_size = 2*cutoff;
int row_lda;
int column_lda;
double upper_boundary;
double lower_boundary;


// Vector of messages, each index is  a direction
// For row layout, 0 is upper direction, 1 is lower
// Still declare these containers statically

message_container_t upper_container;
message_container_t lower_container;
std::vector<message_container_t> message_containers;

// Particle In/Out Vectors
int particle_going_out_upper_num;
std::vector<particle_t> particle_going_out_upper;
int particle_going_out_lower_num;
std::vector<particle_t> particle_going_out_lower;
int particle_coming_in_upper_num;
std::vector<particle_t> particle_coming_in_upper;
int particle_coming_in_lower_num;
std::vector<particle_t> particle_coming_in_lower;

// // Ghost Particle Structures
//int num_ghost_particles_from_upper;
//std::vector<particle_t> ghost_particles_from_upper;
std::unordered_map<int, std::unordered_set<particle_t *>> ghost_particles_upper_bins;

// int num_ghost_particles_from_lower;
// std::vector<particle_t> ghost_particles_from_lower;
std::unordered_map<int, std::unordered_set<particle_t *>> ghost_particles_lower_bins;

// int num_ghost_particles_to_upper;
// std::vector<particle_t> ghost_particles_to_upper;
// int num_ghost_particles_to_lower;
// std::vector<particle_t> ghost_particles_to_lower;

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


void print_map(std::unordered_map<int, std::unordered_set<particle_t *>> const &m)
{
    for (auto const &pair: m) 
	{
        std::cout << "Bin : " << pair.first << "\n";
			for(auto const &parts : pair.second)
			{
			std::cout << "p_id : " << parts->id << " ";
			}
			std::cout << "\n";
	}
}

int calculate_bin_number(double x, double y, double size, double bin_size, int row_lda, int column_lda){
    if ((y < upper_boundary)){
        // Upper neighbor
        return -2;
    }
    if (( y >= lower_boundary)){
        // Lower neighbor
        return -1;
    }
    y = y - upper_boundary;
    double quotient;
    quotient = y/bin_size;
    int row = (int) quotient;

    quotient = x/bin_size;
    int column = (int) quotient;

    return row* column_lda + column;

}

// Helper function that works if particle is in bin
int calculate_if_edge(double y, double size){
    // 0 is upper edge
    // 1 is lower edge
    if ((y - upper_boundary) < cutoff){
        return 0; // In upper edge
    }
    if ((lower_boundary - y) <= cutoff){
        return 1; // In lower edge
    }
    return -1;
}


void update_boundary_particles(double size){
    // reset container for ghost particles
    for (auto container: message_containers) {
        container.ghost_particles_going_out.clear();
        container.num_ghost_particles_going_out = 0;
    }
    //ghost_particles_going_out.clear(); - upper/lower
    //num_ghost_particles_going_out = 0; - upper/lower
    // push back the particles within interaction range of particles in the adjacent cells
    int rows_to_check[] = {0, row_lda - 2, row_lda - 1};
    for (auto row: rows_to_check) {
        for(int i = 0; i < column_lda; i++){
            for (auto it = bins[i + row*column_lda].begin(); it != bins[i].end(); ++it){
                particle_t p = **it;
                int ind = calculate_if_edge(p.y, size);
                if (ind >= 0) {
                    auto container = message_containers.at(ind);
                    container.ghost_particles_going_out.push_back(p);
                    container.num_ghost_particles_going_out++;
		    
		    //ghost_particles_going_out.push_back(p);
		    //num_ghost_particles_going_out++;
                }
            }
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
        return;
    }
    
    //bins[origin_bin].erase(&p);
    /*
    // new_bin == -1, going out from upper boundary
    if (new_bin == -1){
        particle_going_out_upper.push_back(p);
        particle_going_out_upper_num += 1;
    }
    // new_bin == -2, going out from lower boundary
    if (new_bin == -2){
        particle_going_out_lower.push_back(p);
        particle_going_out_lower_num += 1;
    }
    */

    // Convert bin # (which neighbor) to index
    
    auto container = message_containers.at(new_bin + 2);
    container.particle_going_out.push_back(p);
    container.particle_going_out_num++;
    
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

void apply_force_bin_upper_boundary(int row, int column, int row2, int column2){
    // Force every particle in each bin to interact
    if (!check_boundary(row, column) && !check_boundary(row2, column2)){
        return;
    }
    // 
    if(column != 0){
        for (auto it2 = ghost_particles_upper_bins[column-1].begin(); it2 != ghost_particles_upper_bins[column-1].end(); ++it2){
            for (auto it = bins[column+row*column_lda].begin(); it != bins[column+row*column_lda].end(); ++it){
                // Interact particles
                apply_force_bi_direction(**it,**it2);
            }
        }
    }
    //
    for (auto it2 = ghost_particles_upper_bins[column].begin(); it2 != ghost_particles_upper_bins[column].end(); ++it2){
        for (auto it = bins[column+row*column_lda].begin(); it != bins[column+row*column_lda].end(); ++it){
            // Interact particles
            apply_force_bi_direction(**it,**it2);
        }
    }
    //
    if(column != (column_lda -1)){
        for (auto it2 = ghost_particles_upper_bins[column+1].begin(); it2 != ghost_particles_upper_bins[column+1].end(); ++it2){
            for (auto it = bins[column+row*column_lda].begin(); it != bins[column+row*column_lda].end(); ++it){
                // Interact particles
                apply_force_bi_direction(**it,**it2);
            }
        }
    }

}

void apply_force_bin_lower_boundary(int row, int column, int row2, int column2){
    // Force every particle in each bin to interact
    if (!check_boundary(row, column) && !check_boundary(row2, column2)){
        return;
    }
    //
    if (column != 0){
        for (auto it2 = ghost_particles_lower_bins[column-1].begin(); it2 != ghost_particles_lower_bins[column-1].end(); ++it2){
            for (auto it = bins[column+row*column_lda].begin(); it != bins[column+row*column_lda].end(); ++it){
                // Interact particles
                apply_force_bi_direction(**it,**it2);
            }
        }
    }
    //
    for (auto it2 = ghost_particles_lower_bins[column].begin(); it2 != ghost_particles_lower_bins[column].end(); ++it2){
        for (auto it = bins[column+row*column_lda].begin(); it != bins[column+row*column_lda].end(); ++it){
            // Interact particles
            apply_force_bi_direction(**it,**it2);
        }
    }
    //
    if (column != (column_lda -1)){
        for (auto it2 = ghost_particles_lower_bins[column+1].begin(); it2 != ghost_particles_lower_bins[column+1].end(); ++it2){
            for (auto it = bins[column+row*column_lda].begin(); it != bins[column+row*column_lda].end(); ++it){
                // Interact particles
                apply_force_bi_direction(**it,**it2);
            }
        }
    }

}

bool check_if_ranks_fit(int rank, int num_procs) {
    return rank >= 0 && rank < num_procs;
}

void recv_ghost_particles(int rank, int num_procs, int container_ind, int rank_diff) {
    auto container = message_containers.at(container_ind);
    auto other_rank = rank + rank_diff;
    if (!check_if_ranks_fit(other_rank, num_procs)) return;
    MPI_Recv(&container.num_ghost_particles_coming_in, // num_ghost_particles_coming_in,
                1,
                MPI_INT,
                other_rank,
                0,
                MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    //std::cout<<"zone: "<< rank << "received lower" << particle_beyond_lower_boundary_num<< " paticles" <<"\n";
    //ghost_particles_coming_in.resize(container.num_ghost_particles_coming_in); -upper/lower
    container.ghost_particles_coming_in.resize(container.num_ghost_particles_coming_in);
    MPI_Recv(&container.ghost_particles_coming_in[0], // ghost_particles_coming_in[0],
                container.num_ghost_particles_coming_in, //num_ghost_particles_coming_in,
                PARTICLE,
                other_rank,
                0,
                MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
}

void send_ghost_particles(int rank, int num_procs, int container_ind, int rank_diff) {
    auto container = message_containers.at(container_ind);
    int other_rank = rank + rank_diff;

    if (!check_if_ranks_fit(other_rank, num_procs)) return;
    MPI_Send(&container.num_ghost_particles_going_out, // &num_ghost_particles_going_out - upper/lower,
                1,
                MPI_INT,
                other_rank,
                0,
                MPI_COMM_WORLD);
    MPI_Send(&container.ghost_particles_going_out[0], // ghost_particles_going_out[0], - upper/lower
                container.num_ghost_particles_going_out, // num_ghost_particles_going_out
                PARTICLE,
                other_rank,
                0,
                MPI_COMM_WORLD);
}

void send_recv_ghost_particles(int rank , int num_procs){
    // Send upper info to upper neighbor
    send_ghost_particles(rank, num_procs, 0, -1);
    // Receive lower info from lower neighbor
    recv_ghost_particles(rank, num_procs, 1, 1);
    // Send lower info to lower neighbor
    send_ghost_particles(rank, num_procs, 1, 1);
    // Receive upper info from upper neighbor
    recv_ghost_particles(rank, num_procs, 0, -1);
}


void send_particles(int rank, int num_procs, int container_ind, int rank_diff) {
    auto container = message_containers.at(container_ind);
    int other_rank = rank + rank_diff;
    if (!check_if_ranks_fit(other_rank, num_procs)) return;
    MPI_Send(&container.particle_going_out_num, //container.particle_going_out_num,
                1,
                MPI_INT,
                other_rank,
                0,
                MPI_COMM_WORLD);
    MPI_Send(&container.particle_going_out[0],    //particle_going_out[0], -upper/lower
                container.particle_going_out_num, //particle_going_out_num, -upper/lower
                PARTICLE,
                other_rank,
                1,
                MPI_COMM_WORLD);
	//std::cout << "Send - P : " << rank << " incoming " << container.particle_coming_in_num << " from " << other_rank << "\n";
}

void recv_particles(int rank, int num_procs, int container_ind, int rank_diff) {
    auto container = message_containers.at(container_ind);
    int other_rank = rank + rank_diff;
    if (!check_if_ranks_fit(other_rank, num_procs)) return;
    MPI_Recv(&container.particle_coming_in_num, //container.particle_coming_in_num,
                1,
                MPI_INT,
                other_rank,
                0,
                MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    container.particle_coming_in.resize(container.particle_coming_in_num);
    //particle_coming_in.resize(container.particle_coming_in_num);
    MPI_Recv(&container.particle_coming_in[0],    //particle_coming_in[0], -upper/lower
                container.particle_coming_in_num, //particle_coming_in_num, -upper/lower
                PARTICLE,
                other_rank,
                1,
                MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
	//std::cout << "Recv - P : " << rank << " incoming " << container.particle_coming_in_num << " from " << other_rank << "\n";
}

void send_recv_particles(int rank, int num_procs) {
    // Send upper info to upper neighbor
    send_particles(rank, num_procs, 0, -1);
    // Receive lower info from lower neighbor
    recv_particles(rank, num_procs, 1, 1);
    // Send lower info to lower neighbor
    send_particles(rank, num_procs, 1, 1);
    // Receive upper info from upper neighbor
    recv_particles(rank, num_procs, 0, -1);
}

void init_simulation(particle_t* parts, int num_parts, double size, int rank, int num_procs) {
    if ((size / num_procs/bin_size) < 1){
        bin_size = cutoff;
    }
    zone_size = size / num_procs;


    double quotient = zone_size / bin_size;
    row_lda = (int) ceil(quotient);

    quotient = size / bin_size;
    column_lda = (int) ceil(quotient);

    upper_boundary = zone_size * rank;
    lower_boundary = zone_size * (rank + 1);

    const int space = ceil(1.5 * bin_size * bin_size * 1. / density);
    for(int i = 0; i< row_lda*column_lda; ++i){
        bins[i].reserve(space);
    }

    for (int i = 0; i < num_parts; ++i){
        int index;
        index = calculate_bin_number(parts[i].x,parts[i].y, size, bin_size, row_lda, column_lda);
        if (index >= 0){
            bins[index].insert(&parts[i]);
        }
    }

    number_particles_sending = 0;
    for( int i =0 ; i < row_lda; i++){
        for( int j =0; j < column_lda; j++){
            for (auto it = bins[j+i*column_lda].begin(); it != bins[j+i*column_lda].end(); ++it){
                number_particles_sending+=1;
            }
        }
    }
    message_containers.push_back(upper_container);
    message_containers.push_back(lower_container);

    update_boundary_particles(size);
    send_recv_ghost_particles(rank, num_procs);
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
                particle_t* ptr = *it;
                flatten_particles.push_back(ptr);
            }
        }
    }
}
// update GZs
void generate_particle_beyond_boundary_bins(){
    const int space = ceil(1.5 * bin_size * bin_size * 1. / density);
    for(int i = 0; i<column_lda; ++i){
        ghost_particles_upper_bins[i].clear();
        ghost_particles_upper_bins[i].reserve(space);
        ghost_particles_lower_bins[i].clear();
        ghost_particles_lower_bins[i].reserve(space);
    }
    // int sum = 0;
    auto upper_container = message_containers.at(0);
    for (auto it = upper_container.ghost_particles_coming_in.begin(); it != upper_container.ghost_particles_coming_in.end(); ++it){
        int index;
        double quotient = ((*it).x)/bin_size;
        index = int(quotient);
        ghost_particles_upper_bins[index].insert(&(*it));
    }
    // std::cout<<"sum"<<sum<<"\n";
    auto lower_container = message_containers.at(1);
    for (auto it = lower_container.ghost_particles_coming_in.begin(); it != lower_container.ghost_particles_coming_in.end(); ++it){
        int index;
        double quotient = ((*it).x)/bin_size;
        index = int(quotient);
        ghost_particles_lower_bins[index].insert(&(*it));
    }
}


void simulate_one_step(particle_t* parts, int num_parts, double size, int rank, int num_procs) {
    step += 1;
    generate_particle_beyond_boundary_bins();
    // reset messages passing from boundaries  
    for (auto container: message_containers) {
        container.particle_coming_in.clear();
        container.particle_coming_in_num = 0;
        container.particle_going_out.clear();
        container.particle_going_out_num = 0;
    }

    /*
    particle_coming_in.clear();
    particle_coming_in_num = 0;
    particle_going_out.clear();
    particle_going_out_num = 0;
    */
    // reset acceleration to zero
	for( int i =0 ; i < row_lda; i++)
	{
        for( int j =0; j < column_lda; j++)
		{
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
    // particle_going_out.reserve(going_out_space);
    update_flatten_particles();
    for (auto it = flatten_particles.begin(); it != flatten_particles.end(); ++it) {
        move(**it, size);
    }

    flatten_particles.clear();

    // send_recv_particles(rank, num_procs);
    send_recv_particles(rank, num_procs);
    
    for (auto container: message_containers) {
        for (auto it = container.particle_coming_in.begin(); it != container.particle_coming_in.end(); ++it) {
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
        container.particle_coming_in.clear();
        container.particle_coming_in_num = 0;
    }
   
   /*
   // accept the newly incoming particles from upper direction and assign b$
   for (auto it = particle_possible_coming_in_upper.begin(); it != particle)
   {
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
        // std::cout << "rank" << rank << "index "<< index << " step " << st$
        if (index >= 0){
            bins[index].insert(temp);
        }
    }

    // accept the newly incoming particles from lower direction and assign b$
    for (auto it = particle_possible_coming_in_lower.begin(); it != particle)
    {
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
        // std::cout <<  "rank" <<rank <<"index "<< index << " step " << ste$
        if (index >= 0){
            bins[index].insert(temp);
    }
    // reset the container for incoming particles
    particle_possible_coming_in_lower.clear();
    particle_possible_coming_in_lower_num = 0;
    particle_possible_coming_in_upper.clear();
    particle_possible_coming_in_upper_num = 0;
    */

    // update the boundary particles
    update_boundary_particles(size);
    // send/recv the ghost particles
    send_recv_ghost_particles(rank, num_procs);

}

void gather_for_save(particle_t* parts, int num_parts, double size, int rank, int num_procs) {
    // Write this function such that at the end of it, the master (rank == 0)
    // processor has an in-order view of all particles. That is, the array
    // parts is complete and sorted by particle id.

    // std::cout << "123" << "\n";
    number_particles_sending = 0;
    particle_send_gathering.clear();
    //print_map(bins);
    std::cout << " Sizeof bins : " << sizeof(bins) << "\n";
    for( int i =0 ; i < row_lda; i++)
    {
        for( int j =0; j < column_lda; j++)
	{
            for (auto it = bins[j+i*column_lda].begin(); it != bins[j+i*column_lda].end(); ++it)
	    {
                particle_t p = **it;
                particle_send_gathering.push_back(p);
                //std::cout << " Gathered pid : " << (**it).id << "\n";
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
            number_particles_receiving_sum += number_particles_receiving[i];
        }
        std::cout << "--------------------------------\n\n # Particles : " << number_particles_receiving_sum << " step "<< step <<"\n ----------------------------- \n\n";
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
            // set particles[pid-1] = pid (sorting)
	    parts[particle_receive_gathering[i].id-1] = particle_receive_gathering[i];
        }
	/*for( int i = 0; i < num_parts; ++i){
	  std::cout << parts[i].id << " ";
	  if( !( (i + 1) % 100)){std::cout << "\n";}
	}
	std::cout << "\n";*/
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
