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
// Map | particle idx : particle
std::unordered_map<int, std::unordered_set<particle_t *>> bins;
// Map | particle idx : bin idx
std::unordered_map<int, int> particle_to_bin;
// Program constants
double zone_size;
double bin_size = cutoff;
// nr of rows in domain
int    irow_lda;
// nr of cols in domain
int    column_lda;
// y-value of upper and lower boundary
double upper_boundary;
double lower_boundary;

// data containers for boundary interactions
// outgoing particles
int particle_going_out_num;
std::vector<particle_t> particle_going_out;
// outgoing particle from upper boundary
int particle_going_out_upper_num;
std::vector<particle_t> particle_going_out_upper;
// outgoing particles from lower boundary
int particle_going_out_lower_num;
std::vector<particle_t> particle_going_out_lower;
// incoming particles upper boundary
int particle_possible_coming_in_upper_num;
std::vector<particle_t> particle_possible_coming_in_upper;
// incoming particles lower boundary
int particle_possible_coming_in_lower_num;
std::vector<particle_t> particle_possible_coming_in_lower;
// particles gone beyond upper boundary - Upper Ghost Zone (GZ) particles
int particle_beyond_upper_boundary_num;
std::vector<particle_t> particle_beyond_upper_boundary;
// particles gone beyond lower boundary - Lower GZ particles
int particle_beyond_lower_boundary_num;
std::vector<particle_t> particle_beyond_lower_boundary;
// Upper GZ particles
int upper_boundary_particle_num;
std::vector<particle_t> particle_upper_boundary;
// Lower GZ particles
int lower_boudary_particle_num;
std::vector<particle_t> particle_lower_boundary;
// ?
std::vector<particle_t *> flatten_particles;
// data container for gather particles to rank 0 processorgc
// before correctness check
std::vector<particle_t > particle_send_gathering;
// nr partaicles sent from all P with rank != 0
int number_particles_sending;
// nr particles received at P with rank=0gc
// from all P with rank != 0
int number_particles_receiving_sum;
// data container for particle reception at P_rank0
std::vector<particle_t > particle_receive_gathering;
// ?gc
std::vector<int > displacement;
// vector to store number of particles from all bins not contained in P_rank0 ?
std::vector<int> number_particles_receiving;
int step = 0;

int calculate_bin_number(double x, double y, double size, double bin_size, int row_lda, int column_lda){
    if ((y < upper_boundary)){
        return -1;
    }
    if (( y >= lower_boundary)){
        return -2;
    }
    // ?
	y = y - upper_boundary;
    // row index
	int row = int( y / bin_size );
    // col index
	int column = int( x / bin_size );
    return row* column_lda + column;

}

void update_boundary_particles(){
    // reset the boundary containers and particle numbers?
	particle_upper_boundary.clear();
    upper_boundary_particle_num = 0;
    particle_lower_boundary.clear();
    lower_boudary_particle_num = 0;
	// iterate through all the columns in the bins
    for(int i = 0; i < column_lda; i++){
		// iterate through all the rows within each bin
		// add particles to Upper GZ from current bin
        for (auto it = bins[i].begin(); it != bins[i].end(); ++it){
            particle_upper_boundary.push_back(**it);
            upper_boundary_particle_num += 1;
        }
		// add particles to Upper GZ from above bin
        for (auto it = bins[i + column_lda].begin(); it != bins[i+ column_lda].end(); ++it){
            particle_upper_boundary.push_back(**it);
            upper_boundary_particle_num += 1;
        }
		// add particles to Lower GZ from current bin
		// layer in current bin included in GZ of neigh is at i + (irow_lda-2)
        for (auto it = bins[i + (irow_lda-2)*column_lda].begin(); it != bins[i+ (irow_lda-2)*column_lda].end(); ++it){
            particle_lower_boundary.push_back(**it);
            lower_boudary_particle_num += 1;
        }
		// add particles to Lower GZ from below bingc
		// layer in current bin included in GZ of current
        for (auto it = bins[i + (irow_lda-1)*column_lda].begin(); it != bins[i+ (irow_lda-1)*column_lda].end(); ++it){
            particle_lower_boundary.push_back(**it);
            lower_boudary_particle_num += 1;
        }
    }
}


// unidirectional force function
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

// bidirectional force function (stores same but opposite force on neighbor)
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
}

// given row and col indices
// check if at boundary
bool check_boundary(int row , int column){
    if ((row < 0) or (row>=irow_lda)){
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

    // origin of the current bin (same as bin number)
	int origin_bin = calculate_bin_number(p.x,p.y, size,bin_size,irow_lda, column_lda);
    // update velocities
	p.vx += p.ax * dt;
    p.vy += p.ay * dt;
	// update pos
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
    // compute new bin indices from new particle position
    int new_bin = calculate_bin_number(p.x, p.y, size, bin_size, irow_lda, column_lda);
    if(origin_bin == new_bin) return;

    // move particle within processor or bin?
	if (new_bin >= 0 ){
        bins[origin_bin].erase(&p);
        bins[new_bin].insert(&p);
    }
	// move particle from current to upper processor or bin?
    if (new_bin == -1){
        bins[origin_bin].erase(&p);
        particle_going_out_upper.push_back(p);
        particle_going_out_upper_num += 1;
    }
	// move particle from current to lower processor or bin?
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

// why not bidirectional here?
void apply_force_bin_self(int row, int column, int row2, int column2){
    // Force every particle in each bin to interact
    if (!check_boundary(row, column) && !check_boundary(row2, column2)){
        return;
    }
    for (auto it2 = bins[column2+row2*column_lda].begin(); it2 != bins[column2+row2*column_lda].end(); ++it2){
        for (auto it = bins[column+row*column_lda].begin(); it != bins[column+row*column_lda].end(); ++it){
            // Interact particles
            //apply_force_bi_direction(**it, **it2);
			apply_force_one_direction(**it, **it2);
        }
    }
}

void apply_force_bin_upper_boundary(int row, int column, int row2, int column2){
    // Force every particle in each bin to interact
    if (!check_boundary(row, column) && !check_boundary(row2, column2)){
        return;
    }
    for (auto it2 = particle_beyond_upper_boundary.begin(); it2 != particle_beyond_upper_boundary.end(); ++it2){
        for (auto it = bins[column+row*column_lda].begin(); it != bins[column+row*column_lda].end(); ++it){
            // Interact particles
            apply_force_bi_direction(**it,*it2);
        }
    }
}


void apply_force_bin_lower_boundary(int row, int column, int row2, int column2){
    // Force every particle in each bin to interact
    if (!check_boundary(row, column) && !check_boundary(row2, column2)){
        return;
    }
    for (auto it2 = particle_beyond_lower_boundary.begin(); it2 != particle_beyond_lower_boundary.end(); ++it2){
        for (auto it = bins[column+row*column_lda].begin(); it != bins[column+row*column_lda].end(); ++it){
            // Interact particles
            apply_force_bi_direction(**it,*it2);
        }
    }
}

void send_recv_boundary_particles(int rank , int num_procs){
    /*
    MPI_RECV(buf, count, datatype, source, tag, comm, status)
	OUT buf initial address of receive buffer (choice)
	IN count number of elements in receive buffer (non-negative
	integer)
	IN datatype datatype of each receive buffer element (handle)
	IN source rank of source or MPI_ANY_SOURCE (integer)
	IN tag message tag or MPI_ANY_TAG (integer)
	IN comm communicator (handle)
	OUT status status object (status)
	*/
	if(rank != (num_procs -1))
	{
		// send particle number
        MPI_Send(&lower_boudary_particle_num, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
		// send particle from lower boundary
        MPI_Send(&particle_lower_boundary[0], lower_boudary_particle_num, PARTICLE, rank+1, 0, MPI_COMM_WORLD);
    }
    if(rank != 0)
	{
        // recv particle number 
		MPI_Recv(&particle_beyond_upper_boundary_num, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        particle_beyond_upper_boundary.resize(particle_beyond_upper_boundary_num);
        
		MPI_Recv(&particle_beyond_upper_boundary[0],  particle_beyond_upper_boundary_num, 
		         PARTICLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    }
    if(rank != 0)
	{
        MPI_Send(&upper_boundary_particle_num, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD);
        MPI_Send(&particle_upper_boundary[0], upper_boundary_particle_num, PARTICLE, rank-1, 0, MPI_COMM_WORLD);
    }
    if(rank != (num_procs -1))
	{
        MPI_Recv(&particle_beyond_lower_boundary_num, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        particle_beyond_lower_boundary.resize(particle_beyond_lower_boundary_num);
        MPI_Recv(&particle_beyond_lower_boundary[0], particle_beyond_lower_boundary_num, PARTICLE, 
		         rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

}
void init_simulation(particle_t* parts, int num_parts, double size, int rank, int num_procs) {
    zone_size = size / num_procs;
    // nr of bins in the GZ
    double quotient = zone_size / bin_size;
    // row index of the GZ
	irow_lda = (int) ceil(quotient);
    // number of bins in the domain
    quotient = size / bin_size;
	// nr col in domain
    column_lda = (int) ceil(quotient);
    // row index of Upper GZ
    upper_boundary = zone_size * rank;
    // row index of Lower GZ
	lower_boundary = zone_size * (rank + 1);
    // std::cout<<"upper_boundary" <<upper_boundary << "rank" << rank << "\n";
    // std::cout<<"lower_boundary" <<lower_boundary << "rank" << rank << "\n";
    // whats this? allocate space elements in bin
	const int space = ceil(1.5 * bin_size * bin_size * 1. / density);
    for(int i = 0; i< irow_lda*column_lda; ++i){
        bins[i].reserve(space);
    }
    // compute index
    for (int i = 0; i < num_parts; ++i){
        int index;
        index = calculate_bin_number(parts[i].x,parts[i].y, size, bin_size, irow_lda, column_lda);
        if (index >= 0){
            bins[index].insert(&parts[i]);
        }
    }

    number_particles_sending = 0;
    for( int i =0 ; i < irow_lda; i++){
        for( int j =0; j < column_lda; j++){
            for (auto it = bins[j+i*column_lda].begin(); it != bins[j+i*column_lda].end(); ++it){
                number_particles_sending+=1;
            }
        }
    }
    // std::cout << "rank "<<rank<< "num" << number_particles_sending << "\n";

    update_boundary_particles();
    send_recv_boundary_particles(rank, num_procs);
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
    for( int i =0 ; i < irow_lda; i++){
        for( int j =0; j < column_lda; j++){
            for (auto it = bins[j+i*column_lda].begin(); it != bins[j+i*column_lda].end(); ++it){
                // particle_t* ptr = *it;
                flatten_particles.push_back(*it);
            }
        }
    }
}


void simulate_one_step(particle_t* parts, int num_parts, double size, int rank, int num_procs) {
    step += 1;
    particle_possible_coming_in_lower.clear();
    particle_possible_coming_in_upper.clear();
    particle_going_out_lower.clear();
    particle_going_out_lower_num = 0;
    particle_going_out_upper.clear();
    particle_going_out_upper_num = 0;
    // reset particles acceleration
    for( int i =0 ; i < irow_lda; i++)
	{
        for( int j =0; j < column_lda; j++)
		{
            for (auto it = bins[j+i*column_lda].begin(); it != bins[j+i*column_lda].end(); ++it)
			{
                (*it)->ax = 0;
                (*it)->ay = 0;
            }
        }
    }

    for( int i =1 ; i < irow_lda - 1; i++)
	{
        for( int j =0; j < column_lda; j++)
		{
            apply_force_bins(i, j);
        }
    }

    for( int j =0; j < column_lda; j++)
	{
        apply_force_bins_upper_boundary(0,j);
        apply_force_bins_lower_boundary(irow_lda-1,j);
    }


    for( int j =0; j < column_lda; j++)
	{
        apply_force_bin_upper_boundary(1, j, 1, j);
        apply_force_bin_lower_boundary(irow_lda-2, j, irow_lda-2, j);
    }
    // nr outgoing particles
    int going_out_space = ceil(1.5 * bin_size * zone_size * 1. / density);
    particle_going_out.reserve(going_out_space);
    update_flatten_particles();
    for (auto it = flatten_particles.begin(); it != flatten_particles.end(); ++it)
	{
        move(**it, size);
    }
    flatten_particles.clear();
    if (rank != 0)
	{
        MPI_Send(&particle_going_out_upper_num, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD);
		MPI_Recv(&particle_possible_coming_in_upper_num, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Send(&particle_going_out_upper[0], particle_going_out_upper_num, PARTICLE, rank-1, 1, MPI_COMM_WORLD);
		particle_possible_coming_in_upper.resize(particle_possible_coming_in_upper_num);
        MPI_Recv(&particle_possible_coming_in_upper[0], particle_possible_coming_in_upper_num, PARTICLE, rank-1, 1,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
	if(rank != (num_procs -1))
	{
		MPI_Send(&particle_going_out_lower_num, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
        MPI_Send(&particle_going_out_lower[0], particle_going_out_lower_num, PARTICLE, rank+1, 1, MPI_COMM_WORLD);
		particle_possible_coming_in_lower.resize(particle_possible_coming_in_lower_num);
        MPI_Recv(&particle_possible_coming_in_lower[0], particle_possible_coming_in_lower_num, PARTICLE, rank+1, 1,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	for (auto it = particle_possible_coming_in_upper.begin(); it != particle_possible_coming_in_upper.end(); ++it) {
        particle_t temp = *it;
        int index;
        index = calculate_bin_number(temp.x,temp.y, size, bin_size, irow_lda, column_lda);
        if (index >= 0){
            bins[index].insert( &(*it) );
        }
    }

    for (auto it = particle_possible_coming_in_lower.begin(); it != particle_possible_coming_in_lower.end(); ++it) {
        particle_t temp = *it;        
		int index;
        index = calculate_bin_number(temp.x,temp.y, size, bin_size, irow_lda, column_lda);
        if (index >= 0){
            bins[index].insert( &(*it) );
        }
    }


    // clear and reset data containers for update on GZ
    particle_possible_coming_in_lower.clear();
    particle_possible_coming_in_lower_num = 0;
    particle_possible_coming_in_upper.clear();
    particle_possible_coming_in_upper_num = 0;
    update_boundary_particles();
    send_recv_boundary_particles(rank, num_procs);
}




// Write this function such that at the end of it, the master (rank == 0)
// processor has an in-order view of all particles. That is, the array
// parts is complete and sorted by particle id.
void gather_for_save(particle_t* parts, int num_parts, double size, int rank, int num_procs)
{
    number_particles_sending = 0;
    particle_send_gathering.clear();
    for( int i =0 ; i < irow_lda; i++){
        for( int j =0; j < column_lda; j++){
            for (auto it = bins[j+i*column_lda].begin(); it != bins[j+i*column_lda].end(); ++it){
                particle_send_gathering.push_back(**it);
                number_particles_sending+=1;
            }
        }
    }
    number_particles_receiving.clear();
    number_particles_receiving.reserve(num_procs);
    if(rank ==0){
        MPI_Gather(&number_particles_sending, 1, MPI_INT, &number_particles_receiving[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
	else{
		MPI_Gather(&number_particles_sending, 1, MPI_INT, nullptr,                        0, MPI_INT, 0, MPI_COMM_WORLD);
    }

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

        MPI_Gatherv(&particle_send_gathering[0], number_particles_sending, PARTICLE, &particle_receive_gathering[0],
                    &number_particles_receiving[0], &displacement[0],      PARTICLE, 0, MPI_COMM_WORLD);
        for( int i =0; i < num_parts; i++){
            parts[particle_receive_gathering[i].id-1] = particle_receive_gathering[i];
        }
    }
    else{
        MPI_Gatherv(&particle_send_gathering[0], number_particles_sending, PARTICLE, nullptr,
		            nullptr,                     nullptr,                  PARTICLE, 0, MPI_COMM_WORLD);
    }
}
