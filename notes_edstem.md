# Notes from EdStem forum

## Zhanghao Wu:

### Q: What's the conversion between processors and ranks and nodes? If we run with 2 codes and 64 ranks, does that mean we have 128 processors?
- When there are 2 nodes with 64 processes each, there will be 128 ranks.
- Q: I see, so which one corresponds to the p in the O(n/p) scaling relation? 
  A: p should be the total number of processors used (2*64=128 in the example), i.e. the number ranks. 

### Q: What is the point of void gatheri\_for\_save(particle\_t* parts, int num_parts, double size, int rank, int num_procs)? Since parts is never modified, why would it ever go out of order? Why sort it at all?
- A: The function is only used for correctness checks. As each processor is only in charge of the particles in the block assigned to it, the function is used for gathering the particles (with the position informations) from all the processors into the same processor. The gathered particles will not be in the order of particle.id. We need to reorder them and save them back to parts . 

### Q: How to create a communicator for every process within a node? I am thinking of dividing the processors into two communicators, one for each node. I want to batch all communication across the nodes into one big message, but allow point-to-point communication within each node, since I think it is a lot longer in latency to communicate across nodes than communication within a node.
- A: For simplicity, in this homework, it is fine to consider all the processors having the same topology, i.e. not separating the processors in the same node vs different nodes. If you still would like to split the two groups, Professor Grigori's lecture on "Advanced MPI and Collective Communication Algorithms" introduces MPI_Comm_split. I think you're still left with the problem of mapping ranks to nodes, and I think the MPI_Comm_split_type has a way to partition processes based on shared memory using the type MPI_COMM_TYPE_SHARED (see https://www.open-mpi.org/doc/v3.0/man3/MPI_Comm_split_type.3.php). Hopefully, this way of splitting groups is available on Perlmutter.- A Tianyu Liang: Take a look at this Inter-communicator Operations (mpi-forum.org) However, I wouldn't recommend this approach because it overcomplicates things.

### Q: If our implementation change the ordering the particles, is that count as incorrect? 
- For example, if the particles are 1,2,3,4,5,...,10 and we split them into 2 procs, and when we gather them back, the particles are 1,4,6,9,10,2,3,5,7,8 but the positions are still correct. How can we measure correctness in comparison to verf.out?
- A. You need to reorder the particles based the particle.id . Our correctness checker will require the particles to be in the same order as the original parts array. The gather_for_save will not be called in the performance benchmark. Please feel free to reorder the particles in rank 0.

### Q: Can we get some performance target for HW2.2?
- A: We don't know the numbers from last year because this is our first time using Perlmutter. A well-implemented algorithm would have super linear strong scaling from 1 node to 2 nodes.  We ran some code from our GSIs, and the time for two nodes (64 processes each) on 6 million particles is from 9~20 seconds.

### Q: What's the reason that the main.cpp does this: MPI\_Bcast(parts, num\_parts, PARTICLE, 0, MPI\_COMM\_WORLD); 
- Isn't this meaning that we are broadcasting all particles to all MPI processors? This doesn't seems to be memory efficient...
- A: The primary purpose of that line is to distribute all the particles to all the processors, so in  init_simulation you don't have to worry about communication, as the parts array will contain all the particles on all the processors. It is indeed not memory efficient, but we trade it for your convenience to implement the init_simulation.
