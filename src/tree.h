
#pragma once
#include <mpi.h>
#include <unordered_map>

#include "allvars.h"
#include "hilbert.h"

struct topTree {

  static const int THRES = 2500;
  int my_rank, total_ranks;
  std::vector<uint64_t> leaves_to_add;
  std::vector<uint64_t> let;
  std::unordered_map<uint64_t, int> hash_counter;

  Domain domain;

  topTree()
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &total_ranks);

    if (my_rank == 0)
      {
        std::cout << "Constructing a top tree \n";
      }
  }

  ~topTree() 
  {
    leaves_to_add.clear();
    let.clear();
    hash_counter.clear();
  }

  void setDomain(double xmin, double xmax, 
		 double ymin, double ymax, 
		 double zmin, double zmax) 
  {
    domain.Xmin = xmin;
    domain.Xmax = xmax;
    domain.Ymin = ymin;
    domain.Ymax = ymax;
    domain.Zmin = zmin;    
    domain.Zmax = zmax;

    domain.X0[0]= (xmin + xmax)/2;
    double R0   = (xmax-xmin)/2;
    domain.X0[1]= (ymin + ymax)/2;
    R0 = std::max( (ymax-ymin)/2, R0);
    domain.X0[2]= (zmin + zmax)/2;
    R0 = std::max( (zmax-zmin)/2, R0);
    domain.R0 = R0 * 1.001;

    if (my_processor.my_rank == 0)
    {
        std::cout << ">>>>>> X0: [" << domain.X0[0] <<", " 
        << domain.X0[1] << ", " << domain.X0[2] << "]\n>>>>>> R0:" << domain.R0 << std::endl << std::endl;
    }
  }

  void buildTopTree(Particles& particles)
  {
    using namespace hilbert;

    double x_min =  MAXDOUBLE, y_min =  MAXDOUBLE, z_min =  MAXDOUBLE;
    double x_max = -MAXDOUBLE, y_max = -MAXDOUBLE, z_max = -MAXDOUBLE;

    for (auto &p : particles) {
      x_min = std::min(p.X[0], x_min);
      x_max = std::max(p.X[0], x_max);
      y_min = std::min(p.X[1], y_min);
      y_max = std::max(p.X[1], y_max);
      z_min = std::min(p.X[2], z_min);
      z_max = std::max(p.X[2], z_max);
    }

    double x_min_global, y_min_global, z_min_global;
    MPI_Allreduce(&x_min, &x_min_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&y_min, &y_min_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&z_min, &z_min_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    double x_max_global, y_max_global, z_max_global;
    MPI_Allreduce(&x_max, &x_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&y_max, &y_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&z_max, &z_max_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    /*if (my_processor.my_rank == 0)
      {
        std::cout << "Global Domain: \n";
        std::cout << ">>>>>> xmin: " << x_min_global << " xmax: " << x_max_global << std::endl;
        std::cout << ">>>>>> ymin: " << y_min_global << " ymax: " << y_max_global << std::endl;
        std::cout << ">>>>>> zmin: " << z_min_global << " zmax: " << z_max_global << std::endl;
      }
    */

    setDomain(x_min_global, x_max_global, y_min_global, y_max_global, z_min_global, z_max_global);

    hilbert::setDomain(domain.X0, domain.R0);

    for (auto&p :particles)
      {
	std::array<double, 3> coords{p.X[0], p.X[1], p.X[2]};
	p.hp_key= hilbert::getKey(coords);
	p.done=false;
	if (my_processor.my_rank == -1)
	  {
	    std::cout << p.hp_key <<" ==> "<< std::bitset<64>(p.hp_key)<< std::endl;
	  }
      }  
    std::vector<uint64_t> leaves_to_process;

    for(uint64_t j=9; j<73; j++)
      {
	leaves_to_process.push_back(j);
      }

    auto current_level = getLevel(leaves_to_process[0]);

    //std::cout << "Starting from level: " << current_level << std::endl;

    double start_time = MPI_Wtime();

    do{

      // set the hash_counter to each desired leaf cell to be zero
      hash_counter.clear();

      for(size_t j=0; j<leaves_to_process.size(); j++)
	{
	  hash_counter[leaves_to_process[j]] = 0;
	}

      for(auto&p: particles)
	{
	  if (!p.done) 
	    {
	      auto ans_id = getAncestor(p.hp_key, current_level);

	      if (hash_counter.find(ans_id) != hash_counter.end())
		{
		  hash_counter[ans_id]++;
		  //point parent top-node to be the current leave
		  p.top_node = ans_id;            
		}    
	    }
	}

      //this is to store the counts in the root_leaves
      std::vector<uint64_t> counts(leaves_to_process.size());
      for (size_t j=0; j<leaves_to_process.size(); j++)
	{
	  counts[j] = hash_counter[leaves_to_process[j]];
	}

      //find out the global sum 
      std::vector<uint64_t> counts_global(leaves_to_process.size(), 0);
      MPI_Allreduce(counts.data(), counts_global.data(), (int)counts.size(), 
		    MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

      //clear the buffer 
      leaves_to_add.clear();

      for (size_t i=0; i<counts_global.size(); i++)
	{
	  //update the hash_counter to be consistent with a global one
	  hash_counter[leaves_to_process[i]] = counts_global[i];

	  //push the children for large nodes
	  if (counts_global[i] > THRES)
	    {
	      for (int j=0; j<8; j++)
		{
		  leaves_to_add.push_back(getChild(leaves_to_process[i]) + j);
		}
	    }

	  //And only push the non-empty cells to locally essentil tree (let)
	  if (counts_global[i] > 0)
	    {
	      let.push_back(leaves_to_process[i]);
	    }
	}

      for (auto&p : particles)
        {
	  if (!p.done)
            {
	      if (hash_counter[p.top_node] <= THRES)
                {
		  p.done = true;
                }
            }
        }

      //double end_time = MPI_Wtime();
      //std::cout << "LET took " << end_time - start_time << " sec for level " << current_level << std::endl;

      uint64_t leaves_size = leaves_to_add.size();
      uint64_t leaves_tot_size;

      MPI_Allreduce(&leaves_size, &leaves_tot_size, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

      if (leaves_tot_size > 0)
	{
	  leaves_to_process = leaves_to_add;
	  current_level++;
	}
      else
	{
	  //no more leaves to process
	  break;
	}
  
    } while (true);
  
    double end_time = MPI_Wtime();

    if (my_rank == 0)
      {
        std::cout << "LET took " << end_time - start_time << " sec\n";
      }

  }

};
