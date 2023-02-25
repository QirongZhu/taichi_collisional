
#pragma once
#include <mpi.h>
#include <unordered_map>

#include "allvars.h"
#include "hilbert.h"

const int THRES = 1000;

void buildTopTree(Particles& particles)
{
  using namespace hilbert;

  std::vector<uint64_t> leaves_to_process;
  for(uint64_t j=9; j<73; j++)
  {
    leaves_to_process.push_back(j);
  }

  std::vector<uint64_t> leaves_to_add;
  std::vector<uint64_t> let;
  std::unordered_map<uint64_t, int> hash_counter;

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

  std::cout << "LET took " << end_time - start_time << " sec\n";

  if(my_processor.my_rank < 0)
  {
    for (auto& c: let)
    {
        std::cout << " " << c;
    }
    std::cout << std::endl;
  }
}