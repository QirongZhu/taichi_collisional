#ifndef io_h
#define io_h
#include <numeric>      // std::iota
#include <algorithm>    // sort
#include "exafmm.h"
#include "hdf5.h"

namespace exafmm {

void read_snapshot(int snapshot_num) {

  static char dummyString[200];
  strncpy(dummyString, input_fname+9, strlen(input_fname)-5-9);  
  snapnum = atoi(dummyString);

  herr_t      status;
  hid_t file_id, hdf5_headergrp, hdf5_attribute, dataset;
        
  file_id = H5Fopen(input_fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  hdf5_headergrp = H5Gopen(file_id, "/Header");
        
  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "NumPartTotal");
  H5Aread(hdf5_attribute, H5T_NATIVE_UINT, &numBodies);
  H5Aclose(hdf5_attribute);
        
  hdf5_attribute = H5Aopen_name(hdf5_headergrp, "Time");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &t_now);
  H5Aclose(hdf5_attribute);
        
  H5Gclose(hdf5_headergrp);       
        
  mainsys.n   = numBodies;
  mainsys.part=(struct particle*) malloc(numBodies*sizeof(struct particle));
  mainsys.last=&mainsys.part[numBodies-1];        

  double *data = new double[numBodies];

  dataset = H5Dopen(file_id, "/Posx/");
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  status  = H5Dclose(dataset);
  for (unsigned int b=0; b<numBodies; b++) mainsys.part[b].pos[0] = data[b];
        
  dataset = H5Dopen(file_id, "/Posy/");
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  status  = H5Dclose(dataset);
  for (unsigned int b=0; b<numBodies; b++) mainsys.part[b].pos[1] = data[b];
        
  dataset = H5Dopen(file_id, "/Posz/");
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  status  = H5Dclose(dataset);
  for (unsigned int b=0; b<numBodies; b++) mainsys.part[b].pos[2] = data[b];
        
  dataset = H5Dopen(file_id, "/Velx/");
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  status  = H5Dclose(dataset);
  for (unsigned int b=0; b<numBodies; b++) mainsys.part[b].vel[0] = data[b];
        
  dataset = H5Dopen(file_id, "/Vely/");
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  status  = H5Dclose(dataset);
  for (unsigned int b=0; b<numBodies; b++) mainsys.part[b].vel[1] = data[b];
        
  dataset = H5Dopen(file_id, "/Velz/");
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  status  = H5Dclose(dataset);
  for (unsigned int b=0; b<numBodies; b++) mainsys.part[b].vel[2] = data[b];
        
  dataset = H5Dopen(file_id, "/Mass/");
  status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  status  = H5Dclose(dataset);
  for (unsigned int b=0; b<numBodies; b++) mainsys.part[b].mass = data[b];
        
  delete[] data;

  unsigned int *id = new unsigned int[numBodies];
  dataset = H5Dopen(file_id, "/ID/");
  status  = H5Dread(dataset, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
  status  = H5Dclose(dataset);
  for (unsigned int b=0; b<numBodies; b++) mainsys.part[b].id = id[b];
  delete[] id;

  status  = H5Fclose(file_id);

  if(status < 0) printf("Reading snapshot error\n");
      
  printf("Reading IC done with %d particles \n", numBodies);
  fflush(stdout); 
}


void write_snapshot(int snapshot_num, struct sys s) {
  hsize_t     dims[1] = {s.n};
  herr_t      status;
  hid_t       file_id, space_id, dset_id, memspace, handle=0;
  char buf[500];
  sprintf(buf, "snapshot_%d.hdf5", snapshot_num);

  double *data = new double[s.n];

  file_id = H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // first write the header information
  handle = H5Gcreate(file_id, "/Header", 0);

  hid_t hdf5_dataspace, hdf5_attribute;

  printf("Writing %d particles at t=%g \n", numBodies, t_now);
  fflush(stdout);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "NumPartTotal", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, &numBodies);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);

  hdf5_dataspace = H5Screate(H5S_SCALAR);
  hdf5_attribute = H5Acreate(handle, "Time", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &t_now);
  H5Aclose(hdf5_attribute);
  H5Sclose(hdf5_dataspace);
  //Header writing done 

  std::vector<unsigned int> ids(s.n);
  std::iota(ids.begin(), ids.end(), 0);
  sort(ids.begin(), ids.end(),
    [&s](unsigned int i1, unsigned int i2) {return s.part[i1].id<s.part[i2].id;});

  memspace= H5Screate_simple(1, dims, NULL);
  space_id= H5Screate_simple(1, dims, NULL);
  for (unsigned int b=0; b<s.n; b++)
    data[b] = s.part[ids[b]].pos[0];
  dset_id = H5Dcreate(file_id, "Posx", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
  status  = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
  status  = H5Sclose(space_id);
  status  = H5Sclose(memspace);
  status  = H5Dclose(dset_id);

  memspace= H5Screate_simple(1, dims, NULL);
  space_id= H5Screate_simple(1, dims, NULL);
  for (unsigned int b=0; b<s.n; b++)
    data[b] = s.part[ids[b]].pos[1];
  dset_id = H5Dcreate(file_id, "Posy", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
  status  = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
  status  = H5Sclose(space_id);
  status  = H5Sclose(memspace);
  status  = H5Dclose(dset_id);

  memspace= H5Screate_simple(1, dims, NULL);
  space_id= H5Screate_simple(1, dims, NULL);
  for (unsigned int b=0; b<s.n; b++)
    data[b] = s.part[ids[b]].pos[2];
  dset_id = H5Dcreate(file_id, "Posz", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
  status  = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
  status  = H5Sclose(space_id);
  status  = H5Sclose(memspace);
  status  = H5Dclose(dset_id);

  memspace= H5Screate_simple(1, dims, NULL);
  space_id= H5Screate_simple(1, dims, NULL);
  for (unsigned int b=0; b<s.n; b++)
    data[b] = s.part[ids[b]].vel[0];
  dset_id = H5Dcreate(file_id, "Velx", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
  status  = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
  status  = H5Sclose(space_id);
  status  = H5Sclose(memspace);
  status  = H5Dclose(dset_id);

  memspace= H5Screate_simple(1, dims, NULL);
  space_id= H5Screate_simple(1, dims, NULL);
  for (unsigned int b=0; b<s.n; b++)
    data[b] = s.part[ids[b]].vel[1];
  dset_id = H5Dcreate(file_id, "Vely", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
  status  = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
  status  = H5Sclose(space_id);
  status  = H5Sclose(memspace);
  status  = H5Dclose(dset_id);

  memspace= H5Screate_simple(1, dims, NULL);
  space_id= H5Screate_simple(1, dims, NULL);
  for (unsigned int b=0; b<s.n; b++)
    data[b] = s.part[ids[b]].vel[2];
  dset_id = H5Dcreate(file_id, "Velz", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
  status  = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
  status  = H5Sclose(space_id);
  status  = H5Sclose(memspace);
  status  = H5Dclose(dset_id);

  memspace= H5Screate_simple(1, dims, NULL);
  space_id= H5Screate_simple(1, dims, NULL);
  for (unsigned int b=0; b<s.n; b++)
    data[b] = s.part[ids[b]].mass;
  dset_id = H5Dcreate(file_id, "Mass", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
  status  = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
  status  = H5Sclose(space_id);
  status  = H5Sclose(memspace);
  status  = H5Dclose(dset_id);

#ifdef OUTPUTPOT
  memspace= H5Screate_simple(1, dims, NULL);
  space_id= H5Screate_simple(1, dims, NULL);
  for (unsigned int b=0; b<s.n; b++)
    data[b] = s.part[ids[b]].pot;
  dset_id = H5Dcreate(file_id, "Potential", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
  status  = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
  status  = H5Sclose(space_id);
  status  = H5Sclose(memspace);
  status  = H5Dclose(dset_id);
#endif

#ifdef OUTPUTACC
  memspace= H5Screate_simple(1, dims, NULL);
  space_id= H5Screate_simple(1, dims, NULL);
  for (unsigned int b=0; b<s.n;b++)
    data[b] = s.part[ids[b]].acc[0];
  dset_id = H5Dcreate(file_id, "Accx", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
  status  = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
  status  = H5Sclose(space_id);
  status  = H5Sclose(memspace);
  status  = H5Dclose(dset_id);

  memspace= H5Screate_simple(1, dims, NULL);
  space_id= H5Screate_simple(1, dims, NULL);
  for (unsigned int b=0; b<s.n; b++)
    data[b] = s.part[ids[b]].acc[1];
  dset_id = H5Dcreate(file_id, "Accy", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
  status  = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
  status  = H5Sclose(space_id);
  status  = H5Sclose(memspace);
  status  = H5Dclose(dset_id);

  memspace= H5Screate_simple(1, dims, NULL);
  space_id= H5Screate_simple(1, dims, NULL);
  for (unsigned int b=0; b<s.n; b++)
    data[b] = s.part[ids[b]].acc[2];
  dset_id = H5Dcreate(file_id, "Accz", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
  status  = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, space_id, H5P_DEFAULT, data);
  status  = H5Sclose(space_id);
  status  = H5Sclose(memspace);
  status  = H5Dclose(dset_id);
#endif

  delete[] data;

  memspace= H5Screate_simple(1, dims, NULL);
  space_id= H5Screate_simple(1, dims, NULL);
  
  unsigned int *id = new unsigned int[s.n];

  for (unsigned int b=0; b<s.n; b++)
    id[b] = s.part[ids[b]].id;
  dset_id = H5Dcreate(file_id, "ID", H5T_NATIVE_UINT, space_id, H5P_DEFAULT);
  status  = H5Dwrite(dset_id, H5T_NATIVE_UINT, memspace, space_id, H5P_DEFAULT, id);
  status  = H5Sclose(space_id);
  status  = H5Sclose(memspace);
  status  = H5Dclose(dset_id);

  delete[] id;

  status  = H5Gclose(handle);
  status  = H5Fclose(file_id);
  if(status < 0) printf("Writing snapshot error\n");
}
}
#endif