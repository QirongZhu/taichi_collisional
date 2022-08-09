#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void allocate_gadget_memory(void);
void free_gadget_memory(void);
void load_snapshot(char *fname, int files);

struct io_header
{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];	/* fills to 256 Bytes */
} header1;

int NumPart;

struct particle_data
{
  float Pos[3];
  float Vel[3];
  float Mass;
  int Type;
  unsigned int Id;
} *PPP;

double Time;
double Redshift;


/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.  */

void load_snapshot(char *fname, int files)
{
  FILE *fd;
  char buf[200];
  int i, k, dummy, ntot_withmasses;
  int n, pc, pc_new;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

  for(i = 0, pc = 0; i < files; i++, pc = pc_new)
    {
      if(files > 1)
	sprintf(buf, "%s.%d", fname, i);
      else
	sprintf(buf, "%s", fname);

      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s`\n", buf);
	  exit(0);
	}

      //      printf("reading `%s' ...\n", buf);
      //      fflush(stdout);

      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header1, sizeof(header1), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);

      if(files == 1)
	{
	  for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
	    NumPart += header1.npart[k];
	}
      else
	{
	  for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
	    NumPart += header1.npartTotal[k];
	}

      for(k = 0, ntot_withmasses = 0; k < 6; k++)
	{
	  if(header1.mass[k] == 0)
	    ntot_withmasses += header1.npart[k];
	}

      if(i == 0)
	allocate_gadget_memory();

      SKIP;

      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
	      fread(&PPP[pc_new].Pos[0], sizeof(float), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;

      SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
	      fread(&PPP[pc_new].Vel[0], sizeof(float), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;

      SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
	      fread(&PPP[pc_new].Id, sizeof(unsigned int), 1, fd);
	      pc_new++;
	    }
	}
      SKIP;

      if(ntot_withmasses > 0)
	SKIP;

      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
	      PPP[pc_new].Type = k;

	      if(header1.mass[k] == 0)
		fread(&PPP[pc_new].Mass, sizeof(float), 1, fd);
	      else
		PPP[pc_new].Mass = header1.mass[k];
	      pc_new++;
	    }
	}

      if(ntot_withmasses > 0)
	SKIP;

      fclose(fd);
    }

  Time = header1.time;
  Redshift = header1.time;

}

/* this routine allocates the memory for the 
 * particle data.  */
void allocate_gadget_memory(void)
{
  //  printf("allocating memory for N=%d ...\n", NumPart);
  if(!(PPP = (particle_data *) malloc(NumPart * sizeof(struct particle_data))))
    {
      fprintf(stderr, "failed to allocate memory.\n");
      exit(0);
    }
  //  printf("allocating memory...done\n");
}

void free_gadget_memory(void)
{
  free(PPP);
  //  printf("freeing memory...done\n");
}
