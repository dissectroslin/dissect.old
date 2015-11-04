/*****************************************************************************
 *  DISSECT: a tool for performing genomic analysis with large sample sizes
 *  Copyright (C) 2014-2015 Oriol Canela-Xandri and Albert Tenesa
 *                          The Roslin Institute (University of Edinburgh)
 *
 *  This file is part of DISSECT.
 *
 *  DISSECT is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  DISSECT is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with DISSECT.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

#include "communicator.h"
#include "global.h"
#include "misc.h"
#include "options.h"

#include <mpi.h>
#include <omp.h>
#include <cmath>
#include <cstring>
#include <iostream>
#include <string>
#include <stdlib.h>

Communicator::Communicator(int argc, char **argv)
{
  // Initiate MPI
  int tmp;
  tmp = MPI_Init(&argc, &argv);
  if (tmp != MPI_SUCCESS) {
    error("Error: MPI can not be started. Terminating.");
  }
  MPI_Comm_rank(MPI_COMM_WORLD, &this->mpiRank);
  MPI_Comm_size(MPI_COMM_WORLD, &this->mpiNumTasks);
  MPI_Get_processor_name(this->hostName, &this->lenHostName);
  this->setMPIRoot();
  
  // Initiate Cblas context 
  // For the moment I assume that the number of processes have a integer square root 
  //int temp = ceil(sqrt(float(this->mpiNumTasks)));
  //this->nProcCols = temp;
  //this->nProcRows = temp;

  // allow a non-square number of processes
  //this->nProcCols=sqrt(this->mpiNumTasks);
  //this->nProcRows=(this->mpiNumTasks+this->nProcCols-1)/this->nProcCols; //This solution have problemes for some numbers: e.g. 60
  
  this->nProcCols = floor(sqrt(this->mpiNumTasks));
  while(this->mpiNumTasks % this->nProcCols != 0)
  {
    this->nProcCols--;
    if(this->nProcCols <= 1)
    {
      error("Error: The number of MPI tasks can not be a prime number.");
    }
  }
  this->nProcRows = this->mpiNumTasks/this->nProcCols;
  if(this->nProcCols*this->nProcRows != this->mpiNumTasks)
  {
    error("Error: The number of MPI tasks can not be a prime number.");
  }

  //this->nDefaultBlockRows = 3;
  this->nDefaultBlockRows = 64;
  this->nDefaultBlockCols = 64;


  //allow block size to be overwritten by an environment variable
  char* tmpBS;
  tmpBS = getenv ("BLACS_BLOCKSIZE");
  if (tmpBS!=NULL){
    this->nDefaultBlockRows = atoi(tmpBS);
    this->nDefaultBlockCols = atoi(tmpBS);
    if(this->nDefaultBlockRows <= 0 || this->nDefaultBlockCols <= 0)
    {
      error("Error: The environment variable BLACS_BLOCKSIZE is not a valid integer value.");
    }
  }



  Cblacs_pinfo(&this->myId, &this->nProc);
  Cblacs_get(0, 0, &this->context);
  Cblacs_gridinit(&this->context, "Row-major", this->nProcRows, this->nProcCols); //Attention: This must be Row-major for compatibility with MPI IO
  Cblacs_pcoord(this->context, this->myId, &this->myRow, &this->myCol);
  
  if(this->mpiRoot == true && (this->myRow != 0 || this->myCol != 0))
  {
    error("Error: An internal error happened after initializing BLACS grid.");
  }
}

Communicator::~Communicator()
{
  Cblacs_gridexit(this->context);
  MPI_Finalize();
}

std::string Communicator::creationMessage()
{
  std::stringstream messageBuffer;
  messageBuffer << this->nProcRows << "x" << this->nProcCols << " BLACS grid successfully created." << std::endl;
  messageBuffer << "Using " << this->nProcRows*this->nProcCols << " process" << ((this->nProcRows*this->nProcCols>1)?"es":"") << "." << std::endl;
  messageBuffer << "BLACS block size: "<< this->nDefaultBlockRows << "x" << this->nDefaultBlockCols << std::endl;
  return messageBuffer.str();
}

void Communicator::error(std::string e)
{
  std::cerr << "**********************************" << std::endl;
  std::cerr << e << std::endl;
  std::cerr << "**********************************" << std::endl;
  MPI_Abort(MPI_COMM_WORLD, -1);
}

void Communicator::setMPIRoot()
{
  this->mpiRoot = this->mpiRank == 0;
}

void Communicator::barrier()
{
  Cblacs_barrier(communicator->context, "All");
}

void Communicator::broadcast(int * values, int nValues)
{
  MPI_Bcast(values, nValues, MPI_INT, 0, MPI_COMM_WORLD);
}

void Communicator::broadcast(double * values, int nValues)
{
  MPI_Bcast(values, nValues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void Communicator::broadcast(char * values, int nValues)
{
  MPI_Bcast(values, nValues, MPI_CHAR, 0, MPI_COMM_WORLD);
}

void Communicator::broadcast(std::string & s)
{
  int size = s.size();
  broadcast(&size);
  
  char * cs = new char [size + 1];
  if(this->mpiRoot)
  {
    std::strcpy(cs, s.c_str());
  }
  
  broadcast(cs, size + 1);
  
  s = std::string(cs, size);
  
  delete [] cs;
}

void Communicator::broadcast(bool * value)
{
  int sendValue;
  if(this->mpiRoot)
  {
    sendValue = ((*value==true)?1:0);
  }
  MPI_Bcast(&sendValue, 1, MPI_INT, 0, MPI_COMM_WORLD);
  *value = ((sendValue==1)?true:false);
}

int * Communicator::gather(int * values, int nValues)
{
  int success;
  int * result = NULL;
  if(this->mpiRoot)
  {
    result = new int [this->mpiNumTasks*nValues];
  }
  success = MPI_Gather( values, nValues, MPI_INT, result, nValues, MPI_INT, 0, MPI_COMM_WORLD);
  if(success != MPI_SUCCESS)
  {
    misc.error("Error: Error on MPI gather communication. Terminating.", 0);
  }
  return result;
}

double * Communicator::gather(double * values, int nValues)
{
  int success;
  double * result = NULL;
  if(this->mpiRoot)
  {
    result = new double [this->mpiNumTasks*nValues];
  }
  success = MPI_Gather( values, nValues, MPI_DOUBLE, result, nValues, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if(success != MPI_SUCCESS)
  {
    misc.error("Error: Error on MPI gather communication. Terminating.", 0);
  }
  return result;
}

int * Communicator::asymmetricGather(int * values, int nValues, int * totalSize)
{
  int success;
  int * sizes = NULL;
  int * displacements = NULL;
  
  if(this->mpiRoot)
  {
    displacements = new int [this->mpiNumTasks];
  }
  
  //Gather the sizes of each array to gather. Only on root. On ohter processes sizes == NULL.
  sizes = gather(&nValues, 1);

  //Compute the total resulting size and displacements. Only in root.
  *totalSize = 0;  
  if(this->mpiRoot)
  {
    for(int i = 0; i<this->mpiNumTasks; i++)
    {
      displacements[i] = *totalSize;
      *totalSize += sizes[i];
    }
  }
  
  int * result = NULL;
  if(this->mpiRoot)
  {
    result = new int [*totalSize];
  }

  success = MPI_Gatherv(values, nValues, MPI_INT, result, sizes, displacements, MPI_INT, 0, MPI_COMM_WORLD);
  if(success != MPI_SUCCESS)
  {
    //char error_string[BUFSIZ];
    //int length_of_error_string, error_class;
    
    //MPI_Error_class(success, &error_class);
    //MPI_Error_string(error_class, error_string, &length_of_error_string);
    std::cout << success << std::endl;
    misc.error("Error: Error on MPI gather communication. Terminating.", 0);
  }
  
  if(this->mpiRoot)
  {
    delete [] sizes;
    delete [] displacements;
  }
  
  return result;
}

double * Communicator::asymmetricGather(double * values, int nValues, int * totalSize)
{
  int success;
  int * sizes = NULL;
  int * displacements = NULL;
  
  if(this->mpiRoot)
  {
    displacements = new int [this->mpiNumTasks];
  }
  
  //Gather the sizes of each array to gather. Only on root. On ohter processes sizes == NULL.
  sizes = gather(&nValues, 1);
  
  //Compute the total resulting size and displacements. Only in root.
  *totalSize = 0;  
  if(this->mpiRoot)
  {
    for(int i = 0; i<this->mpiNumTasks; i++)
    {
      displacements[i] = *totalSize;
      *totalSize += sizes[i];
    }
  }
  
  double * result = NULL;
  if(this->mpiRoot)
  {
    result = new double [*totalSize];
  }
  
  success = MPI_Gatherv(values, nValues, MPI_DOUBLE, result, sizes, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if(success != MPI_SUCCESS)
  {
    //char error_string[BUFSIZ];
    //int length_of_error_string, error_class;
    
    //MPI_Error_class(success, &error_class);
    //MPI_Error_string(error_class, error_string, &length_of_error_string);
    std::cout << success << std::endl;
    misc.error("Error: Error on MPI gather communication. Terminating.", 0);
  }
  
  if(this->mpiRoot)
  {
    delete [] sizes;
    delete [] displacements;
  }
  
  return result;
}

void Communicator::mpiDebug(int pid)
{
  //Code for attaching gdb to a particular process
  if((communicator->mpiRank == pid || pid < 0) /*&& options.mpiDebug*/){
   int i = 0;
   std::cout << "PID " << communicator->mpiRank << " on " << communicator->hostName << " ready for attach\n" << std::endl;
   fflush(stdout);
   while (0 == i)
   {
    }
  }
}
