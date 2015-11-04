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

#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#include <mpi.h>

#include <sstream>
#include <string>

class Communicator
{
public:
  //MPI vars
  bool mpiOn; 		///<we are using MPI?
  bool mpiRoot; 	///<Is the current process the root process?
  int mpiRank;		///<Which is this process id?
  int mpiNumTasks;
  char hostName[MPI_MAX_PROCESSOR_NAME];
  int lenHostName;
  
  //BLACS vars
  int nProc;
  int nProcRows;
  int nProcCols;
  int context;
  int myId;
  int myRow;		///< The grid row id of the current process
  int myCol;		///< The grid column id of the current process
  
  //Default Matrix parameters;
  int nDefaultBlockRows;	///<The default number of rows per block when creating a distributed matrix.
  int nDefaultBlockCols;	///<The default number of columns per block when creating a distributed matrix.
  
  Communicator(int argc, char **argv);
  ~Communicator();
  
  //returns communicator creation message
  std::string creationMessage();
  
  void error(std::string e);
  
  void setMPIRoot();
  
  void barrier();
  
  void broadcast(int * values, int nValues = 1);
  void broadcast(double * values, int nValues = 1);
  void broadcast(char * values, int nValues = 1);
  void broadcast(bool * value);
  void broadcast(std::string & s);
  
  /**
   * Gather arrays from different nodes.
   * 
   * Gather arrays from different nodes and stores in root. The returned pointer must be free'd after use.
   * 
   * \param values Pointer to the local array of values.
   * \param nValues Number of values in the local array.
   * \return arrays joined together. Size is nValues*mpiNumTasks. The returned pointer must be free'd after use. Valid pointer only on root. NULL on other nodes.
   */
  int * gather(int * values, int nValues);
  double * gather(double * values, int nValues);
  
  /**
   * Gather arrays of different sizes from different nodes.
   * 
   * Gather arrays of different from different nodes and stores in root. The returned pointer must be free'd after use.
   * 
   * \param values Pointer to the local array of values.
   * \param nValues Number of values in the local array.
   * \param totalSize pointer to a variable where the total size of the result will be stored.
   * \return arrays joined together. Size is the value in totalSize. The returned pointer must be free'd after use. Valid pointer only on root. NULL on other nodes.
   */
  int * asymmetricGather(int * values, int nValues, int * totalSize);
  double * asymmetricGather(double * values, int nValues, int * totalSize);
  
  void mpiDebug(int pid = -1);
};

#endif
