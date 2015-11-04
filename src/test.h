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

#ifndef TEST_H
#define TEST_H

#include <mpi.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

#include <unistd.h>

#include <omp.h>

#include "main.h"
#include "global.h"
#include "matrix.h"
#include "options.h"
#include "communicator.h"
#include "misc.h"
#include "genotype.h"
#include "grm.h"
#include "kernel.h"
#include "phenotype.h"
#include "covariate.h"
#include "reml.h"
#include "simulatephenotype.h"
#include "pca.h"
#include "auxiliar.h"
#include "covariancematrix.h"
#include "analysis.h"
#include "blockmatrix.h"

void test1();
void test2();
void test2bis();
void test3();
void test4();
void test5();
void test6(); //Test filtering columns
void test7(); //Test filtering individuals in grm, phenotypes and covariates and genotypes
void test8(); //Test copySubMatrix() method of Matrix class.
void test9(); //Test asymmetric filtering individuals in grm
void test10(); //joinMatrices()
void test11(); //Matrix().invert()
void test12(); //Test define genotype Regions.
void test13(); //Test communicator gather and asymmetricGather
void test14(); //Test grm prunning
void test15(); //Test covariate syncronization
void test16(); //Test matrix eigenvalues and matrix bending/
void test17(); //Test covariance matrix functions
void test18(); //Test joining genotypes.
void test19(); //Test grm addition.
void test20(); //Test grm construction from different sources.
void test20bis(); //Test grm construction from different sources after switching to Kernel class.
void test21(); //Test Matrix write and read
void test22(); //Test loading matrices from different threads.
void test23(); //Test communicator asymmetricGather with doubles
void test24(); //Test parallel genotype loading
void test25(); //Test trace of product of two matrices.
void test26(); //Test new function for filter matrix rows and columns.
void test27(); //Test new matrix write/read functions.
void test28(); //Test matrix pack/unpack functions.
void test29(); //Test MPI read/write
void test29bis(); //Test MPI read/write
void test30(); //Test pdtran/pdlaset/pdlacpy change to _
void test31(); //Test pdtran/pdlaset/pdlacpy change to _
void test32(); //Test MPI/noMPI grm outptut
void test33(); //Test symmetricInvert with float/double cases and the determinant computation.
void test34(); //Test matrix transposition.
void test35(); //Test string broadcast.
void test36(); //Test statistical functions;
void test37(); //Test diagonal() Matrix method;
void test38(); //Test covariate NA dealing;
void test39(); //Test box_muller;

void test40(); //Test generalResorting
void test41(); //Test generalResorting
void test42(double *original, int nr, int nc, int bnr, int bnc, std::map<int, int> & rowsOriginDestination, std::map<int, int> & colsOriginDestination); //Test generalResorting2
void test43ARCHER(); //Test generalResorting performance.

void test44(); //Test diagonalDistribution modifications.
void test45(); //Test multiplication after diagonalDistribution modifications.
void test45bis(); //Same as test45 but with matrices of higher dimension.
void test45bis2(); //Same as test45 but with a scaling factor.
void test45bis3(); //Testing symmetric properties after multiplication.
void test46(); //Test element wise division and multiplication.
void test47(); //Test addition with diagonal distributed matrices.
void test48(); //Test inversion with diagonal matrices.
void test49(); //Test trace, diagonal and other stuff with diagonal distributed matrices.
void test50(); //Test diagonal() and transpose()
void test51(); //joinMatrices with diagonal distributed matrices
void test52(); //Test matrixToStandardVector with diagonal distributed matrices.

void test53(); //Test GRM diagonalization.
void test53bis(); //Test GRM diagonalization after switch to Kernel class.
void test54(); //Test GRM diagonalization read/write.
void test54bis(); //Test GRM diagonalization read/write after switch to Kernel class.
void test55(); //Test GRM copying and some matrix checks with diagonal matrices.
void test55bis(); //Test GRM copiing and some matrix checks with diagonal matrices after switch to Kernel class.

void test56(); //Test Phenotype column counting and testing.

void test57(); //Test GRM substitution by Kernel.
void test57aux(GRM *g1, GRM *g2, Kernel *k1, Kernel *k2, bool symmetrize = true);

void test58(); //Test genotype grouping by byOrderedFixedSize

void test59(); //Test scatter int vector: scatterVectorRet() method
void test60(); //test makeIntersectionMatrix() method
void test61(); //Test diagonal() and diagonal() when the matrix is not a square matrix.

void test62(); //Test get dependent columns;

void test63(); //Test --extract option

void test64(); //Test getGlobalIndexOutsideRange() method
void test65(); //Test getGlobalIndexInsideRange() method

void test66(); //Test GRM check overlaping SNPs 

void test67(); //Test BlockMatrix;
void test67bis(); //Test BlockMatrix and error with mkl scalapack inversion
void test67bis2(int nBlocks, int nRows, int nCols, int defBlockRow, int defBlockCol); //Test BlockMatrix and error with mkl scalapack inversion

void test68(); //Test BlockMatrix multiplications with normal distributed matrices.
void test68bis(int nBlocks, int nRows, int nCols, int defBlockRow, int defBlockCol, int colDimensionTestM); //Test BlockMatrix multiplications with normal distributed matrices.

void test69(); //Test Epistasis Kernel computation.

void test70(); //Get table from file.

void test71(); //Test prepareRAWREML

#endif
