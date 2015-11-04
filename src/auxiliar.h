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

#ifndef AUXILIAR_H
#define AUXILIAR_H

#include "grm.h"
#include "genotype.h"

#include <vector>
#include <map>
#include <string>

class Matrix;
class Genotype;
class Kernel;

/**
 * Computes the intersection between many string vectors
 * 
 * Example of usage: intersectionStringVectors(3, &v1, &v2, &v3)
 * Where v1, v2 , v3 are of type std::vector<std::string>.
 * Returns the intersection between v1, v2, v3
 * 
 * \param count The number of vectors to intersect.
 * \param ... List of pointers to string vectors.
 * \return The intersection between the argument vectors
 */
std::vector<std::string> intersectionStringVectors(int count, ...);

/**
 * Returns the intersection of the elements of a vector and a template ordered according to the template.
 * 
 * returns the intersection of the elements of sourceVector and templateVector in the order that appear in the templateVector
 * 
 * \param templateVector The vector that will define the order
 * \param sourceVector The vector of the elements to be sorted.
 * \return The intersection between templateVector and sourceVector ordered according to templateVector.
 */
std::vector<std::string> orderVectorAsTemplate(std::vector<std::string> & templateVector, std::vector<std::string> & sourceVector);

/**
 * Returns the difference between the elements of two vectors.
 * 
 * Returns the a vectors with the elements in baseVector which not are in vectorToSubstract. The two vectors cannot have repeated elements and all elements in vectorToSubstract must be in baseVector.
 * 
 * \param baseVector The vector where elements will be substracted.
 * \param vectorToSubstract The vector that defines elements to be substracted.
 * \return The elements in baseVector which not are in vectorToSubstract.
 */
std::vector<std::string> differenceBetweenTwoVectors(std::vector<std::string> & baseVector, std::vector<std::string> & vectorToSubstract);

/**
 * Return a random sample of a vector elements.
 * 
 * Return a random sample (without repetition) of a vector of elements. Elements
 * are returned in the same order that in the sourceVector.
 * Uses the c rand() generator. So, it do not have "ideal" random distribution.
 * 
 * \param sourceVector The vector from which elements will be get randomly
 * \param nElements The number of elements to get from sourceVector.
 * \return A random sample of nElements from sourceVector in the same order than in sourceVector.
 */
std::vector<std::string> getRandomSample(std::vector<std::string> & sourceVector, int nElements);

/**
 * Get the values in a map using the keys stored in a vector.
 * 
 * Get the values in a map using the keys stored in a vector. The values are returned in the same order appear in the vector
 * 
 * \param keysVector The keys that will be get from valuesMap
 * \param valuesMap The map from which the values will be got
 * \return The values in valuesMap witth the keys in keysVector. The order of the values is the same than the order in keysVector
 */
std::vector<int> extractMapValues(std::vector<std::string> & keysVector, std::map<std::string, int> & valuesMap);

/**
 * Get first column in a file and return in a vector.
 * 
 * Get first column in a file and return in a vector. If file is empty, raise an error.
 * 
 * \param f the name of the file where lines must be read
 * \return A vector storing the first column on each element. Returns on all processes.
 */
void getListFromFile(std::string f, std::vector<std::string> & list);

/**
 * Get first column in a file and return in a vector (ATTENTION: Untested function).
 * 
 * Get first column in a file and return in a vector. If file is empty, raise an error.
 * 
 * \param f the name of the file where lines must be read
 * \return A vector storing the first column on each element. Returns on all processes.
 */
void getListFromFile(std::string f, std::vector<double> & list);

/**
 * Get first n columns in a file and return in a vector of vectors.
 * 
 * Get forst n columns in a file and return in a vector of vectors. If file is empty, or all rows do not have not enough columns, raise an error.
 * 
 * \param f the name of the file where lines must be read
 * \param nColumns The number of columns to get
 * \return A vector of vectors storing the columns on each element. Returns on all processes.
 */
void getTableFromFile(std::string f, std::vector< std::vector<std::string> > & table, int nColumns);

/**
 * Get the number of file columns based on the columns of first line.
 */
int getNumberOfFileColumns(std::string f);

/**
 * Least Common Multiple
 */
int leastCommonMultiple(int a, int b);

/**
 * Greatest Common Divisor
 */
int greatestCommonDivisor(int a, int b);

/**
 * Return the estimated variance of a vector.
 * 
 * \param m matrix that should have only one column or only one row.
 * \return The standard deviation.
 */
double computeVariance(Matrix * m);

/**
 * Int to std::string
 * 
 * \param i Integer to be converted to string
 */
std::string getString(int i);

/**
 * Int to std::string (other name for getString)
 * 
 * \param i Integer to be converted to string
 */
std::string i2s(int i);

/**
 * Replaces the spaces and tabs of a string with underscores.
 * 
 * \param s the string with spaces and/or tabs
 * \return a string with spaces and tabs replaced.
 */
std::string spacetab2underscore(std::string s);

/**
 * Load Genotype using the program arguments.
 * 
 * \return A pointer to the loaded Genotype.
 */
Genotype * loadGenotypeUsingOptions();

/**
 * Create a new GRM using the program arguments.
 * 
 * If GRM is created from genotype file, a pointer to this genotype could be returned.
 * 
 * \param[in] returnGenotype a pointer to the genotype for computing the grm must be returned? If not, all genotypes used will be cleared.
 * \param[out] returnedGenotype Pointer to a genotype pointer. returnGenotype is true and a genotype file is used for computing the GRM. Then, the genotype pointer is returned on this var. If no genotype is loaded, a NULL pointer is returned.
 * \return A pointer to the created GRM.
 */
Kernel * loadGRMUsingOptions(bool returnGenotype = false, Genotype **returnedGenotype = NULL);

/**
 * Returns a list of phenotypes to analyze based on options.
 */
std::vector<int> getPhenotyesForAnalysis();

#define MBIG 1000000000 // According to Knuth, any large MBIG, and
#define MSEED 161803398 // any smaller (but still large MSEED can
                        // be substituted for the above values.
#define MZ 0           
#define FAC (1.0/MBIG)

/**
 * Returns a uniform random deviate between 0.0 and 1.0.
 * Set idnum to any negative value to initialize or reinitialize the sequence.
 * From numerical recipes. Implemented by Manuel Valenzuela (23 enero 1996).
 * Maybe this should be rewtritten. I am not sure this will be necessary.
 * 
 * \param idnum seed
 */
float ran3(long *idnum);

/**
 * Returns a uniform random variate between [0.0, 1.0].
 * (based on NumRec routine rand3(), based itself on this-or-that
 * from Knuth.  Not a linear congruence generator.
 * To initialize/reinitialize, pass it a negative long int; it has a 
 * memory, so passing it the same initializer multiple times
 * during a run of the program will produce different values.
 * 
 * \param idnum seed
 */
double unif_rand_dbl(long *idum);

/**
 * normal random variate generator
 * 
 * \param m mean
 * \param s standard deviation
 * \param idnum seed
 */
double box_muller(double m, double s, long *idnum);

double chi1_CDF(int df, double x);

double FStatCDF(double df1, double df2, double x);

double tStatCDF(double df, double x);

std::string getFileName(std::string & path);

/////////////////////////////////////////
//Just for gdb and debugging purposes:

int compareMatrices(Matrix * m1, Matrix * m2, double threshold = -1.);
int compareGlobalMatrices(Matrix * m1, Matrix * m2, double threshold = -1.);

GRM * loadGRMUsingOptionsOld(bool returnGenotype = false, Genotype **returnedGenotype = NULL);

#endif
