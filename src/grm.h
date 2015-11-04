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

/***********************************************
 * 
 * ATTENTION: This class is no longer used. It is replaced by Kernel class.
 * 
 */


#ifndef GRM_H
#define GRM_H

#include "matrix.h"
#include "genotype.h"

#include <string>
#include <vector>
#include <map>

class Options;

class GRM
{
public:
  Matrix *grm;					///<The GRM matrix.  Set to NULL when diagonalized == true.
  Matrix *N;					///<Matrix that stores the number of genotypes shared between two individuals.  Set to NULL when diagonalized == true.
  
  Matrix *eigenValues;                          ///<The eigenvalues of the diagonalized GRM. Set to NULL when diagonalized == false.
  Matrix *eigenVectors;                         ///<The eigenvectors after GRM diagonalization. Set to NULL when diagonalized == false.
  
  std::vector<Individual> individuals;  	///<List of the Individuals Ids. Only in the mpiRoot process. Only valid if asymmetric == false.
  std::vector<std::string> individualIds;	///<List of the individual keys: familyId + @ + individualId. Only in the mpiRoot process. Only valid if asymmetric == false.
  std::map<std::string, int> individualIdsIdx;	///<Map between the individual keys: familyId + @ + individualId and their position in the grm. Only in the mpiRoot process. Only valid if asymmetric == false.
  int nIndividuals;				///<Total number of individuals. Only valid if asymmetric == false.
  
  std::vector<std::string> SNPIds;              ///<SNP ids used for computing the grm.
  
  bool normalized;                                      ///<true if GRM is normalized (i.e. each grm element is divided by the corresponding N element).
  
  bool asymmetric;                                      ///<true if individuals in rows and columns are different. false otherwise.
  
  bool diagonalized;                                    ///<If grm diagonalized, grm is stored as their diagonal (in eigenValues variable) an their eigenvectors (in eigenVectors).  grm = NULL and N = NULL.
  
  std::vector<Individual> individualsRows;              ///<List of the Individuals Ids in rows. Only in the mpiRoot process. Only valid if asymmetric == true.
  std::vector<std::string> individualIdsRows;           ///<List of the individual keys in rows: familyId + @ + individualId. Only in the mpiRoot process. Only valid if asymmetric == true.
  std::map<std::string, int> individualIdsIdxRows;      ///<Map between the individual keys in rows: familyId + @ + individualId and their position in the grm rows. Only in the mpiRoot process. Only valid if asymmetric == true.
  int nIndividualsRows;                                 ///<Total number of individuals in rows. Only valid if asymmetric == true.
  
  std::vector<Individual> individualsCols;              ///<List of the Individuals Ids in columns. Only in the mpiRoot process. Only valid if asymmetric == true.
  std::vector<std::string> individualIdsCols;           ///<List of the individual keys in columns: familyId + @ + individualId. Only in the mpiRoot process. Only valid if asymmetric == true.
  std::map<std::string, int> individualIdsIdxCols;      ///<Map between the individual keys in columns: familyId + @ + individualId and their position in the grm columns. Only in the mpiRoot process. Only valid if asymmetric == true.
  int nIndividualsCols;                                 ///<Total number of individuals in columns. Only valid if asymmetric == true.
  
  /**
   * Constructor that creates an empty grm
   */
  GRM();
  
  /**
   * Constructor that generates a grm from a Genotype class
   * 
   * \param genotype pointer to the genotype class
   */
  GRM(Genotype * genotype, bool normalizeGRM = true);
  
  /**
   * Constructor that loads a grm from a file
   * 
   * \param f the base name of the file that contains the grm (without the extensions .grm.id or .grm.dat)
   */
  GRM(std::string f);
  
  /**
   * Constructor that greates a new GRM from an existent GRM
   * 
   * \param grm The grm that will be copied
   */
  GRM(GRM * grm);
  
  /**
   * Destructor
   */
  ~GRM();
  
  /**
   * Copy parameters from existent grm
   * 
   * \param grm GRM from which parameters will be copied.
   */
  void copyParameters(GRM *grm);
  
  /**
   * Normalize GRM
   * 
   * Normalizes the GRM. Namely, element-wise division: grm/N. normalized flag is set to true.
   */
  void normalize();
  
  /**
   * Denormalize GRM
   * 
   * Denormalizes the GRM. Namely, element-wise multiplication: grm*N. normalized flag is set to false.
   */
  void denormalize();
  
  /**
   * Returns the normalized GRM
   * 
   * Returns the normalized this->grm if this->diagonalized == false or this->eigenValues if this->diagonalized == true
   */
  Matrix * getNormalizedGRM();
  
  /**
   * Compute the GRM from genotype data
   * 
   * \param genotype pointer to genotype data used for computation
   */
  void computeGRM(Genotype * genotype);
  
  /**
   * Save GRM to a file
   * 
   * \param fn the base name of the file (without the extensions .grm.id or .grm.dat)
   */
  void writeGRM(std::string fn);
  
  /**
   * Load GRM from a file
   * 
   * \param fn the base name of the file (without the extensions .grm.id or .grm.dat)
   */
  void readGRM(std::string f);
  
  /**
   * Filter individuals from GRM
   * 
   * \param keepIndividualIds List of individuals that will be kept. The order of this list must be in the same of this->individualsIds
   * \param filterN if true, this->N matrix will be filtered too. Otherwise, it will be deleted.
   */
  void filterIndividuals(std::vector<std::string> & keepIndividualIds, bool filterN = true);
  
  /**
   * Filter individuals from GRM. Different individuals in rows and columns.
   * 
   * Filter individuals uwing different individuals in rows and columns. this->asymmetric flag fill be set to true.
   * 
   * \param keepIndividualIdsRows List of individuals that will be kept in rows. The order of this list must be in the same of this->individualsIds
   * \param keepIndividualIdsCols List of individuals that will be kept in columns. The order of this list must be in the same of this->individualsIds
   * \param filterN if true, this->N matrix will be filtered too. Otherwise, it will be deleted.
   */
  void filterIndividualsAsymmetric(std::vector<std::string> & keepIndividualIdsRows, std::vector<std::string> & keepIndividualIdsCols, bool filterN = true);
  
  /**
   * Add two GRMs
   * 
   * Performs the operation this = scalingFactor1*grm1 + scalingFactor2*grm2. Before the sum, the GRMs are denormalized. grm1 and grm2 must be different than this.
   * If grm2 == NULL, then grm1 will be added to this: this = scalingFactor2*this + scalingFactor1*grm1. If scalingFactors are equal, SNPIds from two GRMs are joined, if
   * scaling factors are different, SNPIds are substracted according to scalingFactors signs.
   * 
   * \param scalingFactor1 Scaling factor for grm1. Now it can only be 1. or -1. indicating the GRM is substracted or added.
   * \param grm1 the first GRM for the sum.
   * \param scalingFactor2 Scaling factor for grm2
   * \param grm2 the second GRM for the sum. Now it can only be 1. or -1. indicating the GRM is substracted or added.
   */
  void addGRMs(double scalingFactor1, GRM *grm1, double scalingFactor2 = 1., GRM *grm2 = NULL);
  
  /**
   * Search the ids of the maximum number of individuals can be kept to eliminate all elements in the GRM greather than a particular cutoff.
   * 
   * \param cutoff the threshold over which the GRM elements must be filtered.
   */
  std::vector<std::string> searchNoHighRelatedIndividuals(double cutoff);
  
  /**
   * Prune a grm removing the highest related individuals
   * 
   * Remove the minimum number of individuals such there are not relations with a value above cutoff
   * 
   * \param cutoff threshold above which relations will be pruned.
   */
  void pruneGRM(double cutoff);
  
  /**
   * Randomly filter this grm
   * 
   * \param fraction fraction of individuals that will be kept.
   * \param limit minimum number of individuals that must be kept.
   */
  void randomSubSample(double fraction, int minimum);
  
  /**
   * Diagonalize this grm and store their eigenvalues and their eigenvectors.
   */
  void diagonalizeGRM();
  
  /**
   * Recover original GRM from their eigenvalues and their eigenvectors.
   */
  void recoverGRMFromEigenDecomposition();
  
  /**
   * Function for debugging pourposes thhat prints the current GRM
   */
  void printGRM(bool showSNPs = false);
};

#endif
