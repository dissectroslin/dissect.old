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

#ifndef COVARIANCEMATRIX_H
#define COVARIANCEMATRIX_H

#include "matrix.h"
#include "blockmatrix.h"
#include "kernel.h"

#include <vector>
#include <map>
#include <set>

class Matrix;
class REML;
class Phenotype;
class Covariate;

enum VarianceType
{
  variance,
  covariance,
  standardDeviation
};

enum VarianceTypeEffect
{
  genetic,
  environment
};

enum CovarianceMatrixStatus
{
  covarianceMatrix,
  inverseCovarianceMatrix,
  covarianceMatrixByBlocks,
  inverseCovarianceMatrixByBlocks,
  undefined
};

enum ElementType
{
  grmType,
  environtmentalType
};

/**
 * Element of the sum which forms a total covariance matrix
 * This class contains information of the current element. In particular, how a covariance submatrix is combined with the different variances
 * and how them will be added to the total covariance.
 * An element is formed by the product of a matrix and a set of variances. e.g.  M*v1*v2*v3
 */
class Element
{
public:
  std::string subCovarianceId;                  ///<This element joined together with others can be part of a covariance matrix of the overall covariance. This is their shared name. All elements with same Id, should create a symmetric matrix when joined together.
  std::string name;                             ///<Element name which is used for identifying the element.
  ElementType type;                             ///<The type of the matrix in this element.
  Matrix * m;                                   ///<Pointer to a covarianceSubMatrices of CovarianceMatrix class. This is the sub matrix that forms the element.
  subMatrix outcomeSubMatrix;                   ///<Where in the total covariance matrix this elemenent will be added.
  subMatrix sourceSubMatrix;                    ///<Which part of this element subMatrix will be added to the total covariance matrix.
  std::pair<int, int> blockPosition;            ///<If the covariancematrix is made of blocks, in which block is this element.
  double constantFactor;                        ///<A constant scaling factor that multiply this element.
  std::vector<std::string> nameVariances;       ///<The name of the variances of CovarianceMatrix that multiply this element.
  std::vector<int> iVariances;                  ///<Indexs to the variances of CovarianceMatrix that multiply this element. In the same order than nameVariances.
  std::vector<VarianceType> varianceTypes;      ///<This variance is a variance/covariance or is the square root (i.e. standard deviation)? 
};

class Variance
{
public:
  double variance;                 ///<variance value
  double initialVariance;          ///<the initial value before REML start
  std::string name;                ///<name
  std::string group;               ///<Defines on which group variance belongs.
  VarianceType type;               ///<variance type
  VarianceTypeEffect typeEffect;   ///<Accounts for genetic variance or environment variance
  std::set<int> inElements;        ///<Variable indicating in which elements appears this variance
};

/**
 * Global covariance Matrix formed by sub covariance matrices and their variances.
 * This class manages the covariance matrix and is able to return the actualized matrix after variances change or the derivate of the matrix with respect a particular variance.
 */
class CovarianceMatrix
{
public:
  int dimension;                                                ///<The dimension of the covariance matrix
  
  std::vector<Matrix*> covarianceSubMatrices;                   ///<Array of the covariance sub matrices that after summed together will form the total covariance.
  std::vector<std::string> covarianceSubMatrixNames;            ///<Array of covariance submatrix matrix names.
  int nCovarianceSubMatrices;                                   ///<Number of covariance sub matrices.
  
  Matrix * m;                                                   ///<after calling computeCovariance(), this Matrix stores the covariance.
  CovarianceMatrixStatus mStatus;                               ///<Describes which matrix is stored in m.
  Matrix * mDerivate;                                           ///<after calling computeDerivateCovariance(idx), this matrix stores the derivate of the covariance matrix with respect the variance idx.
  //bool mDerivateAllocated;
  
  std::vector<Variance> variances;                              ///<Array of variances
  std::map<std::string, double> varianceGroupExpectedMagnitudes;          ///<
  
  std::vector<Element> elements;                                ///<Array of the elements that form the covariance matrix sum.
  
  bool byBlocks;                                                ///<If byBlocks is true, all elements have properly defined a valid blockPosition i.e. != std::pair<int, int>(-1, -1)
  bool allCovarianceSubMatricesDiagonal;                        ///<If true, all matrices in covarianceSubMatrices are diagonal.
  int blockDimension;                                           ///<Number of blocks in each dimension of the covariance matrix.
  BlockMatrix mInBlocks;                                        ///<Covariance matrix in blocks. It is only used when byBlocks and allCovarianceSubMatricesDiagonal are true.
  
  CovarianceMatrix(int dimension, bool useDiagonalDistribution = false, int tBlockDimension = -1);
  ~CovarianceMatrix();
  
  /**
   * Insert a new covariance submatrix
   * 
   * This method inserts a covariance submatrix. To this end, the variables covarianceSubMatrices, covarianceSubMatrixNames, and nCovarianceSubMatrices are modified.
   * Once a matrix is inserted, a copy of the reference is added. It should be considered that this reference now it is owned by this class. It will manage and delete it.
   * 
   * \param name The name of the new covariance submatrix. This must be different than the previous inserted names
   * \param m Pointer to the matrix. Only a reference of the matrix is copied.
   */
  void insertCovarianceMatrix(std::string name, Matrix * m);
  
  /**
   * Insert a new covariance submatrix from a Kernel matrix
   * 
   * This method inserts a covariance submatrix. To this end, the variables covarianceSubMatrices, covarianceSubMatrixNames, and nCovarianceSubMatrices are modified.
   * Once a matrix is inserted, a copy of the reference is added. It should be considered that this reference now it is owned by this class. It will manage and delete it.
   * 
   * \param name The name of the new covariance submatrix. This must be different than the previous inserted names
   * \param kernel Pointer to the Kernel. Only a reference of the kernel will be copied.
   */
  void insertCovarianceMatrix(std::string name, Kernel * kernel);
  
  /**
   * Define a new variances group
   * 
   * This method creates a new variances group and defines a expected magnitude for this group.
   * 
   * \param name The grou pname. This must be different than the previous created groups.
   * \param groupExpectedMagnitude Expected magnitude for this group.
   */
  void insertVarianceGroup(std::string name, double groupExpectedMagnitude);
  
  /**
   * Insert a new variance
   * 
   * This method inserts a new variance. To this end, the variables nVariances, variances, varianceNames, varianceTypes, varianceInElements are modified.
   * 
   * \param name Name of the new variance. This must be different than the previous inserted names
   * \param groupName Which group this variance belongs.
   * \param type This is a variance or a covariance?
   * \param initialValue Initial value of the variance.
   */
  void insertVariance(std::string name, std::string groupName, VarianceType type, VarianceTypeEffect typeEffect, double initialValue);
  
  /**
   * Insert a new element to the covariance matrix
   * 
   * This method inserts a new element. Only elements in the upper triangular part of the matrix will be used.
   * 
   * \param elementName The name of the new element. This must be unique.
   * \param type This element type.
   * \param insertedCovarianceMatrixName Name of the subCovariance matrix the will be on the element. This must be inserted before in the current object by using insertCovarianceMatrix() method.
   * \param outcomeSubMatrix define the outcomeSubMatrix of the element.
   * \param sourceSubMatrix define the sourceSubMatrix of the element.
   */
  void insertElement(std::string elementName, ElementType type, std::string insertedCovarianceMatrixName, double constantFactor = 1., subMatrix outcomeSubMatrix = subMatrix(), subMatrix sourceSubMatrix = subMatrix());
  
  /**
   * Insert a new element to the covariance matrix
   * 
   * This method inserts a new element.
   * 
   * \param subCovarianceMatrixId The name of the subcovariance matrix it forms part.
   * \param elementName The name of the new element. This must be unique.
   * \param type This element type.
   * \param insertedCovarianceMatrixName Name of the subCovariance matrix the will be on the element. This must be inserted before in the current object by using insertCovarianceMatrix() method.
   * \param blockPosition If the covariancematrix is made of blocks, in which block is this element. It must be std::pair<int, int>(-1, -1) otherwise.
   * \param outcomeSubMatrix define the outcomeSubMatrix of the element.
   * \param sourceSubMatrix define the sourceSubMatrix of the element.
   */
  void insertElement(std::string subCovarianceMatrixId, std::string elementName, ElementType type, std::string insertedCovarianceMatrixName, double constantFactor = 1., std::pair<int, int> blockPosition = std::pair<int, int>(-1, -1), subMatrix outcomeSubMatrix = subMatrix(), subMatrix sourceSubMatrix = subMatrix());
  
  /**
   * Append a variance to an element.
   * 
   * This function adds a variance in the selected element. Many variances can be added to the same element.
   * 
   * \param elementName The name of the element where the variance will be inserted.
   * \param insertedVarianceName Name of the variance that will be on the element. This must be already defined in the current object by using insertVariance() method.
   * \param insertedVarianceTypes Variance type. Allowed values: variance, covariance, standardDeviation.
   */
  void appendVarianceToElement(std::string elementName, std::string varianceName, VarianceType varianceType);
  
  /**
   * Return index of an element
   * 
   * \param elementName The name of the element to be returned.
   * \return The index to the element with name elementName in the vector this->elements
   */
  int getElement(std::string elementName);
  
  /**
   * Return index of an element
   * 
   * \param varianceName The name of the variance to be returned.
   * \return The index to the variance with name varianceName in the vector this->variances
   */
  int getVariance(std::string varianceName);
  
  /**
   * Reinitialize variances to their initial value.
   */
  void reinitializeVariances();
  
  /**
   * Remove all variances from an Element.
   * 
   * \param elementName Name of the element where the variances will be removed.
   */
  void clearElementVariances(std::string elementName);
  
  /**
   * Change the constant factor of an element.
   * 
   * \param elementName Name of the element where the constantFactor will be changed.
   * \param constantFactor The new constantFactor value.
   */
  void changeElementConstantFactor(std::string elementName, double constantFactor);
  
  /**
   * Remove a variance from the CovarianceMatrix class.
   * 
   * The removed variance must not be in any element.
   * 
   * \param name Name of the variance to eliminate.
   */
  void deleteVariance(std::string name);
  
  /**
   * Return the Covariance derivate with respect a particular variance
   * 
   * The method computes the Covariance derivate with respect a particular variance. Stores it in mDerivate variable and returns the pointer to the matrix.
   * 
   * \param idxDerivateVariance The index of the variance in variances used for deriving the covariance matrix.
   * \return the Covariance derivate with respect the idxDerivateVariance variance.
   */
  Matrix * computeDerivateCovariance(int idxDerivateVariance);
  
  /**
   * Compute the total covariance
   * 
   * Compute the total covariance and stores the result on m.
   */
  void computeCovariance();
  
  /**
   * Compute the total covariance in block form
   */
  void computeBlockCovariance();
  
  /**
   * Compute the variance part of an element
   * 
   * \param element The elemento for which compute the variance part
   */
  double computeElementVariance(Element & element);
  
  /**
   * Compute the variance part of the element which is part of the subcovariance defined by subCovarianceId
   * 
   * If the current subcovariance is formed by more then one element, an error is raised.
   * 
   * \param subCovarianceId A string defining the subcovariance.
   */
  double computeSubCovarianceVariance(std::string subCovarianceId);
  
  /**
   * Compute the derivate of the variance part of an element with respect the variance of index idxDerivateVariance
   * 
   * The element must have this variance between their associated variances. Otherwise an error will be raised.
   * 
   * \param element The elemento for which compute the variance derivate part.
   * \param idxDerivateVariance The index of the variance to derive.
   */
  double computeElementVarianceDerivate(Element & element, int idxDerivateVariance);
  
  /**
   * invert the covariance
   * 
   * Compute the the inverse of the covariance matrix and store in m.
   */
  bool invertCovariance(double * logDeterminant, bool useSinglePrecision = false);
  
  /**
   * Constrain the variances Method 1
   * 
   * If a variance of type variance becomes smaller than 0, then it will be constrained to a positive value.
   * ATENTION: It is based on the GCTA approach. This method probably must be improved since it can create new negative variances?? Now is 'improved', but need checking.
   * 
   * \return The number of constrained variances.
   * 
   */
  int constrainVariancesM1();
  
  /**
   * Constrain the variances Method 2
   * 
   * If a variance of type variance becomes smaller than 0, then it will be constrained to a minimum positive value.
   * This method searches the minimum scaling factor for delta that makes that all the variances are at least these minimum positive value.
   * Then variances are recomputed as variances = oldVariances + delta/scalingFactor;
   * 
   * \return The scaling factor.
   */
  double constrainVariancesM2(std::vector<double> & oldVariances, std::vector<double> & delta);
  
  /**
   * Fix variances to 0
   * 
   * If a variance of type variance multiplied by 1e20 becomes smaller than their group expected Magnitude (defined by varianceGroupExpectedMagnitudes),
   * then it will fixed to 0.
   * 
   * \return The number of variances fixed to 0.
   */
  int fixVariancesToZero();
  
  /**
   * Return a string with the variance names.
   * 
   * \param separator string used to separate the names.
   * \return string with variance names.
   */
  std::string getStringVarianceNames(std::string separator = " ");
  
  /**
   * Return a string with the variance values.
   * 
   * \param separator string used to separate the values.
   * \return string with variance names.
   */
  std::string getStringVarianceValues(std::string separator = " ");
  
  /**
   * Get the name of a matrix in CovarianceMatrix the pointer of the matrix
   */
  std::string getMatrixName(Matrix *m);
  
  /**
   * Broadcast the variance values in the root node to the other nodes.
   * 
   * It is based on communicator->broadcast.
   */
  void broadcastVariances();

  /**
   * Broadcast the initial variance values in the root node to the other nodes.
   * 
   * It is based on communicator->broadcast.
   */
  void broadcastInitialVariances();
  
  /**
   * Return the variance values in a std::vector
   * 
   * \return A vector of variances
   */
  std::vector<double> getVectorVariances();
  
  /**
   * get the variance values from a std::vector and copy to this->variances[].variance
   * 
   * \return A vector of variances
   */
  void storeVectorVariances(std::vector<double> & vv);
  
  /**
   * Return the matrix of element in the proper context.
   * 
   * The matrices in elements are stored as a matrix and their position in covariance matrix. In a single REML, the element matrix have same dimensions than
   * covariance matrix, but in bivariate REML, not: the element matrix could be in a submatrix of covariance matrix. This function returns the element matrix
   * inserted in a matrix of the size of covariance matrix with 0,'s in the components not filled by the element.
   * 
   * \param i the index of the matrix element wants to be returned
   * \param m1 A pointer to a matrix with dimensions of covariance matrix. This will be filled with the element matrix.
   */
  void getElementMatrixInContext(int i, Matrix * m1);
  
  /**
   * Return a sub covariance.
   * 
   * The covariance matrix is formed by the sum of subcovariance matrices which in turn are formed by one or more elements (that share element.subCovarianceId == subCovarianceId).
   * This function returns the selected subcovariance matrix. All these subcovariance matrices are expected to be symmetrical.
   * 
   * \param subCovarianceId The id defining the sub covariance matrix.
   */
  Matrix * getSubCovariance(std::string subCovarianceId);
  
  /**
   * Performs the product between a subcovariance and a matrix.
   * 
   * Performs the product scale*subV*m where subV is a subcovariance defined by subCovarianceId and scale is a scalar factor.
   * 
   * \param subCovarianceId Id of the subcovariance matrix to be multiplied with m.
   * \param m Matrix to be multiplied to the subcovariance
   * \param scale Factor to multiply on the product
   * \return The product
   */
  Matrix * multiply(std::string subCovarianceId, Matrix * m, double scale = 1.);
  
  /**
   * Set the variance values
   * 
   * \param newVariances Set the values on variances vector to values in newVariances. The names on both vectors must be the same.
   */
  void setVarianceInitialValues(std::vector<Variance> newVariances);
  
  /**
   * Debugging method
   */
  void showCovarianceMatrix(bool showMatrices = true);
};

#endif
