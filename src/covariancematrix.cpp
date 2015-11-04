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

#include "matrix.h"
#include "blockmatrix.h"
#include "communicator.h"
#include "misc.h"
#include "global.h"
#include "reml.h"
#include "kernel.h"
#include "options.h"
#include "covariancematrix.h"

#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <iomanip>

CovarianceMatrix::CovarianceMatrix(int dimension, bool useDiagonalDistribution, int tBlockDimension)
{
  this->dimension = dimension;
  
  this->nCovarianceSubMatrices = 0;
  
  this->blockDimension = tBlockDimension;
  if( tBlockDimension >= 1)
  {
    this->byBlocks = true;
  }
  
  this->allCovarianceSubMatricesDiagonal = useDiagonalDistribution;
  
  DistributionType baseMatrixDistribution;
  if( useDiagonalDistribution == true && tBlockDimension == 1 )
  {
    baseMatrixDistribution = diagonalDistribution;
  }
  else
  {
    baseMatrixDistribution = MATRIX_DEFAULT_DISTRIBUTION;
  }
  
  this->m = new Matrix(baseMatrixDistribution, this->dimension, this->dimension, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  this->mStatus = undefined;
  this->mDerivate = new Matrix(baseMatrixDistribution, this->dimension, this->dimension, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  
  this->mInBlocks = BlockMatrix();
}

CovarianceMatrix::~CovarianceMatrix()
{
  delete m;
  
  for(int i = 0; i<this->nCovarianceSubMatrices; i++)
  {
    if(this->covarianceSubMatrices[i] != NULL)
    {
      delete this->covarianceSubMatrices[i];
      this->covarianceSubMatrices[i] = NULL;
    }
  }
  
  delete this->mDerivate;
}

void CovarianceMatrix::insertCovarianceMatrix(std::string name, Matrix * m)
{
  if(name == "")
  {
    misc.error("Error: covariance matrix can not be created with name: '" + name + "'.", 0);
  }
  for(int i = 0; i<this->nCovarianceSubMatrices; i++)
  {
    if(this->covarianceSubMatrixNames[i] == name)
    {
      misc.error("Error: There are two covariance matrices with the same name: '" + name + "'", 0);
    }
  }
  if( (this->allCovarianceSubMatricesDiagonal == true && m->distribution != diagonalDistribution) || 
      (this->allCovarianceSubMatricesDiagonal == false && m->distribution == diagonalDistribution)  )
  {
    misc.error("Error: An internal error was happened. Diagonal matrices cannot be mixed with non diagonal matrices in the covariance matrix.", 0);
  }
  
  if(m->symmetric == true)
  {
    m->symmetrizeTriangularMatrix();
  }
  this->covarianceSubMatrixNames.push_back(name);
  this->covarianceSubMatrices.push_back(m);
  this->nCovarianceSubMatrices = this->covarianceSubMatrices.size();
}

void CovarianceMatrix::insertCovarianceMatrix(std::string name, Kernel * kernel)
{
  insertCovarianceMatrix(name, kernel->getNormalizedKernel());
  
  kernel->kernel = NULL;
  kernel->eigenValues = NULL;
  delete kernel;
}

void CovarianceMatrix::insertVarianceGroup(std::string name, double groupExpectedMagnitude)
{
  if(name == "")
  {
    misc.error("Error: this group can not be created with name: '" + name + "'.", 0);
  }
  if(this->varianceGroupExpectedMagnitudes.count(name) != 0)
  {
    misc.error("Error: There are two groups with the same name: '" + name + "'", 0);
  }
  this->varianceGroupExpectedMagnitudes[name] = groupExpectedMagnitude;
}

void CovarianceMatrix::insertVariance(std::string name, std::string groupName, VarianceType type, VarianceTypeEffect typeEffect, double initialValue)
{
  if(name == "")
  {
    misc.error("Error: Invalid variance name.", 0);
  }
  if(type != variance && type != covariance)
  {
    misc.error("Error: Invalid variance type.", 0);
  }
  for(int i = 0; i<this->variances.size(); i++)
  {
    if(this->variances[i].name == name)
    {
      misc.error("Error: There are two variances with the same name: '" + name + "'", 0);
    }
  }
  if(this->varianceGroupExpectedMagnitudes.count(groupName) != 1)
  {
    misc.error("Error: When inserting a new variance, there are not any group with the name: '" + groupName + "'", 0);
  }
  
  Variance newVariance;
  newVariance.name = name;
  newVariance.group = groupName;
  newVariance.type = type;
  newVariance.typeEffect = typeEffect;
  newVariance.variance = initialValue;
  newVariance.initialVariance = initialValue;
  newVariance.inElements = std::set<int>();

  this->variances.push_back(newVariance);
}

void CovarianceMatrix::insertElement(std::string subCovarianceMatrixId, std::string elementName, ElementType type, std::string insertedCovarianceMatrixName, double constantFactor, std::pair<int, int> blockPosition, subMatrix outcomeSubMatrix, subMatrix sourceSubMatrix)
{
  bool flag;
  Element newElement;
  
  if(getElement(elementName) >= 0)
  {
    misc.error("Error: Invalid element name. The element names must be unique.", 0);
  }
  
  newElement.subCovarianceId = subCovarianceMatrixId;
  newElement.name = elementName;
  newElement.type = type;
  
  flag = false;
  for(int i = 0; i<this->nCovarianceSubMatrices; i++)
  {
    if(this->covarianceSubMatrixNames[i] == insertedCovarianceMatrixName)
    {
      newElement.m = this->covarianceSubMatrices[i];
      flag = true;
    }
  }
  if(!flag)
  {
    misc.error("Error: Invalid covariance matrix name when constructing an Element: '" + insertedCovarianceMatrixName + "'", 0);
  }
  
  newElement.iVariances.clear();
  newElement.nameVariances.clear();
  newElement.varianceTypes.clear();
  
  newElement.constantFactor = constantFactor;

  if( this->byBlocks == true )
  {
    if( blockPosition.first < 0 || blockPosition.second < 0 || blockPosition.first >= this->blockDimension || blockPosition.second >= this->blockDimension )
    {
      misc.error("Error: An internal error was happened. Inserting a matrix with no block position properly specified in a blocked covariance matrix.", 0);
    }
  }
  
  newElement.blockPosition = blockPosition;
  newElement.outcomeSubMatrix = outcomeSubMatrix;
  newElement.sourceSubMatrix = sourceSubMatrix;
  
  this->elements.push_back(newElement);
}

void CovarianceMatrix::appendVarianceToElement(std::string elementName, std::string varianceName, VarianceType varianceType)
{
  int idxElement = getElement(elementName);
  if(idxElement < 0)
  {
    misc.error("Error: Variance cannot be appended to a non existing element.", 0);
  }
  for(int j = 0; j<this->variances.size(); j++)
  {
    if(this->variances[j].name == varianceName)
    {
      this->elements[idxElement].iVariances.push_back(j);
      this->elements[idxElement].nameVariances.push_back(varianceName);
      this->elements[idxElement].varianceTypes.push_back( varianceType );
      this->variances[j].inElements.insert(idxElement);
      return;
    }
  }
  misc.error("Error: Invalid variance name when adding variances to an Element: " + varianceName, 0);
}

int CovarianceMatrix::getElement(std::string elementName)
{
  for(int i = 0; i<this->elements.size(); i++)
  {
    if(this->elements[i].name == elementName)
    {
      return i;
    }
  }
  return -1;
}

int CovarianceMatrix::getVariance(std::string varianceName)
{
  for(int i = 0; i<this->variances.size(); i++)
  {
    if(this->variances[i].name == varianceName)
    {
      return i;
    }
  }
  misc.error("Error: An internal error was happened. Malformed Covariance Matrix. The variance name do no exist.", 0);
  return -1;
}

void CovarianceMatrix::reinitializeVariances()
{
  for(int i = 0; i<this->variances.size(); i++)
  {
    this->variances[i].variance = this->variances[i].initialVariance;
  }
  broadcastVariances();
}

void CovarianceMatrix::clearElementVariances(std::string elementName)
{
  int idxElement = getElement(elementName);
  if(idxElement < 0)
  {
    misc.error("Error: An internal error was happened. Inexistent element when trying to clear element.", 0);
  }
  for(int i = 0; i<this->elements[idxElement].iVariances.size(); i++)
  {
    int idxVariance = this->elements[idxElement].iVariances[i];
    int removed = this->variances[idxVariance].inElements.erase(idxElement);
    if(removed != 1)
    {
      misc.error("Error: An internal error was happened. Malformed Covariance Matrix.", 0);
    }
  }
  this->elements[idxElement].iVariances.clear();
  this->elements[idxElement].nameVariances.clear();
  this->elements[idxElement].varianceTypes.clear();
}

void CovarianceMatrix::changeElementConstantFactor(std::string elementName, double constantFactor)
{
  int idxElement = getElement(elementName);
  if(idxElement < 0)
  {
    misc.error("Error: An internal error was happened. Inexistent element when trying to change the constant factor.", 0);
  }
  this->elements[idxElement].constantFactor = constantFactor;
}

void CovarianceMatrix::deleteVariance(std::string name)
{
  int idxVariance = -1;
  for(int i = 0; i<this->variances.size(); i++)
  {
    if(this->variances[i].name == name)
    {
      idxVariance = i;
    }
  }
  if(idxVariance < 0)
  {
    misc.error("Error: An internal error was happened. Invalid variance name when adding variances to an Element: " + name, 0);
  }
  if(this->variances[idxVariance].inElements.size() != 0)
  {
    misc.error("Error: An internal error was happened. Trying to remove a variance before removing it from elements.", 0);
  }
  
  std::string group = this->variances[idxVariance].group;
  this->variances.erase(this->variances.begin() + idxVariance);
  
  bool flag = false;
  for(int i = 0; i<this->variances.size(); i++)
  {
    if(this->variances[i].group == group)
    {
      flag = true;
    }
  }
  if(flag == false)
  {
    int erased = varianceGroupExpectedMagnitudes.erase(group);
    if(erased != 1)
    {
      misc.error("Error: An internal error was happened. Malformed Covariance Matrix.", 0);
    }
  }
  
  for(int i = 0; i < this->elements.size(); i++)
  {
    for(int j = 0; j < this->elements[i].nameVariances.size(); j++)
    {
      if( (this->elements[i].nameVariances[j] == name) || (this->elements[i].nameVariances.size() != this->elements[i].iVariances.size()) )
      {
        misc.error("Error: An internal error was happened. Malformed Covariance Matrix.", 0);
      }
      this->elements[i].iVariances[j] = getVariance(this->elements[i].nameVariances[j]);
    }
  }
}

void CovarianceMatrix::computeCovariance()
{
  this->m->fillWithConstant(0.);
  this->m->symmetric = true;
  this->m->uplo = 'B';
  
  for(int i=0; i<this->elements.size(); i++)
  {
    Element currentElement = this->elements[i];
    
    //Compute the product of the variance elements
    double currentVariance = computeElementVariance(currentElement);
    
    //Add to the covariance a subcovariance matrix multiplied by the previous product of variances.
    if(currentElement.outcomeSubMatrix.active == false || currentElement.sourceSubMatrix.active == false)
    {
      this->m->add(currentElement.m, (i==0)?0.:1., currentVariance);
    }
    else
    {
      this->m->add(currentElement.m, (i==0)?0.:1., currentVariance, currentElement.outcomeSubMatrix, currentElement.sourceSubMatrix);
    }
  }
  this->m->symmetric = true;
  this->m->uplo = 'U';
  this->mStatus = covarianceMatrix;
}

void CovarianceMatrix::computeBlockCovariance()
{
  if( this->byBlocks == false )
  {
    misc.error("Error: An internal error was happened. Covariance matrix cannot be computed by blocks. Blocks undefined.", 0);
  }
  
  int nrBlocks = 0;
  int ncBlocks = 0;
  
  std::vector< std::vector<Matrix*> > blocks = std::vector< std::vector<Matrix*> >(this->blockDimension, std::vector<Matrix*>(this->blockDimension, NULL));
    
  for(int i=0; i<this->elements.size(); i++)
  {
    Element currentElement = this->elements[i];
    
    //Compute the product of the variance elements
    double currentVariance = computeElementVariance(currentElement);
    
    int br = currentElement.blockPosition.first;
    int bc = currentElement.blockPosition.second;
    
    if( br > bc )
    {
      misc.error("Error: An internal error was happened wen computing the covariance matrix in block form. No elements expected below the diagonal.", 0);
    }
    
    //Add to the covariance a subcovariance matrix multiplied by the previous product of variances.
    if(blocks[br][bc] == NULL)
    {
      Matrix * temp = new Matrix(currentElement.m->distribution, currentElement.m->nGlobRows, currentElement.m->nGlobCols);
      temp->fillWithConstant(0.);
      blocks[br][bc] = temp;
    }
    blocks[br][bc]->add(currentElement.m, 1., currentVariance);
  }
  
  this->mInBlocks.clear();
  
  for(int br=0; br<this->blockDimension; br++)
  {
    std::vector<Matrix*> rowBlocks;
    for(int bc=0; bc<this->blockDimension; bc++)
    {
      if( bc >= br && blocks[br][bc] == NULL )
      {
        misc.error("Error: An internal error was happened wen computing the covariance matrix in block form. Empty block.", 0);
      }
      if( bc < br)
      {
        if(blocks[br][bc] != NULL)
        {
          misc.error("Error: An internal error was happened wen computing the covariance matrix in block form. Empty block.", 0);
        }
        blocks[br][bc] = new Matrix(blocks[bc][br]);
      }
      rowBlocks.push_back(blocks[br][bc]);
    }
    this->mInBlocks.addBlockRow(rowBlocks);
  }
  blocks.clear();
  
  this->mStatus = covarianceMatrixByBlocks;
}

Matrix * CovarianceMatrix::computeDerivateCovariance(int idxDerivateVariance)
{
  bool flag = false;
  
  if(variances[idxDerivateVariance].inElements.size() == 0)
  {
    misc.error("Error: An internal error was happened on method computeDerivateCovariance().", 0);
  }
  
  this->mDerivate->fillWithConstant(0.);
  this->mDerivate->symmetric = true;
  this->mDerivate->uplo = 'B';
  
  //Iterate over the elements that has the variance with index idxDerivateVariance in their components
  int i = 0;
  for (std::set<int>::iterator it=this->variances[idxDerivateVariance].inElements.begin(); it!=this->variances[idxDerivateVariance].inElements.end(); ++it)
  {
    flag = true;
    Element currentElement = this->elements[ *it ];
    
    //Compute the derivate of the current element variances product with respect the variance idxDerivateVariance.
    double currentVariance = computeElementVarianceDerivate(currentElement, idxDerivateVariance);
    
    if(currentElement.outcomeSubMatrix.active == false || currentElement.sourceSubMatrix.active == false)
    {
      this->mDerivate->add(currentElement.m, (i==0)?0.:1., currentVariance);
    }
    else
    {
      this->mDerivate->add(currentElement.m, (i==0)?0.:1., currentVariance, currentElement.outcomeSubMatrix, currentElement.sourceSubMatrix);
    }
    i++;
  }
  
  if(flag == false)
  {
    misc.error("Error: An internal error was happened on method computeDerivateCovariance().", 0);
  }
  
  this->mDerivate->symmetric = true;
  this->mDerivate->uplo = 'U';
  
  return this->mDerivate;
}

double CovarianceMatrix::computeElementVariance(Element & element)
{
  if(element.iVariances.size() == 0)
  {
    misc.error("Error: An internal error was happened. Covariance matrix malformed. There is an element without variances.", 0);
  }

  double currentVariance = 1.;
  for(int j = 0; j<element.iVariances.size(); j++)
  {
    if(element.varianceTypes[j] == variance || element.varianceTypes[j] == covariance)
    {
      currentVariance *= this->variances[element.iVariances[j]].variance;
    }
    else if(element.varianceTypes[j] == standardDeviation)
    {
      currentVariance *= sqrt(this->variances[element.iVariances[j]].variance); //The sqrt if this variance is in fact a standard deviation.
    }
    else
    {
      misc.error("Error: An internal error was happened", 0);
    }
  }
  return currentVariance*element.constantFactor;
}

double CovarianceMatrix::computeSubCovarianceVariance(std::string subCovarianceId)
{
  double result = 0.;
  int nElements = 0;
  for(int i=0; i<this->elements.size(); i++)
  {
    if( this->elements[i].subCovarianceId == subCovarianceId )
    {
      result = this->computeElementVariance(this->elements[i]);
      nElements++;
    }
  }
  if(nElements == 0)
  {
    misc.error("Error: An internal error was happened. Computing the variance of a subCovariance with inexistent name.", 0);
  }
  if(nElements != 1)
  {
    misc.error("Error: An internal error was happened. Computing the variance of a subCovariance with more than one element is not allowed.", 0);
  }
  return result;
}

double CovarianceMatrix::computeElementVarianceDerivate(Element & element, int idxDerivateVariance)
{
  bool flag = false;
  double currentVariance = 1.;
  for(int j = 0; j<element.iVariances.size(); j++) //Creates the derivate of the product of variances with respect the variance idxDerivateVariance
  {
    int idx  = element.iVariances[j];
    if(idx != idxDerivateVariance)
    {
      if(element.varianceTypes[j] == variance || element.varianceTypes[j] == covariance)
      {
        currentVariance *= this->variances[idx].variance;
      }
      else if(element.varianceTypes[j] == standardDeviation)
      {
        currentVariance *= sqrt(this->variances[idx].variance);
      }
      else
      {
        misc.error("Error: An internal error was happened.", 0);
      }
    }
    else
    {
      flag = true; //Just check everything is ok.
      if(element.varianceTypes[j] == variance || element.varianceTypes[j] == covariance)
      {
        currentVariance *= 1.;
      }
      else if(element.varianceTypes[j] == standardDeviation)
      {
        currentVariance *= 0.5/sqrt(this->variances[idx].variance); //If the product depends not on the variance but on the standard deviation, the derivate is 1/(2*std) instead of 1.
      }
      else
      {
        misc.error("Error: An internal error was happened.", 0);
      }
    }
  }
  if(!flag)
  {
    misc.error("Error: An internal error was happened when computing a variance derivate.", 0);
  }
  return currentVariance*element.constantFactor;
}

bool CovarianceMatrix::invertCovariance(double * logDeterminant, bool useSinglePrecision)
{
  if( this->allCovarianceSubMatricesDiagonal == false || this->byBlocks == false || this->blockDimension <= 1)
  {
    if(this->mStatus != covarianceMatrix)
    {
      computeCovariance();
    }
    bool inverted = this->m->symmetricInvert(logDeterminant, useSinglePrecision);
    if(inverted == false)
    {
      computeCovariance();
      inverted = this->m->invert(logDeterminant);
      if(inverted == false)
      {
        misc.message << "Error: V matrix can not be inverted. REML iterations cannot continue." << std::endl;
        return false;
      }
    }
    this->mStatus = inverseCovarianceMatrix;
    return true;
  }
  else //When all matrices are diagonal and block dimension is gretaher than one, use blockMatrix for inversion.
  {
    computeBlockCovariance();
    bool inverted = this->mInBlocks.invert(logDeterminant);
    if(inverted == false)
    {
      misc.message << "Error: V matrix can not be inverted. REML iterations cannot continue." << std::endl;
      return false;
    }
    
    this->mStatus = inverseCovarianceMatrixByBlocks;
    
    return true;
  }
}

int CovarianceMatrix::constrainVariancesM1()
{
  int nVariancesConstrained = 0;
  std::map<std::string, int> nVariancesConstrainedGroup;
  std::map<std::string, double> correctedMagnitude;
  std::map<std::string, double> nonCorrectedMagnitude;
  for(std::map<std::string, double>::iterator it = varianceGroupExpectedMagnitudes.begin(); it != varianceGroupExpectedMagnitudes.end(); ++it)
  {
    nVariancesConstrainedGroup[it->first] = 0;
    correctedMagnitude[it->first] = 0.;
    nonCorrectedMagnitude[it->first] = 0.;
  }
  std::vector<bool> constrained(this->variances.size(), false);
  
  for(int i = 0; i < this->variances.size(); i++)
  {
    std::string idGroup = this->variances[i].group;
    double baseMagnitude = this->varianceGroupExpectedMagnitudes[idGroup];
    
    if(this->variances[i].variance < 0 && this->variances[i].type == variance)
    {
      correctedMagnitude[idGroup] += baseMagnitude*options.varianceConstrainProportion - this->variances[i].variance;
      this->variances[i].variance = baseMagnitude*options.varianceConstrainProportion;
      constrained[i] = true;
      nVariancesConstrainedGroup[idGroup]++;
      nVariancesConstrained++;
    }
    else if(this->variances[i].variance > 0 && this->variances[i].type == variance)
    {
      nonCorrectedMagnitude[idGroup] += this->variances[i].variance;
    }
  }
  
  //Use GCTA methods?
  if(options.remlGCTAMode)
  {
    for(std::map<std::string, double>::iterator it = varianceGroupExpectedMagnitudes.begin(); it != varianceGroupExpectedMagnitudes.end(); ++it)
    {
      if(nVariancesConstrainedGroup[it->first] != 0)
      {
        correctedMagnitude[it->first] /= double(nVariancesConstrainedGroup[it->first]);
      }
    }
    
    for(int i = 0; i < this->variances.size(); i++)
    {
      std::string idGroup = this->variances[i].group;
      double baseMagnitude = this->varianceGroupExpectedMagnitudes[idGroup];
      if(constrained[i] == false && this->variances[i].type == variance)
      {
        this->variances[i].variance -= correctedMagnitude[idGroup];
        /*this->variances[i].variance -= correctedMagnitude[idGroup]*this->variances[i].variance/nonCorrectedMagnitude[idGroup];
        if(this->variances[i].variance < 0.)
        {
          this->variances[i].variance = baseMagnitude*options.varianceConstrainProportion;
        }*/
      }
    }
  }
  
  communicator->broadcast(&nVariancesConstrained, 1); 
  broadcastVariances();
  
  options.varianceConstrainProportion = options.varianceConstrainProportion/10.;
  
  return nVariancesConstrained;
}

double CovarianceMatrix::constrainVariancesM2(std::vector<double> & oldVariances, std::vector<double> & delta)
{
  int nVariancesToConstrain = 0;
  std::vector<double> scalingFactors;
  std::vector<bool> infiniteScalingFactor;
  
  if(delta.size() != this->variances.size())
  {
    misc.error("Error: An internal error was happened. The dimensions of variances and delta differ.", 0);
  }
  
  for(int i = 0; i < this->variances.size(); i++)
  {
    std::string idGroup = this->variances[i].group;
    double baseMagnitude = this->varianceGroupExpectedMagnitudes[idGroup];
    infiniteScalingFactor.push_back(false);
    if((oldVariances[i] + delta[i]) < 0 && this->variances[i].type == variance)
    {
      double tempScalingFactor = -delta[i]/(oldVariances[i] - baseMagnitude*options.varianceConstrainProportion);
      if( isinf(tempScalingFactor) || isnan(tempScalingFactor) )
      {
        infiniteScalingFactor[i] = true;
      }
      else
      {
        scalingFactors.push_back( tempScalingFactor );
      }
      nVariancesToConstrain++;
    }
  }
  
  if(nVariancesToConstrain == 0)
  {
    return 1.;
  }
  
  if(scalingFactors.size() == 0)
  {
    misc.error("Error: When constraining variances, an error was happened. Too many negative variances and a scaling factor can not be computed.", 0);
  }
  
  //Search the minimum scaling factor that don't produce negative variances.
  double minimumScalingFactor = -1.;
  bool flagCheck = false; //Just for checking
  for(int i = 0; i < scalingFactors.size(); i++)
  {
    int flag = true;
    for(int j = 0; j < this->variances.size(); j++)
    {
      std::string idGroup = this->variances[j].group;
      double baseMagnitude = this->varianceGroupExpectedMagnitudes[idGroup];
      if((oldVariances[j] + delta[j]) < 0 && this->variances[j].type == variance)
      {
        double temp = oldVariances[j] + (delta[j]/scalingFactors[i]);
        if(temp < 0. && infiniteScalingFactor[j] == false)
        {
          flag = false;
        }
      }
    }
    if(flag && ( (scalingFactors[i] <= minimumScalingFactor) || minimumScalingFactor < 0. ))
    {
      minimumScalingFactor = scalingFactors[i];
      flagCheck = true;
    }
  }
  if(!flagCheck)
  {
    misc.error("Error: An internal error was happened when constraining variances (Method 2).", 0);
  }
  
  for(int i = 0; i < this->variances.size(); i++)
  {
    std::string idGroup = this->variances[i].group;
    double baseMagnitude = this->varianceGroupExpectedMagnitudes[idGroup];
    this->variances[i].variance = oldVariances[i] + (delta[i]/minimumScalingFactor);
    if(this->variances[i].variance < 0 && this->variances[i].type == variance)
    {
      if(infiniteScalingFactor[i] == true)
      {
        this->variances[i].variance = baseMagnitude*options.varianceConstrainProportion;
      }
      else
      {
        misc.error("Error: An internal error was happened when constraining variances (Method 2). Variances still negative.", 0);
      }
    }
  }
  
  broadcastVariances();
  
  return minimumScalingFactor;
}

int CovarianceMatrix::fixVariancesToZero()
{
  int nVariancesFixed = 0;
  for(int i = 0; i < this->variances.size(); i++)
  {
    std::string idGroup = this->variances[i].group;
    double baseMagnitude = this->varianceGroupExpectedMagnitudes[idGroup];
    if((this->variances[i].variance*(1e20)) < baseMagnitude && this->variances[i].type == variance)
    {
      this->variances[i].variance = 0.;
      nVariancesFixed++;
    }
  }
  
  broadcastVariances();
  communicator->broadcast(&nVariancesFixed);
  
  return nVariancesFixed;
}

std::string CovarianceMatrix::getStringVarianceNames(std::string separator)
{
  std::stringstream result;
  for(int i = 0;  i< this->variances.size(); i++)
  {
    if(options.logFieldWidth > this->variances[i].name.size())
    {
      result << std::setw(options.logFieldWidth) << this->variances[i].name;
    }
    else
    {
      result << (i==0?"":separator) << this->variances[i].name;
    }
  }
  return result.str();
}

std::string CovarianceMatrix::getStringVarianceValues(std::string separator)
{
  std::stringstream result;
  for(int i = 0;  i< this->variances.size(); i++)
  {
    result << std::setprecision(options.logOutputPrecision) << std::setw(options.logFieldWidth) << this->variances[i].variance;
  }
  return result.str();
}

std::string CovarianceMatrix::getMatrixName(Matrix *m)
{
  for(int i = 0; i<this->covarianceSubMatrices.size(); i++)
  {
    if(this->covarianceSubMatrices[i] == m)
    {
      return this->covarianceSubMatrixNames[i];
    }
  }
  misc.error("Error: An internal error was happened. In covariance matrix: Matrix not found.", 0);
}

void CovarianceMatrix::broadcastVariances()
{
  double * vv = new double [this->variances.size()];
  
  for(int i = 0; i < this->variances.size(); i++)
  {
    vv[i] = this->variances[i].variance;
  }
  communicator->broadcast(vv, this->variances.size());
  for(int i = 0; i < this->variances.size(); i++)
  {
    this->variances[i].variance = vv[i];
  }
  
  delete [] vv;
}

void CovarianceMatrix::broadcastInitialVariances()
{
  double * vv = new double [this->variances.size()];
  
  for(int i = 0; i < this->variances.size(); i++)
  {
    vv[i] = this->variances[i].initialVariance;
  }
  communicator->broadcast(vv, this->variances.size());
  for(int i = 0; i < this->variances.size(); i++)
  {
    this->variances[i].initialVariance = vv[i];
  }
  
  delete [] vv;
}

std::vector<double> CovarianceMatrix::getVectorVariances()
{
  std::vector<double> vv;
  for(int i = 0; i < this->variances.size(); i++)
  {
    vv.push_back(this->variances[i].variance);
  }
  return vv;
}

void CovarianceMatrix::storeVectorVariances(std::vector<double> & vv)
{
  if(vv.size() != this->variances.size())
  {
    misc.error("Error: An internal error was happened. The size of vector variances is not equal to the argument vector.", 0);
  }
  for(int i = 0; i < this->variances.size(); i++)
  {
    this->variances[i].variance = vv[i];
  }
}

void CovarianceMatrix::getElementMatrixInContext(int i, Matrix * m1)
{
  m1->fillWithConstant(0.);
  
  Element currentElement = this->elements[i];
  
  //Add to the covariance a subcovariance matrix multiplied by the previous product of variances.
  if(currentElement.outcomeSubMatrix.active == false || currentElement.sourceSubMatrix.active == false)
  {
    m1->add(currentElement.m, 0., 1.);
  }
  else
  {
    m1->add(currentElement.m, 0., 1., currentElement.outcomeSubMatrix, currentElement.sourceSubMatrix);
  }
  
  m1->symmetric = true;
  m1->uplo = 'U';
}

Matrix * CovarianceMatrix::getSubCovariance(std::string subCovarianceId)
{
  Matrix * subCovarianceMatrix = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->dimension, this->dimension, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  
  subCovarianceMatrix->fillWithConstant(0.);
  subCovarianceMatrix->symmetric = true;
  subCovarianceMatrix->uplo = 'B';
  
  bool nameExist = false;
  for(int i=0; i<this->elements.size(); i++)
  {
    if(elements[i].subCovarianceId != subCovarianceId)
    {
      continue;
    }
    nameExist = true;
    
    Element currentElement = this->elements[i];
    
    //Compute the product of the variance elements
    double currentVariance = computeElementVariance(currentElement);
    
    //Add to the covariance a subcovariance matrix multiplied by the previous product of variances.
    if(currentElement.outcomeSubMatrix.active == false || currentElement.sourceSubMatrix.active == false)
    {
      subCovarianceMatrix->add(currentElement.m, (i==0)?0.:1., currentVariance);
    }
    else
    {
      subCovarianceMatrix->add(currentElement.m, (i==0)?0.:1., currentVariance, currentElement.outcomeSubMatrix, currentElement.sourceSubMatrix);
    }
  }
  
  if(nameExist == false)
  {
    misc.error("Error: An internal error was happened. Inexistent subCovarianceId.", 0);
  }
  
  subCovarianceMatrix->symmetric = true;
  subCovarianceMatrix->uplo = 'U';
  
  return subCovarianceMatrix;
}

Matrix * CovarianceMatrix::multiply(std::string subCovarianceId, Matrix * m, double scale)
{
  Matrix * mResult;
  if( this->allCovarianceSubMatricesDiagonal == false )
  {
    mResult = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
    Matrix * subCovariance = getSubCovariance(subCovarianceId);
    mResult->multiply(subCovariance, 'N', m, 'N', scale);
    delete subCovariance;
  }
  else
  {
    mResult = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
    Matrix * subCovariance = getSubCovariance(subCovarianceId);
    mResult->multiply(subCovariance, 'N', m, 'N', scale);
    delete subCovariance;
  }
  return mResult;
}

void CovarianceMatrix::setVarianceInitialValues(std::vector<Variance> newVariances)
{
  int noInitialize;
  if(communicator->mpiRoot && newVariances.size() == 0)
  {
    noInitialize = 1;
  }
  communicator->broadcast(&noInitialize);
  if(noInitialize == 1)
  {
    return;
  }
  
  if(communicator->mpiRoot)
  {
    if(newVariances.size() != this->variances.size())
    {
      misc.error("Error: The defined initial variances do not agree with REML variances.", 0);
    }
    for(int i = 0; i < newVariances.size(); i++)
    {
      int check = 0;
      for(int j = 0; j < this->variances.size(); j++)
      {
        if(newVariances[i].name == this->variances[j].name)
        {
          this->variances[j].initialVariance = newVariances[i].variance;
          check++;
        }
      }
      if(check != 1)
      {
        misc.error("Error: Error when setting initial variances. Variance names do not agree.", 0);
      }
    }
  }
  broadcastInitialVariances();
}

void CovarianceMatrix::showCovarianceMatrix(bool showMatrices)
{
  if(showMatrices)
  {
    Matrix * temp = new Matrix(cyclicDistribution, this->dimension, this->dimension, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
    for(int i=0; i<this->elements.size(); i++)
    {
      Element currentElement = this->elements[i];
      temp->fillWithConstant(0.);
      temp->add(currentElement.m, 0., 1, currentElement.outcomeSubMatrix, currentElement.sourceSubMatrix);
      temp->showGlobal(getMatrixName(currentElement.m));
    }
    delete temp;
  }
 
  misc.message << "*********************************************" << std::endl;
  misc.message << getStringVarianceNames() << std::endl;
  misc.message << getStringVarianceValues() << std::endl << std::endl;
  
  for(int i=0; i<this->elements.size(); i++)
  {
    Element currentElement = this->elements[i];
    
    std::string mod1 = "";
    std::string mod2 = "";
    
    
    misc.message << currentElement.constantFactor << "*M(" << getMatrixName(currentElement.m) << ")M ";
    for(int j = 0; j<currentElement.iVariances.size(); j++)
    {
      if(currentElement.varianceTypes[j] == standardDeviation)
      {
        mod1 = "sqrt(";
        mod2 = ")";
      }
      else if(currentElement.varianceTypes[j] == covariance)
      {
        mod1 = "c(";
        mod2 = ")";
      }
      else
      {
        mod1 = "";
        mod2 = "";
      }
      misc.message << "*" << mod1 << this->variances[currentElement.iVariances[j]].name << mod2;
    }
    misc.message << ((i!=(elements.size()-1))?" + ":"");
  }
  misc.message << std::endl;
  misc.message << "*********************************************" << std::endl;
}
