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
#include "options.h"
#include "communicator.h"
#include "misc.h"
#include "global.h"
#include "reml.h"
#include "grm.h"
#include "covariate.h"
#include "genotype.h"
#include "phenotype.h"
#include "covariancematrix.h"
#include "auxiliar.h"
#include "message.h"

#include <fstream>
#include <iomanip>
#include <cmath>

REML::REML(bool argWriteResults)
{
  this->y = NULL;
  this->X = NULL;
  
  this->V = NULL;
  
  this->ViX = NULL;
  this->XtViX_i = NULL;
  this->P = NULL;
  this->Py = NULL;
  this->subVPy = NULL;
  this->AI = NULL;
  this->yPsubVPy_trPsubV = NULL;
  
  this->success = false;
  this->writeResults = argWriteResults;
  
  this->singlePrecisionInversion = options.allowSinglePrecisionInversion;
  
  this->nPhenotypes = 0;
  
  this->usingDiagonalKernels = false;
}

REML::~REML()
{
  if(this->y != NULL)
  {
    delete this->y;
  }
  if(this->X != NULL)
  {
    delete this->X;
  }
  if(this->V != NULL)
  {
    delete this->V;
  }
  for(std::map<std::string, Genotype *>::iterator it = this->SNPsBLUPGenotypes.begin(); it != this->SNPsBLUPGenotypes.end(); ++it)
  {
    delete it->second;
  }
  this->SNPsBLUPGenotypes.clear();
  
  for(std::map<std::string, Matrix *>::iterator it = this->GRMEigenVectors.begin(); it != this->GRMEigenVectors.end(); ++it)
  {
    delete it->second;
  }
  this->GRMEigenVectors.clear();
  
  deleteIntermediateMatrices();
}

void REML::deleteIntermediateMatrices()
{
  if( this->ViX != NULL )
  {
    delete this->ViX;
    this->ViX = NULL;
  }
  if( this->XtViX_i != NULL )
  {
    delete this->XtViX_i;
    this->XtViX_i = NULL;
  }
  if( this->P != NULL )
  {
    delete this->P;
    this->P = NULL;
  }
  if( this->Py != NULL )
  {
    delete this->Py;
    this->Py = NULL;
  }
  if( this->subVPy != NULL )
  {
    delete this->subVPy;
    this->subVPy = NULL;
  }
  if( this->AI != NULL )
  {
    delete this->AI;
    this->AI = NULL;
  }
  if( this->yPsubVPy_trPsubV != NULL )
  {
    delete this->yPsubVPy_trPsubV;
    this->yPsubVPy_trPsubV = NULL;
  }
}

bool REML::prepare(REMLType type, std::vector<Kernel*> & kernels, std::vector<double> weights, std::vector<int> phenotypeColumns, std::vector<double> heritabilities, std::vector<std::pair<std::string, std::string> > covariateFiles)
{
  //Perform some checks and parameter initialization
  
  if(kernels.size() < 1)
  {
    misc.error("Error: An internal error was happened. No GRMs specified for REML analysis.", 0);
  }
  
  if(phenotypeColumns.size() < 1)
  {
    misc.error("Error: An internal error was happened. No Phenotypes specified for REML analysis.", 0);
  }
  
  if(covariateFiles.size() != phenotypeColumns.size())
  {
    misc.error("Error: An internal error was happened. Different number of covariate files than phenotypes columns to analyze specified for REML analysis.", 0);
  }  

  if(weights.size() == 0)
  {
    double equallyDistributedWeight = 1./double(kernels.size());
    for(int i = 0; i<kernels.size(); i++)
    {
      weights.push_back(equallyDistributedWeight);
    }
  }
  else
  {
    if( kernels.size() != weights.size() )
    {
      misc.error("Error: An internal error was happened. Kernel weights not properly defined.", 0);
    }
  }
  
  if(heritabilities.size() == 0)
  {
    for(int i = 0; i<phenotypeColumns.size(); i++)
    {
      heritabilities.push_back(0.5);
    }
  }
  else
  {
    if( heritabilities.size() != phenotypeColumns.size() )
    {
      misc.error("Error: An internal error was happened. Initial heritabilities not properly defined.", 0);
    }
  }
  
  //Set reml type
  
  this->type = type;
  
  //Set the number of phenotypes
  
  this->nPhenotypes = phenotypeColumns.size();
  
  //Sanitize kernels
  
  for(int i = 0; i<kernels.size(); i++)
  {
    bool sanitized = kernels[i]->sanitizeKernel();
    if( sanitized == false )
    {
      return false;
    }
  }
  
  //Prepare the GRM (covs), covariance (X) and phenotype (y) matrices
  
  std::vector<Phenotype*> phenotypes;
  for(int i = 0; i<this->nPhenotypes; i++)
  {
    Phenotype * temp = new Phenotype(cyclicDistribution, options.phenotypesFile, phenotypeColumns[i]);
    phenotypes.push_back(temp);
  }
  
  this->covariateNames.clear();
  
  std::vector<Covariate*> covariates;
  int previousIndex = 0;
  for(int i = 0; i<covariateFiles.size(); i++)
  {
    if(options.joinCovariatesVertically == true) //This option is untested
    {
      Covariate * tempCovariate = new Covariate(covariateFiles[i].first, covariateFiles[i].second, phenotypes[i]->individualIds, false);
      for(int j = 0; j<covariates.size(); j++)
      {
        covariates[j]->syncronizeDiscreteCovariateCategoriesWith(tempCovariate);
      }
      covariates.push_back(tempCovariate);
    }
    else
    {
      Covariate * tempCovariate = new Covariate(covariateFiles[i].first, covariateFiles[i].second, phenotypes[i]->individualIds);
      covariates.push_back(tempCovariate);
      this->addCovariatesNames(tempCovariate, addSuffix("p" + i2s(i + 1) + "-"), previousIndex);
      previousIndex += tempCovariate->covariates->nGlobCols;
    }
  }
  
  if(options.joinCovariatesVertically == true) //This option is untested
  {
    for(int i = 0; i<covariates.size(); i++)
    {
      covariates[i]->parseRawCovariates(phenotypes[i]->individualIds, covariates.size(), i);
    }
    this->addCovariatesNames(covariates[0]);
  }
  
  //Search for common individuals
  
  std::vector< std::vector<std::string> > commonIndividualsInGRMOrder(this->nPhenotypes, std::vector<std::string>() );
  std::vector< std::vector<int> > idxsKeptKernel0(this->nPhenotypes, std::vector<int>() );
  int totalIndividuals;
  if(communicator->mpiRoot)
  {
    for(int i = 0; i<this->nPhenotypes; i++)
    {
      std::vector<std::string> tempIndividuals = intersectionStringVectors(3, &(kernels[0]->individualIds), &phenotypes[i]->individualIds, &covariates[i]->individualIds);
      for(int j = 1; j<kernels.size(); j++)
      {
        tempIndividuals = intersectionStringVectors(2, &tempIndividuals, &(kernels[j]->individualIds));
      }
      
      std::vector<std::string> tempOrderedIndividuals = orderVectorAsTemplate(kernels[0]->individualIds, tempIndividuals);
      commonIndividualsInGRMOrder[i] = tempOrderedIndividuals;
      for(int j = 1; j<kernels.size(); j++)
      {
        std::vector<std::string> tempTest = orderVectorAsTemplate(kernels[j]->individualIds, tempIndividuals);
        if( tempTest != tempOrderedIndividuals )
        {
          misc.error("Error: Not all GRMs/Kernels have the individuals in the same order. In the current version, DISSECT cannot resort it. Please, resort it or contact us.", 0);
        }
      }
      if( tempOrderedIndividuals.size() < 1)
      {
        misc.error("Error: There are not enough individuals for one of the traits for performing the analysis.", 0);
      }
      misc.message << tempOrderedIndividuals.size() << " individuals are kept for the analysis for phenotype " << i + 1 << "." << std::endl;
      
      std::vector<int> tempIdxsKept = extractMapValues(tempOrderedIndividuals, kernels[0]->individualIdsIdx);
      idxsKeptKernel0[i] = tempIdxsKept;
    }
  }

  int nIndividualsBeforeFilterKernel0 = kernels[0]->nIndividuals;
  std::vector<int> nIndividualsTraits;
  int nTotalIndividuals = 0;
  for(int i = 0; i<this->nPhenotypes; i++)
  {
    int nIndividualsTrait = commonIndividualsInGRMOrder[i].size();
    communicator->broadcast(&nIndividualsTrait, 1);
    nIndividualsTraits.push_back(nIndividualsTrait);
    
    nTotalIndividuals += nIndividualsTrait;
  }
  
  ////////////////////////////////////////////////////
  // Init phenotype matrix
  
  std::vector<double> phenotypeVariances;
  for(int i = 0; i<this->nPhenotypes; i++)
  {
    phenotypes[i]->filterIndividuals(commonIndividualsInGRMOrder[i]);
    double phenotypeVariance = phenotypes[i]->computePhenotypeVariance();
    phenotypeVariances.push_back(phenotypeVariance);
    if( i == 0)
    {
      this->y = new Matrix(phenotypes[i]->phenotypes);
    }
    else
    {
      Matrix * temp = new Matrix(this->y);
      this->y->joinMatricesVertically(temp, phenotypes[i]->phenotypes);
      delete temp;
    }
    delete phenotypes[i];
  }
  phenotypes.clear();
  
  
  ////////////////////////////////////////////////////
  // Check whether Kernel diagonalization is coherent. If all kernels are diagonal, they must have the same individuals and this must be shared with phenotypes and covariates.
  // If only some kernels are diagonal it is an error if options.forceUseDiagonalizedKernels == false (default option). Kernels are de-diagonalized otherwise.
  this->usingDiagonalKernels = false;
  bool allKernelsDiagonal = true;
  bool someKernelsDiagonal = false;
  for(int i = 0; i<kernels.size(); i++)
  {
    if(kernels[i]->diagonalized == true)
    {
      someKernelsDiagonal = true;
    }
    else
    {
      allKernelsDiagonal = false;
    }
  }
  
  if( allKernelsDiagonal == false && someKernelsDiagonal == true )
  {
    if(options.forceUseDiagonalizedKernels == true)
    {
      for(int i = 0; i<kernels.size(); i++)
      {
        kernels[i]->recoverKernelFromEigenDecomposition();
      }
    }
    else
    {
      misc.error("Error: Sorry, this analysis cannot be performed with diagonal GRMs/Kernels. You can use the option --force-use-diag-kernels to force converting diagonal GRMs/Kernels to their non-diagonalized form before starting the analysis. However, this have some limitations: any analysis that depends on the normalization matrix cannot be performed from diagonal Kernels/GRMs (e.g. regional analysis)", 0);
    }
  }
  else if( allKernelsDiagonal == true )
  {
    for(int i = 0; i<kernels.size(); i++)
    {
      if(kernels[i]->individualIds != commonIndividualsInGRMOrder[i] || kernels[i]->individualIds != kernels[0]->individualIds)
      {
        misc.error("Error: Sorry, this analysis cannot be performed with diagonal GRMs/Kernels. All kernels have to share same individuals and these also have to be present on phenotype and covariate data.", 0);
      }
    }
    
    if(kernels.size() == 1 && this->nPhenotypes == 1)
    {
      this->GRMEigenVectors[kernels[0]->name] = new Matrix(kernels[0]->eigenVectors);
      this->usingDiagonalKernels = true;
      
      transformUsingEigenVectors(kernels[0]->name, 'T', &(this->y), this->nPhenotypes);
      
      if(this->SNPsBLUPGenotypes.count(kernels[0]->name) != 0)
      {
        this->SNPsBLUPGenotypes[ kernels[0]->name ]->normalizeGenotypes();
        this->SNPsBLUPGenotypes[ kernels[0]->name ]->filterSNPsAndIndividuals(this->SNPsBLUPGenotypes[ kernels[0]->name ]->SNPIds, kernels[0]->individualIds);
        if( this->SNPsBLUPGenotypes[ kernels[0]->name ]->individualIds != kernels[0]->individualIds )
        {
          misc.error("Error: The order of individuals in the Kernel/GRM is different that the order in the genotypes file. Sorry, at this version DISSECT needs individuals in both files must be in the same order in both files.", 0);
        }
        Matrix *temp = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
        temp->multiply(this->SNPsBLUPGenotypes[ kernels[0]->name ]->genotypes, 'N', kernels[0]->eigenVectors, 'N');
        delete this->SNPsBLUPGenotypes[ kernels[0]->name ]->genotypes;
        this->SNPsBLUPGenotypes[ kernels[0]->name ]->genotypes = temp;
      }
    }
    else if( kernels.size() == 1  && this->nPhenotypes > 1 )
    {
      this->GRMEigenVectors[kernels[0]->name] = new Matrix(kernels[0]->eigenVectors);
      this->usingDiagonalKernels = true;
      
      transformUsingEigenVectors(kernels[0]->name, 'T', &(this->y), this->nPhenotypes);
      
      if(this->SNPsBLUPGenotypes.count(kernels[0]->name) != 0)
      {
        misc.error("Error: An internal error was happened. Unexpected genotype data to be corrected when preparing REML for performing the analysis with diagonal Kernels.", 0);
      }
    }
    else
    {
      misc.error("Error: Sorry, this type of analysis is not implemented, yet.", 0);
    }
    
    this->singlePrecisionInversion = false;
  }
  
  ////////////////////////////////////////////////////
  // Init covariate matrix
  for(int i = 0; i<this->nPhenotypes; i++)
  {
    covariates[i]->filterIndividuals(commonIndividualsInGRMOrder[i]);

    if( kernels.size() == 1 && this->GRMEigenVectors.count(kernels[0]->name) != 0 && this->usingDiagonalKernels == true ) //Correct in case we are using diagonal kernels.
    {
      transformUsingEigenVectors(kernels[0]->name, 'T', &(covariates[i]->covariates), 1);
    }
    
    if( i == 0 )
    {
      this->X = new Matrix(covariates[i]->covariates);
    }
    else
    {
      Matrix * temp = new Matrix(this->X);
      if(options.joinCovariatesVertically == true) //Untested option
      {
        if(temp->nGlobCols != covariates[i]->covariates->nGlobCols)
        {
          misc.error("Error: The covariate or quantitative covariate files of at least two traits do not have the same number of columns.", 0); //.
        }
        this->X->joinMatricesVertically(temp, covariates[i]->covariates);
      }
      else
      {
        subMatrix sm1 = subMatrix(0, 0, temp->nGlobRows, temp->nGlobCols);
        subMatrix sm2 = subMatrix(temp->nGlobRows, temp->nGlobCols, covariates[i]->covariates->nGlobRows, covariates[i]->covariates->nGlobCols);
        this->X->joinMatrices(temp, sm1, covariates[i]->covariates, sm2, 0.);
      }
      delete temp;
    }
    delete covariates[i];
  }
  covariates.clear();
  

  ////////////////////////////////////////////////////
  // Create the covariance matrix
  communicator->broadcast(&nTotalIndividuals, 1);
  
  this->dimension = nTotalIndividuals;
  
  this->V = new CovarianceMatrix(this->dimension, allKernelsDiagonal, this->nPhenotypes);
  
  ////////////////////////////////////////////////////
  // Insert the covariance matrices
  std::vector<std::string> kernelNames; //For storing kernels names before deleting them.
  this->vIndividuals.clear();
  for(int i = 0; i<kernels.size(); i++)
  {
    for(int j = 0; j<this->nPhenotypes; j++)
    {
      kernels[i]->normalize();
      Kernel * filteredKernel = new Kernel(kernels[i]); //Kernel for removing trait individuals
      filteredKernel->filterIndividuals(commonIndividualsInGRMOrder[j], false);
      
      if(i == 0)
      {
        this->vIndividuals.push_back( filteredKernel->individuals );
      }
      
      this->V->insertCovarianceMatrix(kernels[i]->name + "_" + i2s(j+1), filteredKernel);
      
      for(int k = j + 1; k<this->nPhenotypes; k++)
      {
        Kernel * kernelCov = new Kernel(kernels[i]); //Kernel for intersection trait1-trait2
        kernelCov->filterIndividualsAsymmetric(commonIndividualsInGRMOrder[j], commonIndividualsInGRMOrder[k], false);
        
        this->V->insertCovarianceMatrix(kernels[i]->name + "_" + i2s(j+1) + "_" + i2s(k+1), kernelCov);
      }
    }
    
    kernelNames.push_back(kernels[i]->name);
    delete kernels[i];
  }
  
  for(int i = 0; i<this->nPhenotypes; i++)
  {
    Matrix * temp;
    if( allKernelsDiagonal == false )
    {
      temp = new Matrix(cyclicDistribution, nIndividualsTraits[i], nIndividualsTraits[i], communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
    }
    else
    {
      temp = new Matrix(diagonalDistribution, nIndividualsTraits[i], nIndividualsTraits[i], communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
    }
    temp->fillDiagonal(1.);
    
    this->V->insertCovarianceMatrix("E_" + i2s(i+1), temp);
  }
  
  std::map<std::pair<int, int>, bool> computeEnvironmentalCovariances;
  if(options.environmentalCovariance)
  {
    for(int i = 0; i<this->nPhenotypes; i++)
    {
      for(int j = i + 1; j<this->nPhenotypes; j++)
      {
        int commonIndividualsBetweenTraits = ( intersectionStringVectors(2, &(commonIndividualsInGRMOrder[i]), &(commonIndividualsInGRMOrder[j])) ).size();
        int tempBiggestSample = (commonIndividualsInGRMOrder[i].size()>commonIndividualsInGRMOrder[j].size()?commonIndividualsInGRMOrder[i].size():commonIndividualsInGRMOrder[j].size());
        double factor = double(commonIndividualsBetweenTraits)/double(tempBiggestSample);

        if( misc.gt(factor < 0.1) )
        {
          misc.message << "Less than 10% of individuals were measured for both traits. The residual covariance component for traits " << i + 1 << " and " << j + 1 << " is discarded.";
          computeEnvironmentalCovariances[std::pair<int, int>(i, j)] = false;
          continue;
        }
        
        Matrix * covar;
        if( allKernelsDiagonal == false )
        {
          Matrix * temp = new Matrix(cyclicDistribution, nIndividualsBeforeFilterKernel0, nIndividualsBeforeFilterKernel0, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
          temp->fillDiagonal(1.);
          covar = new Matrix(cyclicDistribution);
          temp->filterRowsAndColumns(covar, idxsKeptKernel0[i], idxsKeptKernel0[j]);
          delete temp;
          
          covar->symmetric = false;
          covar->uplo = 'B';
        }
        else
        {
          if( nIndividualsTraits[i] != nIndividualsTraits[j] )
          {
            misc.error("Error: An internal error was happened. traits with different numbers of individuals using diagonal Kernels is not allowed.", 0);
          }
          covar = new Matrix(diagonalDistribution, nIndividualsTraits[i], nIndividualsTraits[j], communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
          covar->fillDiagonal(1.);
        }
  
        this->V->insertCovarianceMatrix("E_" + i2s(i+1) + "_" + i2s(j+1), covar);
        
        computeEnvironmentalCovariances[std::pair<int, int>(i, j)] = true;
      }
    }
  }
  else
  {
    for(int i = 0; i<this->nPhenotypes; i++)
    {
      for(int j = i + 1; j<this->nPhenotypes; j++)
      {
        computeEnvironmentalCovariances[std::pair<int, int>(i, j)] = false;
      }
    }
  }
  
  ////////////////////////////////////////////////////
  // Insert variance groups
  for(int i = 0; i<this->nPhenotypes; i++)
  {
    this->V->insertVarianceGroup("Phenotype_" + i2s(i+1), phenotypeVariances[i]);
    for(int j = i + 1; j<this->nPhenotypes; j++)
    {
      this->V->insertVarianceGroup("Phenotype_" + i2s(i+1) + "_" + i2s(j+1), 0.5*sqrt(phenotypeVariances[i]*phenotypeVariances[j]));
    }
  }
  
  ////////////////////////////////////////////////////
  // Insert variances
  
  for(int i = 0; i<kernels.size(); i++)
  {
    for(int j = 0; j<this->nPhenotypes; j++)
    {
      this->V->insertVariance("Var(" + kernelNames[i] + "_p" + i2s(j+1) + ")", "Phenotype_" + i2s(j+1), variance, genetic, phenotypeVariances[j]*heritabilities[j]*weights[i]);
      
      for(int k = j + 1; k<this->nPhenotypes; k++)
      {
        this->V->insertVariance("Covar(" + kernelNames[i] + "_p" + i2s(j+1) + "-" + i2s(k+1) + ")", "Phenotype_" + i2s(j+1) + "_" + i2s(k+1), covariance, genetic, 0.5*sqrt( phenotypeVariances[j]*heritabilities[j]*weights[i]*phenotypeVariances[k]*heritabilities[k]*weights[i] ) );
      }
    }
  }

  for(int i = 0; i<this->nPhenotypes; i++)
  {
    this->V->insertVariance("Var(E_p" + i2s(i+1) + ")", "Phenotype_" + i2s(i+1), variance, environment, phenotypeVariances[i]*(1.-heritabilities[i]));
    
    for(int j = i + 1; j<this->nPhenotypes; j++)
    {
      if( computeEnvironmentalCovariances[std::pair<int, int>(i, j)] == true )
      {
        this->V->insertVariance("Covar(E_p" + i2s(i+1) + "-" + i2s(j+1) + ")", "Phenotype_" + i2s(i+1) + "_" + i2s(j+1), covariance, environment, 0.5*sqrt(phenotypeVariances[i]*(1.-heritabilities[i])*phenotypeVariances[j]*(1.-heritabilities[j])) );
      }
    }
  }
  
  ////////////////////////////////////////////////////
  // Insert Elements
  
  for(int i = 0; i<kernels.size(); i++)
  {
    int shiftRows = 0;
    int shiftColumns = 0;
    for(int j = 0; j<this->nPhenotypes; j++)
    {
      this->V->insertElement(kernelNames[i], kernelNames[i] + "_" + i2s(j + 1), grmType, kernelNames[i] + "_" + i2s(j + 1), 1., std::pair<int, int>(j, j), subMatrix(shiftRows, shiftColumns, nIndividualsTraits[j], nIndividualsTraits[j]), subMatrix(0, 0, nIndividualsTraits[j], nIndividualsTraits[j]));
      this->V->appendVarianceToElement(kernelNames[i] + "_" + i2s(j + 1), "Var(" + kernelNames[i] + "_p" + i2s(j + 1) + ")", variance);
      
      for(int k = j + 1; k<this->nPhenotypes; k++)
      {
        shiftColumns += nIndividualsTraits[ k - 1 ];
        this->V->insertElement(kernelNames[i], kernelNames[i] + "_" + i2s(j + 1) + "_" + i2s(k + 1), grmType, kernelNames[i] + "_" + i2s(j + 1) + "_" + i2s(k + 1), 1., std::pair<int, int>(j, k), subMatrix(shiftRows, shiftColumns, nIndividualsTraits[j], nIndividualsTraits[k]), subMatrix(0, 0, nIndividualsTraits[j], nIndividualsTraits[k]));
        this->V->appendVarianceToElement(kernelNames[i] + "_" + i2s(j + 1) + "_" + i2s(k + 1), "Covar(" + kernelNames[i] + "_p" + i2s(j + 1) + "-" + i2s(k + 1) + ")", covariance);
      }
      
      shiftRows += nIndividualsTraits[j];
      shiftColumns = shiftRows;
    }
  }

  int shiftRows = 0;
  int shiftColumns = 0;
  for(int i = 0; i<this->nPhenotypes; i++)
  {
    this->V->insertElement("E", "E_" + i2s(i + 1), environtmentalType, "E_" + i2s(i + 1), 1., std::pair<int, int>(i, i), subMatrix(shiftRows, shiftColumns, nIndividualsTraits[i], nIndividualsTraits[i]), subMatrix(0, 0, nIndividualsTraits[i], nIndividualsTraits[i]));
    this->V->appendVarianceToElement("E_" + i2s(i + 1), "Var(E_p" + i2s(i + 1) + ")", variance);
    
    for(int j = i + 1; j<this->nPhenotypes; j++)
    {
      shiftColumns += nIndividualsTraits[ j - 1 ];
      if( computeEnvironmentalCovariances[std::pair<int, int>(i, j)] == true )
      {
        this->V->insertElement("E", "E_" + i2s(i + 1) + "_" + i2s(j + 1), environtmentalType, "E_" + i2s(i + 1) + "_" + i2s(j + 1), 1., std::pair<int, int>(i, j), subMatrix(shiftRows, shiftColumns, nIndividualsTraits[i], nIndividualsTraits[j]), subMatrix(0, 0, nIndividualsTraits[i], nIndividualsTraits[j]));
        this->V->appendVarianceToElement("E_" + i2s(i + 1) + "_" + i2s(j + 1), "Covar(E_p" + i2s(i + 1) + "-" + i2s(j + 1) + ")", covariance);
      }
    }
    
    shiftRows += nIndividualsTraits[i];
    shiftColumns = shiftRows;
  }
  
  kernels.clear();
  
  return true;
}

bool REML::prepare(Matrix* yparam, Matrix* Xparam, std::vector<Matrix*> & kernels, double heritability, std::vector<double> weights)
{
  //Perform some checks and parameter initialization
  
  if(kernels.size() < 1)
  {
    misc.error("Error: An internal error was happened. No GRMs specified for REML analysis.", 0);
  }
  
  if(weights.size() == 0)
  {
    double equallyDistributedWeight = 1./double(kernels.size());
    for(int i = 0; i<kernels.size(); i++)
    {
      weights.push_back(equallyDistributedWeight);
    }
  }
  else
  {
    if( kernels.size() != weights.size() )
    {
      misc.error("Error: An internal error was happened. Kernel weights not properly defined.", 0);
    }
  }
  
  if(heritability < 0)
  {
    heritability = 0.5;
  }
  
  //Set reml type
  
  this->type = rawREMLType;
  
  //Set the number of phenotypes
  
  this->nPhenotypes = 1;
  
  //Check dimensions
  
  for(int i = 0; i<kernels.size(); i++)
  {
    if( kernels[i]->nGlobRows != yparam->nGlobRows || kernels[i]->nGlobCols != yparam->nGlobRows )
    {
      misc.error("Error: An internal error was happened. Different number of individuals in kernels and phenotypes when computing raw REML.", 0);
    }
  }
  if( yparam->nGlobRows != Xparam->nGlobRows )
  {
    misc.error("Error: An internal error was happened. Different number of individuals in covariates and phenotypes when computing raw REML.", 0);
  }
  
  int nTotalIndividuals = yparam->nGlobRows;
  
  ////////////////////////////////////////////////////
  // Set phenotype and covariate matrices
  
  this->y = yparam;
  this->X = Xparam;
  
  ////////////////////////////////////////////////////
  // Diagonal kernels?
  this->usingDiagonalKernels = false;
  bool allKernelsDiagonal = true;
  bool someKernelsDiagonal = false;
  for(int i = 0; i<kernels.size(); i++)
  {
    if(kernels[i]->distribution == diagonalDistribution)
    {
      someKernelsDiagonal = true;
    }
    else
    {
      allKernelsDiagonal = false;
    }
  }
  
  if( allKernelsDiagonal == false && someKernelsDiagonal == true )
  {
    misc.error("Error: An internal error was happened. Diagonal kernels mixed with nondiagonal ones have been found.", 0);
  }
  
  if( allKernelsDiagonal == true )
  {
    this->singlePrecisionInversion = false;
  }
  
  ////////////////////////////////////////////////////
  // Create the covariance matrix
  communicator->broadcast(&nTotalIndividuals, 1);
  
  this->dimension = nTotalIndividuals;
  
  this->V = new CovarianceMatrix(this->dimension, allKernelsDiagonal, this->nPhenotypes);
  
  ////////////////////////////////////////////////////
  // Insert the covariance matrices
  std::vector<std::string> kernelNames; //For storing kernels names before deleting them.
  this->vIndividuals.clear();
  for(int i = 0; i<kernels.size(); i++)
  {
    this->V->insertCovarianceMatrix("K_" + i2s(i + 1), kernels[i]);
    kernelNames.push_back("K_" + i2s(i + 1));
  }
  this->vIndividuals.push_back( std::vector<Individual>() );
  
  Matrix * temp;
  if( allKernelsDiagonal == false )
  {
    temp = new Matrix(cyclicDistribution, nTotalIndividuals, nTotalIndividuals, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  }
  else
  {
    temp = new Matrix(diagonalDistribution, nTotalIndividuals, nTotalIndividuals, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  }
  temp->fillDiagonal(1.);
  
  this->V->insertCovarianceMatrix("E", temp);

  
 ////////////////////////////////////////////////////
  // Insert variance groups
  double phenotypeVariance = computeVariance(this->y);
  this->V->insertVarianceGroup("Phenotype_1", phenotypeVariance);
  
  ////////////////////////////////////////////////////
  // Insert variances
  
  for(int i = 0; i<kernels.size(); i++)
  {
    this->V->insertVariance("Var(" + kernelNames[i] + ")", "Phenotype_1", variance, genetic, phenotypeVariance*heritability*weights[i]);
  }
  this->V->insertVariance("Var(E)", "Phenotype_1", variance, environment, phenotypeVariance*(1.-heritability));
  
  ////////////////////////////////////////////////////
  // Insert Elements
  
  for(int i = 0; i<kernels.size(); i++)
  {
    this->V->insertElement(kernelNames[i], kernelNames[i], grmType, kernelNames[i], 1., std::pair<int, int>(0, 0), subMatrix(0, 0, nTotalIndividuals, nTotalIndividuals), subMatrix(0, 0, nTotalIndividuals, nTotalIndividuals));
    this->V->appendVarianceToElement(kernelNames[i], "Var(" + kernelNames[i] + ")", variance);
  }

  int shiftRows = 0;
  int shiftColumns = 0;

  this->V->insertElement("E", "E", environtmentalType, "E", 1., std::pair<int, int>(0, 0), subMatrix(0, 0, nTotalIndividuals, nTotalIndividuals), subMatrix(0, 0, nTotalIndividuals, nTotalIndividuals));
  this->V->appendVarianceToElement("E", "Var(E)", variance);
 
  kernels.clear();
  
  return true;
}

void REML::createCovarMatrix(Matrix * srcCovarMatrix, Matrix * resultCovarMatrix)
{
  int dimension1 = srcCovarMatrix->nGlobRows;
  int dimension2 = srcCovarMatrix->nGlobCols;
  Matrix * srcCovarMatrixTransp = new Matrix(cyclicDistribution);
  srcCovarMatrixTransp->transpose(srcCovarMatrix);
  resultCovarMatrix->joinMatrices(srcCovarMatrix, subMatrix(0, dimension1, dimension1, dimension2), srcCovarMatrixTransp, subMatrix(dimension1, 0, dimension2, dimension1));
  delete srcCovarMatrixTransp;
}

void REML::transformUsingEigenVectors(std::string name, char t, Matrix ** m, int numberBlocks)
{
  if( this->usingDiagonalKernels == false || this->GRMEigenVectors.count(name) == 0 )
  {
    misc.error("Error: An internal error was happened. The matrix cannot be transformed using eigenvectors. This is not a diagonal analysis or eigenvectors not present.", 0);
  }
  
  Matrix * eigenVectors = this->GRMEigenVectors[name];
  
  if( numberBlocks == 1)
  {
    Matrix * temp = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
    temp->multiply(eigenVectors, t, *m, 'N');
    delete *m;
    *m = temp;
  }
  else
  {
    if( numberBlocks*eigenVectors->nGlobCols != (*m)->nGlobRows )
    {
      misc.error("Error: An internal error was happened. The matrix cannot be transformed using eigenvectors. Dimensions do not match.", 0);
    }
    
    Matrix * temp = new Matrix((*m)->distribution, (*m)->nGlobRows, (*m)->nGlobCols);
    temp->fillWithConstant(0.);
    
    int rowShift = 0;
    for(int ib = 0; ib < numberBlocks; ib++)
    {
      Matrix * mBlock = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, eigenVectors->nGlobRows, (*m)->nGlobCols);
      mBlock->fillWithConstant(0.);
      subMatrix smDest(mBlock);
      subMatrix smSrc(rowShift, 0, eigenVectors->nGlobRows, (*m)->nGlobCols);
      mBlock->add(*m, 1., 1., smDest, smSrc);
      
      Matrix * transformedmBlock = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
      transformedmBlock->multiply(eigenVectors, t, mBlock, 'N');
      
      temp->add(transformedmBlock, 1., 1., smSrc, smDest);
      
      rowShift += eigenVectors->nGlobRows;
    }
    
    delete *m;
    *m = temp;
  }
}

void REML::addCovariatesNames(Covariate * cm, std::string prefix, int shift)
{
  if( communicator->mpiRoot == true )
  {
    CovariateNames tempCovariateNames;
    int idx = shift;
    for( int i = 0; i < cm->meanNames.size(); i++ )
    {
      std::string tempName = prefix + cm->meanNames[i];
      std::pair<std::string, int> tempPair(tempName, idx);
      tempCovariateNames.meanNames.push_back(tempPair);
      idx++;
    }
    for( int i = 0; i < cm->discreteCovarNames.size(); i++ )
    {
      std::string tempName = prefix + cm->discreteCovarNames[i];
      std::pair<std::string, int> tempPair(tempName, idx);
      tempCovariateNames.discreteCovarNames.push_back(tempPair);
      idx++;
    }
    for( int i = 0; i < cm->quantitativeCovarNames.size(); i++ )
    {
      std::string tempName = prefix + cm->quantitativeCovarNames[i];
      std::pair<std::string, int> tempPair(tempName, idx);
      tempCovariateNames.quantitativeCovarNames.push_back(tempPair);
      idx++;
    }
    this->covariateNames.push_back(tempCovariateNames);

    //Check there are not repeated indices.    
    std::set<int> checkIndices;
    for(int idxPheno = 0; idxPheno<this->covariateNames.size(); idxPheno++)
    {
      for( int i = 0; i < this->covariateNames[idxPheno].meanNames.size(); i++ )
      {
        int testIdx = this->covariateNames[idxPheno].meanNames[i].second;
        if( checkIndices.find(testIdx) != checkIndices.end() )
        {
          misc.error("Error: An internal error was happened. When creating covariates names indexes. Repeated index.", 0);
        }
        checkIndices.insert(testIdx);
      }
      for( int i = 0; i < this->covariateNames[idxPheno].discreteCovarNames.size(); i++ )
      {
        int testIdx = this->covariateNames[idxPheno].discreteCovarNames[i].second;
        if( checkIndices.find(testIdx) != checkIndices.end() )
        {
          misc.error("Error: An internal error was happened. When creating covariates names indexes. Repeated index.", 0);
        }
        checkIndices.insert(testIdx);
      }
      for( int i = 0; i < this->covariateNames[idxPheno].quantitativeCovarNames.size(); i++ )
      {
        int testIdx = this->covariateNames[idxPheno].quantitativeCovarNames[i].second;
        if( checkIndices.find(testIdx) != checkIndices.end() )
        {
          misc.error("Error: An internal error was happened. When creating covariates names indexes. Repeated index.", 0);
        }
        checkIndices.insert(testIdx);
      }
    }//End for checking
  } //End mpi root
}

std::string REML::addSuffix(std::string suffix)
{
  if(this->nPhenotypes > 1)
  {
    return suffix;
  }
  return "";
}

void REML::fixCovarianceMatrixVariances()
{
  if(this->type == bivariateREMLType)
  {
    this->V->clearElementVariances("GRM_1_2");
    this->V->deleteVariance("Covar(GRM_p1-2)");
    this->V->changeElementConstantFactor("GRM_1_2", options.fixedCorrelation);
    this->V->appendVarianceToElement("GRM_1_2", "Var(GRM_p1)", standardDeviation);
    this->V->appendVarianceToElement("GRM_1_2", "Var(GRM_p2)", standardDeviation);
  }
  else
  {
    misc.error("Error: Sorry, for this kind of REML, the variance fixation is still not implemented.", 0);
  }
}

double REML::computeREML()
{
  bool iterate = true;
  int nIterations = 0;
  double logLikelihoodDifference;
  
  this->V->reinitializeVariances();
  
  this->logLikelihood = -1.e50;
  
  deleteIntermediateMatrices();
  this->ViX = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  this->XtViX_i = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  this->P = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  this->Py = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  this->subVPy = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->dimension, this->V->variances.size());
  this->AI = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  this->yPsubVPy_trPsubV = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->V->variances.size(), 1);
  
  this->delta = std::vector<double>(this->V->variances.size(), 0.);
  
  misc.message << "Starting REML iterations..." << std::endl;
  misc.message << "Initial variance values:       " << this->V->getStringVarianceValues() << std::endl;
  misc.message << std::setw(9) << "Step  M " << std::setw(3+4) << "Time" << std::setw(options.logFieldWidth+1) << "Prev. LgL " << this->V->getStringVarianceNames() << std::endl;
  
  this->success = true;
  misc.setGetElapsedTime("REMLIteration");
  misc.setGetElapsedTime("REMLAnalysis");
  
  //Start REML iterations
  while(iterate)
  {
    std::stringstream stepResults;
    this->stepModifications = "";
    //REML step
    if(nIterations == 0 || options.REMLMethod == 1)
    {
      stepResults << std::setw(5) << "EM";
      emREMLStep();
    }
    else
    {
      stepResults << std::setw(5) << "AI";
      aiREMLStep();
    }
    if(this->success == false)
    {
      iterate = false;
      break;
    }
    
    //delta loglikelihood
    this->logLikelihoodDifference = computeLogLikelihood();
    
    //Constrain variances
    int nConstrained = this->V->constrainVariancesM1();
    bool variancesConstrained = false;
    std::stringstream constrainedMessage;
    if(nConstrained != 0)
    {
      constrainedMessage << " (" << nConstrained << " of " << this->V->variances.size() << " variances constrained)";
      variancesConstrained = true;
    }
    double test = double(nConstrained)/double(this->V->variances.size());
    if(test > 0.5 && nIterations == 0)
    {
      misc.error("Error: More than half of the variance components are constrained in the first step. REML stopped.", 0);
    }
    else if( test > 0.5 && nIterations != 0 )
    {
      if(options.remlGCTAMode == false)
      {
        this->V->constrainVariancesM2(this->oldVariances, this->delta);
        nConstrained = 0;
        constrainedMessage.str( std::string() );
        constrainedMessage.clear();
        constrainedMessage << " (variances constrained by a scaling factor)";
        variancesConstrained = true;
        
        if( options.allowFixingVariancesToZero == true )
        {
          int nFixedToZero = this->V->fixVariancesToZero();
          if( nFixedToZero != 0 )
          {
            constrainedMessage.str( std::string() );
            constrainedMessage.clear();
            constrainedMessage << " (variances constrained by a scaling factor, " << nFixedToZero << " variances fixed to 0.)";
          }
        }
      }
      else
      {
        misc.error("Error: More than half of the variance components are constrained in the first step. REML stopped.", 0);
      }
    }
    
    //Write REML step results
    if(this->singlePrecisionInversion == true)
    {
      this->stepModifications += "s";
    }
    stepResults << std::setw(3) << this->stepModifications;
    stepResults << std::setprecision(1) << std::setw(8) << misc.setGetElapsedTime("REMLIteration");
    stepResults << std::setprecision(options.logOutputPrecision + 2) << std::setw(options.logFieldWidth) << this->logLikelihood;
    misc.message << stepResults.str() << this->V->getStringVarianceValues() << constrainedMessage.str() << std::endl;
    
    //Other iteration?
    nIterations++;
    if(nIterations>=options.maxREMLIterations)
    {
      this->success = false;
    }
    
    //REML converged?
    //if(varianceDifferenceNorm() < options.varianceConvergenceThreshold && this->logLikelihoodDifference < 1e-4 && this->logLikelihoodDifference > -1e-2)
    if( misc.gt( allVariancesRelativeDifferencesLowerThan(options.varianceConvergenceThreshold) == true  && this->logLikelihoodDifference < 1e-4 && this->logLikelihoodDifference > -1e-2 && variancesConstrained == false) )
    {
      if(this->singlePrecisionInversion == false)
      {
        iterate = false;
      }
      else
      {
        this->singlePrecisionInversion = false;
      }
    }
    
    //Switch to double precision?
    if(this->singlePrecisionInversion == true && allVariancesRelativeDifferencesLowerThan(options.varianceConvergenceThreshold/10.) == true)
    {
      this->singlePrecisionInversion = false;
    }
    
    //If something failed, stop iterations.
    if(this->success == false)
    {
      iterate = false;
    }
  }
  
  //REML converged?
  if(this->success == true && this->writeResults == true && this->type != rawREMLType)
  {
    computeSummary();
    misc.message << "REML analysis has been finished with success after " << misc.setGetElapsedTime("REMLAnalysis", true) << "!" << std::endl;
    
    if(options.computeIndividualsBLUP)
    {
      computeIndividualsBLUP();
    }
    if(options.computeSNPsBLUP)
    {
      computeSNPsBLUP();
    }
    if(options.computeBLUE)
    {
      computeBLUE();
    }
  }
  
  if(this->type != rawREMLType)
  {
    deleteIntermediateMatrices();
  }
  
  if(this->success == false)
  {
    misc.message << "Sorry, REML failed to converge..." << std::endl;
  }
  
  return this->logLikelihood;
}

void REML::computePMatrix()
{
  if( this->V->mStatus != inverseCovarianceMatrix && this->V->mStatus != inverseCovarianceMatrixByBlocks )
  {
    misc.error("Error: An internal error was happened. P matrix cannot be computed if covariance matrix is not inverted.", 0);
  }
  
  if( this->V->mStatus == inverseCovarianceMatrix )
  {
    this->ViX->multiply(this->V->m, 'N', this->X, 'N');
  }
  else
  {
    BlockMatrix ViXtemp;
    ViXtemp.multiply(this->V->mInBlocks, this->X);
    delete this->ViX;
    this->ViX = ViXtemp.block2distributed();
  }
  
  this->XtViX_i->multiply(this->X, 'T', this->ViX, 'N');
  this->XtViX_i->symmetric = true;
  this->XtViX_i->uplo = 'B';
  
  bool inverted = this->XtViX_i->symmetricInvert(&this->logDetXtViX);
  if(inverted == false)
  {
    this->XtViX_i->multiply(this->X, 'T', this->ViX, 'N');
    this->XtViX_i->symmetric = true;
    this->XtViX_i->uplo = 'B';
    inverted = this->XtViX_i->invert(&this->logDetXtViX);
    if(inverted == false)
    {
      misc.message << "Error: XtViX matrix can not be inverted. REML iterations can not continue.";
      this->success = false;
      return;
    }
  }
  
  Matrix * temp = new Matrix(cyclicDistribution);
  temp->multiply(this->ViX, 'N', this->XtViX_i, 'N');
  this->P->multiply(temp, 'N', this->ViX, 'T');
  this->P->symmetric = true;
  this->P->uplo = 'B';
  if( this->V->mStatus == inverseCovarianceMatrix )
  {
    this->P->add(this->V->m, -1., 1.); //These two steps maybe can be optimized to a single one addapting the multiply function.
  }
  else
  {
    Matrix * Vitemp = this->V->mInBlocks.block2distributed();
    this->P->add(Vitemp, -1., 1.); //These two steps maybe can be optimized to a single one addapting the multiply function.
    delete Vitemp;
  }
  
  delete temp;
}

void REML::computePyMatrix()
{
  this->Py->multiply(this->P, 'N', this->y, 'N');
}

void REML::computeSubVPyMatrix()
{
  for(int i = 0; i<this->V->variances.size(); i++)
  {
    Matrix * subV = this->V->computeDerivateCovariance(i);
    this->subVPy->multiply(subV, 'N', this->Py, 'N', 1., subMatrix(0, i, this->dimension, 1));
  }
}


void REML::computeAIMatrix()
{
  Matrix * PsubVPy = new Matrix(cyclicDistribution);

  PsubVPy->multiply(P, 'N', this->subVPy, 'N');
  this->AI->multiply(this->subVPy, 'T', PsubVPy, 'N', 0.5);
  
  this->AI->symmetric = true;
  bool inverted = this->AI->invert();
  if(inverted == false)
  {
    this->success = false;
    return;
  }
  
  delete PsubVPy;
}

void REML::computeyPsubVPy_trPsubVVector(double scale)
{
  Matrix * yPsubVPy = new Matrix(cyclicDistribution);
  yPsubVPy->multiply(this->Py, 'T', this->subVPy, 'N');
  
  double * vectorTraces = new double [this->V->variances.size()];
  Matrix * temp = new Matrix(cyclicDistribution);
  for(int i=0; i<this->V->variances.size(); i++)
  {
    vectorTraces[i] = this->P->traceOfMatrixProduct(this->V->computeDerivateCovariance(i));
  }
  this->yPsubVPy_trPsubV->scatterMatrix(vectorTraces);
  
  this->yPsubVPy_trPsubV->add(yPsubVPy, -scale, scale);

  delete vectorTraces;
  delete temp;  
  delete yPsubVPy;
}

double REML::computeLogLikelihood()
{
  Matrix * temp = new Matrix(cyclicDistribution);
  
  temp->multiply(this->y, 'T', this->Py, 'N');
  double ytPy;
  temp->gatherMatrix(&ytPy);
  communicator->broadcast(&ytPy, 1);
  
  double previuosLogLikelihood = this->logLikelihood;
  double logLikelihoodDifference = this->logLikelihood;
  this->logLikelihood = -0.5*(this->logDetV + this->logDetXtViX + ytPy);
  logLikelihoodDifference = this->logLikelihood - logLikelihoodDifference;
  
  this->logLikelihoodRelativeDifference = logLikelihoodDifference/fabs(previuosLogLikelihood);
  
  delete temp;
  
  return logLikelihoodDifference;
}

void REML::aiREMLStep()
{
  this->V->computeCovariance();
  this->success = this->V->invertCovariance(&this->logDetV, this->singlePrecisionInversion);
  if(this->success == false)
  {
    return;
  }
  
  computePMatrix();
  if(this->success == false)
  {
    return;
  }
  
  computePyMatrix();
  
  computeSubVPyMatrix();
  
  computeAIMatrix();
  if(this->success == false)
  {
    return;
  }
  
  computeyPsubVPy_trPsubVVector(0.5);
  
  Matrix * delta = new Matrix(cyclicDistribution);
  Matrix * newVariances = new Matrix(cyclicDistribution, this->V->variances.size(), 1);
  
  delta->multiply(this->AI, 'N', this->yPsubVPy_trPsubV, 'N');
  
  //delta->showGlobal("delta");
  std::vector<double> vv = this->V->getVectorVariances();
  this->oldVariances = vv;
  
  newVariances->scatterMatrix(&(vv[0]));
  if(this->logLikelihoodRelativeDifference > options.changeAIStepThreshold)
  {
    if( options.allowSwitchFromAItoEM == true ) //Switch to a EM step.
    {
      emPartialREMLStep(newVariances, delta, 0.5);
      this->stepModifications += "e";
    }
    else
    {
      newVariances->add(delta, 1., options.stepWeightingConstant);
      this->stepModifications += "q";
    }
  }
  else
  {
    newVariances->add(delta, 1., 1.);
  }
  
  newVariances->gatherMatrix(&(vv[0]));
  communicator->broadcast(&(vv[0]), vv.size());
  this->V->storeVectorVariances(vv);
  
  delta->gatherMatrix(&this->delta[0]);
  communicator->broadcast(&this->delta[0], vv.size()); 
  
  delete newVariances;
  delete delta;
}

void REML::emREMLStep()
{
  this->V->computeCovariance();
  this->success = this->V->invertCovariance(&this->logDetV, this->singlePrecisionInversion);
  if(this->success == false)
  {
    return;
  }
  
  computePMatrix();
  if(this->success == false)
  {
    return;
  }
  
  computePyMatrix();
  
  computeSubVPyMatrix();
  
  computeyPsubVPy_trPsubVVector();
  
  double * globyPsubVPy_trPsubV = new double [this->V->variances.size()];
  yPsubVPy_trPsubV->gatherMatrix(globyPsubVPy_trPsubV);
  
  communicator->broadcast(globyPsubVPy_trPsubV, this->V->variances.size());
  
  this->oldVariances = this->V->getVectorVariances();
  
  for(int i=0; i<this->V->variances.size(); i++)
  {
    this->V->variances[i].variance = double(this->dimension)*this->V->variances[i].variance + this->V->variances[i].variance*this->V->variances[i].variance*globyPsubVPy_trPsubV[i];
    this->V->variances[i].variance = this->V->variances[i].variance/double(this->dimension);
  }
  
  std::vector<double> vv = this->V->getVectorVariances();
  communicator->broadcast(&(vv[0]), vv.size());
  this->V->storeVectorVariances(vv);
  
  
  delete [] globyPsubVPy_trPsubV;
}

void REML::emPartialREMLStep(Matrix * newVariances, Matrix * delta, double scale)
{
  double * globyPsubVPy_trPsubV = new double [this->V->variances.size()];
  yPsubVPy_trPsubV->gatherMatrix(globyPsubVPy_trPsubV);
  
  communicator->broadcast(globyPsubVPy_trPsubV, this->V->variances.size());
  
  this->oldVariances = this->V->getVectorVariances();
  
  std::vector<double> globalDelta(this->V->variances.size());
  std::vector<double> globalNewVariances(this->V->variances.size());
  for(int i=0; i<this->V->variances.size(); i++)
  {
    globalDelta[i] = this->V->variances[i].variance;
    this->V->variances[i].variance = double(this->dimension)*this->V->variances[i].variance + this->V->variances[i].variance*this->V->variances[i].variance*(globyPsubVPy_trPsubV[i]/scale);
    this->V->variances[i].variance = this->V->variances[i].variance/double(this->dimension);
    globalDelta[i] = this->V->variances[i].variance - globalDelta[i];
    globalNewVariances[i] = this->V->variances[i].variance;
  }
  
  std::vector<double> vv = this->V->getVectorVariances();
  communicator->broadcast(&(vv[0]), vv.size());
  this->V->storeVectorVariances(vv);
  
  if( communicator->mpiRoot && (globalDelta.size() != delta->nGlobRows || this->V->variances.size() != newVariances->nGlobRows) )
  {
    misc.error("Error: An internal error was happened when performing a partial EM step.", 0);
  }
  delta->scatterMatrix(&(globalDelta[0]));
  newVariances->scatterMatrix(&(globalNewVariances[0]));
  
  delete [] globyPsubVPy_trPsubV;
}

double REML::varianceDifferenceNorm()
{
  if(this->V->variances.size() != this->delta.size())
  {
    misc.error("Error: An internal error was happened in varianceDifferenceNorm().", 0);
  }
  
  double result = 0.;
  double norm = 0.;
  for(int i=0; i<this->V->variances.size(); i++)
  {
    double temp = this->V->variances[i].variance - this->oldVariances[i];
    result += temp * temp;
    norm += this->V->variances[i].variance * this->V->variances[i].variance;
  }
  
  result = result/norm;
  communicator->broadcast(&result);
  
  return result;
}

bool REML::allVariancesRelativeDifferencesLowerThan(double threshold)
{
  bool result = true;
  
  for(int i=0; i<this->V->variances.size(); i++)
  {
    double temp = (this->V->variances[i].variance - this->oldVariances[i]);
    temp /= this->oldVariances[i];
    if( fabs(temp) > threshold )
    {
      result = false;
    }
  }
  
  communicator->broadcast(&result);
  
  return result;
}

void REML::computeSummary()
{
  Message message(options.outFile + ".reml");
  std::vector< std::vector<double> > gAI;
  this->AI->symmetrizeTriangularMatrix();
  this->AI->matrixToStandardVector(gAI);
  
  if(options.REMLMethod == 1)
  {
    misc.error("Error: summary for EM method still not implemented.", 0);
  }
  
  if(communicator->mpiRoot)
  {
    message << "Summary results:" << std::endl;
    message << "-----------------------------\n" << std::endl;
    
    //Variances and their errors
    for(int i = 0; i< this->V->variances.size(); i++)
    {
      message << std::setw(0)<< this->V->variances[i].name << std::setw(options.logFieldWidth) << this->V->variances[i].variance << std::setw(options.logFieldWidth) << sqrt(gAI[i][i]) << std::endl;
    }
    
    //Group variances
    std::map<std::string, std::vector<int> > groupedGeneticVariances;
    std::map<std::string, int > groupedEnvironmentVariances;
    for(int i = 0; i< this->V->variances.size(); i++)
    {
      if(this->V->variances[i].type != variance )
      {
        continue;
      }
      VarianceTypeEffect typeEffect = this->V->variances[i].typeEffect;
      if(typeEffect == environment)
      {
        if(groupedEnvironmentVariances.count(this->V->variances[i].group) == 0)
        {
          groupedEnvironmentVariances[ this->V->variances[i].group ] = i;
        }
        else
        {
          misc.error("Error: An internal error was happened when computing the summary statistics. The structure of variances is unexpected. Two genetic variances in the same group.", 0);
        }
      }
      else
      {
        groupedGeneticVariances[ this->V->variances[i].group ].push_back(i);
      }
    }
    
    if( groupedGeneticVariances.size() != groupedEnvironmentVariances.size() )
    {
      misc.error("Error: An internal error was happened when computing the summary statistics. The structure of variances is unexpected. Not all groups have an environment and a genetic variance.", 0);
    }
    
    //Compute heritabilities
    for(std::map<std::string, int >::iterator itEnv = groupedEnvironmentVariances.begin(); itEnv != groupedEnvironmentVariances.end(); ++itEnv)
    {
      if( groupedGeneticVariances.count(itEnv->first) != 1 )
      {
        misc.error("Error: An internal error was happened when computing the summary statistics. The structure of variances is unexpected. The group names differ between environment and genetic variances.", 0);
      }
      if( groupedGeneticVariances[itEnv->first].size() < 1 )
      {
        misc.error("Error: An internal error was happened when computing the summary statistics. The structure of variances is unexpected. The number of genetic variances in the group is 0.", 0);
      }
      message << "\n- " << itEnv->first << ":\n" << std::endl;
      
      int environmentVarianceIdx = itEnv->second;
      double environmentVariance = this->V->variances[environmentVarianceIdx].variance;
      double totalVariance = environmentVariance;
      
      double varianceTotalVariance = gAI[environmentVarianceIdx][environmentVarianceIdx];
      
      for(int idx = 0; idx < groupedGeneticVariances[itEnv->first].size(); idx++)
      {
        int geneticVarianceIdx = groupedGeneticVariances[itEnv->first][idx];
        double geneticVariance = this->V->variances[geneticVarianceIdx].variance;
        
        totalVariance += geneticVariance;
        varianceTotalVariance += gAI[geneticVarianceIdx][environmentVarianceIdx] + gAI[environmentVarianceIdx][geneticVarianceIdx];
        for(int idx2 = 0; idx2 < groupedGeneticVariances[itEnv->first].size(); idx2++)
        {
          int geneticVarianceIdx2 = groupedGeneticVariances[itEnv->first][idx2];
          varianceTotalVariance += gAI[geneticVarianceIdx][geneticVarianceIdx2];
        }
      }
      message << std::setw(0) << "Var(" << itEnv->first << ")" << std::setw(options.logFieldWidth) << totalVariance << std::setw(options.logFieldWidth) << sqrt(varianceTotalVariance) << std::endl;
      
      for(int idx = 0; idx < groupedGeneticVariances[itEnv->first].size(); idx++)
      {
        int geneticVarianceIdx = groupedGeneticVariances[itEnv->first][idx];
        double geneticVariance = this->V->variances[geneticVarianceIdx].variance;
        
        double varianceGeneticVariance = gAI[geneticVarianceIdx][geneticVarianceIdx];
        double cov = gAI[geneticVarianceIdx][environmentVarianceIdx];
        for(int idx2 = 0; idx2 < groupedGeneticVariances[itEnv->first].size(); idx2++)
        {
          int geneticVarianceIdx2 = groupedGeneticVariances[itEnv->first][idx2];
          cov += gAI[geneticVarianceIdx][geneticVarianceIdx2];
        }
        
        double h2 = geneticVariance/totalVariance;
        double varh2;
        varh2 = varianceGeneticVariance/(geneticVariance*geneticVariance);
        varh2 += varianceTotalVariance/(totalVariance*totalVariance);
        varh2 -= 2.*cov/(geneticVariance*totalVariance);
        varh2 *= h2*h2;
        
        message << std::setw(0) << this->V->variances[geneticVarianceIdx].name << "/Var(" << itEnv->first << ")" << std::setw(options.logFieldWidth) << h2 << std::setw(options.logFieldWidth) << sqrt(varh2) << std::endl;
      }
    }
  }
  
  message << std::endl;
  
  if(communicator->mpiRoot)
  {
    int fieldWidth = 15;
    message << "AI Matrix:" << std::endl;
    message << "-----------------------------\n" << std::endl;
    
    message << std::setw(fieldWidth*2);
    for(int i = 0; i< this->V->variances.size(); i++)
    {
      message << this->V->variances[i].name << std::setw(fieldWidth);
    }
    message << std::endl;
    for(int i = 0; i< this->V->variances.size(); i++)
    {
      message << std::setw(fieldWidth) << this->V->variances[i].name;
      for(int j = 0; j< i + 1; j++)
      {
        message << std::setw(fieldWidth) << gAI[i][j];
      }
      message << std::endl;
    }
    message << std::endl;
  }
  
}

void REML::computeBLUE()
{
  misc.message << "Computing BLUEs..." << std::endl;
  
  Matrix * blue = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  Matrix * temp = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  temp->multiply(this->ViX, 'T', this->y, 'N');
  blue->multiply(this->XtViX_i, 'N', temp, 'N');
  
  std::vector<double> globalBLUE;
  blue->matrixToStandardVector(globalBLUE);
  
  std::vector< std::vector<double> > errorBLUE;
  this->XtViX_i->matrixToStandardVector(errorBLUE);
  
  delete blue;
  delete temp;
  
  if(communicator->mpiRoot)
  {
    int nTest = 0;
    for(int idxPheno = 0; idxPheno<this->covariateNames.size(); idxPheno++)
    {
      nTest += (this->covariateNames[idxPheno].meanNames.size() + this->covariateNames[idxPheno].discreteCovarNames.size() + this->covariateNames[idxPheno].quantitativeCovarNames.size());
    }
    if( globalBLUE.size() !=  nTest)
    {
      misc.error("Error: An internal error was happened. The size of covariate names is not of same dimension of BLUE.", 0);
    }

    for(int idxPheno = 0; idxPheno<this->covariateNames.size(); idxPheno++)
    {
      Message message(options.outFile + addSuffix(".pheno" + i2s(idxPheno + 1)) + ".blue.mean");
      for(int i = 0; i < this->covariateNames[idxPheno].meanNames.size(); i++)
      {
        int idx = this->covariateNames[idxPheno].meanNames[i].second;
        message << this->covariateNames[idxPheno].meanNames[i].first << " " << globalBLUE[ idx ] << " (" << sqrt(errorBLUE[ idx ][ idx ]) << ")" << std::endl;
      }
      
      message.redirect(options.outFile + addSuffix(".pheno" + i2s(idxPheno + 1)) + ".blue.discrete");
      for(int i = 0; i < this->covariateNames[idxPheno].discreteCovarNames.size(); i++)
      {
        int idx = this->covariateNames[idxPheno].discreteCovarNames[i].second;
        message << this->covariateNames[idxPheno].discreteCovarNames[i].first << " " << globalBLUE[ idx ] << " (" << sqrt(errorBLUE[ idx ][ idx ]) << ")" << std::endl;
      }
      
      message.redirect(options.outFile + addSuffix(".pheno" + i2s(idxPheno + 1)) + ".blue.quantitative");
      for(int i = 0; i < this->covariateNames[idxPheno].quantitativeCovarNames.size(); i++)
      {
        int idx = this->covariateNames[idxPheno].quantitativeCovarNames[i].second;
        message << this->covariateNames[idxPheno].quantitativeCovarNames[i].first << " " << globalBLUE[ idx ] << " (" << sqrt(errorBLUE[ idx ][ idx ]) << ")" << std::endl;
      }
    } //End for over different phenotypes
  } //End mpiroot if
}

void REML::computeIndividualsBLUP()
{
  for(int idxBLUP = 0; idxBLUP < individualBLUPNames.size(); idxBLUP++)
  {
    std::string name = individualBLUPNames[idxBLUP];
    
    Matrix * temp;
    Matrix * elementMatrix = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->V->dimension, this->V->dimension, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
    
    misc.message << "Computing Individual BLUPs..." << std::endl;
    
    temp = this->V->multiply(name, this->Py);
    if(this->GRMEigenVectors.count(name) != 0)
    {
      transformUsingEigenVectors(name, 'N', &(temp), this->nPhenotypes);
    }
    std::vector<double> blup;
    temp->matrixToStandardVector(blup);
    
    delete temp;
    
    if(communicator->mpiRoot)
    {
      int nTest = 0;
      for(int indPheno = 0; indPheno<this->vIndividuals.size(); indPheno++)
      {
        nTest += this->vIndividuals[indPheno].size();
      }
      
      if(blup.size() == 0)
      {
        misc.error("Error: An internal error was happened. BLUPs array is empty.", 0);
      }
      if(blup.size() != nTest)
      {
        misc.error("Error: An internal error was happened. Dimensions differ when computing individual BLUPs.", 0);
      }
      
      int shift = 0;
      for(int indPheno = 0; indPheno<this->vIndividuals.size(); indPheno++)
      {
        Message message(options.outFile + "." + spacetab2underscore(name) + addSuffix(".pheno" + i2s(indPheno + 1)) + ".blup.indiv");
        message << "FID" << " " << "IID" << " " << "BLUP";

        message << std::endl;
        for(int i = 0; i < vIndividuals[indPheno].size(); i++)
        {
          message << this->vIndividuals[indPheno][i].familyID << " " << this->vIndividuals[indPheno][i].individualID << " " << blup[ shift + i ] << std::endl;
        }
        shift += vIndividuals[indPheno].size();
      } //End for storing different traits
    } //End mpiRoot if
  } //End for over different subCovariances
}

void REML::computeSNPsBLUP()
{
  if( this->nPhenotypes != 1 )
  {
    return;
  }
  
  for(std::map<std::string, Genotype *>::iterator it = this->SNPsBLUPGenotypes.begin(); it != this->SNPsBLUPGenotypes.end(); ++it)
  {
    std::string name = it->first;
    Genotype* genotypes = it->second;
    int phenoIdx = 0;

    misc.message << "Computing SNP BLUPs (" << name << ")..." << std::endl;

    if(genotypes == NULL)
    {
      misc.error("Error: An internal error was happened when computing SNP BLUPs. No genotypes are loaded.", 0);
    }
    
    if( this->mSNPIds.count(name) == 0 )
    {
      misc.error("Error: An internal error was happened. " + name + " is not defined in the mSNPIds on REML.", 0);
    }
    std::multiset<std::string> test1(this->mSNPIds[name].begin(), this->mSNPIds[name].end());
    std::multiset<std::string> test2(genotypes->SNPIds.begin(), genotypes->SNPIds.end());
    if( test1 != test2 )
    {
      misc.message << "WARNING: SNPs in genotype data for computing SNP BLUPs differ from SNPs used for computing GRM." << std::endl;
    }

    std::vector<std::string> tempIndividualIds;
    for(int i = 0; i < this->vIndividuals[phenoIdx].size(); i++)
    {
      tempIndividualIds.push_back(this->vIndividuals[phenoIdx][i].familyID + "@" + this->vIndividuals[phenoIdx][i].individualID);
    }
    
    std::vector< std::string > blupNames;
    Matrix * temp = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);

    double elementVariance = this->V->computeSubCovarianceVariance(name);
    communicator->broadcast(&elementVariance);
    
    genotypes->normalizeGenotypes();
    Genotype * filteredGenotypes =  new Genotype();;
    genotypes->filterSNPsAndIndividuals(genotypes->SNPIds, tempIndividualIds, true, filteredGenotypes);
    if(filteredGenotypes->individualIds != tempIndividualIds)
    {
      misc.error("Error: An error was happened when computing SNP BLUPs. The individuals in the genotype data are not the same or in the same order than in GRMs used for REML analysis.", 0);
    }
    
    Matrix *tempColumnWithOnes = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, filteredGenotypes->nIndividuals, 1);
    Matrix *tempNNonMissings = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
    tempColumnWithOnes->fillWithConstant(1.);
    tempNNonMissings->multiply(filteredGenotypes->missings, 'N', tempColumnWithOnes, 'N');
    std::vector<double> nNonMissings;
    tempNNonMissings->matrixToStandardVector(nNonMissings);
    delete tempColumnWithOnes;
    delete tempNNonMissings;
    
    temp->multiply(filteredGenotypes->genotypes, 'N', this->Py, 'N');
    
    std::vector<double> snpBLUPs;
    temp->matrixToStandardVector(snpBLUPs);

    delete temp;
    
    if(communicator->mpiRoot)
    {
      if( ( snpBLUPs.size() != filteredGenotypes->SNPs.size() ) || ( nNonMissings.size() != snpBLUPs.size() ) )
      {
        misc.error("Error: An internal error was happened when storing SNP BLUPs.", 0);
      }
      
      Message message(options.outFile + "." + spacetab2underscore(name) + ".blup.snps");
      message << "SNP ALLELE BLUP STDEV MEAN NBLUP" << std::endl;
      for(int i = 0; i < snpBLUPs.size(); i++)
      {
        snpBLUPs[i] = snpBLUPs[i]*double(filteredGenotypes->nIndividuals)*elementVariance/(nNonMissings[i]*double(filteredGenotypes->nSNPs));
        
        message << filteredGenotypes->SNPs[i].name << " " << filteredGenotypes->SNPs[i].allele2;
        message << " " << std::setprecision(14) << snpBLUPs[i];
        message << " " << std::setprecision(14) << filteredGenotypes->SNPs[i].standardDev;
        message << " " << std::setprecision(14) << 2.*filteredGenotypes->SNPs[i].p2;
        message << " " << std::setprecision(14) << snpBLUPs[i]/filteredGenotypes->SNPs[i].standardDev;
        message << std::endl;
      }
    }
    
    delete filteredGenotypes;

  } //End different subcovariances iteration.
}
