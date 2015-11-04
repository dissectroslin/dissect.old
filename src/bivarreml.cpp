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

#include "bivarreml.h"
#include "genotype.h"
#include "kernel.h"
#include "covariancematrix.h"
#include "reml.h"
#include "covariate.h"
#include "phenotype.h"
#include "matrix.h"
#include "options.h"
#include "auxiliar.h"
#include "message.h"
#include "misc.h"
#include "global.h"
#include "results.h"

#include <cmath>
#include <set>

BivarREML::BivarREML(bool argWriteResults)
{
  this->writeResults = argWriteResults;
  
  this->reml = NULL;
}

BivarREML::~BivarREML()
{
  if(this->reml != NULL)
  {
    delete this->reml;
    this->reml = NULL;
  }
}

void BivarREML::compute()
{
  if(this->reml != NULL)
  {
    delete this->reml;
  }
  
  Kernel *grm = loadGRMUsingOptions();
  this->reml = new REML(this->writeResults);
  prepare(grm);
  this->reml->computeREML();

  if(options.fixCorrelation == true)
  {
    double baseLogLikelihood = this->reml->logLikelihood;
    double baseNVariances = this->reml->V->variances.size();
    
    this->reml->fixCovarianceMatrixVariances();
    this->reml->writeResults = false;
    this->reml->computeREML();
    compareREMLs(baseLogLikelihood, baseNVariances);
  }
  
  delete this->reml;
  this->reml = NULL;
}

void BivarREML::computeRegional()
{
  if(this->reml != NULL)
  {
    delete this->reml;
  }
  
  //Load genotypes and GRM
  Genotype * genotype = loadGenotypeUsingOptions();
  Kernel * grm;
  if(options.grmFile == "")
  {
    grm = new Kernel(genotype);
  }
  else
  {
    grm = new Kernel(options.grmFile);
    genotype->normalizeGenotypes();
  }
  
  if(grm->individualIds != genotype->individualIds)
  {
    misc.error("Error: The individuals in GRM and genotypes are not the same. Maybe you are not using the proper GRM or genotype file?", 0);
  }

  if( grm->diagonalized == true )
  {
    if(options.forceUseDiagonalizedKernels == true)
    {
      grm->recoverKernelFromEigenDecomposition();
    }
    else
    {
      misc.error("Error: Sorry, this analysis cannot be performed with diagonal GRMs. You can use the option --force-use-diag-grms to force converting diagonal GRMs to their non-diagonalized form before starting the analysis.", 0);
    }
  }
  
  //Create groups
  genotype->groupSNPs(options.regionBy);
  
  //Start iterations over groups
  std::string baseOutFile = options.outFile;
  int nGroups = genotype->groupedSNPs.size();
  communicator->broadcast(&nGroups);
  std::map<std::string, std::set<std::string> >::iterator it = genotype->groupedSNPs.begin();
  for(int i = 0; i < nGroups; i++)
  {
    std::string group = "";
    if(communicator->mpiRoot)
    {
      group = it->first;
      it++;
    }
    communicator->broadcast(group);
    
    options.outFile = baseOutFile + "." + group;
    
    this->reml = new REML();
    
    Genotype *regionalGenotype = new Genotype();
    genotype->genotypeOfSNPsGroup(group, regionalGenotype);
    double proportionRegionalTotal = double(regionalGenotype->nSNPs)/double(genotype->nSNPs);
    Kernel *regionalGRM = new Kernel(regionalGenotype);
    Kernel *globalGRM = new Kernel();
    globalGRM->addKernels(1., grm, -1., regionalGRM);
    delete regionalGenotype;
    
    std::vector<std::string> snpIntersection = orderVectorAsTemplate(regionalGRM->randomVarNames, grm->randomVarNames);
    if( snpIntersection != regionalGRM->randomVarNames )
    {
      misc.error("Error: SNPs in regional GRM are not present in global GRM. This could be due that some SNPs in genotype file are not present in the the GRM (i.e. they are not used for computing the global GRM).", 0);
    }
    
    std::vector<Kernel*> grms;
    grms.push_back(globalGRM);
    grms.push_back(regionalGRM);
    
    std::vector<std::string> names;
    names.push_back("Global GRM");
    names.push_back("Regional GRM");
    
    std::vector<double> weights;
    weights.push_back(1. - proportionRegionalTotal);
    weights.push_back(proportionRegionalTotal);
    
    prepare(grms, names, weights);
    this->reml->computeREML();
    
    delete this->reml;
    this->reml = NULL;
  }
  
  delete grm;
  delete genotype;
}

void BivarREML::prepare(Kernel * grm)
{
  //Perform some checks
  
  if(grm->diagonalized == true)
  {
    if(options.forceUseDiagonalizedKernels == true)
    {
      grm->recoverKernelFromEigenDecomposition();
    }
    else
    {
      misc.error("Error: Sorry, this analysis cannot be performed with diagonal GRMs. You can use the option --force-use-diag-grms to force converting diagonal GRMs to their non-diagonalized form before starting the analysis.", 0);
    }
  }
  
  //Set reml type
  
  this->reml->type = bivariateREMLType;
  
  //Prepare the GRM (covs), covariance (X) and phenotype (y) matrices
  Kernel *grm1 = grm;
  Kernel *grm2 = new Kernel(grm1);
  Kernel *grm12 = new Kernel(grm1);

  this->reml->SNPIds = grm->randomVarNames;
  
  Phenotype * phenotype1 = new Phenotype(cyclicDistribution, options.phenotypesFile, options.phenotypeColumn1);
  Phenotype * phenotype2 = new Phenotype(cyclicDistribution, options.phenotypesFile, options.phenotypeColumn2);
  this->reml->meanNames.clear();
  this->reml->quantitativeCovarNames.clear();
  this->reml->discreteCovarNames.clear();
  Covariate * covariate1;
  Covariate * covariate2;
  if(options.joinCovariatesVertically == true)
  {
    covariate1 = new Covariate(options.covarsFile1, options.qCovarsFile1, phenotype1->individualIds, false);
    covariate2 = new Covariate(options.covarsFile2, options.qCovarsFile2, phenotype2->individualIds, false);
    covariate1->syncronizeDiscreteCovariateCategoriesWith(covariate2);
    covariate1->parseRawCovariates(phenotype1->individualIds, 2, 0);
    covariate2->parseRawCovariates(phenotype2->individualIds, 2, 1);
    this->reml->addCovariatesNames(covariate1);
  }
  else
  {
    covariate1 = new Covariate(options.covarsFile1, options.qCovarsFile1, phenotype1->individualIds);
    covariate2 = new Covariate(options.covarsFile2, options.qCovarsFile2, phenotype2->individualIds);
    this->reml->addCovariatesNames(covariate1, "t1-");
    this->reml->addCovariatesNames(covariate2, "t2-", covariate1->covariates->nGlobCols);
  }
  
  std::vector<std::string> commonIndividuals1, commonIndividuals2, commonIndividualsInGRMOrder1, commonIndividualsInGRMOrder2;
  std::vector<int> idxsKeptGRM1, idxsKeptGRM2;
  int commonIndividualsBetweenTraits;
  int totalIndividuals;
  int nIndividualsBeforeFilter;
  if(communicator->mpiRoot)
  {
    commonIndividuals1 = intersectionStringVectors(3, &grm1->individualIds, &phenotype1->individualIds, &covariate1->individualIds);
    commonIndividualsInGRMOrder1 = orderVectorAsTemplate(grm1->individualIds, commonIndividuals1);
    commonIndividuals2 = intersectionStringVectors(3, &grm2->individualIds, &phenotype2->individualIds, &covariate2->individualIds);
    commonIndividualsInGRMOrder2 = orderVectorAsTemplate(grm2->individualIds, commonIndividuals2);
    
    idxsKeptGRM1 = extractMapValues(commonIndividualsInGRMOrder1, grm1->individualIdsIdx);
    idxsKeptGRM2 = extractMapValues(commonIndividualsInGRMOrder2, grm2->individualIdsIdx);
    nIndividualsBeforeFilter = grm1->individualIds.size(); //The individuals in grm1 and grm2 must be the same at this point.
    
    commonIndividualsBetweenTraits = (intersectionStringVectors(2, &commonIndividuals1, &commonIndividuals2)).size();
    totalIndividuals = (commonIndividuals1.size()>commonIndividuals2.size()?commonIndividuals1.size():commonIndividuals2.size());
  }
  communicator->broadcast(&commonIndividualsBetweenTraits, 1);
  communicator->broadcast(&totalIndividuals, 1);
  communicator->broadcast(&nIndividualsBeforeFilter, 1);
  if(totalIndividuals == 0)
  {
    misc.error("Error: There are not enough individuals to start the analysis.", 0);
  }
  
  bool computeResidualCovariance = options.environmentalCovariance;
  if( (double(commonIndividualsBetweenTraits)/double(totalIndividuals) < 0.1) && computeResidualCovariance == true )
  {
    misc.message << "Less than 10% of individuals were measured for both traits. The residual covariance component is discarded.";
    computeResidualCovariance = false;
  }
  
  misc.message << commonIndividualsInGRMOrder1.size() << " individuals are kept for the analysis for phenotype 1." << std::endl;
  misc.message << commonIndividualsInGRMOrder2.size() << " individuals are kept for the analysis for phenotype 2." << std::endl;
  
  ////////////////////////////////////////////////////
  // Init phenotype matrix
  phenotype1->filterIndividuals(commonIndividualsInGRMOrder1);
  phenotype2->filterIndividuals(commonIndividualsInGRMOrder2);
  double phenotypeVariance1 = phenotype1->computePhenotypeVariance();
  double phenotypeVariance2 = phenotype2->computePhenotypeVariance();
  this->reml->y = new Matrix(cyclicDistribution);
  this->reml->y->joinMatricesVertically(phenotype1->phenotypes, phenotype2->phenotypes);
  delete phenotype1;
  delete phenotype2;
  
  ////////////////////////////////////////////////////
  // Init covariate matrix
  covariate1->filterIndividuals(commonIndividualsInGRMOrder1);
  covariate2->filterIndividuals(commonIndividualsInGRMOrder2);
  this->reml->X = new Matrix(cyclicDistribution);
  if(options.joinCovariatesVertically == true)
  {
    if(covariate1->covariates->nGlobCols != covariate2->covariates->nGlobCols)
    {
      misc.error("Error: The covariate or quantitative covariate files of both traits do not have the same number of columns.", 0);
    }
    this->reml->X->joinMatricesVertically(covariate1->covariates, covariate2->covariates);
  }
  else
  {
    subMatrix sm1 = subMatrix(0, 0, covariate1->covariates->nGlobRows, covariate1->covariates->nGlobCols);
    subMatrix sm2 = subMatrix(covariate1->covariates->nGlobRows, covariate1->covariates->nGlobCols, covariate2->covariates->nGlobRows, covariate2->covariates->nGlobCols);
    this->reml->X->joinMatrices(covariate1->covariates, sm1, covariate2->covariates, sm2, 0.);
  }
  delete covariate1;
  delete covariate2;

  ////////////////////////////////////////////////////
  // Create the covariance matrix
  int nIndividualsTrait1 = commonIndividualsInGRMOrder1.size();
  int nIndividualsTrait2 = commonIndividualsInGRMOrder2.size();
  communicator->broadcast(&nIndividualsTrait1, 1);
  communicator->broadcast(&nIndividualsTrait2, 1);
  this->reml->dimension = nIndividualsTrait1 + nIndividualsTrait2;
  
  this->reml->V = new CovarianceMatrix(this->reml->dimension);
  
  ////////////////////////////////////////////////////
  // Insert the covariance matrices
  grm1->filterIndividuals(commonIndividualsInGRMOrder1, false);
  grm2->filterIndividuals(commonIndividualsInGRMOrder2, false);
  grm12->filterIndividualsAsymmetric(commonIndividualsInGRMOrder1, commonIndividualsInGRMOrder2, false);
  this->reml->individuals.clear();
  this->reml->individuals.reserve( grm1->individuals.size() + grm2->individuals.size() );
  this->reml->individuals.insert( this->reml->individuals.end(), grm1->individuals.begin(), grm1->individuals.end() );
  this->reml->individuals.insert( this->reml->individuals.end(), grm2->individuals.begin(), grm2->individuals.end() );

  
  Matrix * grm12Covar = new Matrix(cyclicDistribution);
  this->reml->createCovarMatrix(grm12->getNormalizedKernel(), grm12Covar);
  delete grm12;
  
  this->reml->V->insertCovarianceMatrix("GRM1", grm1);
  this->reml->V->insertCovarianceMatrix("GRM2", grm2);
  this->reml->V->insertCovarianceMatrix("GRM12", grm12Covar);
  delete grm1;
  delete grm2;
  delete grm12Covar;
  
  Matrix * temp1 = new Matrix(cyclicDistribution, nIndividualsTrait1, nIndividualsTrait1, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  Matrix * temp2 = new Matrix(cyclicDistribution, nIndividualsTrait2, nIndividualsTrait2, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  temp1->fillDiagonal(1.);
  temp2->fillDiagonal(1.);
  this->reml->V->insertCovarianceMatrix("E1", temp1);
  this->reml->V->insertCovarianceMatrix("E2", temp2);
  delete temp1;
  delete temp2;
  
  if(computeResidualCovariance)
  {
    temp1 = new Matrix(cyclicDistribution, nIndividualsBeforeFilter, nIndividualsBeforeFilter, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
    temp1->fillDiagonal(1.);
    
    Matrix * temp2 = new Matrix(cyclicDistribution);
    temp1->filterRowsAndColumns(temp2, idxsKeptGRM1, idxsKeptGRM2);
    delete temp1;
    
    Matrix * tempCovar = new Matrix(cyclicDistribution);
    this->reml->createCovarMatrix(temp2, tempCovar);
    tempCovar->symmetric = true;
    tempCovar->uplo = 'B';
    delete temp2;

    this->reml->V->insertCovarianceMatrix("E12", tempCovar);
    delete tempCovar;
  }
  
  ////////////////////////////////////////////////////
  // Insert variance groups
  this->reml->V->insertVarianceGroup("Phenotype 1", phenotypeVariance1);
  this->reml->V->insertVarianceGroup("Phenotype 2", phenotypeVariance2);
  this->reml->V->insertVarianceGroup("Phenotype 12", 0.5*sqrt(phenotypeVariance1*phenotypeVariance2));
  
  ////////////////////////////////////////////////////
  // Insert variances
  this->reml->V->insertVariance("Var(GRM1)", "Phenotype 1", variance, genetic, phenotypeVariance1*options.initialh2Trait1);
  this->reml->V->insertVariance("Var(GRM2)", "Phenotype 2", variance, genetic, phenotypeVariance2*options.initialh2Trait2);
  this->reml->V->insertVariance("Covar(GRM12)", "Phenotype 12", covariance, genetic, 0.5*sqrt(phenotypeVariance1*options.initialh2Trait1*phenotypeVariance2*options.initialh2Trait2) );
  
  this->reml->V->insertVariance("Var(E1)", "Phenotype 1", variance, environment, phenotypeVariance1*(1.-options.initialh2Trait1));
  this->reml->V->insertVariance("Var(E2)", "Phenotype 2", variance, environment, phenotypeVariance2*(1.-options.initialh2Trait2));
  if(computeResidualCovariance)
  {
    this->reml->V->insertVariance("Covar(E12)", "Phenotype 12", covariance, environment, 0.5*sqrt(phenotypeVariance1*(1.-options.initialh2Trait1)*phenotypeVariance2*(1.-options.initialh2Trait2)));
  }

  ////////////////////////////////////////////////////
  // Insert Elements
  this->reml->V->insertElement("GRM1", grmType, "GRM1", 1., subMatrix(0, 0, nIndividualsTrait1, nIndividualsTrait1), subMatrix(0, 0, nIndividualsTrait1, nIndividualsTrait1));
  this->reml->V->appendVarianceToElement("GRM1", "Var(GRM1)", variance);
  
  this->reml->V->insertElement("GRM2", grmType, "GRM2", 1., subMatrix(nIndividualsTrait1, nIndividualsTrait1, nIndividualsTrait2, nIndividualsTrait2), subMatrix(0, 0, nIndividualsTrait2, nIndividualsTrait2));
  this->reml->V->appendVarianceToElement("GRM2", "Var(GRM2)", variance);
  
  this->reml->V->insertElement("GRM12", grmType, "GRM12", 1., subMatrix(), subMatrix());
  this->reml->V->appendVarianceToElement("GRM12", "Covar(GRM12)", covariance);
  
  this->reml->V->insertElement("E1", environtmentalType, "E1", 1., subMatrix(0, 0, nIndividualsTrait1, nIndividualsTrait1), subMatrix(0, 0, nIndividualsTrait1, nIndividualsTrait1));
  this->reml->V->appendVarianceToElement("E1", "Var(E1)", variance);
  
  this->reml->V->insertElement("E2", environtmentalType, "E2", 1., subMatrix(nIndividualsTrait1, nIndividualsTrait1, nIndividualsTrait2, nIndividualsTrait2), subMatrix(0, 0, nIndividualsTrait2, nIndividualsTrait2));
  this->reml->V->appendVarianceToElement("E2", "Var(E2)", variance);
  
  if(computeResidualCovariance)
  {
    this->reml->V->insertElement("E12", environtmentalType, "E12", 1., subMatrix(), subMatrix());
    this->reml->V->appendVarianceToElement("E12", "Covar(E12)", covariance);
  }
}

void BivarREML::prepare(std::vector<Kernel*> & grms, std::vector<std::string> names, std::vector<double> weights)
{
  //Perform some checks
  
  if(grms.size() < 1)
  {
    misc.error("Error: An internal error was happened. No GRMs for REML analysis.", 0);
  }
  if( grms.size() != names.size() )
  {
    misc.error("Error: An internal error was happened. GRM names not properly defined.", 0);
  }
  if(weights.size() == 0)
  {
    double equallyDistributedWeight = 1./double(grms.size());
    for(int i = 0; i<grms.size(); i++)
    {
      weights.push_back(equallyDistributedWeight);
    }
  }
  else
  {
    if( grms.size() != weights.size() )
    {
      misc.error("Error: An internal error was happened. GRM weights not properly defined.", 0);
    }
  }
  
  for(int i = 0; i<grms.size(); i++)
  {
    if(grms[i]->diagonalized == true)
    {
      if(options.forceUseDiagonalizedKernels == true)
      {
        grms[i]->recoverKernelFromEigenDecomposition();
      }
      else
      {
        misc.error("Error: Sorry, this analysis cannot be performed with diagonal GRMs. You can use the option --force-use-diag-grms to force converting diagonal GRMs to their non-diagonalized form before starting the analysis.", 0);
      }
    }
  }
  
  //Set reml type
  
  this->reml->type = multipleBivarREMLType;
  
  //Prepare SNPs
  
  this->reml->SNPIds.clear();
  
  //Prepare the GRM (covs), covariance (X) and phenotype (y) matrices
  
  Phenotype * phenotype1 = new Phenotype(cyclicDistribution, options.phenotypesFile, options.phenotypeColumn1);
  Phenotype * phenotype2 = new Phenotype(cyclicDistribution, options.phenotypesFile, options.phenotypeColumn2);
  this->reml->meanNames.clear();
  this->reml->quantitativeCovarNames.clear();
  this->reml->discreteCovarNames.clear();
  Covariate * covariate1;
  Covariate * covariate2;
  if(options.joinCovariatesVertically == true)
  {
    covariate1 = new Covariate(options.covarsFile1, options.qCovarsFile1, phenotype1->individualIds, false);
    covariate2 = new Covariate(options.covarsFile2, options.qCovarsFile2, phenotype2->individualIds, false);
    covariate1->syncronizeDiscreteCovariateCategoriesWith(covariate2);
    covariate1->parseRawCovariates(phenotype1->individualIds, 2, 0);
    covariate2->parseRawCovariates(phenotype2->individualIds, 2, 1);
    this->reml->addCovariatesNames(covariate1);
  }
  else
  {
    covariate1 = new Covariate(options.covarsFile1, options.qCovarsFile1, phenotype1->individualIds);
    covariate2 = new Covariate(options.covarsFile2, options.qCovarsFile2, phenotype2->individualIds);
    this->reml->addCovariatesNames(covariate1, "t1-");
    this->reml->addCovariatesNames(covariate2, "t2-", covariate1->covariates->nGlobCols);
  }
  
  std::vector<std::string> commonIndividuals1, commonIndividuals2, commonIndividualsInGRMOrder1, commonIndividualsInGRMOrder2;
  std::vector<int> idxsKeptGRM1, idxsKeptGRM2;
  int commonIndividualsBetweenTraits;
  int totalIndividuals;
  int nIndividualsBeforeFilter;
  if(communicator->mpiRoot)
  {
    commonIndividuals1 = intersectionStringVectors(3, &(grms[0]->individualIds), &phenotype1->individualIds, &covariate1->individualIds);
    commonIndividuals2 = intersectionStringVectors(3, &(grms[0]->individualIds), &phenotype2->individualIds, &covariate2->individualIds);
    for(int i = 1; i<grms.size(); i++)
    {
      commonIndividuals1 = intersectionStringVectors(2, &commonIndividuals1, &(grms[i]->individualIds));
      commonIndividuals2 = intersectionStringVectors(2, &commonIndividuals2, &(grms[i]->individualIds));
    }
    commonIndividualsInGRMOrder1 = orderVectorAsTemplate(grms[0]->individualIds, commonIndividuals1);
    commonIndividualsInGRMOrder2 = orderVectorAsTemplate(grms[0]->individualIds, commonIndividuals2);
    for(int i = 1; i<grms.size(); i++)
    {
      std::vector<std::string> tempTest1 = orderVectorAsTemplate(grms[i]->individualIds, commonIndividuals1);
      std::vector<std::string> tempTest2 = orderVectorAsTemplate(grms[i]->individualIds, commonIndividuals2);
      if( tempTest1 != commonIndividualsInGRMOrder1 || tempTest2 != commonIndividualsInGRMOrder2 )
      {
        misc.error("Error: Not all GRMs have the individuals in the same order. At this moment, DISSECT cannot resort it. Please, resort it or contact us.", 0);
      }
    }

    idxsKeptGRM1 = extractMapValues(commonIndividualsInGRMOrder1, grms[0]->individualIdsIdx);
    idxsKeptGRM2 = extractMapValues(commonIndividualsInGRMOrder2, grms[0]->individualIdsIdx);
    nIndividualsBeforeFilter = grms[0]->individualIds.size();
    
    commonIndividualsBetweenTraits = (intersectionStringVectors(2, &commonIndividuals1, &commonIndividuals2)).size();
    totalIndividuals = (commonIndividuals1.size()>commonIndividuals2.size()?commonIndividuals1.size():commonIndividuals2.size());
  }
  communicator->broadcast(&commonIndividualsBetweenTraits, 1);
  communicator->broadcast(&totalIndividuals, 1);
  communicator->broadcast(&nIndividualsBeforeFilter, 1);
  if(totalIndividuals == 0)
  {
    misc.error("Error: There are not enough individuals to start the analysis.", 0);
  }
  
  bool computeResidualCovariance = options.environmentalCovariance;
  if( (double(commonIndividualsBetweenTraits)/double(totalIndividuals) < 0.1) && computeResidualCovariance == true)
  {
    misc.message << "Less than 10% of individuals were measured for both traits. The residual covariance component is discarded.";
    computeResidualCovariance = false;
  }
  
  misc.message << commonIndividualsInGRMOrder1.size() << " individuals are kept for the analysis for phenotype 1." << std::endl;
  misc.message << commonIndividualsInGRMOrder2.size() << " individuals are kept for the analysis for phenotype 2." << std::endl;
  
  ////////////////////////////////////////////////////
  // Init phenotype matrix
  phenotype1->filterIndividuals(commonIndividualsInGRMOrder1);
  phenotype2->filterIndividuals(commonIndividualsInGRMOrder2);
  double phenotypeVariance1 = phenotype1->computePhenotypeVariance();
  double phenotypeVariance2 = phenotype2->computePhenotypeVariance();
  this->reml->y = new Matrix(cyclicDistribution);
  this->reml->y->joinMatricesVertically(phenotype1->phenotypes, phenotype2->phenotypes);
  delete phenotype1;
  delete phenotype2;
  
  ////////////////////////////////////////////////////
  // Init covariate matrix
  covariate1->filterIndividuals(commonIndividualsInGRMOrder1);
  covariate2->filterIndividuals(commonIndividualsInGRMOrder2);
  this->reml->X = new Matrix(cyclicDistribution);
  if(options.joinCovariatesVertically == true)
  {
    if(covariate1->covariates->nGlobCols != covariate2->covariates->nGlobCols)
    {
      misc.error("Error: The covariate or quantitative covariate files of both traits do not have the same number of columns.", 0); //.
    }
    this->reml->X->joinMatricesVertically(covariate1->covariates, covariate2->covariates);
  }
  else
  {
    subMatrix sm1 = subMatrix(0, 0, covariate1->covariates->nGlobRows, covariate1->covariates->nGlobCols);
    subMatrix sm2 = subMatrix(covariate1->covariates->nGlobRows, covariate1->covariates->nGlobCols, covariate2->covariates->nGlobRows, covariate2->covariates->nGlobCols);
    this->reml->X->joinMatrices(covariate1->covariates, sm1, covariate2->covariates, sm2, 0.);
  }
  delete covariate1;
  delete covariate2;
  
  ////////////////////////////////////////////////////
  // Create the covariance matrix
  int nIndividualsTrait1 = commonIndividualsInGRMOrder1.size();
  int nIndividualsTrait2 = commonIndividualsInGRMOrder2.size();
  communicator->broadcast(&nIndividualsTrait1, 1);
  communicator->broadcast(&nIndividualsTrait2, 1);
  this->reml->dimension = nIndividualsTrait1 + nIndividualsTrait2;
  
  this->reml->V = new CovarianceMatrix(this->reml->dimension);
  
  ////////////////////////////////////////////////////
  // Insert the covariance matrices
  for(int i = 0; i<grms.size(); i++)
  {
    grms[i]->normalize();
    Kernel * grm2 = new Kernel(grms[i]); //GRM for the second trait
    Kernel * grm12 = new Kernel(grms[i]); //GRM for intersection trait1-trait2
    grms[i]->filterIndividuals(commonIndividualsInGRMOrder1, false);
    grm2->filterIndividuals(commonIndividualsInGRMOrder2, false);
    grm12->filterIndividualsAsymmetric(commonIndividualsInGRMOrder1, commonIndividualsInGRMOrder2, false);
    
    if( i == 0 )
    {
      this->reml->individuals.clear();
      this->reml->individuals.reserve( grms[i]->individuals.size() + grm2->individuals.size() );
      this->reml->individuals.insert( this->reml->individuals.end(), grms[i]->individuals.begin(), grms[i]->individuals.end() );
      this->reml->individuals.insert( this->reml->individuals.end(), grm2->individuals.begin(), grm2->individuals.end() );
    }
    
    Matrix * grm12Covar = new Matrix(cyclicDistribution);
    this->reml->createCovarMatrix(grm12->getNormalizedKernel(), grm12Covar);
    delete grm12;
    
    this->reml->V->insertCovarianceMatrix(names[i] + "_1", grms[i]);
    this->reml->V->insertCovarianceMatrix(names[i] + "_2", grm2);
    this->reml->V->insertCovarianceMatrix(names[i] + "_12", grm12Covar);
    delete grms[i];
    delete grm2;
    delete grm12Covar;
  }
 
  Matrix * temp1 = new Matrix(cyclicDistribution, nIndividualsTrait1, nIndividualsTrait1, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  Matrix * temp2 = new Matrix(cyclicDistribution, nIndividualsTrait2, nIndividualsTrait2, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  temp1->fillDiagonal(1.);
  temp2->fillDiagonal(1.);
  this->reml->V->insertCovarianceMatrix("E1", temp1);
  this->reml->V->insertCovarianceMatrix("E2", temp2);
  delete temp1;
  delete temp2;
  
  if(computeResidualCovariance)
  {
    temp1 = new Matrix(cyclicDistribution, nIndividualsBeforeFilter, nIndividualsBeforeFilter, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
    temp1->fillDiagonal(1.);
    
    Matrix * temp2 = new Matrix(cyclicDistribution);
    temp1->filterRowsAndColumns(temp2, idxsKeptGRM1, idxsKeptGRM2);
    delete temp1;
    
    Matrix * tempCovar = new Matrix(cyclicDistribution);
    this->reml->createCovarMatrix(temp2, tempCovar);
    tempCovar->symmetric = true;
    tempCovar->uplo = 'B';
    delete temp2;
    
    this->reml->V->insertCovarianceMatrix("E12", tempCovar);
    delete tempCovar;
  }
  
  ////////////////////////////////////////////////////
  // Insert variance groups
  this->reml->V->insertVarianceGroup("Phenotype 1", phenotypeVariance1);
  this->reml->V->insertVarianceGroup("Phenotype 2", phenotypeVariance2);
  this->reml->V->insertVarianceGroup("Phenotype 12", 0.5*sqrt(phenotypeVariance1*phenotypeVariance2));
  
  ////////////////////////////////////////////////////
  // Insert variances
  
  for(int i = 0; i<grms.size(); i++)
  {
    this->reml->V->insertVariance("Var(" + names[i] + "_1)", "Phenotype 1", variance, genetic, phenotypeVariance1*options.initialh2Trait1*weights[i]);
    this->reml->V->insertVariance("Var(" + names[i] + "_2)", "Phenotype 2", variance, genetic, phenotypeVariance2*options.initialh2Trait2*weights[i]);
    this->reml->V->insertVariance("Covar(" + names[i] + "_12)", "Phenotype 12", covariance, genetic, 0.5*sqrt( phenotypeVariance1*options.initialh2Trait1*weights[i]*phenotypeVariance2*options.initialh2Trait2*weights[i] ) );
  }
  
  this->reml->V->insertVariance("Var(E1)", "Phenotype 1", variance, environment, phenotypeVariance1*(1.-options.initialh2Trait1));
  this->reml->V->insertVariance("Var(E2)", "Phenotype 2", variance, environment, phenotypeVariance2*(1.-options.initialh2Trait2));
  if(computeResidualCovariance)
  {
    this->reml->V->insertVariance("Covar(E12)", "Phenotype 12", covariance, environment, 0.5*sqrt(phenotypeVariance1*(1.-options.initialh2Trait1)*phenotypeVariance2*(1.-options.initialh2Trait2)) );
  }
  
  ////////////////////////////////////////////////////
  // Insert Elements
  
  for(int i = 0; i<grms.size(); i++)
  {
    this->reml->V->insertElement(names[i] + "_1", grmType, names[i] + "_1", 1., subMatrix(0, 0, nIndividualsTrait1, nIndividualsTrait1), subMatrix(0, 0, nIndividualsTrait1, nIndividualsTrait1));
    this->reml->V->appendVarianceToElement(names[i] + "_1", "Var(" + names[i] + "_1)", variance);
    
    this->reml->V->insertElement(names[i] + "_2", grmType, names[i] + "_2", 1., subMatrix(nIndividualsTrait1, nIndividualsTrait1, nIndividualsTrait2, nIndividualsTrait2), subMatrix(0, 0, nIndividualsTrait2, nIndividualsTrait2));
    this->reml->V->appendVarianceToElement(names[i] + "_2", "Var(" + names[i] + "_2)", variance);
    
    this->reml->V->insertElement(names[i] + "_12", grmType, names[i] + "_12", 1., subMatrix(), subMatrix());
    this->reml->V->appendVarianceToElement(names[i] + "_12", "Covar(" + names[i] + "_12)", covariance);
  }

  
  this->reml->V->insertElement("E1", environtmentalType, "E1", 1., subMatrix(0, 0, nIndividualsTrait1, nIndividualsTrait1), subMatrix(0, 0, nIndividualsTrait1, nIndividualsTrait1));
  this->reml->V->appendVarianceToElement("E1", "Var(E1)", variance);
  
  this->reml->V->insertElement("E2", environtmentalType, "E2", 1., subMatrix(nIndividualsTrait1, nIndividualsTrait1, nIndividualsTrait2, nIndividualsTrait2), subMatrix(0, 0, nIndividualsTrait2, nIndividualsTrait2));
  this->reml->V->appendVarianceToElement("E2", "Var(E2)", variance);
  
  if(computeResidualCovariance)
  {
    this->reml->V->insertElement("E12", environtmentalType, "E12", 1., subMatrix(), subMatrix());
    this->reml->V->appendVarianceToElement("E12", "Covar(E12)", covariance);
  }
  
  grms.clear();
}

void BivarREML::compareREMLs(double baseLogLikelihood, double baseNVariances)
{
  if(communicator->mpiRoot)
  {
    double LogRatio = 2.0*(baseLogLikelihood - this->reml->logLikelihood);
    if(LogRatio < 0.)
    {
      LogRatio = 0.;
    }
    int df = baseNVariances - this->reml->V->variances.size();
    
    Message message(options.outFile + ".reml.fixed");
    
    message << "LRT\t" << std::setprecision(3) << LogRatio << std::endl;
    message << "df\t" << std::setprecision(1) << df << std::endl;
    message << "Pval\t" << std::setprecision(4) << 0.5*chi1_CDF(df, LogRatio) << std::endl;
  }
}