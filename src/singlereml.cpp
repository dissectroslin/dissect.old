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

#include "singlereml.h"
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

#include <sstream>
#include <vector>

SingleREML::SingleREML(bool argWriteResults)
{
  this->writeResults = argWriteResults;
  
  this->reml = NULL;
}

SingleREML::~SingleREML()
{
  if(this->reml != NULL)
  {
    delete this->reml;
    this->reml = NULL;
  }
}

void SingleREML::compute()
{
  if(this->reml != NULL)
  {
    delete this->reml;
  }
  
  //Prepare data
  Genotype *genotypes;
  Kernel *grmBase = loadGRMUsingOptions(true, &genotypes);
  precomputeInitialValues(0, grmBase);
  if(options.computeSNPsBLUP)
  {
    if(genotypes == NULL)
    {
      genotypes = loadGenotypeUsingOptions();
    }
  }
  else
  {
    if(genotypes != NULL)
    {
      delete genotypes;
      genotypes = NULL;
    }
  }
  
  Kernel *epistaticGRMBase = NULL;
  if( options.computeEpistasisVariance == true )
  {
    epistaticGRMBase = new Kernel(grmBase, kernelEpistaticGRM, genotypes);
  }
  
  //Compute REML
  std::vector<int> phenotypesForAnalyze = getPhenotyesForAnalysis();
  Kernel * grm;
  Kernel * epistaticGRM;
  std::string baseOutFile = options.outFile;
  for(int i = 0; i<phenotypesForAnalyze.size(); i++)
  {
    if(phenotypesForAnalyze.size() != 1)
    {
      grm = new Kernel(grmBase);
      std::stringstream sstemp;
      sstemp << phenotypesForAnalyze[i];
      options.outFile = baseOutFile + "." + sstemp.str();
      options.phenotypeColumn = phenotypesForAnalyze[i];
      if( options.computeEpistasisVariance == true )
      {
        epistaticGRM = new Kernel(epistaticGRMBase);
      }
    }
    else
    {
      grm = grmBase;
      grmBase = NULL;
      if( options.computeEpistasisVariance == true )
      {
        epistaticGRM = epistaticGRMBase;
        epistaticGRMBase = NULL;
      }
    }
  
    std::vector<Kernel*> grms;
    grm->name = "GRM";
    grms.push_back(grm);
    
    std::vector<double> weights;
    weights.push_back(1.);
    
    std::vector<int> phenotypeColumns;
    phenotypeColumns.push_back(options.phenotypeColumn);

    std::vector<double> heritabilities;
    heritabilities.push_back(options.initialh2Trait);
    
    std::vector<std::pair<std::string, std::string> > covariateFiles;
    covariateFiles.push_back(std::pair<std::string, std::string>(options.covarsFile, options.qCovarsFile));
  
    if( options.computeEpistasisVariance == true )
    {
      epistaticGRM->name = "epi";
      grms.push_back(epistaticGRM);
      weights.clear(); //weights will be computed by prepare function.
    }
    
    this->reml = new REML(this->writeResults);

    if(genotypes != NULL && options.computeSNPsBLUP)
    {
      this->reml->mSNPIds[grm->name] = grm->randomVarNames;
      this->reml->SNPsBLUPGenotypes[grm->name] = genotypes;
    }
    if(options.computeIndividualsBLUP == true)
    {
      this->reml->individualBLUPNames.push_back(grm->name);
    }
    
    bool prepared = this->reml->prepare(singleREMLType, grms, weights, phenotypeColumns, heritabilities, covariateFiles);
    if( prepared == true )
    {
      this->reml->computeREML();
    }
    else
    {
      misc.message << "Sorry, a problem was happened while preparing data for performing REML. The MLM cannot be fitted. Please, check the logs." << std::endl;
    }

    this->reml->SNPsBLUPGenotypes.clear();    
    delete this->reml;
    this->reml = NULL;

  }
  
  if(genotypes != NULL)
  {
    delete genotypes;
  }
  if(grmBase != NULL)
  {
    delete grmBase;
  }
  if(epistaticGRMBase != NULL)
  {
    delete epistaticGRMBase;
  }
}

void SingleREML::computeRegional()
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
    globalGRM->name = "Global GRM";
    grms.push_back(globalGRM);
    regionalGRM->name = "Regional GRM";
    grms.push_back(regionalGRM);
    
    std::vector<double> weights;
    weights.push_back(1. - proportionRegionalTotal);
    weights.push_back(proportionRegionalTotal);
    
    std::vector<int> phenotypeColumns;
    phenotypeColumns.push_back(options.phenotypeColumn);

    std::vector<double> heritabilities;
    heritabilities.push_back(options.initialh2Trait);
    
    std::vector<std::pair<std::string, std::string> > covariateFiles;
    covariateFiles.push_back(std::pair<std::string, std::string>(options.covarsFile, options.qCovarsFile));
    
    if(options.computeIndividualsBLUP == true)
    {
      this->reml->individualBLUPNames.push_back(globalGRM->name);
      this->reml->individualBLUPNames.push_back(regionalGRM->name);
    }
    
    bool prepared = this->reml->prepare(regionalSingleREMLType, grms, weights, phenotypeColumns, heritabilities, covariateFiles);
    if( prepared == true )
    {
      this->reml->computeREML();
    }
    else
    {
      misc.message << "WARNING: Sorry, a problem was happened while preparing data for performing REML on region " << group << ". The MLM cannot be fitted on this region. Please, check the logs." << std::endl;
    }
    
    delete this->reml;
    this->reml = NULL;
  }
  
  delete grm;
  delete genotype;
}

void SingleREML::computeMultipleGroups()
{
  if(this->reml != NULL)
  {
    delete this->reml;
  }
  
  misc.error("Error: Not tested. Probably global GRM should not be specified.", 0);
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
  
  misc.error("Error: Here appropiate checks must be performed. i.e. SNPs if global GRM? Diagonalized global GRM? etc.", 0);
  if(grm->individualIds != genotype->individualIds)
  {
    misc.error("Error: The individuals in GRM and genotypes are not the same. Maybe are you not using the proper GRM or genotype file?", 0);
  }
  
  //Create groups
  genotype->groupSNPs(options.regionBy);
  
  //Start iterations over groups
  misc.message << "Computing all regional GRMs..." << std::endl;
  misc.setGetElapsedTime("ConstructRegionalGRMs");
  
  std::vector<Kernel*> grms;
  std::vector<double> weights;
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
    
    Genotype *regionalGenotype = new Genotype();
    genotype->genotypeOfSNPsGroup(group, regionalGenotype);
    double proportionRegionalTotal = double(regionalGenotype->nSNPs)/double(genotype->nSNPs);
    Kernel *regionalGRM = new Kernel(regionalGenotype);
    delete regionalGenotype;
    
    regionalGRM->name = group;
    grms.push_back(regionalGRM);
    weights.push_back(proportionRegionalTotal);
  }
  delete grm;
  delete genotype;
  
  misc.message << "Regional GRMs computed after " << misc.setGetElapsedTime("ConstructRegionalGRMs", true) << std::endl;
  
  this->reml = new REML();
  
  std::vector<int> phenotypeColumns;
  phenotypeColumns.push_back(options.phenotypeColumn);

  std::vector<double> heritabilities;
  heritabilities.push_back(options.initialh2Trait);
  
  std::vector<std::pair<std::string, std::string> > covariateFiles;
  covariateFiles.push_back(std::pair<std::string, std::string>(options.covarsFile, options.qCovarsFile));
  
  bool prepared = this->reml->prepare(regionalSingleREMLType, grms, weights, phenotypeColumns, heritabilities, covariateFiles);
  if( prepared == true )
  {
    this->reml->computeREML();
  }
  else
  {
    misc.message << "Sorry, a problem was happened while preparing data for performing REML. The MLM cannot be fitted. Please, check the logs." << std::endl;
  }
    
  delete this->reml;
  this->reml = NULL;
}

void SingleREML::computeRaw(Kernel *srcGRM)
{
  std::vector<Kernel*> grms;
  srcGRM->name = "GRM";
  grms.push_back(srcGRM);
  
  std::vector<double> weights;
  weights.push_back(1.);
  
  std::vector<int> phenotypeColumns;
  phenotypeColumns.push_back(options.phenotypeColumn);

  std::vector<double> heritabilities;
  heritabilities.push_back(options.initialh2Trait);
  
  std::vector<std::pair<std::string, std::string> > covariateFiles;
  covariateFiles.push_back(std::pair<std::string, std::string>(options.covarsFile, options.qCovarsFile));
  
  bool prepared = this->reml->prepare(singleREMLType, grms, weights, phenotypeColumns, heritabilities, covariateFiles);
  if( prepared == true )
  {
    this->reml->computeREML();
  }
  else
  {
    misc.message << "Sorry, a problem was happened while preparing data for performing REML. The MLM cannot be fitted. Please, check the logs." << std::endl;
  }
}

/*
void SingleREML::prepare(Kernel * grm, Genotype *genotypes)
{
  //Set reml type
  this->reml->type = singleREMLType;
  
  //Copy SNPs Ids used for computing the GRM
  this->reml->SNPIds = grm->randomVarNames;
  
  //If computing SNP BLUPs, set genotype data.
  if(options.computeSNPsBLUP)
  {
    this->reml->SNPsBLUPGenotypes = genotypes;
  }
  
  //Prepare the GRM (covs), covariance (X) and phenotype (y) matrices
  Phenotype * phenotype = new Phenotype(cyclicDistribution, options.phenotypesFile, options.phenotypeColumn1);
  Covariate * covariate = new Covariate(options.covarsFile1, options.qCovarsFile1, phenotype->individualIds);
  //Covariate * covariate = new Covariate("", options.qCovarsFile1);
  this->reml->meanNames.clear();
  this->reml->quantitativeCovarNames.clear();
  this->reml->discreteCovarNames.clear();
  this->reml->addCovariatesNames(covariate);
  
  std::vector<std::string> commonIndividuals, commonIndividualsInGRMOrder;
  if(communicator->mpiRoot)
  {
    commonIndividuals = intersectionStringVectors(3, &grm->individualIds, &phenotype->individualIds, &covariate->individualIds);
    commonIndividualsInGRMOrder = orderVectorAsTemplate(grm->individualIds, commonIndividuals);
  }
  misc.message << commonIndividualsInGRMOrder.size() << " individuals are kept for the analysis." << std::endl;
  
  phenotype->filterIndividuals(commonIndividualsInGRMOrder);
  covariate->filterIndividuals(commonIndividualsInGRMOrder);
  grm->filterIndividuals(commonIndividualsInGRMOrder, false);
  this->reml->individuals = grm->individuals;
  
  double phenotypeVariance = phenotype->computePhenotypeVariance();
  this->reml->y = new Matrix(phenotype->phenotypes);
  this->reml->X = new Matrix(covariate->covariates);
  delete phenotype;
  delete covariate;
  
  //Start covariance matrix
  int tempDimension = commonIndividualsInGRMOrder.size();
  communicator->broadcast(&tempDimension);
  this->reml->dimension = tempDimension;
  
  this->reml->V = new CovarianceMatrix(this->reml->dimension, grm->diagonalized);
  
  //Correct corresponfing matrices if grm is diagonalized
  
  if( grm->diagonalized == true )
  {
    Matrix * temp;
    
    temp = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
    temp->multiply(grm->eigenVectors, 'T', this->reml->y, 'N');
    delete this->reml->y;
    this->reml->y = temp;
    
    temp = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
    temp->multiply(grm->eigenVectors, 'T', this->reml->X, 'N');
    delete this->reml->X;
    this->reml->X = temp;
    
    if(this->reml->SNPsBLUPGenotypes != NULL)
    {
      this->reml->SNPsBLUPGenotypes->normalizeGenotypes();
      temp = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
      temp->multiply(this->reml->SNPsBLUPGenotypes->genotypes, 'N', grm->eigenVectors, 'N');
      delete this->reml->SNPsBLUPGenotypes->genotypes;
      this->reml->SNPsBLUPGenotypes->genotypes = temp;
    }
    
    this->reml->GRMEigenVectors = new Matrix(grm->eigenVectors);
  }
  
  //Add elements on covariance matrix
  
  if(options.correctLinkageDisequilibrium == true) //Replace current grm matrix by the grm*grm/N product matrix
  {
    double tempFactor;
    if(this->reml->SNPsBLUPGenotypes != NULL) //Replace current genotype matrix by the genotype*grm/N product matrix
    {
      misc.message << "Correcting genotypes for individual BLUPs for LD..." << std::endl;
      this->reml->SNPsBLUPGenotypes->normalizeGenotypes();
      Matrix * temp = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
      grm->normalize();
      tempFactor = double(grm->randomVarNames.size());
      communicator->broadcast(&tempFactor);
      temp->multiply(this->reml->SNPsBLUPGenotypes->genotypes, 'N', grm->getNormalizedKernel(), 'N', tempFactor);
      delete this->reml->SNPsBLUPGenotypes->genotypes;
      this->reml->SNPsBLUPGenotypes->genotypes = temp;
    }
    
    misc.message << "Correcting grm for LD..." << std::endl;
    Matrix * temp = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
    grm->normalize();
    tempFactor = double(grm->randomVarNames.size())/double(grm->nIndividuals);
    communicator->broadcast(&tempFactor);
    temp->multiply(grm->getNormalizedKernel(), 'N', grm->getNormalizedKernel(), 'N', tempFactor);
    this->reml->V->insertCovarianceMatrix("GRM", temp);
    delete temp;
    delete grm;
  }
  else
  {
    this->reml->V->insertCovarianceMatrix("GRM", grm);
    delete grm;
  }
  
  this->reml->V->insertVarianceGroup("Phenotype", phenotypeVariance);
  
  this->reml->V->insertVariance("Var(GRM)", "Phenotype", variance, genetic, phenotypeVariance*options.initialh2Trait1);
  this->reml->V->insertVariance("Var(E)", "Phenotype", variance, environment, phenotypeVariance*(1.-options.initialh2Trait1));
  //this->reml->V->variances[0] = phenotypeVariance*0.5; //0.375;
  //this->reml->V->variances[1] = phenotypeVariance*0.5; //0.375;
  
  this->reml->V->insertElement("GRM", grmType, "GRM", 1., subMatrix(), subMatrix());
  this->reml->V->appendVarianceToElement("GRM", "Var(GRM)", variance);
  
  this->reml->V->insertElement("Environmental", environtmentalType, "Identity", 1., subMatrix(), subMatrix());
  this->reml->V->appendVarianceToElement("Environmental", "Var(E)", variance);
  
  //Set initial values if specified somwhow
  this->reml->V->setVarianceInitialValues(this->subsampleREMLResults.variances);
  
}

void SingleREML::prepare(std::vector<Kernel*> & grms, std::vector<std::string> names, std::vector<double> weights)
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
  
  this->reml->type = multipleSingleREMLType;
  
  //Prepare SNPs
  
  this->reml->SNPIds.clear();
  
  //Prepare, covariance (X) and phenotype (y) matrices
  
  Phenotype * phenotype = new Phenotype(cyclicDistribution, options.phenotypesFile, options.phenotypeColumn1);
  Covariate * covariate = new Covariate(options.covarsFile1, options.qCovarsFile1, phenotype->individualIds);
  this->reml->meanNames.clear();
  this->reml->quantitativeCovarNames.clear();
  this->reml->discreteCovarNames.clear();
  this->reml->addCovariatesNames(covariate);
  
  std::vector<std::string> commonIndividuals, commonIndividualsInGRMOrder;
  if(communicator->mpiRoot)
  {
    commonIndividuals = intersectionStringVectors(3, &(grms[0]->individualIds), &phenotype->individualIds, &covariate->individualIds);
    for(int i = 1; i<grms.size(); i++)
    {
      commonIndividuals = intersectionStringVectors(2, &commonIndividuals, &(grms[i]->individualIds));
    }
    commonIndividualsInGRMOrder = orderVectorAsTemplate(grms[0]->individualIds, commonIndividuals);
    for(int i = 1; i<grms.size(); i++)
    {
      std::vector<std::string> tempTest = orderVectorAsTemplate(grms[i]->individualIds, commonIndividuals);
      if( tempTest != commonIndividualsInGRMOrder )
      {
        misc.error("Error: Not all GRMs have the individuals in the same order. At this moment, DISSECT cannot resort it. Please, resort it or contact us.", 0);
      }
    }
  }
  misc.message << commonIndividualsInGRMOrder.size() << " individuals are kept for the analysis." << std::endl;
  
  phenotype->filterIndividuals(commonIndividualsInGRMOrder);
  covariate->filterIndividuals(commonIndividualsInGRMOrder);
  for(int i = 0; i<grms.size(); i++)
  {
    grms[i]->normalize();
    grms[i]->filterIndividuals(commonIndividualsInGRMOrder, false);
  }
  this->reml->individuals = grms[0]->individuals;
  
  double phenotypeVariance = phenotype->computePhenotypeVariance();
  this->reml->y = new Matrix(phenotype->phenotypes);
  this->reml->X = new Matrix(covariate->covariates);
  delete phenotype;
  delete covariate;
  
  //Start covariance matrix
  int tempDimension = commonIndividualsInGRMOrder.size();
  communicator->broadcast(&tempDimension);
  this->reml->dimension = tempDimension;
  
  this->reml->V = new CovarianceMatrix(this->reml->dimension);
  
  for(int i = 0; i<grms.size(); i++)
  {
    this->reml->V->insertCovarianceMatrix(names[i], grms[i]);
    delete grms[i];
  }
  
  ////////////////////////////////////////////////////
  // Insert variance groups
  this->reml->V->insertVarianceGroup("Phenotype", phenotypeVariance);
  
  ////////////////////////////////////////////////////
  // Insert variances
  for(int i = 0; i<grms.size(); i++)
  {
    std::string tempName = "Var(" + names[i] + ")";
    this->reml->V->insertVariance(tempName, "Phenotype", variance, genetic, phenotypeVariance*options.initialh2Trait1*weights[i]);
  }
  
  this->reml->V->insertVariance("Var(E)", "Phenotype", variance, environment, phenotypeVariance*(1.-options.initialh2Trait1));
  
  ////////////////////////////////////////////////////
  // Insert elements
  for(int i = 0; i<grms.size(); i++)
  {
    this->reml->V->insertElement(names[i], grmType, names[i], 1., subMatrix(), subMatrix());
    std::string tempName = "Var(" + names[i] + ")";
    this->reml->V->appendVarianceToElement(names[i], tempName, variance);
  }
  
  this->reml->V->insertElement("Environmental", environtmentalType, "Identity", 1., subMatrix(), subMatrix());
  this->reml->V->appendVarianceToElement("Environmental", "Var(E)", variance);
  
  grms.clear();
}*/

void SingleREML::precomputeInitialValues(int type, Kernel *srcGRM)
{
  REMLResults results;
  results.variances.clear();
  results.logLikelihood = 0;
  
  if(options.initialVariancesFile != "" )
  {
    misc.message << "Reading initial variance values from file [ " << options.initialVariancesFile << " ]" << std::endl;
    if( communicator->mpiRoot == true )
    {
      std::ifstream file;
      std::string line;
      
      misc.checkFileExists(options.initialVariancesFile);
      file.open(options.initialVariancesFile.c_str());
      
      while(getline(file,line))
      {
        if(!file)
        {
          break;
        }
        
        std::istringstream sstemp(line);
        
        Variance variance;
        
        sstemp >> variance.name;
        if( (sstemp >> variance.variance).fail() )
        {
          misc.error("Error: The variance " + variance.name + " has not a valid value in file [ " + options.initialVariancesFile + " ].", 0);
        }
        results.variances.push_back(variance);
      }
      file.close();
    }
  }
  else if(options.computeREMLInSubsample == true)
  {
    //If number of individuals is not enough, dissable subsampling.
    int avoidPrecomputing = 0;
    if(3*options.minimumSubsample > srcGRM->nIndividuals)
    {
      avoidPrecomputing = 1;
    }
    communicator->broadcast(&avoidPrecomputing);
    if(avoidPrecomputing == 1)
    {
      misc.message << "Random subsampling is disabled due the low number of individuals." << std::endl;
      this->subsampleREMLResults = results;
      return;
    }
    
    misc.message << "Estimating starting variance values from REML on random subsampling." << std::endl;
    misc.message.tab = "  ";
    misc.setGetElapsedTime("REMLSubsampling");
    
    //Start random subsampling
    std::vector<REMLResults> partialREMLResults;
    for(int i = 0; i < options.nSubSampleIterations; i++)
    {
      SingleREML singleREML(false);
      if(type == 0)
      {
        Kernel *grm = new Kernel(srcGRM);
        grm->randomSubSample(options.initialSubsampleFraction, options.minimumSubsample);
        singleREML.computeRaw(grm);
      }
      
      if(singleREML.reml != NULL)
      {
        if(singleREML.reml->success == true) //If REML finished successful, store results.
        {
          REMLResults temp;
          temp.variances = singleREML.reml->V->variances;
          temp.logLikelihood = singleREML.reml->logLikelihood;
          partialREMLResults.push_back(temp);
        }
      }
    }
    
    if(partialREMLResults.size() > 0) //Average results from all REML analysis
    {
      results.variances = partialREMLResults[0].variances;
      results.logLikelihood = partialREMLResults[0].logLikelihood;
      for(int i = 1; i < partialREMLResults.size(); i++) //Average variances
      {
        std::vector<Variance> variancesToAdd = partialREMLResults[i].variances;
        if(variancesToAdd.size() != results.variances.size())
        {
          misc.error("Error: An internal error was happened when estimating initial variances using random subsampling. Different number of variances.", 0);
        }
        for(int j = 0; j < results.variances.size(); j++)
        {
          if(results.variances[j].name != variancesToAdd[j].name)
          {
            misc.error("Error: An internal error was happened when estimating initial variances using random subsampling. Different variance names.", 0);
          }
          results.variances[j].variance += variancesToAdd[j].variance;
        }
      }
      
      for(int j = 0; j < results.variances.size(); j++)
      {
        results.variances[j].variance /= double(partialREMLResults.size());
      }
    } //End of results averaging
    
    misc.message.tab = "";
    misc.message << "Initial variance values have been estimated after " << misc.setGetElapsedTime("REMLSubsampling", true) << std::endl;
  }
  
  this->subsampleREMLResults = results;
}
