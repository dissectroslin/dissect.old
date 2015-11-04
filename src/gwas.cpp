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

#include "gwas.h"
#include "analysis.h"
#include "reml.h"
#include "auxiliar.h"

#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstdio>

GWAS::GWAS()
{
  this->currentFile = "";
  this->currentGRMFile = "";
  this->currentGRMBase = NULL;
  this->currentGRMFiltered = NULL;
  this->accumulatedh2 = options.initialh2Trait;
  this->nTests = 0;
  
  if(options.genotypeFile != "")
  {
    this->genotypeFiles.push_back(options.genotypeFile);
  }
  else if(options.genotypeListFile != "")
  {
    getListFromFile(options.genotypeListFile, this->genotypeFiles);
  }
  else if(options.genotypeAndGRMsListFile != "")
  {
    std::vector< std::vector<std::string> > files;
    getTableFromFile(options.genotypeAndGRMsListFile, files, 2);
    std::map<std::string, std::vector<std::string> > grms2genotypes;
    for(int i = 0; i<files.size(); i++) //We load genotype and grm files in this way, because we want all genotype with same grms together in the list.
    {
      misc.checkFileExists(files[i][0] + ".bed");
      misc.checkFileExists(files[i][0] + ".bim");
      misc.checkFileExists(files[i][0] + ".fam");
      misc.checkFileExists(files[i][1] + ".grm.ids");
      misc.checkFileExists(files[i][1] + ".grm.dat");
      misc.checkFileExists(files[i][1] + ".grm.snps");
      grms2genotypes[ files[i][1] ].push_back( files[i][0] );
    }
    for( std::map<std::string, std::vector<std::string> >::iterator it = grms2genotypes.begin(); it != grms2genotypes.end(); ++it )
    {
      std::vector<std::string> genoFiles = it->second;
      std::string grmFile = it->first;
      for(int i = 0; i<genoFiles.size(); i++)
      {
        if( this->grmFiles.count(genoFiles[i]) != 0 )
        {
          misc.error("Error: There is at least one genotype file repeated in [ " + options.genotypeAndGRMsListFile + " ] file. Genotype files have to be unique.", 0);
        }
        this->genotypeFiles.push_back(genoFiles[i]);
        this->grmFiles[ genoFiles[i] ] = grmFile;
      }
    }
  }
  else
  {
    misc.error("Error: An internal error was happened when loading the genotype for a GWAS analysis. No files specified.", 0);
  }
  
}

GWAS::~GWAS()
{
}

void GWAS::computeGWAS()
{
  std::string backupOutFile = options.outFile;
  if(this->genotypeFiles.size() > 1)
  {
    misc.message.tab = "  ";
  }
  
  for(int i = 0; i<this->genotypeFiles.size(); i++)
  {
    //Load needed data: genotypes, covars, phenotypes and grms
    this->currentFile = this->genotypeFiles[i];
    Genotype * genotype = new Genotype(this->genotypeFiles[i]);
    Phenotype * phenotype = new Phenotype(cyclicDistribution, options.phenotypesFile, options.phenotypeColumn);
    Covariate * covariate = new Covariate(options.covarsFile, options.qCovarsFile, phenotype->individualIds);
  
    genotype->normalizeGenotypes();
    
    std::vector<std::string> commonIndividualsInGenotypeOrder;
    if(this->grmFiles.count( this->currentFile ) == 0)
    {
      //Get shared individuals and SNPs in group and filter
      std::vector<std::string> commonIndividuals = intersectionStringVectors(3, &genotype->individualIds, &phenotype->individualIds, &covariate->individualIds);
      commonIndividualsInGenotypeOrder = orderVectorAsTemplate(genotype->individualIds, commonIndividuals);
    }
    else
    {
      if( this->grmFiles[this->currentFile] != this->currentGRMFile ) //Only load a new GRM if file changes.
      {
        if(this->currentGRMBase != NULL)
        {
          delete this->currentGRMBase;
        }
        this->currentGRMFile = this->grmFiles[this->currentFile];
        this->currentGRMBase = new Kernel(this->currentGRMFile);
        this->currentGRMBase->normalize();
      }
      if( this->currentGRMFiltered != NULL )
      {
        delete this->currentGRMFiltered;
      }
      this->currentGRMFiltered = new Kernel(this->currentGRMBase);
      this->accumulatedh2 = options.initialh2Trait;
      this->nTests = 0;
      
      //Get shared individuals and SNPs in group and filter
      std::vector<std::string> commonIndividuals = intersectionStringVectors(4, &genotype->individualIds, &phenotype->individualIds, &covariate->individualIds, &this->currentGRMFiltered->individualIds);
      commonIndividualsInGenotypeOrder = orderVectorAsTemplate(genotype->individualIds, commonIndividuals);
      std::vector<std::string> commonIndividualsInGRMOrder = orderVectorAsTemplate(this->currentGRMFiltered->individualIds, commonIndividuals);
      if( commonIndividualsInGenotypeOrder != commonIndividualsInGRMOrder )
      {
        misc.error("Error: At this current version, when performing this analysis, individuals in genotype files must be in the same order than in grm files. Aborting.", 0);
      }
      if( this->currentGRMFiltered->diagonalized == true && commonIndividualsInGRMOrder != this->currentGRMFiltered->individualIds )
      {
        misc.error("Error: When using diagonal GRMs, the individuals in phenotype, genotype, covars and GRM files must be all the same without missing phenotypes. Aborting.", 0);
      }
      this->currentGRMFiltered->filterIndividuals(commonIndividualsInGRMOrder, false);
      
      //If GRM is diagonal, correct phenotypes, genotypes and covariates
      if(this->currentGRMFiltered->diagonalized == true)
      {
        Matrix * tempPhenos = new Matrix();
        tempPhenos->multiply(this->currentGRMFiltered->eigenVectors, 'T', phenotype->phenotypes, 'N');
        delete phenotype->phenotypes;
        phenotype->phenotypes = tempPhenos;
        
        Matrix * tempCovars = new Matrix();
        tempCovars->multiply(this->currentGRMFiltered->eigenVectors, 'T', covariate->covariates, 'N');
        delete covariate->covariates;
        covariate->covariates = tempCovars;
        
        Matrix * tempGenos = new Matrix();
        tempGenos->multiply(genotype->genotypes, 'N', this->currentGRMFiltered->eigenVectors, 'N'); //Geno*EigV = (EigVt*Genot)t
        delete genotype->genotypes;
        genotype->genotypes = tempGenos;
      }
    }
    
    phenotype->filterIndividuals(commonIndividualsInGenotypeOrder);
    covariate->filterIndividuals(commonIndividualsInGenotypeOrder);
    genotype->filterSNPsAndIndividuals(genotype->SNPIds, commonIndividualsInGenotypeOrder, false);
    
    if( this->genotypeFiles.size() > 1 )
    {
      options.outFile += "." + getFileName(this->genotypeFiles[i]);
    }
    
    //Start analysis
    if(options.analysis != recursiveGWASAnalysis) //Perform grouped GWAS or standard GWAS
    {
      misc.setGetElapsedTime("GWAS");
      misc.message << "Starting GWAS analysis..." << std::endl;
      if(options.regionBy != ungrouped)
      {
        genotype->groupSNPs(options.regionBy);
        computeGroupedGWAS(genotype, phenotype, covariate);
      }
      else
      {
        computeIndividualGWAS(genotype, phenotype, covariate);
      }
      misc.message << "GWAS analysis finished after " << misc.setGetElapsedTime("GWAS", true) << "." << std::endl;
    }
    else //Perform recursive GWAS
    {
      misc.setGetElapsedTime("recursiveGWAS");
      this->significanceThreshold = options.significanceThresholdFilterSNPs;
      this->significantSNPIds.clear();
      std::vector<std::string> subsetToAnalyzeSNPIds = genotype->SNPIds;
      int iteration = 1;
      std::string backupOutFile2 = options.outFile;
      options.fixedGroupSize = floor(double(genotype->nIndividuals)*options.relationFitSNPsIndividuals);
      if(options.fixedGroupSize < 1)
      {
        misc.error("Error: The resultant SNPs group size for a recursive GWAS is less than 1. This could happen because sample size is too small or because the parameter passed by --rgwas-ratio is too small.", 0);
      }
      misc.message << "Starting recursive GWAS analysis using groups of " << options.fixedGroupSize << " SNPs..." << std::endl;
      misc.message.tab = "  ";
      while( misc.gt( subsetToAnalyzeSNPIds.size() != this->significantSNPIds.size() && (options.recursiveGWASMaxIterations < 1 || options.recursiveGWASMaxIterations >= iteration) ) )
      {
        options.outFile = backupOutFile2 + ".iter" + getString(iteration);
        subsetToAnalyzeSNPIds = this->significantSNPIds;
        misc.message << "Performing iteration " << iteration << "." << std::endl;
        
        genotype->groupSNPs(byOrderedFixedSize, subsetToAnalyzeSNPIds);
        computeGroupedGWAS(genotype, phenotype, covariate);
        
        misc.message << "Analysis finished. " << ((iteration==1)?genotype->nSNPs:subsetToAnalyzeSNPIds.size()) << " SNPs analyzed." << std::endl;
        
        iteration++;
        
        if( misc.gt(this->significantSNPIds.size() == 0) )
        {
          misc.message << "No significant SNPs found. Stopping iterations." << std::endl;
          break;
        }
      }
      options.outFile = backupOutFile2;
      misc.message.tab = "";
      misc.message << "recursive GWAS analysis finished after " << iteration - 1 << " iterations. It needed " << misc.setGetElapsedTime("recursiveGWAS", true) << "." << std::endl;
    }
    
    options.outFile = backupOutFile;
    
    delete genotype;
    delete phenotype;
    delete covariate;
  }
  
  this->currentFile = "";
  
  if(this->currentGRMBase != NULL)
  {
    delete this->currentGRMBase;
    this->currentGRMBase = NULL;
  }
  if( this->currentGRMFiltered != NULL )
  {
    delete this->currentGRMFiltered;
    this->currentGRMFiltered = NULL;
  }
  
  if(this->genotypeFiles.size() > 1)
  {
    misc.message.tab = "";
    misc.message << "GWAS analysis finished on all files." << std::endl;
  }
}

void GWAS::computeGroupedGWAS(Genotype * genotype, Phenotype * phenotype, Covariate * covariate)
{
  std::map<std::string, GLMResults > results;
  std::map<std::string, std::vector<SNP> > resultSNPInfo;
  
  std::map<std::string, std::vector<std::string> > groupsWithErrors;
 
  int nGroups = genotype->groupedSNPs.size();
  communicator->broadcast(&nGroups);
  std::map<std::string, std::set<std::string> >::iterator it = genotype->groupedSNPs.begin();
  
  GLMResults reducedResults;
  Matrix * br = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  computeGLM(phenotype->phenotypes, covariate->covariates, br, reducedResults);
  delete br;
  
  //Start iterations over groups
  for(int i = 0; i < nGroups; i++)
  {
    std::string group = "";
    if(communicator->mpiRoot)
    {
      group = it->first;
      it++;
    }
    
    Genotype *groupGenotype = new Genotype();      
    genotype->genotypeOfSNPsGroup(group, groupGenotype, false);
    
    bool iterate = true;
    while(iterate == true)
    {
      //Set the analysis matrices
      Matrix * gTranspose = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
      gTranspose->transpose(groupGenotype->genotypes);
      Matrix * incidence = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
      incidence->joinMatricesHorizontally(covariate->covariates, gTranspose);

      //Fit GLM
      GLMResults groupResults;
      Matrix * fixedEffects = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
      bool success = computeGLM(phenotype->phenotypes, incidence, fixedEffects, groupResults);
      
      //If success, collect results
      if(success == true)
      {
        bool succesGroupSignificance = computeGroupSignificance(phenotype->phenotypes, reducedResults, groupResults);
        //If group significance can be calculated, collect results. Otherwise, set success to false.
        if(succesGroupSignificance == true)
        {
          if(communicator->mpiRoot)
          {
            results[group] = groupResults;
            resultSNPInfo[group] = groupGenotype->SNPs;
          }
          getLessSignificantCorrelatedSNPs(options.correlatedSNPsThreshold, groupGenotype, groupResults, covariate->nCovariates); //Get SNPs correlated with other SNPs more significant.
          iterate = false;
        }
        else
        {
          success = false;
        }
      }
      
      //If fitting failed, try to filter dependent SNPs in case their are not already filtered.
      if(success == false)
      {
        if( misc.gt( groupsWithErrors.count(group) != 0 ) ) //This group was already filtered for dependent SNPs?
        {
          iterate = false;
        }
        else
        {
          std::vector<int> dependentColumns = gTranspose->getDependentColumns();
          if( misc.gt(dependentColumns.size() == 0) ) //Are there any dependent SNP?
          {
            iterate = false;
          }
          else //Filter dependent SNPs
          {
            std::vector<std::string> dependentSNPs;
            for(int i = 0; i<dependentColumns.size(); i++)
            {
              dependentSNPs.push_back( groupGenotype->SNPIds[ dependentColumns[i] ] );
            }
            misc.message << "Sorry, linear model for group < " << group << " > in file [ " << this->currentFile << " ] cannot be fitted. Maybe there are linear dependent SNPs. Trying to filter dependent SNPs." << std::endl;
            std::vector<std::string> independentSNPs = differenceBetweenTwoVectors(groupGenotype->SNPIds, dependentSNPs);
            groupGenotype->filterSNPsAndIndividuals(independentSNPs, groupGenotype->individualIds);
            groupsWithErrors[group] = dependentSNPs;
          }
          if( misc.gt(groupGenotype->nSNPs == 0) )
          {
            iterate = false;
          }
        }
        if(iterate == false) //No solution with current method. Store all group in the file with unfitted SNPs.
        {
          misc.message << "Sorry, linear model for group < " << group << " > in file [ " << this->currentFile << " ] cannot be fitted even after filtering problematic SNPs. Could there are some linear dependence between covariates?" << std::endl;
          groupsWithErrors[group] = groupGenotype->SNPIds;
        }
        groupResults = GLMResults();
      }
      
      delete gTranspose;
      delete incidence;
      delete fixedEffects;
    }
    
    delete groupGenotype;      
  }
  storeResults(results, covariate, resultSNPInfo);
  
  if(communicator->mpiRoot == true && groupsWithErrors.size() != 0)
  {
    misc.message << "There are " << groupsWithErrors.size() << " groups that cannot be fitted or with dependent SNPs. They are stored in file [ " << (options.outFile + ".multi.gwas.unfitted") << " ]." << std::endl;
    Message message(options.outFile + ".multi.gwas.unfitted");
    for(std::map<std::string, std::vector<std::string> >::iterator it = groupsWithErrors.begin(); it != groupsWithErrors.end(); ++it)
    {
      std::string group = it->first;
      std::vector<std::string> unfittedSNPs = it->second;
      for (int i = 0; i<unfittedSNPs.size(); i++)
      {
          message << group << " " << unfittedSNPs[i] << std::endl;
      }
    }
  }
}

void GWAS::computeIndividualGWAS(Genotype * genotype, Phenotype * phenotype, Covariate * covariate)
{

  std::map<std::string, GLMResults > results;
  std::map<std::string, std::vector<SNP> > resultSNPInfo;
  std::vector<std::string> unfittedSNPs;
  
  Matrix * gTransposed = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  gTransposed->transpose(genotype->genotypes);
  
  for(int isnp = 0; isnp < genotype->nSNPs; isnp++ )
  {
    //Get data for one SNP
    Matrix * snpData = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, genotype->nIndividuals, 1, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
    snpData->fillWithConstant(0.);
    snpData->add(gTransposed, 0., 1., subMatrix(0, 0, snpData->nGlobRows, 1), subMatrix(0, isnp, gTransposed->nGlobRows, 1));
    
    Matrix * incidence = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
    incidence->joinMatricesHorizontally(covariate->covariates, snpData);
    
    GLMResults SNPResults;
    Matrix * fixedEffects = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
    bool success = computeGLM(phenotype->phenotypes, incidence, fixedEffects, SNPResults);

    //If success, collect results
    if(success == false)
    {
      unfittedSNPs.push_back(genotype->SNPs[isnp].name);
    }
    else
    {
      if(communicator->mpiRoot)
      {
        std::vector<SNP> temp;
        temp.push_back(genotype->SNPs[isnp]);
        results[ genotype->SNPs[isnp].name ] = SNPResults;
        resultSNPInfo[ genotype->SNPs[isnp].name ] = temp;
      }
    }
    
    delete snpData;    
    delete incidence;
    delete fixedEffects;
  }
  delete gTransposed;
  
  storeResults(results, covariate, resultSNPInfo);
  
  if(communicator->mpiRoot && unfittedSNPs.size() != 0)
  {
    misc.message << "Warning: There are " << unfittedSNPs.size() << " SNPs which cannot be fitted. They are stored in file [ " << (options.outFile + ".gwas.unfitted") << " ]." << std::endl;
    Message message(options.outFile + ".gwas.unfitted");
    for(int i = 0; i < unfittedSNPs.size(); i++)
    {
      message << unfittedSNPs[i] << std::endl;
    }
  }
}

bool GWAS::computeGLM(Matrix * y, Matrix * X, Matrix * b, GLMResults & results)
{
  if( this->grmFiles.count( this->currentFile ) == 0 )
  {
    return computeGLMWithoutCovariance(y, X, b, results);
  }
  else
  {
    return computeGLMWithCovariance(y, X, b, results);
  }
}

bool GWAS::computeGLMWithoutCovariance(Matrix * y, Matrix * X, Matrix * b, GLMResults & results)
{
  if(X->nGlobCols > X->nGlobRows)
  {
    misc.message << "Warning: There are more fixed effects than data points. I am unable to fit the model." << std::endl;
    return false;
  }
  
  //Adjust linear model
  Matrix * XtX_i = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  XtX_i->multiply(X, 'T', X, 'N');
  
  bool inverted = XtX_i->symmetricInvert();
  if(inverted == false)
  {
    //XtX_i->multiply(X, 'T', X, 'N');
    //inverted = XtX_i->invert();
    //if(inverted == false)
    //{
    //  return false;
    //}
    results = GLMResults();
    return false;
  }
  
  Matrix * Xty = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  Xty->multiply(X, 'T', y, 'N');
  
  b->multiply(XtX_i, 'N', Xty, 'N');
  
  //Compute some matrices for SSE and t-statistics estimation
  Matrix * btXty = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  btXty->multiply(b, 'T', Xty, 'N');
  Matrix * yty = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  yty->multiply(y, 'T', y, 'N');
  
  if(btXty->nGlobRows != 1 || btXty->nGlobCols != 1 || yty->nGlobRows != 1 || yty->nGlobCols != 1)
  {
    misc.error("Error: An internal error was happened when computing GLM.", 0);
  }
  
  //Collect results
  std::vector<double> globalbtXty;
  std::vector<double> globalyty;
  std::vector<double> globalXtX_iDiagonal;
  results.SE.clear();
  results.tStatistics.clear();
  results.tStatisticPValues.clear();
  
  b->matrixToStandardVector(results.b);
  btXty->matrixToStandardVector(globalbtXty);
  yty->matrixToStandardVector(globalyty);
  globalXtX_iDiagonal = XtX_i->diagonal();
  
  if(communicator->mpiRoot)
  {
    double n_q_1 = double(y->nGlobRows) - double(b->nGlobRows);
    
    results.btXty = globalbtXty[0];
    
    results.SSE = globalyty[0] - results.btXty;
    results.MSE = results.SSE/double(n_q_1);

    for(int i = 0; i < results.b.size(); i++)
    {
      double temp = sqrt(results.MSE*globalXtX_iDiagonal[i]);
      results.SE.push_back(temp);
      results.tStatistics.push_back(results.b[i]/temp);
      results.tStatisticPValues.push_back( 2*tStatCDF(n_q_1, fabs(results.b[i]/temp) ) );
    }
  }
  
  delete Xty;
  delete XtX_i;
  delete btXty;
  delete yty;
  
  results.type = OLSModelType;
  results.success = true;
  
  return true;
}

bool GWAS::computeGLMWithCovariance(Matrix * y, Matrix * X, Matrix * b, GLMResults & results)
{
  if(X->nGlobCols > X->nGlobRows)
  {
    misc.message << "Warning: There are more fixed effects than data points. I am unable to fit the model." << std::endl;
    return false;
  }
  if( this->currentGRMFiltered == NULL )
  {
    misc.error("Error: An internal error was happened. Invalid GRM pointer when fittin the liner model in a GWAS.", 0);
  }
  
  REML reml;

  Matrix * remly = new Matrix(y);
  Matrix * remlX = new Matrix(X);
  Matrix * remlKernel = new Matrix(this->currentGRMFiltered->getNormalizedKernel());
  
  std::vector<Matrix*> kernels;
  kernels.push_back(remlKernel);
  
  std::vector<double> weights;
  weights.push_back(1.);
  
  double initialh2;
  if( this->nTests != 0 )
  {
    initialh2 = this->accumulatedh2/double(this->nTests);
    communicator->broadcast(&initialh2);
  }
  else
  {
    initialh2 = options.initialh2Trait;
  }
  
  bool prepared = reml.prepare(remly, remlX, kernels, initialh2,  weights);
  if( prepared == true )
  {
    reml.computeREML();
  }
  else
  {
    misc.message << "Sorry, a problem was happened while preparing data for performing REML. The MLM cannot be fitted. Please, check the logs." << std::endl;
    return false;
  }
  if(reml.success == false)
  {
    return false;
  }
  
  Matrix * temp = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  temp->multiply(reml.ViX, 'T', reml.y, 'N');
  b->multiply(reml.XtViX_i, 'N', temp, 'N');
  delete temp;
  
  if(this->nTests == 0)
  {
    this->accumulatedh2 = reml.V->variances[0].variance/(reml.V->variances[0].variance + reml.V->variances[1].variance);
  }
  else
  {
    this->accumulatedh2 += ( reml.V->variances[0].variance/(reml.V->variances[0].variance + reml.V->variances[1].variance) );
  }
  this->nTests++;
    
  results = GLMResults();
  
  b->matrixToStandardVector(results.b);
  std::vector<double> globalXtViX_iDiagonal = reml.XtViX_i->diagonal();
  
  if(communicator->mpiRoot)
  {
    for(int i = 0; i < results.b.size(); i++)
    {
      double temp = sqrt(globalXtViX_iDiagonal[i]);
      results.SE.push_back(temp);
      double chi2Test = (results.b[i]*results.b[i])/globalXtViX_iDiagonal[i];
      results.chi2Statistics.push_back(chi2Test);
      results.chi2StatisticsPValues.push_back( chi1_CDF(1, chi2Test) );
    }
  }
  
  results.type = REMLModelType;
  results.success = true;
  
  return true;
}

bool GWAS::computeGroupSignificance(Matrix * y, GLMResults & reducedResults, GLMResults & results)
{
  if( results.type != OLSModelType )
  {
    return true;
  }
  
  if(reducedResults.success == false || results.success == false)
  {
    return false;
  }

  if( misc.gt( (results.btXty - reducedResults.btXty) < 0. ) )
  {
    return false;
  }
  
  if(communicator->mpiRoot)
  {
    double h = double(results.b.size()) - double(reducedResults.b.size());
    double n_q_1 = double(y->nGlobRows) - double(results.b.size());
    results.SSR = results.btXty - reducedResults.btXty;
    results.MSR = results.SSR/double(h);
    results.FStatistic = results.MSR / results.MSE;
    results.FStatisticPValue = FStatCDF(double(h), double(n_q_1), results.FStatistic);
  }
  
  return true;
}

void GWAS::storeResults(std::map<std::string, GLMResults > & effects, Covariate * covariate, std::map<std::string, std::vector<SNP> > & effectsSNPs)
{
  this->significantSNPIds.clear();
  
  if(communicator->mpiRoot)
  {
    Message messageMean(options.outFile + ".gwas.mean");
    Message messageDiscrete(options.outFile + ".gwas.discrete");
    Message messageQuantitative(options.outFile + ".gwas.quantitative");
    Message messageSNPs(options.outFile + ".gwas.snps");
    messageMean << "GROUP NAME BETA SE PV" << std::endl;
    messageDiscrete << "GROUP NAME BETA SE PV" << std::endl;
    messageQuantitative << "GROUP NAME BETA SE PV" << std::endl;
    messageSNPs << "GROUP SNP ALLELE MEAN STDEV BETA NBETA SE PV GROUPPV" << std::endl;
    
    for(std::map<std::string, GLMResults >::iterator it = effects.begin(); it != effects.end(); ++it)
    {
      std::string group = it->first;
      GLMResults groupResults = it->second;
      std::vector<double> globalEffects = groupResults.b;
      std::vector<double> globalPValues;
      if( groupResults.type == OLSModelType )
      {
        globalPValues = groupResults.tStatisticPValues;
      }
      else if( groupResults.type == REMLModelType )
      {
        globalPValues = groupResults.chi2StatisticsPValues;
      }
      else
      {
        misc.error("Error: An internal error was happened. Invalid model type.", 0);
      }
      std::vector<SNP> SNPs = effectsSNPs[group];
      
      if( globalEffects.size() != (covariate->meanNames.size() + covariate->discreteCovarNames.size() + covariate->quantitativeCovarNames.size() + SNPs.size()) )
      {
        misc.error("Error: An internal error was happened. The size of covariate names is not of same dimension of fitted model fixed effects.", 0);
      }
      int shift = 0;
      
      for(int i = 0; i < covariate->meanNames.size(); i++)
      {
        messageMean << group << " " << covariate->meanNames[i];
        messageMean << " " << globalEffects[ i + shift ];
        messageMean << " " << groupResults.SE[ i + shift ];
        messageMean << " " << globalPValues[ i + shift ];
        messageMean << std::endl;
      }
      shift += covariate->meanNames.size();
      
      for(int i = 0; i < covariate->discreteCovarNames.size(); i++)
      {
        messageDiscrete << group << " " << covariate->discreteCovarNames[i];
        messageDiscrete << " " << globalEffects[ i + shift ];
        messageDiscrete << " " << groupResults.SE[ i + shift ];
        messageDiscrete << " " << globalPValues[ i + shift ];
        messageDiscrete << std::endl;
      }
      shift += covariate->discreteCovarNames.size();
      
      for(int i = 0; i < covariate->quantitativeCovarNames.size(); i++)
      {
        messageQuantitative << group << " " << covariate->quantitativeCovarNames[i];
        messageQuantitative << " " << globalEffects[ i + shift ];
        messageQuantitative << " " << groupResults.SE[ i + shift ];
        messageQuantitative << " " << globalPValues[ i + shift ];
        messageQuantitative << std::endl;
      }
      shift += covariate->quantitativeCovarNames.size();
      
      for(int i = 0; i < SNPs.size(); i++)
      {
        messageSNPs << group;
        messageSNPs << " " << SNPs[i].name;
        messageSNPs << " " << SNPs[i].allele2;
        messageSNPs << " " << std::setprecision(3) << 2.*SNPs[i].p2;
        messageSNPs << " " << std::setprecision(3) << SNPs[i].standardDev;
        messageSNPs << " " << globalEffects[ i + shift ];
        messageSNPs << " " << std::setprecision(5) << globalEffects[ i + shift ]/SNPs[i].standardDev;
        messageSNPs << " " << groupResults.SE[ i + shift ];
        messageSNPs << " " << globalPValues[ i + shift ];
        messageSNPs << " " << groupResults.FStatisticPValue;
        messageSNPs << std::endl;
        
        if(globalPValues[ i + shift ] < this->significanceThreshold)
        {
          this->significantSNPIds.push_back(SNPs[i].name);
        }
      }
    }
  }

  if(communicator->mpiRoot)
  {
    if(this->correlatedSNPIds.size() != 0)
    {
      std::vector<std::string> temp(this->correlatedSNPIds.begin(), this->correlatedSNPIds.end());
      std::vector<std::string> intersection = intersectionStringVectors(2, &(this->significantSNPIds), &temp);
      if(intersection.size() != 0)
      {
        std::vector<std::string> uncorrelatedSignifcantSNPIds = differenceBetweenTwoVectors(this->significantSNPIds, intersection);
        misc.message << this->significantSNPIds.size() - uncorrelatedSignifcantSNPIds.size() << " correlated SNPs removed." << std::endl;
        this->significantSNPIds = uncorrelatedSignifcantSNPIds;
        
        Message messageCorrelated(options.outFile + ".gwas.correlatedSNPs");
        for(int i = 0; i < intersection.size(); i++)
        {
          messageCorrelated << intersection[i] << std::endl;
        }
      }
      this->correlatedSNPIds.clear();
    }
  }
}

void GWAS::getLessSignificantCorrelatedSNPs(double threshold, Genotype* genotypes, GLMResults results, int shift)
{
  if(threshold <= 0.)
  {
    return;
  }
  
  Matrix * correlations = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  Matrix * normalization = new Matrix(MATRIX_DEFAULT_DISTRIBUTION);
  correlations->multiply(genotypes->genotypes, 'N', genotypes->genotypes, 'T');
  normalization->multiply(genotypes->missings, 'N', genotypes->missings, 'T');
  correlations->elementWiseDivision(normalization);
  delete normalization;
  
  std::vector<int> idxSNPs1;
  std::vector<int> idxSNPs2;
  correlations->getGlobalIndexOutsideRange(-threshold, threshold, idxSNPs1, idxSNPs2);

  std::vector<std::string> correlatedSNPs;
  if(communicator->mpiRoot == true)
  {
    if( (correlations->nGlobRows + shift) != results.b.size())
    {
      misc.error("Error: An internal error was happened. Discordant number of elements when searching for uncorrelated SNPs.", 0);
    }
    std::set<int> idxsToRemove;
    for(int i = 0; i<idxSNPs1.size(); i++)
    {
      int idxS1 = idxSNPs1[i];
      int idxS2 = idxSNPs2[i];
      
      if( (idxS1 + shift) >= results.tStatisticPValues.size() || (idxS2 + shift) >= results.tStatisticPValues.size() )
      {
        misc.error("Error: An internal error was happened. Incorrect indices.", 0);
      }
      
      if( idxS1 > idxS2 ) //Because the matrix is symmetric.
      {
        if( results.tStatisticPValues[ idxS1 + shift ] < results.tStatisticPValues[ idxS2 + shift ] )
        {
          idxsToRemove.insert(idxS2);
        }
        else
        {
          idxsToRemove.insert(idxS1);
        }
      }
    }
    
    for (std::set<int>::iterator it = idxsToRemove.begin(); it != idxsToRemove.end(); ++it)
    {
      this->correlatedSNPIds.insert( genotypes->SNPIds[*it] );
    }
  }
  
  delete correlations;
}

void GWAS::debugWrite(Matrix * y, Matrix * X, Matrix * b)
{
  if(options.debug)
  {
    y->debugWrite(options.outFile + ".y");
    X->debugWrite(options.outFile + ".X");
    b->debugWrite(options.outFile + ".b");
  }
}