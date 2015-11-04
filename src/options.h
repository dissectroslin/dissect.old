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

#ifndef OPTIONS_H
#define OPTIONS_H

#include "range.h"
#include "genotype.h"

#include <string>
#include <sstream>
#include <map>

enum AnalysisToPerform
{
  withoutAnalysis,                         ///<No analysis specified
  justCheckAnalysis,                       ///<Just check if the options can be parsed properly
  showHelpAnalysis,                        ///<Show help?
  makeGRMAnalysis,                         ///<Compute the GRM?
  REMLAnalysis,                            ///<Perform REML analysis?
  GWASAnalysis,                            ///<Perform GWAS analysis?
  recursiveGWASAnalysis,                   ///<Perform a recurive GWAS analysis.
  bivarREMLAnalysis,                       ///<Perform bivariate REML analysis?
  multivarREMLAnalysis,                    ///<Perform multi REML analysis?
  simulatePhenotypeAnalysis,               ///<SimulatePhenotypes?
  predictPhenotypeAnalysis,                ///<PredictPhenotypes?
  PCAAnalysis                              ///<Perform a PCA?
};

class Options
{
public:
  ////////////////////////////
  // Actions
  AnalysisToPerform analysis;
  
  ////////////////////////////
  // Base input files
  std::string genotypeFile;             ///<Genotype file
  std::string genotypeListFile;         ///<File with a list of genotype files
  std::string grmFile;                  ///<GRM file
  std::string phenotypesFile;           ///<Phenotypes file
  std::string covarsFile;               ///<Covariates file
  std::vector<std::string> covarsFiles; ///<Covariates files list (for bivariate multi REML analysis)
  std::string qCovarsFile;              ///<Quantitative covariates file
  std::vector<std::string> qCovarsFiles;///<Quantitative covariates files list (for bivariate multi REML analysis)
  
  ////////////////////////////
  // Output files
  std::string outFile;                  ///<Output file
  
  ////////////////////////////
  // GRM
  int GRMJoinMethod;                    ///<How GRM will be constructed from multiple genotype files. 0: computing grm from each file and then add all together. 1: join all genotypes, then compute grm.
  int useMPIForWriting;                 ///<GRM writing will be serialized if this is false. It could be useful for some file systems.
  bool diagonalizeGRM;                  ///<Diagonalize computed GRMs?
  bool writeAlsoNormalGRM;              ///<If diagonalizeGRMs == true, store also undiagonalized GRM
  double minimumNumberOverlapingSNPs;   ///<The minimum number of overlapping SNPs for considering correct the correlation between two individuals.
  double maximumProportionOfElimitaedIndividuals; ///<The maximum proportion of individuals allowed to be filtered due to not enouth SNP overlapping before disregarding a GRM.
  
  ////////////////////////////
  // REML
  bool computeBLUE;                     ///<Compute the best linear unbiased estimator?
  bool computeIndividualsBLUP;          ///<Compute the best linear unbiased predictor for the total genetic effect?
  bool computeSNPsBLUP;                 ///<Compute the best linear unbiased predictor for each individual SNP?
  bool fixCorrelation;                  ///<Repeat REML with a fixed correlation?
  double fixedCorrelation;              ///<Value of the fixed correlation.
  int phenotypeColumn;                  ///<Column of the phenotype used in RMEL in file phenotypesFile.
  std::vector<int> phenotypeColumns;    ///<Columns of the phenotypes in file phenotypesFile to analyze i nbiver/multi REML.
  bool analyzeAllPhenos;                ///<Iteratively analyze all phenotypes.
  int REMLMethod;                       ///<REML method 0->AI, 1->EM.
  bool environmentalCovariance;         ///<Environmental covariance must be computed on bivariate REML?
  bool computeREMLInSubsample;          ///<Compute REML in a subsample before REML with all samples.
  int nSubSampleIterations;             ///<Number of REML subsample computations.
  double initialSubsampleFraction;      ///<Proportion of individuals that will be used as a first REML approximation.
  int minimumSubsample;                 ///<Minimum number of individuals on random subsample.
  double initialh2Trait;                ///<The initial estimated heretability used for estimating initial variance values for trait in REML
  std::vector<double> initialh2Traits;  ///<The initial estimated heretabilities used for estimating initial variance values for traits when computing bivar/multi REML;
  double varianceConvergenceThreshold;  ///<The threshold for variance convergence.
  double changeAIStepThreshold;         ///<Threshold for switching from AI to short AI or EM
  double allowSwitchFromAItoEM;         ///<When loglikelihood relative difference > changeAIStepThreshold, a switch from AI to EM is performed when this is true. (needs testing)
  std::string initialVariancesFile;     ///<Define a file with initial variances in.
  bool correctLinkageDisequilibrium;    ///<Correct for linkage disequillibrium in single REML?
  bool joinCovariatesVertically;        ///<If true, on bivar REML analysis, joinCovariates vertically.
  bool forceUseDiagonalizedKernels;     ///<Force the use of diagonalized GRMs for analysis that do not allow them. GRMs will be undiagonalized before use.
  bool computeEpistasisVariance;        ///<Compute the variance due to epistatic effects?
  double stepWeightingConstant;         ///<Constant used for weighting AI REML step when previous change in likelihood is large.
  
  ////////////////////////////
  // GWAS
  bool gwasWithMixedModel;              ///<Perform GWAS using a mixed model correcting for population structure.
  double correlatedSNPsThreshold;       ///<Threshold above which two SNPs are considered to be correlated enough for filtering one of two.
  std::string genotypeAndGRMsListFile;  ///<File with a list of genotype files and GRM files used for GWAS correcting for population structure.
  
  ////////////////////////////
  // Recursive GWAS
  double relationFitSNPsIndividuals;      ///<When performing a grouped recursive GWAS, the ratio of the number of SNPs to fit related to the number of samples.
  double significanceThresholdFilterSNPs; ///<Significance threshold for filtering SNPs in a recursive GWAS.
  int recursiveGWASMaxIterations;         ///<The maximum number of iterations for recursive GWAS. If negative, there are no maximum.
  
  ////////////////////////////
  // Limits
  int maxREMLIterations;                ///<Max number of REML iterations
  double varianceConstrainProportion;   ///<If a variance is negative in a REML step. It will be constrained to a base variance multiplied by this variable.
  long randomSeed;                      ///<Seed used for random number generation.
  int logOutputPrecision;
  int resultsOutputPrecision;
  int logFieldWidth;
  int resultsFieldWidth;
  double grmCutoff;                     ///<Threshold over which individuals will be pruned in the grm
  bool pruneGRM;                        ///<Prune GRMs?
  
  ////////////////////////////
  // Phenotype simulation vars
  std::string effectsSizeFile;          ///<file of SNP effects
  std::string adjustEffectsGenotypeListFile; ///<genotype file for adjusting SNP effects using their frequencies.
  double simulationh2;                  ///<h square for the simulation.
  double prevalence;                    ///<Prevalence of the binary trait
  bool simulateBinaryTrait;             ///<The trait is binary?
  bool simulateQuantitativeTrait;       ///<The trait is quantitative?
  
  ////////////////////////////
  // Phenotype prediction vars
  std::string snpEffectsFile;           ///<file with SNP effects
  
  ////////////////////////////
  // PCA
  int nEigenValues;                     ///<The number of eigenvalues will be stored.
  
  ////////////////////////////
  // Regional analysis
  bool regionalAnalysis;                ///<Perform regional analysis?
  GroupBy regionBy;                     ///<How regions will be defined?
  std::string regionsFile;              ///<SNPs are grouped by groups defined in this file. This could be used in REML regional analysis or GWAS.
  int regionSize;                       ///<In a regional analysis, the regions (defined by position) are of this size.
  int regionOverlapping;                ///<In a regional analysis, the regions (defined by position) have this overlapping.
  int minSNPsInRegions;                 ///<The minimum SNPs ina region for computing the region.
  int fixedGroupSize;                   ///<The size of the group when SNPs grouped using ordered fixed sized groups.
  
  ////////////////////////////
  // Divers
  bool removeMissings;                  ///<If true, individuals with a missing value on any fixed effect will be removed.
  bool remlGCTAMode;
  int defaultBlockSize;                 ///<Set the default matrix distribution block sizes.
  bool allowSinglePrecisionInversion;   ///<Allow single precision on matrix inversions.
  std::string fileSNPsToKeep;           ///<When loading a genotype file, keep only SNPs in this file.
  bool allowFixingVariancesToZero;      ///<Allow variances to be fixed to zero when performing a REML.
  
  ////////////////////////////
  // Debug options
  bool mpiDebug;
  bool debug;
  
  ////////////////////////////
  // Parsed options
  std::map<std::string, std::string> parsedOptions;     ///<Parsed options
  
  Options();
  Options(int argc, char **argv);
  ~Options();
  
  void showOptions();
  
  void defaultOptions();
  void parseOptions(int argc, char **argv);
  void checkIncompatibilities();
  void fillMissingOptions();
  
  void setAnalysis(AnalysisToPerform analysisToSet);
  std::string getString(int argc, char **argv, int i);
  std::string getFileName(int argc, char **argv, int i);
  int getInt(int argc, char **argv, int i, Range range = Range());
  double getDouble(int argc, char **argv, int i, Range range = Range());
  std::vector<std::string> getStringList(int argc, char **argv, int i);
  std::vector<std::string> getFileNameList(int argc, char **argv, int i);
  std::vector<int> getIntList(int argc, char **argv, int i, Range range);
  std::vector<double> getDoubleList(int argc, char **argv, int i, Range range);
  
  void appendParsedOptions(std::string option, std::string parameters);
  void showParsedOptions();
  
  /**
   * Check if passed GRMs are Diagonalize
   * 
   * \param fn file name of the grm
   * \return true if grm is diagonal, false otherwise.
   */
  bool isGRMDiagonal(std::string fn);
};

#endif
