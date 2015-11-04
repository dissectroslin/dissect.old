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

#include "analysis.h"
#include "singlereml.h"
#include "multireml.h"

#include <cstdlib>

Analysis::Analysis()
{
  unsigned int seed = options.randomSeed;
  std::srand ( seed );
}

Analysis::~Analysis()
{
}

void Analysis::makeGRM()
{
  if(options.regionalAnalysis == false)
  {
    Kernel *grm = loadGRMUsingOptions();
    bool sanitized = grm->sanitizeKernel();
    if( sanitized == true )
    {
      if( options.diagonalizeGRM == true )
      {
        if( options.writeAlsoNormalGRM == true )
        {
          grm->writeKernel(options.outFile + ".nondiagonal");
        }
        grm->diagonalizeKernel();
      }
      grm->writeKernel(options.outFile);
    }
    else
    {
      misc.message << "Sorry, a problem was happened while computing the GRM. Please, check the logs." << std::endl;
    }
    delete grm;
  }
  else
  {
    Genotype * genotype = loadGenotypeUsingOptions();
    genotype->groupSNPs(options.regionBy);
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
      Genotype *regionalGenotype = new Genotype();
      genotype->genotypeOfSNPsGroup(group, regionalGenotype);
      Kernel *grm = new Kernel(regionalGenotype);
      bool sanitized = grm->sanitizeKernel();
      if( sanitized == true )
      {
        if( options.diagonalizeGRM == true )
        {
          if( options.writeAlsoNormalGRM == true )
          {
            grm->writeKernel(options.outFile + ".nondiagonal");
          }
          grm->diagonalizeKernel();
        }
        grm->writeKernel(options.outFile);
      }
      else
      {
        misc.message << "WARNING: Sorry, a problem was happened while computing the GRM on region " << group << ". Please, check the logs." << std::endl;
      }
      delete grm;
      delete regionalGenotype;
    }
    options.outFile = baseOutFile;
  }
}

void Analysis::makeREML()
{
  if(options.regionalAnalysis == false)
  {
    SingleREML singleREML;
    singleREML.compute();
  }
  else
  {
    SingleREML singleREML;
    singleREML.computeRegional();
  }
}

void Analysis::makeMultivarREML()
{
  if(options.regionalAnalysis == false) //Regional analysis?
  {
    MultiREML multiREML;
    multiREML.compute();
  }
  else
  {
    MultiREML multiREML;
    multiREML.computeRegional();
  }
}

void Analysis::makeSimulatePhenotype()
{
  Genotype * genotype = loadGenotypeUsingOptions();
  Genotype * adjustEffectsGenotype = NULL;
  if( options.adjustEffectsGenotypeListFile != "" )
  {
    adjustEffectsGenotype = new Genotype();
    adjustEffectsGenotype->loadList(options.adjustEffectsGenotypeListFile);
  }
  SimulatePhenotype simulatePhenotype(genotype, options.effectsSizeFile, adjustEffectsGenotype);
  simulatePhenotype.simulatePhenotypes();
}

void Analysis::makePredictPhenotype()
{
  Genotype *genotype;
  if(options.genotypeFile != "") //Read single genotype file
  {
    genotype = new Genotype(options.genotypeFile);
    PredictPhenotype predictPhenotype(genotype, options.snpEffectsFile);
    predictPhenotype.predictPhenotypes();
    predictPhenotype.storePredictions();
  }
  else if (options.genotypeListFile != "") //Read multiple genotype files
  {
    misc.message << "Predicting phenotypes from genotypes specified in [ " << options.genotypeListFile << " ]..." << std::endl;
    misc.message.tab = "  ";
    misc.setGetElapsedTime("AddingEffects");
    std::vector<std::string> fileList;
    getListFromFile(options.genotypeListFile, fileList);

    genotype = new Genotype(fileList[0]);
    PredictPhenotype predictPhenotype(genotype, options.snpEffectsFile, false);
    predictPhenotype.predictPhenotypes();

    for(int i = 1; i < fileList.size(); i++)
    {
      genotype = new Genotype(fileList[i]);
      PredictPhenotype phenotypeToAdd(genotype, options.snpEffectsFile, false);
      phenotypeToAdd.predictPhenotypes();
      predictPhenotype.addMoreEffects(&phenotypeToAdd);
    }
    misc.message.tab = "";
    misc.message << "Genotype effects added successfully after " << misc.setGetElapsedTime("AddingEffects", true) << ". Storing results..." << std::endl;
    predictPhenotype.storePredictions();
  }
  else
  {
    misc.error("Error: No genotype file(s) specified.", 0);
  }
}

void Analysis::makePCA()
{
  Kernel *grm = loadGRMUsingOptions();
  bool sanitized = grm->sanitizeKernel();
  if( sanitized != true )
  {
    misc.message << "WARNING: Sorry, a problem was happened while computing GRM for performing a PCA. Please, check the logs." << std::endl;
  }
  PCA pca(grm);
  delete grm;
}

void Analysis::makeGWAS()
{
  GWAS gwas;
  gwas.computeGWAS();
}

void Analysis::makeRecursiveGWAS()
{
  GWAS gwas;
  gwas.computeGWAS();
}


