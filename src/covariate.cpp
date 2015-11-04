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

#include "covariate.h"
#include "matrix.h"
#include "misc.h"
#include "global.h"
#include "options.h"
#include "auxiliar.h"
#include "communicator.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <set>

Covariate::Covariate(std::string fcov, std::string fqcov, std::vector<std::string> emptyIndividualIds, bool constructMatrix)
{
  this->covariates = NULL;
  
  this->nIndividuals = 0;
  this->nCovariates = 0;
  this->nQuantitativeCovariates = 0;
  this->nDiscreteCovariates = 0;
  
  readRawCovariate(fqcov, this->rawQCovars);
  readRawCovariate(fcov, this->rawCovars);
  this->rawCleared = false;
  getDiscreteCovariateCategories(this->rawCovars);
  
  if(constructMatrix == true)
  {
    parseRawCovariates(emptyIndividualIds);
  }
}

Covariate::~Covariate()
{
  if(this->covariates!=NULL)
  {
    delete this->covariates;
    this->covariates = NULL;
  }
}

void Covariate::parseRawCovariates(std::vector<std::string> emptyIndividualIds, int nMeans, int idxThisMean)
{
  if(this->rawCleared == true)
  {
    misc.error("Error: An internal error was happened. The covariates can not be parsed when raw vectors are empty.", 0);
  }

  std::vector< std::vector<double> > qCovars;
  std::vector< std::vector<double> > Covars;
  
  reestructureQuantitativeCovariate(this->rawQCovars, qCovars);
  reestructureDiscreteCovariateUsingDifferences(this->rawCovars, Covars);
  
  //Compute the total number of covariates and individuals
  //this->nIndividuals = Covars.size()!=0?Covars.size():qCovars.size();
  this->nCovariates = 0;
  this->nCovariates +=nMeans;
  this->nCovariates += this->nDiscreteCovariates;
  this->nCovariates += this->nQuantitativeCovariates;
  
  this->meanNames.clear();
  for(int i = 0; i < nMeans; i++)
  {
    std::stringstream ss;
    ss << "mean-" << i;
    this->meanNames.push_back(ss.str());
  }
  
  //If no covariates are specified, then, the number of individuals is 0.
  //However the mean matrix must be created. To this end, the individuals in emptyIndividuals are used.
  if(this->nIndividuals == 0)
  {
    if(communicator->mpiRoot == true)
    {
      this->individualIds = emptyIndividualIds;
      for(int i = 0; i < this->individualIds.size(); i++)
      {
        individualIdsIdx[ this->individualIds[i] ] = i;
      }
      this->nIndividuals = this->individualIds.size();
    }
    communicator->broadcast(&this->nIndividuals, 1);
  }
  
  //Create the covariance matrix
  constructCovariateMatrix(Covars, qCovars, nMeans, idxThisMean);
  
  clearRaw();
}

void Covariate::clearRaw()
{
  this->rawCovars.clear();
  this->rawQCovars.clear();
  this->rawCleared = true;
}

void Covariate::readRawCovariate(std::string f, std::vector< std::vector<std::string> > & covars)
{
  std::ifstream file;
  std::string line;
  
  if(f == "")
  {
    covars.clear();
    return;
  }
  
  if(communicator->mpiRoot)
  {
    misc.message << "Reading covariates data from file [ " << f << " ] ..." << std::endl;
    
    misc.checkFileExists(f);
    file.open(f.c_str());
    
    covars.clear();
    int nCovars = -1;
    int idx = 0;
    while(getline(file,line))
    {
      if(!file)
      {
	break;
      }
      
      std::istringstream sstemp(line);
      
      std::string familyId, individualId;
      sstemp >> familyId >> individualId;
      std::string key = familyId + "@" + individualId;
      if(this->nIndividuals != 0) //If a file is readed before, nIndividuals should be != 0 and the individualIds between two files must be the same.
      {
	if(idx>=this->individualIds.size())
	{
	  misc.error("Error: There are not the same individuals in covariate and quantitative covariate files.", 0);
	}
	if(this->individualIds[idx] != familyId + "@" + individualId)
	{
	  misc.error("Error: There are not the same individuals in covariate and quantitative covariate files.", 0);
	}
      }
      else
      {
	
	this->individualIds.push_back(key);
	if(this->individualIdsIdx.count(key) != 0)
	{
	  misc.error("Error: The individual with family Id: " + familyId + " and individual Id: " + individualId + " appears more than one time in file [ " + f + " ].", 0);
	}
	this->individualIdsIdx[key] = idx;
      }

      std::vector<std::string> covarsLine;
      std::string sCovar;
      
      //Read all de covariates for the current individual
      while(sstemp >> sCovar)
      {
	covarsLine.push_back(sCovar);
        if(sCovar == "NA")
        {
          this->individualIdsWithMissingData.insert(key);
        }
      }
      
      //Check that each individual has the same number of covariates
      if(nCovars == -1)
      {
	nCovars = covarsLine.size();
      }
      if(nCovars != covarsLine.size())
      {
	misc.error("Error: The covariables file is not properly formated. The individual " + familyId + " " + individualId + " has less/more covariables than previous individuals.", 0);
      }

      covars.push_back(covarsLine);
      idx++;
    }
    file.close();
    
    if(covars.size()==0)
    {
      misc.error("Error: The covariates file [ " + f + " ] is empty.", 0);
    }
    misc.message << covars[0].size() << " covariates have been read for " << covars.size() << " individuals." << std::endl;
    
    if(this->nIndividuals != 0 && this->nIndividuals != this->individualIds.size())
    {
      misc.error("Error: There are not the same individuals in covariate and quantitative covariate files.", 0);
    }
    this->nIndividuals = this->individualIds.size();
  }
  communicator->broadcast(&this->nIndividuals, 1);
}

void Covariate::getDiscreteCovariateCategories(std::vector< std::vector<std::string> > & rawCovars)
{
  if(this->rawCleared == true)
  {
    misc.error("Error: An internal error was happened. The covariate categories can not be parsed when raw vectors are empty.", 0);
  }
  
  //Get all the posible values of the discrete covars.
  if(communicator->mpiRoot)
  {
    if(rawCovars.size() == 0)
    {
      return;
    }
    
    int nTotalValuesCovar = 0;
    this->covarCategories.clear();
    for(int j=0; j<rawCovars[0].size(); j++)
    {
      std::map<std::string, int> temp;
      int nValuesCurrentCovar = 0;
      for(int i=0; i<rawCovars.size(); i++)
      {
        if(temp.count(rawCovars[i][j]) != 0)
        {
          continue;
        }
        temp[ rawCovars[i][j] ] = nValuesCurrentCovar;
        nValuesCurrentCovar++;
        nTotalValuesCovar++;
      }
      if(double(temp.size()) > double(this->nIndividuals)*0.1)
      {
        misc.error("Error: There are too many categories in the " + i2s( j + 1 ) + " column of covariates file.", 0);
      }
      this->covarCategories.push_back(temp);
    }
    this->nDiscreteCovariates = nTotalValuesCovar;
  }
}

void Covariate::syncronizeDiscreteCovariateCategoriesWith(Covariate * covariate)
{
  if(this->rawCleared == true)
  {
    misc.error("Error: An internal error was happened. Covariates can not be syncronized when raw vectors are empty.", 0);
  }
  
  if(communicator->mpiRoot)
  {
    std::vector< std::map<std::string, int> > combinedCovarCategories;
    if( this->covarCategories.size() != covariate->covarCategories.size() )
    {
      misc.error("Error: Error when combining two covariate files. It seems that their columns are not the same or not represent the same thing.", 0);
    }
    
    int nTotalValuesCovar = 0;
    for(int i = 0; i < this->covarCategories.size(); i++)
    {
      //get the categories in both Covariate objects
      std::set<std::string> categories;
      for(std::map<std::string, int>::iterator it = this->covarCategories[i].begin(); it != this->covarCategories[i].end(); ++it)
      {
        categories.insert(it->first);
      }
      for(std::map<std::string, int>::iterator it = covariate->covarCategories[i].begin(); it != covariate->covarCategories[i].end(); ++it)
      {
        categories.insert(it->first);
      }
      //index the categories
      std::map<std::string, int> indexedCategories;
      int idx = 0;
      for(std::set<std::string>::iterator it = categories.begin(); it != categories.end(); ++it)
      {
        indexedCategories[*it] = idx;
        idx++;
        nTotalValuesCovar++;
      }
      //Append to the new vector of covarCategories;
      combinedCovarCategories.push_back(indexedCategories);
    }
    
    this->covarCategories = combinedCovarCategories;
    covariate->covarCategories = combinedCovarCategories;
    this->nDiscreteCovariates = nTotalValuesCovar;
    covariate->nDiscreteCovariates = nTotalValuesCovar;
  }
}

void Covariate::reestructureQuantitativeCovariate(std::vector< std::vector<std::string> > & rawCovars, std::vector< std::vector<double> > & covars)
{
  covars.clear();
  if(communicator->mpiRoot)
  {
    this->quantitativeCovarNames.clear();
    if(rawCovars.size() != 0)
    {
      this->nQuantitativeCovariates = rawCovars[0].size();
    }
    else
    {
      return;
    }
    
    //Parse data in rawCovars
    std::vector<double> means(rawCovars[0].size(), 0.);
    std::vector<double> nNonMissing(rawCovars[0].size(), 0.);
    std::vector<int> allMissing(rawCovars[0].size(), 1);
    for(int i=0; i<rawCovars.size(); i++)
    {
      std::vector<double> covarsLine;
      for(int j=0; j<rawCovars[i].size(); j++)
      {
        double covar;
        if(rawCovars[i][j] != "NA")
        {
          std::istringstream sstemp2(rawCovars[i][j]);
          if( (sstemp2 >> covar).fail() )
          {
            misc.error("Error: The individual " + this->individualIds[i] + " has not a valid value in one of their covariables: " + rawCovars[i][j], 0);
          }
          covarsLine.push_back(covar);
          means[j] += covar;
          nNonMissing[j] += 1;
          allMissing[j] = 0;
        }
        else
        {
          covarsLine.push_back(0.);
        }
      }
      covars.push_back(covarsLine);
    }
    
    //Check there are no columns with all missings.
    for(int j=0; j<rawCovars[0].size(); j++)
    {
      if(allMissing[j] == 1)
      {
        misc.error("Error: A column in the quantitative covariates file have all values missing.", 0);
      }
    }
    
    //Replace NA values with the mean
    for(int i=0; i<rawCovars.size(); i++)
    {
      for(int j=0; j<rawCovars[i].size(); j++)
      {
        if(rawCovars[i][j] == "NA")
        {
          covars[i][j] = means[j]/nNonMissing[j];
        }
      }
    }
    
    //Store names in order
    for(int j=0; j<rawCovars[0].size(); j++)
    {
      std::stringstream ss;
      ss << "q-col-" << j;
      this->quantitativeCovarNames.push_back(ss.str());
    }
  }
}


void Covariate::reestructureDiscreteCovariate(std::vector< std::vector<std::string> > & rawCovars, std::vector< std::vector<double> > & covars)
{
  covars.clear();
  if(communicator->mpiRoot)
  {
    if(rawCovars.size() == 0)
    {
      return;
    }
   
    //Generate the matrix of dummy variables (each value for each covariate has a column. Each individual has a 1 on the columns that match the their value and category, 0 otherwise)
    for(int i=0; i<rawCovars.size(); i++)
    {
      std::vector<double> temp(this->nDiscreteCovariates, 0);
      int shift = 0;
      for(int j=0; j<rawCovars[0].size(); j++)
      {
	if(this->covarCategories[j].count(rawCovars[i][j]) == 0)
	{
	  misc.error("Error: An internal error was happened.", 0);
	}
	temp[ shift + this->covarCategories[j][ rawCovars[i][j] ] ] = 1.;
	shift += this->covarCategories[j].size();
      }
      covars.push_back(temp);
    }
    
    //Store names in order
    int namesShift = 0;
    this->discreteCovarNames = std::vector<std::string>(this->nDiscreteCovariates, "");
    for(int j=0; j<rawCovars[0].size(); j++)
    {
      for(std::map<std::string, int>::iterator it = this->covarCategories[j].begin(); it != this->covarCategories[j].end(); ++it)
      {
        std::stringstream ss;
        ss << "col-" << j << "(" << it->first << ")";
        this->discreteCovarNames[ namesShift + it->second ] = ss.str();
      }
      namesShift += this->covarCategories[j].size();
    }
  }
}

void Covariate::reestructureDiscreteCovariateUsingDifferences(std::vector< std::vector<std::string> > & rawCovars, std::vector< std::vector<double> > & covars)
{
  covars.clear();
  if(communicator->mpiRoot)
  {
    if(rawCovars.size() == 0)
    {
      return;
    }
    
    this->nDiscreteCovariates = this->nDiscreteCovariates - rawCovars[0].size(); //The first category of each covariate is removed.
    
    //Generate the matrix of dummy variables (each value for each covariate has a column. Each individual has a 1 on the columns that match the their value and category, 0 otherwise)
    for(int i=0; i<rawCovars.size(); i++)
    {
      std::vector<double> temp(this->nDiscreteCovariates, 0);
      int shift = 0;
      for(int j=0; j<rawCovars[0].size(); j++)
      {
        if(this->covarCategories[j][ rawCovars[i][j] ] == 0) //The first category will not be included.
        {
          shift += this->covarCategories[j].size() - 1;
          continue;
        }
        if(this->covarCategories[j].count(rawCovars[i][j]) == 0)
        {
          misc.error("Error: An internal error was happened.", 0);
        }
        temp[ shift + this->covarCategories[j][ rawCovars[i][j] ] - 1 ] = 1.;
        shift += this->covarCategories[j].size() - 1;
      }
      covars.push_back(temp);
    }
    
    //Store names in order
    int namesShift = 0;
    this->discreteCovarNames = std::vector<std::string>(this->nDiscreteCovariates, "");
    for(int j=0; j<rawCovars[0].size(); j++)
    {
      std::string baseCategoryName = "";
      for(std::map<std::string, int>::iterator it = this->covarCategories[j].begin(); it != this->covarCategories[j].end(); ++it)
      {
        if(it->second == 0)
        {
          baseCategoryName = it->first;
        }
      }
      for(std::map<std::string, int>::iterator it = this->covarCategories[j].begin(); it != this->covarCategories[j].end(); ++it)
      {
        if(it->second == 0)
        {
          continue;
        }
        std::stringstream ss;
        ss << "col-" << j << "(" << it->first << "/" << baseCategoryName << ")";
        this->discreteCovarNames[ namesShift + it->second - 1 ] = ss.str();
      }
      namesShift += this->covarCategories[j].size() - 1;
    }
  }
}

void Covariate::constructCovariateMatrix(std::vector< std::vector<double> > & Covars, std::vector< std::vector<double> > & qCovars, int nMeans, int idxThisMean)
{
  double * covariatesTempMatrix;
  
  if (communicator->mpiRoot) {
    if(Covars.size()!=0 && qCovars.size()!=0 && Covars.size() != qCovars.size())
    {
      misc.error("Error: The number of individuals in the covariates file differ from the number of individuals in the quantitative covariate file.", 0);
    }
    
    //Copy Covars and qCovars to covariatesTempMatrix
    covariatesTempMatrix = new double [this->nCovariates*this->nIndividuals];
    
    int columnShift = 0;
    for(int i = 0; i < nMeans; i++)
    {
      double value = ((i == idxThisMean)?1.:0.);
      for(int r = 0; r < this->nIndividuals; r++)
      {
        for(int c = 0; c < 1; c++)
        {
          covariatesTempMatrix[(c+columnShift)*this->nIndividuals + r] = value;
        }
      }
      columnShift++;
    }
    if(Covars.size() != 0)
    {
      for(int r = 0; r < Covars.size(); r++)
      {
	for(int c = 0; c < Covars[r].size(); c++)
	{
	  covariatesTempMatrix[(c+columnShift)*this->nIndividuals + r] = Covars[r][c];
	}
      }
      columnShift += Covars[0].size();
    }
    if(qCovars.size() != 0)
    {
      for(int r = 0; r < qCovars.size(); r++)
      {
	for(int c = 0; c < qCovars[r].size(); c++)
	{
	  covariatesTempMatrix[(c+columnShift)*this->nIndividuals + r] = qCovars[r][c];
	}
      }
      columnShift += qCovars[0].size();
    }
  }
  
  communicator->broadcast(&this->nIndividuals, 1);
  communicator->broadcast(&this->nCovariates, 1);
  
  //Distribute covariatesTempMatrix between nodes.
  if(this->covariates!=NULL)
  {
    delete this->covariates;
    this->covariates = NULL;
  }
  
  this->covariates = new Matrix(MATRIX_DEFAULT_DISTRIBUTION, this->nIndividuals, this->nCovariates, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  this->covariates->scatterMatrix(covariatesTempMatrix);
  
  if (communicator->mpiRoot) {
    delete [] covariatesTempMatrix;
  }
}

void Covariate::filterIndividuals(std::vector<std::string> keepIndividualIds)
{
  double * originalGlobMatrix;
  double * filteredGlobMatrix;
  if(communicator->mpiRoot)
  {
    originalGlobMatrix = new double [this->covariates->nGlobRows*this->covariates->nGlobCols];
    filteredGlobMatrix = new double [keepIndividualIds.size()*this->nCovariates];
  }
  
  this->covariates->gatherMatrix(originalGlobMatrix);
  
  if(communicator->mpiRoot)
  {
    std::map<std::string, int> newIndividualIdsIdx;
    std::set<std::string> newIndividualIdsWithMissingData;
    
    for(int r=0; r<keepIndividualIds.size(); r++) //Is it better access rows contiguously and searxh in map each step? Or inverse?
    {
      if(this->individualIdsIdx.count(keepIndividualIds[r]) == 0)
      {
	misc.error("Error: An internal error was happened. The individual '" + keepIndividualIds[r]  + "' is not in the covariance matrix.", 0);
      }
      if(newIndividualIdsIdx.count(keepIndividualIds[r]) != 0)
      {
	misc.error("Error: An internal error was happened. The individual '" + keepIndividualIds[r]  + "' is repeated.", 0);
      }
      if(this->individualIdsWithMissingData.find(keepIndividualIds[r]) != this->individualIdsWithMissingData.end())
      {
        newIndividualIdsWithMissingData.insert(keepIndividualIds[r]);
      }
      
      int oldRowIdx = this->individualIdsIdx[keepIndividualIds[r]];
      for(int c=0; c<this->nCovariates; c++)
      {
	filteredGlobMatrix[c*keepIndividualIds.size() + r] = originalGlobMatrix[c*this->nIndividuals + oldRowIdx];
      }
      newIndividualIdsIdx[keepIndividualIds[r]] = r;
    }
    this->nIndividuals = keepIndividualIds.size();
    this->individualIds = keepIndividualIds;
    this->individualIdsIdx = newIndividualIdsIdx;
    
    this->individualIdsWithMissingData = newIndividualIdsWithMissingData;
  }
  communicator->broadcast(&this->nIndividuals, 1);
  communicator->barrier();
  
  this->covariates->initParameters(this->nIndividuals, this->nCovariates, communicator->nDefaultBlockRows, communicator->nDefaultBlockCols);
  this->covariates->scatterMatrix(filteredGlobMatrix);
  
  if(communicator->mpiRoot)
  {
    delete [] originalGlobMatrix;
    delete [] filteredGlobMatrix;
  }
}

void Covariate::printCovariate(int nSetw)
{
  double *cov;
  
  if (communicator->mpiRoot) {
    cov = new double [this->nIndividuals*this->nCovariates];
  }
  
  this->covariates->gatherMatrix(cov);
  if (communicator->mpiRoot) {
    misc.message << "Covariates Matrix:\n";
    misc.message << this->individualIds.size() << " " << this->individualIdsIdx.size() << std::endl ;
    for (int r = 0; r < this->nIndividuals; r++)
    {
      misc.message << this->individualIds[r] << " " << this->individualIdsIdx[this->individualIds[r]] << ": " ;
      for (int c = 0; c < this->nCovariates; c++)
      {
	misc.message << cov[c*this->nIndividuals + r] << " ";
      }
      misc.message << std::endl;
    }
    misc.message << std::endl;
    
    for (int r = 0; r < this->nIndividuals; r++)
    {
      misc.message << this->individualIds[r] << " " << this->individualIdsIdx[this->individualIds[r]] << ": mu " ;
      int shift = 0;
      for (int c = 0; c < this->meanNames.size(); c++)
      {
        misc.message << std::setw(nSetw) << cov[(c + shift)*this->nIndividuals + r] << " ";
      }
      shift += this->meanNames.size();
      for (int c = 0; c < this->nDiscreteCovariates; c++)
      {
        if(cov[(c + shift)*this->nIndividuals + r] == 1)
        {
          misc.message << std::setw(nSetw) << discreteCovarNames[c];
        }
        else
        {
          misc.message << std::setw(nSetw) << "-";
        }
      }
      shift += this->nDiscreteCovariates;
      for (int c = 0; c < this->nQuantitativeCovariates; c++)
      {
        misc.message << std::setw(nSetw) << cov[(c + shift)*this->nIndividuals + r] << " ";
      }
      misc.message << std::endl;
    }
    misc.message << std::endl;
    
    misc.message.flush();
    
    misc.message << "With missings:" << std::endl;
    for(std::set<std::string>::iterator it = this->individualIdsWithMissingData.begin(); it != this->individualIdsWithMissingData.end(); ++it)
    {
      misc.message << *it << std::endl;
    }
  }
  
  if (communicator->mpiRoot) {
    delete [] cov;
  }
}
