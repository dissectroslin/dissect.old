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

#include "auxiliar.h"
#include "global.h"
#include "misc.h"
#include "options.h"
#include "matrix.h"
#include "kernel.h"
#include "genotype.h"
#include "phenotype.h"

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <cstdarg>
#include <cmath>

#ifdef BOOSTLIB
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#endif

std::vector<std::string> intersectionStringVectors(int count, ...)
{
  std::vector<std::string> vResult;
  
  va_list vectors;
  va_start(vectors, count);
  
  vResult = *va_arg(vectors, std::vector<std::string>*);
  std::sort(vResult.begin(), vResult.end());
  for(int i = 1; i < count; i++)
  {
    std::vector<std::string> vTempResult;
    std::vector<std::string> vAdd = *va_arg(vectors, std::vector<std::string>*);
    std::sort(vAdd.begin(), vAdd.end());
    
    std::set_intersection(vResult.begin(), vResult.end(), vAdd.begin(), vAdd.end(), std::back_inserter(vTempResult));
    vResult = vTempResult;
  }
  
  return vResult;
}

std::vector<std::string> orderVectorAsTemplate(std::vector<std::string> & templateVector, std::vector<std::string> & sourceVector)
{
  std::set<std::string> sourceSet(sourceVector.begin(), sourceVector.end());
  std::vector<std::string> resultVector;
  for(int i = 0; i<templateVector.size(); i++)
  {
    if(sourceSet.find(templateVector[i]) != sourceSet.end())
    {
      resultVector.push_back(templateVector[i]);
    }
  }
  
  return resultVector;
}

std::vector<std::string> differenceBetweenTwoVectors(std::vector<std::string> & baseVector, std::vector<std::string> & vectorToSubstract)
{
  std::vector<std::string> test = orderVectorAsTemplate(vectorToSubstract, baseVector);
  if(test != vectorToSubstract)
  {
    misc.error("Error: Invalid operation. Two vectors cannot be substracted if all the elements of one are not contained into the other.", 0);
  }
  
  std::set<std::string> substractSet(vectorToSubstract.begin(), vectorToSubstract.end());
  std::set<std::string> baseSet(baseVector.begin(), baseVector.end());
  if(substractSet.size() != vectorToSubstract.size() || baseSet.size() != baseVector.size() )
  {
    misc.error("Error: Invalid operation. Two vectors cannot be substracted if one of them have repeated elements.", 0);
  }
  std::vector<std::string> resultVector;
  for(int i = 0; i<baseVector.size(); i++)
  {
    if(substractSet.find(baseVector[i]) == substractSet.end())
    {
      resultVector.push_back(baseVector[i]);
    }
  }
  
  return resultVector;
}

std::vector<std::string> getRandomSample(std::vector<std::string> & sourceVector, int nElements)
{
  if(nElements >= sourceVector.size() || nElements < 1)
  {
    misc.error("Error: An internal error was happened. Wrong number of elements for a random sample.", 0);
  }
  
  std::vector<std::string> shuffledVector = sourceVector;
  
  std::random_shuffle ( shuffledVector.begin(), shuffledVector.end() );
  
  std::vector<std::string> result(shuffledVector.begin(), shuffledVector.begin() + nElements);
  return orderVectorAsTemplate(sourceVector, result);
}

std::vector<int> extractMapValues(std::vector<std::string> & keysVector, std::map<std::string, int> & valuesMap)
{
  std::vector<int> result;
  
  for(int i = 0; i < keysVector.size(); i++)
  {
    if(valuesMap.count(keysVector[i]) == 0)
    {
      misc.error("Error: An internal error was happened. The key is not in the map.", 0);
    }
    result.push_back( valuesMap[ keysVector[i] ] );
  }
  
  return result;
}

void getListFromFile(std::string f, std::vector<std::string> & list)
{
  list.clear();
  std::ifstream file;
  misc.checkFileExists(f);
  file.open(f.c_str());
  
  std::string line;
  while(getline(file,line))
  {
    std::istringstream sstemp(line);
    std::string element;
    sstemp >> element;
    if(element != "")
    {
      list.push_back(element);
    }
  }

  file.close();
  
  if(list.size() == 0)
  {
    misc.error("Error: File [ " + f + " ] is empty.", 0);
  }
}

void getListFromFile(std::string f, std::vector<double> & list)
{
  std::vector<std::string> temp;
  
  getListFromFile(f, temp);
  
  list.clear();
  
  for(int i = 0; i<temp.size(); i++)
  {
    double value;
    std::stringstream sstemp(temp[i]);
    if( (sstemp >> value).fail() )
    {
      misc.error("Error: The value " + temp[i] + " in file [ " + f + " ] is not a valid number.", 0);
    }
    list.push_back(value);
  }
}

void getTableFromFile(std::string f, std::vector< std::vector<std::string> > & table, int nColumns)
{
  table.clear();
  std::ifstream file;
  misc.checkFileExists(f);
  file.open(f.c_str());
  
  std::string line;
  while(getline(file,line))
  {
    std::istringstream sstemp(line);
    std::string element;
    std::vector<std::string> row;
    int nElements = 0;
    while((sstemp >> element) && nElements < nColumns)
    {
      if(element != "")
      {
        row.push_back(element);
        nElements++;
      }
    }
    if(row.size() ==  0)
    {
      continue;
    }
    if(row.size() != nColumns)
    {
      misc.error(" Error: Line: \n " + line + "\n in file [ " + f + " ] does not have at least " + i2s(nColumns) + " elements.", 0);
    }
    table.push_back(row);
  }

  file.close();
  
  if(table.size() == 0)
  {
    misc.error("Error: File [ " + f + " ] is empty.", 0);
  }
}

int getNumberOfFileColumns(std::string f)
{
  std::ifstream file;
  std::string line;

  misc.checkFileExists(f);
  file.open(f.c_str());
  
  getline(file,line);
  if(!file)
  {
    misc.error("Error: File [ " + f + " ] is empty.", 0);
  }
  
  std::istringstream sstemp(line);
  std::string dummy;
  int nColumns = 0;
  while(sstemp >> dummy)
  { 
    nColumns++;
  }
    
  if(nColumns < 1)
  {
    misc.error("Error: Phenotypes file [ " + f + " ] is not properly formated.", 0);
  }
  
  file.close();

  return nColumns;
}

int leastCommonMultiple(int a, int b)
{
  if(a < 0 || b < 0)
  {
    misc.error("Error: An internal error was happened while computing the Least Common Multiple", 0);
  }
  return (a*b)/greatestCommonDivisor(a,b);
}

int greatestCommonDivisor(int a, int b)
{
  if(a < 0 || b < 0)
  {
    misc.error("Error: An internal error was happened while computing the Greatest Common Divisor", 0);
  }
  if(b == 0)
  {
    return a;
  }
  else
  {
    return greatestCommonDivisor(b, a % b);
  }
}

double computeVariance(Matrix * m)
{
  int dimension;
  double *globalM;
  double variance;
  
  if(m->nGlobRows != 1 && m->nGlobCols != 1)
  {
    misc.error("Error: Unable to compute the variance of a matrix with more than 1 rows and columns.", 0);
  }
  
  if( (m->nGlobRows == 1 && m->nGlobCols <= 1) || m->nGlobRows <= 1 && m->nGlobCols == 1) 
  {
    misc.error("Error: Unable to compute the phenotypic variance of less than two elements.", 0);
  }
  
  if(m->nGlobRows > m->nGlobCols)
  {
    dimension = m->nGlobRows;
  }
  else
  {
    dimension = m->nGlobCols;
  }
  
  if (communicator->mpiRoot) {
    globalM = new double [dimension];
  }
  
  m->gatherMatrix(globalM);
  
  if (communicator->mpiRoot)
  {
    double avg = 0.;
    for(int r = 0; r<dimension; r++)
    {
      avg += globalM[r];
    }
    avg /= double(dimension);
    variance = 0.;
    for(int r = 0; r<dimension; r++)
    {
      double temp = globalM[r] - avg;
      variance += temp * temp;
    }
    variance /= (double(dimension) - 1.);
  }
  
  communicator->broadcast(&variance, 1);
  
  if (communicator->mpiRoot) {
    delete [] globalM;
  }
  
  return variance;
}

std::string getString(int i)
{
  std::stringstream ss;
  ss << i;
  return ss.str();
}

std::string i2s(int i)
{
  return getString(i);
}

std::string spacetab2underscore(std::string s)
{
  std::string text = s;
  std::replace(text.begin(), text.end(), ' ', '_');
  std::replace(text.begin(), text.end(), '\t', '_');
  return text;
}

Genotype * loadGenotypeUsingOptions()
{
  Genotype *genotype;
  if(options.genotypeFile != "") //Read single genotype file
  {
    genotype = new Genotype(options.genotypeFile);
  }
  else if (options.genotypeListFile != "") //Read multiple genotype files
  {
    genotype = new Genotype();
    genotype->loadList(options.genotypeListFile);
  }
  else
  {
    misc.error("Error: No genotype file(s) specified.", 0);
  }
  return genotype;
}

Kernel * loadGRMUsingOptions(bool returnGenotype, Genotype **returnedGenotype)
{
  Kernel *grm = NULL;
  
  if(returnGenotype == true)
  {
    if(returnedGenotype == NULL)
    {
      misc.error("Error: An internal error was happened when loading the GRM. Wrong pointer for returning a genotype.", 0);
    }
    *returnedGenotype = NULL;
  }
  
  if( returnGenotype == true && options.GRMJoinMethod == 0 )
  {
    misc.message << "Changing GRM joining method to method 1." << std::endl;
    options.GRMJoinMethod = 1;
  }
  
  if(options.grmFile == "")
  {
    if(options.genotypeListFile == "") //Create GRM from single genotype file
    {
      Genotype *genotype = new Genotype(options.genotypeFile);
      grm = new Kernel(genotype);
      if(returnGenotype == false)
      {
        delete genotype;
      }
      else
      {
        *returnedGenotype = genotype;
      }
    }
    else //Create GRM from multiple genotype files
    {
      if(options.GRMJoinMethod == 0) //Create GRM using method 0: First compute GRMs for each genotype, then add GRMs.
      {
        misc.message << "Computing GRM using a list of genotype files. Method 0: First compute GRMs for each genotype, then add GRMs." << std::endl;
        misc.message.tab = "  ";
        misc.setGetElapsedTime("GRMPartComputation");
        std::ifstream file;
        misc.checkFileExists(options.genotypeListFile);
        file.open(options.genotypeListFile.c_str());
        
        int idx = 0;
        std::string line;
        std::vector<std::string> loadedSNPs;
        while(getline(file,line))
        {
          std::istringstream sstemp(line);
          std::string partialFileName;
          sstemp >> partialFileName;
          
          if( idx == 0 )
          {
            Genotype *genotype = new Genotype(line);
            loadedSNPs = genotype->SNPIds;
            
            grm = new Kernel(genotype, false);
            if(returnGenotype == false)
            {
              delete genotype;
            }
            else
            {
              misc.error("Error: An internal error was happened. GRM cannot be loaded using method 0 when genotypes must also be loaded.", 0);
            }
          }
          else
          {
            Genotype *genotype = new Genotype(line);
            std::vector<std::string> SNPsIntersection = intersectionStringVectors(2, &loadedSNPs, &genotype->SNPIds); //Check there are not intersections between SNPs in files.
            if(SNPsIntersection.size() != 0)
            {
              misc.error("Error: When computing a GRM from multiple genotype files. There are SNPs repeated in different files.", 0);
            }
            loadedSNPs.insert( loadedSNPs.end(), genotype->SNPIds.begin(), genotype->SNPIds.end() );
            
            Kernel *grmToAdd = new Kernel(genotype, false);
            delete genotype;
            misc.message << "Adding GRM..." << std::endl;
            grm->addKernels(1., grmToAdd);
            delete grmToAdd;
          }
          
          idx++;
        }
        file.close();
        misc.message.tab = "";
        grm->normalize();
        misc.message << "GRM computation finished after " << misc.setGetElapsedTime("GRMPartComputation", true) << std::endl;
      }
      else  //Create GRM using method 1: First join genotypes, then compute GRM.
      {
        misc.message << "Computing GRM using a list of genotype files. Method 1: First join genotypes, then compute GRM." << std::endl;
        misc.message.tab = "  ";
        misc.setGetElapsedTime("GRMPartComputation");
        Genotype *genotype = new Genotype();
        genotype->loadList(options.genotypeListFile);
        grm = new Kernel(genotype);
        misc.message.tab = "";
        misc.message << "GRM computation finished after " << misc.setGetElapsedTime("GRMPartComputation", true) << std::endl;
        if(returnGenotype == false)
        {
          delete genotype;
        }
        else
        {
          *returnedGenotype = genotype;
        }
      }
    }
  }
  else
  {
    grm = new Kernel(options.grmFile);
  }
  return grm;
}

std::vector<int> getPhenotyesForAnalysis()
{
  std::vector<int> phenotypesForAnalyze;
  if(options.analyzeAllPhenos == false)
  {
    phenotypesForAnalyze.push_back(options.phenotypeColumn);
  }
  else
  {
    int nMaxPhenoCol = getNumberOfFileColumns(options.phenotypesFile) - 2;
    for(int i = 0; i<nMaxPhenoCol; i++)
    {
      phenotypesForAnalyze.push_back(i + 1);
    }
  }
  return phenotypesForAnalyze;
}

float ran3(long *idnum)
{
  static int inext, inextp;
  static long ma[56];   /* The value 56 (range ma[1..55]) is special */
  /* and should not be modified; see Knuth.    */
  static int iff=0;
  long mj, mk;
  int i, ii, k;
  
  if (*idnum<0 || iff==0) {    /* Initialization */
    iff = 1;
    /* Initialize ma[55] using the seed idnum and the large number MSEED */
    mj = MSEED - (*idnum<0 ? -*idnum : *idnum);
    mj %= MBIG;
    ma[55] = mj;
    mk = 1;
    /* Now initizalize the rest of the table, in a slightly       */
    /* random order, with numbers that are not especially random. */
    for (i=1; i<=54; i++) {
      ii = (21*i) % 55;
      ma[ii] = mk;
      mk = mj - mk;
      if (mk < MZ) mk += MBIG;
      mj = ma[ii];
    }
    /* We randomize them by "warming up the generator." */
    for (k=1; k<=4; k++) 
      for (i=1; i<=55; i++) {
        ma[i] -= ma[1+(i+30) % 55];
        if (ma[i]<MZ) ma[i] += MBIG;
      }
      inext = 0;     /* Prepare indices for our first generated number. */
      inextp = 31;   /* The constant 31 is special; see Knuth */
      *idnum = 1;
  }
  /* Here is where we start, except on initialization */
  if (++inext==56) inext = 1;   /* Initizalize inext and inextp, wrapping    */
    if (++inextp==56) inextp = 1; /* arround 56 to 1.                          */
      mj = ma[inext] - ma[inextp];  /* Generate new random number substractively */
      if (mj<MZ) mj +=MBIG;         /* Make sure that it is in range.            */
        ma[inext] = mj;               /* Store it                                  */
        
        
        return mj*FAC;                /* and output the derived uniform deviate.   */
}

double unif_rand_dbl(long *idnum)
{
  double highorder = (double) ran3(idnum);
  double loworder = (double) ran3(idnum);
  
  return highorder + loworder*FAC;
}

double box_muller(double m, double s, long *idnum)     
{                                       
  double x1, x2, w, y1;
  static double y2;
  static int use_last = 0;
  
  int nRepetitions = 0;
  if (use_last)
  {
    y1 = y2;
    use_last = 0;
  }
  else
  {
    do {
      x1 = 2.0 * unif_rand_dbl(idnum) - 1.0;
      x2 = 2.0 * unif_rand_dbl(idnum) - 1.0;
      w = x1 * x1 + x2 * x2;
      
      nRepetitions++;
      if(nRepetitions > 100000)
      {
        misc.error("Error: Sorry, an internal error was happened in the random number generator.", 0);
      }
    } while ( w >= 1.0 );
    
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w;
    use_last = 1;
  }
  return( m + y1 * s );
}

double chi1_CDF(int df, double x)
{
  #ifdef BOOSTLIB
    return boost::math::cdf(boost::math::complement(boost::math::chi_squared(df), x));
  #else
    return -1.;
  #endif
}

double FStatCDF(double df1, double df2, double x)
{
#ifdef BOOSTLIB
  return boost::math::cdf(boost::math::complement(boost::math::fisher_f(df1, df2), x));
#else
  return -1.;
#endif
}

double tStatCDF(double df, double x)
{
#ifdef BOOSTLIB
  return boost::math::cdf(boost::math::complement(boost::math::students_t(df), x));
#else
  return -1.;
#endif
}

std::string getFileName(std::string & path)
{
  return path.substr(path.find_last_of("/\\") + 1);
}

/////////////////////////////////////////
//Just for gdb and debugging purposes:

std::string& SSS (const char* s)
{
  return *(new std::string(s));
}

int compareMatrices(Matrix * m1, Matrix * m2, double threshold)
{
  if( (m1->nRows != m2->nRows) || (m1->nCols != m2->nCols) )
  {
    misc.error("Error: Comparing matrices with different dimensions.", 0);
    return 0;
  }
  if( (m1->nBlockRows != m2->nBlockRows) || (m1->nBlockCols != m2->nBlockCols) )
  {
    misc.error("Error: Comparing matrices with different structure.", 0);
    return 0;
  }
  
  int result = 1;
  
  #pragma omp parallel for
  for(int c = 0; c<m1->nCols; c++)
  {
    for(int r = 0; r<m1->nRows; r++)
    {
      if(m1->m[c*m1->nRows + r] != m2->m[c*m1->nRows + r] && threshold < 0.)
      {
        result = 0;
        std::cout << r << " " << c << " " << m1->m[c*m1->nRows + r] << " " << m2->m[c*m1->nRows + r] << " " <<  std::endl;
      }
      if( ( fabs(m1->m[c*m1->nRows + r] - m2->m[c*m1->nRows + r])/std::min(fabs(m1->m[c*m1->nRows + r]), fabs(m2->m[c*m1->nRows + r])) ) > threshold && threshold >= 0. )
      {
        result = 0;
        std::cout << r << " " << c << " " << m1->m[c*m1->nRows + r] << " " << m2->m[c*m1->nRows + r] << " " <<  std::endl;
      }
    }
  }
  
  if(result == 0)
  {
    std::cout << "\n*****************************\nMatrices differ\n*****************************\n" << std::endl;
    return 0;
  }
  return 1;
  
//   int * results = NULL;
//   //Gather the sizes of each array to gather. Only on root. On ohter processes results == NULL.
//   results = gather(&result, 1);
//   
//   if(communicator->mpiRoot)
//   {
//     for(int i = 0; i<communicator->mpiNumTasks; i++)
//     {
//       if(results[i] == 0)
//       {
//         result = 0;
//       }
//     }
//     delete [] results;
//   }
//   communicator->broadcast(result);
//   return result;
}

int compareGlobalMatrices(Matrix * m1, Matrix * m2, double threshold)
{
  std::vector< std::vector<double> > gm1;
  std::vector< std::vector<double> > gm2;
  
  m1->matrixToStandardVector(gm1);
  m2->matrixToStandardVector(gm2);
  
  if( (m1->nGlobRows != m2->nGlobRows) || (m1->nGlobCols != m2->nGlobCols) )
  {
    misc.error("Error: Comparing matrices with different dimensions.", 0);
    return 0;
  }
  
  int result = 1;
  
  if(communicator->mpiRoot)
  {
    for(int c = 0; c<m1->nGlobCols; c++)
    {
      for(int r = 0; r<m1->nGlobRows; r++)
      {
        if(gm1[r][c] != gm2[r][c] && threshold < 0.)
        {
          result = 0;
          std::cout << r << " " << c << " " << gm1[r][c] << " " << gm2[r][c] << " " <<  std::endl;
        }
        if( ( fabs(gm1[r][c] - gm2[r][c])/std::min(fabs(gm1[r][c]), fabs(gm2[r][c])) ) > threshold && threshold >= 0. )
        {
          result = 0;
          std::cout << r << " " << c << " " << gm1[r][c] << " " << gm2[r][c] << " " <<  std::endl;
        }
      }
    }
    
    if(result == 0)
    {
      std::cout << "\n\n\n*****************************\nMatrices differ\n*****************************\n\n\n" << std::endl;
      misc.error("ERROR", 0);
      return 0;
    }
  }
  return 1;
}

GRM * loadGRMUsingOptionsOld(bool returnGenotype, Genotype **returnedGenotype)
{
  GRM *grm = NULL;
  
  if(returnGenotype == true)
  {
    if(returnedGenotype == NULL)
    {
      misc.error("Error: An internal error was happened when loading the GRM. Wrong pointer for returning a genotype.", 0);
    }
    *returnedGenotype = NULL;
  }
  
  if( returnGenotype == true && options.GRMJoinMethod == 0 )
  {
    misc.message << "Changing GRM joining method to method 1." << std::endl;
    options.GRMJoinMethod = 1;
  }
  
  if(options.grmFile == "")
  {
    if(options.genotypeListFile == "") //Create GRM from single genotype file
    {
      Genotype *genotype = new Genotype(options.genotypeFile);
      grm = new GRM(genotype);
      if(returnGenotype == false)
      {
        delete genotype;
      }
      else
      {
        *returnedGenotype = genotype;
      }
    }
    else //Create GRM from multiple genotype files
    {
      if(options.GRMJoinMethod == 0) //Create GRM using method 0: First compute GRMs for each genotype, then add GRMs.
      {
        misc.message << "Computing GRM using a list of genotype files. Method 0: First compute GRMs for each genotype, then add GRMs." << std::endl;
        misc.message.tab = "  ";
        misc.setGetElapsedTime("GRMPartComputation");
        std::ifstream file;
        misc.checkFileExists(options.genotypeListFile);
        file.open(options.genotypeListFile.c_str());
        
        int idx = 0;
        std::string line;
        std::vector<std::string> loadedSNPs;
        while(getline(file,line))
        {
          std::istringstream sstemp(line);
          std::string partialFileName;
          sstemp >> partialFileName;
          
          if( idx == 0 )
          {
            Genotype *genotype = new Genotype(line);
            loadedSNPs = genotype->SNPIds;
            
            grm = new GRM(genotype, false);
            if(returnGenotype == false)
            {
              delete genotype;
            }
            else
            {
              misc.error("Error: An internal error was happened. GRM cannot be loaded using method 0 when genotypes must also be loaded.", 0);
            }
          }
          else
          {
            Genotype *genotype = new Genotype(line);
            std::vector<std::string> SNPsIntersection = intersectionStringVectors(2, &loadedSNPs, &genotype->SNPIds); //Check there are not intersections between SNPs in files.
            if(SNPsIntersection.size() != 0)
            {
              misc.error("Error: When computing a GRM from multiple genotype files. There are SNPs repeated in different files.", 0);
            }
            loadedSNPs.insert( loadedSNPs.end(), genotype->SNPIds.begin(), genotype->SNPIds.end() );
            
            GRM *grmToAdd = new GRM(genotype, false);
            delete genotype;
            misc.message << "Adding GRM..." << std::endl;
            grm->addGRMs(1., grmToAdd);
            delete grmToAdd;
          }
          
          idx++;
        }
        file.close();
        misc.message.tab = "";
        grm->normalize();
        misc.message << "GRM computation finished after " << misc.setGetElapsedTime("GRMPartComputation", true) << std::endl;
      }
      else  //Create GRM using method 1: First join genotypes, then compute GRM.
      {
        misc.message << "Computing GRM using a list of genotype files. Method 1: First join genotypes, then compute GRM." << std::endl;
        misc.message.tab = "  ";
        misc.setGetElapsedTime("GRMPartComputation");
        Genotype *genotype = new Genotype();
        genotype->loadList(options.genotypeListFile);
        grm = new GRM(genotype);
        misc.message.tab = "";
        misc.message << "GRM computation finished after " << misc.setGetElapsedTime("GRMPartComputation", true) << std::endl;
        if(returnGenotype == false)
        {
          delete genotype;
        }
        else
        {
          *returnedGenotype = genotype;
        }
      }
    }
  }
  else
  {
    grm = new GRM(options.grmFile);
  }
  return grm;
}

