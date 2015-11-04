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

#ifndef BIVARREML_H
#define BIVARREML_H

#include "matrix.h"
#include "genotype.h"
#include "covariate.h"
#include "reml.h"
#include "covariancematrix.h"
#include "kernel.h"

#include <vector>
#include <set>
#include <utility>
#include <fstream>

class BivarREML
{
public:
  REML * reml;
  
  bool writeResults;
  
  REMLResults subsampleREMLResults;
  
  BivarREML(bool argWriteResults = true);
  ~BivarREML();
  
  /**
   * Initialize REML parameters
   * 
   * \param grm The grm that will be used for preparing REML. It will be deleted in this function.
   */
  void prepare(Kernel * grm);
  
  /**
   * Initialize REML parameters assuming multiple GRMs
   * 
   * \param grms Pointers to the grms that will be used for preparing REML. All these grms will be deleted and this vector cleared.
   * \param names Names that will be used for grms.
   * \param weights weights for the initial variances of each grm in grms. If empty, 1./grms.size() will be used for each grm.
   */
  void prepare(std::vector<Kernel*> & grms, std::vector<std::string> names, std::vector<double> weights = std::vector<double>());
  
  /**
   * Compute REML
   */
  void compute();
  
  /**
   * Compute regional REML
   * 
   * Compute regional REML. Each computation uses a regional GRM plus a global GRM with the regional part substracted.
   */
  void computeRegional();
  
  /**
   * Compare the results of two REMLs and store a p-value.
   */
  void compareREMLs(double baseLogLikelihood, double baseNVariances);
};

#endif