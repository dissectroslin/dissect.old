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

#include "message.h"
#include "communicator.h"
#include "global.h"
#include "misc.h"

#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

Message::Message()
{
  this->output = &std::cout;
  this->log = NULL;
  this->logActive = false;
  this->outputIsFile = false;
  
  this->tab = "";
  this->newLine = false;
}

Message::Message(std::ostream * out)
{
  this->output = out;
  this->log = NULL;
  this->logActive = false;
  this->outputIsFile = false;
  
  this->tab = "";
  this->newLine = false;
}

Message::Message(std::ostream * out, std::ostream * logstream)
{
  this->output = out;
  this->log = logstream;
  this->logActive = true;
  this->outputIsFile = false;
  
  this->tab = "";
  this->newLine = false;
}

Message::Message(std::string fname)
{
  std::ofstream * file = new std::ofstream();
  misc.openOutputFileInRoot(*file, fname);
  this->output = file;
  this->log = NULL;
  this->logActive = false;
  this->outputIsFile = true;
  
  this->tab = "";
  this->newLine = false;
}

Message::~Message()
{
  if(this->outputIsFile == true)
  {
    if(communicator->mpiRoot)
    {
      ((std::ofstream*)this->output)->close();
    }
  }
}

void Message::redirect(std::string fname)
{
  if(this->outputIsFile == true)
  {
    if(communicator->mpiRoot)
    {
      ((std::ofstream*)this->output)->close();
    }
    delete this->output;
    std::ofstream * file = new std::ofstream();
    misc.openOutputFileInRoot(*file, fname);
    this->output = file;
  }
  else
  {
    misc.error("Error: An internal error was happened. This message can not be reconverted to a file type.", 0);
  }
}

void Message::setWidth(int width)
{
  this->output->width(width);
  if(this->logActive == true)
  {
    this->log->width(width);
  }
}

void Message::flush()
{
  this->output->flush();
  if(this->logActive)
  {
    this->log->flush();
  }
}

Message& Message::operator<<(std::ostream& (*fn)(std::ostream&))
{
  this->newLine = true;
  if(communicator->mpiRoot)
  {
    *(this->output) << fn;
    if(this->logActive)
    {
      *(this->log) << fn;
    }
  }
  return *this;
}

Message& Message::operator<<(std::_Setw manip)
{
  if(communicator->mpiRoot)
  {
    *(this->output) << manip;
    if(this->logActive)
    {
      *(this->log) << manip;
    }
  }
  return *this;
}

Message& Message::operator<<(std::_Setprecision manip)
{
  if(communicator->mpiRoot)
  {
    *(this->output) << manip;
    if(this->logActive)
    {
      *(this->log) << manip;
    }
  }
  return *this;
}
