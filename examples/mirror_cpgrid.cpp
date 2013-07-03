/*
  Copyright 2009, 2010, 2013 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010, 2013 Statoil ASA.

  This file is part of The Open Porous Media project  (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#if HAVE_CONFIG_H
#include "config.h"
#endif

//#include <dune/grid/CpGrid.hpp>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>


/**
 * @file mirror_cpgrid.cpp
 * @brief Mirror corner-point grid taken from grdecl file
 * 
 * Results in a periodic grid
 *
 * @author Håvard Berland <havb@statoil.com>
 * @author Atgeirr F Rasmussen <atgeirr@sintef.no>
 * @author Bård Skaflestad     <bard.skaflestad@sintef.no>
 *
 */


int main(int argc, char** argv)
{

    // Process input parameters
    if (argc != 3) {
        std::cout << "Usage: mirror_cpgrid filename.grdecl direction" << std::endl;
        exit(1);
    }
    const char* eclipsefilename = argv[1];
    const char* direction = argv[2];  
    
    bool convert_to_SI = false;
    Opm::EclipseGridParser eclParser(eclipsefilename, convert_to_SI);
    
}