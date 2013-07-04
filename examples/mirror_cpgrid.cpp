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

#include <iostream>

#include <opm/core/io/eclipse/EclipseGridParser.hpp>


/**
 * @file mirror_cpgrid.cpp
 * @brief Mirror corner-point grid taken from grdecl file
 * 
 * Results in a periodic grid
 *
 * @author Lars Vingli Odsæter <laods@statoil.com>
 *
 */

using namespace Opm;
using namespace std;


/// Print init message in new grid filename
void printInitMessage(ofstream& out, const char* origfilename, string direction) {
    ifstream infile;
    infile.open(origfilename, ios::in);
    if (!infile) {
        cerr << "Can't open input file " << string(origfilename) << endl;
        exit(1);
    }
    // Print init message and copy comments from original grid file
    out << "-- This grdecl file is generated from OPM application 'mirror_cpgrid.cpp' (Author: Lars Vingli Odsæter)." << endl
        << "-- The grid '" << string(origfilename) << "' is mirrored around itself in the " << direction << " direction." << endl
        << "-- Thus, the resulting grid should be periodic in the " << direction << "-direction." << endl
        << "-- Original comments taken from '" << string(origfilename) << "':" << endl;
    string nextInLine;
    while (getline(infile, nextInLine)) {
        if (nextInLine.substr(0,2) == "--") {
            out << nextInLine << endl;
        }
        else {
            break;
        }
    }
    out << endl;   
}

/// Write keyword values to file
template <class T>
void printKeywordValues(ofstream& out, string keyword, vector<T> values, int nCols) {
    out << keyword << endl;
    int col = 0;
    typename std::vector<T>::iterator iter;
    for (iter = values.begin(); iter != values.end(); ++iter) {
        out << *iter << " ";            
        ++col;
        // Break line for every nCols entry.
        if (col == nCols) {
            out << endl;
            col = 0;
        }
    }
    if (col != 0)
        out << endl;
    out << "/" << endl << endl;
}

void mirror_mapaxes(EclipseGridParser parser, string direction, ofstream& out) {
    // Assumes axis aligned with x/y-direction
    cout << "Warning: Keyword MAPAXES not fully understood. Result should be verified manually." << endl;
    if (parser.hasField("MAPAXES")) {
        vector<double> mapaxes = parser.getFloatingPointValue("MAPAXES");
        vector<double> mapaxes_mirrored = mapaxes;
        if (direction == "x") {
            mapaxes_mirrored[4] = (mapaxes[4]-mapaxes[2])*2 + mapaxes[2];
        }
        else if (direction == "y") {
            mapaxes_mirrored[1] = (mapaxes[1]-mapaxes[3])*2 + mapaxes[3];
        }
        printKeywordValues(out, "MAPAXES", mapaxes_mirrored, 2);
    }
}

void mirror_specgrid(EclipseGridParser parser, string direction, ofstream& out) {
    SPECGRID specgrid = parser.getSPECGRID();
    vector<int> dim = specgrid.dimensions;
    if (direction == "x")      {dim[0] *= 2;}
    else if (direction == "y") {dim[1] *= 2;}
    else                       {cerr << "Direction should be either x or y" << endl; exit(1);}
    out << "SPECGRID" << endl << dim[0] << " " << dim[1] << " " << dim[2] << " "
        << specgrid.numres << " " << specgrid.qrdial << endl << "/" << endl << endl;
}

void mirror_coord(EclipseGridParser parser, string direction, ofstream& out) {
    
}

void mirror_zcorn(EclipseGridParser parser, string direction, ofstream& out) {
    
}

void mirror_celldata(string keyword, EclipseGridParser parser, string direction, ofstream& out) {
    // Handle ACTNUM and SATNUM specially.
}
    
    
int main(int argc, char** argv)
{
    // Set output precision
    int decimals = 4;
    
    // Process input parameters
    if (argc != 3) {
        cout << "Usage: mirror_cpgrid filename.grdecl direction" << endl;
        exit(1);
    }
    const char* eclipsefilename = argv[1];
    string direction = string(argv[2]);
    if ( ! ((direction == "x") || (direction == "y")) ) {
        cerr << "Unrecognized input parameter for direction: '" << direction 
             << "'. Should be either x or y (maybe also z later)." << endl;
        exit(1);
    }
    
    // Parse grdecl file
    cout << "Parsing grid file '" << eclipsefilename << "' ..." << endl;
    bool convert_to_SI = false;
    EclipseGridParser eclParser(eclipsefilename, convert_to_SI);
    if ( ! (eclParser.hasField("SPECGRID") && eclParser.hasField("COORD") && eclParser.hasField("ZCORN")) ) {
        cerr << "Grid file " << eclipsefilename << "are missing keywords SPECGRID, COORD or ZCORN!" << endl;
        exit(1);
    }
    
    // Create new grid file
    string mirrored_eclipsefilename = string(eclipsefilename);
    string::size_type last_dot = mirrored_eclipsefilename.find_last_of('.');
    mirrored_eclipsefilename = mirrored_eclipsefilename.substr(0, last_dot) + "_mirrored-" + direction + ".grdecl";
    ofstream outfile;
    outfile.open(mirrored_eclipsefilename, ios::out | ios::trunc);
    if (!outfile) {
        cerr << "Can't open output file " << mirrored_eclipsefilename << endl;
        exit(1);
    }
    outfile.precision(decimals);
    outfile.setf(ios::fixed);

    // Print init message
    printInitMessage(outfile, eclipsefilename, direction);
    
    // Mirror keywords    
    mirror_mapaxes(eclParser, direction, outfile);
    mirror_specgrid(eclParser, direction, outfile);
    mirror_coord(eclParser, direction, outfile);
    mirror_zcorn(eclParser, direction, outfile);
    mirror_celldata("ACTNUM", eclParser, direction, outfile);
    mirror_celldata("PERMX", eclParser, direction, outfile);
    mirror_celldata("PORO", eclParser, direction, outfile);
    mirror_celldata("SATNUM", eclParser, direction, outfile);

}