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
 * The input grid is mirrored in either the x- or y-direction, resulting in a periodic grid.
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

/// Mirror keyword MAPAXES in parser
void mirror_mapaxes(EclipseGridParser parser, string direction, ofstream& out) {
    // Assumes axis aligned with x/y-direction
    cout << "Warning: Keyword MAPAXES not fully understood. Result should be verified manually." << endl;
    if (parser.hasField("MAPAXES")) {
        vector<double> mapaxes = parser.getFloatingPointValue("MAPAXES");
        vector<double> mapaxes_mirrored = mapaxes;
        // Double the length of the coordinate axis
        if (direction == "x") {
            mapaxes_mirrored[4] = (mapaxes[4]-mapaxes[2])*2 + mapaxes[2];
        }
        else if (direction == "y") {
            mapaxes_mirrored[1] = (mapaxes[1]-mapaxes[3])*2 + mapaxes[3];
        }
        printKeywordValues(out, "MAPAXES", mapaxes_mirrored, 2);
    }
}

/// Mirror keyword SPECGRID in parser
void mirror_specgrid(EclipseGridParser parser, string direction, ofstream& out) {
    // We only need to multiply the dimension by 2 in the correct direction.
    SPECGRID specgrid = parser.getSPECGRID();
    vector<int> dim = specgrid.dimensions;
    if (direction == "x")      {dim[0] *= 2;}
    else if (direction == "y") {dim[1] *= 2;}
    else                       {cerr << "Direction should be either x or y" << endl; exit(1);}
    out << "SPECGRID" << endl << dim[0] << " " << dim[1] << " " << dim[2] << " "
        << specgrid.numres << " " << specgrid.qrdial << endl << "/" << endl << endl;
}

/// Mirror keyword COORD in parser
void mirror_coord(EclipseGridParser parser, string direction, ofstream& out) {
    // We assume uniform spacing in x and y directions and parallel top and bottom faces
    vector<int> dimensions = parser.getSPECGRID().dimensions;
    vector<double> coord = parser.getFloatingPointValue("COORD");
    const int entries_per_pillar = 6;
    vector<double> coord_mirrored;
    // Handle the two directions differently due to ordering of the pillars.
    if (direction == "x") {
        // Total entries in mirrored ZCORN. Number of pillars times 6
        const int entries = (2*dimensions[0] + 1) * (dimensions[1] + 1) * entries_per_pillar;
        // Entries per line in x-direction. Number of pillars in x-direction times 6
        const int entries_per_line = entries_per_pillar*(dimensions[0] + 1);
        coord_mirrored.assign(entries, 0.0);
        // Distance between pillars in x-directiion
        const double spacing = coord[entries_per_pillar]-coord[0];
        vector<double>::iterator it_new = coord_mirrored.begin();
        vector<double>::iterator it_orig;
        // Loop through each pillar line in the x-direction
        for (it_orig = coord.begin(); it_orig != coord.end(); it_orig += entries_per_line) {
            // Copy old pillars
            copy(it_orig, it_orig + entries_per_line, it_new);
            // Add new pillars in between
            it_new += entries_per_line;
            vector<double> next_vec(it_orig + entries_per_line - entries_per_pillar, it_orig + entries_per_line);
            for (int r=0; r < dimensions[0]; ++r) {
                next_vec[0] += spacing;
                next_vec[3] += spacing;
                copy(next_vec.begin(), next_vec.end(), it_new);
                it_new += entries_per_pillar;
            }
        }
    }
    else if (direction == "y") {
        // Total entries in mirrored ZCORN. Number of pillars times 6
        const int entries = (dimensions[0] + 1) * (2*dimensions[1] + 1) * entries_per_pillar;
        // Entries per line in y-direction. Number of pillars in y-direction times 6
        const int entries_per_line = entries_per_pillar*(dimensions[1] + 1);
        coord_mirrored.assign(entries, 0.0);
        // Distance between pillars in y-directiion
        const double spacing = coord[entries_per_line + 1]-coord[1];
        vector<double>::iterator it_new = coord_mirrored.begin();
        // Copy old pillars
        copy(coord.begin(), coord.end(), it_new);
        // Add new pillars at the end
        it_new += coord.size();
        vector<double> next_vec(coord.end() - entries_per_line, coord.end());
        for ( ; it_new != coord_mirrored.end(); it_new += entries_per_line) {
            for (int i = 1; i < entries_per_line; i += 3) {
                next_vec[i] += spacing;
            }
            copy(next_vec.begin(), next_vec.end(), it_new);
        }
    }
    else {
        cerr << "Direction should be either x or y" << endl;
        exit(1);
    }
    // Write new COORD values to output file
    printKeywordValues(out, "COORD", coord_mirrored, 6);
}

/// Mirror keyword ZCORN in parser
void mirror_zcorn(EclipseGridParser parser, string direction, ofstream& out) {
    vector<int> dimensions = parser.getSPECGRID().dimensions;
    vector<double> zcorn = parser.getFloatingPointValue("ZCORN");
    vector<double> zcorn_mirrored;
    // Handle the two directions differently due to ordering of the pillars.
    if (direction == "x") {
        // Total entries in mirrored ZCORN. Eight corners per cell.
        const int entries = dimensions[0]*2*dimensions[1]*dimensions[2]*8;
        zcorn_mirrored.assign(entries, 0.0);
        // Entries per line in x-direction. Two for each cell.
        const int entries_per_line = dimensions[0]*2;
        vector<double>::iterator it_new = zcorn_mirrored.begin();
        vector<double>::iterator it_orig = zcorn.begin();
        // Loop through each line and copy old corner-points and add new (which are the old reversed)
        for ( ; it_orig != zcorn.end(); it_orig += entries_per_line) {
            vector<double> next_vec(it_orig, it_orig + entries_per_line);
            vector<double> next_reversed = next_vec;
            reverse(next_reversed.begin(), next_reversed.end());
            // Copy old corner-points
            copy(it_orig, it_orig + entries_per_line, it_new);
            it_new += entries_per_line;
            // Add new corner-points
            copy(next_reversed.begin(), next_reversed.end(), it_new);
            it_new += entries_per_line;
        }
    }
    else if (direction == "y") {
        // Total entries in mirrored ZCORN. Eight corners per cell.
        const int entries = dimensions[0]*dimensions[1]*2*dimensions[2]*8;
        zcorn_mirrored.assign(entries, 0.0);
        // Entries per line in x-direction. Two for each cell.
        const int entries_per_line_x = dimensions[0]*2;
        // Entries per layer of corner-points. Four for each cell
        const int entries_per_layer = dimensions[0]*dimensions[1]*4;
        vector<double>::iterator it_new = zcorn_mirrored.begin();
        vector<double>::iterator it_orig = zcorn.begin();
        // Loop through each layer and copy old corner-points and add new (which are the old reordered) 
        for ( ; it_orig != zcorn.end(); it_orig += entries_per_layer) {
            // Copy old corner-points
            copy(it_orig, it_orig + entries_per_layer, it_new);
            it_new += entries_per_layer;
            // Add new corner-points
            vector<double> next_vec(it_orig, it_orig + entries_per_layer);
            vector<double> next_reordered(entries_per_layer, 0.0);
            vector<double>::iterator it_next = next_vec.end();
            vector<double>::iterator it_reordered = next_reordered.begin();
            // Reorder next entries
            for ( ; it_reordered != next_reordered.end(); it_reordered += entries_per_line_x) {
                copy(it_next - entries_per_line_x, it_next, it_reordered);
                it_next -= entries_per_line_x;
            }
            copy(next_reordered.begin(), next_reordered.end(), it_new);
            it_new += entries_per_layer;
        }
    }
    else {
        cerr << "Direction should be either x or y" << endl;
        exit(1);
    }
    // Write new ZCORN values to output file
    printKeywordValues(out, "ZCORN", zcorn_mirrored, 8);
}

vector<int> getKeywordValues(string keyword, EclipseGridParser parser, int dummy) {
    return parser.getIntegerValue(keyword);
}

vector<double> getKeywordValues(string keyword, EclipseGridParser parser, double dummy) {
    return parser.getFloatingPointValue(keyword);
}

/// Mirror keywords that have one value for each cell
template <class T>
void mirror_celldata(string keyword, EclipseGridParser parser, string direction, ofstream& out) {
    if ( ! parser.hasField(keyword)) {
        cout << "Ignoring keyword " << keyword << " as it was not found." << endl;
        return;
    }
    // Get data from eclipse parser
    vector<int> dimensions = parser.getSPECGRID().dimensions;
    vector<T> values = getKeywordValues(keyword, parser, T(0.0));
    vector<T> values_mirrored(2*dimensions[0]*dimensions[1]*dimensions[2], 0.0);
    // Handle the two directions differently due to ordering of the pillars.
    if (direction == "x") {
        typename vector<T>::iterator it_orig = values.begin();
        typename vector<T>::iterator it_new = values_mirrored.begin();
        // Loop through each line and copy old cell data and add new (which are the old reversed)
        for ( ; it_orig != values.end(); it_orig += dimensions[0]) {
            // Copy old cell data
            copy(it_orig, it_orig + dimensions[0], it_new);
            it_new += dimensions[0];
            // Add new cell data
            vector<double> next_vec(it_orig, it_orig + dimensions[0]);
            vector<double> next_reversed = next_vec;
            reverse(next_reversed.begin(), next_reversed.end());
            copy(next_reversed.begin(), next_reversed.end(), it_new);
            it_new += dimensions[0];
        }
    }
    else if (direction =="y") {
        typename vector<T>::iterator it_orig = values.begin();
        typename vector<T>::iterator it_new = values_mirrored.begin();
        // Entries per layer
        const int entries_per_layer = dimensions[0]*dimensions[1];
        // Loop through each layer and copy old cell data and add new (which are the old reordered) 
        for ( ; it_orig != values.end(); it_orig += entries_per_layer) {
            // Copy old cell data
            copy(it_orig, it_orig + entries_per_layer, it_new);
            it_new += entries_per_layer;
            // Add new cell data
            vector<T> next_vec(it_orig, it_orig + entries_per_layer);
            vector<T> next_reordered(entries_per_layer, 0.0);
            typename vector<T>::iterator it_next = next_vec.end();
            typename vector<T>::iterator it_reordered = next_reordered.begin();
            // Reorder next entries
            for ( ; it_reordered != next_reordered.end(); it_reordered += dimensions[0]) {
                copy(it_next - dimensions[0], it_next, it_reordered);
                it_next -= dimensions[0];
            }
            copy(next_reordered.begin(), next_reordered.end(), it_new);
            it_new += entries_per_layer;
        }
    }
    else {
        cerr << "Direction should be either x or y" << endl;
        exit(1);
    }
    // Write new keyword values to output file
    printKeywordValues(out, keyword, values_mirrored, 8);
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
    mirror_celldata<int>("ACTNUM", eclParser, direction, outfile);
    mirror_celldata<double>("PERMX", eclParser, direction, outfile);
    mirror_celldata<double>("PERMY", eclParser, direction, outfile);
    mirror_celldata<double>("PERMZ", eclParser, direction, outfile);
    mirror_celldata<double>("PORO", eclParser, direction, outfile);
    mirror_celldata<int>("SATNUM", eclParser, direction, outfile);

}