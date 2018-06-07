/**
 * GeoDa TM, Copyright (C) 2011-2015 by Luc Anselin - all rights reserved
 *
 * This file is part of GeoDa.
 * 
 * GeoDa is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GeoDa is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __GEODA_CENTER_DISTANCES_CALC_H__
#define __GEODA_CENTER_DISTANCES_CALC_H__

#include <map>
#include <boost/bimap.hpp>
#include <wx/string.h>

/** We ultimately need all distance pairs, sorted by distance. */

struct UnOrdIntPair {
	UnOrdIntPair();
	UnOrdIntPair(int i, int j);
	bool operator<(const UnOrdIntPair& s) const;
	bool operator==(const UnOrdIntPair& s) const;
	int i;
	int j;
	wxString toStr();
};

typedef boost::bimap<int, UnOrdIntPair> pairs_bimap_type;
typedef std::map<UnOrdIntPair, double> dist_map_type;

class DistancesCalc {
public:
	
	
};

#endif
