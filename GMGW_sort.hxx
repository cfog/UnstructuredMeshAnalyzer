/*
 * GMGW_sort.hxx
 *
 *  Created on: Nov 18, 2016
 *      Author: cfog
 */

/*
    Copyright (C) The University of British Columbia, 2018.

    This file is part of UnstructuredMeshAnalyzer.

    UnstructuredMeshAnalyzer is free software: you can redistribute it
    and/or modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation, either version 3 of
    the License, or (at your option) any later version.

    UnstructuredMeshAnalyzer is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with UnstructuredMeshAnalyzer.  If not, see
    <https://www.gnu.org/licenses/>.
*/

#ifndef GMGW_SORT_HXX_
#define GMGW_SORT_HXX_


#define min3(a, b, c) (a < b ? (a < c ? a : c) : (b < c ? b : c))
#define max3(a,b,c) (a > b ? (a > c ? a : c) : (b > c ? b : c))

bool triCompare(const struct triSort& triA, const struct triSort& triB);
bool quadCompare(const struct quadSort& quadA, const struct quadSort& quadB);

struct triSort {
    int m_cell;
    unsigned int m_verts[3];
    void set(const int cell, const unsigned int v0,
             const unsigned int v1, const unsigned int v2)
    {
        m_cell = cell;
        // Sort the verts in ascending order.
        m_verts[0] = min3(v0, v1, v2);
        m_verts[2] = max3(v0, v1, v2);
        m_verts[1] = (v0 + v1 + v2) - m_verts[0] - m_verts[2];
    }
    bool operator==(const triSort& t) const {
        return (m_verts[0] == t.m_verts[0] &&
                m_verts[1] == t.m_verts[1] &&
                m_verts[2] == t.m_verts[2]);
    }
    friend bool operator<(const triSort& triA, const triSort& triB)
    {
        int diff = triA.m_verts[0] - triB.m_verts[0];
        if (diff == 0) {
            diff = triA.m_verts[1] - triB.m_verts[1];
            if (diff == 0) {
                diff = triA.m_verts[2] - triB.m_verts[2];
            }
        }
        return diff < 0;
    }
};

struct triHash {
    size_t operator()(const triSort& ts) const {
        // Hash it!
        static const size_t init =   0xcbf29ce484222325ULL;
        static const size_t factor = 0x00000100000001b3ULL;
        size_t result = init;
        result ^= ts.m_verts[0];
        result *= factor;
        result ^= ts.m_verts[1];
        result *= factor;
        result ^= ts.m_verts[2];
        result *= factor;
        return result;
    }
};

bool triCompare(const triSort& triA, const triSort& triB)
{
    int diff = triA.m_verts[0] - triB.m_verts[0];
    if (diff == 0) {
        diff = triA.m_verts[1] - triB.m_verts[1];
        if (diff == 0) {
            diff = triA.m_verts[2] - triB.m_verts[2];
            if (diff == 0) {
                diff = triA.m_cell - triB.m_cell;
            }
        }
    }
    return diff < 0;
}

struct quadSort {
    int m_cell;
    unsigned int m_verts[4];
    void set(const int cell, const unsigned int v0,
             const unsigned int v1, const unsigned int v2,
             const unsigned int v3)
    {
        m_cell = cell;
        m_verts[0] = v0;
        m_verts[1] = v1;
        m_verts[2] = v2;
        m_verts[3] = v3;
        std::sort(m_verts, m_verts+4);
    }
    bool operator==(const quadSort& q) const {
        return (m_verts[0] == q.m_verts[0] &&
                m_verts[1] == q.m_verts[1] &&
                m_verts[2] == q.m_verts[2] &&
                m_verts[3] == q.m_verts[3]);
    }
    friend bool operator<(const quadSort& quadA, const quadSort& quadB)
    {
        int diff = quadA.m_verts[0] - quadB.m_verts[0];
        if (diff == 0) {
            diff = quadA.m_verts[1] - quadB.m_verts[1];
            if (diff == 0) {
                diff = quadA.m_verts[2] - quadB.m_verts[2];
                if (diff == 0) {
                    diff = quadA.m_verts[3] - quadB.m_verts[3];
                }
            }
        }
        return diff < 0;
    }
};

struct quadHash {
    size_t operator()(const quadSort& qs) const {
        // Hash it!
        static const size_t init =   0xcbf29ce484222325ULL;
        static const size_t factor = 0x00000100000001b3ULL;
        size_t result = init;
        result ^= qs.m_verts[0];
        result *= factor;
        result ^= qs.m_verts[1];
        result *= factor;
        result ^= qs.m_verts[2];
        result *= factor;
        result ^= qs.m_verts[3];
        result *= factor;
        return result;
    }
};

bool quadCompare(const quadSort& quadA, const quadSort& quadB)
{
    int diff = quadA.m_verts[0] - quadB.m_verts[0];
    if (diff == 0) {
        diff = quadA.m_verts[1] - quadB.m_verts[1];
        if (diff == 0) {
            diff = quadA.m_verts[2] - quadB.m_verts[2];
            if (diff == 0) {
                diff = quadA.m_verts[3] - quadB.m_verts[3];
                if (diff == 0) {
                    diff = quadA.m_cell - quadB.m_cell;
                }
            }
        }
    }
    return diff < 0;
}




#endif /* GMGW_SORT_HXX_ */
