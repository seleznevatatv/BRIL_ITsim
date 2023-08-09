/**
 * @file ClusterModel.cc
 * @author Max Bensink (maxbensink@outlook.com)
 * @version 1.0
 * @date 2023-08-01
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "ClusterModel.h"

#include <algorithm>

ClusterModel::ClusterModel(int col, int row)
    : col(col), row(row), ncols(1), nrows(1) {}

void ClusterModel::add_hit(int local_col, int local_row)
{
    if (local_col >= SIZE_QCORE_HORIZONTAL * (ncols - 1) && local_col < SIZE_QCORE_HORIZONTAL * ncols &&
        local_row >= SIZE_QCORE_VERTICAL * (nrows - 1) && local_row < SIZE_QCORE_VERTICAL * nrows)
    {
        hit_map.push_back(std::make_pair(local_col, local_row));
    }
    else
    {
        throw std::runtime_error(std::string("ERROR wrong col row for the cluster model ") + std::to_string(local_col) + " " + std::to_string(local_row) + " " + std::to_string(ncols) + " " + std::to_string(nrows));
    }
}

void ClusterModel::next_col()
{
    ncols++;
}

void ClusterModel::next_row()
{
    nrows++;
}

bool ClusterModel::is_touching_prev_col() const
{
    for (auto hit : hit_map)
    {
        if (hit.first == 0)
        {
            return true;
        }
    }
    return false;
}

bool ClusterModel::is_touching_next_col() const
{
    for (auto hit : hit_map)
    {
        if (hit.first == SIZE_QCORE_HORIZONTAL * ncols - 1)
        {
            return true;
        }
    }
    return false;
}

bool ClusterModel::is_touching_prev_row() const
{
    for (auto hit : hit_map)
    {
        if (hit.second == 0)
        {
            return true;
        }
    }
    return false;
}

bool ClusterModel::is_touching_next_row() const
{
    for (auto hit : hit_map)
    {
        if (hit.second == SIZE_QCORE_VERTICAL * nrows - 1)
        {
            return true;
        }
    }
    return false;
}

std::pair<int, int> ClusterModel::get_hit_global_coordinate(const std::pair<int, int> &hit) const
{
    if (std::find(hit_map.begin(), hit_map.end(), hit) == hit_map.end())
    {
        throw std::runtime_error("ERROR did not find this hit in the hit map");
    }
    else
    {
        int col_global = hit.first + col * SIZE_QCORE_HORIZONTAL;
        int row_global = hit.second + row * SIZE_QCORE_VERTICAL;
        return std::make_pair(col_global, row_global);
    }
}

std::string ClusterModel::to_string() const
{
    std::string res = "Cluster, col = " + std::to_string(col) + ", ncols = " + std::to_string(ncols) +
                      ", row = " + std::to_string(row) + ", nrows = " + std::to_string(nrows) + "\n";
    res += "\tHits LOCAL coordinates: ";
    for (auto hit : hit_map)
    {
        res += "(" + std::to_string(hit.first) + "," + std::to_string(hit.second) + "), ";
    }
    res += "\n";
    res += "\tHits GLOBAL coordinates: ";
    for (auto hit : hit_map)
    {
        res += "(" + std::to_string(hit.first + col * SIZE_QCORE_HORIZONTAL) + "," +
               std::to_string(hit.second + row * SIZE_QCORE_VERTICAL) + "), ";
    }
    res += "\n";
    return res;
}