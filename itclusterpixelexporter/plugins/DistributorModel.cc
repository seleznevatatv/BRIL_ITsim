/**
 * @file DistributorModel.cc
 * @author Max Bensink (maxbensink@outlook.com)
 * @version 1.0
 * @date 2023-08-01
 *
 * @copyright Copyright (c) 2023
 *
 */


#include <iostream>
#include <vector>

#include "ClusterModel.h"

DistributorModel::DistributorModel(const std::vector<ProcessedQuarterCoreModel> &qcore_buf) : qcore_buf(qcore_buf) {}

std::tuple<std::vector<ClusterModel>, std::vector<ProcessedQuarterCoreModel>> DistributorModel::run()
{
    std::vector<ProcessedQuarterCoreModel> standalone_buf;
    std::vector<ProcessedQuarterCoreModel> row_merger_buf;

    for (const auto &qcore : qcore_buf)
    {
        // Ready to distribute
        if (qcore.is_last_in_event)
        {
            row_merger_buf.push_back(qcore);
        }
        else if (qcore.is_last && !qcore.is_neighbour)
        {
            if (!qcore.is_touching_prev_col && !qcore.is_touching_next_col)
            {
                standalone_buf.push_back(qcore);
            }
            else
            {
                row_merger_buf.push_back(qcore);
            }
        }
        else if (qcore.is_last)
        {
            if (!qcore.is_touching_prev_row && !qcore.is_touching_prev_col && !qcore.is_touching_next_col)
            {
                standalone_buf.push_back(qcore);
            }
            else
            {
                row_merger_buf.push_back(qcore);
            }
        }
        else if (!qcore.is_neighbour)
        {
            if (!qcore.is_touching_next_row && !qcore.is_touching_prev_col && !qcore.is_touching_next_col)
            {
                standalone_buf.push_back(qcore);
            }
            else
            {
                row_merger_buf.push_back(qcore);
            }
        }
        else if (qcore.col == 0)
        {
            if (qcore.row == 0)
            {
                if (!qcore.is_touching_next_col && !qcore.is_touching_next_row)
                {
                    standalone_buf.push_back(qcore);
                }
                else
                {
                    row_merger_buf.push_back(qcore);
                }
            }
            else
            { // Any row
                if (!qcore.is_touching_prev_row && !qcore.is_touching_next_col && !qcore.is_touching_next_row && !qcore.is_touching_prev_col)
                {
                    standalone_buf.push_back(qcore);
                }
                else
                {
                    row_merger_buf.push_back(qcore);
                }
            }
        }
        else
        {
            row_merger_buf.push_back(qcore);
        }
    }

    // Process standalone cluster list
    std::vector<ClusterModel> standalone_cluster_list;
    for (const auto &qcore : standalone_buf)
    {
        for (const auto &cluster : qcore.cluster_list)
        {
            standalone_cluster_list.push_back(cluster);
        }
    }

    // Done
    return std::make_tuple(standalone_cluster_list, row_merger_buf);
}