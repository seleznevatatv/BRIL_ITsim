/**
 * @file QuarterCoreProcessor.cc
 * @author Max Bensink (maxbensink@outlook.com)
 * @version 1.0
 * @date 2023-08-01
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "ClusterModel.h"

#include <algorithm>
#include <ranges>

ProcessedQuarterCoreModel::ProcessedQuarterCoreModel(const QuarterCore &qcore)
{
    // just make a model out of the original quarter core
    hit_map = qcore.get_unsparsified_hit_map();
    col = qcore.col;
    row = qcore.row;
    is_last = qcore.is_last;
    is_neighbour = qcore.is_neighbour;
    is_last_in_event = qcore.is_last_in_event;

    // get touch info
    unpacked_hitmap = qcore.get_unsparsified_hit_map();

    is_touching_prev_col = *std::max_element(unpacked_hitmap[0].begin(), unpacked_hitmap[0].end()) > 0 ? 1 : 0;
    is_touching_next_col = *std::max_element(unpacked_hitmap[SIZE_QCORE_HORIZONTAL - 1].begin(), unpacked_hitmap[SIZE_QCORE_HORIZONTAL - 1].end()) > 0 ? 1 : 0;

    // Check if there is any non-zero value in the first column of unpacked_hitmap
    bool touching_prev_row = std::any_of(unpacked_hitmap.begin(), unpacked_hitmap.end(),
                                         [](const std::vector<int> &row)
                                         { return row[0] > 0; });

    is_touching_prev_row = touching_prev_row ? 1 : 0;

    // Check if there is any non-zero value in the last column of unpacked_hitmap
    bool touching_next_row = std::any_of(unpacked_hitmap.begin(), unpacked_hitmap.end(),
                                         [](const std::vector<int> &row)
                                         {
                                             if (row.size() > SIZE_QCORE_HORIZONTAL)
                                                 return false;
                                             else
                                                 return true;
                                         });

    is_touching_next_row = touching_next_row ? 1 : 0;

    // get for quarter quarter cores
    std::vector<std::vector<std::vector<int>>> quarter_quarter_core_list;

    if (SIZE_QCORE_VERTICAL == 2 && SIZE_QCORE_HORIZONTAL == 8)
    {
        quarter_quarter_core_list.push_back({hit_map[0], hit_map[1]});
        quarter_quarter_core_list.push_back({hit_map[2], hit_map[3]});
        quarter_quarter_core_list.push_back({hit_map[4], hit_map[5]});
        quarter_quarter_core_list.push_back({hit_map[6], hit_map[7]});
    }
    else if (SIZE_QCORE_VERTICAL == 4 && SIZE_QCORE_HORIZONTAL == 4)
    {
        quarter_quarter_core_list.push_back({{hit_map[0][0], hit_map[0][1]}, {hit_map[1][0], hit_map[1][1]}});
        quarter_quarter_core_list.push_back({{hit_map[2][0], hit_map[2][1]}, {hit_map[3][0], hit_map[3][1]}});
        quarter_quarter_core_list.push_back({{hit_map[0][2], hit_map[0][3]}, {hit_map[1][2], hit_map[1][3]}});
        quarter_quarter_core_list.push_back({{hit_map[2][2], hit_map[2][3]}, {hit_map[3][2], hit_map[3][3]}});
    }
    else
    {
        std::cout << "ERROR wrong size of the qcore" << std::endl;
        exit(1);
    }

    // stage 0 (does not depend on the qcore size)

    std::vector<int> has_cluster(4, 0);
    for (size_t qqcore_id = 0; qqcore_id < quarter_quarter_core_list.size(); ++qqcore_id)
    {
        const std::vector<std::vector<int>> &quarter_quarter_core = quarter_quarter_core_list[qqcore_id];
        for (const auto &col_vec : quarter_quarter_core)
        {
            if (std::any_of(col_vec.begin(), col_vec.end(), [](int val)
                            { return val == 1; }))
            {
                has_cluster[qqcore_id] = 1;
                break;
            }
        }
    }

    // stage 1
    std::vector<ClusterModel> cluster_list_0;
    if (has_cluster[0] && has_cluster[1] && (*std::max_element(quarter_quarter_core_list[0][1].begin(), quarter_quarter_core_list[0][1].end()) > 0) && (*std::max_element(quarter_quarter_core_list[1][0].begin(), quarter_quarter_core_list[1][0].end()) > 0))
    {
        ClusterModel temp_cluster(col, row);
        for (int qqcore_id = 0; qqcore_id < 2; ++qqcore_id)
        {
            int qcore_offset_col = qqcore_id * 2;
            int qcore_offset_row = 0;
            for (int col_id = 0; col_id < 2; ++col_id)
            {
                for (int row_id = 0; row_id < 2; ++row_id)
                {
                    if (quarter_quarter_core_list[qqcore_id][col_id][row_id] > 0)
                    {
                        temp_cluster.add_hit(qcore_offset_col + col_id, qcore_offset_row + row_id);
                    }
                }
            }
        }
        cluster_list_0.push_back(temp_cluster);
    }
    else
    {
        for (int qqcore_id = 0; qqcore_id < 2; ++qqcore_id)
        {
            int qcore_offset_col = qqcore_id * 2;
            int qcore_offset_row = 0;
            if (has_cluster[qqcore_id])
            {
                ClusterModel temp_cluster(col, row);
                for (int col_id = 0; col_id < 2; ++col_id)
                {
                    for (int row_id = 0; row_id < 2; ++row_id)
                    {
                        if (quarter_quarter_core_list[qqcore_id][col_id][row_id] > 0)
                        {
                            temp_cluster.add_hit(qcore_offset_col + col_id, qcore_offset_row + row_id);
                        }
                    }
                }
                cluster_list_0.push_back(temp_cluster);
            }
        }
    }

    std::vector<ClusterModel> cluster_list_1;
    if (has_cluster[2] && has_cluster[3] && (*std::max_element(quarter_quarter_core_list[2][1].begin(), quarter_quarter_core_list[2][1].end()) > 0) && (*std::max_element(quarter_quarter_core_list[3][0].begin(), quarter_quarter_core_list[3][0].end()) > 0))
    {
        ClusterModel temp_cluster(col, row);
        for (int qqcore_id = 2; qqcore_id < 4; ++qqcore_id)
        {
            int qcore_offset_col = (SIZE_QCORE_VERTICAL == 2 && SIZE_QCORE_HORIZONTAL == 8) ? qqcore_id * 2 : (qqcore_id - 2) * 2;
            int qcore_offset_row = (SIZE_QCORE_VERTICAL == 2 && SIZE_QCORE_HORIZONTAL == 8) ? 0 : 2;
            for (int col_id = 0; col_id < 2; ++col_id)
            {
                for (int row_id = 0; row_id < 2; ++row_id)
                {
                    if (quarter_quarter_core_list[qqcore_id][col_id][row_id] > 0)
                    {
                        temp_cluster.add_hit(qcore_offset_col + col_id, qcore_offset_row + row_id);
                    }
                }
            }
        }
        cluster_list_1.push_back(temp_cluster);
    }
    else
    {
        for (int qqcore_id = 2; qqcore_id < 4; ++qqcore_id)
        {
            int qcore_offset_col = (SIZE_QCORE_VERTICAL == 2 && SIZE_QCORE_HORIZONTAL == 8) ? qqcore_id * 2 : (qqcore_id - 2) * 2;
            int qcore_offset_row = (SIZE_QCORE_VERTICAL == 2 && SIZE_QCORE_HORIZONTAL == 8) ? 0 : 2;
            if (has_cluster[qqcore_id])
            {
                ClusterModel temp_cluster(col, row);
                for (int col_id = 0; col_id < 2; ++col_id)
                {
                    for (int row_id = 0; row_id < 2; ++row_id)
                    {
                        if (quarter_quarter_core_list[qqcore_id][col_id][row_id] > 0)
                        {
                            temp_cluster.add_hit(qcore_offset_col + col_id, qcore_offset_row + row_id);
                        }
                    }
                }
                cluster_list_1.push_back(temp_cluster);
            }
        }
    }

    // final cluster list
    cluster_list.reserve(cluster_list_0.size() + cluster_list_1.size());

    // stage 3
    std::vector<int> cluster_was_used_0(2, 0);
    std::vector<int> cluster_was_used_1(2, 0);
    for (size_t cluster_0_id = 0; cluster_0_id < cluster_list_0.size(); ++cluster_0_id)
    {
        ClusterModel &cluster_0 = cluster_list_0[cluster_0_id];
        for (size_t cluster_1_id = 0; cluster_1_id < cluster_list_1.size(); ++cluster_1_id)
        {
            ClusterModel &cluster_1 = cluster_list_1[cluster_1_id];
            bool is_merged = false;
            for (const auto &hit_0 : cluster_0.hit_map)
            {
                for (const auto &hit_1 : cluster_1.hit_map)
                {
                    if (abs(hit_0.first - hit_1.first) <= 1 && abs(hit_0.second - hit_1.second) <= 1)
                    {
                        is_merged = true;
                        break;
                    }
                }
                if (is_merged)
                {
                    break;
                }
            }
            if (is_merged)
            {
                cluster_was_used_0[cluster_0_id] = 1;
                cluster_was_used_1[cluster_1_id] = 1;
                ClusterModel new_cluster(col, row);
                for (const auto &hit_1 : cluster_1.hit_map)
                {
                    new_cluster.add_hit(hit_1.first, hit_1.second);
                }
                cluster_list.push_back(new_cluster);
                break;
            }
        }
    }

    // append the ones not merged
    for (size_t cluster_0_id = 0; cluster_0_id < cluster_list_0.size(); ++cluster_0_id)
    {
        if (cluster_was_used_0[cluster_0_id] == 0)
        {
            cluster_list.push_back(cluster_list_0[cluster_0_id]);
        }
    }

    for (size_t cluster_1_id = 0; cluster_1_id < cluster_list_1.size(); ++cluster_1_id)
    {
        if (cluster_was_used_1[cluster_1_id] == 0)
        {
            cluster_list.push_back(cluster_list_1[cluster_1_id]);
        }
    }
}