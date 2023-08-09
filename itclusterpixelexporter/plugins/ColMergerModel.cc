/**
 * @file ColMergerModel.cc
 * @author Max Bensink (maxbensink@outlook.com)
 * @version 1.0
 * @date 2023-08-01
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "ClusterModel.h"

ColMergerModel::ColMergerModel(const std::vector<ClusterModel> &col_merger_buf) : col_merger_buf(col_merger_buf), _current_state("Unknown"), not_adjacent(0), already_used(0), not_found(0), rhs_must_be_bigger(0)
{
    _state_map["FillColumn1"] = [this]()
    { _FillColumn1(); };
    _state_map["FillColumn2"] = [this]()
    { _FillColumn2(); };
    _state_map["MergeColumns"] = [this]()
    { _MergeColumns(); };
    _state_map["Done"] = [this]()
    { _Done(); };
}

std::vector<ClusterModel> ColMergerModel::run()
{
    _current_state = "FillColumn1";
    processed_cluster_list.clear();

    while (true)
    {
        // Run state machine
        auto func = _state_map.find(_current_state);
        if (func != _state_map.end())
        {
            func->second();
        }
        else
        {
            // Unknown state
            break;
        }

        // Check if we are done
        if (_current_state == "Done")
        {
            break;
        }
    }

    // Return
    return processed_cluster_list;
}

void ColMergerModel::_FillColumn1()
{
    column_buf_1.clear();

    if (!col_merger_buf.empty())
    {
        int prev_col = col_merger_buf[0].col;

        while (true)
        {
            if (!col_merger_buf.empty())
            {
                ClusterModel cluster = col_merger_buf[0];
                if (cluster.col == prev_col)
                {
                    // Once more - now pop
                    cluster = col_merger_buf.front();
                    col_merger_buf.erase(col_merger_buf.begin());

                    if (cluster.is_touching_next_col())
                    {
                        column_buf_1.push_back(cluster);
                    }
                    else
                    {
                        processed_cluster_list.push_back(cluster);
                    }
                }
                else
                {
                    if (!column_buf_1.empty())
                    {
                        // Switch to the next col
                        _current_state = "FillColumn2";
                    }
                    else
                    {
                        // The cluster was not touching the next column - go back to the same state
                        _current_state = "FillColumn1";
                    }
                    break;
                }
            }
            else
            {
                _empty_column_1();
                _current_state = "Done";
                break;
            }
        }
    }
    else
    {
        _current_state = "Done";
    }
}

void ColMergerModel::_empty_column_1()
{
    for (const auto &cluster : column_buf_1)
    {
        processed_cluster_list.push_back(cluster);
    }
    column_buf_1.clear();
}

void ColMergerModel::_FillColumn2()
{
    column_buf_2.clear();

    if (!col_merger_buf.empty())
    {
        // Set the column that we expect
        int current_col = column_buf_1[0].col + column_buf_1[0].ncols;

        while (true)
        {
            // Until it is empty
            if (!col_merger_buf.empty())
            {
                ClusterModel cluster = col_merger_buf[0];
                if (cluster.col == current_col)
                {
                    cluster = col_merger_buf.front();
                    col_merger_buf.erase(col_merger_buf.begin());

                    column_buf_2.push_back(cluster);
                }
                else
                {
                    if (!column_buf_2.empty())
                    {
                        _current_state = "MergeColumns";
                    }
                    else
                    {
                        _empty_column_1();
                        _current_state = "FillColumn1";
                    }
                    break;
                }
            }
            else
            {
                // Check that we managed to put something in the buffer
                if (!column_buf_2.empty())
                {
                    _current_state = "MergeColumns";
                }
                else
                {
                    _empty_column_1();
                    _current_state = "Done";
                }
                break;
            }
        }
    }
    else
    {
        _empty_column_1();
        _current_state = "Done";
    }
}

ClusterModel ColMergerModel::_merge_two_clusters(const ClusterModel &cluster1, const ClusterModel &cluster2)
{
    // Check the assumption
    if (cluster2.col < cluster1.col)
    {
#ifdef THROW_ERROR
        throw std::runtime_error("ERROR (column): This method assumes that cluster2.col is always bigger");
#else
        std::cerr << "ERROR (column): This method assumes that cluster2.col is always bigger" << std::endl;

        rhs_must_be_bigger++;
#endif
    }

#ifdef THROW_ERROR
    const ClusterModel &lcluster = cluster1;
    const ClusterModel &rcluster = cluster2;
#else
    const ClusterModel &lcluster = cluster1.col < cluster2.col ? cluster1 : cluster2;
    const ClusterModel &rcluster = cluster1.col < cluster2.col ? cluster2 : cluster1;
#endif

    ClusterModel new_cluster(lcluster);

    if (lcluster.row <= rcluster.row)
    {
        new_cluster.ncols += 1;
        for (const auto &hit : rcluster.hit_map)
        {
            int rel_pos_hit_col = hit.first + SIZE_QCORE_HORIZONTAL * (rcluster.col - lcluster.col); // Just move columns
            int rel_pos_hit_row = hit.second + SIZE_QCORE_VERTICAL * (rcluster.row - lcluster.row);  // NOTE the difference is here
            new_cluster.hit_map.push_back(std::make_pair(rel_pos_hit_col, rel_pos_hit_row));
            int nrows_calc = (rcluster.row + rcluster.nrows) - (lcluster.row + lcluster.nrows);
            if (nrows_calc >= 0)
            {
                new_cluster.nrows += nrows_calc;
            }
        }
    }
    else
    {
        new_cluster = ClusterModel(lcluster.col, rcluster.row);

        new_cluster.nrows = rcluster.nrows;

        int nrows_calc = (lcluster.row + lcluster.nrows) - (rcluster.row + rcluster.nrows);
        if (nrows_calc >= 0)
        {
            new_cluster.nrows += nrows_calc;
        }
        new_cluster.ncols = lcluster.ncols + 1;

        // Add hits now
        for (const auto &hit : lcluster.hit_map)
        {
            int rel_pos_hit_col = hit.first;
            int rel_pos_hit_row = hit.second + SIZE_QCORE_VERTICAL * (lcluster.row - rcluster.row); // NOTE the difference is here
            new_cluster.hit_map.push_back(std::make_pair(rel_pos_hit_col, rel_pos_hit_row));
        }

        for (const auto &hit : rcluster.hit_map)
        {
            int rel_pos_hit_col = hit.first + SIZE_QCORE_HORIZONTAL * (rcluster.col - lcluster.col); // Just move columns
            int rel_pos_hit_row = hit.second;
            new_cluster.hit_map.push_back(std::make_pair(rel_pos_hit_col, rel_pos_hit_row));
        }
    }

    // Return
    return new_cluster;
}

void ColMergerModel::_MergeColumns()
{
    // Time for processing
    std::vector<ClusterModel> updated_processor_buf;
    std::vector<int> cluster2_was_used(column_buf_2.size(), 0);

    for (size_t i = 0; i < column_buf_1.size(); i++)
    {
        const auto &cluster1 = column_buf_1[i];

        int cluster1_was_used = 0;
        int cluster2_id = 0;

        for (size_t j = 0; j < column_buf_2.size(); j++)
        {
            const auto &cluster2 = column_buf_2[j];

            // Just in case - this should not happen
            if (cluster1.col + cluster1.ncols != cluster2.col)
            {
#ifdef THROW_ERROR
                throw std::runtime_error("in merge columns, not adjacent columns");
#else
                std::cerr << "ERROR (COL_MERGER): in merge columns, not adjacent columns" << std::endl;
                not_adjacent++;
#endif
            }
            else
            {
                int is_touching = 0;

                // Relative position of row
                if (cluster1.row <= cluster2.row)
                {
                    for (const auto &hit1 : cluster1.hit_map)
                    {
                        for (const auto &hit2 : cluster2.hit_map)
                        {
                            int rel_pos_hit1_col = hit1.first;
                            int rel_pos_hit2_col = hit2.first + SIZE_QCORE_HORIZONTAL * cluster1.ncols; // Just move columns
                            int rel_pos_hit1_row = hit1.second;
                            int rel_pos_hit2_row = hit2.second + SIZE_QCORE_VERTICAL * (cluster2.row - cluster1.row); // NOTE the difference is here
                            if (abs(rel_pos_hit2_col - rel_pos_hit1_col) <= 1 && abs(rel_pos_hit2_row - rel_pos_hit1_row) <= 1)
                            {
                                is_touching = 1;
                                break;
                            }
                        }
                        if (is_touching == 1)
                        {
                            break;
                        }
                    }
                }
                else
                {
                    for (const auto &hit1 : cluster1.hit_map)
                    {
                        for (const auto &hit2 : cluster2.hit_map)
                        {
                            int rel_pos_hit1_col = hit1.first;
                            int rel_pos_hit2_col = hit2.first + SIZE_QCORE_HORIZONTAL * cluster1.ncols;               // Just move columns
                            int rel_pos_hit1_row = hit1.second + SIZE_QCORE_VERTICAL * (cluster1.row - cluster2.row); // NOTE the difference is here
                            int rel_pos_hit2_row = hit2.second;
                            if (abs(rel_pos_hit2_col - rel_pos_hit1_col) <= 1 && abs(rel_pos_hit2_row - rel_pos_hit1_row) <= 1)
                            {
                                is_touching = 1;
                                break;
                            }
                        }
                        if (is_touching == 1)
                        {
                            break;
                        }
                    }
                }

                // Now check if it's touching
                if (is_touching == 1)
                {
                    // Check reuse
                    if (cluster1_was_used > 0)
                    {

                        if (cluster2_was_used[cluster2_id] > 0)
                        {
#ifdef THROW_ERROR
                            throw std::runtime_error("ERROR (COL_MERGER): Impossible both clusters were already used");
#else
                            std::cerr << "ERROR (COL_MERGER): Impossible both clusters were already used" << std::endl;
                            already_used++;
#endif
                        }
                        else
                        {
                            cluster2_was_used[cluster2_id]++;
                        }

                        // Find clusters to merge
                        int found_id = -1;
                        for (size_t cluster_id = 0; cluster_id < updated_processor_buf.size(); cluster_id++)
                        {
                            const ClusterModel &current_cluster = updated_processor_buf[cluster_id];
                            for (const auto &hit : current_cluster.hit_map)
                            {
                                if (compare_vertexes(current_cluster.get_hit_global_coordinate(hit), cluster1.get_hit_global_coordinate(cluster1.hit_map[0])))
                                {
                                    found_id = static_cast<int>(cluster_id);
                                    break;
                                }
                            }
                            if (found_id >= 0)
                            {
                                break;
                            }
                        }

                        if (found_id < 0)
                        {
#ifdef THROW_ERROR
                            throw std::runtime_error("ERROR (COL_MERGER): Did not find matching cluster");
#else
                            std::cerr << "ERROR (COL_MERGER): Did not find matching cluster" << std::endl;
                            not_found++;
#endif
                        }
                        else
                        {
                            const ClusterModel &cluster_to_merge_1 = updated_processor_buf[found_id];
                            const ClusterModel &cluster_to_merge_2 = cluster2;
                            updated_processor_buf.push_back(_merge_two_clusters(cluster_to_merge_1, cluster_to_merge_2));
                        }
                    }
                    else if (cluster2_was_used[cluster2_id] > 0)
                    {

                        if (cluster1_was_used > 0)
                        {
#ifdef THROW_ERROR
                            throw std::runtime_error("ERROR (COL_MERGER): Impossible both clusters were already used");
#else
                            std::cerr << "ERROR (COL_MERGER): Impossible both clusters were already used" << std::endl;
                            already_used++;
#endif
                        }
                        else
                        {
                            cluster1_was_used++;
                        }

                        // Find clusters to merge
                        int found_id = -1;
                        for (size_t cluster_id = 0; cluster_id < updated_processor_buf.size(); cluster_id++)
                        {
                            const ClusterModel &current_cluster = updated_processor_buf[cluster_id];
                            for (const auto &hit : current_cluster.hit_map)
                            {
                                if (compare_vertexes(current_cluster.get_hit_global_coordinate(hit), cluster2.get_hit_global_coordinate(cluster2.hit_map[0])))
                                {
                                    found_id = static_cast<int>(cluster_id);
                                    break;
                                }
                            }
                            if (found_id >= 0)
                            {
                                break;
                            }
                        }

                        if (found_id < 0)
                        {
                            std::cerr << "Needed: " << std::endl;
                            std::cerr << cluster2.to_string() << std::endl;
                            std::cerr << "List: " << std::endl;
                            for (const auto &cluster : updated_processor_buf)
                            {
                                std::cerr << cluster.to_string() << std::endl;
                            }

                            throw std::runtime_error("Did not find matching cluster");
                        }
                        else
                        {
                            const ClusterModel &cluster_to_merge_2 = updated_processor_buf[found_id];
                            const ClusterModel &cluster_to_merge_1 = cluster1;
                            updated_processor_buf.push_back(_merge_two_clusters(cluster_to_merge_1, cluster_to_merge_2));
                        }
                    }
                    else
                    {
                        // Increment cluster use
                        cluster1_was_used++;
                        cluster2_was_used[cluster2_id]++;

                        const ClusterModel &cluster_to_merge_1 = cluster1;
                        const ClusterModel &cluster_to_merge_2 = cluster2;
                        updated_processor_buf.push_back(_merge_two_clusters(cluster_to_merge_1, cluster_to_merge_2));
                    }
                }

                // Cluster2 id increment
                cluster2_id++;
            }
        }

        // Check if cluster1 was used at least once
        if (cluster1_was_used == 0)
        {
            processed_cluster_list.push_back(cluster1);
        }
    }

    // The iteration is over, check if any of the cluster2 were not used
    for (size_t cluster2_id = 0; cluster2_id < column_buf_2.size(); cluster2_id++)
    {
        if (cluster2_was_used[cluster2_id] == 0)
        {
            if (column_buf_2[cluster2_id].is_touching_next_col())
            {
                updated_processor_buf.push_back(column_buf_2[cluster2_id]);
            }
            else
            {
                processed_cluster_list.push_back(column_buf_2[cluster2_id]);
            }
        }
    }

    // Filter clusters
    column_buf_1.clear();
    for (const auto &cluster : updated_processor_buf)
    {
        if (cluster.is_touching_next_col())
        {
            column_buf_1.push_back(cluster);
        }
        else
        {
            processed_cluster_list.push_back(cluster);
        }
    }

    // And now check if we have something to process still
    if (!column_buf_1.empty())
    {
        _current_state = "FillColumn2";
    }
    else
    {
        _current_state = "FillColumn1";
    }
}

void ColMergerModel::_Done() const
{
    // No way to get out of this state
}

unsigned int ColMergerModel::get_not_adjacent() const
{
    return this->not_adjacent;
}

unsigned int ColMergerModel::get_already_used() const
{
    return this->already_used;
}

unsigned int ColMergerModel::get_not_found() const
{
    return this->not_found;
}

unsigned int ColMergerModel::get_rhs_must_be_bigger() const
{
    return this->rhs_must_be_bigger;
}