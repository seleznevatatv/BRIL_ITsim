/**
 * @file RowMergerModel.cc
 * @author Max Bensink (maxbensink@outlook.com)
 * @version 1.0
 * @date 2023-08-01
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "ClusterModel.h"

RowMergerModel::RowMergerModel(const std::vector<ProcessedQuarterCoreModel> &row_merger_buf) : row_merger_buf(row_merger_buf), qcore_id(0), processor_buffer_not_empty(0), already_used(0), not_found(0)
{
    _state_map["FillBuffer"] = [this]()
    { _FillBuffer(); };
    _state_map["AppendQcore"] = [this]()
    { _AppendQcore(); };
    _state_map["Done"] = [this]()
    { _Done(); };
}

std::tuple<std::vector<ClusterModel>, std::vector<ClusterModel>> RowMergerModel::run()
{
    _current_state = "FillBuffer";
    processor_buf.clear();
    col_merger_buf.clear();
    processed_cluster_list.clear();
    qcore_id = 0;

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
    return std::make_tuple(processed_cluster_list, col_merger_buf);
}

void RowMergerModel::_FillBuffer()
{
    if (processor_buf.empty())
    {
        if (qcore_id < row_merger_buf.size())
        {
            processor_buf = row_merger_buf[qcore_id].cluster_list;
            qcore_id++;
            _current_state = "AppendQcore";
        }
        else
        {
            _current_state = "Done";
        }
    }
    else
    {
// Error: Processor buffer should be empty at this point
// Handle the error as per your requirement (e.g., throw an exception or log an error and exit the program).
// For demonstration purposes, we are just printing an error message and exiting the program.
#ifdef THROW_ERROR

        throw std::runtime_error("Error (ROW_MERGER): Processor buffer is not empty when filling buffer!");
#else
        std::cerr << "ERROR (ROW_MERGER): Processor buffer is not empty when filling buffer!" << std::endl;
        processor_buffer_not_empty++;
#endif
    }
}

void RowMergerModel::_AppendQcore()
{
    // Check if it's not the last cluster
    if (qcore_id < row_merger_buf.size())
    {
        const ProcessedQuarterCoreModel &qcore_to_append = row_merger_buf[qcore_id];
        const std::vector<ClusterModel> &cluster2_list = qcore_to_append.cluster_list;

        // We are doing all the processing within a single column
        if (qcore_to_append.col == processor_buf[0].col)
        {
            // Create a new processor buffer
            std::vector<ClusterModel> updated_processor_buf;
            std::vector<int> cluster2_was_used(cluster2_list.size(), 0);

            for (size_t i = 0; i < processor_buf.size(); i++)
            {
                const auto &cluster1 = processor_buf[i];

                int cluster1_was_used = 0;

                // Check if this cluster is big enough to match with the appended qcore
                if (cluster1.row + cluster1.nrows == qcore_to_append.row)
                {
                    // Iterate over cluster2 list
                    for (size_t cluster2_id = 0; cluster2_id < cluster2_list.size(); cluster2_id++)
                    {
                        // Finally time to check all matching clusters
                        int is_touching = 0;
                        for (const auto &hit1 : cluster1.hit_map)
                        {
                            for (const auto &hit2 : cluster2_list[cluster2_id].hit_map)
                            {
                                if (std::abs(hit2.first - hit1.first) <= 1 &&
                                    std::abs(hit2.second + SIZE_QCORE_VERTICAL * cluster1.nrows - hit1.second) <= 1)
                                {
                                    is_touching = 1;
                                    break;
                                }
                            }
                            if (is_touching)
                            {
                                break;
                            }
                        }

                        // Now check if it was used
                        if (is_touching)
                        {
                            std::pair<int, int> prcl, cl1;

                            if (cluster1_was_used > 0)
                            {
                                if (cluster2_was_used[cluster2_id] > 0)
                                {
// Error: Impossible that cluster2 was used together with cluster1
#ifdef THROW_ERROR
                                    throw std::runtime_error("ERROR: Impossible cluster2 was used together with cluster1");
#else
                                    std::cerr << "ERROR: Impossible cluster2 was used together with cluster1" << std::endl;
                                    already_used++;
#endif
                                }
                                else
                                {
                                    cluster2_was_used[cluster2_id]++;
                                }

                                int found_id = -1;
                                for (size_t processed_cluster_id = 0; processed_cluster_id < updated_processor_buf.size(); processed_cluster_id++)
                                {
                                    const ClusterModel &processed_cluster = updated_processor_buf[processed_cluster_id];
                                    for (const auto &hit : processed_cluster.hit_map)
                                    {

                                        prcl = processed_cluster.get_hit_global_coordinate(hit);
                                        cl1 = cluster1.get_hit_global_coordinate(cluster1.hit_map[0]);

                                        if (compare_vertexes(processed_cluster.get_hit_global_coordinate(hit), cluster1.get_hit_global_coordinate(cluster1.hit_map[0])))
                                        {
                                            found_id = static_cast<int>(processed_cluster_id);
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
                                    // Error: Did not find used cluster
#ifdef THROW_ERROR
                                    throw std::runtime_error("ERROR: Did not find used cluster");
#else
                                    std::cerr << "ERROR (ROW_MERGER): Did not find used cluster" << std::endl;
#endif
                                }
                                else
                                {
                                    ClusterModel current_cluster = updated_processor_buf[static_cast<size_t>(found_id)];
                                    for (const auto &hit_2 : cluster2_list[cluster2_id].hit_map)
                                    {
                                        current_cluster.hit_map.push_back(std::make_pair(hit_2.first, hit_2.second + SIZE_QCORE_VERTICAL * (current_cluster.nrows - 1)));
                                    }
                                    updated_processor_buf[static_cast<size_t>(found_id)] = current_cluster;
                                }
                            }
                            else if (cluster2_was_used[cluster2_id] > 0)
                            {
                                if (cluster1_was_used > 0)
                                {
// Error: Impossible that cluster2 was used together with cluster1
#ifdef THROW_ERROR
                                    throw std::runtime_error("ERROR (ROW_MERGER): Impossible cluster2 was used together with cluster1");
#else
                                    std::cerr << "ERROR (ROW_MERGER): Impossible cluster2 was used together with cluster1" << std::endl;
#endif
                                }
                                else
                                {
                                    cluster1_was_used++;
                                }

                                // std::cout << "NOTE: Cluster2 already used before - appending to it" << std::endl;

                                int found_id = 0;
                                for (size_t processed_cluster_id = 0; processed_cluster_id < updated_processor_buf.size(); processed_cluster_id++)
                                {
                                    const ClusterModel &processed_cluster = updated_processor_buf[processed_cluster_id];
                                    for (const auto &hit : processed_cluster.hit_map)
                                    {
                                        if (compare_vertexes(processed_cluster.get_hit_global_coordinate(hit), cluster2_list[cluster2_id].get_hit_global_coordinate(cluster2_list[cluster2_id].hit_map[0])))
                                        {
                                            found_id = static_cast<int>(processed_cluster_id);
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
// Error: Did not find used cluster
#ifdef THROW_ERROR
                                    throw std::runtime_error("ERROR: Did not find used cluster");
#else
                                    std::cerr << "ERROR (ROW_MERGER): Impossible cluster2 was used together with cluster1" << std::endl;
#endif
                                }
                                else
                                {
                                    ClusterModel current_cluster = updated_processor_buf[static_cast<size_t>(found_id)];
                                    for (const auto &hit_1 : cluster1.hit_map)
                                    {
                                        current_cluster.hit_map.push_back(std::make_pair(hit_1.first, hit_1.second));
                                    }
                                    updated_processor_buf[static_cast<size_t>(found_id)] = current_cluster;
                                }
                            }
                            else
                            {
                                ClusterModel new_cluster = cluster1;
                                new_cluster.next_row();
                                for (const auto &hit : cluster2_list[cluster2_id].hit_map)
                                {
                                    new_cluster.add_hit(hit.first, hit.second + SIZE_QCORE_VERTICAL * (new_cluster.nrows - 1));
                                }
                                updated_processor_buf.push_back(new_cluster);

                                cluster1_was_used++;
                                cluster2_was_used[cluster2_id]++;
                            }
                        }

                        // Increment cluster2 id
                        cluster2_id++;
                    }

                    // Distribute cluster1 if not used
                    if (cluster1_was_used == 0)
                    {
                        distribute_one_cluster(cluster1);
                    }
                }
                else
                {
                    // This cluster1 goes away
                    distribute_one_cluster(cluster1);
                }
            }

            // Cluster2 are appended to the processor buffer
            for (size_t cluster2_id = 0; cluster2_id < cluster2_list.size(); cluster2_id++)
            {
                if (cluster2_was_used[cluster2_id] == 0)
                {
                    if (cluster2_list[cluster2_id].is_touching_next_row())
                    {
                        updated_processor_buf.push_back(cluster2_list[cluster2_id]);
                    }
                    else
                    {
                        distribute_one_cluster(cluster2_list[cluster2_id]);
                    }
                }
            }

            // Assign for processing
            processor_buf.clear();
            for (const auto &cluster : updated_processor_buf)
            {
                if (cluster.is_touching_next_row())
                {
                    processor_buf.push_back(cluster);
                }
                else
                {
                    distribute_one_cluster(cluster);
                }
            }

            // Increment qcore id
            qcore_id++;

            // Check if we need to change state
            if (processor_buf.empty())
            {
                _current_state = "FillBuffer";
            }
        }
        else
        {
            // In this case, we have to empty all clusters in the buffer and start processing this new column
            // State remains the same
            distribute_all_clusters();
            processor_buf = cluster2_list;
            qcore_id++;
        }
    }
    else
    {
        // Seems that we are done
        // All done
        distribute_all_clusters();
        _current_state = "Done";
    }
}

void RowMergerModel::distribute_one_cluster(const ClusterModel &cluster)
{
    if (cluster.is_touching_prev_col() || cluster.is_touching_next_col())
    {
        col_merger_buf.push_back(cluster);
    }
    else
    {
        processed_cluster_list.push_back(cluster);
    }
}

void RowMergerModel::distribute_all_clusters()
{
    for (const auto &cluster : processor_buf)
    {
        distribute_one_cluster(cluster);
    }
    processor_buf.clear();
}

void RowMergerModel::_Done() const
{
    // No way to get out of this state
    // The function does nothing here, as it is the last state.
}

unsigned int RowMergerModel::get_processor_buffer_not_empty() const
{
    return this->processor_buffer_not_empty;
}

unsigned int RowMergerModel::get_already_used() const
{
    return this->already_used;
}

unsigned int RowMergerModel::get_not_found() const
{
    return this->not_found;
}