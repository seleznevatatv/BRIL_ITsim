/**
 * @file DetermineMetadeta.cc
 * @author Max Bensink (maxbensink@outlook.com)
 * @version 1.0
 * @date 2023-08-01
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <vector>
#include <map>
#include <algorithm>

#include "ClusterModel.h"

typedef std::pair<int, int> cluster_key_t;

std::vector<QuarterCore> determine_metadeta(const std::vector<QuarterCore> &qcores)
{
    std::vector<cluster_key_t> keys;
    std::map<cluster_key_t, QuarterCore> qcore_mapping;
    std::map<int, std::vector<int>> pos_mapping;
    std::vector<QuarterCore> qcores_width_metadata;

    for (auto &qcore : qcores)
        pos_mapping[qcore.col].push_back(qcore.row);

    for (auto &col : pos_mapping)
        std::sort(col.second.begin(), col.second.end());

    for (const auto &col_data : pos_mapping)
    {
        auto &col = col_data.first;

        for (const auto &row : col_data.second)
        {
            keys.push_back(std::make_pair(col, row));
            qcore_mapping[std::make_pair(col, row)] = *std::find_if(qcores.begin(), qcores.end(), [&col, &row](const QuarterCore &qcore)
                                                                    { return col == qcore.col && row == qcore.row; });
        }
    }

    for (size_t i = 0; i < keys.size(); i++)
    {

        const auto &current_key = keys[i];
        auto &current_qcore = qcore_mapping[current_key];

        if (i == keys.size() - 1)
        {
            current_qcore.is_last = 1;
            current_qcore.is_last_in_event = 1;
        }
        else
        {
            current_qcore.is_last_in_event = 0;

            auto next_key = keys[i + 1];

            current_qcore.is_last = current_key.first == next_key.first ? 0 : 1;
        }

        if (i == 0)
            current_qcore.is_neighbour = 0;
        else
        {
            auto prev_key = keys[i - 1];

            current_qcore.is_neighbour = ((current_key.first == prev_key.first) && (current_key.second == (prev_key.second + 1))) ? 1 : 0;
        }

    }

    for (const auto &col_data : pos_mapping)
    {
        auto &col = col_data.first;

        for (const auto &row : col_data.second)
        {
            qcores_width_metadata.push_back(qcore_mapping[std::make_pair(col, row)]);
        }
    }

    return qcores_width_metadata;
}

bool compare_vertexes(std::pair<int, int> a, std::pair<int, int> b)
{
    return a.first == b.first && a.second == b.second;
}