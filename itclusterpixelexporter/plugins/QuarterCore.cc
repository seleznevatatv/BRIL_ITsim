/**
 * @file QuarterCore.cc
 * @author Max Bensink (maxbensink@outlook.com)
 * @version 1.0
 * @date 2023-08-01
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <iostream>

#include "ClusterModel.h"

QuarterCore::QuarterCore() : col(0), row(0), is_last(0), is_neighbour(0), is_last_in_event(0) {}
QuarterCore::QuarterCore(int qcol, int qrow) : col(qcol), row(qrow), is_last(0), is_neighbour(0), is_last_in_event(0) {}

void QuarterCore::add_hit(int col, int row)
{
    if (col >= 0 && col < SIZE_QCORE_HORIZONTAL && row >= 0 && row < SIZE_QCORE_VERTICAL)
    {
        hit_map.push_back(std::make_pair(col, row));
    }
    else
    {
        throw std::runtime_error("Error col row out of range");
    }
}

std::vector<std::vector<int>> QuarterCore::get_unsparsified_hit_map() const
{
    std::vector<std::vector<int>> unpacked_hitmap(SIZE_QCORE_HORIZONTAL, std::vector<int>(SIZE_QCORE_VERTICAL, 0));
    for (const auto &hit : hit_map)
    {
        unpacked_hitmap[hit.first][hit.second] = 1;
    }
    return unpacked_hitmap;
}

int QuarterCore::get_hit_4x4(int id) const
{
    int rem8 = id % 8;
    int row = (id - rem8) / 4 + (id % 2);
    int col = (rem8 - rem8 % 2) / 2;

    auto point = std::make_pair(col, row);

    return std::find(hit_map.begin(), hit_map.end(), point) != hit_map.end() ? 1 : 0;
}

std::pair<int, int> QuarterCore::get_binary_tree() const
{
    int binary_tree = 0;
    int binary_tree_length = 0;
    std::vector<int> pair_hitor(8, 0);
    std::vector<int> quarter_hitor(4, 0);
    std::vector<int> half_hitor(2, 0);
    int pair_id = 0; // start from right bottom corner

    while (pair_id < 8)
    {
        // misc
        int quarter_id = (pair_id - (pair_id % 2)) / 2;
        int half_id = (quarter_id - (quarter_id % 2)) / 2;

        // get info
        int hit_bottom = get_hit_4x4((8 - pair_id) * 2 - 1);
        int hit_top = get_hit_4x4((8 - pair_id) * 2 - 2);

        // assign hit pair or
        pair_hitor[pair_id] = hit_top || hit_bottom;

        // encode if something is there
        if (!hit_top && hit_bottom)
        { // 01
            // will be actually coded as 0
            binary_tree |= (0x0 << binary_tree_length);
            binary_tree_length += 1;
        }
        else if (hit_top && !hit_bottom)
        { // 10
            binary_tree |= (0x2 << binary_tree_length);
            binary_tree_length += 2;
        }
        else if (hit_top && hit_bottom)
        { // 11
            binary_tree |= (0x3 << binary_tree_length);
            binary_tree_length += 2;
        }

        // every second pair - send the s3 info
        if (pair_id % 2 == 1)
        {
            quarter_hitor[quarter_id] = pair_hitor[pair_id] || pair_hitor[pair_id - 1];

            // within this context - every second quarter
            if (quarter_id % 2 == 1)
            {
                half_hitor[half_id] = quarter_hitor[quarter_id] || quarter_hitor[quarter_id - 1];

                // s3 right stored here
                if (!pair_hitor[pair_id - 2] && pair_hitor[pair_id - 3])
                { // 01
                    // will be actually coded as 0
                    binary_tree |= (0x0 << binary_tree_length);
                    binary_tree_length += 1;
                }
                else if (pair_hitor[pair_id - 2] && !pair_hitor[pair_id - 3])
                { // 10
                    binary_tree |= (0x2 << binary_tree_length);
                    binary_tree_length += 2;
                }
                else if (pair_hitor[pair_id - 2] && pair_hitor[pair_id - 3])
                { // 11
                    binary_tree |= (0x3 << binary_tree_length);
                    binary_tree_length += 2;
                }

                // s3 left stored here
                if (!pair_hitor[pair_id] && pair_hitor[pair_id - 1])
                { // 01
                    // will be actually coded as 0
                    binary_tree |= (0x0 << binary_tree_length);
                    binary_tree_length += 1;
                }
                else if (pair_hitor[pair_id] && !pair_hitor[pair_id - 1])
                { // 10
                    binary_tree |= (0x2 << binary_tree_length);
                    binary_tree_length += 2;
                }
                else if (pair_hitor[pair_id] && pair_hitor[pair_id - 1])
                { // 11
                    binary_tree |= (0x3 << binary_tree_length);
                    binary_tree_length += 2;
                }

                // s2 stored here
                if (!quarter_hitor[quarter_id] && quarter_hitor[quarter_id - 1])
                { // 01
                    // will be actually coded as 0
                    binary_tree |= (0x0 << binary_tree_length);
                    binary_tree_length += 1;
                }
                else if (quarter_hitor[quarter_id] && !quarter_hitor[quarter_id - 1])
                { // 10
                    binary_tree |= (0x2 << binary_tree_length);
                    binary_tree_length += 2;
                }
                else if (quarter_hitor[quarter_id] && quarter_hitor[quarter_id - 1])
                { // 11
                    binary_tree |= (0x3 << binary_tree_length);
                    binary_tree_length += 2;
                }
            }
        }

        // now increment pair id
        pair_id++;
    }

    // now we can add s1 info
    if (!half_hitor[1] && half_hitor[0])
    { // 01
        // will be actually coded as 0
        binary_tree |= (0x0 << binary_tree_length);
        binary_tree_length += 1;
    }
    else if (half_hitor[1] && !half_hitor[0])
    { // 10
        binary_tree |= (0x2 << binary_tree_length);
        binary_tree_length += 2;
    }
    else if (half_hitor[1] && half_hitor[0])
    { // 11
        binary_tree |= (0x3 << binary_tree_length);
        binary_tree_length += 2;
    }

    return std::make_pair(binary_tree, binary_tree_length);
}

std::string QuarterCore::to_string() const
{
    std::string temp = "\tqcol = " + std::to_string(col) + ", qrow = " + std::to_string(row);
    temp += "\n\tis_last_in_event = " + std::to_string(is_last_in_event) + ", is_last = " + std::to_string(is_last) + ", is_neighbour = " + std::to_string(is_neighbour);
    temp += "\n\thit map:\n";
    std::vector<std::vector<int>> unpacked_hitmap = get_unsparsified_hit_map();
    for (int row = 0; row < SIZE_QCORE_VERTICAL; ++row)
    {
        temp += "\t\t";
        for (int col = 0; col < SIZE_QCORE_HORIZONTAL; ++col)
        {
            temp += (std::to_string(unpacked_hitmap[col][row]) + " ");
        }
        temp += "\n";
    }
    temp += "\tbinarytree:\n";
    std::pair<int, int> binary_tree_info = get_binary_tree();
    temp += ("\t\t" + std::to_string(binary_tree_info.first) + " (Length: " + std::to_string(binary_tree_info.second) + ")\n\n");
    return temp;
}