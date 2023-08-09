/**
 * @file ClusterModel.h
 * @author your name (you@domain.com)
 * @version 0.1
 * @date 2023-08-01
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <iostream>
#include <vector>
#include <functional>
#include <unordered_map>

// #define THROW_ERROR

/// @brief The horizontal size of the QuarterCore.
constexpr int SIZE_QCORE_HORIZONTAL = 4;

/// @brief The vertical size of the QuarterCore.
constexpr int SIZE_QCORE_VERTICAL = 4;

/// @struct ClusterModel
/// @brief Represents a cluster in the QCore data.
/// @details A cluster is a collection of hits in the same column and adjacent rows.
struct ClusterModel
{
    /// @brief Constructs a ClusterModel with the given column and row.
    /// @param col The column of the cluster.
    /// @param row The row of the cluster.
    ClusterModel(int col, int row);

    /// @brief Adds a hit to the cluster with the specified local column and row.
    /// @param local_col The local column of the hit within the cluster.
    /// @param local_row The local row of the hit within the cluster.
    void add_hit(int local_col, int local_row);

    /// @brief Moves to the next column in the cluster.
    void next_col();

    /// @brief Moves to the next row in the cluster.
    void next_row();

    /// @brief Checks if the cluster is touching the previous column.
    /// @return True if the cluster is touching the previous column, false otherwise.
    bool is_touching_prev_col() const;

    /// @brief Checks if the cluster is touching the next column.
    /// @return True if the cluster is touching the next column, false otherwise.
    bool is_touching_next_col() const;

    /// @brief Checks if the cluster is touching the previous row.
    /// @return True if the cluster is touching the previous row, false otherwise.
    bool is_touching_prev_row() const;

    /// @brief Checks if the cluster is touching the next row.
    /// @return True if the cluster is touching the next row, false otherwise.
    bool is_touching_next_row() const;

    /// @brief Returns the global coordinate of the specified hit.
    /// @param hit The local coordinate of the hit within the cluster.
    /// @return The global coordinate (column, row) of the hit.
    std::pair<int, int> get_hit_global_coordinate(const std::pair<int, int> &hit) const;

    /// @brief Converts the cluster to a string representation.
    /// @return The string representation of the cluster.
    std::string to_string() const;

    std::vector<std::pair<int, int>> hit_map; ///< The list of hits in the cluster.
    int col;                                  ///< The column of the cluster.
    int row;                                  ///< The row of the cluster.
    int ncols;                                ///< The number of columns in the cluster.
    int nrows;                                ///< The number of rows in the cluster.
};

/// @struct QuarterCore
/// @brief Represents a QuarterCore in the QCore data.
struct QuarterCore
{
    /// @brief Default constructor for QuarterCore.
    QuarterCore();

    /// @brief Constructs a QuarterCore with the given column and row.
    /// @param qcol The column of the QuarterCore.
    /// @param qrow The row of the QuarterCore.
    QuarterCore(int qcol, int qrow);

    /// @brief Adds a hit to the QuarterCore with the specified column and row.
    /// @param col The column of the hit.
    /// @param row The row of the hit.
    void add_hit(int col, int row);

    /// @brief Returns the 4x4 hit at the specified index in the hit map.
    /// @param id The index of the hit in the hit map.
    /// @return The 4x4 hit value at the specified index.
    int get_hit_4x4(int id) const;

    /// @brief Returns the unsparsified hit map of the QuarterCore.
    /// @return The unsparsified hit map as a 2D vector of integers.
    std::vector<std::vector<int>> get_unsparsified_hit_map() const;

    /// @brief Returns the sparsified hit map of the QuarterCore.
    /// @return The sparsified hit map as a list of global coordinates (column, row).
    std::vector<std::pair<int, int>> get_sparsified_hit_map() const;

    /// @brief Returns the binary tree representation of the QuarterCore.
    /// @return The binary tree representation as a pair of integers.
    std::pair<int, int> get_binary_tree() const;

    /// @brief Converts the QuarterCore to a string representation.
    /// @return The string representation of the QuarterCore.
    std::string to_string() const;

    std::vector<std::pair<int, int>> hit_map; ///< The list of hits in the QuarterCore.
    int col;                                  ///< The column of the QuarterCore.
    int row;                                  ///< The row of the QuarterCore.
    int is_last;                              ///< Indicator for the last QuarterCore in the data.
    int is_neighbour;                         ///< Indicator for the neighboring QuarterCore.
    int is_last_in_event;                     ///< Indicator for the last QuarterCore in the event.
};

/// @brief Determines metadata for QuarterCores in the given list.
/// @param qcores The list of QuarterCores.
/// @return A vector of QuarterCores with metadata.
std::vector<QuarterCore> determine_metadeta(const std::vector<QuarterCore> &qcores);

/// @struct ProcessedQuarterCoreModel
/// @brief Represents a processed QuarterCore with metadata.
struct ProcessedQuarterCoreModel
{
    /// @brief Constructs a ProcessedQuarterCoreModel using the given QuarterCore.
    /// @param qcore The QuarterCore to be processed.
    ProcessedQuarterCoreModel(const QuarterCore &qcore);

    std::vector<std::vector<int>> hit_map, unpacked_hitmap; ///< The hit map and unpacked hit map of the processed QuarterCore.
    int col;                                                ///< The column of the processed QuarterCore.
    int row;                                                ///< The row of the processed QuarterCore.
    int is_last;                                            ///< Indicator for the last QuarterCore in the data.
    int is_neighbour;                                       ///< Indicator for the neighboring QuarterCore.
    int is_last_in_event;                                   ///< Indicator for the last QuarterCore in the event.
    int is_touching_prev_col;                               ///< Indicator for touching the previous column.
    int is_touching_next_col;                               ///< Indicator for touching the next column.
    int is_touching_prev_row;                               ///< Indicator for touching the previous row.
    int is_touching_next_row;                               ///< Indicator for touching the next row.
    std::vector<ClusterModel> cluster_list;                 ///< List of clusters in the processed QuarterCore.
};

/// @class DistributorModel
/// @brief Represents a model for distributing QuarterCores.
class DistributorModel
{
public:
    /// @brief Constructs a DistributorModel with the given QuarterCore buffer.
    /// @param qcore_buf The buffer of processed QuarterCores.
    DistributorModel(const std::vector<ProcessedQuarterCoreModel> &qcore_buf);

    /// @brief Runs the DistributorModel to distribute QuarterCores and clusters.
    /// @return A tuple containing the list of clusters and the list of processed QuarterCores.
    std::tuple<std::vector<ClusterModel>, std::vector<ProcessedQuarterCoreModel>> run();

private:
    std::vector<ProcessedQuarterCoreModel> qcore_buf; ///< The buffer of processed QuarterCores.
};

/// @class RowMergerModel
/// @brief Represents a model for merging rows of QuarterCores.
class RowMergerModel
{
public:
    /// @brief Constructs a RowMergerModel with the given buffer of processed QuarterCores.
    /// @param row_merger_buf The buffer of processed QuarterCores to be merged.
    RowMergerModel(const std::vector<ProcessedQuarterCoreModel> &row_merger_buf);

    /// @brief Runs the RowMergerModel to merge rows of QuarterCores and clusters.
    /// @return A tuple containing the list of clusters and the list of remaining processed QuarterCores.
    std::tuple<std::vector<ClusterModel>, std::vector<ClusterModel>> run();

    unsigned int get_processor_buffer_not_empty() const;
    unsigned int get_already_used() const;
    unsigned int get_not_found() const;

private:
    /// @brief Fills the buffer for processing.
    void _FillBuffer();

    /// @brief Appends the current QuarterCore to the processor buffer.
    void _AppendQcore();

    /// @brief Handles the "Done" state of the RowMergerModel.
    void _Done() const;

    /// @brief Distributes one cluster to the processor buffer or the Column Merger buffer.
    /// @param cluster The cluster to be distributed.
    void distribute_one_cluster(const ClusterModel &cluster);

    /// @brief Distributes all clusters from the processor buffer to the Column Merger buffer.
    void distribute_all_clusters();

    std::vector<ProcessedQuarterCoreModel> row_merger_buf;             ///< The buffer of processed QuarterCores to be merged.
    std::vector<ClusterModel> processor_buf;                           ///< The buffer for processing clusters.
    std::vector<ClusterModel> col_merger_buf;                          ///< The buffer for Column Merger clusters.
    std::vector<ClusterModel> processed_cluster_list;                  ///< The list of processed clusters.
    size_t qcore_id;                                                   ///< The current QuarterCore ID in the buffer.
    std::string _current_state;                                        ///< The current state of the RowMergerModel.
    std::unordered_map<std::string, std::function<void()>> _state_map; ///< The state map for RowMergerModel.

    unsigned int processor_buffer_not_empty;
    unsigned int already_used;
    unsigned int not_found;
};

/// @class ColMergerModel
/// @brief Represents a model for merging columns of clusters.
class ColMergerModel
{
public:
    /// @brief Constructs a ColMergerModel with the given buffer of clusters.
    /// @param col_merger_buf The buffer of clusters to be merged.
    ColMergerModel(const std::vector<ClusterModel> &col_merger_buf);

    /// @brief Runs the ColMergerModel to merge columns of clusters.
    /// @return The list of processed clusters after column merging.
    std::vector<ClusterModel> run();

    unsigned int get_not_adjacent() const;
    unsigned int get_already_used() const;
    unsigned int get_not_found() const;
    unsigned int get_rhs_must_be_bigger() const;

private:
    /// @brief Fills Column 1 with clusters.
    void _FillColumn1();

    /// @brief Empties Column 1 and processes its clusters.
    void _empty_column_1();

    /// @brief Fills Column 2 with clusters.
    void _FillColumn2();

    /// @brief Merges two clusters together.
    /// @param cluster1 The first cluster to merge.
    /// @param cluster2 The second cluster to merge.
    /// @return The merged cluster resulting from the merge.
    ClusterModel _merge_two_clusters(const ClusterModel &cluster1, const ClusterModel &cluster2);

    /// @brief Merges the clusters from Column 1 and Column 2.
    void _MergeColumns();

    /// @brief Handles the "Done" state of the ColMergerModel.
    void _Done() const;

    std::vector<ClusterModel> col_merger_buf;                          ///< The buffer of clusters to be merged.
    std::vector<ClusterModel> processed_cluster_list;                  ///< The list of processed clusters after column merging.
    std::vector<ClusterModel> column_buf_1;                            ///< The buffer for Column 1 clusters.
    std::vector<ClusterModel> column_buf_2;                            ///< The buffer for Column 2 clusters.
    std::string _current_state;                                        ///< The current state of the ColMergerModel.
    std::unordered_map<std::string, std::function<void()>> _state_map; ///< The state map for ColMergerModel.

    unsigned int not_adjacent;
    unsigned int already_used;
    unsigned int not_found;
    unsigned int rhs_must_be_bigger;
};

/**
 * @brief compair two vertices to check if they have the same value
 *
 * @param a first vertex
 * @param b second vertex
 * @return true
 * @return false
 */
bool compare_vertexes(std::pair<int, int> a, std::pair<int, int> b);
