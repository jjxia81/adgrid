#pragma once

#include <array>
#include <functional>
#include <span>

#include "indirect_value.hpp"
#include "../strong_type/strong_type.hpp"


namespace mtet {

using Scalar = double;
class MTetMeshImpl;

using VertexId =
    strong::type<uint64_t, struct VertexId_, strong::equality, strong::default_constructible>;
using TetId = strong::type<uint64_t, struct TetId_, strong::default_constructible>;
using EdgeId = strong::type<uint64_t, struct EdgeId_, strong::default_constructible>;

bool operator==(TetId t0, TetId t1);
bool operator==(EdgeId e0, EdgeId e1);

class MTetMesh
{
public:
    MTetMesh();
    ~MTetMesh();

    MTetMesh(const MTetMesh&) noexcept;
    MTetMesh& operator=(const MTetMesh&) noexcept;

    MTetMesh(MTetMesh&&) noexcept;
    MTetMesh& operator=(MTetMesh&&) noexcept;

public:
    VertexId add_vertex(Scalar x, Scalar y, Scalar z);
    TetId add_tet(VertexId v0, VertexId v1, VertexId v2, VertexId v3);
    void initialize_connectivity();

public:
    bool has_vertex(VertexId vertex_id) const;
    bool has_tet(TetId tet_id) const;
    bool has_edge(EdgeId edge_id) const;

    std::span<Scalar, 3> get_vertex(VertexId vertex_id);
    std::span<VertexId, 4> get_tet(TetId tet_id);

    std::span<const Scalar, 3> get_vertex(VertexId vertex_id) const;
    std::span<const VertexId, 4> get_tet(TetId tet_id) const;

    EdgeId get_edge(TetId tet_id, uint8_t local_index) const;
    std::array<VertexId, 2> get_edge_vertices(EdgeId edge_id) const;

    /**
     * Get a tet that contains the given edge.
     */
    TetId get_edge_tet(EdgeId edge_id) const;

    /**
     * Get the mirror of a given tet with `tet_id` across its local face indexed by `local_index`.
     */
    TetId get_mirror(TetId tet_id, uint8_t local_index) const;

    size_t get_num_vertices() const;
    size_t get_num_tets() const;

public:
    /**
     * Split the edge of the given tet with the given local edge id.
     *
     * The split will insert a new vertex at the mid point of the specified edge.
     *
     * @param edge_id        The id of the edge to split.
     *
     * @return [vid, e0_id, e1_id]
     *         vid: The id of the new vertex.
     *         e0_id: The id of the first half of the split edge.
     *         e1_id: The id of the second half of the split edge.
     */
    std::tuple<VertexId, EdgeId, EdgeId> split_edge(EdgeId edge_id);

    /**
     * Split the edge of the given tet with the given local edge id.
     *
     * The split will insert a new vertex at the mid point of the specified edge.
     *
     * @param tet_id        The id of the tet to split.
     * @param local_edge_id The local edge id of the edge to split.
     *
     * @return [vid, e0_id, e1_id]
     *         vid: The id of the new vertex.
     *         e0_id: The id of the first half of the split edge.
     *         e1_id: The id of the second half of the split edge.
     *
     * Let the oriented tet be [v0, v1, v2, v3]. The local edges order are the following:
     *   0: [v0, v1],
     *   1: [v1, v2],
     *   2: [v2, v0],
     *   3: [v0, v3],
     *   4: [v1, v3],
     *   5: [v2, v3]
     *
     * The indices are illustrated in the ASCII drawing below.
     *
     *               v2
     *             ,/|`\
     *           ,/  |  `\
     *         ,2    |.   `1
     *       ,/       5     `\
     *     ,/         |       `\
     *   v0-------0-- |. -------v1
     *     `\.         |      ,/
     *        `\.      |    ,4
     *           `3.   |. ,/
     *              `\. |/
     *                 `v3
     */
    std::tuple<VertexId, EdgeId, EdgeId> split_edge(TetId tet_id, uint8_t local_edge_id);

public:
    void par_foreach_vertex(
        const std::function<void(VertexId, std::span<const Scalar, 3>)>& callback) const;
    void seq_foreach_vertex(
        const std::function<void(VertexId, std::span<const Scalar, 3>)>& callback) const;
    void par_foreach_tet(
        const std::function<void(TetId, std::span<const VertexId, 4>)>& callback) const;
    void seq_foreach_tet(
        const std::function<void(TetId, std::span<const VertexId, 4>)>& callback) const;
    void foreach_edge_in_tet(
        TetId tet_id,
        const std::function<void(EdgeId, VertexId, VertexId)>& callback);
    void foreach_tet_around_edge(EdgeId edge_id, const std::function<void(TetId)>& callback) const;

private:
    nonstd::indirect_value<MTetMeshImpl> m_impl;
};


} // namespace mtet
