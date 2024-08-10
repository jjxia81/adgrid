//
//  mesh_refine.cpp
//  adaptive_mesh_refinement
//
//  Created by Yiwen Ju on 8/4/24.
//

#include "grid_refine.h"

/// hash for mounting a boolean that represents the activeness to a tet
/// since the tetid isn't const during the process, mount the boolean using vertexids of 4 corners.
uint64_t vertexHash(std::span<VertexId, 4>& x)
{
    ankerl::unordered_dense::hash<uint64_t> hash_fn;
    return hash_fn(value_of(x[0])) + hash_fn(value_of(x[1])) + hash_fn(value_of(x[2])) + hash_fn(value_of(x[3]));
}

bool gridRefine(
                const int mode,
                const bool curve_network,
                const double threshold,
                const double alpha,
                const int max_elements,
                const size_t funcNum,
                const std::function<llvm_vecsmall::SmallVector<Eigen::RowVector4d, 20>(std::span<const Scalar, 3>, size_t)> func,
                const std::function<std::pair<std::array<double, 2>, llvm_vecsmall::SmallVector<int, 20>>(llvm_vecsmall::SmallVector<std::array<double, 2>, 20>)> csg_func,
                const bool discretize_later,
                mtet::MTetMesh &grid,
                tet_metric &metric_list,
                std::array<double, timer_amount> profileTimer
                )
{
    /// Precomputing active multiples' indices:
    init_multi(funcNum, mode);
    
    /// Tet Metric
    int sub_call_two = 0;
    int sub_call_three = 0;

    /// initialize vertex map: vertex index -> {{f_i, gx, gy, gz} | for all f_i in the function}
    
    using IndexMap = ankerl::unordered_dense::map<uint64_t, llvm_vecsmall::SmallVector<Eigen::RowVector4d, 20>>;
    IndexMap vertex_func_grad_map;
    vertex_func_grad_map.reserve(grid.get_num_vertices());

    ///initialize activeness map: four vertexids (v0, v1, v2, v3) -> hash(v0, v1, v2, v3) -> active boolean
    using activeMap = ankerl::unordered_dense::map<uint64_t, bool>;
    activeMap vertex_active_map;
    vertex_active_map.reserve(grid.get_num_tets());

    grid.seq_foreach_vertex([&](VertexId vid, std::span<const Scalar, 3> data)
                            {vertex_func_grad_map[value_of(vid)] = func(data, funcNum);});

    auto comp = [](std::pair<mtet::Scalar, mtet::EdgeId> e0,
                   std::pair<mtet::Scalar, mtet::EdgeId> e1)
    { return e0.first < e1.first; };
    std::vector<std::pair<mtet::Scalar, mtet::EdgeId>> Q;

    Eigen::Matrix<double, 4, 3> pts;
    std::array<llvm_vecsmall::SmallVector<Eigen::RowVector4d, 20>,4> tet_info;
    for (size_t i = 0; i < tet_info.size(); i++){
        tet_info[i].resize(funcNum);
    }
    double activeTet = 0;
    auto push_longest_edge = [&](mtet::TetId tid)
    {
        std::span<VertexId, 4> vs = grid.get_tet(tid);
        {
            //Timer eval_timer(evaluation, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            for (int i = 0; i < 4; ++i)
            {
                auto vid = vs[i];
                auto coords = grid.get_vertex(vid);
                pts.row(i) = Eigen::RowVector3d({coords[0], coords[1], coords[2]});
                llvm_vecsmall::SmallVector<std::array<double, 4>, 20> func_gradList(funcNum);
                
                llvm_vecsmall::SmallVector<Eigen::RowVector4d, 20> func_info(funcNum);
                std::array<double, 4> func_grad;
                if (!vertex_func_grad_map.contains(value_of(vid))) {
                    func_info = func(coords, funcNum);
                    vertex_func_grad_map[value_of(vid)] = func_info;
                }
                else {
                    //func_gradList = vertex_func_grad_map[value_of(vid)];
                    func_info = vertex_func_grad_map[value_of(vid)];
                }
                tet_info[i] = func_info;
            }
            //eval_timer.Stop();
        }
        bool isActive = 0;
        bool subResult;
        {
            //Timer sub_timer(subdivision, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            switch (mode){
                case IA:
                    subResult = critIA(pts, tet_info, funcNum, threshold, curve_network, isActive, sub_call_two, sub_call_three);
                    break;
                case MI:
                    subResult = critMI(pts, tet_info, funcNum, threshold, curve_network, isActive, sub_call_two, sub_call_three);
                    break;
                case CSG:
                    subResult = critCSG(pts, tet_info, funcNum, csg_func, threshold, curve_network, isActive, sub_call_two, sub_call_three);
                    break;
                default:
                    throw std::runtime_error("no implicit complexes specified");
                    return false;
            }
            //sub_timer.Stop();
        }
        vertex_active_map[vertexHash(vs)] = isActive;
        //Timer eval_timer(evaluation, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
        if (subResult)
        {
            mtet::EdgeId longest_edge;
            mtet::Scalar longest_edge_length = 0;
            grid.foreach_edge_in_tet(tid, [&](mtet::EdgeId eid, mtet::VertexId v0, mtet::VertexId v1)
                                     {
                auto p0 = grid.get_vertex(v0);
                auto p1 = grid.get_vertex(v1);
                mtet::Scalar l = (p0[0] - p1[0]) * (p0[0] - p1[0]) + (p0[1] - p1[1]) * (p0[1] - p1[1]) +
                (p0[2] - p1[2]) * (p0[2] - p1[2]);
                if (l > longest_edge_length) {
                    longest_edge_length = l;
                    longest_edge = eid;
                } });
            Q.emplace_back(longest_edge_length, longest_edge);
            //eval_timer.Stop();
            return true;
        }
        //eval_timer.Stop();
        return false;
    };


    {
        Timer timer(total_time, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
        
        // Initialize priority queue.
        grid.seq_foreach_tet([&](mtet::TetId tid, [[maybe_unused]] std::span<const mtet::VertexId, 4> vs)
                             { push_longest_edge(tid); });
        std::make_heap(Q.begin(), Q.end(), comp);
        
        // Keep splitting the longest edge
        while (!Q.empty())
        {
            std::pop_heap(Q.begin(), Q.end(), comp);
            auto [edge_length, eid] = Q.back();
            if (!grid.has_edge(eid)){
                Q.pop_back();
                continue;
            }
            //implement alpha value:
            mtet::Scalar comp_edge_length = alpha * edge_length;
            bool addedActive = false;
            grid.foreach_tet_around_edge(eid,[&](mtet::TetId tid){
                std::span<VertexId, 4> vs = grid.get_tet(tid);
                if(vertex_active_map.contains(vertexHash(vs))){
                    if (vertex_active_map[vertexHash(vs)]){
                        mtet::EdgeId longest_edge;
                        mtet::Scalar longest_edge_length = 0;
                        grid.foreach_edge_in_tet(tid, [&](mtet::EdgeId eid_active, mtet::VertexId v0, mtet::VertexId v1)
                                                 {
                            auto p0 = grid.get_vertex(v0);
                            auto p1 = grid.get_vertex(v1);
                            mtet::Scalar l = (p0[0] - p1[0]) * (p0[0] - p1[0]) + (p0[1] - p1[1]) * (p0[1] - p1[1]) +
                            (p0[2] - p1[2]) * (p0[2] - p1[2]);
                            if (l > longest_edge_length) {
                                longest_edge_length = l;
                                longest_edge = eid_active;
                            }
                        });
                        if (longest_edge_length > comp_edge_length) {
                            Q.emplace_back(longest_edge_length, longest_edge);
                            addedActive = true;
                        }
                    }
                }
            });
            if(addedActive){
                std::push_heap(Q.begin(), Q.end(), comp);
                continue;
            }
            Q.pop_back();
            std::array<VertexId, 2> vs_old = grid.get_edge_vertices(eid);
            //Timer split_timer(splitting, [&](auto profileResult){profileTimer = combine_timer(profileTimer, profileResult);});
            auto [vid, eid0, eid1] = grid.split_edge(eid);
            //split_timer.Stop();
            //std::cout << "Number of elements: " << mesh.get_num_tets() << std::endl;
            if (grid.get_num_tets() > max_elements) {
                break;
            }
            grid.foreach_tet_around_edge(eid0, [&](mtet::TetId tid)
                                         {
                if (push_longest_edge(tid)) {
                    std::push_heap(Q.begin(), Q.end(), comp);
                } });
            grid.foreach_tet_around_edge(eid1, [&](mtet::TetId tid)
                                         {
                if (push_longest_edge(tid)) {
                    std::push_heap(Q.begin(), Q.end(), comp);
                } });
        }
        timer.Stop();
    }
    
    std::vector<mtet::TetId> activeTetId;
    grid.seq_foreach_tet([&](mtet::TetId tid, std::span<const VertexId, 4> data) {
        std::span<VertexId, 4> vs = grid.get_tet(tid);
        std::array<std::valarray<double>,4> vallPoints;
        for (int i = 0; i < 4; i++){
            vallPoints[i] = {0.0,0.0,0.0};
        }
        for (int i = 0; i < 4; i++){
            VertexId vid = vs[i];
            std::span<Scalar, 3> coords = grid.get_vertex(vid);
            vallPoints[i][0] = coords[0];
            vallPoints[i][1] = coords[1];
            vallPoints[i][2] = coords[2];
        }
        double ratio = tet_radius_ratio(vallPoints);
        if (ratio < metric_list.min_radius_ratio){
            metric_list.min_radius_ratio = ratio;
        }
        if(vertex_active_map.contains(vertexHash(vs))){
            if (vertex_active_map[vertexHash(vs)]){
                metric_list.active_tet++;
                activeTetId.push_back(tid);
                if (ratio < metric_list.active_radius_ratio){
                    metric_list.active_radius_ratio = ratio;
                }
            }
        }
    });
    metric_list.total_tet = grid.get_num_tets();
    metric_list.two_func_check = sub_call_two;
    metric_list.three_func_check = sub_call_three;
    //profiled time(see details in time.h) and profiled number of calls to zero
    std::cout << time_label[0] << ": " << profileTimer[0] << std::endl;

    if (discretize_later){
        /// save the grid output for discretization tool
        save_mesh_json("grid.json", grid);
        /// save the mesh output for isosurfacing tool
        save_function_json("function_value.json", grid, vertex_func_grad_map, funcNum);
    }
    // write mesh and active tets
    mtet::save_mesh("tet_grid.msh", grid);
    mtet::save_mesh("active_tets.msh", grid, std::span<mtet::TetId>(activeTetId));
    return true;
}
