#include "iostream"
#include "fstream"
#include "Eigen/Dense"
#include "iomanip"
#include "backend/solve.h"
#include "time_check.h"

using namespace std;
//以一定的格式写入CSV表格
void WriteToCSV(std::string name, Eigen::MatrixXd matrix){
    std::ofstream f(name.c_str());
    const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");
    f << matrix.format(CSVFormat);
}

namespace whitewizard_vio{
namespace backend{
    solve::solve(ProblemType problemtype) :
        problemType_(problemtype){
            verticies_marg_.clear();
        }

    solve::~solve(){
        global_vertex_id = 0;
        // std::cout << "solver is deleted" << std::endl;
    }
    // 加入顶点
    bool solve::AddVertex(std::shared_ptr<Vertex> vertex) {
        if(verticies_.find(vertex->Id()) != verticies_.end()){
        // std::cout << "this vertex" << vertex->Id() << " has been added" << std::endl;
        return false;
    }
        else{
        verticies_.insert(pair<unsigned long, shared_ptr<Vertex>>(vertex->Id(), vertex));
        //不是重复的顶点插入
    }

    if(problemType_ == ProblemType::SLAM_PROBLEM)//如果是slam问题要以稀疏方式存储，hessian矩阵要resize
        if (IsPoseVertex(vertex)) {
            ResizePoseHessiansWhenAddingPose(vertex);
            }
    return true;
    }
    
    void solve::AddOrderingSLAM(std::shared_ptr<whitewizard_vio::backend::Vertex> v){
        if(IsPoseVertex(v)) {
            v->SetOrderingId(ordering_poses_);
            idx_pose_vertices_.insert(pair<ulong, std::shared_ptr<Vertex>>(v->Id(), v));
            ordering_poses_ += v -> LocalDimension();  
        }else if(IsLandmarkVertex(v)){
            v->SetOrderingId(ordering_landmarks_);
            ordering_landmarks_ += v->LocalDimension();
            idx_landmark_vertices_.insert(pair<ulong, std::shared_ptr<Vertex>>(v->Id(), v));
        }
    }

    void solve::ResizePoseHessiansWhenAddingPose(shared_ptr<Vertex> v){
        int size = H_prior_.rows() + v->LocalDimension();
        H_prior_.conservativeResize(size, size);
        b_prior_.conservativeResize(size);

        b_prior_.tail(v->LocalDimension()).setZero();
        H_prior_.rightCols(v->LocalDimension()).setZero();
        H_prior_.bottomRows(v->LocalDimension()).setZero();
    } 

    void solve::ExtendHessiansPriorSize(int dim)
    {
        int size = H_prior_.rows() + dim;
        H_prior_.conservativeResize(size, size);
        b_prior_.conservativeResize(size);

        b_prior_.tail(dim).setZero();
        H_prior_.rightCols(dim).setZero();
        H_prior_.bottomRows(dim).setZero();
    }

    bool solve::IsPoseVertex(std::shared_ptr<whitewizard_vio::backend::Vertex> v){
        std::string type = v->TypeInfo();
        return type == std::string("vertex pose") || type == std::string("vertex inverse depth");
    }

    bool solve::IsLandmarkVertex(std::shared_ptr<whitewizard_vio::backend::Vertex> v){
        std::string type = v->TypeInfo();
        return type == string("vertex point XYZ") || type == string("vertex inverse depth");
    }

    bool solve::AddEdge(std::shared_ptr<edge> edges){
        if(edges_.find(edges->Id()) == edges_.end()){
            edges_.insert(pair<ulong, std::shared_ptr<edge>>(edges->Id(), edges));
        } else {
            // std::cout << "this edge " << edges->Id() << "has been added" << std::endl;
            return false;
        }

        for(auto &vertex : edges->Verticies() ){
            vertexToEdge_.insert(pair<ulong, shared_ptr<edge>>(vertex->Id(), edges));
        }
        return true;
    }

    std::vector<shared_ptr<edge>> solve::GetConnectedEdges(std::shared_ptr<Vertex> vertex){
        std::vector<shared_ptr<edge>> edges;
        auto range = vertexToEdge_.equal_range(vertex->Id());

        for (auto i = range.first; i != range.second; i++)
        {
            if(edges_.find(i->second->Id()) == edges_.end())
                continue;//这个边还需要存在

            edges.emplace_back(i->second);
        }
        return edges;
    }

    bool solve::RemoveVertex(std::shared_ptr<Vertex> vertex){
        //先找到这个顶点
        if(verticies_.find(vertex->Id()) == verticies_.end()){
            std::cout << "this vertex " << vertex->Id() << " is not in the queue" << std::endl;
            return false;
        }
        std::vector<shared_ptr<edge>> remove_edge = GetConnectedEdges(vertex);
        for (size_t i = 0; i < remove_edge.size(); i++)
        {
            RemoveEdge(remove_edge[i]);
        }

        if(IsPoseVertex(vertex))
            idx_pose_vertices_.erase(vertex->Id());
        else
            idx_landmark_vertices_.erase(vertex->Id());

        vertex->SetOrderingId(-1);
        verticies_.erase(vertex->Id());
        vertexToEdge_.erase(vertex->Id());
        
        return true;
    }

    bool solve::RemoveEdge(std::shared_ptr<edge> edges){
        //先找到这个边
        if(edges_.find(edges->Id()) == edges_.end() ){
            std::cout << "this edge " << edges->Id() << " is not in the queue" << std::endl;
            return false;
        }

        edges_.erase(edges->Id());
        return true;
    }

    bool solve::Solve(int iterations){
        if(edges_.size() == 0 || verticies_.size() == 0){
            return false;
        }

        timecheck t_solve;
        SetOrdering();//统计维数，准备构建H矩阵
        MakeHessian();//遍历edge，构建H矩阵
        ComputeLambdaInitLM();//LM初始化

        //LM迭代——————————————————————————————————————
        bool stop = false;
        int iter = 0;
        double last_chi_ = 1e20;
        while (!stop && (iter < iterations))
        {
            std::cout << "iter: " << iter << " , chi= " << currentChi_ << " , Lambda= " << currentLambda_ << std::endl;
            bool calSuccess = false;
            int false_cnt = 0;
            while (!calSuccess && false_cnt < 10 )
            {
                //  解线性方程
                SolveLinearSystem();
                //  状态更新
                UpdateStates();
                //  判断当前是否可行，判断LM 的阻尼因子lambda更新策略
                calSuccess = IsGoodStepInLM();
                //  如果计算出来了，构建hessian矩阵
                if(calSuccess){
                    MakeHessian();
                    false_cnt = 0;
                } else{
                    false_cnt++;
                    RollbackStates();//回滚，误差没有下降
                }
            }
            iter++;
            //当误差和第一次相比数量级小于1e-5则退出优化
            if(last_chi_ - currentChi_ < 1e-5)
            {
                stop = true;
            }
            last_chi_ = currentChi_;
        }
        t_hessian_cost_ = 0;
        return true;
    }

    bool solve::SolveGenericProblem(int iterations) {
    return true;
    }

    void solve::SetOrdering(){
        //重新记次
        ordering_poses_ = 0;
        ordering_generic_ = 0;
        ordering_landmarks_ = 0;
        for (auto vertex: verticies_) {
            ordering_generic_ += vertex.second->LocalDimension();  // 所有的优化变量总维数

            if (problemType_ == ProblemType::SLAM_PROBLEM)    // 如果是 slam 问题，还要分别统计 pose 和 landmark 的维数，后面会对他们进行排序
            {
                AddOrderingSLAM(vertex.second);
            }

        }

        if (problemType_ == ProblemType::SLAM_PROBLEM) {
            // 这里要把 landmark 的 ordering 加上 pose 的数量，就保持了 landmark 在后,而 pose 在前
            ulong all_pose_dimension = ordering_poses_;
            for (auto landmarkVertex : idx_landmark_vertices_) {
                landmarkVertex.second->SetOrderingId(
                    landmarkVertex.second->OrderingId() + all_pose_dimension
                );
            }
        }

    }

    bool solve::CheckOrdering() {
        if (problemType_ == ProblemType::SLAM_PROBLEM) {
            int current_ordering = 0;
            for (auto v: idx_pose_vertices_) {
                assert(v.second->OrderingId() == current_ordering);
                current_ordering += v.second->LocalDimension();
            }

            for (auto v: idx_landmark_vertices_) {
                assert(v.second->OrderingId() == current_ordering);
                current_ordering += v.second->LocalDimension();
            }
        }
        return true;
    }

    void solve::MakeHessian() {
        timecheck t_h;
        // 直接构造大的 H 矩阵
        ulong size = ordering_generic_;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(size, size));
        Eigen::Matrix<double, Eigen::Dynamic, 1> b(Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(size));

        for (auto &edge: edges_) {

            edge.second->ComputeResidual();
            edge.second->ComputeJacobians();

            auto jacobians = edge.second->Jacobians();
            auto verticies = edge.second->Verticies();
            assert(jacobians.size() == verticies.size());
            for (size_t i = 0; i < verticies.size(); ++i) {
                auto v_i = verticies[i];
                if (v_i->IsFixed()) continue;    // Hessian 里不需要添加它的信息，也就是它的雅克比为 0

                auto jacobian_i = jacobians[i];
                ulong index_i = v_i->OrderingId();
                ulong dim_i = v_i->LocalDimension();

                // 鲁棒核函数会修改残差和信息矩阵，如果没有设置 robust cost function，就会返回原来的
                double drho;
                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> robustInfo(edge.second->Information().rows(),edge.second->Information().cols());
                edge.second->RobustInfo(drho,robustInfo);

                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> JtW = jacobian_i.transpose() * robustInfo;
                for (size_t j = i; j < verticies.size(); ++j) {
                    auto v_j = verticies[j];

                    if (v_j->IsFixed()) continue;

                    auto jacobian_j = jacobians[j];
                    ulong index_j = v_j->OrderingId();
                    ulong dim_j = v_j->LocalDimension();

                    assert(v_j->OrderingId() != -1);
                    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> hessian = JtW * jacobian_j;

                    // 所有的信息矩阵叠加起来
                    H.block(index_i, index_j, dim_i, dim_j).noalias() += hessian;
                    if (j != i) {
                        // 对称的下三角
                        H.block(index_j, index_i, dim_j, dim_i).noalias() += hessian.transpose();

                        }
                    }
                    b.segment(index_i, dim_i).noalias() -= drho * jacobian_i.transpose()* edge.second->Information() * edge.second->Residual();
                }

            }
        Hessian_ = H;
        b_ = b;
        t_hessian_cost_ += t_h.sum_time();

        if(H_prior_.rows() > 0)
        {
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H_prior_tmp = H_prior_;
            Eigen::Matrix<double, Eigen::Dynamic, 1> b_prior_tmp = b_prior_;

            /// 遍历所有 POSE 顶点，然后设置相应的先验维度为 0 .  fix 外参数, SET PRIOR TO ZERO
            /// landmark 没有先验
            for (auto vertex: verticies_) {
                if (IsPoseVertex(vertex.second) && vertex.second->IsFixed() ) {
                    int idx = vertex.second->OrderingId();
                    int dim = vertex.second->LocalDimension();
                    H_prior_tmp.block(idx,0, dim, H_prior_tmp.cols()).setZero();
                    H_prior_tmp.block(0,idx, H_prior_tmp.rows(), dim).setZero();
                    b_prior_tmp.segment(idx,dim).setZero();
    //                std::cout << " fixed prior, set the Hprior and bprior part to zero, idx: "<<idx <<" dim: "<<dim<<std::endl;
                }
            }
            Hessian_.topLeftCorner(ordering_poses_, ordering_poses_) += H_prior_tmp;
            b_.head(ordering_poses_) += b_prior_tmp;
        }

        delta_x_ = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(size);  // initial delta_x = 0_n;


    }

    /*
    解算 Hx = b,可以用PCG
    */
    void solve::SolveLinearSystem() {


        if (problemType_ == ProblemType::GENERIC_PROBLEM) {
            // PCG solver
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H = Hessian_;
            for (size_t i = 0; i < Hessian_.cols(); ++i) {
                H(i, i) += currentLambda_;
            }
            // delta_x_ = PCGSolver(H, b_, H.rows() * 2);
            delta_x_ = H.ldlt().solve(b_);

        } else {

            int reserve_size = ordering_poses_;
            int marg_size = ordering_landmarks_;
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Hmm = Hessian_.block(reserve_size, reserve_size, marg_size, marg_size);
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Hpm = Hessian_.block(0, reserve_size, reserve_size, marg_size);
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Hmp = Hessian_.block(reserve_size, 0, marg_size, reserve_size);
            Eigen::Matrix<double, Eigen::Dynamic, 1> bpp = b_.segment(0, reserve_size);
            Eigen::Matrix<double, Eigen::Dynamic, 1> bmm = b_.segment(reserve_size, marg_size);

            // Hmm 是对角线矩阵，它的求逆可以直接为对角线块分别求逆，如果是逆深度，对角线块为1维的，则直接为对角线的倒数，这里可以加速
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Hmm_inv(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(marg_size, marg_size));

            for (auto landmarkVertex : idx_landmark_vertices_) {
                int idx = landmarkVertex.second->OrderingId() - reserve_size;
                int size = landmarkVertex.second->LocalDimension();
                Hmm_inv.block(idx, idx, size, size) = Hmm.block(idx, idx, size, size).inverse();
            }

            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> tempH = Hpm * Hmm_inv;
            H_pp_schur_ = Hessian_.block(0, 0, ordering_poses_, ordering_poses_) - tempH * Hmp;
            b_pp_schur_ = bpp - tempH * bmm;

            Eigen::Matrix<double, Eigen::Dynamic, 1> delta_x_pp(Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(reserve_size));

            for (ulong i = 0; i < ordering_poses_; ++i) {
                H_pp_schur_(i, i) += currentLambda_;              // LM Method
            }

            // timecheck t_linearsolver;
            delta_x_pp =  H_pp_schur_.ldlt().solve(b_pp_schur_);//  SVec.asDiagonal() * svd.matrixV() * Ub;    
            delta_x_.head(reserve_size) = delta_x_pp;
            // std::cout << " Linear Solver Time Cost: " << t_linearsolver.sum_time() << std::endl;

            // step3: solve Hmm * delta_x = bmm - Hmp * delta_x_pp;
            Eigen::Matrix<double, Eigen::Dynamic, 1> delta_x_ll(marg_size);
            delta_x_ll = Hmm_inv * (bmm - Hmp * delta_x_pp);
            delta_x_.tail(marg_size) = delta_x_ll;

    //        std::cout << "schur time cost: "<< t_Hmminv.sum_time()<<std::endl;
        }

    }

    void solve::UpdateStates() {

        // update vertex
        for (auto vertex: verticies_) {
            vertex.second->BackUpParameters();    // 保存上次的估计值

            ulong idx = vertex.second->OrderingId();
            ulong dim = vertex.second->LocalDimension();
            Eigen::Matrix<double, Eigen::Dynamic, 1> delta = delta_x_.segment(idx, dim);
            vertex.second->Plus(delta);
        }

        // update prior
        if (err_prior_.rows() > 0) {
            // BACK UP b_prior_
            b_prior_backup_ = b_prior_;
            err_prior_backup_ = err_prior_;

            /// update with first order Taylor, b' = b + \frac{\delta b}{\delta x} * \delta x
            /// \delta x = Computes the linearized deviation from the references (linearization points)
            b_prior_ -= H_prior_ * delta_x_.head(ordering_poses_);       // update the error_prior
            err_prior_ = -Jt_prior_inv_ * b_prior_.head(ordering_poses_ - 15);

        }

    }

    void solve::RollbackStates() {

        // update vertex
        for (auto vertex: verticies_) {
            vertex.second->RollBackParameters();
        }

        // Roll back prior_
        if (err_prior_.rows() > 0) {
            b_prior_ = b_prior_backup_;
            err_prior_ = err_prior_backup_;
        }
    }

    /// LM
    void solve::ComputeLambdaInitLM() {
        ni_ = 2.;
        currentLambda_ = -1.;
        currentChi_ = 0.0;

        for (auto edge: edges_) {
            currentChi_ += edge.second->RobustChi2();
        }
        if (err_prior_.rows() > 0)
            currentChi_ += err_prior_.norm();
        currentChi_ *= 0.5;

        stopThresholdLM_ = 1e-10 * currentChi_;          // 迭代条件为 误差下降 1e-10 倍

        double maxDiagonal = 0;
        ulong size = Hessian_.cols();
        assert(Hessian_.rows() == Hessian_.cols() && "Hessian is not square");
        for (ulong i = 0; i < size; ++i) {
            maxDiagonal = std::max(fabs(Hessian_(i, i)), maxDiagonal);
        }

        maxDiagonal = std::min(5e10, maxDiagonal);
        double tau = 1e-5;  // 1e-5
        currentLambda_ = tau * maxDiagonal;
    //        std::cout << "currentLamba_: "<<maxDiagonal<<" "<<currentLambda_<<std::endl;
    }

    void solve::AddLambdatoHessianLM() {
        ulong size = Hessian_.cols();
        assert(Hessian_.rows() == Hessian_.cols() && "Hessian is not square");
        for (ulong i = 0; i < size; ++i) {
            Hessian_(i, i) += currentLambda_;
        }
    }

    void solve::RemoveLambdaHessianLM() {
        ulong size = Hessian_.cols();
        assert(Hessian_.rows() == Hessian_.cols() && "Hessian is not square");
        for (ulong i = 0; i < size; ++i) {
            Hessian_(i, i) -= currentLambda_;
        }
    }

    bool solve::IsGoodStepInLM() {
        double scale = 0;
    //    scale = 0.5 * delta_x_.transpose() * (currentLambda_ * delta_x_ + b_);
    //    scale += 1e-3;    // 确保这个值不为零
        scale = 0.5* delta_x_.transpose() * (currentLambda_ * delta_x_ + b_);
        scale += 1e-6;    // 确保这个值不为零

        // recompute residuals after update state
        double tempChi = 0.0;
        for (auto edge: edges_) {
            edge.second->ComputeResidual();
            tempChi += edge.second->RobustChi2();
        }
        if (err_prior_.size() > 0)
            tempChi += err_prior_.norm();
        tempChi *= 0.5;          // 1/2 * err^2

        double rho = (currentChi_ - tempChi) / scale;
        if (rho > 0 && isfinite(tempChi))   // last step was good, 误差在下降
        {
            double alpha = 1. - pow((2 * rho - 1), 3);
            alpha = std::min(alpha, 2. / 3.);
            double scaleFactor = (std::max)(1. / 3., alpha);
            currentLambda_ *= scaleFactor;
            ni_ = 2;
            currentChi_ = tempChi;
            return true;
        } else {
            currentLambda_ *= ni_;
            ni_ *= 2;
            return false;
        }
    }


    /// PCG雅可比
    Eigen::Matrix<double, Eigen::Dynamic, 1> solve::PCGSolver(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &A, const Eigen::Matrix<double, Eigen::Dynamic, 1> &b, int maxIter = -1) {
        assert(A.rows() == A.cols() && "PCG solver ERROR: A is not square matrix");
        int rows = b.rows();
        int n = maxIter < 0 ? rows : maxIter;
        Eigen::Matrix<double, Eigen::Dynamic, 1> x(Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(rows));
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> M_inv = A.diagonal().asDiagonal().inverse();
        Eigen::Matrix<double, Eigen::Dynamic, 1> r0(b);  // initial r = b - A*0 = b
        Eigen::Matrix<double, Eigen::Dynamic, 1> z0 = M_inv * r0;
        Eigen::Matrix<double, Eigen::Dynamic, 1> p(z0);
        Eigen::Matrix<double, Eigen::Dynamic, 1> w = A * p;
        double r0z0 = r0.dot(z0);
        double alpha = r0z0 / p.dot(w);
        Eigen::Matrix<double, Eigen::Dynamic, 1> r1 = r0 - alpha * w;
        int i = 0;
        double threshold = 1e-6 * r0.norm();
        while (r1.norm() > threshold && i < n) {
            i++;
            Eigen::Matrix<double, Eigen::Dynamic, 1> z1 = M_inv * r1;
            double r1z1 = r1.dot(z1);
            double belta = r1z1 / r0z0;
            z0 = z1;
            r0z0 = r1z1;
            r0 = r1;
            p = belta * p + z1;
            w = A * p;
            alpha = r1z1 / p.dot(w);
            x += alpha * p;
            r1 -= alpha * w;
        }
        return x;
    }


    //marg 所有和 frame 相连的 edge: imu factor, projection factor
    //如果某个landmark和该frame相连，但是又不想加入marg, 那就把改edge先去掉

    bool solve::Marginalize(const std::vector<std::shared_ptr<Vertex> > margVertexs, int pose_dim) {

        SetOrdering();
        /// 找到需要 marg 的 edge, margVertexs[0] is frame, its edge contained pre-intergration
        std::vector<shared_ptr<edge>> marg_edges = GetConnectedEdges(margVertexs[0]);

        std::unordered_map<int, shared_ptr<Vertex>> margLandmark;
        // 构建 Hessian 的时候 pose 的顺序不变，landmark的顺序要重新设定
        int marg_landmark_size = 0;
    //    std::cout << "\n marg edge 1st id: "<< marg_edges.front()->Id() << " end id: "<<marg_edges.back()->Id()<<std::endl;
        for (size_t i = 0; i < marg_edges.size(); ++i) {
    //        std::cout << "marg edge id: "<< marg_edges[i]->Id() <<std::endl;
            auto verticies = marg_edges[i]->Verticies();
            for (auto iter : verticies) {
                if (IsLandmarkVertex(iter) && margLandmark.find(iter->Id()) == margLandmark.end()) {
                    iter->SetOrderingId(pose_dim + marg_landmark_size);
                    margLandmark.insert(make_pair(iter->Id(), iter));
                    marg_landmark_size += iter->LocalDimension();
                }
            }
        }
    //    std::cout << "pose dim: " << pose_dim <<std::endl;
        int cols = pose_dim + marg_landmark_size;
        /// 构建误差 H 矩阵 H = H_marg + H_pp_prior
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H_marg(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(cols, cols));
        Eigen::Matrix<double, Eigen::Dynamic, 1> b_marg(Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(cols));
        int ii = 0;
        for (auto edge: marg_edges) {
            edge->ComputeResidual();
            edge->ComputeJacobians();
            auto jacobians = edge->Jacobians();
            auto verticies = edge->Verticies();
            ii++;

            assert(jacobians.size() == verticies.size());
            for (size_t i = 0; i < verticies.size(); ++i) {
                auto v_i = verticies[i];
                auto jacobian_i = jacobians[i];
                ulong index_i = v_i->OrderingId();
                ulong dim_i = v_i->LocalDimension();

                double drho;
                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> robustInfo(edge->Information().rows(),edge->Information().cols());
                edge->RobustInfo(drho,robustInfo);

                for (size_t j = i; j < verticies.size(); ++j) {
                    auto v_j = verticies[j];
                    auto jacobian_j = jacobians[j];
                    ulong index_j = v_j->OrderingId();
                    ulong dim_j = v_j->LocalDimension();

                    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> hessian = jacobian_i.transpose() * robustInfo * jacobian_j;

                    assert(hessian.rows() == v_i->LocalDimension() && hessian.cols() == v_j->LocalDimension());
                    // 所有的信息矩阵叠加起来
                    H_marg.block(index_i, index_j, dim_i, dim_j) += hessian;
                    if (j != i) {
                        // 对称的下三角
                        H_marg.block(index_j, index_i, dim_j, dim_i) += hessian.transpose();
                    }
                }
                b_marg.segment(index_i, dim_i) -= drho * jacobian_i.transpose() * edge->Information() * edge->Residual();
            }

        }
            std::cout << "edge factor cnt: " << ii <<std::endl;

        /// marg landmark
        int reserve_size = pose_dim;
        if (marg_landmark_size > 0) {
            int marg_size = marg_landmark_size;
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Hmm = H_marg.block(reserve_size, reserve_size, marg_size, marg_size);
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Hpm = H_marg.block(0, reserve_size, reserve_size, marg_size);
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Hmp = H_marg.block(reserve_size, 0, marg_size, reserve_size);
            Eigen::Matrix<double, Eigen::Dynamic, 1> bpp = b_marg.segment(0, reserve_size);
            Eigen::Matrix<double, Eigen::Dynamic, 1> bmm = b_marg.segment(reserve_size, marg_size);

            // Hmm 是对角线矩阵，它的求逆可以直接为对角线块分别求逆，如果是逆深度，对角线块为1维的，则直接为对角线的倒数，这里可以加速
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Hmm_inv(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(marg_size, marg_size));
            for (auto iter: margLandmark) {
                int idx = iter.second->OrderingId() - reserve_size;
                int size = iter.second->LocalDimension();
                Hmm_inv.block(idx, idx, size, size) = Hmm.block(idx, idx, size, size).inverse();
            }

            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> tempH = Hpm * Hmm_inv;
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Hpp = H_marg.block(0, 0, reserve_size, reserve_size) - tempH * Hmp;
            bpp = bpp - tempH * bmm;
            H_marg = Hpp;
            b_marg = bpp;
        }

        Eigen::Matrix<double, Eigen::Dynamic, 1> b_prior_before = b_prior_;
        if(H_prior_.rows() > 0)
        {
            H_marg += H_prior_;
            b_marg += b_prior_;
        }

        /// marg frame and speedbias
        int marg_dim = 0;

        // index 大的先移动
        for (int k = margVertexs.size() -1 ; k >= 0; --k)
        {

            int idx = margVertexs[k]->OrderingId();
            int dim = margVertexs[k]->LocalDimension();
    //        std::cout << k << " "<<idx << std::endl;
            marg_dim += dim;
            // 将 row i 移动矩阵最下面
            Eigen::MatrixXd temp_rows = H_marg.block(idx, 0, dim, reserve_size);
            Eigen::MatrixXd temp_botRows = H_marg.block(idx + dim, 0, reserve_size - idx - dim, reserve_size);
            H_marg.block(idx, 0, reserve_size - idx - dim, reserve_size) = temp_botRows;
            H_marg.block(reserve_size - dim, 0, dim, reserve_size) = temp_rows;

            // 将 col i 移动矩阵最右边
            Eigen::MatrixXd temp_cols = H_marg.block(0, idx, reserve_size, dim);
            Eigen::MatrixXd temp_rightCols = H_marg.block(0, idx + dim, reserve_size, reserve_size - idx - dim);
            H_marg.block(0, idx, reserve_size, reserve_size - idx - dim) = temp_rightCols;
            H_marg.block(0, reserve_size - dim, reserve_size, dim) = temp_cols;

            Eigen::VectorXd temp_b = b_marg.segment(idx, dim);
            Eigen::VectorXd temp_btail = b_marg.segment(idx + dim, reserve_size - idx - dim);
            b_marg.segment(idx, reserve_size - idx - dim) = temp_btail;
            b_marg.segment(reserve_size - dim, dim) = temp_b;
        }

        double eps = 1e-8;
        int m2 = marg_dim;
        int n2 = reserve_size - marg_dim;   // marg pose
        Eigen::MatrixXd Amm = 0.5 * (H_marg.block(n2, n2, m2, m2) + H_marg.block(n2, n2, m2, m2).transpose());

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes(Amm);
        Eigen::MatrixXd Amm_inv = saes.eigenvectors() * Eigen::VectorXd(
                (saes.eigenvalues().array() > eps).select(saes.eigenvalues().array().inverse(), 0)).asDiagonal() *
                                saes.eigenvectors().transpose();

        Eigen::VectorXd bmm2 = b_marg.segment(n2, m2);
        Eigen::MatrixXd Arm = H_marg.block(0, n2, n2, m2);
        Eigen::MatrixXd Amr = H_marg.block(n2, 0, m2, n2);
        Eigen::MatrixXd Arr = H_marg.block(0, 0, n2, n2);
        Eigen::VectorXd brr = b_marg.segment(0, n2);
        Eigen::MatrixXd tempB = Arm * Amm_inv;
        H_prior_ = Arr - tempB * Amr;
        b_prior_ = brr - tempB * bmm2;

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes2(H_prior_);
        Eigen::VectorXd S = Eigen::VectorXd((saes2.eigenvalues().array() > eps).select(saes2.eigenvalues().array(), 0));
        Eigen::VectorXd S_inv = Eigen::VectorXd(
                (saes2.eigenvalues().array() > eps).select(saes2.eigenvalues().array().inverse(), 0));

        Eigen::VectorXd S_sqrt = S.cwiseSqrt();
        Eigen::VectorXd S_inv_sqrt = S_inv.cwiseSqrt();
        Jt_prior_inv_ = S_inv_sqrt.asDiagonal() * saes2.eigenvectors().transpose();
        err_prior_ = -Jt_prior_inv_ * b_prior_;

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> J = S_sqrt.asDiagonal() * saes2.eigenvectors().transpose();
        H_prior_ = J.transpose() * J;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> tmp_h = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>( (H_prior_.array().abs() > 1e-9).select(H_prior_.array(),0) );
        H_prior_ = tmp_h;

        // remove vertex and remove edge
        for (size_t k = 0; k < margVertexs.size(); ++k) {
            RemoveVertex(margVertexs[k]);
        }

        for (auto landmarkVertex: margLandmark) {
            RemoveVertex(landmarkVertex.second);
        }

        return true;

    }

}
}






