#ifndef WHITEWIZARD_BACKEND_PROBLEM_H
#define WHITEWIZARD_BACKEND_PROBLEM_H

#include <unordered_map>
#include <map>
#include <memory>
#include "edge.h"
#include "vertex.h"

namespace whitewizard_vio{
namespace backend{
    typedef unsigned long ulong;
    //hash元素定义
    typedef std::map<unsigned long, std::shared_ptr<Vertex>> HashVertex;
    typedef std::unordered_map<unsigned long, std::shared_ptr<edge>> HashEdge;
    typedef std::unordered_multimap<unsigned long, std::shared_ptr<edge>> HashVertexIdToEdge;

    class solve
    {
    private:
        
        // Solve的实现，解通用问题
        bool SolveGenericProblem(int iterations);

        // Solve的实现，解SLAM问题
        bool SolveSLAMProblem(int iterations);

        // 设置各顶点的ordering_index
        void SetOrdering();

        // 对新的slam顶点集合进行排序
        void AddOrderingSLAM(std::shared_ptr<Vertex> v);

        // 构造大H矩阵
        void MakeHessian();

        // schur求解SBA
        void SchurSBA();

        // 解线性方程
        void SolveLinearSystem();

        // 更新状态变量
        void UpdateStates();

        void RollbackStates(); // 有时候 update 后残差会变大，需要退回去，重来

        // 计算并更新Prior部分
        void ComputePrior();

        // 判断一个顶点是否为Pose顶点
        bool IsPoseVertex(std::shared_ptr<Vertex> v);

        // 判断一个顶点是否为landmark顶点
        bool IsLandmarkVertex(std::shared_ptr<Vertex> v);

        // 在新增顶点后，需要调整几个hessian的大小
        void ResizePoseHessiansWhenAddingPose(std::shared_ptr<Vertex> v);

        // 检查ordering是否正确
        bool CheckOrdering();

        // 获取某个顶点连接到的边
        std::vector<std::shared_ptr<edge>> GetConnectedEdges(std::shared_ptr<Vertex> vertex);

        // Levenberg
        // 计算LM算法的初始Lambda
        void ComputeLambdaInitLM();

        // Hessian 对角线加上或者减去  Lambda
        void AddLambdatoHessianLM();

        void RemoveLambdaHessianLM();

        // LM 算法中用于判断 Lambda 在上次迭代中是否可以，以及Lambda怎么缩放
        bool IsGoodStepInLM();

        // PCG 迭代线性求解器
        Eigen::Matrix<double, Eigen::Dynamic, 1> PCGSolver(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &A, const Eigen::Matrix<double, Eigen::Dynamic, 1> &b, int maxIter);

        double currentLambda_;
        double currentChi_;
        double stopThresholdLM_;    // LM 迭代退出阈值条件
        double ni_;                 //控制 Lambda 缩放大小



        // 整个信息矩阵
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Hessian_;
        Eigen::Matrix<double, Eigen::Dynamic, 1> b_;
        Eigen::Matrix<double, Eigen::Dynamic, 1> delta_x_;

        // 先验部分信息
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H_prior_;
        Eigen::Matrix<double, Eigen::Dynamic, 1> b_prior_;
        Eigen::Matrix<double, Eigen::Dynamic, 1> b_prior_backup_;
        Eigen::Matrix<double, Eigen::Dynamic, 1> err_prior_backup_;

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Jt_prior_inv_;
        Eigen::Matrix<double, Eigen::Dynamic, 1> err_prior_;

        // SBA的Pose部分
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H_pp_schur_;
        Eigen::Matrix<double, Eigen::Dynamic, 1> b_pp_schur_;
        // Heesian 的 Landmark 和 pose 部分
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H_pp_;
        Eigen::Matrix<double, Eigen::Dynamic, 1> b_pp_;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H_ll_;
        Eigen::Matrix<double, Eigen::Dynamic, 1> b_ll_;

        // all vertices
        HashVertex verticies_;

        // all edges
        HashEdge edges_;

        // 由vertex id查询edge
        HashVertexIdToEdge vertexToEdge_;

        // Ordering related
        ulong ordering_poses_ = 0;
        ulong ordering_landmarks_ = 0;
        ulong ordering_generic_ = 0;
        std::map<unsigned long, std::shared_ptr<Vertex>> idx_pose_vertices_;        // 以ordering排序的pose顶点
        std::map<unsigned long, std::shared_ptr<Vertex>> idx_landmark_vertices_;    // 以ordering排序的landmark顶点

        // verticies need to marg. <Ordering_id_, Vertex>
        HashVertex verticies_marg_;

        bool bDebug = false;
        double t_hessian_cost_ = 0.0;
        double t_PCGsovle_cost_ = 0.0;

    public:

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        /*
        SLAM问题还是通用的问题
        如果是SLAM问题那么pose和landmark是区分开的，Hessian以稀疏方式存储
        SLAM问题只接受一些特定的Vertex和Edge
        如果是通用问题那么hessian是稠密的，除非用户设定某些vertex为marginalized
        */
        enum class ProblemType {
            SLAM_PROBLEM,
            GENERIC_PROBLEM
        };
        ProblemType problemType_;
        solve(ProblemType problemtype);
        ~solve();
        //vertex增减
        bool AddVertex(std::shared_ptr<Vertex> vertex);

        bool RemoveVertex(std::shared_ptr<Vertex> vertex);
        //edge增减
        bool AddEdge(std::shared_ptr<edge> edge);

        bool RemoveEdge(std::shared_ptr<edge> edge);

        /**
         * 取得在优化中被判断为outlier部分的边，方便前端去除outlier
         * @param outlier_edges
         */
        void GetOutlierEdges(std::vector<std::shared_ptr<edge>> &outlier_edges);

        /**
         * 求解此问题
         * @param iterations
         * @return
         */
        bool Solve(int iterations = 10);

        /// 边缘化一个frame和以它为host的landmark
        bool Marginalize(std::shared_ptr<Vertex> frameVertex,
                const std::vector<std::shared_ptr<Vertex>> &landmarkVerticies);

        bool Marginalize(const std::shared_ptr<Vertex> frameVertex);
        bool Marginalize(const std::vector<std::shared_ptr<Vertex> > frameVertex,int pose_dim);

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> GetHessianPrior(){ return H_prior_;}
        Eigen::Matrix<double, Eigen::Dynamic, 1> GetbPrior(){ return b_prior_;}
        Eigen::Matrix<double, Eigen::Dynamic, 1> GetErrPrior(){ return err_prior_;}
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> GetJtPrior(){ return Jt_prior_inv_;}

        void SetHessianPrior(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& H){H_prior_ = H;}
        void SetbPrior(const Eigen::Matrix<double, Eigen::Dynamic, 1>& b){b_prior_ = b;}
        void SetErrPrior(const Eigen::Matrix<double, Eigen::Dynamic, 1>& b){err_prior_ = b;}
        void SetJtPrior(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& J){Jt_prior_inv_ = J;}

        void ExtendHessiansPriorSize(int dim);

        //test compute prior
        void TestComputePrior();

    };
    
}
}
#endif