//顶点属性配置
#ifndef WHITEWIZARD_BACKEND_VERTEX_H
#define WHITEWIZARD_BACKEND_VERTEX_H


#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <map>

namespace whitewizard_vio{
namespace backend{
    extern unsigned long global_vertex_id;

    class Vertex
    {
    protected:
        Eigen::Matrix<double, Eigen::Dynamic, 1> parameters_;   // 实际存储的变量值
        Eigen::Matrix<double, Eigen::Dynamic, 1> parameters_backup_; // 每次迭代优化中对参数进行备份，用于回滚
        int local_dimension_;   // 局部参数化维度
        unsigned long id_;  // 顶点的id，自动生成
        unsigned long ordering_id_ = 0;

        bool fixed_ = false;    // 是否固定
    public:
        explicit Vertex(int num_dimension, int local_dimension = -1);
        virtual ~Vertex();

        // 返回变量维度
        int Dimension() const;

        /// 返回变量本地维度
        int LocalDimension() const;

        // 该顶点的id
        unsigned long Id() const { return id_; }

        // 返回参数值
        Eigen::Matrix<double, Eigen::Dynamic, 1> Parameters() const { return parameters_; }

        // 返回参数值的引用
        Eigen::Matrix<double, Eigen::Dynamic, 1> &Parameters() { return parameters_; }

        // 设置参数值
        void SetParameters(const Eigen::Matrix<double, Eigen::Dynamic, 1> &params) { parameters_ = params; }

        // 备份和回滚参数，用于丢弃一些迭代过程中不好的估计
        void BackUpParameters() { parameters_backup_ = parameters_; }
        void RollBackParameters() { parameters_ = parameters_backup_; }

        // 加法，可重定义
        // 默认是向量加
        virtual void Plus(const Eigen::Matrix<double, Eigen::Dynamic, 1> &delta);

        // 返回顶点的名称，在子类中实现
        virtual std::string TypeInfo() const = 0;

        int OrderingId() const { return ordering_id_; }

        void SetOrderingId(unsigned long id) { ordering_id_ = id; };

        // 固定该点的估计值
        void SetFixed(bool fixed = true) {
            fixed_ = fixed;
        }

        // 测试该点是否被固定
        bool IsFixed() const { return fixed_; }
    };
     
}
}
#endif