//whitewizard created on Sep 5th
#ifndef WHITEWIZARD_BACKEND_EDGE_H
#define WHITEWIZARD_BACKEND_EDGE_H

#include <memory>
#include <string>
#include <Eigen/Dense>
#include "loss_function.h"

namespace whitewizard_vio{
    namespace backend{
        /*
        边计算残差，预测-观测值，在构造函数中定义
        代价函数：残差*数据*残差，数值在backend中求和后最小化
        */
        class Vertex;
        class edge{
            protected:
            unsigned long id_;  // edge id
            int ordering_id_;   //edge id in problem

            std::vector<std::string> verticies_types_;  // 各顶点类型信息，用于debug
            std::vector<std::shared_ptr<Vertex>> verticies_; // 该边对应的顶点
            Eigen::Matrix<double, Eigen::Dynamic, 1> residual_;                 // 残差
            std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> jacobians_;  // 雅可比，每个雅可比维度是 residual x vertex[i]
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> info_matrix_;             // 信息矩阵
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> sqrt_info_matrix_;
            Eigen::Matrix<double, Eigen::Dynamic, 1> observation_;              // 观测信息

            LossFunction *lossfunction_;

            public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
            /*
            构造函数，配置雅可比空间
            residual_dim 残差维度
            num_verticies 定点数量
            verticies_types 定点类型
            */
            explicit edge(int residual_dim, int num_verticies, const std::vector<std::string> &verticies_types);

            virtual ~edge();

            //get id
            unsigned long Id() const {return id_;}

            /*
            设置一个顶点
            vertex 顶点对应的ertex对象
             */
            bool AddVertex(std::shared_ptr<Vertex> vertex) {
                verticies_.emplace_back(vertex);
                return true;
            }
            /*
            创建一个顶点
            vertices 按引用顺序排列
            */
            bool SetVertex(const std::vector<std::shared_ptr<Vertex>> &vertices){
                verticies_ = vertices;
                return true;
            }

            /*
            返回一个顶点
            */
            std::shared_ptr<Vertex> GetVertex(int i) {
                return verticies_[i];
            }
            
            /*
            返回所有顶点
            */
            std::vector<std::shared_ptr<Vertex>> Verticies() const {
                return verticies_;
            }

            /*
            返回关联顶点个数
            */
            size_t NumVertices() const { return verticies_.size(); }

            /*
            计算雅可比
            本后端不支持自动求导，需要实现每个子类的雅可比计算方法
            后面的残差计算也都实现在子类中
            */

            virtual void ComputeJacobians() = 0;
            /*
             返回边的类型信息
            */
            virtual std::string TypeInfo() const = 0;

            /*
            计算残差
            */
            virtual void ComputeResidual() = 0;

            // 计算平方误差，会乘以信息矩阵
            double Chi2() const;
            double RobustChi2() const;

            // 返回残差
            Eigen::Matrix<double, Eigen::Dynamic, 1> Residual() const { return residual_; }

            // 返回雅可比
            std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> Jacobians() const { return jacobians_; }

            // 设置信息矩阵
            void SetInformation(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &information) {
                info_matrix_ = information;
                sqrt_info_matrix_ = Eigen::LLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(info_matrix_).matrixL().transpose();
            }

            // 返回信息矩阵
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Information() const {
                return info_matrix_;
            }

            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> SqrtInformation() const {
                return sqrt_info_matrix_;
            }

            void SetLossFunction(LossFunction* ptr){ lossfunction_ = ptr; }
            LossFunction* GetLossFunction(){ return lossfunction_;}
            void RobustInfo(double& drho, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& info) const;

            // 设置观测信息
            void SetObservation(const Eigen::Matrix<double, Eigen::Dynamic, 1> &observation) {
                observation_ = observation;
            }

            // 返回观测信息
            Eigen::Matrix<double, Eigen::Dynamic, 1> Observation() const { return observation_; }

            // 检查边的信息是否全部设置
            bool CheckValid();

            int OrderingId() const { return ordering_id_; }

            void SetOrderingId(int id) { ordering_id_ = id; };
        };
    }
}
#endif