//先验edge的边是一元的，顶点有Ti
#ifndef WHITEWIZARD_BACKEND_EDGE_PRIOR_H
#define WHITEWIZARD_BACKEND_EDGE_PRIOR_H

#include <memory>
#include <string>

#include <Eigen/Dense>

#include "edge.h"

namespace whitewizard_vio{
namespace backend{
    class edgepiror : public edge{
        private:
            Eigen::Matrix<double, 3, 1> Pp_;   // pose prior
            Eigen::Quaterniond   Qp_;   // Rotation prior
        
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

            edgepiror(const Eigen::Matrix<double, 3, 1> &p, const Eigen::Quaterniond &q):
                edge(6, 1, std::vector<std::string>{"VertexPose"}),
                Pp_(p), Qp_(q) {}
            
            // 返回边的类型信息
            virtual std::string TypeInfo() const override { return "edgepiror"; }

            // 计算残差
            virtual void ComputeResidual() override;

            // 计算雅可比
            virtual void ComputeJacobians() override;
        };
    }
}
#endif