//IMU误差，边为四元边,相连的顶点有 Pi Mi Pj Mj
#ifndef WHITEWIZARD_BACKEND_IMU_EDGE_H
#define WHITEWIZARD_BACKEND_IMU_EDGE_H

#include <memory>
#include <string>
#include "../thirdparty/Sophus/sophus/se3.hpp"

#include "backend/edge.h"
#include "../integration_base.h"

namespace whitewizard_vio{
namespace backend{
    class edgeimu : public edge{
        private:
            enum StateOrder
        {
            O_P = 0,
            O_R = 3,
            O_V = 6,
            O_BA = 9,
            O_BG = 12
        };
            IntegrationBase* pre_integration_;
            static Eigen::Matrix<double, 3, 1> gravity_;

            Eigen::Matrix<double, 3, 3> dp_dba_ = Eigen::Matrix<double, 3, 3>::Zero();
            Eigen::Matrix<double, 3, 3> dp_dbg_ = Eigen::Matrix<double, 3, 3>::Zero();
            Eigen::Matrix<double, 3, 3> dr_dbg_ = Eigen::Matrix<double, 3, 3>::Zero();
            Eigen::Matrix<double, 3, 3> dv_dba_ = Eigen::Matrix<double, 3, 3>::Zero();
            Eigen::Matrix<double, 3, 3> dv_dbg_ = Eigen::Matrix<double, 3, 3>::Zero();
        

        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

            explicit edgeimu(IntegrationBase* _pre_integration):pre_integration_(_pre_integration),
                            edge(15, 4, std::vector<std::string>{"VertexPose", "VertexSpeedBias", "VertexPose", "VertexSpeedBias"}){}
        
            // 返回边的类型信息
            virtual std::string TypeInfo() const override { return "edgeimu"; }

            // 计算残差
            virtual void ComputeResidual() override;

            // 计算雅可比
            virtual void ComputeJacobians() override;
        };
    }
}
#endif
