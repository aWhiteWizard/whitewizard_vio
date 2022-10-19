#include "backend/edge_imu.h"
#include "backend/vertex_pose.h"
#include "backend/vertex_speedbias.h"

#include <iostream>

namespace whitewizard_vio{
namespace backend{
using Sophus::SO3d;

    Eigen::Matrix<double, 3, 1> edgeimu::gravity_ = Eigen::Matrix<double, 3, 1> (0, 0, 9.8);
    void edgeimu::ComputeResidual(){
        Eigen::Matrix<double, Eigen::Dynamic, 1> param_0 = verticies_[0]->Parameters();
        Eigen::Quaterniond Qi(param_0[6], param_0[3], param_0[4], param_0[5]);
        Eigen::Matrix<double, 3, 1> Pi = param_0.head<3>();

        Eigen::Matrix<double, Eigen::Dynamic, 1> param_1 = verticies_[1]->Parameters();
        Eigen::Matrix<double, 3, 1> Vi = param_1.head<3>();
        Eigen::Matrix<double, 3, 1> Bai = param_1.segment(3, 3);
        Eigen::Matrix<double, 3, 1> Bgi = param_1.tail<3>();

        Eigen::Matrix<double, Eigen::Dynamic, 1> param_2 = verticies_[2]->Parameters();
        Eigen::Quaterniond Qj(param_2[6], param_2[3], param_2[4], param_2[5]);
        Eigen::Matrix<double, 3, 1> Pj = param_2.head<3>();

        Eigen::Matrix<double, Eigen::Dynamic, 1> param_3 = verticies_[3]->Parameters();
        Eigen::Matrix<double, 3, 1> Vj = param_3.head<3>();
        Eigen::Matrix<double, 3, 1> Baj = param_3.segment(3, 3);
        Eigen::Matrix<double, 3, 1> Bgj = param_3.tail<3>();

        residual_ = pre_integration_->evaluate(Pi, Qi, Vi, Bai, Bgi,
                                Pj, Qj, Vj, Baj, Bgj);
        SetInformation(pre_integration_->covariance.inverse());
    }

    void edgeimu::ComputeJacobians() {

        Eigen::Matrix<double, Eigen::Dynamic, 1> param_0 = verticies_[0]->Parameters();
        Eigen::Quaterniond Qi(param_0[6], param_0[3], param_0[4], param_0[5]);
        Eigen::Matrix<double, 3, 1> Pi = param_0.head<3>();

        Eigen::Matrix<double, Eigen::Dynamic, 1> param_1 = verticies_[1]->Parameters();
        Eigen::Matrix<double, 3, 1> Vi = param_1.head<3>();
        Eigen::Matrix<double, 3, 1> Bai = param_1.segment(3, 3);
        Eigen::Matrix<double, 3, 1> Bgi = param_1.tail<3>();

        Eigen::Matrix<double, Eigen::Dynamic, 1> param_2 = verticies_[2]->Parameters();
        Eigen::Quaterniond Qj(param_2[6], param_2[3], param_2[4], param_2[5]);
        Eigen::Matrix<double, 3, 1> Pj = param_2.head<3>();

        Eigen::Matrix<double, Eigen::Dynamic, 1> param_3 = verticies_[3]->Parameters();
        Eigen::Matrix<double, 3, 1> Vj = param_3.head<3>();
        Eigen::Matrix<double, 3, 1> Baj = param_3.segment(3, 3);
        Eigen::Matrix<double, 3, 1> Bgj = param_3.tail<3>();

        double sum_dt = pre_integration_->sum_dt;
        Eigen::Matrix3d dp_dba = pre_integration_->jacobian.template block<3, 3>(O_P, O_BA);
        Eigen::Matrix3d dp_dbg = pre_integration_->jacobian.template block<3, 3>(O_P, O_BG);

        Eigen::Matrix3d dq_dbg = pre_integration_->jacobian.template block<3, 3>(O_R, O_BG);

        Eigen::Matrix3d dv_dba = pre_integration_->jacobian.template block<3, 3>(O_V, O_BA);
        Eigen::Matrix3d dv_dbg = pre_integration_->jacobian.template block<3, 3>(O_V, O_BG);

        if (pre_integration_->jacobian.maxCoeff() > 1e8 || pre_integration_->jacobian.minCoeff() < -1e8)
        {
            std::cout << "numerical unstable in preintegration" << std::endl;
        }

            Eigen::Matrix<double, 15, 6, Eigen::RowMajor> jacobian_pose_i;
            jacobian_pose_i.setZero();

            jacobian_pose_i.block<3, 3>(O_P, O_P) = -Qi.inverse().toRotationMatrix();
            jacobian_pose_i.block<3, 3>(O_P, O_R) = Utility::skewSymmetric(Qi.inverse() * (0.5 * G * sum_dt * sum_dt + Pj - Pi - Vi * sum_dt));

            Eigen::Quaterniond corrected_delta_q1 = pre_integration_->delta_q * Utility::deltaQ(dq_dbg * (Bgi - pre_integration_->linearized_bg));
            jacobian_pose_i.block<3, 3>(O_R, O_R) = -(Utility::Qleft(Qj.inverse() * Qi) * Utility::Qright(corrected_delta_q1)).bottomRightCorner<3, 3>();

            jacobian_pose_i.block<3, 3>(O_V, O_R) = Utility::skewSymmetric(Qi.inverse() * (G * sum_dt + Vj - Vi));

            if (jacobian_pose_i.maxCoeff() > 1e8 || jacobian_pose_i.minCoeff() < -1e8)
            {
            std::cout << "numerical unstable in preintegration" << std::endl;
            }
            jacobians_[0] = jacobian_pose_i;
    //jacobians[1]

            Eigen::Matrix<double, 15, 9, Eigen::RowMajor> jacobian_speedbias_i;
            jacobian_speedbias_i.setZero();
            jacobian_speedbias_i.block<3, 3>(O_P, O_V - O_V) = -Qi.inverse().toRotationMatrix() * sum_dt;
            jacobian_speedbias_i.block<3, 3>(O_P, O_BA - O_V) = -dp_dba;
            jacobian_speedbias_i.block<3, 3>(O_P, O_BG - O_V) = -dp_dbg;
            jacobian_speedbias_i.block<3, 3>(O_R, O_BG - O_V) = -Utility::Qleft(Qj.inverse() * Qi * pre_integration_->delta_q).bottomRightCorner<3, 3>() * dq_dbg;


            jacobian_speedbias_i.block<3, 3>(O_V, O_V - O_V) = -Qi.inverse().toRotationMatrix();
            jacobian_speedbias_i.block<3, 3>(O_V, O_BA - O_V) = -dv_dba;
            jacobian_speedbias_i.block<3, 3>(O_V, O_BG - O_V) = -dv_dbg;

            jacobian_speedbias_i.block<3, 3>(O_BA, O_BA - O_V) = -Eigen::Matrix3d::Identity();

            jacobian_speedbias_i.block<3, 3>(O_BG, O_BG - O_V) = -Eigen::Matrix3d::Identity();

            jacobians_[1] = jacobian_speedbias_i;

    //jacobians[2]
            Eigen::Matrix<double, 15, 6, Eigen::RowMajor> jacobian_pose_j;
            jacobian_pose_j.setZero();

            jacobian_pose_j.block<3, 3>(O_P, O_P) = Qi.inverse().toRotationMatrix();

            Eigen::Quaterniond corrected_delta_q2 = pre_integration_->delta_q * Utility::deltaQ(dq_dbg * (Bgi - pre_integration_->linearized_bg));
            jacobian_pose_j.block<3, 3>(O_R, O_R) = Utility::Qleft(corrected_delta_q2.inverse() * Qi.inverse() * Qj).bottomRightCorner<3, 3>();

            jacobians_[2] = jacobian_pose_j;

    //jacobians[3]
            Eigen::Matrix<double, 15, 9, Eigen::RowMajor> jacobian_speedbias_j;
            jacobian_speedbias_j.setZero();

            jacobian_speedbias_j.block<3, 3>(O_V, O_V - O_V) = Qi.inverse().toRotationMatrix();

            jacobian_speedbias_j.block<3, 3>(O_BA, O_BA - O_V) = Eigen::Matrix3d::Identity();

            jacobian_speedbias_j.block<3, 3>(O_BG, O_BG - O_V) = Eigen::Matrix3d::Identity();

            jacobians_[3] = jacobian_speedbias_j;
    }
}
}