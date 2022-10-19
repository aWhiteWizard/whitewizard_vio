
#ifndef WHITEWIZARD_BACKEND_REPRO_EDGE_H
#define WHITEWIZARD_BACKEND_REPRO_EDGE_H


#include <memory>
#include <string>

#include <Eigen/Dense>

#include "edge.h"

namespace whitewizard_vio{
namespace backend{
    /*
    视觉重投影误差的边，对于逆深度，三元边，相连的顶点有
    1、路标点的逆深度InveseDepth
    2、第一次观测到该路标点的source Camera的位姿T_World_From_Body1，
    3、观测到该路标点的mearsurement Camera位姿T_World_From_Body2。
    并且，verticies_顶点顺序必须为InveseDepth、T_World_From_Body1、T_World_From_Body2。
    */
    class edgereprojection : public edge
    {
    private:
        Eigen::Matrix<double, 3, 1> pts_i_, pts_j_;
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    edgereprojection(const Eigen::Matrix<double, 3, 1> &pts_i, const Eigen::Matrix<double, 3, 1> &pts_j)
        : edge(2, 4, std::vector<std::string>{"VertexInverseDepth", "VertexPose", "VertexPose", "VertexPose"}) {
        pts_i_ = pts_i;
        pts_j_ = pts_j;
    }

    /// 返回边的类型信息
    virtual std::string TypeInfo() const override { return "EdgeReprojection"; }

    /// 计算残差
    virtual void ComputeResidual() override;

    /// 计算雅可比
    virtual void ComputeJacobians() override;

    };
    /*
    视觉重投影误差的边，对于世界坐标系XYZ，二元边，相连的顶点有
    1、路标点的世界坐标系XYZ
    2、观测到该路标点的 Camera 的位姿T_World_From_Body1。
    并且，verticies_顶点顺序必须为 XYZ、T_World_From_Body1。
    */
   class edgereprojectionXYZ : public edge
   {
   private:
        //从IMU转换到cam
        Eigen::Quaterniond qic;
        Eigen::Matrix<double, 3, 1> tic;

        //测量值
        Eigen::Matrix<double, 3, 1> obs_;
   public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        edgereprojectionXYZ(const Eigen::Matrix<double, 3, 1> &pts_i)
        : edge(2, 2, std::vector<std::string>{"VertexXYZ", "VertexPose"}) {
        obs_ = pts_i;
    }

    /// 返回边的类型信息
    virtual std::string TypeInfo() const override { return "edgereprojectionXYZ"; }

    /// 计算残差
    virtual void ComputeResidual() override;

    /// 计算雅可比
    virtual void ComputeJacobians() override;

    void SetTranslationImuFromCamera(Eigen::Quaterniond &qic_, Eigen::Matrix<double, 3, 1> &tic_);
   };
    /*
    仅计算重投影pose
    */
   class edgereprojectionpose : public edge
   {
   private:
    Eigen::Matrix<double, 3, 1> landmark_world_;
    Eigen::Matrix<double, 3, 3> K_;
   public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    edgereprojectionpose(const Eigen::Matrix<double, 3, 1> &landmark_world, const Eigen::Matrix<double, 3, 3> &K) :
        edge(2, 1, std::vector<std::string>{"Vertexpose"}),
        landmark_world_(landmark_world), K_(K) {}

    /// 返回边的类型信息
    virtual std::string TypeInfo() const override { return "pose"; }

    /// 计算残差
    virtual void ComputeResidual() override;

    /// 计算雅可比
    virtual void ComputeJacobians() override;
   };
}
}

#endif