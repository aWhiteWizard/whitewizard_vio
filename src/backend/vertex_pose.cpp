#include "../../include/backend/vertex_pose.h"
#include "../thirdparty/Sophus/sophus/se3.hpp"


namespace whitewizard_vio{
namespace backend{
    
    void VertexPose::Plus(const Eigen::Matrix<double, Eigen::Dynamic, 1> &delta) {
        Eigen::Matrix<double, Eigen::Dynamic, 1> &parameters = Parameters();
        parameters.head<3>() += delta.head<3>();
        Eigen::Quaterniond q(parameters[6], parameters[3], parameters[4], parameters[5]);
        q = q * Sophus::SO3d::exp(Eigen::Matrix<double, 3, 1> (delta[3], delta[4], delta[5])).unit_quaternion();  // right multiplication with so3
        q.normalized();
        parameters[3] = q.x();
        parameters[4] = q.y();
        parameters[5] = q.z();
        parameters[6] = q.w();
    }
}
}