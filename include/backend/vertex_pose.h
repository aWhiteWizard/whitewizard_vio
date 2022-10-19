//vertex姿态数据类型,xyZ和四元数，7个维度。

/*
  需要维护[H|b]矩阵中的如下数据块
  p: pose, m:mappoint
 
     Hp1_p2    
     Hp2_p2    Hp2_m1    Hp2_m2    Hp2_m3     |    bp2
                         
                         Hm2_m2               |    bm2
                                   Hm2_m3     |    bm3
  1. 若该Camera为source camera，则维护vHessionSourceCamera；
  2. 若该Camera为measurement camera, 则维护vHessionMeasurementCamera；
  3. 并一直维护m_HessionDiagonal；
 */
#ifndef WHITEWIZARD_BACKEND_VERTEXPOSE_H
#define WHITEWIZARD_BACKEND_VERTEXPOSE_H

#include <memory>
#include "vertex.h"

namespace whitewizard_vio{
namespace backend{

    class VertexPose : public Vertex
    {
    private:
        /* data */
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        VertexPose() : Vertex(7,6) {};
        virtual void Plus(const Eigen::Matrix<double, Eigen::Dynamic, 1> &delta) override;
        std::string TypeInfo() const{
            return "vertex pose";
        }
    };

}
}

#endif
