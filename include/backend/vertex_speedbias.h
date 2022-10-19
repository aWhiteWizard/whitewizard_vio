//vertex速度偏置，九个维度
#ifndef WHITEWIZARD_BACKEND_VERTEXSPEEDBIAS_H
#define WHITEWIZARD_BACKEND_VERTEXSPEEDBIAS_H

#include <memory>
#include "vertex.h"

namespace whitewizard_vio{
namespace backend{
    class VertexSpeedBias : public Vertex {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        VertexSpeedBias() : Vertex(9) {}

        std::string TypeInfo() const {
            return "VertexSpeedBias";
        }

    };
}
}

#endif