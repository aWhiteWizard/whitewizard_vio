//xyz形式的顶点的vertex
#ifndef WHITEWIZARD_BACKEND_VERTEXPOINT_H
#define WHITEWIZARD_BACKEND_VERTEXPOINT_H

#include "vertex.h"

namespace whitewizard_vio{
namespace backend{
    class VertexPointXYZ : public Vertex {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        VertexPointXYZ() : Vertex(3) {}

        std::string TypeInfo() const { return "vertex point XYZ"; }
};
}
}
#endif