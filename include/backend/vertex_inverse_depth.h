//逆深度形式的顶点

#ifndef WHITEWIZARD_BACKEND_VERTEX_INVERSE_H
#define WHITEWIZARD_BACKEND_VERTEX_INVERSE_H
#include "vertex.h"
namespace whitewizard_vio{
namespace backend{
    class VertexInverseDepth : public Vertex
    {
    private:
        /* data */
    public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    VertexInverseDepth() : Vertex(1) {}

    virtual std::string TypeInfo() const { return "vertex inverse depth"; }
    };
    
}
}
#endif