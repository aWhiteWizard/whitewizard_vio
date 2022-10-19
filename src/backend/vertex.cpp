#include "backend/vertex.h"
#include <iostream>

namespace whitewizard_vio{
namespace backend{

    unsigned long global_vertex_id = 0;

    Vertex::Vertex(int num_dimension, int local_dimension) {
        parameters_.resize(num_dimension, 1);
        local_dimension_ = local_dimension > 0 ? local_dimension : num_dimension;
        id_ = global_vertex_id++;

    }

    Vertex::~Vertex() {}

    int Vertex::Dimension() const {
        return parameters_.rows();
    }

    int Vertex::LocalDimension() const {
        return local_dimension_;
    }

    void Vertex::Plus(const Eigen::Matrix<double, Eigen::Dynamic, 1> &delta) {
        parameters_ += delta;
    }

}
}