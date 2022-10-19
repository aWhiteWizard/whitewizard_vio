#include "backend/vertex.h"
#include "backend/edge.h"
#include <iostream>

using namespace std;

namespace whitewizard_vio{
namespace backend{
    unsigned long global_edge_id = 0;

    edge::edge(int residual_dim, int num_verticies, const std::vector<std::string> &verticies_types){
        residual_.resize(residual_dim, 1);
        if(!verticies_types.empty())
            verticies_types_ = verticies_types;
        jacobians_.resize(num_verticies);
        id_ = global_edge_id++;

        Eigen::MatrixXd info_matrix(residual_dim, residual_dim);
        info_matrix.setIdentity();
        info_matrix_ = info_matrix;
        lossfunction_ = NULL;
    }

    edge::~edge(){}

    double edge::Chi2() const{
        return residual_.transpose() * info_matrix_ * residual_;
    }

    double edge::RobustChi2() const{
        double e2 = this->Chi2();
        if(lossfunction_)
        {
            Eigen::Vector3d rho;
            lossfunction_->Compute(e2,rho);
            e2 = rho[0];
        }
        return e2;
    }
    void edge::RobustInfo(double &drho, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &info)
    const {
        if(lossfunction_)
        {
            double e2 = this -> Chi2();
            Eigen::Vector3d rho;
            lossfunction_->Compute(e2, rho);
            Eigen::Matrix<double, Eigen::Dynamic, 1> weight_err = sqrt_info_matrix_ * residual_;

            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> robust_info(info_matrix_.rows(), info_matrix_.cols());
            robust_info.setIdentity();
            robust_info *= rho[1];
            if(rho[1] + 2*rho[2]*e2 > 0.)
            {
                robust_info += 2 * rho[2] * weight_err * weight_err.transpose();
            }
            info = robust_info * info_matrix_;
            drho = rho[1];
        }else
        {
            drho = 1.0;
            info = info_matrix_;
        }
    }
    bool edge::CheckValid() {
        if (!verticies_types_.empty()) {
            // check type info
            for (size_t i = 0; i < verticies_.size(); ++i) {
                if (verticies_types_[i] != verticies_[i]->TypeInfo()) {
                    cout << "Vertex type does not match, should be " << verticies_types_[i] <<
                        ", but set to " << verticies_[i]->TypeInfo() << endl;
                    return false;
                }
            }
        }
        return true;
    }
}
}