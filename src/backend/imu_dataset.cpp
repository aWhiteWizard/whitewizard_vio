#include "backend/imu_dataset.h"

using Sophus::SO3d;

namespace whitewizard_vio{
namespace backend{

    void imu_dataset::Propagate(double dt, const Eigen::Matrix<double, 3, 1> &acc, const Eigen::Matrix<double, 3, 1> &gyr){
        dt_buf_.emplace_back(dt);
        acc_buf_.emplace_back(acc);
        gyr_buf_.emplace_back(gyr);

        Sophus::SO3d dR = Sophus::SO3d::exp((gyr - bg_) * dt);
        delta_r_ = delta_r_ * dR;
        delta_v_ += delta_r_ * (acc - ba_) * dt;
        delta_p_ += delta_v_ * dt + 0.5 * (delta_r_ * (acc - ba_) * dt * dt);
        sum_dt_ += dt;

        // update jacobians w.r.t. bg and ba
        dr_dbg_ -= delta_r_.inverse().matrix() * SO3d::JacobianR(((gyr - bg_) * dt)) * dt;
        dv_dba_ -= delta_r_.matrix() * dt;
        dv_dbg_ -= delta_r_.matrix() * SO3d::hat(acc - ba_) * dr_dbg_ * dt;
        dp_dba_ += dv_dba_ * dt - 0.5 * delta_r_.matrix() * dt * dt;
        dp_dbg_ += dv_dbg_ * dt - 0.5 * delta_r_.matrix() * SO3d::hat(acc - ba_) * dr_dbg_ * dt * dt;

        // propagate noise
        A_.block<3, 3>(0, 0) = dR.inverse().matrix();
        B_.block<3, 3>(0, 0) = SO3d::JacobianR(dR.log());

        A_.block<3, 3>(3, 0) = -delta_r_.matrix() * SO3d::hat(acc - ba_) * dt;
        A_.block<3, 3>(3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
        B_.block<3, 3>(3, 3) = delta_r_.matrix() * dt;

        A_.block<3, 3>(6, 0) = -0.5 * delta_r_.matrix() * SO3d::hat(acc - ba_) * dt * dt;
        A_.block<3, 3>(6, 3) = Eigen::Matrix<double, 3, 3>::Identity() * dt;
        A_.block<3, 3>(6, 6) = Eigen::Matrix<double, 3, 3>::Identity();
        B_.block<3, 3>(6, 3) = 0.5 * delta_r_.matrix() * dt * dt;

        covariance_measurement_ = A_ * covariance_measurement_ * A_.transpose() + B_ * noise_measurement_ * B_.transpose();
    }

    void imu_dataset::Repropagate(){
        auto dt = dt_buf_;
        auto acc_buf = acc_buf_;
        auto gyr_buf = gyr_buf_;
        Reset();

        for (size_t i = 0; i < dt.size(); ++i) {
            Propagate(dt[i], acc_buf[i], gyr_buf[i]);
        }
    }

    void imu_dataset::Correct(const Eigen::Matrix<double, 3, 1> &delta_ba, const Eigen::Matrix<double, 3, 1> &delta_bg){
        delta_r_ = delta_r_ * SO3d::exp(dr_dbg_ * delta_bg);
        delta_v_ += dv_dba_ * delta_ba + dv_dbg_ * delta_bg;
        delta_p_ += dp_dba_ * delta_ba + dp_dbg_ * delta_bg;
    }
}
}