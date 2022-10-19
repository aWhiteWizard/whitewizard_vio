#pragma once


//对IMU的数据进行整合处理，包括噪声 ，数据发布
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <map>
#include "../thirdparty/Sophus/sophus/se3.hpp"

namespace whitewizard_vio{
namespace backend{
    class imu_dataset
    {
    private:
        // raw data from IMU
        std::vector<double> dt_buf_;
        std::vector<Eigen::Matrix<double, 3, 1>, Eigen::aligned_allocator<Eigen::Matrix<double, 3, 1>>>  acc_buf_;
        std::vector<Eigen::Matrix<double, 3, 1>, Eigen::aligned_allocator<Eigen::Matrix<double, 3, 1>>>  gyr_buf_;

        // pre-integrated IMU measurements
        double sum_dt_ = 0;
        Sophus::SO3d delta_r_;  // dR
        Eigen::Matrix<double, 3, 1> delta_v_ = Eigen::Matrix<double, 3, 1>::Zero();    // dv
        Eigen::Matrix<double, 3, 1> delta_p_ = Eigen::Matrix<double, 3, 1>::Zero();    // dp

        // gravity, biases
        static Eigen::Matrix<double, 3, 1> gravity_;
        Eigen::Matrix<double, 3, 1> bg_ = Eigen::Matrix<double, 3, 1>::Zero();    // initial bias of gyro
        Eigen::Matrix<double, 3, 1> ba_ = Eigen::Matrix<double, 3, 1>::Zero();    // initial bias of accelerator

        // jacobian w.r.t bg and ba
        Eigen::Matrix<double, 3, 3> dr_dbg_ = Eigen::Matrix<double, 3, 3>::Zero();
        Eigen::Matrix<double, 3, 3> dv_dbg_ = Eigen::Matrix<double, 3, 3>::Zero();
        Eigen::Matrix<double, 3, 3> dv_dba_ = Eigen::Matrix<double, 3, 3>::Zero();
        Eigen::Matrix<double, 3, 3> dp_dbg_ = Eigen::Matrix<double, 3, 3>::Zero();
        Eigen::Matrix<double, 3, 3> dp_dba_ = Eigen::Matrix<double, 3, 3>::Zero();

        // noise propagation
        Eigen::Matrix<double, 9, 9> covariance_measurement_ = Eigen::Matrix<double, 9, 9>::Zero();
        Eigen::Matrix<double, 6, 6> covariance_random_walk_ = Eigen::Matrix<double, 6, 6>::Zero();
        Eigen::Matrix<double, 9, 9> A_ = Eigen::Matrix<double, 9, 9>::Zero();
        Eigen::Matrix<double, 9, 6> B_ = Eigen::Matrix<double, 9, 6>::Zero();

        // raw noise of imu measurement
        Eigen::Matrix<double, 6, 6> noise_measurement_ = Eigen::Matrix<double, 6, 6>::Identity();
        Eigen::Matrix<double, 6, 6> noise_random_walk_ = Eigen::Matrix<double, 6, 6>::Identity();

        // accelerometer measurement noise standard deviation
        constexpr static double acc_noise_ = 0.2;
        // gyroscope measurement noise standard deviation
        constexpr static double gyr_noise_ = 0.02;
        // accelerometer bias random walk noise standard deviation
        constexpr static double acc_random_walk_ = 0.0002;
        // gyroscope bias random walk noise standard deviation
        constexpr static double gyr_random_walk_ = 2.0e-5;

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        /**
         * constructor, with initial bias a and bias g
         * @param ba
         * @param bg
         */
        explicit imu_dataset(const Eigen::Matrix<double, 3, 1> &ba, const Eigen::Matrix<double, 3, 1> &bg) : ba_(ba), bg_(bg) {
            const Eigen::Matrix<double, 3, 3> i3 = Eigen::Matrix<double, 3, 3>::Identity();
            noise_measurement_.block<3, 3>(0, 0) = (acc_noise_ * acc_noise_) * i3;
            noise_measurement_.block<3, 3>(3, 3) = (gyr_noise_ * gyr_noise_) * i3;
            noise_random_walk_.block<3, 3>(0, 0) = (acc_random_walk_ * acc_random_walk_) * i3;
            noise_random_walk_.block<3, 3>(3, 3) = (gyr_random_walk_ * gyr_random_walk_) * i3;
        }

        ~imu_dataset() {}

        /**
         * propage pre-integrated measurements using raw IMU data
         * @param dt
         * @param acc
         * @param gyr_1
         */
        void Propagate(double dt, const Eigen::Matrix<double, 3, 1> &acc, const Eigen::Matrix<double, 3, 1> &gyr);

        /**
         * according to pre-integration, when bias is updated, pre-integration should also be updated using
         * first-order expansion of ba and bg
         *
         * @param delta_ba
         * @param delta_bg
         */
        void Correct(const Eigen::Matrix<double, 3, 1> &delta_ba, const Eigen::Matrix<double, 3, 1> &delta_bg);

        void SetBiasG(const Eigen::Matrix<double, 3, 1> &bg) { bg_ = bg; }

        void SetBiasA(const Eigen::Matrix<double, 3, 1> &ba) { ba_ = ba; }

        /// if bias is update by a large value, redo the propagation
        void Repropagate();

        /// reset measurements
        /// NOTE ba and bg will not be reset, only measurements and jacobians will be reset!
        void Reset() {
            sum_dt_ = 0;
            delta_r_ = Sophus::SO3d();  // dR
            delta_v_ = Eigen::Matrix<double, 3, 1>::Zero();    // dv
            delta_p_ = Eigen::Matrix<double, 3, 1>::Zero();    // dp

            // jacobian w.r.t bg and ba
            dr_dbg_ = Eigen::Matrix<double, 3, 3>::Zero();
            dv_dbg_ = Eigen::Matrix<double, 3, 3>::Zero();
            dv_dba_ = Eigen::Matrix<double, 3, 3>::Zero();
            dp_dbg_ = Eigen::Matrix<double, 3, 3>::Zero();
            dp_dba_ = Eigen::Matrix<double, 3, 3>::Zero();

            // noise propagation
            covariance_measurement_ = Eigen::Matrix<double, 9, 9>::Zero();
            covariance_random_walk_ = Eigen::Matrix<double, 6, 6>::Zero();
            A_ = Eigen::Matrix<double, 9, 9>::Zero();
            B_ = Eigen::Matrix<double, 9, 6>::Zero();
        }

        /**
         * get the jacobians from r,v,p w.r.t. biases
         * @param _dr_dbg
         * @param _dv_dbg
         * @param _dv_dba
         * @param _dp_dbg
         * @param _dp_dba
         */
        void GetJacobians(Eigen::Matrix<double, 3, 3> &dr_dbg, Eigen::Matrix<double, 3, 3> &dv_dbg, Eigen::Matrix<double, 3, 3> &dv_dba, Eigen::Matrix<double, 3, 3> &dp_dbg, Eigen::Matrix<double, 3, 3> &dp_dba) const {
            dr_dbg = dr_dbg_;
            dv_dbg = dv_dbg_;
            dv_dba = dv_dba_;
            dp_dbg = dp_dbg_;
            dp_dba = dp_dba_;
        }

        Eigen::Matrix<double, 3, 3> GetDrDbg() const { return dr_dbg_; }

        // 得到噪声的协方差
        Eigen::Matrix<double, 9, 9> GetCovarianceMeasurement() const {
            return covariance_measurement_;
        }

        // 得到随机游走的协方差
        Eigen::Matrix<double, 6, 6> GetCovarianceRandomWalk() const {
            return noise_random_walk_ * sum_dt_;
        }

        // 得到时间
        double GetSumDt() const {
            return sum_dt_;
        }

        // 得到数据

        void GetDeltaRVP(Sophus::SO3d &delta_r, Eigen::Matrix<double, 3, 1> &delta_v, Eigen::Matrix<double, 3, 1> &delta_p) const {
            delta_r = delta_r_;
            delta_v = delta_v_;
            delta_p = delta_p_;
        }

        Eigen::Matrix<double, 3, 1> GetDv() const { return delta_v_; }

        Eigen::Matrix<double, 3, 1> GetDp() const { return delta_p_; }

        Sophus::SO3d GetDr() const { return delta_r_; }
    };
}
}