#pragma once

//whitewizard create at sep 11th
#include <stdio.h>
#include <queue>
#include <map>
#include <mutex>
#include <thread>
#include <fstream>
#include <condition_variable>

#include <pangolin/pangolin.h>

#include <Eigen/Dense>
#include <Eigen/Core>

#include "opencv4/opencv2/opencv.hpp"
#include "opencv4/opencv2/highgui/highgui.hpp"
#include "opencv4/opencv2/imgproc/imgproc_c.h"
#include "opencv4/opencv2/core/core_c.h"
#include "opencv4/opencv2/core/types_c.h"

#include "estimator.h"
#include "parameter.h"
#include "feature_track.h"
//imu
struct IMU_MESSAGE
{
    double header;
    Eigen::Vector3d linear_acc;
    Eigen::Vector3d angular_vec;    
};
typedef std::shared_ptr<IMU_MESSAGE const> ImuConstPtr;
//这是关于typedef的用法，定义新类型，在IMU_MESSAGE类中不占据内存，只是定义了一个新类型，封装在了类中，便于使用。

//image
struct IMG_MESSAGE{
    double header;
    std::vector<Eigen::Vector3d> points;
    std::vector<int> point_id;
    std::vector<float> point_u;
    std::vector<float> point_v;
    std::vector<float> point_vel_x;
    std::vector<float> point_vel_y;
};
typedef std::shared_ptr<IMG_MESSAGE const> ImgConstPtr;

class System
{
    public:
        System(std::string sConfig_files);
        ~System();
        //发布数据
        void PubImgData(double dStampSec, cv::Mat &img);
        void PubImuData(double StampSec, const Eigen::Vector3d &Gyr, const Eigen::Vector3d &Acc);
        //后端线程
        void BackEnd();
        void Draw();
        void bPubImgData(double StampSec, const std::vector<cv::Point2f> &FeaturePoints);

        pangolin::OpenGlRenderState s_cam;
        pangolin::View d_cam;
    
    private:
        //特征跟踪
        std::vector<uchar> track_status;
        std::vector<float> track_error;

        FeatureTrack trackerData[NUM_OF_CAM];
        double first_image_time;
        int pub_count = 1;
        bool first_image_flag = true;
        double last_image_time = 0;
        bool init_pub = 0;

        //估计器
        Estimator estimator;
        std::condition_variable con;
        double current_time = -1;
        std::queue<ImuConstPtr> imu_buf;
        std::queue<ImgConstPtr> feature_buf;

        int sum_of_wait = 0;
        std::mutex m_buf;
        std::mutex m_state;
        std::mutex i_buf;
        std::mutex m_estimator;

        double latest_time;
        Eigen::Vector3d tmp_P;
        Eigen::Quaterniond tmp_Q;
        Eigen::Vector3d tmp_V;
        Eigen::Vector3d tmp_Ba;
        Eigen::Vector3d tmp_Bg;
        Eigen::Vector3d acc_0;
        Eigen::Vector3d gyr_0;

        bool init_feature = 0;
        bool init_imu = 1;
        double last_imu_t = 0;
        std::ofstream ofs_pose;
        std::vector<Eigen::Vector3d> vPath_to_draw;
        bool bStart_backend;
        std::vector<std::pair<std::vector<ImuConstPtr>, ImgConstPtr>> getMeasurements();

};

