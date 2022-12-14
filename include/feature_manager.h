#ifndef WHITEWIZARD_FEATURE_MANAGER_H
#define WHITEWIZARD_FEATURE_MANAGER_H

#include <list>
#include <algorithm>
#include <vector>
#include <numeric>
#include <map>
#include <Eigen/Dense>

#include "parameter.h"

using namespace Eigen;
using namespace std;

class FeatureFrame
{
    private:

    public:
        double cur_td;
        Vector3d point;
        Vector2d uv;
        Vector2d velocity;
        double z;
        bool is_used;
        double parallax;
        MatrixXd A;
        VectorXd b;
        double dep_gradient;
        FeatureFrame(const Eigen::Matrix<double, 7, 1> &_point, double td)
        {
            point.x() = _point(0);
            point.y() = _point(1);
            point.z() = _point(2);
            uv.x() = _point(3);
            uv.y() = _point(4);
            velocity.x() = _point(5);
            velocity.y() = _point(6);
            cur_td = td;
        }
};

class FeatureId
{
    private:

    public:
        const int feature_id;
        int start_frame;
        vector<FeatureFrame> feature_per_frame;

        int used_num;
        bool is_outlier;
        bool is_margin;
        double estimated_depth;
        int solve_flag; // 0 haven't solve yet; 1 solve succ; 2 solve fail;

        Vector3d gt_p;

        FeatureId(int _feature_id, int _start_frame)
            : feature_id(_feature_id), start_frame(_start_frame),
            used_num(0), estimated_depth(-1.0), solve_flag(0)
        {
        }

        int endFrame();
};

class FeatureManager
{
    private:
        double compensatedParallax2(const FeatureId &it_per_id, int frame_count);
        const Matrix3d *Rs;
        Matrix3d ric[NUM_OF_CAM];
    public:
        FeatureManager(Matrix3d _Rs[]);

        void setRic(Matrix3d _ric[]);

        void clearState();

        int getFeatureCount();

        bool addFeatureCheckParallax(int frame_count, const map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> &image, double td);
        void debugShow();
        vector<pair<Vector3d, Vector3d>> getCorresponding(int frame_count_l, int frame_count_r);

        void setDepth(const VectorXd &x);
        void removeFailures();
        void clearDepth(const VectorXd &x);
        VectorXd getDepthVector();
        void triangulate(Vector3d Ps[], Vector3d tic[], Matrix3d ric[]);
        void removeBackShiftDepth(Eigen::Matrix3d marg_R, Eigen::Vector3d marg_P, Eigen::Matrix3d new_R, Eigen::Vector3d new_P);
        void removeBack();
        void removeFront(int frame_count);
        void removeOutlier();
        list<FeatureId> feature;
        int last_track_num;
};



#endif