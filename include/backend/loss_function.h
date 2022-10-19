//提供各个损失函数的基类
#ifndef WHITEWIZARD_LOSS_FUNCTION_H
#define WHITEWIZARD_LOSS_FUNCTION_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <map>

namespace whitewizard_vio{
    namespace backend{
        /*
            计算error的比例因子:
            error是e^T
            输出：
            rho[0]: 实际的误差值
            rho[0]: 一阶导数
            rho[0]: 二阶导数
        */
       class LossFunction
       {
       private:
        /* data */
       public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        virtual ~LossFunction() {}

        virtual void Compute(double err2, Eigen::Vector3d& rho) const=0;
       };

       /*
            一般的损失函数
            使用空指针和loss function效果相同,override
       */
      class NormalLoss : public LossFunction{
        virtual void Compute(double err2, Eigen::Vector3d& rho) const override 
        {
            // TODO:: whether multiply 1/2
            rho[0] = err2;
            rho[1] = 1;
            rho[2] = 0;
        }
      };

      /*
            Huber损失函数
            Huber(e) = e^2                  if e <= delta
            Huber(e) = delta*(2*e - delta)  if e > delta
      */
        class HuberLoss : public LossFunction
        {
        private:
            double delta_;
        public:
            explicit HuberLoss(double delta) : delta_(delta) {}

            virtual void Compute(double err2, Eigen::Vector3d& rho) const override;
        };
        /*
            柯西损失函数
        */
        class CauchyLoss : public LossFunction
        {
            private:
                double delta_;
            public:
                explicit CauchyLoss(double delta) : delta_(delta) {}

                virtual void Compute(double err2, Eigen::Vector3d& rho) const override;
        };
        /*
            Tukey损失函数
        */
        class TukeyLoss : public LossFunction
        {
            private:
                double delta_;
            public:
                explicit TukeyLoss(double delta) : delta_(delta) {}

                virtual void Compute(double err2, Eigen::Vector3d& rho) const override;
        };       
    }
}

#endif