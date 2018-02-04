#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

double UKF::NormalizeAngle(double angle) {
    while(angle > M_PI) angle -= 2.0*M_PI;
    while(angle <-M_PI) angle += 2.0*M_PI;
    return angle;
}

UKF::UKF() {
  // 是否使用lidar/radar测量数据
  use_laser_ = true;
  use_radar_ = true;

  // 状态方程初始化
  x_ = VectorXd::Zero(5);

  // initial covariance matrix
  P_ = MatrixXd::Zero(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.7;
  
  /* --------------------------传感器内置参数 DO NOT MODIFY ----------------------------------*/
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  /*-----------------------------------------------------------------------------------------*/


  time_us_ = 0.0;

  // 状态方程的维度
  n_x_ = 5;

  // 增强后状态方程的维度
  n_aug_ = 7;
  
  // for reuse, set variable
  n_sig = 2*n_aug_+1;
  
  // 扩散参数
  lambda_ = 3 - n_aug_;
  
  // 权重向量
  weights_ = VectorXd::Zero(n_sig);
  weights_.fill(0.5 / (n_aug_ + lambda_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  // 预测的σ点矩阵
  Xsig_pred_ = MatrixXd::Zero(n_x_, n_sig);

  // 测量噪声的协方差矩阵-radar / lidar
  int n_z = 3;
  R_radar = MatrixXd::Zero(n_z,n_z);
  R_radar << std_radr_ * std_radr_, 0, 0,
             0, std_radphi_ * std_radphi_, 0,
             0, 0, std_radrd_ * std_radrd_;

  R_laser = MatrixXd::Zero(2, 2);
  R_laser << std_laspx_ * std_laspx_, 0,
             0, std_laspy_ * std_laspy_;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  // 初始化
  if (!is_initialized_) {
    // LADAR 测量数据
    if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
      double p_x = meas_package.raw_measurements_(0);
      double p_y = meas_package.raw_measurements_(1);
    x_ << p_x, p_y, 0, 0, 0;
    }
    // RADAR 测量数据
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
      double rho = meas_package.raw_measurements_(0);
      double angle = meas_package.raw_measurements_(1);
      double p_x = rho * cos(angle);
      double p_y = rho * sin(angle);

      x_ << p_x, p_y, 0, 0, 0;
    }

    P_ << 1, 0., 0., 0., 0.,
          0., 1, 0., 0., 0.,
          0., 0., 1, 0., 0.,
          0., 0., 0., 1, 0.,
          0., 0., 0., 0., 1;
  
    // 前一时间步
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;

    return;
  }

  // 预测
  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0; // 单位: s
  time_us_ = meas_package.timestamp_;
  Prediction(delta_t);

  // 更新
  if ((meas_package.sensor_type_ == MeasurementPackage::LASER) && use_laser_ ) {
    UpdateLidar(meas_package);
  } else if ((meas_package.sensor_type_ == MeasurementPackage::RADAR) && use_radar_) {
    UpdateRadar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 */
void UKF::Prediction(double delta_t) {
  // 估计目标位置
  /*---------------------------------------------------------------*/
  /*                           1.生成σ点                            */
  /*---------------------------------------------------------------*/

  // 创建增强后的状态矩阵
  VectorXd X_aug = VectorXd::Zero(n_aug_);
  X_aug.head(5) = x_;
  X_aug(5) = 0;
  X_aug(6) = 0;

  // 增强的协方差矩阵
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);


  MatrixXd Q = MatrixXd::Zero(2,2);
  Q << std_a_ * std_a_, 0,
       0, std_yawdd_ * std_yawdd_;

  P_aug.topLeftCorner(5,5) = P_;
  P_aug.bottomRightCorner(2,2) = Q;

  // 矩阵平方根
  MatrixXd L = P_aug.llt().matrixL();

  // σ噪点矩阵
  MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, n_sig);
  Xsig_aug.col(0) = X_aug;

  for(int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i+1) = X_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(n_aug_+i+1) = X_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  /*---------------------------------------------------------------*/
  /*              2.根据增强的σ点计算预测的状态矩阵                   */
  /*---------------------------------------------------------------*/

  // 预测的σ点集
  for (int i = 0; i < n_sig; i++)
  {
    // 提取状态向量中的变量
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double u = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yaw_d = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yaw_dd = Xsig_aug(6,i);

    double p_x_pred, p_y_pred;

    // 注意除法
    if (fabs(yaw_d) > 0.0001) {
      p_x_pred = p_x + u / yaw_d * ( sin (yaw + yaw_d * delta_t) - sin(yaw));
      p_y_pred = p_y + u / yaw_d * ( cos(yaw) - cos(yaw + yaw_d * delta_t) );
    }
    else {
      p_x_pred = p_x + u * delta_t * cos(yaw);
      p_y_pred = p_y + u * delta_t * sin(yaw);
    }

    double u_pred = u;
    double yaw_pred = yaw + yaw_d * delta_t;
    double yaw_d_pred = yaw_d;

    // 添加噪声
    p_x_pred += 0.5 * nu_a * (delta_t * delta_t) * cos(yaw);
    p_y_pred += 0.5 * nu_a * (delta_t * delta_t) * sin(yaw);
    u_pred += nu_a * delta_t;
    yaw_pred += 0.5 * nu_yaw_dd * (delta_t * delta_t);
    yaw_d_pred += nu_yaw_dd * delta_t;

    // 得到添加噪声的预测状态矩阵
    Xsig_pred_(0, i) = p_x_pred;
    Xsig_pred_(1, i) = p_y_pred;
    Xsig_pred_(2, i) = u_pred;
    Xsig_pred_(3, i) = yaw_pred;
    Xsig_pred_(4, i) = yaw_d_pred;
  }

  /*---------------------------------------------------------------*/
  /*              3.计算预测的平均值与协方差                         */
  /*---------------------------------------------------------------*/

  // 预测状态向量
  x_.fill(0.0);
  x_ = Xsig_pred_ * weights_;

  // 状态协方差矩阵
  P_.fill(0.0);
  
  for (int i = 0; i < n_sig; i++) {
    // 状态矩阵的差分
    VectorXd diff = Xsig_pred_.col(i) - x_;

    // 极坐标角度的正则化
    diff(3) = NormalizeAngle(diff(3));

    P_ += weights_(i) * diff * diff.transpose();
  }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  
  // LIDAR的测量值
  VectorXd z = meas_package.raw_measurements_;
  
  MatrixXd H_laser = MatrixXd::Zero(2, 5);
  H_laser << 1, 0, 0, 0, 0,
             0, 1, 0, 0, 0;

  // 计算卡尔曼增益矩阵
  VectorXd z_pred = H_laser * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_laser.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = H_laser * PHt + R_laser;
  MatrixXd Si = S.inverse();
  MatrixXd K =  PHt * Si;

  // 更新x_ 与 P_
  x_ += (K * y);
  // P_ = (I - K * H_laser) * P_ I为单位矩阵
  P_ -= K * H_laser * P_;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  // 测量数据的维度 曲率半径，角度，角速度
  int n_z = 3;

  /*---------------------------------------------------------------*/
  /*              1. 结合噪声矩阵进行极坐标预测                      */
  /*---------------------------------------------------------------*/
  // 创建噪声矩阵
  MatrixXd Zsig = MatrixXd::Zero(n_z, 2 * n_aug_ + 1);

  // σ点矩阵转换
  for (int i = 0; i < n_sig; i++) {  // 2n+1个σ点
    double p_x  = Xsig_pred_(0,i);
    double p_y  = Xsig_pred_(1,i);
    double u   = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double u_x = u * cos(yaw);
    double u_y = u * sin(yaw);

    // 测量矩阵转换
    Zsig(0,i) = sqrt(p_x * p_x + p_y * p_y);  //曲率半径
    if(Zsig(0,i)<0.0001)
      Zsig(0,i) = 0.0001;
    
    if ((p_y == 0.0) && (p_x == 0.0))
      Zsig(1,i) = atan2(0.0001, 0.0);
    else
      Zsig(1,i) = atan2(p_y,p_x);   // 角度
    
    Zsig(2,i) = (p_x*u_x + p_y*u_y ) / Zsig(0,i);
  }

  // 预测的测量向量均值
  VectorXd z_pred = VectorXd::Zero(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < n_sig; i++) {
      z_pred += weights_(i) * Zsig.col(i);
  }

  // 测量的协方差矩阵 S
  MatrixXd S = MatrixXd::Zero(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < n_sig; i++) {  // 2n+1σ点
    // 测量残差
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = NormalizeAngle(z_diff(1));

    S += weights_(i) * z_diff * z_diff.transpose();
  }
  
  S += R_radar;

  /*---------------------------------------------------------------*/
  /*       2. 根据交叉协方差矩阵与卡尔曼增益进行状态更新               */
  /*---------------------------------------------------------------*/
  //交叉协方差矩阵 Tc
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);

  // 当前RADAR测量向量
  VectorXd z = VectorXd::Zero(n_z);
  z = meas_package.raw_measurements_;

  // 计算Tc，与计算S过程大致相同
  Tc.fill(0.0);
  for (int i = 0; i < n_sig; i++) {

    // 计算测量-预测残差，并对角度进行正则化
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = NormalizeAngle(z_diff(1));

    // 状态矩阵残差
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = NormalizeAngle(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // 卡尔曼增益矩阵
  MatrixXd K = Tc * S.inverse();

  // 残差
  VectorXd z_diff = z - z_pred;

  // 角度的正则化
  z_diff(1) = NormalizeAngle(z_diff(1));

  // 状态更新
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
}