#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
    * Calculate the RMSE here.
  */
    VectorXd rmse(4);
    rmse << 0.0, 0.0, 0.0, 0.0;
    // 输入有效性验证
    if (estimations.size() != ground_truth.size() || estimations.size()==0 ){
        cout << "Invalid input data" << endl;
        return rmse;
    }
	//accumulate squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i){
    	VectorXd res = estimations[i] - ground_truth[i];
		// coefficient-wise 乘法
        // more to see: https://eigen.tuxfamily.org/dox/group__TutorialArrayClass.html
		res = res.array()*res.array();
		rmse += res;
	}
	// 平均值
	rmse = rmse/estimations.size();
	// 计算平方根
	rmse = rmse.array().sqrt();

	return rmse;
}
