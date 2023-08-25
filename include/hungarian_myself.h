#include<iostream>
#include<vector>
#include <chrono>
#include<Eigen/Dense>

using namespace std;
class Hungarian 
{
public:
	void solve(vector<vector<int>> CostMat, vector<int>& assignment);
private:
	int cost;
	int nRows, nCols, n;
	int uzr, uzc, step;
	Eigen::MatrixXd minR, minC;
	Eigen::MatrixXd dMat;
	Eigen::MatrixXd startZ;
	Eigen::MatrixXd coverColumn, coverRow, primeZ, cR, cC;
	vector<int> ridx, cidx;
	void init(vector<vector<int>> CostMat);
	void step1();
	void step2();
	void step3();
	void step4();
	void step5();
	void step6();
	Eigen::MatrixXd bsxfun(Eigen::MatrixXd dm, Eigen::MatrixXd mr, string mode);
	Eigen::MatrixXd isequal(Eigen::MatrixXd a, Eigen::MatrixXd b);
	Eigen::MatrixXd findx(Eigen::MatrixXd coverMat);
	void outerplus(vector<int>& rdx, vector<int>& cdx, int& minval);
	void findx2(Eigen::MatrixXd& tminr, Eigen::MatrixXd& tminc, Eigen::MatrixXd& dtMat);
};