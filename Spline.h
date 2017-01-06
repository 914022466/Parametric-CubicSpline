#ifndef CUBICSPLINE
#define CUBICSPLINE

#include <vector>
#include <cmath>


//自定义的点类型
struct MYPOINT
{
	double x, y ,z;
	
};

//用于生成参数化三次插值样条曲线
class Spline {

public:

	//生成的曲线点
	std::vector<double> spline_u;
	
	MYPOINT begin_k1;//起点导矢
	MYPOINT end_k1;//终点导矢
	int interval;//生成曲线点的数量
public:
	//空的构造函数，你可以稍后调用setVertex来设置插值点
	Spline(){}
	//设置插值点
	void setVertex(std::vector<MYPOINT> point, int interval_ );

	//曲线点导矢
	std::vector<MYPOINT> daoshi;
	//获得点
	std::vector<MYPOINT> get_point();

	std::vector<MYPOINT> spline_create(int i = 2);
	
private:
	std::vector<MYPOINT> m_aVertex;//插值点,数目为num
	std::vector<double> distance;//弦长(num-1)
	std::vector<double> distance_sum_n;
	std::vector<double> canshu_u;//插值点的参数值(num)
	std::vector<double> delta;//delta(num-1)
	std::vector<MYPOINT> all_point;//曲线点

	std::vector<std::vector<double>>  matrix_m;
	std::vector<std::vector<double>>  matrix_right;

	//插值点导矢
	std::vector<std::vector<double>>  daoshi_vec;
	std::vector<MYPOINT> _p;

	std::vector<double> m_Cx;//每个插值点的导矢的x
	std::vector<double> m_Cy;//每个插值点的导矢的y
	std::vector<double> m_Cz;//每个插值点的导矢的z

private:
	//规范积累弦长参数化
	void parameterization();

	void setMatrix(int condition);

	void getDerivative();

	void tridia(int n, int k, std::vector<std::vector<double>>&dia, std::vector<std::vector<double>>&t,
		int iflag);

	void getSplineDerivat(int i, double u, int j, MYPOINT &p);

	void getCubicCuvValue(int i, std::vector<double> &tx, std::vector<double> &ty, std::vector<double> &tz,double u,
		int j, MYPOINT &p);

};

#endif