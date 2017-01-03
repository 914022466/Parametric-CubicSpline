#pragma once
#include <vector>
#include <cmath>
#define INTERVAL 200

struct MYPOINT
{
	double x, y ,z;
	int run,wait;
};


//
std::vector<MYPOINT> convert2d(double temp_ry,double temp_rx,std::vector<MYPOINT> point,double DSPL_OVERRIDE,double ox,double oy);



class SPLINE {

public:
	std::vector<double> spline_u;
	//std::vector<std::vector<MYPOINT>>  spine_point_all;//曲线点
	MYPOINT begin_k1;//起点导矢
	MYPOINT end_k1;//终点导矢
	//空的构造函数，你可以稍后调用setVertex来设置插值点
	SPLINE()
	{

	}
	//设置插值点
	void setVertex(std::vector<MYPOINT> point, int interval_ = INTERVAL);
	//point_x插值点x坐标，point_y插值点y坐标

	//曲线点导矢
	std::vector<MYPOINT> daoshi;
	//获得点
	std::vector<MYPOINT> get_point();
	//移动插值点
	std::vector<MYPOINT> spine_modify(int begin, int end, int choice, MYPOINT new_position);
	std::vector<MYPOINT> spline_create(int i = 2);
	int interval;
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
	void canshuhua();

	void settheM(int condition);

	void getdaoshi();

	void Tridia(int n, int k, std::vector<std::vector<double>>&dia, std::vector<std::vector<double>>&t,
		int iflag);

	//////////////////////////////////////////////////////////////////////////
	//下列函数由于经常调用，故在此写作内联函数，牺牲空间换时间
	//////////////////////////////////////////////////////////////////////////


	void GetSplineDerivat(int i, double u, int j, MYPOINT &p)
	{
		int n = m_aVertex.size();
		double u1 = u;
		std::vector<double> tx(4);
		std::vector<double> ty(4);
		std::vector<double> tz(4);

		tx[0] = m_aVertex[i].x;
		tx[1] = m_Cx[i];
		tx[2] = 2 / delta[i] * (3 / delta[i] * (m_aVertex[i + 1].x - tx[0]) - 2 * tx[1] - m_Cx[i + 1]);
		tx[3] = 6 / delta[i] / delta[i] * (-2 / delta[i] * (m_aVertex[i + 1].x - tx[0]) + tx[1] + m_Cx[i + 1]);


		ty[0] = m_aVertex[i].y;
		ty[1] = m_Cy[i];
		ty[2] = 2 / delta[i] * (3 / delta[i] * (m_aVertex[i + 1].y - ty[0]) - 2 * ty[1] - m_Cy[i + 1]);
		ty[3] = 6 / delta[i] / delta[i] * (-2 / delta[i] * (m_aVertex[i + 1].y - ty[0]) + ty[1] + m_Cy[i + 1]);

		tz[0] = m_aVertex[i].z;
		tz[1] = m_Cz[i];
		tz[2] = 2 / delta[i] * (3 / delta[i] * (m_aVertex[i + 1].z - tz[0]) - 2 * tz[1] - m_Cz[i + 1]);
		tz[3] = 6 / delta[i] / delta[i] * (-2 / delta[i] * (m_aVertex[i + 1].z - tz[0]) + tz[1] + m_Cz[i + 1]);

		GetCubicCuvValue(i, tx, ty,tz, u1, j, p);
	}


	//////////////////////////////////////////////////////////////////////////
	//下列函数由于经常调用，故在此写作内联函数，牺牲空间换时间
	//////////////////////////////////////////////////////////////////////////
	void GetCubicCuvValue(int i, std::vector<double> &tx, std::vector<double> &ty, std::vector<double> &tz,double u,
		int j, MYPOINT &p)
	{
		double pjx = 0, pjy = 0,pjz=0;
		double du = u - canshu_u[i];
		double fj = 4 - j;
		for (int k = 3; k >= j; k--)
		{
			pjx = (pjx / fj)*du + tx[k];
			pjy = (pjy / fj)*du + ty[k];
			pjz= (pjz / fj)*du + tz[k];
			fj = fj - 1;
		}
		p.x = pjx;
		p.y = pjy;
		p.z = pjz;
	}

	//SPLINE(std::vector<double> point_x, std::vector<double> point_y, int interval_ = INTERVAL);

};

