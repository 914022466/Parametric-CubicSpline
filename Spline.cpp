#include "spline.h"




std::vector<MYPOINT> Spline::spline_create(int i)
{
	parameterization();
	setMatrix(i);
	getDerivative();
	 	for (double  u=0;u<1.0-0.8/(interval-1);u+=1.0/(interval-1))
	 	{
			int j=0;
			for (;j<canshu_u.size()-1;j++)
			{
				if (u>=canshu_u[j]&&u<canshu_u[j+1])
				{
					break;
				}
			}
	 		MYPOINT c;
	 		getSplineDerivat(j,u,0,c);
	 		MYPOINT daoshi_temp;
	 		getSplineDerivat(j,u,1,daoshi_temp);
	 		all_point.push_back(c);
			daoshi.push_back(daoshi_temp);
			spline_u.push_back(u);
	 	}
		MYPOINT c;
		getSplineDerivat(canshu_u.size()-2,1,0,c);
		MYPOINT daoshi_temp;
		getSplineDerivat(canshu_u.size()-2,1,1,daoshi_temp);
		all_point.push_back(c);
		daoshi.push_back(daoshi_temp);
		spline_u.push_back(1);
		return all_point;
}

void Spline::getDerivative()
{
	std::vector<std::vector<double>> daoshi = matrix_right;
	std::vector<std::vector<double>> m_tmep = matrix_m;
	tridia(m_aVertex.size() - 1, 3, m_tmep, daoshi, 0);
	daoshi_vec = daoshi;
	for (int i = 0; i<m_aVertex.size(); i++)
	{
		m_Cx.push_back(daoshi_vec[i][0]);
		m_Cy.push_back(daoshi_vec[i][1]);
		m_Cz.push_back(daoshi_vec[i][2]);
	}
}

void Spline::tridia(int n, int k, std::vector<std::vector<double>>&dia, std::vector<std::vector<double>>&t,
	int iflag)
{
	if (iflag == 0)
	{
		for (int i = 1; i <= n; i++)
		{
			dia[i][0] = -(dia[i][0] / dia[i - 1][1]);
			dia[i][1] = dia[i][1] + dia[i - 1][2] * dia[i][0];
		}
		for (int i = 1; i <= n; i++)
			for (int lk = 0; lk<k; lk++)
				t[i][lk] = t[i][lk] + t[i - 1][lk] * dia[i][0];
		for (int lk = 0; lk<k; lk++)
			t[n][lk] = t[n][lk] / dia[n][1];
		for (int i = n - 1; i >= 0; i--)
			for (int lk = 0; lk<k; lk++)
				t[i][lk] = (t[i][lk] - dia[i][2] * t[i + 1][lk]) / dia[i][1];
	}
	return;
}





void Spline::setMatrix(int condition)
{
	//condition=1，按导矢,condition=2，按自然边界
	matrix_m.resize(m_aVertex.size());
	matrix_right.resize(m_aVertex.size());
	//从1开始到倒数第二个
	for (int i = 1; i<m_aVertex.size() - 1; i++)
	{
		matrix_m[i].push_back(delta[i]);//ai=delta[i]
		matrix_m[i].push_back(2 * (delta[i] + delta[i - 1]));//bi=2*(ai+ci)
		matrix_m[i].push_back(delta[i - 1]);//ci=delta[i-1]
		matrix_right[i].push_back(3 * (delta[i] * _p[i - 1].x / delta[i - 1] + delta[i - 1] * _p[i].x / delta[i]));
		matrix_right[i].push_back(3 * (delta[i] * _p[i - 1].y / delta[i - 1] + delta[i - 1] * _p[i].y / delta[i]));
		matrix_right[i].push_back(3 * (delta[i] * _p[i - 1].z / delta[i - 1] + delta[i - 1] * _p[i].z / delta[i]));
	}
	if (condition == 1)
	{
		matrix_m[0].push_back(0); matrix_right[0].push_back(begin_k1.x);
		matrix_m[0].push_back(1); matrix_right[0].push_back(begin_k1.y);
		matrix_m[0].push_back(0); matrix_right[0].push_back(begin_k1.z);
		matrix_m[m_aVertex.size() - 1].push_back(0); matrix_right[m_aVertex.size() - 1].push_back(end_k1.x);
		matrix_m[m_aVertex.size() - 1].push_back(1); matrix_right[m_aVertex.size() - 1].push_back(end_k1.y);
		matrix_m[m_aVertex.size() - 1].push_back(0); matrix_right[m_aVertex.size() - 1].push_back(end_k1.z);
	}
	if (condition == 2)
	{
		matrix_m[0].push_back(0); matrix_right[0].push_back(3 * _p[0].x / delta[0]);
		matrix_m[0].push_back(2); matrix_right[0].push_back(3 * _p[0].y / delta[0]);
		matrix_m[0].push_back(1); matrix_right[0].push_back(3 * _p[0].z / delta[0]);
		matrix_m[m_aVertex.size() - 1].push_back(2); matrix_right[m_aVertex.size() - 1].push_back(3 * _p[m_aVertex.size() - 2].x / delta[m_aVertex.size() - 2]);
		matrix_m[m_aVertex.size() - 1].push_back(1); matrix_right[m_aVertex.size() - 1].push_back(_p[m_aVertex.size() - 2].y / delta[m_aVertex.size() - 2]);
		matrix_m[m_aVertex.size() - 1].push_back(0); matrix_right[m_aVertex.size() - 1].push_back(_p[m_aVertex.size() - 2].z / delta[m_aVertex.size() - 2]);


	}


}

void Spline::parameterization()
{
	double distance_sum = 0;//总弦长

	for (int i = 0; i<m_aVertex.size() - 1; i++)
	{
		double temp_dis_power = (m_aVertex[i].x - m_aVertex[i + 1].x)*(m_aVertex[i].x - m_aVertex[i + 1].x) + (m_aVertex[i].y - m_aVertex[i + 1].y)*(m_aVertex[i].y - m_aVertex[i + 1].y)+ (m_aVertex[i].z - m_aVertex[i + 1].z)*(m_aVertex[i].z - m_aVertex[i + 1].z);
		MYPOINT temp_point = { m_aVertex[i + 1].x - m_aVertex[i].x,m_aVertex[i + 1].y - m_aVertex[i].y,m_aVertex[i + 1].z - m_aVertex[i].z };
		_p.push_back(temp_point);
		double temp_dis = std::sqrt(temp_dis_power);
		distance_sum += temp_dis;
		distance.push_back(temp_dis);
		distance_sum_n.push_back(distance_sum);
	}
	canshu_u.push_back(0);//第一个为0
	for (int i = 0; i<distance_sum_n.size(); i++)
	{
		canshu_u.push_back(distance_sum_n[i] / distance_sum);
	}
	for (int i = 0; i<canshu_u.size() - 1; i++)
	{
		delta.push_back(canshu_u[i + 1] - canshu_u[i]);
	}
}
std::vector<MYPOINT> Spline::get_point()
{
	return all_point;
}

void Spline::setVertex(std::vector<MYPOINT> point, int interval_)
{
	for (int i = 0; i<point.size(); i++)
	{
		if (m_aVertex.size() > 0)
			if (m_aVertex[m_aVertex.size() - 1].x == point[i].x && m_aVertex[m_aVertex.size() - 1].y == point[i].y&&m_aVertex[m_aVertex.size() - 1].z == point[i].z)
				continue;
		m_aVertex.push_back(point[i]);
	}
	interval = interval_;
}


//////////////////////////////////////////////////////////////////////////
//下列函数由于经常调用，故在此写作内联函数，牺牲空间换时间
//////////////////////////////////////////////////////////////////////////

inline void Spline::getSplineDerivat(int i, double u, int j, MYPOINT &p)
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

	getCubicCuvValue(i, tx, ty,tz, u1, j, p);
}

inline void Spline::getCubicCuvValue(int i, std::vector<double> &tx, std::vector<double> &ty, std::vector<double> &tz,double u,
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


//SPLINE::SPLINE(std::vector<double> point_x, std::vector<double> point_y, int interval_)
//{
//	//不会在m_aVertex中出现两个相邻点相同
//	for (int i = 0; i<point_x.size(); i++)
//	{
//		if (m_aVertex.size()>0)
//			if (m_aVertex[m_aVertex.size() - 1].x == point_x[i] && m_aVertex[m_aVertex.size() - 1].y == point_y[i])//如果是同一个点就跳过
//				continue;
//		MYPOINT a = { point_x[i],point_y[i] };
//		m_aVertex.push_back(a);
//	}
//	interval = interval_;
//
//
//}