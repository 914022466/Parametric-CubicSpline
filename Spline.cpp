#include "spline.h"

#include <iostream>
using  std::cout;
std::vector<MYPOINT> SPLINE::spline_create(int i)
{
	canshuhua();
	settheM(i);
	getdaoshi();
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
	 		GetSplineDerivat(j,u,0,c);
	 		MYPOINT daoshi_temp;
	 		GetSplineDerivat(j,u,1,daoshi_temp);
	 		all_point.push_back(c);
			daoshi.push_back(daoshi_temp);
			spline_u.push_back(u);
	 	}
		MYPOINT c;
		GetSplineDerivat(canshu_u.size()-2,1,0,c);
		MYPOINT daoshi_temp;
		GetSplineDerivat(canshu_u.size()-2,1,1,daoshi_temp);
		all_point.push_back(c);
		daoshi.push_back(daoshi_temp);
		spline_u.push_back(1);
		return all_point;
}

void SPLINE::getdaoshi()
{
	std::vector<std::vector<double>> daoshi = matrix_right;
	std::vector<std::vector<double>> m_tmep = matrix_m;
	Tridia(m_aVertex.size() - 1, 3, m_tmep, daoshi, 0);
	daoshi_vec = daoshi;
	for (int i = 0; i<m_aVertex.size(); i++)
	{
		m_Cx.push_back(daoshi_vec[i][0]);
		m_Cy.push_back(daoshi_vec[i][1]);
		m_Cz.push_back(daoshi_vec[i][2]);
	}
}

void SPLINE::Tridia(int n, int k, std::vector<std::vector<double>>&dia, std::vector<std::vector<double>>&t,
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



//保证all_point,daoshi正确

std::vector<MYPOINT> SPLINE::spine_modify(int begin, int end, int choice, MYPOINT new_position)
{
	//newpointlist下只有三个点
	std::vector<MYPOINT> newpointlist;
	newpointlist.push_back(all_point[begin]);
	newpointlist.push_back(new_position);
	newpointlist.push_back(all_point[end]);
	std::vector<MYPOINT> newsplinepoint,newtempdaoshi;
	SPLINE z;
	z.setVertex(newpointlist,(end-begin+1));
	MYPOINT bigi =daoshi[begin];
	MYPOINT en = daoshi[end];
	z.begin_k1 = bigi;
	z.end_k1 = en;
	newsplinepoint=z.spline_create(1);//生成局部曲线z
	newtempdaoshi=z.daoshi;

	double u;int j=0;
	for (int i = begin;i<=end; i++)
	{
		all_point[i]=newsplinepoint[j];
		daoshi[i]=newtempdaoshi[j];
		j++;
	}
	return all_point;

}


void SPLINE::settheM(int condition)
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

void SPLINE::canshuhua()
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
std::vector<MYPOINT> SPLINE::get_point()
{
	return all_point;
}

void SPLINE::setVertex(std::vector<MYPOINT> point, int interval_)
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



std::vector<MYPOINT> convert2d(double temp_ry,double temp_rx,std::vector<MYPOINT> point,double DSPL_OVERRIDE,double ox,double oy)
{
	std::vector<double> x,y,z;
	for (int i=0;i<point.size();i++)
	{
		x.push_back(point[i].x);
		y.push_back(point[i].y);
		z.push_back(point[i].z);
	}
	std::vector<MYPOINT> temp_p,p_screen;temp_p.resize(point.size());p_screen.resize(point.size());
	std::vector<double> d_x,d_y;d_x.resize(point.size());d_y.resize(point.size());
	for (int i=0;i<point.size();i++)
	{
		temp_p[i].y = y[i] * cos(temp_rx) - x[i] * sin(temp_rx);
		temp_p[i].x = y[i] * sin(temp_rx) + x[i] * cos(temp_rx);
		temp_p[i].z = z[i];
		p_screen[i].z = temp_p[i].z * cos(temp_ry) - temp_p[i].y * sin(temp_ry);
		p_screen[i].y = temp_p[i].z * sin(temp_ry) + temp_p[i].y * cos(temp_ry);
		p_screen[i].x = temp_p[i].x;

		d_x[i] = p_screen[i].x * DSPL_OVERRIDE + ox;
		/*DSPL_OVERRIDE 是缩放变量，ox和			 			oy是平面显示中心点坐标*/
		d_y[i] = oy - p_screen[i].y * DSPL_OVERRIDE;
	}
	std::vector<MYPOINT> result=point;
	for (int i=0;i<point.size();i++)
	{
		result[i].x=d_x[i];
		result[i].y=d_y[i];
	}
	return result;
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