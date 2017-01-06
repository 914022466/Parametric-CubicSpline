#include"spline.h"
#include<iostream>
using std::cout;
using std::vector;
int main()
{
	vector<MYPOINT> a;vector<MYPOINT>  ds;
	MYPOINT  ax1 = { 1,2,3 }; a.push_back(ax1);
	MYPOINT  ax2 = { 2,3,4 }; a.push_back(ax2);
	MYPOINT  ax3 = { 3,2,4 }; a.push_back(ax3);
	//MYPOINT  ax4 = { 3,2,8 }; a.push_back(ax4);
	//MYPOINT  ax5 = { 3,3,8 }; a.push_back(ax5);
	//MYPOINT  ax6 = { 4,4,9 }; a.push_back(ax6);

	Spline s;
	vector<MYPOINT>  re,re2,re3;
	s.setVertex(a,200);
	re=s.spline_create();
	

	return 0;
}