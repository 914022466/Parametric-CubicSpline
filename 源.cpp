#include"spline.h"
#include<iostream>
using namespace std;
int main()
{
	vector<MYPOINT> a;vector<MYPOINT>  ds;
	MYPOINT  ax1 = { 1,2,3 }; a.push_back(ax1);
	MYPOINT  ax2 = { 2,3,4 }; a.push_back(ax2);
	MYPOINT  ax3 = { 3,2,4 }; a.push_back(ax3);
	//MYPOINT  ax4 = { 3,2,8 }; a.push_back(ax4);
	//MYPOINT  ax5 = { 3,3,8 }; a.push_back(ax5);
	//MYPOINT  ax6 = { 4,4,9 }; a.push_back(ax6);


	ds=convert2d(1.57,0,a,1,0,0);
	SPLINE s;
	vector<MYPOINT>  re,re2,re3;
	s.setVertex(a,97);
	re=s.spline_create();
	MYPOINT  ax7 = { 4,2,9 }; 
	re2=s.spine_modify(2,5,10,ax7);
	MYPOINT  ax8 = { 3,4,10 }; 
	re3=s.spine_modify(2,7,45,ax8);

	return 0;
}