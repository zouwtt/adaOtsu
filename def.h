#ifndef DEF_H //就是头文件名（全大写后加个_H

#define DEF_H

#include <iostream>
#include <opencv.hpp>
#include <stack>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <fstream>
#include <string>
#define DTYPE double
struct PIX {
	int p;
	DTYPE Diff;
};
struct derivative {
	int p;
	DTYPE der;
};

struct NODE {
	double O;            /* old amount of obj. function up to prev. class */
	struct NODE* p_back; /* back pointer */
};

typedef struct _Edge
{
	double dDist;//Euclidean distance between two pixel points pt1 and pt2
	CvPoint pt1;
	CvPoint pt2;
	_Edge(double _dDist = 0, CvPoint _pt1 = cvPoint(0, 0), CvPoint _pt2 = cvPoint(0, 0))
		: dDist(_dDist), pt1(_pt1), pt2(_pt2) {}
}Edge;
struct CompareEdges {
	bool operator()(const Edge& edge1, const Edge& edge2) const {
		return edge1.dDist < edge2.dDist;
	}
};
struct ROWELEMENT {
	unsigned int col;
	struct ROWELEMENT* p_lastnext;
};
struct NODE** trellis;      /* trellis structure MxL */
DTYPE* N;                   /* numerator of objective function length: L */
DTYPE* D;                   /* denominator of objective function length: L */
DTYPE  Otemp, Ntemp, Dtemp;  /* temporary files */
int    stage;               /* variable for current stage */
unsigned int* maxcols;     /* needed by the smawkalgorithm to store the maxima positions */
unsigned int c;             /* current class */
using namespace std;
using namespace cv;
const int dx[4] = { 1, -1, 0, 0 }; // 4邻域的x方向偏移
const int dy[4] = { 0, 0, 1, -1 }; // 4邻域的y方向偏移
// 定义一个结构体来存储连通区域的信息
struct ConnectedComponent {
	cv::Mat mask;  // 连通区域的掩码
	int size;      // 连通区域的大小
};

#endif