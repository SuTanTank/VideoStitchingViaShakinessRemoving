#include <opencv2/opencv.hpp>
#include <iostream>
#include <cassert>
#include <vector>

#include "maxflow/graph.h"

#define INT_MAX 1000000


using namespace std;
using namespace cv;

typedef Graph<int,int,int> GraphType;

void getEdgeMap(int edgeThre, Mat &rgb, Mat &edge)
{
	cvtColor(rgb, edge, CV_BGR2GRAY);
	// ver 1
	//blur(edge, edge, Size(3, 3));
	//Canny(edge, edge, edgeThre, edgeThre * 3);
	//blur(edge, edge, Size(30, 30));

	// ver 2
	blur(edge, edge, Size(3, 3));
	Canny(edge, edge, edgeThre, edgeThre * 3);
	blur(edge, edge, Size(30, 30));

	return;
}

bool isBlack(Vec3b pixel)
{
	if (pixel[0] == 0 && pixel[1] == 0 && pixel[2] == 0)
		return true;
	else
		return false;
}

Mat FillHoles(const Mat& src)
{
	//assume input is uint8 B & W (0 or 1)
	//this function imitates imfill(image,'hole')
	Mat src1 = src.clone();
	
	Mat dst = Mat::zeros(src.rows, src.cols, CV_8UC1);

	std::vector<std::vector<Point> > contours;
	std::vector<Vec4i> hierarchy;

	findContours(src1, contours, hierarchy, CV_RETR_CCOMP, CV_CHAIN_APPROX_SIMPLE);

	drawContours(dst, contours, 0, Scalar::all(255), CV_FILLED);

	return dst;
}

Mat FindBoundaries(Mat src)
{
	Mat src1 = src.clone();
	Mat dst = Mat::zeros(src.size(), CV_8UC1);
	Mat eroded = Mat::zeros(src.size(), CV_8UC1);

	erode(src1, eroded,Mat());
	absdiff(src1, eroded, dst);

	return dst;	
}

void process_err(char *errInfo)
{
	printf_s("%s\n", errInfo);
}

void findSeamAndBlend(int stable, int start, int end, int ab, float scale)
{
	char srcFolder[] = ".";
	char trgFolder[] = ".";
	char tmpFolder[] = ".";


	char filename[255];
	int edgeThre = 50;
	Mat A0, B0, A1, B1;
	Mat A0s, B0s, A1s, B1s;
	Mat A0e, B0e, A1e, B1e;
	Mat edge0, edge1;
	Mat A0sMask, B0sMask, A1sMask, B1sMask;

	sprintf(filename, "%s/A%d.png", srcFolder, start);
	A1 = imread(filename);
	resize(A1, A1s, Size(), scale, scale);
	getEdgeMap(50, A1s, A1e);
	sprintf(filename, "%s/B%d.png", srcFolder, start);
	B1 = imread(filename);
	resize(B1, B1s, Size(), scale, scale);
	getEdgeMap(50, B1s, B1e);

	Mat A1sMasktmp = Mat::zeros(A1s.size(), CV_8UC1);
	Mat B1sMasktmp = Mat::zeros(B1s.size(), CV_8UC1);

	cvtColor(A1s, A1sMasktmp, CV_BGR2GRAY);
	cvtColor(B1s, B1sMasktmp, CV_BGR2GRAY);

	A1sMask = FillHoles(A1sMasktmp);
	B1sMask = FillHoles(B1sMasktmp);

	int width = (int)A1s.cols;
	int height = (int)A1s.rows;
	//Mat mask(width, height, CV_8U, 0);
	Mat mask = Mat::zeros(height, width, CV_8UC1);
	Mat saveMask = Mat::zeros(height, width, CV_8UC1);
	int est_nodes = height * width * 2;
	int est_edges = est_nodes * 5;
	//GraphType g0(est_nodes, est_edges);

	const int PRESET_WEIGHT = 1;

	for (int frameIndex = start; frameIndex < end; frameIndex++)
	{
		A0 = A1;
		A0s = A1s;
		A0e = A1e;
		B0 = B1;
		B0s = B1s;
		B0e = B1e;

		A0sMask = A1sMask;
		B0sMask = B1sMask;

		sprintf(filename, "%s/A%d.png", srcFolder, frameIndex + 1);
		A1 = imread(filename);
		resize(A1, A1s, Size(), scale, scale);
		getEdgeMap(50, A1s, A1e);
		sprintf(filename, "%s/B%d.png", srcFolder, frameIndex + 1);
		B1 = imread(filename);
		resize(B1, B1s, Size(), scale, scale);
		getEdgeMap(50, B1s, B1e);

		A1sMasktmp = Mat::zeros(A1s.size(), CV_8UC1);
		B1sMasktmp = Mat::zeros(B1s.size(), CV_8UC1);

		cvtColor(A1s, A1sMasktmp, CV_BGR2GRAY);
		cvtColor(B1s, B1sMasktmp, CV_BGR2GRAY);

		A1sMask = FillHoles(A1sMasktmp);
		B1sMask = FillHoles(B1sMasktmp);


		Mat imgAll = Mat::zeros(A0sMask.size(), CV_8UC1);
		bitwise_or(A0sMask, B0sMask, imgAll);

		//imwrite("D:/test.png", imgAll);


		edge0 = A0e + B0e;
		edge1 = A1e + B1e;

		//sprintf_s(filename, "%s/edge0.png", tmpFolder);
		//imwrite(filename, edge0);
		//sprintf_s(filename, "%s/edge1.png", tmpFolder);
		//imwrite(filename, edge1);


		Mat graphcut;
		Mat halfA, halfB;
		Mat halfA_mask, halfB_mask;

		GraphType g(est_nodes, est_edges);

		for (int i = 0; i < est_nodes; i++) {
			g.add_node();
		}

		int maxy = 0, miny = height, maxx = 0, minx = width;

		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				int idx0 = y*width + x;
				int idx1 = y*width + x + height * width;

				int a0m = (int)A0sMask.at<uchar>(y, x);
				int b0m = (int)B0sMask.at<uchar>(y, x);
				int a1m = (int)A1sMask.at<uchar>(y, x);
				int b1m = (int)B1sMask.at<uchar>(y, x);

				if (a0m < 128 && b0m >= 128)
					g.add_tweights(idx0, 0, INT_MAX);
				if (a0m >= 128 && b0m < 128)
					g.add_tweights(idx0, INT_MAX, 0);

				if (a1m < 128 && b1m >= 128)
					g.add_tweights(idx1, 0, INT_MAX);
				if (a1m >= 128 && b1m < 128)
					g.add_tweights(idx1, INT_MAX, 0);

				if (frameIndex > start)
				{
					int inLastFrame = (int)imgAll.at<uchar>(y, x);

					if (inLastFrame > 128)
					{
						// use preset edge
						if (mask.at<uchar>(y, x) == 0) {
							g.add_tweights(idx0, INT_MAX, 0);    //INT_MAX
						}
						else {
							g.add_tweights(idx0, 0, INT_MAX);
						}
					}
					
				}
			}
		}


		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				int idx0 = y*width + x;
				int idx1 = y*width + x + height * width;

				Vec3b a0 = A0s.at<Vec3b>(y, x);
				Vec3b b0 = B0s.at<Vec3b>(y, x);
				Vec3b a1 = A1s.at<Vec3b>(y, x);
				Vec3b b1 = B1s.at<Vec3b>(y, x);

				float weight = 1;
				//weight = pow(1 - r1, 4) * 5 + 1;
				//float weight = 1;

				double cap00 = weight * (norm(a0, b0) + edge0.at<uchar>(y, x));
				double cap10 = weight * (norm(a1, b1) + edge1.at<uchar>(y, x));

				// Add right edge
				if (x + 1 < width) {
					Vec3b a01 = A0s.at<Vec3b>(y, x + 1);
					Vec3b b01 = B0s.at<Vec3b>(y, x + 1);
					double cap01 = weight * (norm(a01, b01) + edge0.at<uchar>(y, x + 1));
					g.add_edge(idx0, idx0 + 1, (int)(cap00 + cap01), (int)(cap00 + cap01));

					Vec3b a11 = A1s.at<Vec3b>(y, x + 1);
					Vec3b b11 = B1s.at<Vec3b>(y, x + 1);
					double cap11 = weight * norm(a11, b11) + edge1.at<uchar>(y, x + 1);
					g.add_edge(idx1, idx1 + 1, (int)(cap10 + cap11), (int)(cap10 + cap11));
				}

				// Add bottom edge
				if (y + 1 < height) {
					Vec3b a02 = A0s.at<Vec3b>(y + 1, x);
					Vec3b b02 = B0s.at<Vec3b>(y + 1, x);
					double cap02 = weight * norm(a02, b02) + edge0.at<uchar>(y + 1, x);
					g.add_edge(idx0, idx0 + width, (int)(cap00 + cap02), (int)(cap00 + cap02));

					Vec3b a12 = A1s.at<Vec3b>(y + 1, x);
					Vec3b b12 = B1s.at<Vec3b>(y + 1, x);
					double cap12 = weight * norm(a12, b12) + edge1.at<uchar>(y + 1, x);
					g.add_edge(idx1, idx1 + width, (int)(cap10 + cap12), (int)(cap10 + cap12));
				}

				// Add back edge
				//norm(a0, b1) + norm(a1, b0) + 
				double cap3 = weight * (norm(a0, a1) + A0e.at<uchar>(y, x) + norm(b0, b1) + B0e.at<uchar>(y, x) + norm(a0, b1) + A1e.at<uchar>(y, x) + norm(a1, b0) + B1e.at<uchar>(y, x)) * 0.01;
				g.add_edge(idx0, idx1, (int)(weight * (norm(a0, b1) + A0e.at<uchar>(y, x) + B1e.at<uchar>(y, x)) * 0.02 * stable), (int)(weight * (norm(a1, b0) + A1e.at<uchar>(y, x) + B0e.at<uchar>(y, x)) * 0.02 * stable));

				//g.add_edge(idx0, idx1, (int)cap3 * 1000, (int)cap3 * 1000);
			}
		}
		int flow = g.maxflow();
		std::cout << frameIndex << "\tmax flow: " << flow << endl;


		

		graphcut = A0.clone();
		halfA = A0.clone();
		halfB = B0.clone();
		halfA_mask = Mat::zeros(A0.rows, A0.cols, CV_8U);
		halfB_mask = halfA_mask.clone();
		
		namedWindow("Result", WINDOW_AUTOSIZE);

		if (frameIndex == start)
		{
			for (int y = 0; y < A0.rows; y++) {
				for (int x = 0; x < A0.cols; x++) {
					Vec3b a0 = A0.at<Vec3b>(y, x);
					Vec3b b0 = B0.at<Vec3b>(y, x);
					Vec3b a1 = A1.at<Vec3b>(y, x);
					Vec3b b1 = B1.at<Vec3b>(y, x);
					int xx = (int)x * scale;
					int yy = (int)y * scale;
					int idx = yy*width + xx;
					if (isBlack(a0) && !isBlack(b0)){
						graphcut.at<Vec3b>(y, x) = B0.at<Vec3b>(y, x);
						halfA.at<Vec3b>(y, x) = Vec3b(0, 0, 0);
						halfB_mask.at<uchar>(y, x) = 255;
					}
					else
					if (!isBlack(a0) && isBlack(b0)){
						graphcut.at<Vec3b>(y, x) = A0.at<Vec3b>(y, x);
						halfB.at<Vec3b>(y, x) = Vec3b(0, 0, 0);
						halfA_mask.at<uchar>(y, x) = 255;
					}
					else
					if (g.what_segment(idx) == GraphType::SOURCE) {
						//graphcut.at<Vec3b>(y, x) = A0.at<Vec3b>(y, x);
						graphcut.at<Vec3b>(y, x) = Vec3b(0, 0, 0);
						halfA_mask.at<uchar>(y, x) = 255;
						halfB.at<Vec3b>(y, x) = Vec3b(0, 0, 0);
					}
					else {
						graphcut.at<Vec3b>(y, x) = Vec3b(0, 0, 0);
						halfB_mask.at<uchar>(y, x) = 255;
						halfA.at<Vec3b>(y, x) = Vec3b(0, 0, 0);
					}
				}
			}
			//sprintf(filename, "%s/C%d.png", trgFolder, frameIndex);
			//imwrite(filename, graphcut);
			sprintf(filename, "%s/C%dA.png", trgFolder, frameIndex);
			imwrite(filename, halfA);
			sprintf(filename, "%s/C%dB.png", trgFolder, frameIndex);
			imwrite(filename, halfB);
			/*sprintf(filename, "../C%dA_m.png", frameIndex);
			imwrite(filename, halfA_mask);
			sprintf(filename, "../C%dB_m.png", frameIndex);
			imwrite(filename, halfB_mask);*/
			Ptr<detail::Blender> blender;
			blender = detail::Blender::createDefault(detail::Blender::MULTI_BAND, true);
			detail::MultiBandBlender* mb = dynamic_cast<detail::MultiBandBlender*>(blender.get());
			//mb->setNumBands(7);
			blender->prepare(Rect(0, 0, A0.cols, A0.rows));
			blender->feed(A0, halfA_mask, Point(0, 0));
			blender->feed(B0, halfB_mask, Point(0, 0));
			Mat res, res_mask;
			blender->blend(res, res_mask);
			res.convertTo(res, CV_8UC3);
			std::sprintf(filename, "%s/D%d.png", trgFolder, frameIndex);
			imwrite(filename, res);
			imshow("result", res); waitKey(1);
		}

		mask = Mat::zeros(height, width, CV_8UC1);
		saveMask = Mat::zeros(height, width, CV_8UC1);

		halfA = A1.clone();
		halfB = B1.clone();
		halfA_mask = Mat::zeros(A0.rows, A0.cols, CV_8U);
		halfB_mask = halfA_mask.clone();


		for (int y = 0; y < A1.rows; y++) {
			for (int x = 0; x < A1.cols; x++) {
				int xx = (int)x * scale;
				int yy = (int)y * scale;
				int idx = yy*width + xx + height * width;
				if (g.what_segment(idx) == GraphType::SOURCE) {
					graphcut.at<Vec3b>(y, x) = A1.at<Vec3b>(y, x);
					mask.at<uchar>(yy, xx) = 0;
					saveMask.at<uchar>(yy, xx) = 0;
					halfA_mask.at<uchar>(y, x) = 255;
					halfB.at<Vec3b>(y, x) = Vec3b(0, 0, 0);

				}
				else {
					graphcut.at<Vec3b>(y, x) = B1.at<Vec3b>(y, x);
					mask.at<uchar>(yy, xx) = 1;
					saveMask.at<uchar>(yy, xx) = 255;
					halfB_mask.at<uchar>(y, x) = 255;
					halfA.at<Vec3b>(y, x) = Vec3b(0, 0, 0);
				}
			}
		}
		//sprintf(filename, "%s/mask%d.png", trgFolder, frameIndex);
		//imwrite(filename, saveMask);


		//sprintf(filename, "%s/C%d.png", trgFolder, frameIndex + 1);
		//imwrite(filename, graphcut);
		sprintf(filename, "%s/C%dA.png", trgFolder, frameIndex + 1);
		imwrite(filename, halfA);
		sprintf(filename, "%s/C%dB.png", trgFolder, frameIndex + 1);
		imwrite(filename, halfB);
		//sprintf(filename, "../C%dA_m.png", frameIndex + 1);
		//imwrite(filename, halfA_mask);
		//sprintf(filename, "../C%dB_m.png", frameIndex + 1);
		//imwrite(filename, halfB_mask);
		Ptr<detail::Blender> blender;
		blender = detail::Blender::createDefault(detail::Blender::MULTI_BAND, true);
		detail::MultiBandBlender* mb = dynamic_cast<detail::MultiBandBlender*>(blender.get());
		//mb->setNumBands(7);
		blender->prepare(Rect(0, 0, A0.cols, A0.rows));
		blender->feed(A0, halfA_mask, Point(0, 0));
		blender->feed(B0, halfB_mask, Point(0, 0));
		Mat res, res_mask;
		blender->blend(res, res_mask);
		res.convertTo(res, CV_8UC3);
		sprintf(filename, "%s/D%d.png", trgFolder, frameIndex + 1);
		imwrite(filename, res);
		cv::imshow("result", res); waitKey(1);
	}
	return;
}

int main(int argc, char **argv)
{
	if (argc == 6){
		// full-auto mode
		int stable = atoi(argv[1]);    // 1 - 5
		int start = atoi(argv[2]);
		int end = atoi(argv[3]);
		int ab = atoi(argv[4]);        // {0, 1}
		float scale = atof(argv[5]);   // 0.2
		findSeamAndBlend(stable, start, end, ab, scale);
	}
	else
	{
		cout << "Usage: ./SeamCut [1-5] start end [01] scale" << endl;
	}
	return 0;
	
}
