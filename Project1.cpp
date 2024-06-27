#ifdef _DEBUG
#pragma comment(lib, "opencv_world300d.lib")
#else
#pragma comment(lib, "opencv_world300.lib")
#endif // DEBUG
#include "def.h"
#include "ExtractContour.h"
#include "compareDifferences.h"
#include "otsuExhaust.h"
#include "otsuOrigDP.h"
#include "otsuDivDP.h"
#include "ExtractConnectedRegion.h"
#include "otsuSmawkDP.h"
#include "GenerateAndShowHistogram.h"
#include "ShowHistogram.h"
#include "PerformImageSegmentation.h"
#include "region.h"
#include "ConvertToLineArt.h"
#include "pca.h"
#include "UM.h"
#include "saveBestThresholds.h"
void main()
{
	string inputDirectory = "E:/thresholding/Project1/image/BSDS300/test/";  
	string outputDirectory = "E:/thresholding/Project1/output/200-80/"; 
	std::string txtFilePath = "E:/thresholding/Project1/image/BSDS300/iids_test.txt"; 
	std::ifstream inputFile(txtFilePath);
	if (!inputFile.is_open()) {
		std::cerr << "Error: Could not open the input file" << std::endl;
		return;
	}
	std::vector<std::string> imageNames;
	std::string imageName;
	while (std::getline(inputFile, imageName)) {
		imageNames.push_back(imageName);
	}
	inputFile.close();
	int L = 256;
	double p[256];
	int i;
	int newL = 80;
	int M;
	int bestM;
	DTYPE Omax;
	int* candidate = (int*)malloc(newL * sizeof(int));
	for (int i = 0; i < L; i++) p[i] = 0;
	auto start_time = std::chrono::high_resolution_clock::now();
	for (const std::string& imageName : imageNames) {
		auto start_time = std::chrono::high_resolution_clock::now();
		DTYPE bestOmax = 0.0;
		int* bestTH = NULL;
		//string inputDirectory = inputDirectory + imageName + ".jpg";
		string inputFileName = inputDirectory + imageName + ".jpg";
		cv::Mat image = cv::imread(inputFileName, cv::IMREAD_GRAYSCALE);
		if (image.empty()) {
			cerr << "Error: Could not read image " << imageName << endl;
			continue; 
		}
		cv::Mat thresholding = image.clone();
		//image = pcaimage(image);
		//ShowHistogram(image,0);
		cv::Mat Hist;
		Hist = ShowHistogram(image, 0);
		cv::Mat smoothedHist;
		smoothedHist = ShowHistogram(image, 1);
		cv::Mat img = image.clone();
		//cv::Mat thresholding = image.clone();
		cv::Size imgSize = image.size();
		GenerateAndShowHistogram(Hist, smoothedHist, candidate, newL, p);
		for (int M = 11; M <= 11; M++) {
			int* th = (int*)malloc((M - 1) * sizeof(int));/*存放阈值*/
			Omax = otsuSmawkDP(M, newL, p, candidate, th);
			if (Omax > bestOmax) {
				bestOmax = Omax;
				bestM = M;
				if (bestTH != NULL) {
					free(bestTH);
				}
				bestTH = th;
			}
			else {
				free(th);
			}
		}
		
		free(bestTH);
		auto end_time = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
		std::cout << "代码运行时间: " << duration.count() << " 毫秒" << std::endl;
	}
	auto end_time = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
	std::cout << "代码运行时间: " << duration.count() << " 毫秒" << std::endl;
	return;
}