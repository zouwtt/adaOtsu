#ifndef GENERATEANDSHOWHISTOGRAM_H //就是头文件名（全大写后加个_H

#define GENERATEANDSHOWHISTOGRAM_H

void saveThresholds(const string& outputPath, int* TH, int can);
int compareI(const void* a, const void* b);
void sort(DTYPE* op, PIX* pixelDifferences, int* cur, int newL, int win);
void GenerateAndShowHistogram(cv::Mat& Hist,cv::Mat& smoothedHist, int* candidate, int newL, double* p) {
	int L = 256;
	// 获取总像素数
    double op[256];
    for (int i = 0; i < L; i++) op[i] = 0;
	int total = 0;
	for (int i = 0; i < L; ++i) {
		total += smoothedHist.at<float>(i);
	}
	// 计算每个像素值的概率并填充数组 p
	for (int i = 0; i < L; ++i) {
		p[i] = smoothedHist.at<float>(i) / total;
        op[i] = smoothedHist.at<float>(i);
	}
	PIX* pixelDifferences;
	if ((pixelDifferences = (PIX*)malloc(L * sizeof(PIX))) == NULL) {
		printf("malloc for pixel differences array failed\n");
	}
	for (int i = 0; i < L; i++) {
		pixelDifferences[i].Diff = 0;
		pixelDifferences[i].p = 0;
	}
	sort(op, pixelDifferences, candidate, newL, 1);
	qsort(candidate, newL, sizeof(int), compareI);
	//for (int i = 0; i < newL; i++) {
	//	printf("%d\n ", candidate[i]);
	//}
}
int compareI(const void* a, const void* b) {
	int* p1 = (int*)a;
	int* p2 = (int*)b;
	if (*p1 > *p2)
	{
		return 1;
	}
	else if (*p1 == *p2)
	{
		return 0;
	}
	else if (*p1 < *p2)
	{
		return -1;
	}
	return 0;
}
void sort(DTYPE* op, PIX* pixelDifferences, int* cur, int newL, int win) {
    int i;
    int L = 256;
    int num = 0; // Counter for elements with positive differences

    // Calculate pixel differences (left and right neighbors)
    for (i = win; i < L; i++) {
        DTYPE leftNeighbor = (i > win - 1) ? op[i - win] - op[i] : 0;
        DTYPE rightNeighbor = (i < L - win) ? op[i + win] - op[i] : 0;
        pixelDifferences[i].Diff = leftNeighbor + rightNeighbor;
        pixelDifferences[i].p = i;

        if (pixelDifferences[i].Diff > 0) {
            num++; // Increment counter for positive differences
        }
    }
  //  printf("num=%d\n", num);
    // Create a temporary array to store only elements with positive differences
    PIX* tempArray = (PIX*)malloc(num * sizeof(PIX));
    if (tempArray == NULL) {
        // Handle memory allocation error
        return;
    }

    // Copy elements with positive differences to the temporary array
    int tempIndex = 0;
    for (i = win; i < L; i++) {
        if (pixelDifferences[i].Diff > 0) {
            tempArray[tempIndex++] = pixelDifferences[i];
        }
    }

    // Calculate mean and standard deviation of differences in the temporary array
    DTYPE sum = 0.0;
    for (i = 0; i < num; i++) {
        sum += tempArray[i].Diff;
    }
    DTYPE mean = sum / num;

    DTYPE sum_squared_diff = 0.0;
    for (i = 0; i < num; i++) {
        sum_squared_diff += pow(tempArray[i].Diff - mean, 2);
    }
    DTYPE std_dev = sqrt(sum_squared_diff / num);

    // Count elements satisfying the condition tempArray[i].Diff > (mean + 1.645 * std_dev)
    int count = 0;
    int k = 0.7;
    for (i = 0; i < num; i++) {
        while (mean - k * std_dev<0) {
            k= k - 0.1;

        }
        if (tempArray[i].Diff > (mean - k * std_dev)) {
            count++;
        }
    }

    // Sort the temporary array based on differences
    qsort(tempArray, num, sizeof(PIX), compareDifferences);
    // Copy sorted elements back to the original array based on the condition
    for (i = 0; i < newL; i++) {
            cur[i] = tempArray[i].p;
    }
    saveThresholds("E:/灰度图像/thresholds100.txt", cur, 30);
    // Free memory allocated for the temporary array
    free(tempArray);
}
void saveThresholds(const string& outputPath, int* TH, int can) {
    ofstream outputFile(outputPath, ios::app);
    if (!outputFile.is_open()) {
        cerr << "Error: Could not open the output file" << endl;
        return;
    }
    for (int i = 0; i < can - 1; i++) {
        outputFile << TH[i] << " ";
    }
    outputFile << "\n";

    outputFile.close();
}
#endif
