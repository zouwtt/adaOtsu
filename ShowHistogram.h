#ifndef SHOWHISTOGRAM_H //就是头文件名（全大写后加个_H

#define  SHOWHISTOGRAM_H
void SmoothHistogram(cv::Mat& hist) {
	int smoothRange = 3; // 平滑窗口的范围，你可以根据需要调整

	for (int i = smoothRange; i < hist.rows - smoothRange; ++i) {
		float sum = 0.0f;
		for (int j = -smoothRange; j <= smoothRange; ++j) {
			sum += hist.at<float>(i + j);
		}
		hist.at<float>(i) = sum / (2 * smoothRange + 1);
	}
}
cv::Mat ShowHistogram(cv::Mat& image, int flag) {
	// 创建直方图存储
	cv::Mat hist;
	int histSize = 256;  // 直方图的 bin 数量
	float range[] = { 0, 256 };  // 像素值范围
	const float* ranges[] = { range };
	cv::calcHist(&image, 1, 0, cv::Mat(), hist, 1, &histSize, ranges, true, false);
	if (flag == 1) {
		SmoothHistogram(hist);
	};
	// 获取直方图的最大值以便进行归一化
	double max_value;
	cv::minMaxLoc(hist, 0, &max_value);

	// 获取窗口的大小
	cv::Size windowSize(800, 600);  // 自定义窗口大小，您可以根据实际需要修改

	// 计算直方图的宽度
	int histWidth = windowSize.width;

	// 创建直方图绘图对象
	cv::Mat histImage = cv::Mat::zeros(windowSize, CV_8UC3);
	histImage.setTo(cv::Scalar(240, 240, 240)); // 设置直方图的背景颜色为灰色

	// 使用渐变颜色绘制直方图
	for (int i = 0; i < histSize; i++) {
		float binValue = hist.at<float>(i);
		int intensity = cvRound(binValue * windowSize.height / max_value); // 根据最大值进行归一化

		cv::Scalar color(255, 140, 30); // 渐变色

		int barHeight = std::min(intensity, windowSize.height);
		cv::rectangle(histImage, cv::Point(i * histWidth / histSize, windowSize.height), cv::Point((i + 1) * histWidth / histSize, windowSize.height - barHeight), color, -1);
		cv::rectangle(histImage, cv::Point(i * histWidth / histSize, windowSize.height), cv::Point((i + 1) * histWidth / histSize, windowSize.height - barHeight), cv::Scalar(0, 0, 0), 1);
	}

	// 添加刻度标签到X轴
	cv::Point textPos;
	cv::Scalar textColor(0, 0, 0); // 字体颜色为黑色
	cv::Point axisStart(0, windowSize.height);
	cv::Point axisEnd(windowSize.width, windowSize.height);
	int textOffset = 20;

	for (int i = 0; i <= histSize; i += 20) {
		std::string text = std::to_string(i * 256 / histSize);
		textPos.x = i * histWidth / histSize - textOffset;
		textPos.y = windowSize.height + 20;
		cv::putText(histImage, text, textPos, cv::FONT_HERSHEY_SIMPLEX, 0.5, textColor);
	}
	// 绘制坐标轴
	//cv::rectangle(histImage, axisStart, axisEnd, cv::Scalar(255, 255, 255), 3);
	//cv::namedWindow("Histogram with Axis and Style", cv::WINDOW_NORMAL);
	//cv::imshow("Histogram with Axis and Style", histImage);
	//cv::waitKey(0);
	//cv::destroyWindow("Histogram with Axis and Style");
	return hist;
}
#endif

