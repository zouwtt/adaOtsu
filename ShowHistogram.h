#ifndef SHOWHISTOGRAM_H //����ͷ�ļ�����ȫ��д��Ӹ�_H

#define  SHOWHISTOGRAM_H
void SmoothHistogram(cv::Mat& hist) {
	int smoothRange = 3; // ƽ�����ڵķ�Χ������Ը�����Ҫ����

	for (int i = smoothRange; i < hist.rows - smoothRange; ++i) {
		float sum = 0.0f;
		for (int j = -smoothRange; j <= smoothRange; ++j) {
			sum += hist.at<float>(i + j);
		}
		hist.at<float>(i) = sum / (2 * smoothRange + 1);
	}
}
cv::Mat ShowHistogram(cv::Mat& image, int flag) {
	// ����ֱ��ͼ�洢
	cv::Mat hist;
	int histSize = 256;  // ֱ��ͼ�� bin ����
	float range[] = { 0, 256 };  // ����ֵ��Χ
	const float* ranges[] = { range };
	cv::calcHist(&image, 1, 0, cv::Mat(), hist, 1, &histSize, ranges, true, false);
	if (flag == 1) {
		SmoothHistogram(hist);
	};
	// ��ȡֱ��ͼ�����ֵ�Ա���й�һ��
	double max_value;
	cv::minMaxLoc(hist, 0, &max_value);

	// ��ȡ���ڵĴ�С
	cv::Size windowSize(800, 600);  // �Զ��崰�ڴ�С�������Ը���ʵ����Ҫ�޸�

	// ����ֱ��ͼ�Ŀ��
	int histWidth = windowSize.width;

	// ����ֱ��ͼ��ͼ����
	cv::Mat histImage = cv::Mat::zeros(windowSize, CV_8UC3);
	histImage.setTo(cv::Scalar(240, 240, 240)); // ����ֱ��ͼ�ı�����ɫΪ��ɫ

	// ʹ�ý�����ɫ����ֱ��ͼ
	for (int i = 0; i < histSize; i++) {
		float binValue = hist.at<float>(i);
		int intensity = cvRound(binValue * windowSize.height / max_value); // �������ֵ���й�һ��

		cv::Scalar color(255, 140, 30); // ����ɫ

		int barHeight = std::min(intensity, windowSize.height);
		cv::rectangle(histImage, cv::Point(i * histWidth / histSize, windowSize.height), cv::Point((i + 1) * histWidth / histSize, windowSize.height - barHeight), color, -1);
		cv::rectangle(histImage, cv::Point(i * histWidth / histSize, windowSize.height), cv::Point((i + 1) * histWidth / histSize, windowSize.height - barHeight), cv::Scalar(0, 0, 0), 1);
	}

	// ��ӿ̶ȱ�ǩ��X��
	cv::Point textPos;
	cv::Scalar textColor(0, 0, 0); // ������ɫΪ��ɫ
	cv::Point axisStart(0, windowSize.height);
	cv::Point axisEnd(windowSize.width, windowSize.height);
	int textOffset = 20;

	for (int i = 0; i <= histSize; i += 20) {
		std::string text = std::to_string(i * 256 / histSize);
		textPos.x = i * histWidth / histSize - textOffset;
		textPos.y = windowSize.height + 20;
		cv::putText(histImage, text, textPos, cv::FONT_HERSHEY_SIMPLEX, 0.5, textColor);
	}
	// ����������
	//cv::rectangle(histImage, axisStart, axisEnd, cv::Scalar(255, 255, 255), 3);
	//cv::namedWindow("Histogram with Axis and Style", cv::WINDOW_NORMAL);
	//cv::imshow("Histogram with Axis and Style", histImage);
	//cv::waitKey(0);
	//cv::destroyWindow("Histogram with Axis and Style");
	return hist;
}
#endif

