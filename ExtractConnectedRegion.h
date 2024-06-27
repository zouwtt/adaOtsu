#ifndef EXTRACTCONNECTEDREGION_H //就是头文件名（全大写后加个_H

#define EXTRACTCONNECTEDREGION_H

int CalcDiffSquare(int pixel1, int pixel2);
int ExtractConnectedRegion(vector<int>& vecLabel, const vector<vector<int>>& vecImgData, const CvSize& ImgSize, const int nMinRegion)
{
	int nRegionNum(0);
	if (vecLabel.empty() || vecImgData.empty() || (vecLabel.size() != ImgSize.width * ImgSize.height) || nMinRegion < 1)
		return nRegionNum;

	CvMat* pMat = cvCreateMat(ImgSize.height, ImgSize.width, CV_32FC1);
	int nStep = pMat->step / sizeof(pMat->data.fl[0]);
	float* pData = pMat->data.fl;

	for (int i = 0, y = 0; y < ImgSize.height; y++)
	{
		for (int x = 0; x < ImgSize.width; x++, i++)
		{
			pData[x] = -vecLabel[i];
		}
		pData += nStep;
}

	// arching the connected regions
	pData = pMat->data.fl;
	for (int y = 0; y < ImgSize.height; y++)
	{
		for (int x = 0; x < ImgSize.width; x++)
		{
			if (pData[x] <= 0)
			{
				nRegionNum++;
				cvFloodFill(pMat, cvPoint(x, y), cvScalar(nRegionNum));
			}
		}
		pData += nStep;
	}

	//calculating the size of each region
	vector<int> vecRegionNum(nRegionNum);
	pData = pMat->data.fl;
	for (int y = 0; y < ImgSize.height; y++)
	{
		for (int x = 0; x < ImgSize.width; x++)
		{
			int nIdx = pData[x] - 1;
			vecRegionNum[nIdx]++;
		}
		pData += nStep;
	}

	//calculating the edge
	vector<Edge> vecConEdge;
	pData = pMat->data.fl;
	int nChannels = vecImgData.size();
	for (int i = 0, y = 0; y < ImgSize.height; y++)
	{
		for (int x = 0; x < ImgSize.width; x++, i++)
		{
			int nCurIdx = (int)pData[x] - 1;//current pixel: (x,y)
			Edge edge;
			edge.pt1 = cvPoint(x, y);

			if (x < ImgSize.width - 1)//left neighbor pixel: (x + 1, y)
			{
				int nLeftIdx = (int)pData[x + 1] - 1;
				if ((nCurIdx != nLeftIdx) && (vecRegionNum[nCurIdx] < nMinRegion || vecRegionNum[nLeftIdx] < nMinRegion))
				{
					edge.pt2 = cvPoint(x + 1, y);
					edge.dDist = 0;
					for (int k = 0; k < nChannels; k++)
					{
						edge.dDist += CalcDiffSquare(vecImgData[k][i], vecImgData[k][i + 1]);
					}

					vecConEdge.push_back(edge);
				}
			}

			if (y < ImgSize.height - 1)//bottom neighbor pixel: (x, y + 1)
			{
				int nBottomIdx = (int)pData[x + ImgSize.width] - 1;
				if ((nCurIdx != nBottomIdx) && (vecRegionNum[nCurIdx] < nMinRegion || vecRegionNum[nBottomIdx] < nMinRegion))
				{
					edge.pt2 = cvPoint(x, y + 1);
					edge.dDist = 0;
					// TODO DEBUG
					for (int k = 0; k < nChannels; k++)
					{
						edge.dDist += CalcDiffSquare(vecImgData[k][i], vecImgData[k][ImgSize.width + i]);
						vecConEdge.push_back(edge);
					}
				}
			}
		}
		pData += nStep;
	}
	//sort(vecConEdge.begin(), vecConEdge.end());
	std::sort(vecConEdge.begin(), vecConEdge.end(), CompareEdges());
	//merging the small region into its neighbor
	// TODO DEBUG 584~610
	int nEdgeNum = vecConEdge.size();
	for (int i = 0; i < nEdgeNum; i++)
	{
		int nLabel1 = cvmGet(pMat, vecConEdge[i].pt1.y, vecConEdge[i].pt1.x) - 1;
		int nLabel2 = cvmGet(pMat, vecConEdge[i].pt2.y, vecConEdge[i].pt2.x) - 1;
		if (nLabel1 != nLabel2 && (vecRegionNum[nLabel1] < nMinRegion || vecRegionNum[nLabel2] < nMinRegion))
		{
			if (vecRegionNum[nLabel2] > vecRegionNum[nLabel1])
				cvFloodFill(pMat, vecConEdge[i].pt1, cvScalar(nLabel2 + 1));
			else
				cvFloodFill(pMat, vecConEdge[i].pt2, cvScalar(nLabel1 + 1));

			int nNum = vecRegionNum[nLabel1];
			vecRegionNum[nLabel1] += vecRegionNum[nLabel2];
			vecRegionNum[nLabel2] += nNum;
		}
	}

	pData = pMat->data.fl;
	for (int y = 0; y < ImgSize.height; y++)
	{
		for (int x = 0; x < ImgSize.width; x++)
		{
			pData[x] = -pData[x];
		}
		pData += nStep;
	}

	//finding the connected regions again
	nRegionNum = 0;
	pData = pMat->data.fl;
	for (int y = 0; y < ImgSize.height; y++)
	{
		for (int x = 0; x < ImgSize.width; x++)
		{
			if (pData[x] <= 0)
			{
				nRegionNum++;
				cvFloodFill(pMat, cvPoint(x, y), cvScalar(nRegionNum));
			}
		}
		pData += nStep;
	}

	pData = pMat->data.fl;
	for (int i = 0, y = 0; y < ImgSize.height; y++)
	{
		for (int x = 0; x < ImgSize.width; x++, i++)
		{
			vecLabel[i] = pData[x] - 1; // arting from zero
		}
		pData += nStep;
	}
	return nRegionNum;
}
int CalcDiffSquare(int pixel1, int pixel2) {
	int diff = pixel1 - pixel2;
	int diffSquare = diff * diff;
	return diffSquare;
}
#endif
