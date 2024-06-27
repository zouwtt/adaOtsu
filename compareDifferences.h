#ifndef COMPAREDIFFERENCES_H //����ͷ�ļ�����ȫ��д��Ӹ�_H

#define COMPAREDIFFERENCES_H
int compareDifferences(const void* a, const void* b) {
	struct PIX* p1 = (PIX*)a;
	struct PIX* p2 = (PIX*)b;
	if (p1->Diff < p2->Diff)
	{
		return 1;
	}
	else if (p1->Diff == p2->Diff)
	{
		return 0;
	}
	else if (p1->Diff > p2->Diff)
	{
		return -1;
	}
	return 0;
}
#endif
