#ifndef OTSUSMAWKDP_H //就是头文件名（全大写后加个_H

#define OTSUSMAWKDP_H

void mfill(unsigned int m, unsigned int rowM, unsigned int rowO, struct ROWELEMENT* reducedMatrix, int* candidate);
struct ROWELEMENT* reduce(unsigned int m, unsigned int p, unsigned int rowM, unsigned int rowO, struct ROWELEMENT* myMatrix, int* candidate);
void msearch(unsigned int m, unsigned int p, unsigned int rowM, unsigned int rowO,
	struct ROWELEMENT* myMatrix, struct ROWELEMENT* lastMatrix, int* candidate);
double otsuSmawkDP(int M, int L, DTYPE* p, int* candidate, int* thresholds)
{
	unsigned int i, j;
	unsigned int Q = L - M + 1;          /* nodes per class */
	struct NODE* searchNODE;             /* used for the backtracking */
	DTYPE Omax = 0;
	int* index;
	unsigned int   Omaxpos = 0;
	unsigned int thisSize;
	unsigned int colmemSize;
	struct ROWELEMENT* colMem;           /* memory for the smawk algorithm */

#ifdef TIMING
	struct timeval tStart;      /* start time */
	struct timeval tEnd;        /* end time   */
	struct timeval tDelta;      /* time difference */
#endif

#ifdef NCALC
	nc = 0;
#endif


	/* calculate the size of the memory needed by the smawk algorithm */
	thisSize = Q;
	colmemSize = Q + 1;
	while (thisSize >= 2) {
		colmemSize += thisSize + 1;
		thisSize = thisSize / 2;
	}
	/* allocale the memory for the smawk algorithm */
	if ((colMem = (ROWELEMENT*)malloc(colmemSize * sizeof(struct ROWELEMENT))) == NULL) {
		printf("malloc for colMem failed\n");
	}

	if ((maxcols = (unsigned int*)malloc(Q * sizeof(unsigned int))) == NULL) {
		printf("malloc for maxcols failed\n");
	}
	if ((index = (int*)malloc(M * sizeof(int))) == NULL) {
		printf("malloc for index array failed\n");
	}

	/* allocate memory for the trellis structure */
	if ((trellis = (NODE**)malloc(M * sizeof(struct NODE*))) == NULL) {
		printf("malloc for trellis M-dim falied\n");
	}
	for (i = 0; i < M; i++) {
		if ((trellis[i] = (NODE*)malloc(L * sizeof(struct NODE))) == NULL) {
			printf("malloc for trellis L-dim failed\n");
		}
		memset(trellis[i], 0, L * sizeof(struct NODE));
	}

	/* allocate memory for the N and D array */
	if ((N = (DTYPE*)malloc(256 * sizeof(DTYPE))) == NULL) {
		printf("malloc for N array failed\n");
	}
	if ((D = (DTYPE*)malloc(256 * sizeof(DTYPE))) == NULL) {
		printf("malloc for D array failed\n");
	}

#ifdef TIMING
	/* start timing here */
	gettimeofday(&tStart, NULL);
#endif

	/* initialize N and D array*/
	N[0] = 0;
	D[0] = p[0];
	for (i = 1; i < 256; i++) {
		N[i] = N[i - 1] + i * p[i];
		D[i] = D[i - 1] + p[i];
	}
	/* initialize the trellis structure for the lowest level (class = 0) */
	if (D[candidate[0]] == 0)            /* to avoid division by zero */
		trellis[0][0].O = 0;
	else
		trellis[0][0].O = N[candidate[0]] * N[candidate[0]] / D[candidate[0]];
	for (i = 1; i < Q; i++) {
		if (D[candidate[i]] == 0)
			trellis[0][i].O = trellis[0][i - 1].O;
		else
			trellis[0][i].O = N[candidate[i]] * N[candidate[i]] / D[candidate[i]];
	}

	/* begin shortest (here longest) path algorithm */
	for (c = 1; c < (M - 1); c++) {
		msearch(Q, Q, 1, 0, colMem, NULL, candidate);
	}

	c = M - 1; i = L - 1; Omax = 0;
	for (j = M - 2; j < (L - 1); j++) {
		//Ntemp = N[i] - N[j];
		//Dtemp = D[i] - D[j];
		Ntemp = N[candidate[i]] - N[candidate[j]];
		Dtemp = D[candidate[i]] - D[candidate[j]];
		if (Dtemp == 0)
			Otemp = trellis[c - 1][j].O;
		else
			Otemp = trellis[c - 1][j].O + Ntemp * Ntemp / Dtemp;
		if (Otemp > Omax) {
			Omax = Otemp;
			Omaxpos = j;
		}
	}
	trellis[c][i].O = Omax;
	trellis[c][i].p_back = &trellis[c - 1][Omaxpos];
	        

	searchNODE = &(trellis[M - 1][L - 1]);
	for (i = M - 1; i > 0; i--) {
		searchNODE = searchNODE->p_back;
		index[i - 1] = (searchNODE - trellis[i - 1]);
	}
	for (i = M - 1; i > 0; i--) {

		thresholds[i - 1] = candidate[index[i - 1]];
	}
#ifdef TIMING
	/* end timing, calculate time difference */
	gettimeofday(&tEnd, NULL);
	timevalSubtract(&tDelta, &tEnd, &tStart);
#ifdef NCALC
	printf("%f\t", (double)(tDelta.tv_sec) + 1.0E-6 * (double)(tDelta.tv_usec));
#else
	printf("%f\n", (double)(tDelta.tv_sec) + 1.0E-6 * (double)(tDelta.tv_usec));
#endif
#endif

#ifdef NCALC
	nc += 2 * (L - M + 1);
	printf("%u\n", nc);
#endif
	DTYPE s = 0.0f;
	for (int j = 0; j < M - 1; j++) {
		s += p[thresholds[j]];
	}
	Omax = (1.0 - s) * Omax;

	free(colMem);
	free(maxcols);
	for (i = 0; i < M; i++) {
		free(trellis[i]);
	}
	free(trellis);
	free(N);
	free(D);

	return Omax;
}
void mfill(unsigned int m, unsigned int rowM, unsigned int rowO, struct ROWELEMENT* reducedMatrix, int* candidate)
{
	struct ROWELEMENT* pCol;
	unsigned int evenMaxPos;
	unsigned int maxCol;
	unsigned int row;
	DTYPE max;
	DTYPE temp;

#ifdef NCALC
	nc += m; /* matrix is square, one calculation per column */
#endif

	if ((m % 2) == 0) {
		/* even number of rows, the last row of the matrix is even, we know the position of this maxima */
		row = rowO + (m - 2) * rowM;
		evenMaxPos = maxcols[row + rowM];
		pCol = reducedMatrix;
		while (pCol->col > evenMaxPos) {
			pCol = pCol->p_lastnext;
		}
	}
	else {
		/* odd number off rows, the last row of the matrix is odd */
		row = rowO + (m - 1) * rowM;
		pCol = reducedMatrix;
	}

	/* pCol points now to the column where the first maximum in an odd row could be */
	/* do this for all odd rows, except for the first one                           */
	while (row > rowO) {
		evenMaxPos = maxcols[row - rowM];
		max = 0;
		while (pCol->col > evenMaxPos) {
			if (pCol->col > row) {
				temp = 0;
			}
			else {
				Dtemp = D[candidate[c + row]] - D[candidate[c - 1 + pCol->col]];
				if (Dtemp == 0) {
					temp = trellis[c - 1][c - 1 + pCol->col].O;
				}
				else {
					Ntemp = N[candidate[c + row]] - N[candidate[c - 1 + pCol->col]];
					temp = trellis[c - 1][c - 1 + pCol->col].O + Ntemp * Ntemp / Dtemp;
				}
			}
			if (temp >= max) {
				max = temp;
				maxCol = pCol->col;
			}
			pCol = pCol->p_lastnext;
		}
		if (pCol->col > row) {
			temp = 0;
		}
		else {
			Dtemp = D[candidate[c + row]] - D[candidate[c - 1 + pCol->col]];
			if (Dtemp == 0) {
				temp = trellis[c - 1][c - 1 + pCol->col].O;
			}
			else {
				Ntemp = N[candidate[c + row]] - N[candidate[c - 1 + pCol->col]];
				temp = trellis[c - 1][c - 1 + pCol->col].O + Ntemp * Ntemp / Dtemp;
			}
		}
		if (temp >= max) {
			maxcols[row] = pCol->col;
			trellis[c][c + row].O = temp;
			trellis[c][c + row].p_back = &trellis[c - 1][c - 1 + pCol->col];
		}
		else {
			maxcols[row] = maxCol;
			trellis[c][c + row].O = max;
			trellis[c][c + row].p_back = &trellis[c - 1][c - 1 + maxCol];
		}
		row = row - 2 * rowM;
	}
	/* for the first row */
	max = 0;
	while (pCol->p_lastnext != NULL) {
		if (pCol->col > row) {
			temp = 0;
		}
		else {
			Dtemp = D[candidate[c + row]] - D[candidate[c - 1 + pCol->col]];
			if (Dtemp == 0) {
				temp = trellis[c - 1][c - 1 + pCol->col].O;
			}
			else {
				Ntemp = N[candidate[c + row]] - N[candidate[c - 1 + pCol->col]];
				temp = trellis[c - 1][c - 1 + pCol->col].O + Ntemp * Ntemp / Dtemp;
			}
		}
		if (temp >= max) {
			max = temp;
			maxCol = pCol->col;
		}
		pCol = pCol->p_lastnext;
	}
	if (pCol->col > row) {
		temp = 0;
	}
	else {
		Dtemp = D[candidate[c + row]] - D[candidate[c - 1 + pCol->col]];
		if (Dtemp == 0) {
			temp = trellis[c - 1][c - 1 + pCol->col].O;
		}
		else {
			Ntemp = N[candidate[c + row]] - N[candidate[c - 1 + pCol->col]];
			temp = trellis[c - 1][c - 1 + pCol->col].O + Ntemp * Ntemp / Dtemp;
		}
	}
	if (temp >= max) {
		maxcols[row] = pCol->col;
		trellis[c][c + row].O = temp;
		trellis[c][c + row].p_back = &trellis[c - 1][c - 1 + pCol->col];
	}
	else {
		maxcols[row] = maxCol;
		trellis[c][c + row].O = max;
		trellis[c][c + row].p_back = &trellis[c - 1][c - 1 + maxCol];
	}
	return;
}
struct ROWELEMENT* reduce(unsigned int m, unsigned int p, unsigned int rowM, unsigned int rowO, struct ROWELEMENT* myMatrix, int* candidate)
{
	unsigned int j = 1;
	unsigned int k = 0;
	unsigned int x;
	unsigned int y;
	unsigned int ncols = p;
	DTYPE val1;
	DTYPE val2;

	/* go trough matrix columns from left to right, delete columns until matrix is square */
	while (ncols > m) {
		x = rowO + k * rowM;
		y = (myMatrix[j].p_lastnext)->col;
		if (y > x) {
			val1 = 0;
		}
		else {
#ifdef NCALC
			nc++;
#endif 
			Dtemp = D[candidate[c + x]] - D[candidate[c - 1 + y]];
			if (Dtemp == 0) {
				val1 = trellis[c - 1][c - 1 + y].O;
			}
			else {
				Ntemp = N[candidate[c + x]] - N[candidate[c - 1 + y]];
				val1 = trellis[c - 1][c - 1 + y].O + Ntemp * Ntemp / Dtemp;
			}
		}
		y = myMatrix[j].col;
		if (y > x) {
			val2 = 0;
		}
		else {
#ifdef NCALC
			nc++;
#endif
			Dtemp = D[candidate[c + x]] - D[candidate[c - 1 + y]];
			if (Dtemp == 0) {
				val2 = trellis[c - 1][c - 1 + y].O;
			}
			else {
				Ntemp = N[candidate[c + x]] - N[candidate[c - 1 + y]];
				val2 = trellis[c - 1][c - 1 + y].O + Ntemp * Ntemp / Dtemp;
			}
		}

		if (val1 >= val2) {
			if (k < (m - 1)) {
				k++;
			}
			else {
				/* delete column j */
				myMatrix[j + 1].p_lastnext = myMatrix[j].p_lastnext;
				ncols--;
			}
			j++;
		}
		else {
			/* delete column j - 1 */
			myMatrix[j].p_lastnext = (myMatrix[j].p_lastnext)->p_lastnext;
			ncols--;
			if (k != 0)
				k--;
			else
				j++;
		}
	}
	/* the rightmost element (dummy) in the column array points to the rightmost              */
	/* column element of our matrix, return this pointer.                                     */
	return myMatrix[p].p_lastnext;
}
void msearch(unsigned int m, unsigned int p, unsigned int rowM, unsigned int rowO,
	struct ROWELEMENT* myMatrix, struct ROWELEMENT* lastMatrix, int* candidate)
{
	struct ROWELEMENT* temp;
	struct ROWELEMENT* reducedMatrix;
	unsigned int i;

	/* if lastMatrix is NULL, we are not in a recursion, the col-indices are initialized with global values */
	/* there is one column element more than the matrix has rows, this is needed for reduce                 */
	if (lastMatrix == NULL) {
		myMatrix[0].col = 0;
		myMatrix[0].p_lastnext = NULL;
		for (i = 1; i < p + 1; i++) {
			myMatrix[i].col = i;
			myMatrix[i].p_lastnext = &myMatrix[i - 1];
		}
	}
	/* otherwise we have to initialize our matrix with the values of the last matrix */
	/* lastMatrix points to the rightmost column element of the previous matrix      */
	/* this means we fill our matrix from right to left                              */
	else {
		/* rightmost element, dummy for reduce */
		myMatrix[p].col = 0;
		myMatrix[p].p_lastnext = &myMatrix[p - 1];
		temp = lastMatrix;
		i = p - 1;
		do {
			/* traversing backward, matrix has always more than one column */
			myMatrix[i].col = temp->col;
			myMatrix[i].p_lastnext = &myMatrix[i - 1];
			temp = temp->p_lastnext;
			i--;
		} while (temp->p_lastnext != NULL);
		myMatrix[0].col = temp->col;
		myMatrix[0].p_lastnext = NULL;
	}
	/* reduce the matrix to a m*m matrix */
	reducedMatrix = reduce(m, p, rowM, rowO, myMatrix, candidate);
	/* check if t\ he   reduced matrix has size 1x1, if yes we are at the lowest point of the recursion */
	/* fill in the column number & the value in the result array, then return                       */
	if (reducedMatrix->p_lastnext == NULL) {
		maxcols[rowO] = reducedMatrix->col;
		trellis[c][c + rowO].p_back = &trellis[c - 1][c - 1 + reducedMatrix->col];
#ifdef NCALC
		nc++;
#endif
		Dtemp = D[candidate[c + rowO]] - D[candidate[c - 1 + reducedMatrix->col]];
		if (Dtemp == 0)
			trellis[c][c + rowO].O = trellis[c - 1][c - 1 + reducedMatrix->col].O;
		else {
			Ntemp = N[candidate[c + rowO]] - N[candidate[c - 1 + reducedMatrix->col]];
			trellis[c][c + rowO].O = trellis[c - 1][c - 1 + reducedMatrix->col].O + Ntemp * Ntemp / Dtemp;
		}
		/* this recursion is finished, return */
		return;
	}
	/* call msearch recursively, reduced matrix is m * m, we pass only even rows */
	msearch(m / 2, m, 2 * rowM, rowM + rowO, myMatrix + p + 1, reducedMatrix, candidate);
	/* call mfill to fill in the values of the odd rows */
	mfill(m, rowM, rowO, reducedMatrix, candidate);
	return;
}

#endif
