/*
 * =====================================================================================
 *
 *       Filename:  space.cpp
 *
 *    Description:  计算节点序数、坐标等网格参数
 *
 *        Version:  1.0
 *        Created:  2015年11月30日 23时39分42秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Duan Junming (DSEC), duanjm@pku.edu.cn
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef _SPACE_H__
#define _SPACE_H__

#include "vvector.h"
#include "MAT.h"

typedef MAT<double> EN;
typedef MAT<double> CD;

void triangle_en(int N, EN& en) {
	en.resize(3, N*N*2);
	for(int j = 0; j < N; ++j) {
		for(int i = 0; i < N; ++i) {
			en[0][j*N*2+i*2] = (N+1)*j + i;
			en[1][j*N*2+i*2] = (N+1)*j + i + N + 2;
			en[2][j*N*2+i*2] = (N+1)*j + i + N + 1;
			en[0][j*N*2+i*2+1] = (N+1)*j + i;
			en[1][j*N*2+i*2+1] = (N+1)*j + i + 1;
			en[2][j*N*2+i*2+1] = (N+1)*j + i + N + 2;
		}
	}
}

void triangle_cd(int N, CD& cd) {
	cd.resize(2, (N+1)*(N+1));
	for(int j = 0; j <= N; ++j) {
		for(int i = 0; i <= N; ++i) {
			cd[0][(N+1)*j+i] = i*1./(N);
			cd[1][(N+1)*j+i] = j*1./(N);
		}
	}
}

void rectangle_en(int N, EN& en) {
	en.resize(4, N*N);
	for(int j = 0; j < N; ++j) {
		for(int i = 0; i < N; ++i) {
			en[0][j*N+i] = (N+1)*j + i;
			en[1][j*N+i] = (N+1)*j + i + 1;
			en[2][j*N+i] = (N+1)*j + i + N + 2;
			en[3][j*N+i] = (N+1)*j + i + N + 1;
		}
	}
}

void rectangle_cd(int N, CD& cd) {
	cd.resize(2, (N+1)*(N+1));
	for(int j = 0; j <= N; ++j) {
		for(int i = 0; i <= N; ++i) {
			cd[0][(N+1)*j+i] = i*1./(N);
			cd[1][(N+1)*j+i] = j*1./(N);
		}
	}
}

#endif
