/* 
 * File:   bwtSW.h
 * Author: daniel
 *
 * Created on 1 de abril de 2016, 19:05
 */
#include "fmi.h"
#include "ranges.h"
#include "alignment.h"
#include "read.h"
#include "dawg.h"
#ifndef BWTSW_H
#define	BWTSW_H

typedef struct BwtSW_
{
    Fmi * fmi;
    long long nNuc;
   
} BwtSW;

typedef struct RowDP_
{
	int * G;
	int * D;
	int * I;
	
}RowDP;

BwtSW * createBwtSW(char * file);
void removeBwtSW(BwtSW * bwtSW) ;
void saveBwtSW(char * file , BwtSW *bwtSW);
BwtSW * loadBwtSW(char * file);
void alignReadBwtSW(BwtSW * bwtSW, Read * read);

void alignReadBwtRec(BwtSW * bwtSW, Dawg * dawg, long long sp, long long ep, int ** I, int ** D, int **G,int depth, long long * sa, int * maxV);
RowDP * createRowDP(int length);
void removeRowDP(RowDP * r,int depth);
void removeRowDPTree(RowDP *r);


#endif	/* BWTSW_H */

