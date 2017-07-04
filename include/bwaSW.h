/* 
 * File:   bwaSW.h
 * Author: daniel
 *
 * Created on 1 de abril de 2016, 19:05
 */
#include "fmi.h"
#include "ranges.h"
#include "alignment.h"
#include "read.h"
#include "dawg.h"
#ifndef BWASH_H
#define	BWASH_H

typedef struct BwaSW_
{
    Fmi * fmi;
    long long nNuc;
   
} BwaSW;
/*
typedef struct ColumnDP_
{
	ColumnDP_ ** childs;
	int * G;
	int * D;
	int * I;
	
}ColumnDP;

BwaSW * createBwaSW(char * file);
void removeBwaSW(BwaSW * bwaSW) ;
void saveBwaSW(char * file , BwaSW *bwaSW);
BwaSW * loadBwaSW(char * file);
void alignReadBwaSW(BwaSW * bwaSW, Read * read);

void alignReadBwaSWv2(BwaSW * bwaSW, Read * read);

ColumnDP * createComlumnDP(int length);
void removeColumnDP(ColumnDP * c);
void removeColumnDPTree(ColumnDP *c);

*/
#endif	/* BWASW_H */

