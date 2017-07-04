/* 
 * File:   fmi.h
 * Author: daniel
 *
 * Created on 29 de octubre de 2015, 21:52
 */
#include "ranges.h"
#include "byte.h"
#ifndef FMI_H
#define	FMI_H

#define pCheckPoint 128
#define pSuffixArray 32
typedef struct Fmi_
{
   long long * ranking;
   long long ** checkpoints;
   long long positionDollar;
   u_int8_t * bwt;
   long long * positions;
   long long  nNuc;
   
} Fmi;
long long * setPositions(long long * ranking, long long nNuc);
Fmi * createFMI(const char *text,long long nNuc);
void updateCheckpoints(char c, long long pos,long long ** checkpoints);
void updateRanking(char c,long long * ranking);
long long getMinRank(char c,long long * ranking);
long long getMaxRank(char c,long long * ranking);
char getCharacterBWT(Fmi fmi, long long i);
long long getCheckpoints(Fmi * fmi, long long pos, char c);
//void findSubString(FMI *fmi,const char * subS);
void pathString(long long pos,Fmi * fmi,int length,char *chain);
long long getPosition(Fmi fmi, long long i);
void printPosition(Fmi fmi, long long i, FILE * f);

void removeFMI(Fmi * fmi);
void saveFMI(FILE * f , Fmi *fmi);
Fmi * loadFMI(FILE * f );


#endif	/* FMI_H */

