/* 
 * File:   bwa.h
 * Author: daniel
 *
 * Created on 31 de enero de 2016, 19:16
 */
#include "fmi.h"
 #include "alignmentRange.h"
 #include "read.h"
#ifndef BWA_H
#define	BWA_H
typedef struct BwaFmi_
{
    Fmi * forwardFmi;
    Fmi * mirrorFmi;
    long long nNuc;
	int gap;   
	
} BwaFmi;
void saveBwaFMI(char * file , BwaFmi *bwaFmi);
BwaFmi * loadBwaFMI(char * file, int gap);
Ranges * inexactMatchingRec(const char * pattern, int i,int * qualities,int q, int z, long long k2, long long l2, Fmi *fmi, int * d,char * patternAux, int gap);
void inexactMatching(Read read, int z, BwaFmi bwaFmi);
BwaFmi * createBwaFMI(char * file, int gap);
void removeBwaFmi(BwaFmi * bwaFmi) ;
int * calculateD(char * pattern, BwaFmi  bwaFmi);

#endif	/* BWA_H */

