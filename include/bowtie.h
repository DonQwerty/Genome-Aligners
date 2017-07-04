/* 
 * File:   bowtie.h
 * Author: daniel
 *
 * Created on 31 de enero de 2016, 19:05
 */
#include "fmi.h"
#include "ranges.h"
#include "alignment.h"
 #include "read.h"
#ifndef BOWTIE_H
#define	BOWTIE_H
#define MAX_QUALITY 70


typedef struct BowtieFmi_
{
    Fmi * forwardFmi;
    Fmi * mirrorFmi;
    long long nNuc;
   
} BowtieFmi;

void saveBowtieFMI(char * file , BowtieFmi *bowtieFmi);
BowtieFmi * loadBowtieFMI(char * file,int seed);
BowtieFmi * createBowtieFMI(char * file, int seed);
void updateSpEp(Fmi *fmi, std::string pattern,int length,long long *sp,long long *ep);
void exactMatchBowtie(BowtieFmi *bowtieFMI,std::string pattern);
Ranges * caseRecurrent(Fmi *fmi, std::string pattern,int * qualitySeq, int pos,int minPos, int mismatches, int quality, long long sp, long long ep, int flagExact);
void case1(BowtieFmi *fmiBow, std::string pattern,Read read, int mismatches,AlignmentSet * set);
void case2(BowtieFmi *fmiBow, std::string pattern, Read read, int mismatches,AlignmentSet * set);
void case3(BowtieFmi *fmiBow, std::string pattern, Read read,int mismatches,AlignmentSet * set);
void case4(BowtieFmi *fmiBow, std::string pattern, Read read,int mismatches,AlignmentSet * set);
Alignment * createAlignmentBowtie(Fmi *fmi, std::string pattern,int * qualitySeq, std::string text, long long pos, int quality,long long startAlign,char * id);
void inexactMatchBowtie(BowtieFmi *bowtieFMI, Read read, int mismatch);
void removeBowtieFMI(BowtieFmi * bowtieFMI);

#endif	/* BOWTIE_H */


