/* 
 * File:   ranges.h
 * Author: daniel
 *
 * Created on 4 de diciembre de 2015, 21:36
 */

#ifndef RANGES_H
#define	RANGES_H
#include "alignmentRange.h"

typedef struct Ranges_
{
   long numRanges;
   AlignmentRange ** alignmentsR;   
} Ranges;

Ranges * createRanges(long long sp, long long ep,int quality, const char * pattern);
Ranges * joinRanges(Ranges * r1, Ranges * r2);
void destroyRanges(Ranges *r);
#ifdef	__cplusplus
extern "C" {
#endif




#ifdef	__cplusplus
}
#endif

#endif	/* RANGES_H */

