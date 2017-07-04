/* 
 * File:   alignmentRange.h
 * Author: daniel
 *
 * Created on 5 de febrero de 2016, 20:55
 */

#ifndef ALIGNMENTRANGE_H
#define	ALIGNMENTRANGE_H
typedef struct AlignmentRange_
{
   long long sp;
   long long ep;
   int qualityScore;
   char * pattern ;
   
} AlignmentRange;

AlignmentRange * createAlignmentRange(long long sp,long long ep, int quality, const char * pattern);
void destroyAlignmentRange(AlignmentRange * alignment);

#endif	/* ALIGNMENTRANGE_H */

