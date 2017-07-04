/* 
 * File:   alignment.h
 * Author: daniel
 *
 * Created on 5 de febrero de 2016, 19:07
 */

#ifndef ALIGNMENT_H
#define	ALIGNMENT_H



typedef struct Alignment_
{
	char * id;
   long long pos;
   int qualityScore;
   char * pattern;
   
} Alignment;

typedef struct AlignmentSet_
{
   int nAlingment;
   Alignment *  alignments[500];
   
} AlignmentSet;

Alignment * createAlignment(long long pos, int quality,const char * pattern, char * id);
void destroyAlignment(Alignment * alignment);
void printAlignment(Alignment * alignment);

AlignmentSet * createAlignmentSet();
void addAlignment(AlignmentSet * set, Alignment * align);
void destroyAlignmentSet(AlignmentSet * set);
int compareAlign(const void *arg1, const void *arg2);
void sortAlignmentSet(AlignmentSet * set);
void fprintAlignmentSet(FILE * f,AlignmentSet * set);


#endif	/* ALIGNMENT_H */

