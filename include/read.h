/* 
 * File:   read.h
 * Author: daniel
 *
 * Created on 6 de febrero de 2016, 16:47
 */

#ifndef READ_H
#define	READ_H
typedef struct Read_
{
    char * id;
    char * sequence;
    char * sequenceReverse;
    int * quality;
    int * qualityReverse;
   
} Read;

Read * createRead(const char * id,const char * sequence,const char * quality);

void destroyRead(Read * read);

int * getQuality(Read read);

char * getSequence(Read read);

int * getQualityRev(Read read);

char * getSequenceRev(Read read);


char * getId(Read read);


#endif	/* READ_H */

