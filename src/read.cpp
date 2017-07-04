
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "read.h"

Read * createRead(const char * id, const char * sequence, const char * quality) {
    unsigned int i;
    int j;
    Read * read = (Read*) malloc(sizeof (Read));
    read->id = (char*) malloc(strlen(id) + 1);
    read->sequence = (char*) malloc(strlen(sequence) + 1);
    read->quality = (int*) malloc(sizeof(int)*strlen(quality));
    read->sequenceReverse = (char*) malloc(strlen(sequence) + 1);
    read->qualityReverse = (int*) malloc(sizeof(int)*strlen(quality));
    strcpy(read->id, id);
    strcpy(read->sequence, sequence);
    j = strlen(quality) - 1;
    for (i = 0; i < strlen(quality); i++) {
        read->quality[i] = ((int) quality[i]) - 33;
        read->qualityReverse[j] = ((int) quality[i]) - 33;
        j--;
    }
    for (j = strlen(sequence) - 1; j >= 0; j--) {
        read->sequenceReverse[j] = sequence[strlen(quality) - 1 - j];
    }
    read->sequenceReverse[strlen(sequence)] = '\0';
    return read;
}

void destroyRead(Read * read) {
	free(read->id);
    free(read->sequence);
    free(read->quality);
    free(read->sequenceReverse);
    free(read->qualityReverse);
    free(read);
}

int * getQuality(Read read) {
    return read.quality;
}

char * getSequence(Read read) {
    return read.sequence;
}

int * getQualityRev(Read read) {
    return read.qualityReverse;
}

char * getSequenceRev(Read read) {
    return read.sequenceReverse;
}

char * getId(Read read) {
    return read.id;
}



