
#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include "alignment.h"

Alignment * createAlignment(long long pos, int quality, const char * pattern, char * id) {
    Alignment * alignment = (Alignment*) malloc(sizeof (Alignment));
    alignment->id = (char*) malloc(strlen(id) + 1);
    strcpy(alignment->id, id);
    alignment->pos = pos;
    alignment->qualityScore = quality;
    alignment->pattern = (char*) malloc(strlen(pattern) + 1);
    strcpy(alignment->pattern, pattern);
    return alignment;

}

void destroyAlignment(Alignment * alignment) {
    free(alignment->id);
    free(alignment->pattern);
    free(alignment);
}

void printAlignment(Alignment * alignment) {
    printf(" ID: %-*s\tPosition: %-*lld\tQuality: %-*d\t%s\n", 6,alignment->id,9,alignment->pos,2, alignment->qualityScore, alignment->pattern);
}

AlignmentSet * createAlignmentSet() {
    AlignmentSet * set = (AlignmentSet *) malloc(sizeof (AlignmentSet));
    set->nAlingment = 0;
    return set;
}

void addAlignment(AlignmentSet * set, Alignment * align) {
    if (!set) {
        set = (AlignmentSet *) malloc(sizeof (AlignmentSet));
        set->nAlingment = 0;
    }
    set->alignments[set->nAlingment] = align;
    set->nAlingment++;
}

void destroyAlignmentSet(AlignmentSet * set) {
    int i;
    if (!set)
        return;
    for (i = 0; i < set->nAlingment; ++i) {
        destroyAlignment(set->alignments[i]);
    }
    free(set);
}

int compareAlign(const void *arg1, const void *arg2) {
    Alignment ** align1 = (Alignment **) arg1;
    Alignment ** align2 = (Alignment **) arg2;
    return ( (*align1)->qualityScore - (*align2)->qualityScore);
}

void sortAlignmentSet(AlignmentSet * set) {
    qsort(set->alignments, set->nAlingment, sizeof (Alignment*), compareAlign);
}

void fprintAlignmentSet(FILE * f, AlignmentSet * set) {
    int i;
    if(set->nAlingment==0)
        return;
    sortAlignmentSet(set);
    for (i = 0; i < set->nAlingment; ++i) {
        printAlignment(set->alignments[i]);
        break;
    }
}