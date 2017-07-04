
#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include "alignmentRange.h"

AlignmentRange * createAlignmentRange(long long sp, long long ep, int quality, const char * pattern) {
    AlignmentRange * alignmentR = (AlignmentRange*) malloc(sizeof (AlignmentRange));
    alignmentR->sp = sp;
    alignmentR->ep = ep;
    alignmentR->qualityScore = quality;
    if (pattern) {
        alignmentR->pattern = (char*) malloc(strlen(pattern) + 1);
        strcpy(alignmentR->pattern, pattern);
    } else
        alignmentR->pattern = NULL;
    return alignmentR;

}

void destroyAlignmentRange(AlignmentRange * alignmentR) {
    if (alignmentR->pattern)
        free(alignmentR->pattern);
    free(alignmentR);
}

//void printAlignmentBowtie(ALIGNMENT * alignment){
//    printf("Position:%lld   Quality:%d",alignment->position,alignment->qualityScore);
//}
