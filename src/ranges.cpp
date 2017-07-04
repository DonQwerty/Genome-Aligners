/* 
 * File:   ranges.cpp
 * Author: daniel
 *
 * Created on 3 de diciembre de 2015, 19:04
 */

#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include "ranges.h"

Ranges * createRanges(long long sp, long long ep, int quality, const char * pattern) {
    Ranges * ranges = (Ranges*) malloc(sizeof (Ranges));
    ranges->alignmentsR = (AlignmentRange**) malloc(sizeof (AlignmentRange *));
    ranges->numRanges = 1;
    ranges->alignmentsR[0] = createAlignmentRange(sp, ep, quality, pattern);
    return ranges;
}

Ranges * joinRanges(Ranges * r1, Ranges * r2) {
    if (!r1 && !r2)
        return NULL;
    else if (!r1) {
        return r2;
    } else if (!r2) {
        return r1;
    }
    Ranges * ranges = (Ranges*) malloc(sizeof (Ranges));
    int numRanges = r1->numRanges + r2->numRanges;
    ranges->numRanges = r1->numRanges;
    ranges->alignmentsR = (AlignmentRange **) malloc(numRanges * sizeof (AlignmentRange*));
   long long sp;
   long long ep;
   int quality;
   int i,j;
    for ( i = 0; i < r1->numRanges; i++) {
        ranges->alignmentsR[i] = r1->alignmentsR[i];
    }
    for ( i = 0; i < r2->numRanges; i++) {
        sp = r2->alignmentsR[i]->sp;
        ep = r2->alignmentsR[i]->ep;
        quality = r2->alignmentsR[i]->qualityScore;
        for ( j = 0; j < r1->numRanges; j++) {
            if(sp < ranges->alignmentsR[j]->sp){
                if(ep > ranges->alignmentsR[j]->sp){
                    if(ranges->alignmentsR[j]->qualityScore > quality){
                        destroyAlignmentRange(ranges->alignmentsR[j]);
                        ranges->alignmentsR[j] = r2->alignmentsR[i];
                        break;
                    }
                    else{
                        destroyAlignmentRange(r2->alignmentsR[i]);
                        break;
                    }
                }
            }
            else if(sp < ranges->alignmentsR[j]->ep){
                if(ranges->alignmentsR[j]->qualityScore > quality){
                        destroyAlignmentRange(ranges->alignmentsR[j]);
                        ranges->alignmentsR[j] = r2->alignmentsR[i];
                        break;
                }
                else{
                    destroyAlignmentRange(r2->alignmentsR[i]);
                    break;
                }
            }
        }
        if(j == r1->numRanges){
            ranges->alignmentsR[ranges->numRanges] = r2->alignmentsR[i];
            ranges->numRanges++;
        }
    }
    
    free(r1->alignmentsR);
    free(r2->alignmentsR);
    free(r1);
    free(r2);
    return ranges;

}

void destroyRanges(Ranges *r) {
    if (r != NULL) {
        for (int i = 0; i < r->numRanges; i++) {
            destroyAlignmentRange(r->alignmentsR[i]);
        }
        free(r->alignmentsR);
        free(r);
    }
}
