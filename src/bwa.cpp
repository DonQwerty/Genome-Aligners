#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <algorithm>
#include "bwa.h"
#include "alignment.h"


using namespace std;

unsigned int SEED ;
void saveBwaFMI(char * file , BwaFmi *bwaFmi){
    FILE * f;
    f = fopen(file,"w");
    saveFMI(f,bwaFmi->forwardFmi);
    saveFMI(f , bwaFmi->mirrorFmi);
    fclose(f);
}

BwaFmi * loadBwaFMI(char * file, int gap ){
    FILE * f;
    f = fopen(file,"r");
    BwaFmi * bwaFmi = (BwaFmi *)malloc(sizeof(BwaFmi));
    bwaFmi->gap = 0; 
	if(gap){
		bwaFmi->gap = gap;	
	}
    bwaFmi->forwardFmi = loadFMI(f);
    bwaFmi->mirrorFmi = loadFMI(f);
    bwaFmi->nNuc = bwaFmi->forwardFmi->nNuc;
    fclose(f);
    return bwaFmi;
}
Ranges * inexactMatchingRec(const char * pattern, int i,int * qualities,int q, int z, long long k2, long long l2, Fmi *fmi, int * d,char * patternAux, int gap) {
    
    if(i>=0)
     if (z < d[i])
        return NULL;
    if (q > 70)
        return NULL;
    if (i < 0) {
        return createRanges(k2, l2, q, patternAux);
    }
    
    Ranges * inter = NULL;
    if(gap){
		inter = joinRanges(inexactMatchingRec(pattern, i - 1,qualities, q+30, z - 1, k2, l2, fmi, d,patternAux,gap), inter);
		}
    char cadena[] = "GCAT";
    for (int j = 0; j < 4; j++) {
        char c = cadena[j];
        long long k = getMinRank(c, fmi->ranking) + getCheckpoints(fmi, k2, c) - 1;
        long long l = getMinRank(c, fmi->ranking) + getCheckpoints(fmi, l2, c) - 1;
        if (k < l) {
			if(gap){
           		Ranges * interDeletion = inexactMatchingRec(pattern, i,qualities, q+30, z - 1, k, l, fmi, d, patternAux,gap);
           		inter = joinRanges(inter, interDeletion);
			}
            if (c == pattern[i]) {
                Ranges * interWithoutMis = inexactMatchingRec(pattern, i - 1,qualities, q , z, k, l, fmi, d,patternAux, gap);
                inter = joinRanges(inter, interWithoutMis);
            } else {
                patternAux[i] = c;
                Ranges * interWithMis = inexactMatchingRec(pattern, i - 1,qualities, q + qualities[i], z - 1, k, l, fmi, d, patternAux, gap);
                patternAux[i] = pattern[i];
                inter = joinRanges(inter, interWithMis);
            }

        }
    }
    return inter;


}

void inexactMatching(Read read, int z, BwaFmi bwaFmi) {
    long long i,j, ep, sp;
    int *d;
    int q = 0;
    Alignment * align;
    AlignmentSet *set = createAlignmentSet();
    d = calculateD(read.sequenceReverse, bwaFmi);
    char * patternAux=strdup(read.sequence);
    Ranges * inters = inexactMatchingRec(read.sequence, strlen(read.sequence) - 1,read.quality,q, z, 0, bwaFmi.forwardFmi->nNuc, bwaFmi.forwardFmi, d,patternAux, bwaFmi.gap);
    if (inters != NULL) {
        for ( j = 0; j < inters->numRanges; j++) {
            sp = inters->alignmentsR[j]->sp+1;
            ep = inters->alignmentsR[j]->ep+1;
            q = inters->alignmentsR[j]->qualityScore;
            for (i = sp; i < ep; i++) {
                align = createAlignment(getPosition(*(bwaFmi.forwardFmi), i), q, inters->alignmentsR[j]->pattern,read.id);
                addAlignment(set, align);
            }

        }
    }
    fprintAlignmentSet(NULL, set);
    destroyAlignmentSet(set);
    free(d);
    free(patternAux);
    destroyRanges(inters);
}

BwaFmi * createBwaFMI(char * file, int gap) {
    BwaFmi * bwaFmi = (BwaFmi *) malloc(sizeof (BwaFmi));
    fflush(stdout);
    string genomeString;
    ifstream f;
    f.open(file);
    string line;
    getline(f, line);
    fflush(stdout);
    bwaFmi->nNuc = std::stoll(line.c_str(), NULL, 0);
    fflush(stdout);
    while (getline(f, line)) {
        genomeString.append(line);

    }
    genomeString.append("$");

    bwaFmi->forwardFmi = createFMI(genomeString.c_str(), bwaFmi->nNuc);
    reverse(genomeString.begin(), genomeString.end() - 1);
    bwaFmi->mirrorFmi = createFMI(genomeString.c_str(), bwaFmi->nNuc);
	if(gap){
		bwaFmi->gap = gap;	
	}
    //if(gap)
      //  printf("holaaa\n");
    f.close();
    return bwaFmi;

}

void removeBwaFmi(BwaFmi * bwaFmi) {
    removeFMI(bwaFmi->forwardFmi);
    removeFMI(bwaFmi->mirrorFmi);
    free(bwaFmi);
}

int * calculateD(char * pattern, BwaFmi bwaFmi) {
    long long k, l;
    char c;
    int z = 0,flag;
    unsigned int i;
    int * d = (int *) malloc(strlen(pattern) * sizeof (int));
    k = 0;
    l = bwaFmi.nNuc;
    for (i = 0; i <= (strlen(pattern) - 1); i++) {
        flag = 0;
        c = pattern[i];
        if(c=='N') {
            flag = 1;
        }
        else{
            k = getMinRank(c, bwaFmi.mirrorFmi->ranking) + getCheckpoints(bwaFmi.mirrorFmi, k, c) - 1;
            l = getMinRank(c, bwaFmi.mirrorFmi->ranking) + getCheckpoints(bwaFmi.mirrorFmi, l, c) - 1;
        }
        if (k > l || flag ==1) {
            k = 0;
            l = bwaFmi.nNuc;
            z++;
        }
        d[i] = z;
    }
    return d;
}
