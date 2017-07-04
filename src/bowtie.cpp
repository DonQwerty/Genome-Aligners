/* 
 * File:   bowtie.cpp
 * Author: daniel
 *
 * Created on 31 de enero de 2016, 19:04
 */

#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <algorithm>
#include <math.h>
#include <iostream>

#include "bowtie.h"

using namespace std;
unsigned int SEED = 28;
void saveBowtieFMI(char * file , BowtieFmi *bowtieFmi){
    FILE * f;
    f = fopen(file,"w");
    saveFMI(f,bowtieFmi->forwardFmi);
    saveFMI(f , bowtieFmi->mirrorFmi);
    fclose(f);
}

BowtieFmi * loadBowtieFMI(char * file,int seed){
    FILE * f;
    f = fopen(file,"r");
    BowtieFmi * bowtieFmi = (BowtieFmi *)malloc(sizeof(BowtieFmi));
    bowtieFmi->forwardFmi = loadFMI(f);
    bowtieFmi->mirrorFmi = loadFMI(f);
    if(seed){
        SEED = seed;
    }
    bowtieFmi->nNuc = bowtieFmi->forwardFmi->nNuc;
    fclose(f);
    return bowtieFmi;
}

BowtieFmi * createBowtieFMI(char * file,int seed) {
    BowtieFmi * bowtieFMI = (BowtieFmi *) malloc(sizeof (BowtieFmi));
    string genomeString;
    ifstream f;
    f.open(file);
    string line;
    getline(f, line);
    if(seed){
        SEED = seed;   
    }
    bowtieFMI->nNuc = std::stoll(line.c_str(), NULL, 0);
    fflush(stdout);
    while (getline(f, line)) {
        genomeString.append(line);

    }
    genomeString.append("$");

    bowtieFMI->forwardFmi = createFMI(genomeString.c_str(), bowtieFMI->nNuc);
    std::reverse(genomeString.begin(), genomeString.end() - 1);
    bowtieFMI->mirrorFmi = createFMI(genomeString.c_str(), bowtieFMI->nNuc);
    f.close();
    return bowtieFMI;

}

void updateSpEp(Fmi *fmi, std::string pattern, int length, long long *sp, long long *ep) {
    int i = pattern.size() - 1;
    char c = pattern.at(i);
    if(c == 'N'){
            *sp = 0;
            *ep = 0;
            return;
        }
    long long spAux = getMinRank(c, fmi->ranking) - 1;
    long long epAux = getMaxRank(c, fmi->ranking) - 1;
    i--;
    while (spAux < epAux && i >= length) {
        c = pattern.at(i);
        if(c == 'N'){
            spAux = 0;
            epAux = 0;
            break;
        }
        spAux = getMinRank(c, fmi->ranking) + getCheckpoints(fmi, spAux, c) - 1;
        epAux = getMinRank(c, fmi->ranking) + getCheckpoints(fmi, epAux, c) - 1;
        i--;
    }
    *sp = spAux;
    *ep = epAux;
}

void exactMatchBowtie(BowtieFmi *bowtieFMI, std::string pattern) {
    FILE * f = fopen("salida2.txt", "w");
    long long i = 0;
    long long sp;
    long long ep;
    updateSpEp(bowtieFMI->forwardFmi, pattern, 0, &sp, &ep);
    for (i = sp + 1; i < ep + 1; i++) {
        printPosition(*(bowtieFMI->forwardFmi), i, f);

    }
    fclose(f);
}

char * aleatCombination() {
    time_t t;
    srand((unsigned) time(&t));
    int i;
    char * chain = (char*) malloc(4);
    int nuc[4] = {0, 1, 2, 3};
    int r;
    for (i = 0; i <= 3; i++) {
        r = rand() % (4 - i);
        chain[i] = getCharNuc(nuc[r + i]);
        nuc[r + i] = nuc[i];
    }
    return chain;

}

Ranges * caseRecurrent(Fmi *fmi, string pattern, int * qualitySeq, int pos, int minPos, int mismatches, int quality, long long sp, long long ep, int flagExact) {
    if (pos <= minPos) {
        if ((mismatches == 0 && flagExact == -1) || (mismatches >= 0 && flagExact == 0) || (mismatches < flagExact && flagExact > 0)) {
            Ranges * range = createRanges(sp + 1, ep + 1, quality, pattern.substr(pattern.length() - SEED).c_str());
            return range;
        } else
            return NULL;
    }
    Ranges * inter = NULL;
    if (mismatches > 0) {
        char c = pattern.at(pos);
        string patternAux(pattern);
        Ranges * aux = NULL;
        long long spAux;
        long long epAux;
        int qualityAux;
        char * chain;
        chain = aleatCombination();
        char nuc;
        for (int i = 0; i < 4; i++) {
            nuc = chain[i];
            spAux = getMinRank(nuc, fmi->ranking) + getCheckpoints(fmi, sp, nuc) - 1;
            epAux = getMinRank(nuc, fmi->ranking) + getCheckpoints(fmi, ep, nuc) - 1;
            aux = NULL;
            if (spAux < epAux) {
                if (nuc == c) {
                    aux = caseRecurrent(fmi, pattern, qualitySeq, pos - 1, minPos, mismatches, quality, spAux, epAux, flagExact);
                } else {
                    qualityAux = quality + qualitySeq[pos];
                    if (qualityAux < MAX_QUALITY) {
                        patternAux.erase(pos, 1);
                        patternAux.insert(pos, 1, nuc);
                        aux = caseRecurrent(fmi, patternAux, qualitySeq, pos - 1, minPos, mismatches - 1, qualityAux, spAux, epAux, flagExact);
                    }
                }
                inter = joinRanges(aux, inter);
            }
        }
        free(chain);
    } else {
        char c = pattern.at(pos);
        if(c != 'N'){
            sp = getMinRank(c, fmi->ranking) + getCheckpoints(fmi, sp, c) - 1;
            ep = getMinRank(c, fmi->ranking) + getCheckpoints(fmi, ep, c) - 1;
            if (sp < ep)
                inter = caseRecurrent(fmi, pattern, qualitySeq, pos - 1, minPos, mismatches, quality, sp, ep, flagExact);
        }
    }
    return inter;
}

Alignment * createAlignmentBowtie(Fmi *fmi, std::string pattern, int * qualitySeq, std::string text, long long pos, int quality, long long startAlign, char * id) {
    int flag = 0;
    if (SEED != pattern.length()) {
        char * suffix = (char *) malloc((int) pattern.length() - SEED);
        pathString(pos, fmi, pattern.length() - SEED, suffix);
        for (int k = pattern.length() - SEED - 1; k >= 0; k--) {
            if (suffix[k] == pattern.at(k)) {
                text.push_back(pattern.at(k));
            } else {
                quality += qualitySeq[k];
                if (quality > MAX_QUALITY) {
                    flag = 1;
                    break;
                }
                text.push_back(suffix[k]);
            }
        }
        free(suffix);
    }
    if (flag == 0) {
        Alignment * align = createAlignment(startAlign, quality, text.c_str(), id);
        return align;
    } else
        return NULL;
}

void case1(BowtieFmi *fmiBow, std::string pattern, Read read, int mismatches, AlignmentSet * set) {
  //  printf("Case 1\n");
    Fmi * fmi = fmiBow->mirrorFmi;
    long long sp;
    long long ep;
    updateSpEp(fmi, pattern, pattern.length() - SEED / 2, &sp, &ep);
    if(sp>=ep)
        return;
    Ranges * ranges = caseRecurrent(fmi, pattern, getQualityRev(read), pattern.length() - SEED / 2 - 1, pattern.length() - SEED - 1, mismatches, 0, sp, ep, 0);
    if (!ranges)
        return;
    for (int i = 0; i < ranges->numRanges; i++) {
        sp = ranges->alignmentsR[i]->sp;
        ep = ranges->alignmentsR[i]->ep;
        for (long long j = sp; j < ep; j++) {
            int quality = ranges->alignmentsR[i]->qualityScore;
            long long posSuff = llabs(getPosition(*fmi, j) - fmiBow->nNuc + 1);

            if ((fmiBow->nNuc - posSuff) >= (int) pattern.length() - SEED + 1) {
                string text = "";
                text.append(ranges->alignmentsR[i]->pattern);
                reverse(text.begin(), text.end());
                Alignment * align = createAlignmentBowtie(fmi, pattern, getQualityRev(read), text, j, quality, posSuff - SEED + 1,read.id);
                if (!align)
                    continue;
                //printAlignment(align);
                addAlignment(set, align);
            }
        }
    }
    destroyRanges(ranges);

    return;
}

void case2(BowtieFmi *fmiBow, std::string pattern, Read read, int mismatches, AlignmentSet * set) {
 //   printf("Case 2\n");
    Fmi * fmi = fmiBow->forwardFmi;
    long long sp;
    long long ep;
    updateSpEp(fmi, pattern.substr(0, SEED), SEED / 2, &sp, &ep);
    if(sp>=ep)
        return;
    Ranges * ranges = caseRecurrent(fmi, pattern.substr(0, SEED), getQuality(read), SEED / 2 - 1, -1, mismatches, 0, sp, ep, mismatches);
    if (!ranges)
        return;
    for (int i = 0; i < ranges->numRanges; i++) {
        string seedPortion(ranges->alignmentsR[i]->pattern);
        reverse(seedPortion.begin(), seedPortion.end());
        updateSpEp(fmiBow->mirrorFmi, seedPortion, 0, &sp, &ep);
        reverse(pattern.begin(), pattern.end());
        for (long long j = sp + 1; j < ep + 1; j++) {

            int quality = ranges->alignmentsR[i]->qualityScore;
            long long posForw = llabs(getPosition(*fmiBow->mirrorFmi, j) - fmiBow->nNuc + 1);

            if ((fmiBow->nNuc - posForw - 1) >= (int) pattern.length() - SEED) {
                string text = "";
                text.append(ranges->alignmentsR[i]->pattern);
                Alignment * align = createAlignmentBowtie(fmiBow->mirrorFmi, pattern, getQualityRev(read), text, j, quality, posForw - SEED + 1,read.id);
                if (!align)
                    continue;
              //  printAlignment(align);
                addAlignment(set, align);
            }
        }
    }
    destroyRanges(ranges);

    return;
}

void case3(BowtieFmi *fmiBow, std::string pattern, Read read, int mismatches, AlignmentSet * set) {
  //  printf("Case 3\n");
    Fmi * fmi = fmiBow->mirrorFmi;
    long long sp = 0;
    long long ep = fmiBow->nNuc;
    int quality;
    Ranges * ranges1 = caseRecurrent(fmi, pattern, getQualityRev(read), pattern.length() - 1, pattern.length() - SEED / 2 - 1, 1, 0, sp, ep, -1);
    if (!ranges1)
        return;
    for (int i = 0; i < ranges1->numRanges; i++) {
        sp = ranges1->alignmentsR[i]->sp - 1;
        ep = ranges1->alignmentsR[i]->ep - 1;
        if(sp>=ep)
            break;
        quality = ranges1->alignmentsR[i]->qualityScore;
        string seed(ranges1->alignmentsR[i]->pattern);
        pattern.replace(pattern.length() - SEED, SEED, seed);
        Ranges * ranges2 = caseRecurrent(fmi, pattern, getQualityRev(read), pattern.length() - SEED / 2 - 1, pattern.length() - SEED - 1, mismatches - 1, quality, sp, ep, mismatches - 1);
        if (!ranges2)
            continue;
        for (int j = 0; j < ranges2->numRanges; j++) {
            sp = ranges2->alignmentsR[j]->sp;
            ep = ranges2->alignmentsR[j]->ep;
            for (long long k = sp; k < ep; k++) {

                quality = ranges2->alignmentsR[j]->qualityScore;
                long long posSuff = llabs(getPosition(*fmi, k) - fmiBow->nNuc + 1);

                if ((fmiBow->nNuc - posSuff) >= (int) pattern.length() - SEED + 1) {
                    string text = "";
                    text.append(ranges2->alignmentsR[j]->pattern);
                    reverse(text.begin(), text.end());
                    Alignment * align = createAlignmentBowtie(fmiBow->mirrorFmi, pattern, getQualityRev(read), text, k, quality, posSuff - SEED + 1,read.id);
                    if (!align)
                        continue;
                  //  printAlignment(align);
                    addAlignment(set, align);
                }
            }
        }
        destroyRanges(ranges2);
    }
    destroyRanges(ranges1);

    return;
}

void case4(BowtieFmi *fmiBow, std::string pattern, Read read, int mismatches, AlignmentSet * set) {
  //  printf("Case 4\n");
    Fmi * fmi = fmiBow->forwardFmi;
    long long sp = 0;
    long long ep = fmiBow->nNuc;
    int quality;
    Ranges * ranges1 = caseRecurrent(fmi, pattern.substr(0, SEED), getQuality(read), SEED - 1, SEED / 2 - 1, 1, 0, sp, ep, -1);
    if (!ranges1)
        return;
    reverse(pattern.begin(), pattern.end());
    for (int i = 0; i < ranges1->numRanges; i++) {
        sp = ranges1->alignmentsR[i]->sp - 1;
        ep = ranges1->alignmentsR[i]->ep - 1;
        quality = ranges1->alignmentsR[i]->qualityScore;
        string seed(ranges1->alignmentsR[i]->pattern);
        Ranges * ranges2 = caseRecurrent(fmi, seed, getQuality(read), SEED / 2 - 1, -1, 2, quality, sp, ep, -1);
        if (!ranges2)
            continue;
        for (int j = 0; j < ranges2->numRanges; j++) {
            sp = ranges2->alignmentsR[j]->sp;
            ep = ranges2->alignmentsR[j]->ep;
            string seedPortion(ranges2->alignmentsR[j]->pattern);
            reverse(seedPortion.begin(), seedPortion.end());
            updateSpEp(fmiBow->mirrorFmi, seedPortion, 0, &sp, &ep);
            for (long long k = sp + 1; k < ep + 1; k++) {

                quality = ranges2->alignmentsR[j]->qualityScore;
                long long posForw = llabs(getPosition(*fmiBow->mirrorFmi, k) - fmiBow->nNuc + 1);

                if ((fmiBow->nNuc - posForw - 1) >= (int) pattern.length() - SEED) {
                    string text = "";
                    text.append(ranges2->alignmentsR[j]->pattern);
                    Alignment * align = createAlignmentBowtie(fmiBow->mirrorFmi, pattern, getQualityRev(read), text, k, quality, posForw - SEED + 1,read.id);
                    if (!align)
                        continue;
                 //   printAlignment(align);
                    addAlignment(set, align);
                }
            }
        }
        destroyRanges(ranges2);
    }
    destroyRanges(ranges1);

    return;
}

void inexactMatchBowtie(BowtieFmi *bowtieFMI, Read read, int mismatch) {
    //FILE * f = fopen("salida2.txt", "w");
    //    long long i = 0;
    // long long j = 0;
    //  RANGES * ranges = NULL;
    //  int numAlignments = 0;
    string pattern;
    pattern=string(read.sequence);
    string patternReverse(pattern);
    AlignmentSet *set = createAlignmentSet();


    reverse(patternReverse.begin(), patternReverse.end());
    switch (mismatch) {
        case 0:

            case1(bowtieFMI, patternReverse, read, 0, set);
            break;
        case 1:
            case1(bowtieFMI, patternReverse, read, 1, set);
            case2(bowtieFMI, pattern, read, 1, set);
            break;
        case 2:
            case1(bowtieFMI, patternReverse, read, 2, set);
            case2(bowtieFMI, pattern, read, 2, set);
            case3(bowtieFMI, patternReverse, read, 2, set);
            break;
        case 3:
            case1(bowtieFMI, patternReverse, read, 3, set);
            case2(bowtieFMI, pattern, read, 3, set);
            case3(bowtieFMI, patternReverse, read, 3, set);
            case4(bowtieFMI, pattern, read, 3, set);
            break;

    }
    fprintAlignmentSet(NULL, set);
    destroyAlignmentSet(set);

    //free(range);
    //for (i = 0; i < numAlignments; i++) {
    // printAlignmentBowtie(alignments[i]);
    //destroyAlignmentBowtie(alignments[i]);
    // }
    // free(alignments);
    //fclose(f);


}

void removeBowtieFMI(BowtieFmi * bowtieFMI) {
    removeFMI(bowtieFMI->forwardFmi);
    removeFMI(bowtieFMI->mirrorFmi);
    free(bowtieFMI);
}

