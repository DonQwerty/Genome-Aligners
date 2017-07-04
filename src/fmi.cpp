
/* 
 * File:   fmi.cpp
 * Author: daniel
 *
 * Created on 20 de octubre de 2015, 19:04
 */

#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>


#include "fmi.h"

const char * originalChain;

int compare(const void *arg1, const void *arg2) {
    return (strcmp(&originalChain[*(long long *) arg1], &originalChain[*(long long *) arg2]));
}

long long * setPositions(long long * positions, long long nNuc) {
    long long newpos = 0;
    long long * fmiPositions = (long long*) malloc(sizeof (long long)*(nNuc / pSuffixArray + 1));
    int i;
    for (i = 0; i < nNuc; i++) {
        if (i % pSuffixArray == 0) {
            newpos = i / pSuffixArray;
            fmiPositions[newpos] = positions[i];
        }

    }
    return fmiPositions;

}

void updateCheckpoints(long long * check, long long pos, long long ** checkpoints) {
    int i;
    long long newpos = pos / pCheckPoint;

    for (i = 0; i < 4; i++) {
        checkpoints[i][newpos] = check[i];

    }


}

void updateRanking(char c, long long * ranking) {
    ranking[getIntegerNuc(c)]++;
}

long long getMinRank(char c, long long * ranking) {
    long long res = 0;
    switch (c) {
        case 'A':
            res = 1;
            break;
        case 'C':
            res = 1 + ranking[0];
            break;
        case 'G':
            res = 1 + ranking[0] + ranking[1];
            break;
        default:
            res = 1 + ranking[0] + ranking[1] + ranking[2];
            break;
    }
    return res;
}

long long getMaxRank(char c, long long * ranking) {

    return getMinRank(c, ranking) + ranking[getIntegerNuc(c)];
}

char getCharacterBWT(Fmi fmi, long long i) {
    if (i == fmi.positionDollar)
        return '$';
    long long numByte, positionInByte;
    i = (i < fmi.positionDollar) ? i : i - 1;
    numByte = i / 4.0;
    positionInByte = i % 4;
    return getReverseConversion(fmi.bwt[numByte], positionInByte);

}

long long * sortPositions(long long nNuc) {
    long long * positions = (long long*) malloc((nNuc + 1) * sizeof (long long));

    for (long long i = 0; i < nNuc + 1; i++) {
        positions[i] = i;
    }
    qsort(positions, nNuc + 1, sizeof (long long), compare);
    return positions;
}

Fmi * createFMI(const char *text, long long nNuc) {
    Fmi * fmi = (Fmi *) malloc(sizeof (Fmi));
    originalChain = text;
    long long *positions = sortPositions(nNuc);

    fmi->ranking = (long long*) malloc(sizeof (long long)*4);
    for (int i = 0; i < 4; i++) {
        fmi->ranking[i] = 0;
    }
    fmi->positionDollar = 0;
    fmi->checkpoints = (long long**) malloc(sizeof (long long*)*4);
    for (int i = 0; i < 4; i++) {
        fmi->checkpoints[i] = (long long *) malloc(sizeof (long long)*((nNuc + 1) / pCheckPoint + 1));
    }
    fmi->positions = setPositions(positions, nNuc);
    fmi->nNuc = nNuc + 1;
    long long * check = (long long *) malloc(4 * sizeof (long long));
    for (int i = 0; i < 4; i++) {
        check[i] = 0;
    }

    long long i, cont = 0, numByte = 0;
    long long size;
    char c;
    u_int8_t byte;
    int positionInByte;
    size = (long long) ceil(nNuc / 4.0);
    fmi->bwt = (u_int8_t *) malloc(size);
    memset(fmi->bwt, 0, size);
    for (i = 0; i < nNuc + 1; i++) {
        if (positions[i] != 0) {
            c = text[positions[i] - 1];
            updateRanking(c, fmi->ranking);
            check[getIntegerNuc(c)]++;
            if (i % pCheckPoint == 0)
                updateCheckpoints(check, i, fmi->checkpoints);
            numByte = cont / 4.0;
            positionInByte = cont % 4;
            byte = getConversion(c, positionInByte);
            fmi->bwt[numByte] = fmi->bwt[numByte] | byte;
            cont++;
        } else {
            if (i % pCheckPoint == 0)
                updateCheckpoints(check, i, fmi->checkpoints);
            fmi->positionDollar = i;
        }
    }
    free(check);
    free(positions);
    return fmi;
}

void pathString(long long pos, Fmi * fmi, int length, char *chain) {

    char c = getCharacterBWT(*fmi, pos);
    if (length == 1) {
        chain[length - 1] = c;
        return;
    } else {
        chain[length - 1] = c;
        pathString(getCheckpoints(fmi, pos, c) + getMinRank(c, fmi->ranking) - 1, fmi, length - 1, chain);
        return;
    }
}

//void findSubString(FMI *fmi, const char * subS) {
//    long long i;
//    char c;
//    c = subS[0];
//    long long min = getMinRank(c, fmi->ranking);
//    long long max = getMaxRank(c, fmi->ranking);
//    for (i = min; i < max; i++) {
//        pathString(i, fmi, &subS[1]);
//    }
//}

long long getCheckpoints(Fmi * fmi, long long pos, char c) {
    int n = 0;
    long long numByte, positionInByte;
    long long newpos = pos / pCheckPoint;
    if(pos<0)
        return 0;
    if (pos % pCheckPoint == 0)
        return fmi->checkpoints[getIntegerNuc(c)][newpos];
    if(pos >= fmi->positionDollar)
        pos--;
    long long startPos = newpos * pCheckPoint + 1;
    if (startPos > fmi->positionDollar) {
        startPos--;
    }

    for (long long i = (startPos); i <= pos; i++) {
        numByte = i / 4.0;
        positionInByte = i % 4;
        char c2 = getReverseConversion(fmi->bwt[numByte], positionInByte);
        n += (c2 == c) ? 1 : 0;
    }
    return fmi->checkpoints[getIntegerNuc(c)][newpos] + n;
}

long long getPosition(Fmi fmi, long long i) {

    if (i % pSuffixArray == 0) {
        return fmi.positions[i / pSuffixArray];
    }
    char c = getCharacterBWT(fmi, i);
    if (c == '$')
        return 0;
    i = getMinRank(c, fmi.ranking) + getCheckpoints(&fmi, i, c) - 1;
    return getPosition(fmi, i) + 1;


}

void printPosition(Fmi fmi, long long i, FILE * f) {
    fprintf(f, "Posicion: %lld\n", getPosition(fmi, i));

}

void removeFMI(Fmi * fmi) {
    free(fmi->bwt);
    free(fmi->ranking);
    for (int i = 0; i < 4; i++) {
        free(fmi->checkpoints[i]);
    }
    free(fmi->checkpoints);
    free(fmi->positions);
    free(fmi);
}

void saveFMI(FILE * f , Fmi *fmi){
    int i;
    fwrite(fmi, sizeof(Fmi),1,f);
    fwrite(fmi->ranking,sizeof(long long),4,f);
    for ( i = 0; i < 4; i++){
        fwrite(fmi->checkpoints[i],sizeof(long long),((fmi->nNuc + 1) / pCheckPoint + 1),f);
    }
    fwrite(fmi->bwt,sizeof(u_int8_t),(long long) ceil(fmi->nNuc / 4.0),f);
    fwrite(fmi->positions,sizeof(long long), (fmi->nNuc / pSuffixArray + 1), f);
}

Fmi * loadFMI(FILE * f ){
        int i;
    Fmi * fmi = (Fmi *)malloc(sizeof(Fmi));

    fread(fmi, sizeof(Fmi),1,f);
    fmi->ranking = (long long *)malloc(sizeof(long long)*4);
    fread(fmi->ranking,sizeof(long long),4,f);
    fmi->checkpoints = (long long **)malloc(sizeof(long long*)*4);
    for ( i = 0; i < 4; i++){
        fmi->checkpoints[i] = (long long *) malloc(sizeof (long long)*((fmi->nNuc + 1) / pCheckPoint + 1));
        fread(fmi->checkpoints[i],sizeof(long long),((fmi->nNuc + 1) / pCheckPoint + 1),f);
    }
    fmi->bwt = (u_int8_t *)malloc(sizeof(u_int8_t)*(long long) ceil(fmi->nNuc / 4.0));
    fread(fmi->bwt,sizeof(u_int8_t),(long long) ceil(fmi->nNuc / 4.0),f);
    fmi->positions = (long long *)malloc(sizeof(long long)*(fmi->nNuc / pSuffixArray + 1));
    fread(fmi->positions,sizeof(long long), (fmi->nNuc / pSuffixArray + 1), f);
    return fmi;
}
