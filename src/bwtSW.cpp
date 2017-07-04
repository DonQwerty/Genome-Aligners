
#include <iostream>

#include <math.h>
#include <algorithm>
#include "bwtSW.h"
#include "dawg.h"
#include "alignment.h"
#include <algorithm>
#include <stack>          // std::stack


using namespace std;
int s;
BwtSW * createBwtSW(char * file){

	BwtSW * bwtSW = (BwtSW *) malloc(sizeof (BwtSW));
    fflush(stdout);
    string genomeString;
    ifstream f;
    f.open(file);
    string line;
    getline(f, line);
    fflush(stdout);
    bwtSW->nNuc = std::stoll(line.c_str(), NULL, 0);
    fflush(stdout);
    while (getline(f, line)) {
        genomeString.append(line);

    }
    genomeString.append("$");

    bwtSW->fmi = createFMI(genomeString.c_str(), bwtSW->nNuc);
    
    f.close();
    return bwtSW;

}

void removeBwtSW(BwtSW * bwtSW) {
    removeFMI(bwtSW->fmi);
    free(bwtSW);
}

void saveBwtSW(char * file , BwtSW *bwtSW){
    FILE * f;
    f = fopen(file,"w");
    saveFMI(f,bwtSW->fmi);
    fclose(f);
}

BwtSW * loadBwtSW(char * file){
    FILE * f;
    f = fopen(file,"r");
    BwtSW * bwtSW = (BwtSW *)malloc(sizeof(BwtSW));
	
    bwtSW->fmi = loadFMI(f);
    bwtSW->nNuc = bwtSW->fmi->nNuc;
    fclose(f);
    return bwtSW;
}

RowDP * createRowDP(int length, int depth){
	RowDP * r = (RowDP *) malloc(sizeof(RowDP));
	if(!r)
		printf("ERROR\n");

	r->G = (int *) malloc(sizeof(int)*length);
	if(!r->G)
		printf("ERROR\n");
	r->D = (int *) malloc(sizeof(int)*length);
	if(!r->D)
		printf("ERROR\n");
	r->I = (int *) malloc(sizeof(int)*length);
	if(!r->I)
		printf("ERROR\n");
	r->G[0]=-5 + depth*-2;
	r->D[0]=-999;
	return r;
}



void removeRowDP(RowDP * c){
	printf("Eliminando\n");
	if(c){
		free(c->G);
		free(c->D);
		free(c->I);
		free(c);
	}
}
/*
void removeRowDPTree(RowDP *c){
	int i;
	for(i = 0; i< 4; i++){
		if(c->childs[i]!=NULL) {
			printf("holis\n");
			removeRowDPTree(c->childs[i]);
		}
	}
	removeRowDP(c);
}
*/


void alignReadBwtRec(BwtSW * bwtSW, Dawg * dawg, long long sp, long long ep, int ** I, int ** D, int **G,int depth, long long * sa, int * maxV){
	int value,pos,maxim;
	long long spAux, epAux;
	int i,n,j,k;
	int flag;
	RowDP *r ;
	char nuc[]="ACGT";
	for(i = 0; i < 4;i++){
		if(depth == 1){
			spAux = getMinRank(nuc[i], bwtSW->fmi->ranking) - 1;
			epAux = getMaxRank(nuc[i], bwtSW->fmi->ranking) - 1;
		}
		else{
			spAux = getMinRank(nuc[i], bwtSW->fmi->ranking) + getCheckpoints(bwtSW->fmi, sp, nuc[i]) -1;
        	epAux = getMinRank(nuc[i], bwtSW->fmi->ranking) + getCheckpoints(bwtSW->fmi, ep, nuc[i]) -1;
		}

		if(spAux < epAux){
			flag = 0;
			for(j = 1; j < dawg->length ; j++){
    			maxim = -999;
    			for(k = 0; k < dawg->arrayNodes[j]->NFathers; k++){
    				n = dawg->arrayNodes[j]->fathers[k];

    				value = G[depth-1][n];
    				if(nuc[i] == dawg->arrayNodes[j]->word )
    					value+=1;
    				else
    					value+=-3;
    				//value = max(value , I[depth-1][n] );
    				//value = max(value , D[depth-1][n]);

    				if(value >= maxim){
    					pos = n;
    					maxim = value;
    				}
    			}
    			//printf("%d nodo: %d\n",maxim,pos );
				if(I[depth-1][j] > 0 && G[depth-1][j]>0)
					I[depth][j]= max(I[depth-1][j],G[depth-1][j]-5) -2;
				else if(I[depth-1][j] > 0)
					I[depth][j]= I[depth-1][j]-2;
				else if(G[depth-1][j]>0)
					I[depth][j]= G[depth-1][j]-2 -5;
				else
					I[depth][j]= -999;

				if(D[depth][pos] > 0 && G[depth][pos]>0)
					D[depth][j]= max(D[depth][pos],G[depth][pos]-5) -2;
				else if(D[depth][pos] > 0)
					D[depth][j]= D[depth][pos]-2;
				else if(G[depth][pos]>0)
					D[depth][j]= G[depth][pos]-2 -5;
				else
					D[depth][j]= -999;
				maxim = max(maxim,I[depth][j]);

				maxim = max(maxim,D[depth][j]);
				if(maxim>*maxV){
					*maxV=maxim;
					sa[0] = spAux;
					sa[1] = epAux;
				}
				G[depth][j] = maxim;
				//G[depth][j] = (G[depth-1][pos] >= 0) ? maxim : -999;
				//printf("G: %d   I: %d     D:%d   \n",G[depth][j],I[depth][j],D[depth][j]);

				if(maxim>0){
					flag = 1;
				}  	/*
				if(maxim > 25){
					for(l=spAux+1 ; l<=epAux;l++){
    					printf("Posicion:%lld Valor:%d\n",getPosition(*(bwtSW->fmi),l ),maxim);
    				}
				}	*/	  			
    		}
    		if(flag ==1){
    			r = createRowDP(dawg->length,1);
    			G[depth+1] = r->G;
				I[depth+1] = r->I;
				D[depth+1] = r->D;
				free(r);
				alignReadBwtRec(bwtSW, dawg, spAux, epAux, I, D, G,depth+1,sa,maxV);
    		}
    	}
	}
	free(I[depth]);
    free(D[depth]);
    free(G[depth]);


}


void alignReadBwtSW(BwtSW * bwtSW, Read * read){
	Dawg * dawg;
	int **I,**D,**G;
	int i;
	int max= -999;
	int length;
	long long sp= 0;
	long long ep =0;
	long long l;
	dawg = createDawg(read->sequence);
	length = strlen(read->sequence);
	long long * sa = (long long *) malloc(sizeof(long long)*2);
	I = (int **)malloc(sizeof(int*)*2*length);
	D = (int **)malloc(sizeof(int*)*2*length);
	G = (int **)malloc(sizeof(int*)*2*length);
	
	I[0]=(int *)malloc(sizeof(int)* dawg->length);
	D[0]=(int *)malloc(sizeof(int)* dawg->length);
	G[0]=(int *)malloc(sizeof(int)* dawg->length);

	

	for(i = 0; i< dawg->length;i++){

		I[0][i] = -999;
		G[0][i] = 0;
		D[0][i] = -999;
	}
	RowDP *r = createRowDP(dawg->length,1);
	G[1] = r->G;
	I[1] = r->I;
	D[1] = r->D;
	free(r);
	alignReadBwtRec(bwtSW, dawg, sp, ep, I, D, G,1,sa,&max);
	for(l=sa[0]+1 ; l<=sa[1];l++){
		Alignment *a=createAlignment(getPosition(*(bwtSW->fmi),l )+1,max,"-",read->id);
		printAlignment(a);
		destroyAlignment(a);
    }

    free(I[0]);
    free(D[0]);
    free(G[0]);
    free(I);
    free(D);
    free(G);
    free(sa);
    destroyDawg(dawg);
    fflush(stdout);

}
