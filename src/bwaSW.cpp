
#include <iostream>

#include <math.h>
#include <algorithm>
#include "bwaSW.h"
#include "dawg.h"
#include <algorithm>
#include <stack>          // std::stack


using namespace std;
/*
int s;

BwaSW * createBwaSW(char * file){

	BwaSW * bwaSW = (BwaSW *) malloc(sizeof (BwaSW));
    fflush(stdout);
    string genomeString;
    ifstream f;
    f.open(file);
    string line;
    getline(f, line);
    fflush(stdout);
    bwaSW->nNuc = std::stoll(line.c_str(), NULL, 0);
    fflush(stdout);
    while (getline(f, line)) {
        genomeString.append(line);

    }
    genomeString.append("$");

    bwaSW->fmi = createFMI(genomeString.c_str(), bwaSW->nNuc);
    
    f.close();
    return bwaSW;

}

void removeBwaSW(BwaSW * bwaSW) {
    removeFMI(bwaSW->fmi);
    free(bwaSW);
}

void saveBwaSW(char * file , BwaSW *bwaSW){
    FILE * f;
    f = fopen(file,"w");
    saveFMI(f,bwaSW->fmi);
    fclose(f);
}

BwaSW * loadBwaSW(char * file){
    FILE * f;
    f = fopen(file,"r");
    BwaSW * bwaSW = (BwaSW *)malloc(sizeof(BwaSW));
	
    bwaSW->fmi = loadFMI(f);
    bwaSW->nNuc = bwaSW->fmi->nNuc;
    fclose(f);
    return bwaSW;
}

ColumnDP * createComlumnDP(int length){
	int i;
	ColumnDP * c = (ColumnDP *) malloc(sizeof(ColumnDP));
	if(!c)
		printf("ERROOOR\n");
	c->childs = (ColumnDP **) malloc(sizeof(ColumnDP*)*4);
	if(!c->childs)
		printf("ERROOOR\n");
	for(i = 0; i<4;i++){
		c->childs[i] = NULL;
	}
	c->G = (int *) malloc(sizeof(int)*length);
	if(!c->G)
		printf("ERROOOR\n");
	c->D = (int *) malloc(sizeof(int)*length);
	if(!c->D)
		printf("ERROOOR\n");
	c->I = (int *) malloc(sizeof(int)*length);
	if(!c->I)
		printf("ERROOOR\n");
	for(i = 0; i<1 ;i++){
		c->G[i] = 0;
		c->D[i] = 0;
		c->I[i] = 0;
	}
	return c;
}

void removeColumnDP(ColumnDP * c){
	printf("Eliminando\n");
	if(c){
		free(c->childs);
		free(c->G);
		free(c->D);
		free(c->I);
		free(c);
	}
}
void removeColumnDPTree(ColumnDP *c){
	int i;
	for(i = 0; i< 4; i++){
		if(c->childs[i]!=NULL) {
			printf("holis\n");
			removeColumnDPTree(c->childs[i]);
		}
	}
	removeColumnDP(c);
}

void alignReadBwaSWAuxv2(BwaSW * bwaSW, Dawg * dawg, long long sp, long long ep, ColumnDP * father ,ColumnDP *dp, int deep, char c, ColumnDP **maxNode){
	int k,n,i;
	int value,pos,maxim;
	long long spAux, epAux;
	char nuc[]="ACGT";
	s++;
	printf("Bajo: %c\n", c);

	//value = (dawg->arrayNodes[deep]->word == c) ? 1 : -3;
	maxim = -99999;
	for(k = 0; k < dawg->arrayNodes[deep]->NFathers; k++){
		n = dawg->arrayNodes[deep]->fathers[k];
		printf("Padre: %d\n", n);
		value = (dawg->arrayNodes[deep]->word == c) ? 1 : -3;
		value = max(father->G[n]+value, father->I[n]);
		value = max(value,father->D[n]);
		printf("%d\n",value );
		if(value >= maxim){
			pos = n;
			maxim = value;
		}
	}
	if(maxim > 0){
		dp->I[deep] = max(dp->I[pos] , dp->G[pos] -5) -2;
		dp->D[deep] = max(father->D[deep] , father->G[deep] -5) -2;
		dp->G[deep] = maxim;
	}
	else{
		dp->I[deep] = -99999;
		dp->D[deep] = -99999;
		dp->G[deep] = -99999;
	}
	printf("%d\n", maxim);
	if(maxim>11){
		for(i=sp+1 ; i<=ep;i++){
			printf("Posicion:%lld  %d\n",getPosition(*(bwaSW->fmi),i),maxim);
		}
		
	}
	if(dp->G[deep] >= 0){
		for(k = 0 ; k < 4; k++){
			spAux = getMinRank(nuc[k], bwaSW->fmi->ranking) + getCheckpoints(bwaSW->fmi, sp, nuc[k]) -1;
        	epAux = getMinRank(nuc[k], bwaSW->fmi->ranking) + getCheckpoints(bwaSW->fmi, ep, nuc[k]) -1;
        	if(spAux < epAux){
        		if(dp->childs[k] == NULL){
        			dp->childs[k] = createComlumnDP(dawg->length);
        		}
        		//printf("Nodo: %d %d\n",dp->G[1],dp );
        		alignReadBwaSWAuxv2(bwaSW,dawg,spAux,epAux,dp,dp->childs[k],deep,nuc[k],maxNode);
        	}

		}
	}
	printf("Subo\n");
}

void alignReadBwaSWv2(BwaSW * bwaSW, Read * read){
	Dawg * dawg;
	int i,j;
	int length;
	long long sp, ep;
	ColumnDP * dp;
	ColumnDP ** pointer;
	printf("Entro\n");
	dawg = createDawg(read->sequence);
	printf("%s\n",read->sequence );
	length = dawg->length;
	printf("Dawg\n");
	s=0;
	ColumnDP * root = createComlumnDP(length);

	root->D[0] = 0;
	root->G[0] = 0;
	root->I[0] = 0;
	for(i = 1; i<length ;i++){
		root->G[i] = 0;
		root->D[i] = 0;
		root->I[i] = 0;
	}
	for(i = 1; i<length; i++){
		printf("Nivel %d\n",i );
		char nuc[]="ACGT";
		for(j = 0; j < 4;j++){
			sp = getMinRank(nuc[j], bwaSW->fmi->ranking) - 1;
	    	ep = getMaxRank(nuc[j], bwaSW->fmi->ranking) - 1;
	    	if(sp < ep){
	    		if(root->childs[j] == NULL){
	    			dp = createComlumnDP(length);
	    			root->childs[j] = dp;
	    		}
	    		alignReadBwaSWAuxv2(bwaSW,dawg,sp,ep,root,root->childs[j],i,nuc[j],pointer);

	    	}
	    	//printf("%lld--%lld\n", sp,ep);

		}
	}
	printf("antes de liberar\n");
	//removeColumnDPTree(root);
    destroyDawg(dawg);
    printf("liberado\n");

}


void alignReadBwaSWAux(BwaSW * bwaSW, Dawg * dawg, long long sp, long long ep, int ** I, int ** D, int **G,int deep){
	int value,pos,maxim;
	long long spAux, epAux;
	int i,n,j,k;
	long long l;

	I[deep]=(int *)malloc(sizeof(int)*dawg->length);
	D[deep]=(int *)malloc(sizeof(int)*dawg->length);
	G[deep]=(int *)malloc(sizeof(int)*dawg->length);
	I[deep][0] = -99999;
	D[deep][0] = -99999;
	G[deep][0] = -99999;

	char nuc[]="ACGT";
	for(i = 0; i < 4;i++){
		spAux = getMinRank(nuc[i], bwaSW->fmi->ranking) + getCheckpoints(bwaSW->fmi, sp, nuc[i]) -1;
        epAux = getMinRank(nuc[i], bwaSW->fmi->ranking) + getCheckpoints(bwaSW->fmi, ep, nuc[i]) -1;
        if(spAux < epAux){
    		for(j = 1; j < dawg->length ; j++){
    			maxim = -99999;
    			for(k = 0; k < dawg->arrayNodes[j]->NFathers; k++){
    				n = dawg->arrayNodes[j]->fathers[k];
    				value = G[deep-1][n];
    				if(nuc[i] == dawg->arrayNodes[j]->word )
    					value+=1;
    				else
    					value+=-3;
    				value = max(value , I[deep-1][n] );
    				value = max(value , D[deep-1][n]);

    				if(value >= maxim){
    					pos = n;
    					maxim = value;
    				}
    			}
    			if(maxim > 0){
    				I[deep][j] = max(I[deep][pos],G[deep][pos] -5) -2;
					D[deep][j] = max(D[deep-1][j],G[deep-1][j]-5) -2;
					G[deep][j] = maxim;
					//flag =
    			}
    			else{
    				I[deep][j] = -99999;
					D[deep][j] = -99999;
					G[deep][j] = -99999;
    			}
    			if(maxim>20){
    				for(l=sp+1 ; l<=ep;l++){
    					printf("Posicion:%lld\n",getPosition(*(bwaSW->fmi),l ));
    				}
    				
    			}
				
    		}
    		alignReadBwaSWAux(bwaSW, dawg, spAux, epAux, I, D, G,deep +1);	
    	}
	}
	free(I[deep]);
    free(D[deep]);
    free(G[deep]);


}

void alignReadBwaSW(BwaSW * bwaSW, Read * read){
	Dawg * dawg;
	int **I,**D,**G;
	int i,n,j,k;
	int length;
	long long sp, ep;
	int value,maxim;
	dawg = createDawg(read->sequence);
	length = strlen(read->sequence);

	I = (int **)malloc(sizeof(int*)*4*length);
	D = (int **)malloc(sizeof(int*)*4*length);
	G = (int **)malloc(sizeof(int*)*4*length);
	
	I[0]=(int *)malloc(sizeof(int)* dawg->length);
	D[0]=(int *)malloc(sizeof(int)* dawg->length);
	G[0]=(int *)malloc(sizeof(int)* dawg->length);
	I[1]=(int *)malloc(sizeof(int)* dawg->length);
	D[1]=(int *)malloc(sizeof(int)* dawg->length);
	G[1]=(int *)malloc(sizeof(int)* dawg->length);

	I[0][0] = 0;
	D[0][0] = 0;
	G[0][0] = 0;
	I[1][0] = -99999;
	D[1][0] = -99999;
	G[1][0] = -99999;

	for(i = 1; i< dawg->length;i++){

		I[0][i] = -99999;
		D[0][i] = -99999;
		G[0][i] = -99999;
	}
	char nuc[]="ACGT";
	for(i = 0; i < 4;i++){
		sp = getMinRank(nuc[i], bwaSW->fmi->ranking) - 1;
    	ep = getMaxRank(nuc[i], bwaSW->fmi->ranking) - 1;
    	if(sp < ep){
    		for(j = 1; j < dawg->length ; j++){
    			maxim = -99999;
    			for(k = 0; k < dawg->arrayNodes[j]->NFathers; k++){
    				n = dawg->arrayNodes[j]->fathers[k];
    				value = G[0][n];
    				if(nuc[i] == dawg->arrayNodes[j]->word )
    					value+=1;
    				else
    					value+=-3;
    				
    				value = max(value , I[0][n]);
    				value = max(value , D[0][n]);


    				if(value >= maxim){
    					maxim = value;
    				}
    			}
    			if(maxim > 0){
    				I[1][j] = -99999;
					D[1][j] = -99999;
					G[1][j] = maxim;
    			}
    			else{
    				I[1][j] = -99999;
					D[1][j] = -99999;
					G[1][j] = -99999;
    			}
				
    		}
    		alignReadBwaSWAux(bwaSW, dawg, sp, ep, I, D, G,2);
    		
    	}

	}
	free(I[1]);
    free(D[1]);
    free(G[1]);
    free(I[0]);
    free(D[0]);
    free(G[0]);
    free(I);
    free(D);
    free(G);
    destroyDawg(dawg);

}
*/