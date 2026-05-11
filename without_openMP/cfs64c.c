/**************************************************************************************
*
* CdL Magistrale in Ingegneria Informatica
* Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2020/21
* 
* Progetto dell'algoritmo Attention Mechanism 221 231 a
* in linguaggio assembly x86-64 + SSE
* 
* Fabrizio Angiulli, aprile 2019
* 
**************************************************************************************/

/*
* 
* Software necessario per l'esecuzione:
* 
*    NASM (www.nasm.us)
*    GCC (gcc.gnu.org)
* 
* entrambi sono disponibili come pacchetti software 
* installabili mediante il packaging tool del sistema 
* operativo; per esempio, su Ubuntu, mediante i comandi:
* 
*    sudo apt-get install nasm
*    sudo apt-get install gcc
* 
* potrebbe essere necessario installare le seguenti librerie:
* 
*    sudo apt-get install lib64gcc-4.8-dev (o altra versione)
*    sudo apt-get install libc6-dev-i386
* 
* Per generare il file eseguibile:
* 
* nasm -f elf64 fss64.nasm && gcc -m64 -msse -O0 -no-pie sseutils64.o fss64.o fss64c.c -o fss64c -lm && ./fss64c $pars
* 
* oppure
* 
* ./runfss64
* 
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>

#define	type		double
#define	MATRIX		type*
#define	VECTOR		type*

typedef struct {
	MATRIX ds; 		// dataset
	VECTOR labels; 	// etichette
	int* out;		// vettore contenente risultato dim=k
	type sc;		// score dell'insieme di features risultato
	int k;			// numero di features da estrarre
	int N;			// numero di righe del dataset
	int d;			// numero di colonne/feature del dataset
	int display;
	int silent;
} params;

typedef struct {
	type x;			//somma semplice
	type x_quad;	//somma dei quadrati
	int flag;
	type* xy;		//verrote delle somme che coinvolgono la colonna x
	int* flag_valori;
} somme_col;

/*
* 
*	Le funzioni sono state scritte assumento che le matrici siano memorizzate 
* 	mediante un array (float*), in modo da occupare un unico blocco
* 	di memoria, ma a scelta del candidato possono essere 
* 	memorizzate mediante array di array (float**).
* 
* 	In entrambi i casi il candidato dovr� inoltre scegliere se memorizzare le
* 	matrici per righe (row-major order) o per colonne (column major-order).
*
* 	L'assunzione corrente � che le matrici siano in row-major order.
* 
*/

void* get_block(int size, int elements) { 
	return _mm_malloc(elements*size,32); 
}

void free_block(void* p) { 
	_mm_free(p);
}

MATRIX alloc_matrix(int rows, int cols) {
	return (MATRIX) get_block(sizeof(type),rows*cols);
}

int* alloc_int_matrix(int rows, int cols) {
	return (int*) get_block(sizeof(int),rows*cols);
}

void dealloc_matrix(void* mat) {
	free_block(mat);
}

/*
* 
* 	load_data
* 	=========
* 
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
* 
* 	Codifica del file:
* 	primi 4 byte: numero di colonne (N) --> numero intero
* 	successivi 4 byte: numero di righe (M) --> numero intero
* 	successivi N*M*8 byte: matrix data in row-major order --> numeri floating-point a precisione doppia
* 
*****************************************************************************
*	Se lo si ritiene opportuno, � possibile cambiare la codifica in memoria
* 	della matrice. 
*****************************************************************************
* 
*/
MATRIX load_data(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols, status, i;
	
	fp = fopen(filename, "rb");
	
	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}
	
	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);
	
	MATRIX data = alloc_matrix(rows,cols);
	status = fread(data, sizeof(type), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*k = cols;
	
	if (cols != 1) {
		MATRIX data2 = alloc_matrix(rows, cols);
		for (int j = 0; j < cols; j++) {
			for (int i = 0; i < rows; i++) {
				data2[(j * rows) + i] = data[j + (i * cols)];
			}
		}
		return data2;
	}

	return data;
}

/*
* 	save_data
* 	=========
* 
*	Salva su file un array lineare in row-major order
*	come matrice di N righe e M colonne
* 
* 	Codifica del file:
* 	primi 4 byte: numero di colonne (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di righe (M) --> numero intero a 32 bit
* 	successivi N*M*8 byte: matrix data in row-major order --> numeri interi o floating-point a precisione doppia
*/
void save_data(char* filename, void* X, int n, int k) {
	FILE* fp;
	int i;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&k, 4, 1, fp);
		fwrite(&n, 4, 1, fp);
		for (i = 0; i < n; i++) {
			fwrite(X, sizeof(type), k, fp);
			//printf("%i %i\n", ((int*)X)[0], ((int*)X)[1]);
			X += sizeof(type)*k;
		}
	}
	else{
		int x = 0;
		fwrite(&x, 4, 1, fp);
		fwrite(&x, 4, 1, fp);
	}
	fclose(fp);
}

/*
* 	save_out
* 	=========
* 
*	Salva su file un array lineare composto da k+1 elementi.
* 
* 	Codifica del file:
* 	primi 4 byte: contenenti il numero di elementi (k+1)		--> numero intero a 32 bit
* 	successivi 4 byte: numero di righe (1) 						--> numero intero a 32 bit
* 	successivi byte: elementi del vettore 		--> 1 numero floating-point a precisione doppia e k interi
*/
void save_out(char* filename, type sc, int* X, int k) {
	FILE* fp;
	int i;
	int n = 1;
	k++;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&n, 4, 1, fp);
		fwrite(&k, 4, 1, fp);
		fwrite(&sc, sizeof(type), 1, fp);
		fwrite(X, sizeof(int), k, fp);
		//printf("%i %i\n", ((int*)X)[0], ((int*)X)[1]);
	}
	fclose(fp);
}

// PROCEDURE ASSEMBLY

// PROCEDURE ASSEMBLY
extern void prova(params *input);
//extern void calcoloSommatoriaQuadrato(params *input, type media, int k, type *risultato);
//extern void calcoloSommatoria(params *input, type mediaColonna1, type mediaColonna2, int indiceColonna1, int indiceColonna2, type *risultato);
//extern void calcoloMediaAssembly(params *input, int j, type *risultato);
//extern void calcoloValori_n0n1(params* input, unsigned long* n0, unsigned long* n1);
extern void calcoloMediaGruppoAssembly(params *input, int j, int n0, int n1, type* differenzaMedia);
//extern void calcoloPearsonCoefficient2(valori *valore, type* sommaX, type* sommaY, type* sommaXY, type* sommaX2, type* sommaY2);
extern void calcoloPearsonCoefficient2(params *input, int indiceColonna1, int indiceColonna2, type* sommaX, type* sommaX2, type* divisione);
extern void calcoloDeviazioneStandard2(params *input, int j, type* somma_sx, type* somma_dx);

extern type sqrtAssembly(type valore);

//type calcoloQuadrato(type valore);
type calcoloPearsonCoefficient(params *input, int indiceColonna1, int indiceColonna2);


type calcoloScore(int k, int j, params *input);
type calcoloScoreIniziale(int indice_feature_corrente, int indice_colonna, params* input);

unsigned long n0 = 0;
unsigned long n1 = 0;
unsigned long n = 0;
type fat1 = 0.0;
type fat2 = 0.0;
type fattore2 = 0.0;

type rff_passato; //serve per evitare di fare il doppio for per il calcolo dell'rff medio, all'interno della funzione calcolo score. La modifica di questo valore deve avvenire dentro l'algoritmo di CFS perchè noi vogliamo memorizzare l'rff passato corrispondente alla colonna che ha generato lo score massimo. Dunque per questo motivo la signature della funzione calcolo score prenderà un ulteriore parametro tramite cui passare questo valore. Questa variabile verrà inizializzata a 0 nel calcolo score iniziale
type rff_passato_momentaneo; //stesso concetto di prima ma questa variabile la usiamo solo per passare il valore del rff appena calcolato dalla funzione di calcolo score alla funzione CFS perchè non possiamo subito sovrascrivere l'rff_passato, lo sovrascriviamo solo se è quello calcolato in corrispondenza dello score massimo

type rcf_assoluto;
type rcf_passato;
type* mediaGruppo;

somme_col* vettore_somme; //vettore contenente le somme delle colonne
/*
Autore: Giovanni
Ultima modifica: 02/01/2024
Legge funzione: Per ogni colonna, calcola il valore di merito dell'insieme S unito quella colonna.
                Alla fine all'insieme iniziale S aggiunge solo la colonna che ha fatto ottenere
                il valore di merito maggiore.
                In poche parole, è una funzione che calcola il massimo della funzione di merito.
                Ciò viene fatto per "k" colonne.
Funzioni Esterne usate: -calcoloScoreIniziale()
                        -calcoloScore()
*/

void inizializzazione_struttura_di_supporto(params *input)  //viene richiamata nella funzione CorrelationFeatureSelection(...)
{
	vettore_somme = (somme_col*) malloc(input->d * sizeof(somme_col));

	for(int i=0;i<input->d;i++)
	{
		vettore_somme[i].flag=0;	//con questa flag dico che la somma di x e x_quad ancora non è stata calcoalta
		vettore_somme[i].xy = (type*) calloc(sizeof(type), input->d);			//allocco il vettore di somme che coinvolgono la colonna x
		vettore_somme[i].flag_valori = (int*) calloc(sizeof(int), input->d);	//allocco il vettore di somme che coinvolgono la colonna x
		//for(int j=0;j<input->d;j++)
		//	vettore_somme[i].flag_valori[j] = 0;    //con questa flag dico che la somma tra la colonna x e la colonna j ancora non è stata calcolata
	}
}

void CorrelationFeatureSelection(params *input)
{
	inizializzazione_struttura_di_supporto(input);


	//Ignorare questa riga, ritornarci dopo aver letto l'intera funzione
    //visto che alla fine sommiamo progressivamente lo score, lo inizializziamo a 0, che è il valore invariante per la somma
    input->sc=0;
	//calcoloValori_n0n1(input, &n0, &n1);
	for (int i = 0; i < input->N; i++){
		if (input->labels[i] == 0){
			n0++;
		}else{
			n1++;
		}
	}
	type* vettore = (type*)calloc(sizeof(type), input->d);
	mediaGruppo = (type*)calloc(sizeof(type), input->d);
	
	n = n0 + n1;
	fat1 = n0 / (type) n;
	fat2 = n1 / (type) n;
	fattore2 = sqrt(fat1 * fat2);
	//Fattore2 = sqrtAssembly(fat1 * fat2);

    //non ci interessa inizializzare out a tutti "-1" perchè di default viene allocato con "k" posizioni disponibili, quindi le sovrascriveremo tutte
	for(int i=0;i<input->k;i++) {  //Eseguiamo la legge della funzione per "k" volte, quindi vogliamo selezionare "k" colonne
    
        //l'indice "i" farà riferimento al vettore "out" (vettore degli indici delle colonne selezionate)
        //indice_feature_corrente=i (quindi indica la cardinalità dell'insieme S prima di aggiungere la colonna corrente)
        //l'indice "j" farà riferimento alla colonna su cui calcoliamo lo score
        //indice_colonna=i

        //visto che in pratica stiamo calcolando una funzione di massimo, inizializziamo i valori del confronto per capire qual'è il massimo inserendo dei valori fantoccio
        int indice_max=-1;
        type score_max=-1.0;
		type rff_max=0.0;  //va inizializzato a 0.0 così da sovrascrivere l'rff_passato ed inizializzarlo alla prima iterazione
		type rcf_max=0.0;

        //in questo for calcoliamo il massimo merito vero e proprio
        for(int j=0;j<input->d;j++) {
            int flag=1; //la flag viene impostata a true di default perchè si presuppone che il "j" corrente non sia una colonna già inserita. Se notiamo che è già stata inserita allora la impostiamo a false (0)
          	if(vettore[j] == 1)
                flag = 0;
            if(i==0) { //se l'insieme S è vuoto, non possiamo calcolare la correlazione delle colonne fra di loro, quindi calcoliamo solo la correlazione tra la singola colonna e la variabile dicotomica "c"
            
                type scoreCalcolato=calcoloScoreIniziale(i,j,input);

                if(score_max<scoreCalcolato) {
                    score_max=scoreCalcolato;
                    indice_max=j;
					rcf_max=rcf_assoluto;
                }
            } else {   //se l'insieme S contiene almeno una colonna, calcoliamo lo score prendendo anche in considerazione il coefficiente di correlazione di pearson
                //if la colonna corrente non appartiene già all'insieme S
                if(flag) { //questo controllo nell'if(i==0) non lo facciamo perchè ridondante
                    type scoreCalcolato=calcoloScore(i,j,input); //passiamo per riferimento questa variabile perchè ci serve solo per farci restituire un valore
                    if(score_max<scoreCalcolato) {
                        score_max=scoreCalcolato;
                        indice_max=j;
						rff_max=rff_passato_momentaneo;
						rcf_max=rcf_assoluto;
                    }
                }
            }//if else
        }//for su j

        //Alla fine del for, aggiungiamo la colonna con lo score massimo all'insieme S; ovvero, aggiungiamo il suo indice al vettore out
        input->out[i]=indice_max;
		vettore[indice_max] = 1;
        //inoltre aggiorniamo lo score, in quanto lo calcoliamo in maniera progressiva
        input->sc=score_max;  
		//aggironiamo l'rff passato con quello calcolato in corrispondenza dello score massimo
		rff_passato=rff_max;

		rcf_passato=rcf_max;
    }//for su i
}

/*
Autore: Giovanni
Ultima modifica: 31/12/2023
Legge funzione: Effettua il calcolo del coefficiente di correlazione di Pearson
                tra due colonne della matrice "ds"
Funzioni Esterne usate: -calcoloSommatoria()
                        -calcoloSommatoriaQuadrato()
                        -calcoloMedia()
*/
type calcoloPearsonCoefficient(params *input, int indiceColonna1, int indiceColonna2)
{
	type somma_X = 0.0;
	//type somma_Y = 0.0;
	type somma_XY = 0.0;
	type somma_X2 = 0.0;
	//type somma_Y2 = 0.0;
	type divisione = 0.0;

	if(vettore_somme[indiceColonna1].flag_valori[indiceColonna2]!=0 && vettore_somme[indiceColonna1].flag!=0 && vettore_somme[indiceColonna2].flag!=0)
	{
		//printf("colonna[1]: %d\ncolonna[2]: %d",indiceColonna1, indiceColonna2);
		somma_X = vettore_somme[indiceColonna1].x;
		somma_X2 = vettore_somme[indiceColonna1].x_quad;
		//somma_Y = vettore_somme[indiceColonna2].x;
		//somma_Y2 = vettore_somme[indiceColonna2].x_quad;
		divisione = vettore_somme[indiceColonna1].xy[indiceColonna2];
	}else{
		calcoloPearsonCoefficient2(input, indiceColonna1, indiceColonna2, &somma_X, &somma_X2, &divisione);
		//aggiorniamo la struttura dati

		vettore_somme[indiceColonna1].x=somma_X;
		vettore_somme[indiceColonna1].x_quad=somma_X2;
		//vettore_somme[indiceColonna2].x=somma_Y;
		//vettore_somme[indiceColonna2].x_quad=somma_Y2;
		//la somma di coppia va aggiornata su tutti e due
		vettore_somme[indiceColonna1].xy[indiceColonna2]=divisione;
		vettore_somme[indiceColonna2].xy[indiceColonna1]=divisione;
		//aggironiamo le flag
		vettore_somme[indiceColonna1].flag_valori[indiceColonna2]=1;
		vettore_somme[indiceColonna2].flag_valori[indiceColonna1]=1;
		vettore_somme[indiceColonna1].flag=1;
		vettore_somme[indiceColonna2].flag=1;
	}
	//type numeratore = input->N*somma_XY - (somma_X * somma_Y);
	//type denominatore = sqrt(((input->N*somma_X2) - (somma_X*somma_X)) *  ((input->N*somma_Y2) - (somma_Y*somma_Y)));
	return divisione;

}

type pointBiserialCorrelationCoefficient(params *input, int j)
{
	
	type differenzaMedia=0.0;
	
	if(mediaGruppo[j] == 0){
		
		calcoloMediaGruppoAssembly(input, j, n0, n1, &differenzaMedia);
		mediaGruppo[j] = differenzaMedia;
	} else {
		differenzaMedia = mediaGruppo[j];
	}
	

	//COLCOLO DEVIAZIONE STANDARD
	type somma_sx = 0.0;
	type somma_dx = 0.0;
	//type denominatore = 0.0;
	if(vettore_somme[j].flag!=0)
	{
		somma_dx=vettore_somme[j].x;
		somma_sx=vettore_somme[j].x_quad;
		//denominatore = vettore_somme[j].x;
	}else{
		calcoloDeviazioneStandard2(input, j, &somma_sx, &somma_dx);
		//calcoloDeviazioneStandard2(input, j, &denominatore);
		//aggiorniamo la struttura dati
		vettore_somme[j].x=somma_dx;
		vettore_somme[j].x_quad=somma_sx;
		//la somma di coppia non possiamo aggiornarla
		//aggironiamo le flag
		vettore_somme[j].flag=1;
	}
	type denominatore = sqrt((somma_sx - (somma_dx * somma_dx/n))/(n-1));

	type fattore1 =  differenzaMedia / denominatore;

	return fattore1 * fattore2;
}

//funzione merit
type calcoloScore(int indice_feature_corrente, int indice_colonna, params* input)
{
    //aggiorniamo l'insieme di colonne S
    int k=indice_feature_corrente+1;  //+1 perchè la cardinalità dell'insieme S è pari a k, ma dopo che aggiungiamo la colonna corrente per calcolare il nuovo score diventa k+1
    input->out[indice_feature_corrente]=indice_colonna;  //S unione f_i (inseriamo la colonna corrente)
    
    //calcoliamo prima il numeratore
    //calcoliamo il point biserial correlation coefficient

    //type rcf_assoluto=0.0;
    //for(int i=0;i<k;i++) {  //calcoliamo la media dei valori assoluti
    //    rcf_assoluto+=valoreAssoluto(pointBiserialCorrelationCoefficient(input, input->out[i]));
	//}

	type numero = pointBiserialCorrelationCoefficient(input, indice_colonna);
	if (numero < 0){
		rcf_assoluto = rcf_passato + (numero * -1);
	} else {
		rcf_assoluto = rcf_passato + numero ;
	}


    type numeratore=rcf_assoluto;  //semplifichiamo il fatto che la media sia divisa per k, e si moltiplica per k
    
    //calcoliamo ora il denominatore
    //calcoliamo il coefficiente di correlazione di pearson tra la colonna corrente e le k-1 colonne già presenti
    type rff_assoluto=0.0;  
	for(int i=0;i<k-1;i++) {//calcoliamo la media dei valori assoluti (a sto giro i si ferma prima perchè non vogliamo calcolare il coefficiente di pearson tra la colonna appena inserita e se stessa, ovviamente la colonna appena inserita è in coda)
		//printf("IC1: %d - ID2: %d \n",input->out[i], input->out[indice_feature_corrente]);
		type numero = calcoloPearsonCoefficient(input, input->out[i], input->out[indice_feature_corrente]);
		if (numero < 0){
			rff_assoluto += numero * -1;
		} else {
			rff_assoluto += numero;
		}
	}
	//invece di calcolare l'rff tra tutte le coppie di colonne selezionate, lo calcoliamo tra le k-1 già presenti e la k-esima appena inserita, e poi ne sommiamo l'rff passato che è una variabile globale che tiene traccia del calcolo tra tutte le possibili coppie che ci manca
	rff_assoluto+=rff_passato;


	//aggiorniamo l'rff passato per la prossima volta che questa funzione verrà richiamata (in realtà poi deciderà la funzione chiamante se renderlo permanente o meno in base al fatto se lo score calcolato adesso è il max tra tutti)
	rff_passato_momentaneo=rff_assoluto;

	//ora serve calcolare la media, dunque ci serve sapere il numero per cui dividere, tale numero lo calcoliamo tramite la serie triangolare (esplicitata in funzione di k) (chiedere a Giovanni come funziona)
	type den=((k-1)*k)/2;
    type denominatore=(rff_assoluto/den)*(k)*(k-1); //stessa semplificazione di (k-1)/(k-1) di prima

	denominatore+=k;

    //eliminiamo la colonna dall'insieme
    input->out[indice_feature_corrente]=-1;

    //assembliamo
    return numeratore/sqrt(denominatore);
}

/*
Autore: Giovanni
Ultima modifica: 02/01/2024
Legge funzione: Implementa la funzione di calcolo del merito dell'insieme di colonne selezionate.
                Questa versione della funzione serve solo quando la cardinalità dell'insieme S è pari a 0.
                Serve solo per calcolare la correlazione della singola colonna rispetto alla variabile
                dicotomica "c".
Funzioni Esterne usate: -valoreAssoluto()
                        -pointBiserialCorrelationCoefficient()
*/
type calcoloScoreIniziale(int indice_feature_corrente, int indice_colonna, params* input)
{
    //aggiorniamo l'insieme di colonne S
    int k=indice_feature_corrente+1;  //+1 perchè la cardinalità dell'insieme S è pari a k, ma dopo che aggiungiamo la colonna corrente per calcolare il nuovo score diventa k+1
    input->out[indice_feature_corrente]=indice_colonna;  //S unione f_i (inseriamo la colonna corrente)
    
    //calcoliamo prima il numeratore
    //calcoliamo il point biserial correlation coefficient
    //type rcf_assoluto=0.0;
    //rcf_assoluto+=valoreAssoluto(pointBiserialCorrelationCoefficient(input, indice_colonna));

	type numero = pointBiserialCorrelationCoefficient(input, indice_colonna);

	if (numero < 0){
		rcf_assoluto = -1.0 * numero;
	}
	

    type numeratore=rcf_assoluto;  //semplifichiamo il fatto che la media sia divisa per k, e si moltiplica per k
    
    //calcoliamo ora il denominatore
    //calcoliamo il coefficiente di correlazione di pearson ad 1, perchè abbiamo una sola colonna nell'insieme S
    
    //l'unica differenza con l'altra funzione è in queste 2 righe di codice
    type rff_assoluto=1.0;     
    type denominatore=rff_assoluto*(k)*(k-1); //non si semplifica il (k-1) perchè non c'è più il concetto di media, in quando il coefficiente di pearson è proprio 1
    
    denominatore+=k;

    //eliminiamo la colonna dall'insieme
    input->out[indice_feature_corrente]=-1;
	
	//inzializziamo l'rff del passato
	rff_passato=0.0;
	
    //assembliamo
    return numeratore/sqrt(denominatore);
}


void cfs(params *input)
{
	if(input->k > input->d){
		printf("Attenzione, numero di featureas da indovinare maggiore del numero di features totali. Riprovare. \n");
		return;
	}
	// Scrittura scheletro algorirmo 2
	CorrelationFeatureSelection(input);
}

int main(int argc, char** argv) {

	char fname[256];
	char* dsfilename = NULL;
	char* labelsfilename = NULL;
	clock_t t;
	float time;
	
	//
	// Imposta i valori di default dei parametri
	//

	params* input = malloc(sizeof(params));

	input->ds = NULL;
	input->labels = NULL;
	input->k = -1;
	input->sc = -1;

	input->silent = 0;
	input->display = 0;

	//printf("%i\n", sizeof(double));

	//
	// Visualizza la sintassi del passaggio dei parametri da riga comandi
	//

	if(argc <= 1){
		printf("%s -ds <DS> -labels <LABELS> -k <K> [-s] [-d]\n", argv[0]);
		printf("\nParameters:\n");
		printf("\tDS: il nome del file ds2 contenente il dataset\n");
		printf("\tLABELS: il nome del file ds2 contenente le etichette\n");
		printf("\tk: numero di features da estrarre\n");
		printf("\nOptions:\n");
		printf("\t-s: modo silenzioso, nessuna stampa, default 0 - false\n");
		printf("\t-d: stampa a video i risultati, default 0 - false\n");
		exit(0);
	}

	//
	// Legge i valori dei parametri da riga comandi
	//

	int par = 1;
	while (par < argc) {
		if (strcmp(argv[par],"-s") == 0) {
			input->silent = 1;
			par++;
		} else if (strcmp(argv[par],"-d") == 0) {
			input->display = 1;
			par++;
		} else if (strcmp(argv[par],"-ds") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing dataset file name!\n");
				exit(1);
			}
			dsfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-labels") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing labels file name!\n");
				exit(1);
			}
			labelsfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-k") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing k value!\n");
				exit(1);
			}
			input->k = atoi(argv[par]);
			par++;
		} else{
			printf("WARNING: unrecognized parameter '%s'!\n",argv[par]);
			par++;
		}
	}

	//
	// Legge i dati e verifica la correttezza dei parametri
	//

	if(dsfilename == NULL || strlen(dsfilename) == 0){
		printf("Missing ds file name!\n");
		exit(1);
	}

	if(labelsfilename == NULL || strlen(labelsfilename) == 0){
		printf("Missing labels file name!\n");
		exit(1);
	}


	input->ds = load_data(dsfilename, &input->N, &input->d);

	int nl, dl;
	input->labels = load_data(labelsfilename, &nl, &dl);
	
	if(nl != input->N || dl != 1){
		printf("Invalid size of labels file, should be %ix1!\n", input->N);
		exit(1);
	} 

	if(input->k <= 0){
		printf("Invalid value of k parameter!\n");
		exit(1);
	}

	input->out = alloc_int_matrix(input->k, 1);

	//
	// Visualizza il valore dei parametri
	//

	if(!input->silent){
		printf("Dataset file name: '%s'\n", dsfilename);
		printf("Labels file name: '%s'\n", labelsfilename);
		printf("Dataset row number: %d\n", input->N);
		printf("Dataset column number: %d\n", input->d);
		printf("Number of features to extract: %d\n", input->k);
	}

	// COMMENTARE QUESTA RIGA!
	//prova(input);
	//

	//
	// Correlation Features Selection
	//

	t = clock();
	cfs(input);
	t = clock() - t;
	time = ((float)t)/CLOCKS_PER_SEC;

	if(!input->silent)
		printf("CFS time = %.3f secs\n", time);
	else
		printf("%.3f\n", time);

	//
	// Salva il risultato
	//
	sprintf(fname, "out64_%d_%d_%d.ds2", input->N, input->d, input->k);
	save_out(fname, input->sc, input->out, input->k);
	if(input->display){
		if(input->out == NULL)
			printf("out: NULL\n");
		else{
			int i,j;
			printf("sc: %lf, out: [", input->sc);
			for(i=0; i<input->k; i++){
				printf("%i,", input->out[i]);
			}
			printf("]\n");
		}
	}

	if(!input->silent)
		printf("\nDone.\n");

	dealloc_matrix(input->ds);
	dealloc_matrix(input->labels);
	dealloc_matrix(input->out);
	free(input);
	
	return 0;
}
