; ---------------------------------------------------------
; Regressione con istruzioni SSE a 32 bit
; ---------------------------------------------------------
; F. Angiulli
; 23/11/2017
;

;
; Software necessario per l'esecuzione:
;
;     NASM (www.nasm.us)
;     GCC (gcc.gnu.org)
;
; entrambi sono disponibili come pacchetti software 
; installabili mediante il packaging tool del sistema 
; operativo; per esempio, su Ubuntu, mediante i comandi:
;
;     sudo apt-get install nasm
;     sudo apt-get install gcc
;
; potrebbe essere necessario installare le seguenti librerie:
;
;     sudo apt-get install lib32gcc-4.8-dev (o altra versione)
;     sudo apt-get install libc6-dev-i386
;
; Per generare file oggetto:
;
;     nasm -f elf32 fss32.nasm 
;
%include "sseutils32.nasm"

section .data			; Sezione contenente dati inizializzati
	input equ 8 ; puntatore a float, occupa 32 bit (4 bytes)
	mediaValoriColonna equ 12 ; float a 32 bit
	indiceColonna equ 16 ; intero a 32 bit
	risultato equ 20 ; intero a 32 bit
	j equ 12
	indiceColonna1 equ 20 ; intero a 32 bit
	indiceColonna2 equ 24 ; intero a 32 bit
	indiceColonna11 equ 12
	indiceColonna22 equ 16
	risultatoX equ 20
	risultatoY equ 24
	risultatoXY equ 28
	risultatoXX equ 32
	risultatoYY equ 36
	somma_sx equ 16
	somma_dx equ 20
	_one dd 1.0
section .bss			; Sezione contenente dati non inizializzati
	alignb 16
	sc		resd		1

section .text			; Sezione contenente il codice macchina


; ----------------------------------------------------------
; macro per l'allocazione dinamica della memoria
;
;	getmem	<size>,<elements>
;
; alloca un'area di memoria di <size>*<elements> bytes
; (allineata a 16 bytes) e restituisce in EAX
; l'indirizzo del primo bytes del blocco allocato
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)
;
;	fremem	<address>
;
; dealloca l'area di memoria che ha inizio dall'indirizzo
; <address> precedentemente allocata con getmem
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)

extern get_block
extern free_block

%macro	getmem	2
	mov	eax, %1
	push	eax
	mov	eax, %2
	push	eax
	call	get_block
	add	esp, 8
%endmacro

%macro	fremem	1
	push	%1
	call	free_block
	add	esp, 4
%endmacro

; ------------------------------------------------------------
; Funzioni
; ------------------------------------------------------------


;global sqrtAssembly
;global calcoloSommatoriaQuadrato
;global calcoloMediaAssembly
;global calcoloSommatoria
global prova
global calcoloValori_n0n1
global calcoloMediaGruppoAssembly
global calcoloPearsonCoefficient2
global calcoloDeviazioneStandard2

msg	db	'sc:',32,0
nl	db	10,0

prova:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		ebp		; salva il Base Pointer
		mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
		push		ebx		; salva i registri da preservare
		push		esi
		push		edi
		; ------------------------------------------------------------
		; legge i parametri dal Record di Attivazione corrente
		; ------------------------------------------------------------

		; elaborazione
		
		; esempio: stampa input->sc
		mov EAX, [EBP+input]	; indirizzo della struttura contenente i parametri
        ; [EAX] input->ds; 			// dataset
		; [EAX + 4] input->labels; 	// etichette
		; [EAX + 8] input->out;	// vettore contenente risultato dim=(k+1)
		; [EAX + 12] input->sc;		// score dell'insieme di features risultato
		; [EAX + 16] input->k; 		// numero di features da estrarre
		; [EAX + 20] input->N;		// numero di righe del dataset
		; [EAX + 24] input->d;		// numero di colonne/feature del dataset
		; [EAX + 28] input->display;
		; [EAX + 32] input->silent;
		MOVSS XMM0, [EAX+12]
		MOVSS [sc], XMM0 
		prints msg
		printss sc     
		prints nl
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------

		pop	edi		; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp	; ripristina lo Stack Pointer
		pop	ebp		; ripristina il Base Pointer
		ret			; torna alla funzione C chiamante


;type calcoloDeviazioneStandard(params *input, int j, type* somma_sx, type* somma_dx)
calcoloDeviazioneStandard2:
	push		ebp		; salva il Base Pointer
	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push		ebx		; salva i registri da preservare
	push		esi
	push		edi

	; mi salvo i valori di utilità
	mov EAX, [EBP + input]	        ; mi prendo il parametro dell'indirizzo della struct contenente i parametri
	mov EBX, [EAX + 20]     		; mi salvo in EBX il valore di 'N'
	mov ECX, [EAX]                  ; metto in ECX il primo indirizzo del vettore ds  

	; in ECX va il seguente calcolo = indirizzo di memoria in cui comincia la colonna da valutare 
	; "indiceColonna" = indirizzoBaseVettoreDS+(indiceColonna*numeroRighe*4)
	mov EDX, [EBP + indiceColonna11] ; in EDX metto l'indice di colonna da valutare
	imul EDX, 4
	imul EDX, EBX
	add ECX, EDX

	; implementazione della sommatoria
	; siccome bisogna sommare valori al quadrato allora la strategia è quella di salvare in due registri xmm
	; gli stessi valori e poi moltiplicarli	
	xorps XMM6, XMM6   ; azzeramento registro per somma cumulata
	xorps XMM7, XMM7
	;da ora implemento il ciclo
	mov EDI, 0
		
	; di seguito si vuole calcolare (numeroValori)/4
	mov EDX, EBX  ; salvo in EDX il numero di valori contenuti nella colonna
	shr EDX, 4    ; shifting sinistro == divisione per 4. Ottengo tutte le quadruple che devo valutare nel ciclo  
	mov EAX, 0
	CDSCiclo: 
		cmp EDI, EDX ; ripetiamo il ciclo numeroValori/4
		jge CDSCicloFine
		mov ESI, EAX
		imul ESI, 16

		movaps XMM0, [ECX + ESI] ; prendo i successivi quattro valori INDICECOLONNA
		addps XMM7, XMM0 ;somma_dx += xi;
		mulps XMM0, XMM0 ;xi*xi
		addps XMM6, XMM0 ;somma_sx += xi*xi;
		

		add ESI, 16
		movaps XMM0, [ECX + ESI] ; prendo i successivi quattro valori INDICECOLONNA
		addps XMM7, XMM0 ;somma_dx += xi;
		mulps XMM0, XMM0 ;xi*xi
		addps XMM6, XMM0 ;somma_sx += xi*xi;
		
		add ESI, 16
		movaps XMM0, [ECX + ESI] ; prendo i successivi quattro valori INDICECOLONNA
		addps XMM7, XMM0 ;somma_dx += xi;
		mulps XMM0, XMM0 ;xi*xi
		addps XMM6, XMM0 ;somma_sx += xi*xi;

		add ESI, 16
		movaps XMM0, [ECX + ESI] ; prendo i successivi quattro valori INDICECOLONNA
		addps XMM7, XMM0 ;somma_dx += xi;
		mulps XMM0, XMM0 ;xi*xi
		addps XMM6, XMM0 ;somma_sx += xi*xi;

		add EAX, 4

		inc EDI 
		jmp CDSCiclo
		
		;dobbiamo vedere se devo eseguire il ciclo del residuo
	CDSCicloFine:
		imul EDI, 16 ; mi sto già spostando in avanti di EDI*4 valori. Questi valori sono quelli letti dal ciclo
		
		; effettuo le due somme parziali per sommare i 4 valodi del registro XMM3
		haddps XMM6, XMM6 
		haddps XMM6, XMM6
		haddps XMM7, XMM7 
		haddps XMM7, XMM7

	CDSCiclo2:
		cmp EDI, EBX
		jge CDSCicloFine2
		mov ESI, EDI
		imul ESI, 4
		movss XMM0, [ECX + ESI] ; prendo i successivi quattro valori INDICECOLONNA1
		addps XMM7, XMM0 ;somma_dx += xi;
		mulps XMM0, XMM0 ;xi*xi
		addps XMM6, XMM0 ;somma_sx += xi*xi;
		inc EDI
		jmp CDSCiclo2

	CDSCicloFine2:
		; serve per la return della sommatoria
		mov EAX, [EBP + somma_dx] ; prendo l'indirizzo di memoria su cui salvare il valore di somma_X += xi;
		extractps [EAX], XMM7, 0   ; salvo il valore in posizione 0 di xmm
		mov EAX, [EBP + somma_sx] ; prendo l'indirizzo di memoria su cui salvare il valore di somma_Y += yi;
		extractps [EAX], XMM6, 0   ; salvo il valore in posizione 0 di xmm

		; Sequenza di uscita dalla funzione
		pop	edi		; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp	; ripristina lo Stack Pointer
		pop	ebp		; ripristina il Base Pointer
		ret			; torna alla funzione C chiamante

;type calcoloPearsonCoefficient(params *input, int indiceColonna1, int indiceColonna2)
calcoloPearsonCoefficient2:
	push		ebp		; salva il Base Pointer
	mov		ebp, esp	; il Base Pointer punta al Record di Attivazione corrente
	push		ebx		; salva i registri da preservare
	push		esi
	push		edi

	; mi salvo i valori di utilità
	mov EAX, [EBP + input]	        ; mi prendo il parametro dell'indirizzo della struct contenente i parametri
	mov EBX, [EAX + 20]     		; mi salvo in EBX il valore di 'N'
	mov ECX, [EAX]                  ; metto in ECX il primo indirizzo del vettore ds  

	; in ECX va il seguente calcolo = indirizzo di memoria in cui comincia la colonna da valutare 
	; "indiceColonna" = indirizzoBaseVettoreDS+(indiceColonna*numeroRighe*4)
	mov EDX, [EBP + indiceColonna11] ; in EDX metto l'indice di colonna da valutare
	imul EDX, 4
	imul EDX, EBX
	add ECX, EDX

	mov ESI, [EAX] 

	mov EDX, [EBP + indiceColonna22] ; in EDX metto l'indice di colonna da valutare
	imul EDX, 4
	imul EDX, EBX
	add ESI, EDX

	mov EAX, ESI
	
	; implementazione della sommatoria
	; siccome bisogna sommare valori al quadrato allora la strategia è quella di salvare in due registri xmm
	; gli stessi valori e poi moltiplicarli	
	xorps XMM6, XMM6   ; azzeramento registro per somma cumulata
	xorps XMM7, XMM7
	xorps XMM5, XMM5
	xorps XMM4, XMM4
	xorps XMM3, XMM3
	;da ora implemento il ciclo
	mov EDI, 0
		
	; di seguito si vuole calcolare (numeroValori)/4
	mov EDX, EBX  ; salvo in EDX il numero di valori contenuti nella colonna
	shr EDX, 4    ; shifting sinistro == divisione per 4. Ottengo tutte le quadruple che devo valutare nel ciclo  

	mov ESI, 0
	CPCCiclo: 
		cmp EDI, EDX ; ripetiamo il ciclo numeroValori/4
		jge CPCCicloFine
		;mov ESI, EBX
		imul ESI, 16

		movaps XMM0, [ECX + ESI] ; prendo i successivi quattro valori INDICECOLONNA1
		movaps XMM1, [EAX + ESI] ; prendo i successivi quattro valori INDICECOLONNA2
		movaps XMM2, XMM0        ; copia valori INDICECOLONNA1
		
		addps XMM6, XMM0 ;somma_X += xi;
		addps XMM7, XMM1 ;somma_Y += yi;

		mulps XMM0, XMM1 ;xi*yi
		addps XMM5, XMM0 ;somma_XY += xi*yi;

		mulps XMM2, XMM2 ;xi*xi;
		addps XMM3, XMM2 ;somma_X2 += xi*xi;
		mulps XMM1, XMM1 ;yi*yi;
		addps XMM4, XMM1 ;somma_Y2 += yi*yi;
		

		add ESI, 16
		movaps XMM0, [ECX + ESI] ; prendo i successivi quattro valori INDICECOLONNA1
		movaps XMM1, [EAX + ESI] ; prendo i successivi quattro valori INDICECOLONNA2
		movaps XMM2, XMM0        ; copia valori INDICECOLONNA1

		addps XMM6, XMM0 ;somma_X += xi;
		addps XMM7, XMM1 ;somma_Y += yi;

		mulps XMM0, XMM1 ;xi*yi
		addps XMM5, XMM0 ;somma_XY += xi*yi;

		mulps XMM2, XMM2 ;xi*xi;
		addps XMM3, XMM2 ;somma_X2 += xi*xi;
		mulps XMM1, XMM1 ;yi*yi;
		addps XMM4, XMM1 ;somma_Y2 += yi*yi;

		add ESI, 16
		movaps XMM0, [ECX + ESI] ; prendo i successivi quattro valori INDICECOLONNA1
		movaps XMM1, [EAX + ESI] ; prendo i successivi quattro valori INDICECOLONNA2
		movaps XMM2, XMM0        ; copia valori INDICECOLONNA1

		addps XMM6, XMM0 ;somma_X += xi;
		addps XMM7, XMM1 ;somma_Y += yi;

		mulps XMM0, XMM1 ;xi*yi
		addps XMM5, XMM0 ;somma_XY += xi*yi;

		mulps XMM2, XMM2 ;xi*xi;
		addps XMM3, XMM2 ;somma_X2 += xi*xi;
		mulps XMM1, XMM1 ;yi*yi;
		addps XMM4, XMM1 ;somma_Y2 += yi*yi;

		add ESI, 16
		movaps XMM0, [ECX + ESI] ; prendo i successivi quattro valori INDICECOLONNA1
		movaps XMM1, [EAX + ESI] ; prendo i successivi quattro valori INDICECOLONNA2
		movaps XMM2, XMM0        ; copia valori INDICECOLONNA1

		addps XMM6, XMM0 ;somma_X += xi;
		addps XMM7, XMM1 ;somma_Y += yi;

		mulps XMM0, XMM1 ;xi*yi
		addps XMM5, XMM0 ;somma_XY += xi*yi;

		mulps XMM2, XMM2 ;xi*xi;
		addps XMM3, XMM2 ;somma_X2 += xi*xi;
		mulps XMM1, XMM1 ;yi*yi;
		addps XMM4, XMM1 ;somma_Y2 += yi*yi;
		inc EDI
		imul ESI, EDI, 4
		 
		jmp CPCCiclo
		
		;dobbiamo vedere se devo eseguire il ciclo del residuo
	CPCCicloFine:
		imul EDI, 16 ; mi sto già spostando in avanti di EDI*4 valori. Questi valori sono quelli letti dal ciclo
		
		; effettuo le due somme parziali per sommare i 4 valodi del registro XMM3
		haddps XMM6, XMM6 
		haddps XMM6, XMM6
		haddps XMM7, XMM7 
		haddps XMM7, XMM7
		haddps XMM5, XMM5 
		haddps XMM5, XMM5
		haddps XMM4, XMM4 
		haddps XMM4, XMM4
		haddps XMM3, XMM3 
		haddps XMM3, XMM3

	CPCCiclo2:
		cmp EDI, EBX
		jge CPCCicloFine2
		mov ESI, EDI
		imul ESI, 4
		movss XMM0, [ECX + ESI] ; prendo i successivi quattro valori INDICECOLONNA1
		movss XMM1, [EAX + ESI] ; prendo i successivi quattro valori INDICECOLONNA2
		movaps XMM2, XMM0        ; copia valori INDICECOLONNA1

		addps XMM6, XMM0 ;somma_X += xi;
		addps XMM7, XMM1 ;somma_Y += yi;

		mulps XMM0, XMM1 ;xi*yi
		addps XMM5, XMM0 ;somma_XY += xi*yi;

		mulps XMM2, XMM2 ;xi*xi;
		addps XMM3, XMM2 ;somma_X2 += xi*xi;
		mulps XMM1, XMM1 ;yi*yi;
		addps XMM4, XMM1 ;somma_Y2 += yi*yi;

		inc EDI
		jmp CPCCiclo2

	CPCCicloFine2:

		; serve per la return della sommatoria
		mov EAX, [EBP + risultatoY] ; prendo l'indirizzo di memoria su cui salvare il valore di somma_X += xi;
		extractps [EAX], XMM7, 0   ; salvo il valore in posizione 0 di xmm
		mov EAX, [EBP + risultatoX] ; prendo l'indirizzo di memoria su cui salvare il valore di somma_Y += yi;
		extractps [EAX], XMM6, 0   ; salvo il valore in posizione 0 di xmm
		mov EAX, [EBP + risultatoXY] ; prendo l'indirizzo di memoria su cui salvare il valore di somma_XY += xi*yi;
		extractps [EAX], XMM5, 0   ; salvo il valore in posizione 0 di xmm
		mov EAX, [EBP + risultatoXX] ; prendo l'indirizzo di memoria su cui salvare il valore di somma_X2 += xi*xi;
		extractps [EAX], XMM3, 0   ; salvo il valore in posizione 0 di xmm
		mov EAX, [EBP + risultatoYY] ; prendo l'indirizzo di memoria su cui salvare il valore di somma_Y2 += yi*yi;
		extractps [EAX], XMM4, 0   ; salvo il valore in posizione 0 di xmm
		; Sequenza di uscita dalla funzione
		pop	edi		; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp	; ripristina lo Stack Pointer
		pop	ebp		; ripristina il Base Pointer
		ret			; torna alla funzione C chiamante

    
calcoloValori_n0n1:
    ; salvataggio dei registri generali
    push ebp
    mov ebp, esp
    push ebx
    push esi
    push edi
    
    ; caricamento dei parametri
    mov eax, [ebp + input] ;input
    mov ecx, [eax + 4] 	   ;labels
    
    ; calcolo il numero di iterazioni del ciclo - Effettuo la divisione per 4 in quanto voglio prendere il numero totale di quartetti
    mov ebx, [eax+20]
    shr ebx,2
    
    ; inizializzo un contatore
    mov edi, 0
    mov edx, 0
    
    ; Azzerro il registro XMM0 e XMM1. XMM0 conterrà man mano i valori caricati invece XMM1 conterrà le somme
    ;pxor XMM0,XMM0
    pxor XMM1,XMM1
    
    ;Ad ogni iterazione del ciclo verranno caricati 4 valori sul registro XMM0. Si effettuerà la somma con XMM1 che ogni volta avrà nelle varie celle di memoria le somme. 
    ;Dunque, ad ogni passaggio effettueremo solo la somma degli 1 e non degli 0 ritrovandoci alla fine in XMM1 tutte le somme parziali degli 1. Effettuando le dovute HADDPS avremo in prima posizione in XMM1 le somme degli 1
    ;per calcolare il quantitativo di 0 sarà sufficiente fare il TOTALE-QUANTITA_UNO
    
    ciclo_calcoloValorin0n1:
    	cmp edi,ebx
    	JGE residui_calcoloValorin0n1
    	
    	;caricamento dei valori in XMM0
    	MOV edx, edi
    	IMUL edx, 16
    	MOVAPS XMM0, [ecx + edx]
    	ADDPS XMM1,XMM0
    	inc edi
    	
    	JMP ciclo_calcoloValorin0n1
    
    residui_calcoloValorin0n1:
    	HADDPS XMM1,XMM1
    	HADDPS XMM1,XMM1
    	
    	;verifico se devo calcolare i residui
    	shl ebx,2
    	mov edx, [eax+20] ;valore di N
    	sub edx,ebx ;Verifica dei residui TOTALE-VALORI CONSIDERATI
    	
    	cmp edx,0
    	JNE calcolo_residuiCalcoloValorin0n1
    	
    calcolo_residuiCalcoloValorin0n1:
    	mov edi,0 ;contatore
    		
    	;inizio del ciclo
    ciclo_residuiCalcoloValorin0n1:
    	cmp edi,edx
    	JGE fine_calcoloValorin0n1

    	; Calcolo dell'offset nel registro esi
    	MOV ebx, edi
    	IMUL ebx, 16
    	ADD ebx, esi

		; Effettua la somma degli ultimi elementi
		ADDPS XMM1, [ebx]
		INC edi
		jmp ciclo_residuiCalcoloValorin0n1
		 
    fine_calcoloValorin0n1:
    	mov esi, [ebp + 12]  ; Carica il puntatore al risultato di n0
    	mov ebx, [ebp + 16]  ; Carica il puntatore al risultato di n1
    	CVTTSS2SI ecx, XMM1  ;Converte il float in intero
    	mov [ebx],ecx ; Salva il risultato di n1 sul puntatore associato        
    	  
	    ;calcoliamo il valore n1
        mov edx, [eax+20] ;valore di N
    	sub edx,ecx ; Totale valori-Totale valori di n1
    	mov [esi],edx
    
	    pop edi
	    pop esi
	    pop ebx
	    mov esp, ebp
	    pop ebp
	    ret


calcoloMediaGruppoAssembly:
	; salvataggio dei registri generali
    push ebp
    mov ebp, esp
    push ebx
    push esi
    push edi

    ; caricamento dei parametri
    mov eax, [ebp + input]
    mov ecx, [ebp + j] ; carico j in ecx
 
    ; caricamento del valore di N
    mov	ebx, [eax + 20]

    ; caricamento del valore di Ds - Si moltiplica per spostarsi all'inizio della colonna facendo NUMERO_COLONNA*LUNGHEZZA_COLONNA*GRANDEZZA_VALORI
    mov edx, [eax] ; carico in edx il valore del dataset
    imul ecx,4
    imul ecx,ebx 
    add edx,ecx

    ;inizializzazione dei vettori che si occuperanno di tenere le somme.
    XORPS XMM3,XMM3 ;somma parziali degli 1
    XORPS XMM4,XMM4 ;somma parziali degli 0
	movddup XMM5, [_one] ; mi copio quattro volte il valore mediaValoriColonna
    shufps XMM5, XMM5, 0 ;Serve per effettuare lo xor con le labels in modo tale da ottenere il vettore al contrario e poter considerare gli 0 

    ;inizializzazione di un contatore
    mov edi,0
    ;calcolo del numero di iterazioni del for
    shr ebx,2
    ;carico il valore di labels
    mov ecx, [eax+4]
    ;Ad ogni iterazione del ciclo vengono caricati i valori di label in XMM0 e del dataset in XMM1. Tramite la moltiplicazione verranno caricati i valori corrispondenti al valore di 1 nel labels considerato, e andranno in XMM2.
    ;Verranno successivamente sommati in XMM3. Si effettua nuovamente la moltiplicazione solo che si fa prima lo xor tra le labels e il vettore di 1 in modo tale che si invertano gli 1 e gli 0. Infine in XMM3 e XMM4 ci saranno 
    ;le somme parziali che verranno sommate tramite HADDPS
    ciclo_mediaGruppo:
   	cmp edi,ebx
		jge residui_mediaGruppo
		MOV esi, edi
		IMUL esi, 16
		MOVAPS XMM0, [ecx + esi] ;carico il quartetto associato alle Labels

		MOVAPS XMM1, [edx + esi] ;carico il quartetto associato al Dataset | usato per gruppo 0
		MOVAPS XMM2, XMM1 ;carico il quartetto associato al Dataset | usato per gruppo 1
		MULPS XMM2, XMM0 ;in XMM2 ci sono i valori corrispondenti della maschera  GRUPPO 1 * MASCHERA
		XORPS XMM0, XMM5 ;Inversione della maschera
		MULPS XMM1, XMM0 ;in XMM1 ci sono i valori corrispondenti alla maschera GRUPPO 0 * MASCHERA
		ADDPS XMM3, XMM2 ;sommo i valori del gruppo 1 0|0|0|0  +  num1|0|num2|0 + 0|num3|0|num4 =
		ADDPS XMM4, XMM1 ;sommo i valori del gruppo 0
		inc edi
		jmp ciclo_mediaGruppo

    residui_mediaGruppo:
    	HADDPS XMM3,XMM3
    	HADDPS XMM3,XMM3
    	HADDPS XMM4,XMM4
    	HADDPS XMM4,XMM4
    	
    	imul edi,4
    	mov ebx, [eax + 20]
    	
    	residui_calcoloMediaGruppo:
    		cmp edi,ebx
    		jge fine_mediaGruppo
    		
    		mov esi, edi
		imul esi, 4 
    		movss xmm0,[ecx+esi]
    		movss xmm1,[edx+esi]
    		movss xmm2,[edx+esi]
    		
    		mulss xmm1,xmm0
    		pxor xmm0,xmm5
    		mulss xmm2,xmm0
    		
    		addss xmm3,xmm1
    		addss xmm4,xmm2
    		
    		inc edi
    		jmp residui_calcoloMediaGruppo
    	
    	
    	fine_mediaGruppo:
    		mov eax, [ebp+16] ;valore di n0
		mov ebx, [ebp+20] ;valore di n1
		CVTSI2SS XMM1, eax ;valore di n0 in XMM1
		CVTSI2SS XMM2, ebx ;valore di n1 in XMM2
		DIVSS XMM3, XMM2
		DIVSS XMM4, XMM1
		SUBSS XMM4,XMM3
		mov eax, [ebp+24]
		movss [eax],XMM4

    pop edi
    pop esi
    pop ebx
    mov esp, ebp
    pop ebp
    ret