; ---------------------------------------------------------
; Regression con istruzioni AVX a 64 bit
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
;     nasm -f elf64 regression64.nasm
;

%include "sseutils64.nasm"

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
    _ones dq 1.0
section .bss			; Sezione contenente dati non inizializzati
	alignb 32
	sc		resq		1

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
	mov	rdi, %1
	mov	rsi, %2
	call	get_block
%endmacro

%macro	fremem	1
	mov	rdi, %1
	call	free_block
%endmacro

; ------------------------------------------------------------
; Funzione prova
; ------------------------------------------------------------
global sqrtAssembly
;global calcoloSommatoriaQuadrato
;global calcoloMediaAssembly
;global calcoloSommatoria
global prova
global calcoloValori_n0n1
global calcoloMediaGruppoAssembly
global calcoloPearsonCoefficient2
global calcoloDeviazioneStandard2


msg	db 'sc:',32,0
nl	db 10,0

prova:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq						; salva i registri generali

		; ------------------------------------------------------------
		; I parametri sono passati nei registri
		; ------------------------------------------------------------
		; rdi = indirizzo della struct input
		
		; esempio: stampa input->sc
        ; [RDI] input->ds; 			// dataset
		; [RDI + 8] input->labels; 	// etichette
		; [RDI + 16] input->out;	// vettore contenente risultato dim=(k+1)
		; [RDI + 24] input->sc;		// score dell'insieme di features risultato
		; [RDI + 32] input->k; 		// numero di features da estrarre
		; [RDI + 36] input->N;		// numero di righe del dataset
		; [RDI + 40] input->d;		// numero di colonne/feature del dataset
		; [RDI + 44] input->display;
		; [RDI + 48] input->silent;
		VMOVSD		XMM0, [RDI+24]
		VMOVSD		[sc], XMM0
		prints 		msg
		printsd		sc
		prints 		nl
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		
		popaq				; ripristina i registri generali
		mov		rsp, rbp	; ripristina lo Stack Pointer
		pop		rbp		; ripristina il Base Pointer
		ret				; torna alla funzione C chiamante

;											RDI			    RSI           RDX           RCX             R8             R9
;extern void calcoloPearsonCoefficient2(valori *input, type* sommaX, type* sommaY, type* sommaXY, type* sommaX2, type* sommaY2);
calcoloPearsonCoefficient2:
	push		rbp				; salva il Base Pointer
	mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
	pushaq						; salva i registri generali

	; in ECX va il seguente calcolo = indirizzo di memoria in cui comincia la colonna da valutare 
	; "indiceColonna" = indirizzoBaseVettoreDS+(indiceColonna*numeroRighe*4)
	mov R15, [RDI + 36]
	mov R10, [RDI]
	;mov EBX, [RDI+8]
	imul ESI, 8
	imul ESI, r15d
	add R10, RSI ;COL1

	mov R14, [RDI]
	;mov EBX, [RDI + 12]
	imul EDX, 8
	imul EDX, r15d
	add R14, RDX ;COL2
	
	; implementazione della sommatoria
	; siccome bisogna sommare valori al quadrato allora la strategia è quella di salvare in due registri xmm
	; gli stessi valori e poi moltiplicarli	
	vxorps YMM6, YMM6   ; azzeramento registro per somma cumulata
	vxorps YMM7, YMM7
	vxorps YMM5, YMM5
	vxorps YMM4, YMM4
	vxorps YMM3, YMM3
	;da ora implemento il ciclo
	mov R11, 0
	mov R13, 0
		
	; di seguito si vuole calcolare (numeroValori)/4
	shr R15D, 4
	CPCCiclo: 
		cmp R11, R15 ; ripetiamo il ciclo numeroValori/4
		jge CPCCicloFine
		mov R12, R13
		imul R12, 32
		
		vmovapd YMM0, [R10 + R12] ; prendo i successivi quattro valori INDICECOLONNA1
		vmovapd YMM1, [R14 + R12] ; prendo i successivi quattro valori INDICECOLONNA2
		vmovapd YMM2, YMM0        ; copia valori INDICECOLONNA1
		
		vaddpd YMM6, YMM0 ;somma_X += xi;
		vaddpd YMM7, YMM1 ;somma_Y += yi;

		vmulpd YMM0, YMM1 ;xi*yi
		vaddpd YMM5, YMM0 ;somma_XY += xi*yi;

		vmulpd YMM2, YMM2 ;xi*xi;
		vaddpd YMM3, YMM2 ;somma_X2 += xi*xi;
		vmulpd YMM1, YMM1 ;yi*yi;
		vaddpd YMM4, YMM1 ;somma_Y2 += yi*yi;
		

		add R12, 32
		vmovapd YMM0, [R10 + R12] ; prendo i successivi quattro valori INDICECOLONNA1
		vmovapd YMM1, [R14 + R12] ; prendo i successivi quattro valori INDICECOLONNA2
		vmovapd YMM2, YMM0        ; copia valori INDICECOLONNA1

		vaddpd YMM6, YMM0 ;somma_X += xi;
		vaddpd YMM7, YMM1 ;somma_Y += yi;

		vmulpd YMM0, YMM1 ;xi*yi
		vaddpd YMM5, YMM0 ;somma_XY += xi*yi;

		vmulpd YMM2, YMM2 ;xi*xi;
		vaddpd YMM3, YMM2 ;somma_X2 += xi*xi;
		vmulpd YMM1, YMM1 ;yi*yi;
		vaddpd YMM4, YMM1 ;somma_Y2 += yi*yi;

		add R12, 32
		vmovapd YMM0, [R10 + R12] ; prendo i successivi quattro valori INDICECOLONNA1
		vmovapd YMM1, [R14 + R12] ; prendo i successivi quattro valori INDICECOLONNA2
		vmovapd YMM2, YMM0        ; copia valori INDICECOLONNA1

		vaddpd YMM6, YMM0 ;somma_X += xi;
		vaddpd YMM7, YMM1 ;somma_Y += yi;

		vmulpd YMM0, YMM1 ;xi*yi
		vaddpd YMM5, YMM0 ;somma_XY += xi*yi;

		vmulpd YMM2, YMM2 ;xi*xi;
		vaddpd YMM3, YMM2 ;somma_X2 += xi*xi;
		vmulpd YMM1, YMM1 ;yi*yi;
		vaddpd YMM4, YMM1 ;somma_Y2 += yi*yi;

		add R12, 32
		vmovapd YMM0, [R10 + R12] ; prendo i successivi quattro valori INDICECOLONNA1
		vmovapd YMM1, [R14 + R12] ; prendo i successivi quattro valori INDICECOLONNA2
		vmovapd YMM2, YMM0        ; copia valori INDICECOLONNA1

		vaddpd YMM6, YMM0 ;somma_X += xi;
		vaddpd YMM7, YMM1 ;somma_Y += yi;

		vmulpd YMM0, YMM1 ;xi*yi
		vaddpd YMM5, YMM0 ;somma_XY += xi*yi;

		vmulpd YMM2, YMM2 ;xi*xi;
		vaddpd YMM3, YMM2 ;somma_X2 += xi*xi;
		vmulpd YMM1, YMM1 ;yi*yi;
		vaddpd YMM4, YMM1 ;somma_Y2 += yi*yi;
		
		add R13, 4
		inc R11
		jmp CPCCiclo
		
		;dobbiamo vedere se devo eseguire il ciclo del residuo
	CPCCicloFine:
		
		imul R11, 16 ; mi sto già spostando in avanti di EDI*4 valori. Questi valori sono quelli letti dal ciclo
	
		; effettuo le due somme parziali per sommare i 4 valodi del registro XMM3
		vhaddpd YMM6,YMM6 ;haddps XMM6, XMM6 haddps XMM6, XMM6 ;somma_X += xi;
    	vextractf128 xmm8,ymm6,0
    	vextractf128 xmm9,ymm6,1
    	vaddpd xmm8,xmm9 ;in prima posizione ci sarà la somma del gruppo 1

		vhaddpd YMM7,YMM7 ;haddps XMM7, XMM7 haddps XMM7, XMM7 ;somma_Y += yi;
    	vextractf128 xmm10,ymm7,0
    	vextractf128 xmm11,ymm7,1
    	vaddpd xmm10,xmm11 ;in prima posizione ci sarà la somma del gruppo 0
		
		vhaddpd YMM5,YMM5 ;haddps XMM5, XMM5 haddps XMM5, XMM5 ;somma_XY += xi*yi;
    	vextractf128 xmm12,ymm5,0
    	vextractf128 xmm13,ymm5,1
    	vaddpd xmm12,xmm13 ;in prima posizione ci sarà la somma del gruppo 0
		
		vhaddpd YMM4,YMM4 ;haddps XMM4, XMM4 haddps XMM4, XMM4 ;somma_Y2 += yi*yi;
    	vextractf128 xmm14,ymm4,0
    	vextractf128 xmm15,ymm4,1
    	vaddpd xmm14,xmm15 ;in prima posizione ci sarà la somma del gruppo 0

		vxorps YMM6, YMM6
		vxorps YMM7, YMM7
		vhaddpd YMM3,YMM3 ;haddps XMM3, XMM3 haddps XMM3, XMM3 ;somma_X2 += xi*xi;
    	vextractf128 xmm6,ymm3,0
    	vextractf128 xmm7,ymm3,1
    	vaddpd xmm6,xmm7 ;in prima posizione ci sarà la somma del gruppo 0
		
		vxorps YMM0, YMM0
		vxorps YMM1, YMM1
		mov r15d, [RDI+36]

	CPCCiclo2:
		cmp R11, r15
		jge CPCCicloFine2
		mov R13, R11
		imul R13, 8

		vmovsd XMM0, [R10 + R13] ; prendo i successivi quattro valori INDICECOLONNA1
		vmovsd XMM1, [R14 + R13] ; prendo i successivi quattro valori INDICECOLONNA2
		vmovapd XMM2, XMM0        ; copia valori INDICECOLONNA1

		vaddsd XMM8, XMM0 ;somma_X += xi;
		vaddsd XMM10, XMM1 ;somma_Y += yi;

		vmulsd XMM0, XMM1 ;xi*yi
		vaddsd XMM12, XMM0 ;somma_XY += xi*yi;

		vmulsd XMM2, XMM2 ;xi*xi;
		vaddsd XMM6, XMM2 ;somma_X2 += xi*xi;
		vmulsd XMM1, XMM1 ;yi*yi;
		vaddsd XMM14, XMM1 ;somma_Y2 += yi*yi;

		inc R11
		jmp CPCCiclo2

	CPCCicloFine2:
	
	;type numeratore = input->N*somma_XY - (somma_X * somma_Y);
		vmovsd XMM2, XMM6
		vmovsd XMM3, XMM8 
		vmovsd XMM0, XMM8 ;somma_X += xi;
		vmovsd XMM1, XMM10 ;somma_Y += yi;
		
		cvtsi2sd  xmm13, r15
		
		vmulsd XMM12, xmm13 ;input->N*somma_XY
		vmulsd XMM0, XMM1   ;(somma_X * somma_Y)
		vsubsd xmm12, xmm0  ;in xmm12 c'è il numeratore
	;type denominatore = sqrt(((input->N*somma_X2) - (somma_X*somma_X)) *  ((input->N*somma_Y2) - (somma_Y*somma_Y)));

		vmulsd XMM6, xmm13    ;(input->N*somma_X2)
		vmulsd XMM8, XMM8   ;(somma_X*somma_X)
		vmulsd XMM14, xmm13   ;(input->N*somma_Y2)
		vmulsd XMM10, XMM10 ;(somma_Y*somma_Y)
	    vsubsd XMM6, XMM8   ;((input->N*somma_X2) - (somma_X*somma_X))
		vsubsd XMM14, XMM10 ;((input->N*somma_Y2) - (somma_Y*somma_Y))
		vmulsd XMM6, XMM14  ;(somma_X*somma_X)) *  ((input->N*somma_Y2)
		sqrtsd XMM6, XMM6
		vdivsd XMM12, XMM6

		movsd [RCX], xmm3 ;somma_XY += xi*yi;
		movsd [R8], xmm2   ;somma_X2 += xi*xi;
		movsd [R9], XMM12  ;somma_Y2 += yi*yi;

		popaq				; ripristina i registri generali
		mov		rsp, rbp	; ripristina lo Stack Pointer
		pop		rbp		; ripristina il Base Pointer
		ret				; torna alla funzione C chiamante



;type calcoloDeviazioneStandard(params *input, int j, type* somma_sx, type* somma_dx)
calcoloDeviazioneStandard2:
	push		rbp				; salva il Base Pointer
	mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
	pushaq						; salva i registri generali

	
	; mi salvo i valori di utilità
	mov R15, [RDI + 36]     		; mi salvo in EBX il valore di 'N'
	mov RAX, [RDI]                  ; metto in ECX il primo indirizzo del vettore ds  

	; in ECX va il seguente calcolo = indirizzo di memoria in cui comincia la colonna da valutare 
	; "indiceColonna" = indirizzoBaseVettoreDS+(indiceColonna*numeroRighe*4)
	imul ESI, 8
	imul ESI, R15d
	add RAX, RSI

	; implementazione della sommatoria
	; siccome bisogna sommare valori al quadrato allora la strategia è quella di salvare in due registri xmm
	; gli stessi valori e poi moltiplicarli	
	vxorps YMM6, YMM6   ; azzeramento registro per somma cumulata
	vxorps YMM7, YMM7
	;da ora implemento il ciclo
	mov R11,0
	mov R10, 0
	; di seguito si vuole calcolare (numeroValori)/4
	shr R15D, 4
	CDSCiclo: 
		cmp R11, R15 ; ripetiamo il ciclo numeroValori/4
		jge CDSCicloFine
		mov R12, R10
		imul R12, 32

		vmovapd  YMM0, [RAX + R12] ; prendo i successivi quattro valori INDICECOLONNA
		vaddpd  YMM7, YMM0 ;somma_dx += xi;
		vmulpd  YMM0, YMM0 ;xi*xi
		vaddpd  YMM6, YMM0 ;somma_sx += xi*xi;

		add R12, 32
		vmovapd  YMM0, [RAX + R12] ; prendo i successivi quattro valori INDICECOLONNA
		vaddpd  YMM7, YMM0 ;somma_dx += xi;
		vmulpd  YMM0, YMM0 ;xi*xi
		vaddpd  YMM6, YMM0 ;somma_sx += xi*xi;

		add R12, 32
		vmovapd  YMM0, [RAX + R12] ; prendo i successivi quattro valori INDICECOLONNA
		vaddpd  YMM7, YMM0 ;somma_dx += xi;
		vmulpd  YMM0, YMM0 ;xi*xi
		vaddpd  YMM6, YMM0 ;somma_sx += xi*xi;

		add R12, 32
		vmovapd  YMM0, [RAX + R12] ; prendo i successivi quattro valori INDICECOLONNA
		vaddpd  YMM7, YMM0 ;somma_dx += xi;
		vmulpd  YMM0, YMM0 ;xi*xi
		vaddpd  YMM6, YMM0 ;somma_sx += xi*xi;

		add R10, 4

		inc R11 
		jmp CDSCiclo
		
		;dobbiamo vedere se devo eseguire il ciclo del residuo
	CDSCicloFine:
		imul R11, 16 ; mi sto già spostando in avanti di EDI*4 valori. Questi valori sono quelli letti dal ciclo
		
		; effettuo le due somme parziali per sommare i 4 valodi del registro XMM3
		vhaddpd YMM6,YMM6          ;somma_sx += xi*xi;
    	vextractf128 xmm8,ymm6,0
    	vextractf128 xmm9,ymm6,1
    	vaddpd xmm8,xmm9 ;in prima posizione ci sarà la somma del gruppo 1
    			
    	vhaddpd YMM7,YMM7          ;somma_dx += xi;
    	vextractf128 xmm10,ymm7,0
    	vextractf128 xmm11,ymm7,1
    	vaddpd xmm10,xmm11 ;in prima posizione ci sarà la somma del gruppo 0

		vxorps YMM0, YMM0
		mov r15d,[rdi+36]

	CDSCiclo2:
		cmp R11, r15
		jge CDSCicloFine2
		mov R13, R11
		imul R13, 8
		vmovsd XMM0, [RAX + R13] ; prendo i successivi quattro valori INDICECOLONNA1
		vaddsd XMM10, XMM0 ;somma_dx += xi;
		vmulsd XMM0, XMM0 ;xi*xi
		vaddsd XMM8, XMM0 ;somma_sx += xi*xi;
		inc R11
		jmp CDSCiclo2

	CDSCicloFine2:
		;sqrt((somma_sx - (somma_dx * somma_dx/n))/(n-1));
		;vmovsd XMM0, xmm10
		;mov ECX, [RDI + 36]
		;cvtsi2sd xmm15, rcx
		;vdivsd XMM10, XMM15 ;somma_dx/n
		;vmulsd XMM0, XMM10  ;(somma_dx * somma_dx/n)+
		;vsubsd XMM8, XMM0
		;dec ECX
		;cvtsi2sd xmm15, RCX
		;vdivsd XMM8, XMM15
		;sqrtsd XMM8, XMM8

		movsd [RDX], xmm8
		movsd [RCX], xmm10
		
		; Sequenza di uscita dalla funzione
		popaq				; ripristina i registri generali
		mov		rsp, rbp	; ripristina lo Stack Pointer
		pop		rbp		; ripristina il Base Pointer
		ret				; torna alla funzione C chiamante


calcoloMediaGruppoAssembly:
		push rbp	; salva il Base Pointer
		mov  rbp, rsp	; il Base Pointer punta al Record di Attivazione corrente
		pushaq		; salva i registri generali
		
		;caricamento dei valori 
		;in RDI è salvato il paarmetro input
		;in RSI è salvato il valore j
		;in RDX è salvato il valore di n0
		;in RCX è salvato il valore di n1
		;in R8 è salvato il puntatore di differenzaMedia
		mov rbx,[rdi] ;caricamento dataset
		
		mov r10, [rdi+8]  ;caricamento dei labels
		
		;posizionamento del puntatore nel dataset
		imul esi,8
		imul esi,[RDI+36]
		add rbx,rsi
			
		;calcolo del numero di iterazioni del ciclo
		mov r15,0
		mov r15,[rdi+36] ;valore di N
		shr r15d,2
		
		;prepariamo il contatore
		mov r11,0
		
		;azzeriamo i registri
		vxorps YMM0,YMM0
		vxorps YMM1,YMM1
		vxorps YMM2,YMM2
		vxorps YMM3,YMM3
		vxorps YMM4,YMM4
		vxorps YMM5,YMM5
		vmovddup ymm5, [_ones]
		vbroadcastsd ymm5, xmm5
		
	ciclo_mediaGruppo:
		cmp r11,r15
		jge residui_mediaGruppo
		MOV r12,r11
		IMUL r12,32
		vmovapd YMM0,[r10+r12]   ;caricamento valori label
		vmovapd YMM1,[rbx+r12]  ;caricamento dei dati del dataset
		vmovapd YMM2,YMM1
    			
    	vmulpd YMM1,YMM0
		vxorpd YMM0,YMM5
    	vmulpd YMM2,YMM0
			    			
		vaddpd YMM3,YMM1 ;somma gruppo 1
		vaddpd YMM4,YMM2 ;somma gruppo 0
    			
		inc r11
		jmp ciclo_mediaGruppo
    			
    residui_mediaGruppo:
    	vhaddpd YMM3,YMM3
    	vextractf128 xmm8,ymm3,0
    	vextractf128 xmm9,ymm3,1
    	vaddpd xmm8,xmm9 ;in prima posizione ci sarà la somma del gruppo 1
    			
    	vhaddpd YMM4,YMM4
    	vextractf128 xmm10,ymm4,0
    	vextractf128 xmm11,ymm4,1
    	vaddpd xmm10,xmm11 ;in prima posizione ci sarà la somma del gruppo 0
    			
    	;verifichiamo se abbiamo dei residui
    	imul r11,4    	
		mov r15d,[rdi+36]
		    	
	vxorps ymm2, ymm2
	vxorps ymm0, ymm0
	vxorps ymm1, ymm1
	
	ciclo_calcoloMediaGruppoResidui:
		cmp r11, r15
		jge fine_MediaGruppo
		mov r14, r11
		imul r14, 8
		
		vmovsd xmm1, [rbx + r14] ; prendo il successivo valore
		vaddsd xmm2, xmm1
		
		vmovsd xmm0, [r10+r14]   ;caricamento valori label
		vmovsd xmm1, [rbx+r14]  ;caricamento dei dati del dataset
		vmovsd xmm2, xmm1
		
		vmulsd xmm1,xmm0
		pxor xmm0,xmm5
    	vmulsd xmm2,xmm0
			    			
		vaddsd xmm8,xmm1 ;somma gruppo 1
		vaddpd xmm10,xmm2 ;somma gruppo 0
		
		inc r11
		jmp ciclo_calcoloMediaGruppoResidui
			
	fine_MediaGruppo:	
		cvtsi2sd xmm1, edx ;valore di n0
		vbroadcastsd ymm1, xmm1
				
		vcvtsi2sd xmm2, ecx ;valore di n1
		vbroadcastsd ymm2, xmm2
			
		vdivpd xmm12, xmm8, xmm2 ;GRUPPO 1
		vdivpd xmm13, xmm10, xmm1 ;GRUPPO 0
	   	vsubpd xmm12,xmm13 ;DifferenzaMedia
	    		
	   	vextractf128 xmm0,ymm12,0
		movsd [R8],xmm12 
		
		popaq				; ripristina i registri generali
		mov		rsp, rbp	; ripristina lo Stack Pointer
		pop		rbp		; ripristina il Base Pointer
		ret				; torna alla funzione C chiamante		
   
   
;extern void calcoloValori_n0n1(params* input, int* n0, int* n1);
calcoloValori_n0n1:
	push rbp	; salva il Base Pointer
	mov  rbp, rsp	; il Base Pointer punta al Record di Attivazione corrente
	pushaq		; salva i registri generali

    ; lettura dei parametri dal Recordi di Attivazione
    mov eax, [rdi+36]  ;estraiamo il numero di righe del dataset
    mov rbx, [rdi+8] ;estraiamo l'indirizzo delle labels
	mov r8, rsi	;indirizzo n0
	mov r9, rdx	;indirizzo n1
    
    ; corpo della funzione
    vpxor ymm0, ymm0 ;azzeriamo il registro che conterrà la somma delle labels
        ;l'obiettivo è sommare le labels così da ottenere n1
        ;poi per ottenere n0 si fa #labels-n1
        ;dunque occorre definire un puntatore che tenga il
        ;numero delle labels, ovvero rax

        ;serve, poi, la variabile che terrà conto
        ;delle iterazioni correntemente eseguite
        ;nel ciclo
    xor edx, edx
	xor r10, r10   ;r10 è esattamente come edx solo che serve per fare i calcoli a 64 bit, perchè il contatore ci serve sia a 32 (per la compare con eax) che a 64 (per l'accesso in memoria con rbx)
        ;in fine eseguiamo il ciclo
    ciclo_calcoloValoren0n1:
        cmp edx, eax    
        JGE fine_ciclo_calcoloValoren0n1
        vmovapd ymm1, [rbx+r10*8]    ;prendiamo 4 labels da 64 bit, quindi il prossimo indirizzo utile sta a 8 byte di distanza
        vaddpd ymm0, ymm1   ;le sommiamo alle vecchie
        add edx, 4  ;abbiamo appena preso 4 valori da 64 bit, quindi incrementiamo
                    ;il contatore di 4.    
		add r10, 4 
        JMP ciclo_calcoloValoren0n1  ;chiudiamo il for

    fine_ciclo_calcoloValoren0n1:
	
    ;usiamo la riduzione prima di calcolare i residui
    vhaddpd ymm0, ymm0
    ;ora dovremmo sommare i 64 bit più significativi con i 64 meno significiativi,
    ;per farlo procediamo come dis eguito

    ; Estrai i 128 bit più significativi da ymm0
    vextractf128 xmm1, ymm0, 1

    ;sommiamo ciò detto prima
    addpd xmm0, xmm1

    ;potrebbero esserci dei residui da calcolare
    cmp edx, eax
    JE fine_residui_calcoloValoren0n1   ;se r8==rax allora non possiamo avere residui
    sub edx, 3  ;ripristiniamo il contatore
	sub r10, 3

    residui_calcoloValoren0n1:  ;esattamente come il for di sopra ma prende elementi da 64 bit singoli
        CMP edx, eax
        JGE fine_residui_calcoloValoren0n1
        movsd xmm1, [rbx+r10*8]
        addsd xmm0, xmm1
        add edx, 1
		add r10, 1
        JMP residui_calcoloValoren0n1

    fine_residui_calcoloValoren0n1:
    cvtsd2si ecx, xmm0 ;converto il valore di n1 in intero
    sub eax, ecx ;calcoliamo il valore di n0 (rax>=rcx sempre)

    ;restituiamo i risultati
    mov [r8], eax
    mov [r9], ecx
	
    popaq
    mov rsp,rbp
    pop rbp
    ret

sqrtAssembly:
    ; sequenza di ingresso nella funzione
    push rbp ; salva il Base Pointer
    mov rbp, rsp ; il Base Pointer punta al Record di Attivazione corrente

    ; lettura dei parametri dal Recordi di Attivazione
	;in xmm0 c'è il parametro di tipo double
    
	; corpo della funzione
    sqrtsd xmm0, xmm0 ;calcoliamo la radice quadrata
		;il valore viene restituito nei registri xmm0 e xmm1

    ; sequenza di uscita dalla funzione
    mov rsp, rbp ; ripristina lo Stack Pointer
    pop rbp ; ripristina il Base Pointer
    ret ; ritorna alla funzione chiamant