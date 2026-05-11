# Feature-Selection-algorithm
Sviluppo di un sistema di Correlation Feature Selection (CFS) in C, ottimizzato per calcolo ad alte prestazioni. Ho implementato algoritmi per dataset complessi utilizzando Assembly x86-64 (AVX) e x86-32 (SSE) per massimizzare l'efficienza hardware. Integrazione di OpenMP per il parallelismo e gestione di metriche statistiche avanzate.

# How to run
Software required to run:
 - NASM (www.nasm.us)
 - GCC (gcc.gnu.org)
Both are disponibili as software package installs with the packaging tool of the operative system; for example, on Ubuntu, to run the following command
 - `sudo apt-get install nasm`
 - `sudo apt-get install gcc`

you may need to install the followinf library:
 - `sudo apt-get install lib64gcc-4.8-dev` (or other version)
 - `sudo apt-get install libc6-dev-i386`
 -

To fenerate the excutable file:
 - `nasm -f elf64 cfs64.nasm && gcc -m64 -msse -O0 -no-pie sseutils64.o cfs64.o cfs64c.c -o file64c -lm && ./file64c $pars`
 - `nasm -f elf32 cfs32.nasm && gcc -m32 -msse -O0 -no-pie sseutils32.o cfs32.o cfs32c.c -o file32c -lm && ./file32c $pars`
