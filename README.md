# Feature-Selection-algorithm
This project involves the development of a Correlation Feature Selection (CFS) system in C, optimized for high-performance computing. I and my colleagues implemented algorithms for large database using Assembly x86-64 (AVX) and x86-32 (SSE) to maximize hardware efficiency. Additionally, we integrated OpenMP for parallelism and the management of advanced statistical metrics

# How to run
**Prerequisites**
To run this project, you need the following software:
 - [NASM](www.nasm.us)
 - [GCC](gcc.gnu.org)
 - 
Both are available as packages through your operating system's package manager. For example, on Ubuntu, run: `sudo apt-get install nasm` or/and `sudo apt-get install gcc`

you may also need to install the following compatibility libraries:
 - `sudo apt-get install lib64gcc-4.8-dev` (or other version)
 - `sudo apt-get install libc6-dev-i386`

**Generating the Executable**
To compile and run the 64-bit version (AVX): `nasm -f elf64 cfs64.nasm && gcc -m64 -msse -O0 -no-pie sseutils64.o cfs64.o cfs64c.c -o file64c -lm && ./file64c $pars`

To compile and run the 32-bit version (SSE): `nasm -f elf32 cfs32.nasm && gcc -m32 -msse -O0 -no-pie sseutils32.o cfs32.o cfs32c.c -o file32c -lm && ./file32c $pars`
