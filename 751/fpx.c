/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: core functions over GF(p) and GF(p^2)
*********************************************************************************************/
#include <arm_neon.h>
#include <stdio.h>
//test

__inline void fpcopy(const felm_t a, felm_t c)
{ // Copy a field element, c = a.
    unsigned int i;

    for (i = 0; i < NWORDS_FIELD; i++)
        c[i] = a[i];
}


__inline void fpzero(felm_t a)
{ // Zero a field element, a = 0.
    unsigned int i;

    for (i = 0; i < NWORDS_FIELD; i++)
        a[i] = 0;
}


void to_mont(const felm_t a, felm_t mc)
{ // Conversion to Montgomery representation,
  // mc = a*R^2*R^(-1) mod p = a*R mod p, where a in [0, p-1].
  // The Montgomery constant R^2 mod p is the global value "Montgomery_R2". 

    fpmul_mont(a, (digit_t*)&Montgomery_R2, mc);
}


void from_mont(const felm_t ma, felm_t c)
{ // Conversion from Montgomery representation to standard representation,
  // c = ma*R^(-1) mod p = a mod p, where ma in [0, p-1].

    digit_t one[NWORDS_FIELD] = {0};
    
    one[0] = 1;



    fpmul_mont(ma, one, c);



    fpcorrection(c);
}


void copy_words(const digit_t* a, digit_t* c, const unsigned int nwords)
{ // Copy wordsize digits, c = a, where lng(a) = nwords.
    unsigned int i;
        
    for (i = 0; i < nwords; i++) {                      
        c[i] = a[i];
    }
}


void MUL512(const digit_t* a, const digit_t* b, digit_t* c){
		asm(
			"push  {r0-r11,lr}			\n\t"
			"vpush {q4-q7}				\n\t"
			
			"sub sp, #4 * 2 * 12		\n\t"	//384 * 2 (B384, C384) == +24
			///////////////////////////////////////////////////////////////
			
			//1 2 3 4 5 6 7 8
			//9 10 11 12
			//768 / 32 = 24
			//1-1
			"ldr r1, [r0, #4 * 0]			\n\t"//A0
			"ldr r2, [r0, #4 * 1]			\n\t"//A1
			"ldr r3, [r0, #4 * 2]			\n\t"//A2
			"ldr r4, [r0, #4 * 3]			\n\t"//A3
			"ldr r5, [r0, #4 * 4]			\n\t"//A3
			"ldr r6, [r0, #4 * 5]			\n\t"//A3

			"ldr r7, [r0, #4 * 12]			\n\t"//A12
			"ldr r8, [r0, #4 * 13]			\n\t"//A13
			"ldr r9, [r0, #4 * 14]			\n\t"//A14
			"ldr r10, [r0, #4 * 15]			\n\t"//A15
			"ldr r11, [r0, #4 * 16]			\n\t"//A15
			"ldr r12, [r0, #4 * 17]			\n\t"//A15

			"subs r1, r1, r7				\n\t"
			"sbcs r2, r2, r8				\n\t"
			"sbcs r3, r3, r9				\n\t"
			"sbcs r4, r4, r10				\n\t"
			"sbcs r5, r5, r11				\n\t"
			"sbcs r6, r6, r12				\n\t"

			"vmov d0, r1, r2 				\n\t"
			"vmov d1, r3, r4 				\n\t"
			"vmov d2, r5, r6 				\n\t"

			"ldr r1, [r0, #4 * 6]			\n\t"//A0
			"ldr r2, [r0, #4 * 7]			\n\t"//A1
			"ldr r3, [r0, #4 * 8]			\n\t"//A2
			"ldr r4, [r0, #4 * 9]			\n\t"//A3
			"ldr r5, [r0, #4 * 10]			\n\t"//A3
			"ldr r6, [r0, #4 * 11]			\n\t"//A3

			"ldr r7, [r0, #4 * 18]			\n\t"//A12
			"ldr r8, [r0, #4 * 19]			\n\t"//A13
			"ldr r9, [r0, #4 * 20]			\n\t"//A14
			"ldr r10, [r0, #4 * 21]			\n\t"//A15
			"ldr r11, [r0, #4 * 22]			\n\t"//A15
			"ldr r12, [r0, #4 * 23]			\n\t"//A15

			"sbcs r1, r1, r7			\n\t"
			"sbcs r2, r2, r8			\n\t"
			"sbcs r3, r3, r9			\n\t"
			"sbcs r4, r4, r10			\n\t"
			"sbcs r5, r5, r11			\n\t"
			"sbcs r6, r6, r12			\n\t"

			"sbcs r14, r14, r14			\n\t"

			"vmov d3, r1, r2 			\n\t"
			"vmov d4, r3, r4 			\n\t"
			"vmov d5, r5, r6 			\n\t"

			//eor
			"vmov r1, r2, d0 			\n\t"
			"vmov r3, r4, d1 			\n\t"
			"vmov r5, r6, d2 			\n\t"

			"vmov r7, r8, d3 			\n\t"
			"vmov r9, r10, d4 			\n\t"
			"vmov r11, r12, d5 			\n\t"

			"eor r1, r1, r14			\n\t"
			"eor r2, r2, r14			\n\t"
			"eor r3, r3, r14			\n\t"
			"eor r4, r4, r14			\n\t"
			"eor r5, r5, r14			\n\t"
			"eor r6, r6, r14			\n\t"

			"eor r7, r7, r14			\n\t"
			"eor r8, r8, r14			\n\t"
			"eor r9, r9, r14			\n\t"
			"eor r10, r10, r14			\n\t"
			"eor r11, r11, r14			\n\t"
			"eor r12, r12, r14			\n\t"
			

			"adds r1, r1, r14, LSR #31	\n\t"
			"adcs r2, r2, #0			\n\t"
			"adcs r3, r3, #0			\n\t"
			"adcs r4, r4, #0			\n\t"
			"adcs r5, r5, #0			\n\t"
			"adcs r6, r6, #0			\n\t"

			"adcs r7, r7, #0			\n\t"
			"adcs r8, r8, #0			\n\t"
			"adcs r9, r9, #0			\n\t"
			"adcs r10, r10, #0			\n\t"
			"adcs r11, r11, #0			\n\t"
			"adcs r12, r12, #0			\n\t"


			"vmov d0, r1, r2 			\n\t"
			"vmov d1, r3, r4 			\n\t"
			"vmov d2, r5, r6 			\n\t"			
			"vmov d3, r7, r8 			\n\t"
			"vmov d4, r9, r10 			\n\t"
			"vmov d5, r11, r12 			\n\t"


			//carry bit: r14
			"vmov d12[0], r14			\n\t"

			//second round
			"ldr r0, [sp, #4 * 41]		\n\t"//loading r1

			//1-1
			"ldr r1, [r0, #4 * 0]			\n\t"//A0
			"ldr r2, [r0, #4 * 1]			\n\t"//A1
			"ldr r3, [r0, #4 * 2]			\n\t"//A2
			"ldr r4, [r0, #4 * 3]			\n\t"//A3
			"ldr r5, [r0, #4 * 4]			\n\t"//A3
			"ldr r6, [r0, #4 * 5]			\n\t"//A3

			"ldr r7, [r0, #4 * 12]			\n\t"//A12
			"ldr r8, [r0, #4 * 13]			\n\t"//A13
			"ldr r9, [r0, #4 * 14]			\n\t"//A14
			"ldr r10, [r0, #4 * 15]			\n\t"//A15
			"ldr r11, [r0, #4 * 16]			\n\t"//A15
			"ldr r12, [r0, #4 * 17]			\n\t"//A15

			"subs r1, r1, r7			\n\t"
			"sbcs r2, r2, r8			\n\t"
			"sbcs r3, r3, r9			\n\t"
			"sbcs r4, r4, r10			\n\t"
			"sbcs r5, r5, r11			\n\t"
			"sbcs r6, r6, r12			\n\t"

			"vmov d6, r1, r2 			\n\t"
			"vmov d7, r3, r4 			\n\t"
			"vmov d8, r5, r6 			\n\t"

			"ldr r1, [r0, #4 * 6]			\n\t"//A0
			"ldr r2, [r0, #4 * 7]			\n\t"//A1
			"ldr r3, [r0, #4 * 8]			\n\t"//A2
			"ldr r4, [r0, #4 * 9]			\n\t"//A3
			"ldr r5, [r0, #4 * 10]			\n\t"//A3
			"ldr r6, [r0, #4 * 11]			\n\t"//A3

			"ldr r7, [r0, #4 * 18]			\n\t"//A12
			"ldr r8, [r0, #4 * 19]			\n\t"//A13
			"ldr r9, [r0, #4 * 20]			\n\t"//A14
			"ldr r10, [r0, #4 * 21]			\n\t"//A15
			"ldr r11, [r0, #4 * 22]			\n\t"//A15
			"ldr r12, [r0, #4 * 23]			\n\t"//A15

			"sbcs r1, r1, r7			\n\t"
			"sbcs r2, r2, r8			\n\t"
			"sbcs r3, r3, r9			\n\t"
			"sbcs r4, r4, r10			\n\t"
			"sbcs r5, r5, r11			\n\t"
			"sbcs r6, r6, r12			\n\t"

			"sbcs r14, r14, r14			\n\t"

			"vmov d9, r1, r2 			\n\t"
			"vmov d10, r3, r4 			\n\t"
			"vmov d11, r5, r6 			\n\t"

			//eor
			"vmov r1, r2, d6 			\n\t"
			"vmov r3, r4, d7 			\n\t"
			"vmov r5, r6, d8 			\n\t"

			"vmov r7, r8, d9 			\n\t"
			"vmov r9, r10, d10 			\n\t"
			"vmov r11, r12, d11			\n\t"

			"eor r1, r1, r14			\n\t"
			"eor r2, r2, r14			\n\t"
			"eor r3, r3, r14			\n\t"
			"eor r4, r4, r14			\n\t"
			"eor r5, r5, r14			\n\t"
			"eor r6, r6, r14			\n\t"

			"eor r7, r7, r14			\n\t"
			"eor r8, r8, r14			\n\t"
			"eor r9, r9, r14			\n\t"
			"eor r10, r10, r14			\n\t"
			"eor r11, r11, r14			\n\t"
			"eor r12, r12, r14			\n\t"
			

			"adds r1, r1, r14, LSR #31	\n\t"
			"adcs r2, r2, #0			\n\t"
			"adcs r3, r3, #0			\n\t"
			"adcs r4, r4, #0			\n\t"
			"adcs r5, r5, #0			\n\t"
			"adcs r6, r6, #0			\n\t"

			"adcs r7, r7, #0			\n\t"
			"adcs r8, r8, #0			\n\t"
			"adcs r9, r9, #0			\n\t"
			"adcs r10, r10, #0			\n\t"
			"adcs r11, r11, #0			\n\t"
			"adcs r12, r12, #0			\n\t"


			"str r1, [sp, #4 * 0]			\n\t"//loading result
			"str r2, [sp, #4 * 1]			\n\t"//loading result
			"str r3, [sp, #4 * 2]			\n\t"//loading result
			"str r4, [sp, #4 * 3]			\n\t"//loading result
			"str r5, [sp, #4 * 4]			\n\t"//loading result
			"str r6, [sp, #4 * 5]			\n\t"//loading result

			"str r7, [sp, #4 * 6]			\n\t"//loading result
			"str r8, [sp, #4 * 7]			\n\t"//loading result
			"str r9, [sp, #4 * 8]			\n\t"//loading result
			"str r10, [sp, #4 * 9]			\n\t"//loading result
			"str r11, [sp, #4 * 10]			\n\t"//loading result
			"str r12, [sp, #4 * 11]			\n\t"//loading result





/*
			///test


			"ldr r2, [sp, #4 * 42]			\n\t"//loading result

			"vldr.64 d0, [sp, #4*0]			\n\t"
			"vldr.64 d1, [sp, #4*2]			\n\t"
			"vldr.64 d2, [sp, #4*4]			\n\t"
			"vldr.64 d3, [sp, #4*6]			\n\t"
			"vldr.64 d4, [sp, #4*8]			\n\t"
			"vldr.64 d5, [sp, #4*10]			\n\t"

			"vstr.64 d0, [r2, #4*0]			\n\t"
			"vstr.64 d1, [r2, #4*2]			\n\t"
			"vstr.64 d2, [r2, #4*4]			\n\t"
			"vstr.64 d3, [r2, #4*6]			\n\t"
			"vstr.64 d4, [r2, #4*8]			\n\t"
			"vstr.64 d5, [r2, #4*10]			\n\t"
*/		



			//carry bit: r12
			"vmov r12, d12[0]			\n\t"

			//carry bit: r14
			"eor r12, r12, r14			\n\t"//carry integration
			"mvn r12, r12				\n\t"

			"vdup.32 d7, r12	 		\n\t"//duplication

			////////////////////////////////
			"ldr r2, [sp, #4 * 42]			\n\t"//loading result
			"ldr r0, [sp, #4 * 40]			\n\t"//loading op1
			"ldr r1, [sp, #4 * 41]			\n\t"//loading op2			
			////////////////////////////////




			///

			//NEON
			
			"ldr r3, [r1, #4 * 9]			\n\t"//A9
			"ldr r4, [r1, #4 * 10]			\n\t"//A10
			"vtrn.32 d0, d3 				\n\t"
			"ldr r5, [r1, #4 * 11]			\n\t"//A11
			"ldr r6, [r0, #4 * 0]			\n\t"//B0
			"vtrn.32 d1, d4 				\n\t"
			"ldr r7, [r0, #4 * 1]			\n\t"//B1
			"ldr r8, [r0, #4 * 2]			\n\t"//B2
			"vtrn.32 d2, d5 				\n\t"
			"umull r9, r10, r3, r6			\n\t"//90
			"str r9, [r2, #4 * 9]			\n\t"//C9
			"vldr.64 d6, [sp, #4*0]			\n\t"
			"mov r9, #0						\n\t"
			"umaal r9, r10, r4, r6			\n\t"//100
			"vmull.u32 q15, d0, d6[0] 		\n\t"	//
			"mov r11, #0					\n\t"
			"umaal r9, r11, r3, r7			\n\t"//91
			"vmull.u32 q14, d3, d6[0] 		\n\t"	//
			"str r9, [r2, #4 * 10]			\n\t"//C10
			"mov r12, #0					\n\t"
			"vmull.u32 q13, d1, d6[0] 		\n\t"	//
			"umaal r12, r10, r5, r6			\n\t"//110
			"mov r9, #0						\n\t"
			"vmull.u32 q12, d4, d6[0] 		\n\t"	//
			"umaal r9, r11, r4, r7			\n\t"//101
			"umaal r12, r9, r3, r8			\n\t"//92
			"vmull.u32 q11, d2, d6[0] 		\n\t"	//
			"str r12, [r2, #4 * 11]			\n\t"//C11
			"umaal r9, r10, r5, r7			\n\t"//111
			"vmull.u32 q10, d5, d6[0] 		\n\t"	//
			"umaal r9, r11, r4, r8			\n\t"//102
			"str r9, [r2, #4 * 24]			\n\t"//C12
			"veor q4,q4,q4	 				\n\t"
			"umaal r10, r11, r5, r8			\n\t"//112
			"str r10, [r2, #4 * 25]			\n\t"//C13
			"veor q5,q5,q5	 				\n\t"
			"str r11, [r2, #4 * 26]			\n\t"//C14
			"ldr r3, [r1, #4 * 0]			\n\t"//A0
			"veor q6,q6,q6	 				\n\t"
			"ldr r4, [r1, #4 * 1]			\n\t"//A1
			"ldr r5, [r1, #4 * 2]			\n\t"//A2
			"veor q7,q7,q7	 				\n\t"
			"umull r9, r10, r3, r6			\n\t"//00
			"str r9, [r2, #4 * 0]			\n\t"//C0
			"veor q8,q8,q8	 				\n\t"
			"mov r9, #0					\n\t"
			"umaal r9, r10, r4, r6			\n\t"//10
			"veor q9,q9,q9	 				\n\t"
			"mov r11, #0					\n\t"
			"umaal r9, r11, r3, r7			\n\t"//01
			"vtrn.32 q15, q4 				\n\t"
			"str r9, [r2, #4 * 1]			\n\t"//C1
			"mov r12, #0					\n\t"
			"vtrn.32 q14, q5 				\n\t"
			"umaal r12, r10, r5, r6			\n\t"//20
			"mov r9, #0					\n\t"
			"vtrn.32 q13, q6 				\n\t"
			"umaal r9, r11, r4, r7			\n\t"//11
			"umaal r12, r9, r3, r8			\n\t"//02
			"vtrn.32 q12, q7 				\n\t"
			"str r12, [r2, #4 * 2]			\n\t"//C2
			"mov r12, #0					\n\t"
			"vtrn.32 q11, q8 				\n\t"
			"ldr r6, [r0, #4 * 3]			\n\t"//B3
			"umaal r9, r10, r5, r7			\n\t"//21
			"vtrn.32 q10, q9 				\n\t"
			"umaal r9, r11, r4, r8			\n\t"//12
			"umaal r9, r12, r3, r6			\n\t"//03
			"vadd.i64  q14, q14, q4 		\n\t"	//q
			"str r9, [r2, #4 * 3]			\n\t"//C3
			"mov r9, #0					\n\t"
			"vadd.i64  q13, q13, q5 		\n\t"	//q
			"ldr r7, [r0, #4 * 4]			\n\t"//B4
			"umaal r9, r10, r5, r8			\n\t"//22
			"vadd.i64  q12, q12, q6 		\n\t"	//q
			"umaal r9, r11, r4, r6			\n\t"//13
			"umaal r9, r12, r3, r7			\n\t"//04
			"vadd.i64  q11, q11, q7 		\n\t"	//q
			"str r9, [r2, #4 * 4]			\n\t"//C4
			"mov r9, #0					\n\t"			
			"vadd.i64  q10, q10, q8 		\n\t"	//q
			"ldr r8, [r0, #4 * 5]			\n\t"//B5
			"umaal r9, r10, r5, r6			\n\t"//23
			"vqadd.u64 d18, d18, d31 		\n\t"	//d
			"umaal r9, r11, r4, r7			\n\t"//14
			"umaal r9, r12, r3, r8			\n\t"//05
			"vext.8 d6, d6, d30, #4 		\n\t"
			"str r9, [r2, #4 * 5]			\n\t"//C5
			"mov r9, #0					\n\t"			
			"vmlal.u32 q14, d0, d6[0] 		\n\t"	//
			"ldr r6, [r0, #4 * 6]			\n\t"//B6
			"umaal r9, r10, r5, r7			\n\t"//24
			"vmlal.u32 q13, d3, d6[0] 		\n\t"	//
			"umaal r9, r11, r4, r8			\n\t"//15
			"umaal r9, r12, r3, r6			\n\t"//06
			"vmlal.u32 q12, d1, d6[0] 		\n\t"	//
			"str r9, [r2, #4 * 6]			\n\t"//C6
			"mov r9, #0					\n\t"			
			"vmlal.u32 q11, d4, d6[0] 		\n\t"	//
			"ldr r7, [r0, #4 * 7]			\n\t"//B7
			"umaal r9, r10, r5, r8			\n\t"//25
			"vmlal.u32 q10, d2, d6[0] 		\n\t"	//
			"umaal r9, r11, r4, r6			\n\t"//16
			"umaal r9, r12, r3, r7			\n\t"//07
			"vmlal.u32 q9 , d5, d6[0] 		\n\t"	//
			"str r9, [r2, #4 * 7]			\n\t"//C7
			"mov r9, #0					\n\t"			
			"veor q4,q4,q4	 				\n\t"
			"ldr r8, [r0, #4 * 8]			\n\t"//B8
			"umaal r9, r10, r5, r6			\n\t"//26
			"veor q5,q5,q5	 				\n\t"
			"umaal r9, r11, r4, r7			\n\t"//17
			"umaal r9, r12, r3, r8			\n\t"//08
			"veor q6,q6,q6	 				\n\t"
			"str r9, [r2, #4 * 8]			\n\t"//C8
			"ldr r9, [r2, #4 * 9]			\n\t"//C9
			"veor q7,q7,q7	 				\n\t"
			"ldr r6, [r0, #4 * 9]			\n\t"//B9
			"umaal r9, r10, r5, r7			\n\t"//27
			"veor q8,q8,q8	 				\n\t"
			"umaal r9, r11, r4, r8			\n\t"//18
			"umaal r9, r12, r3, r6			\n\t"//09
			"veor q15,q15,q15 				\n\t"
			"str r9, [r2, #4 * 9]			\n\t"//C9
			"ldr r9, [r2, #4 * 10]			\n\t"//C10
			"vtrn.32 q14, q4 				\n\t"
			"ldr r7, [r0, #4 * 10]			\n\t"//B10
			"umaal r9, r10, r5, r8			\n\t"//28
			"vtrn.32 q13, q5 				\n\t"
			"umaal r9, r11, r4, r6			\n\t"//19
			"umaal r9, r12, r3, r7			\n\t"//010
			"vtrn.32 q12, q6 				\n\t"
			"str r9, [r2, #4 * 10]			\n\t"//C10
			"ldr r9, [r2, #4 * 11]			\n\t"//C11
			"vtrn.32 q11, q7 				\n\t"
			"ldr r8, [r0, #4 * 11]			\n\t"//B11
			"umaal r9, r10, r5, r6			\n\t"//29
			"vtrn.32 q10, q8 				\n\t"
			"umaal r9, r11, r4, r7			\n\t"//110
			"umaal r9, r12, r3, r8			\n\t"//011
			"vtrn.32 q9 , q15 				\n\t"
			"str r9, [r2, #4 * 11]			\n\t"//C11
			"ldr r9, [r2, #4 * 24]			\n\t"//C12
			"vadd.i64  q13, q13, q4 		\n\t"	//q
			"ldr r3, [r1, #4 * 3]			\n\t"//A3
			"umaal r9, r10, r3, r6			\n\t"//39
			"vadd.i64  q12, q12, q5 		\n\t"	//q
			"umaal r9, r11, r5, r7			\n\t"//210
			"umaal r9, r12, r4, r8			\n\t"//111
			"vadd.i64  q11, q11, q6 		\n\t"	//q
			"str r9, [r2, #4 * 24]			\n\t"//C12
			"ldr r9, [r2, #4 * 25]			\n\t"//C13
			"vadd.i64  q10, q10, q7 		\n\t"	//q
			"ldr r4, [r1, #4 * 4]			\n\t"//A4
			"umaal r9, r10, r4, r6			\n\t"//49
			"vadd.i64  q9 , q9 , q8 		\n\t"	//q
			"umaal r9, r11, r3, r7			\n\t"//310
			"umaal r9, r12, r5, r8			\n\t"//211
			"vqadd.u64 d30, d30, d29 		\n\t"	//d
			"str r9, [r2, #4 * 25]			\n\t"//C13
			"ldr r9, [r2, #4 * 26]			\n\t"//C14
			"vext.8 d6, d6, d28, #4 		\n\t"
			"ldr r5, [r1, #4 * 5]			\n\t"//A5
			"umaal r9, r10, r5, r6			\n\t"//59
			"veor d6, d6, d7				\n\t"
			"umaal r9, r11, r4, r7			\n\t"//410
			"umaal r9, r12, r3, r8			\n\t"//311
			"vstr.64 d6, [sp, #4*12]			\n\t"
			"str r9, [r2, #4 * 26]			\n\t"//C14
			"umaal r10, r11, r5, r7			\n\t"//510
			"vldr.64 d6, [sp, #4*2]			\n\t"
			"umaal r10, r12, r4, r8			\n\t"//411
			"str r10, [r2, #4 * 27]			\n\t"//C15
			"vmlal.u32 q13, d0, d6[0] 		\n\t"	//
			"umaal r11, r12, r5, r8			\n\t"//511
			"str r11, [r2, #4 * 28]			\n\t"//C16
			"vmlal.u32 q12, d3, d6[0] 		\n\t"	//
			"str r12, [r2, #4 * 29]			\n\t"//C17
			"ldr r6, [r0, #4 * 0]			\n\t"//B0
			"vmlal.u32 q11, d1, d6[0] 		\n\t"	//
			"ldr r7, [r0, #4 * 1]			\n\t"//B1
			"ldr r8, [r0, #4 * 2]			\n\t"//B2
			"vmlal.u32 q10, d4, d6[0] 		\n\t"	//
			"ldr r9, [r2, #4 * 3]			\n\t"//C3
			"mov r10, #0					\n\t"
			"vmlal.u32 q9 , d2, d6[0] 		\n\t"	//
			"umaal r9, r10, r3, r6			\n\t"//30
			"str r9, [r2, #4 * 3]			\n\t"//C3
			"vmlal.u32 q15, d5, d6[0] 		\n\t"	//
			"ldr r9, [r2, #4 * 4]			\n\t"//C4
			"mov r11, #0					\n\t"
			"veor q4,q4,q4	 				\n\t"
			"umaal r9, r10, r4, r6			\n\t"//40
			"umaal r9, r11, r3, r7			\n\t"//31
			"veor q5,q5,q5	 				\n\t"
			"str r9, [r2, #4 * 4]			\n\t"//C4
			"ldr r9, [r2, #4 * 5]			\n\t"//C5
			"veor q6,q6,q6	 				\n\t"
			"mov r12, #0					\n\t"
			"umaal r9, r10, r5, r6			\n\t"//50
			"veor q7,q7,q7	 				\n\t"
			"umaal r9, r11, r4, r7			\n\t"//41
			"umaal r9, r12, r3, r8			\n\t"//32
			"veor q8,q8,q8	 				\n\t"
			"str r9, [r2, #4 * 5]			\n\t"//C5
			"ldr r9, [r2, #4 * 6]			\n\t"//C6
			"veor q14,q14,q14 				\n\t"
			"ldr r6, [r0, #4 * 3]			\n\t"//B3
			"umaal r9, r10, r5, r7			\n\t"//51
			"vtrn.32 q13, q4 				\n\t"
			"umaal r9, r11, r4, r8			\n\t"//42
			"umaal r9, r12, r3, r6			\n\t"//33
			"vtrn.32 q12, q5 				\n\t"
			"str r9, [r2, #4 * 6]			\n\t"//C6
			"ldr r9, [r2, #4 * 7]			\n\t"//C7
			"vtrn.32 q11, q6 				\n\t"
			"ldr r7, [r0, #4 * 4]			\n\t"//B4
			"umaal r9, r10, r5, r8			\n\t"//52
			"vtrn.32 q10, q7 				\n\t"
			"umaal r9, r11, r4, r6			\n\t"//43
			"umaal r9, r12, r3, r7			\n\t"//34
			"vtrn.32 q9 , q8 				\n\t"
			"str r9, [r2, #4 * 7]			\n\t"//C7
			"ldr r9, [r2, #4 * 8]			\n\t"//C8
			"vtrn.32 q15, q14 				\n\t"
			"ldr r8, [r0, #4 * 5]			\n\t"//B5
			"umaal r9, r10, r5, r6			\n\t"//53
			"vadd.i64  q12, q12, q4 		\n\t"	//q
			"umaal r9, r11, r4, r7			\n\t"//44
			"umaal r9, r12, r3, r8			\n\t"//35
			"vadd.i64  q11, q11, q5 		\n\t"	//q
			"str r9, [r2, #4 * 8]			\n\t"//C8
			"ldr r9, [r2, #4 * 9]			\n\t"//C9
			"vadd.i64  q10, q10, q6 		\n\t"	//q
			"ldr r6, [r0, #4 * 6]			\n\t"//B6
			"umaal r9, r10, r5, r7			\n\t"//54
			"vadd.i64  q9 , q9 , q7 		\n\t"	//q
			"umaal r9, r11, r4, r8			\n\t"//45
			"umaal r9, r12, r3, r6			\n\t"//36
			"vadd.i64  q15, q15, q8		\n\t"	//q
			"str r9, [r2, #4 * 9]			\n\t"//C9
			"ldr r9, [r2, #4 * 10]			\n\t"//C10
			"vqadd.u64 d28, d28, d27 		\n\t"	//d
			"ldr r7, [r0, #4 * 7]			\n\t"//B7
			"umaal r9, r10, r5, r8			\n\t"//55
			"vext.8 d6, d6, d26, #4 		\n\t"
			"umaal r9, r11, r4, r6			\n\t"//46
			"umaal r9, r12, r3, r7			\n\t"//37
			"vmlal.u32 q12, d0, d6[0] 		\n\t"	//
			"str r9, [r2, #4 * 10]			\n\t"//C10
			"ldr r9, [r2, #4 * 11]			\n\t"//C11
			"vmlal.u32 q11, d3, d6[0] 		\n\t"	//
			"ldr r8, [r0, #4 * 8]			\n\t"//B8
			"umaal r9, r10, r5, r6			\n\t"//56
			"vmlal.u32 q10, d1, d6[0] 		\n\t"	//
			"umaal r9, r11, r4, r7			\n\t"//47
			"umaal r9, r12, r3, r8			\n\t"//38
			"vmlal.u32 q9 , d4, d6[0] 		\n\t"	//
			"str r9, [r2, #4 * 11]			\n\t"//C11
			"ldr r9, [r2, #4 * 24]			\n\t"//C12
			"vmlal.u32 q15, d2, d6[0] 		\n\t"	//
			"ldr r3, [r1, #4 * 6]			\n\t"//A6
			"umaal r9, r10, r3, r6			\n\t"//66
			"vmlal.u32 q14, d5, d6[0] 		\n\t"	//
			"umaal r9, r11, r5, r7			\n\t"//57
			"umaal r9, r12, r4, r8			\n\t"//48
			"veor q4,q4,q4	 				\n\t"
			"str r9, [r2, #4 * 24]			\n\t"//C12
			"ldr r9, [r2, #4 * 25]			\n\t"//C13
			"veor q5,q5,q5	 				\n\t"
			"ldr r4, [r1, #4 * 7]			\n\t"//A7
			"umaal r9, r10, r4, r6			\n\t"//76
			"veor q6,q6,q6	 				\n\t"
			"umaal r9, r11, r3, r7			\n\t"//67
			"umaal r9, r12, r5, r8			\n\t"//58
			"veor q7,q7,q7	 				\n\t"
			"str r9, [r2, #4 * 25]			\n\t"//C13
			"ldr r9, [r2, #4 * 26]			\n\t"//C14
			"veor q8,q8,q8	 				\n\t"
			"ldr r5, [r1, #4 * 8]			\n\t"//A8
			"umaal r9, r10, r5, r6			\n\t"//86
			"veor q13,q13,q13 				\n\t"
			"umaal r9, r11, r4, r7			\n\t"//77
			"umaal r9, r12, r3, r8			\n\t"//68
			"vtrn.32 q12, q4 				\n\t"
			"str r9, [r2, #4 * 26]			\n\t"//C14
			"ldr r9, [r2, #4 * 27]			\n\t"//C15
			"vtrn.32 q11, q5 				\n\t"
			"ldr r6, [r0, #4 * 9]			\n\t"//B9
			"umaal r9, r10, r5, r7			\n\t"//87
			"vtrn.32 q10, q6 				\n\t"
			"umaal r9, r11, r4, r8			\n\t"//78
			"umaal r9, r12, r3, r6			\n\t"//69
			"vtrn.32 q9 , q7 				\n\t"
			"str r9, [r2, #4 * 27]			\n\t"//C15
			"ldr r9, [r2, #4 * 28]			\n\t"//C16
			"vtrn.32 q15, q8 				\n\t"
			"ldr r7, [r0, #4 * 10]			\n\t"//B10
			"umaal r9, r10, r5, r8			\n\t"//88
			"vtrn.32 q14, q13 				\n\t"
			"umaal r9, r11, r4, r6			\n\t"//79
			"umaal r9, r12, r3, r7			\n\t"//610
			"vadd.i64  q11, q11, q4 		\n\t"	//q
			"str r9, [r2, #4 * 28]			\n\t"//C16
			"ldr r9, [r2, #4 * 29]			\n\t"//C17
			"vadd.i64  q10, q10, q5 		\n\t"	//q
			"ldr r8, [r0, #4 * 11]			\n\t"//B11
			"umaal r9, r10, r5, r6			\n\t"//89
			"vadd.i64  q9 , q9 , q6 		\n\t"	//q
			"umaal r9, r11, r4, r7			\n\t"//710
			"umaal r9, r12, r3, r8			\n\t"//611
			"vadd.i64  q15, q15, q7 		\n\t"	//q
			"str r9, [r2, #4 * 29]			\n\t"//C17
			"umaal r10, r11, r5, r7			\n\t"//810
			"vadd.i64  q14, q14, q8 		\n\t"	//q
			"umaal r10, r12, r4, r8			\n\t"//711
			"str r10, [r2, #4 * 30]			\n\t"//C18
			"vqadd.u64 d26, d26, d25 		\n\t"	//d
			"umaal r11, r12, r5, r8			\n\t"//811
			"str r11, [r2, #4 * 31]			\n\t"//C19
			"vext.8 d6, d6, d24, #4 		\n\t"
			"str r12, [r2, #4 * 32]			\n\t"//C20
			"ldr r6, [r0, #4 * 0]			\n\t"//B0
			"veor d6, d6, d7				\n\t"
			"ldr r7, [r0, #4 * 1]			\n\t"//B1
			"ldr r8, [r0, #4 * 2]			\n\t"//B2
			"vstr.64 d6, [sp, #4*14]			\n\t"
			"ldr r9, [r2, #4 * 6]			\n\t"//C6
			"mov r10, #0					\n\t"
			"vldr.64 d6, [sp, #4*4]			\n\t"
			"umaal r9, r10, r3, r6			\n\t"//60
			"str r9, [r2, #4 * 6]			\n\t"//C6
			"vmlal.u32 q11, d0, d6[0] 		\n\t"	//
			"ldr r9, [r2, #4 * 7]			\n\t"//C7
			"mov r11, #0					\n\t"
			"vmlal.u32 q10, d3, d6[0] 		\n\t"	//
			"umaal r9, r10, r4, r6			\n\t"//70
			"umaal r9, r11, r3, r7			\n\t"//61
			"vmlal.u32 q9 , d1, d6[0] 		\n\t"	//
			"str r9, [r2, #4 * 7]			\n\t"//C7
			"ldr r9, [r2, #4 * 8]			\n\t"//C8
			"vmlal.u32 q15, d4, d6[0] 		\n\t"	//
			"mov r12, #0					\n\t"
			"umaal r9, r10, r5, r6			\n\t"//80
			"vmlal.u32 q14, d2, d6[0] 		\n\t"	//
			"umaal r9, r11, r4, r7			\n\t"//71
			"umaal r9, r12, r3, r8			\n\t"//62
			"vmlal.u32 q13, d5, d6[0] 		\n\t"	//
			"str r9, [r2, #4 * 8]			\n\t"//C8
			"ldr r9, [r2, #4 * 9]			\n\t"//C9
			"veor q4,q4,q4	 				\n\t"
			"ldr r6, [r0, #4 * 3]			\n\t"//B3
			"umaal r9, r10, r5, r7			\n\t"//81
			"veor q5,q5,q5	 				\n\t"
			"umaal r9, r11, r4, r8			\n\t"//72
			"umaal r9, r12, r3, r6			\n\t"//63
			"veor q6,q6,q6	 				\n\t"
			"str r9, [r2, #4 * 9]			\n\t"//C9
			"ldr r9, [r2, #4 * 10]			\n\t"//C10
			"veor q7,q7,q7	 				\n\t"
			"ldr r7, [r0, #4 * 4]			\n\t"//B4
			"umaal r9, r10, r5, r8			\n\t"//82
			"veor q8,q8,q8	 				\n\t"
			"umaal r9, r11, r4, r6			\n\t"//73
			"umaal r9, r12, r3, r7			\n\t"//64
			"veor q12,q12,q12 				\n\t"
			"str r9, [r2, #4 * 10]			\n\t"//C10
			"ldr r9, [r2, #4 * 11]			\n\t"//C11
			"vtrn.32 q11, q4 				\n\t"
			"ldr r8, [r0, #4 * 5]			\n\t"//B5
			"umaal r9, r10, r5, r6			\n\t"//83
			"vtrn.32 q10, q5 				\n\t"
			"umaal r9, r11, r4, r7			\n\t"//74
			"umaal r9, r12, r3, r8			\n\t"//65
			"vtrn.32 q9 , q6 				\n\t"
			"str r9, [r2, #4 * 11]			\n\t"//C11
			"ldr r9, [r2, #4 * 24]			\n\t"//C12
			"vtrn.32 q15, q7 				\n\t"
			"ldr r3, [r1, #4 * 9]			\n\t"//A9
			"umaal r9, r10, r3, r6			\n\t"//93
			"vtrn.32 q14, q8 				\n\t"
			"umaal r9, r11, r5, r7			\n\t"//84
			"umaal r9, r12, r4, r8			\n\t"//75
			"vtrn.32 q13, q12 				\n\t"
			"str r9, [r2, #4 * 24]			\n\t"//C12
			"ldr r9, [r2, #4 * 25]			\n\t"//C13
			"vadd.i64  q10, q10, q4 		\n\t"	//q
			"ldr r4, [r1, #4 * 10]			\n\t"//A10
			"umaal r9, r10, r4, r6			\n\t"//103
			"vadd.i64  q9 , q9 , q5 		\n\t"	//q
			"umaal r9, r11, r3, r7			\n\t"//94
			"umaal r9, r12, r5, r8			\n\t"//85
			"vadd.i64  q15, q15, q6 		\n\t"	//q
			"str r9, [r2, #4 * 25]			\n\t"//C13			//4-8
			"ldr r9, [r2, #4 * 26]			\n\t"//C14
			"vadd.i64  q14, q14, q7 		\n\t"	//q
			"ldr r5, [r1, #4 * 11]			\n\t"//A11
			"umaal r9, r10, r5, r6			\n\t"//113
			"vadd.i64  q13, q13, q8 		\n\t"	//q
			"umaal r9, r11, r4, r7			\n\t"//104
			"umaal r9, r12, r3, r8			\n\t"//95
			"vqadd.u64 d24, d24, d23 		\n\t"	//d
			"str r9, [r2, #4 * 26]			\n\t"//C14
			"ldr r9, [r2, #4 * 27]			\n\t"//C15
			"vext.8 d6, d6, d22, #4 		\n\t"
			"ldr r6, [r0, #4 * 6]			\n\t"//B6
			"umaal r9, r10, r5, r7			\n\t"//114
			"vmlal.u32 q10, d0, d6[0] 		\n\t"	//
			"umaal r9, r11, r4, r8			\n\t"//105
			"umaal r9, r12, r3, r6			\n\t"//96
			"vmlal.u32 q9 , d3, d6[0] 		\n\t"	//
			"str r9, [r2, #4 * 27]			\n\t"//C15
			"ldr r9, [r2, #4 * 28]			\n\t"//C16
			"vmlal.u32 q15, d1, d6[0] 		\n\t"	//
			"ldr r7, [r0, #4 * 7]			\n\t"//B7
			"umaal r9, r10, r5, r8			\n\t"//115
			"vmlal.u32 q14, d4, d6[0] 		\n\t"	//
			"umaal r9, r11, r4, r6			\n\t"//106
			"umaal r9, r12, r3, r7			\n\t"//97
			"vmlal.u32 q13, d2, d6[0] 		\n\t"	//
			"str r9, [r2, #4 * 28]			\n\t"//C16
			"ldr r9, [r2, #4 * 29]			\n\t"//C17
			"vmlal.u32 q12, d5, d6[0] 		\n\t"	//
			"ldr r8, [r0, #4 * 8]			\n\t"//B8
			"umaal r9, r10, r5, r6			\n\t"//116
			"veor q4,q4,q4	 				\n\t"
			"umaal r9, r11, r4, r7			\n\t"//107
			"umaal r9, r12, r3, r8			\n\t"//98
			"veor q5,q5,q5	 				\n\t"
			"str r9, [r2, #4 * 29]			\n\t"//C17
			"ldr r9, [r2, #4 * 30]			\n\t"//C18
			"veor q6,q6,q6	 				\n\t"
			"ldr r6, [r0, #4 * 9]			\n\t"//B9
			"umaal r9, r10, r5, r7			\n\t"//117
			"veor q7,q7,q7	 				\n\t"
			"umaal r9, r11, r4, r8			\n\t"//108
			"umaal r9, r12, r3, r6			\n\t"//99
			"veor q8,q8,q8	 				\n\t"
			"str r9, [r2, #4 * 30]			\n\t"//C18
			"ldr r9, [r2, #4 * 31]			\n\t"//C19
			"veor q11,q11,q11 				\n\t"
			"ldr r7, [r0, #4 * 10]			\n\t"//B10
			"umaal r9, r10, r5, r8			\n\t"//118
			"vtrn.32 q10, q4 				\n\t"
			"umaal r9, r11, r4, r6			\n\t"//109
			"umaal r9, r12, r3, r7			\n\t"//910
			"vtrn.32 q9 , q5 				\n\t"
			"str r9, [r2, #4 * 31]			\n\t"//C19
			"ldr r9, [r2, #4 * 32]			\n\t"//C20
			"vtrn.32 q15, q6 				\n\t"
			"ldr r8, [r0, #4 * 11]			\n\t"//B11
			"umaal r9, r10, r5, r6			\n\t"//119
			"vtrn.32 q14, q7 				\n\t"
			"umaal r9, r11, r4, r7			\n\t"//1010
			"umaal r9, r12, r3, r8			\n\t"//9111
			"vtrn.32 q13, q8 				\n\t"
			"str r9, [r2, #4 * 32]			\n\t"//C20
			"umaal r10, r11, r5, r7			\n\t"//1110
			"vtrn.32 q12, q11 				\n\t"
			"umaal r10, r12, r4, r8			\n\t"//1011
			"str r10, [r2, #4 * 33]			\n\t"//C21
			"vadd.i64  q9 , q9 , q4 		\n\t"	//q
			"umaal r11, r12, r5, r8			\n\t"//1111
			"str r11, [r2, #4 * 34]			\n\t"//C22
			"vadd.i64  q15, q15, q5 		\n\t"	//q
			"str r12, [r2, #4 * 35]			\n\t"//C23
			"ldr r3, [r1, #4 * 21]			\n\t"//A9
			"vadd.i64  q14, q14, q6 		\n\t"	//q
			"ldr r4, [r1, #4 * 22]			\n\t"//A10
			"ldr r5, [r1, #4 * 23]			\n\t"//A11
			"vadd.i64  q13, q13, q7 		\n\t"	//q
			"ldr r6, [r0, #4 * 12]			\n\t"//B0
			"ldr r7, [r0, #4 * 13]			\n\t"//B1
			"vadd.i64  q12, q12, q8 		\n\t"	//q
			"ldr r8, [r0, #4 * 14]			\n\t"//B2
			"ldr r9, [r2, #4 * 33]			\n\t"//C9
			"vqadd.u64 d22, d22, d21 		\n\t"	//d
			"mov r10, #0					\n\t"
			"umaal r9, r10, r3, r6			\n\t"//90
			"vext.8 d6, d6, d20, #4 		\n\t"
			"str r9, [r2, #4 * 33]			\n\t"//C9
			"ldr r9, [r2, #4 * 34]			\n\t"//C10
			"veor d6, d6, d7				\n\t"
			"umaal r9, r10, r4, r6			\n\t"//100
			"mov r11, #0					\n\t"
			"vstr.64 d6, [sp, #4*16]			\n\t"
			"umaal r9, r11, r3, r7			\n\t"//91
			"str r9, [r2, #4 * 34]			\n\t"//C10
			"vldr.64 d6, [sp, #4*6]			\n\t"
			"ldr r12, [r2, #4 * 35]			\n\t"//C11
			"umaal r12, r10, r5, r6			\n\t"//110
			"vmlal.u32 q9 , d0, d6[0] 		\n\t"	//
			"mov r9, #0						\n\t"
			"umaal r9, r11, r4, r7			\n\t"//101
			"vmlal.u32 q15, d3, d6[0] 		\n\t"	//
			"umaal r12, r9, r3, r8			\n\t"//92
			"str r12, [r2, #4 * 35]			\n\t"//C11
			"vmlal.u32 q14, d1, d6[0] 		\n\t"	//
			"umaal r9, r10, r5, r7			\n\t"//111
			"umaal r9, r11, r4, r8			\n\t"//102
			"vmlal.u32 q13, d4, d6[0] 		\n\t"	//
			"str r9, [r2, #4 * 36]			\n\t"//C12
			"umaal r10, r11, r5, r8			\n\t"//112
			"vmlal.u32 q12, d2, d6[0] 		\n\t"	//
			"str r10, [r2, #4 * 37]			\n\t"//C13
			"str r11, [r2, #4 * 38]			\n\t"//C14
			"vmlal.u32 q11, d5, d6[0] 		\n\t"	//
			"ldr r3, [r1, #4 * 12]			\n\t"//A0
			"ldr r4, [r1, #4 * 13]			\n\t"//A1
			"veor q4,q4,q4	 				\n\t"
			"ldr r5, [r1, #4 * 14]			\n\t"//A2
			"ldr r9, [r2, #4 * 24]			\n\t"//C0
			"veor q5,q5,q5	 				\n\t"
			"mov r10, #0					\n\t"
			"umaal r9, r10, r3, r6			\n\t"//00
			"veor q6,q6,q6	 				\n\t"
			"str r9, [r2, #4 * 24]			\n\t"//C0
			"ldr r9, [r2, #4 * 25]			\n\t"//C1
			"veor q7,q7,q7	 				\n\t"
			"umaal r9, r10, r4, r6			\n\t"//10
			"mov r11, #0					\n\t"
			"veor q8,q8,q8	 				\n\t"
			"umaal r9, r11, r3, r7			\n\t"//01
			"str r9, [r2, #4 * 25]			\n\t"//C1
			"veor q10,q10,q10 				\n\t"
			"ldr r12, [r2, #4 * 26]			\n\t"//C2
			"umaal r12, r10, r5, r6			\n\t"//20
			"vtrn.32 q9 , q4 				\n\t"
			"mov r9, #0					\n\t"
			"umaal r9, r11, r4, r7			\n\t"//11
			"vtrn.32 q15, q5 				\n\t"
			"umaal r12, r9, r3, r8			\n\t"//02
			"str r12, [r2, #4 * 26]			\n\t"//C2
			"vtrn.32 q14, q6 				\n\t"
			"ldr r12, [r2, #4 * 27]			\n\t"//C3
			"ldr r6, [r0, #4 * 15]			\n\t"//B3
			"vtrn.32 q13, q7 				\n\t"
			"umaal r9, r10, r5, r7			\n\t"//21
			"umaal r9, r11, r4, r8			\n\t"//12
			"vtrn.32 q12, q8 				\n\t"
			"umaal r9, r12, r3, r6			\n\t"//03
			"str r9, [r2, #4 * 27]			\n\t"//C3
			"vtrn.32 q11, q10 				\n\t"	
			"ldr r9, [r2, #4 * 28]			\n\t"//C4
			"ldr r7, [r0, #4 * 16]			\n\t"//B4
			"vadd.i64  q15, q15, q4 		\n\t"	//q
			"umaal r9, r10, r5, r8			\n\t"//22
			"umaal r9, r11, r4, r6			\n\t"//13
			"vadd.i64  q14, q14, q5 		\n\t"	//q
			"umaal r9, r12, r3, r7			\n\t"//04
			"str r9, [r2, #4 * 28]			\n\t"//C4
			"vadd.i64  q13, q13, q6 		\n\t"	//q
			"ldr r9, [r2, #4 * 29]			\n\t"//C5	
			"ldr r8, [r0, #4 * 17]			\n\t"//B5
			"vadd.i64  q12, q12, q7 		\n\t"	//q
			"umaal r9, r10, r5, r6			\n\t"//23
			"umaal r9, r11, r4, r7			\n\t"//14
			"vadd.i64  q11, q11, q8 		\n\t"	//q
			"umaal r9, r12, r3, r8			\n\t"//05
			"str r9, [r2, #4 * 29]			\n\t"//C5
			"vqadd.u64 d20, d20, d19 		\n\t"	//d
			"ldr r9, [r2, #4 * 30]			\n\t"//C6	
			"ldr r6, [r0, #4 * 18]			\n\t"//B6
			"vext.8 d6, d6, d18, #4 		\n\t"
			"umaal r9, r10, r5, r7			\n\t"//24
			"umaal r9, r11, r4, r8			\n\t"//15
			"vmlal.u32 q15, d0, d6[0] 		\n\t"	//
			"umaal r9, r12, r3, r6			\n\t"//06
			"str r9, [r2, #4 * 30]			\n\t"//C6
			"vmlal.u32 q14, d3, d6[0] 		\n\t"	//
			"ldr r9, [r2, #4 * 31]			\n\t"//C7	
			"ldr r7, [r0, #4 * 19]			\n\t"//B7
			"vmlal.u32 q13, d1, d6[0] 		\n\t"	//
			"umaal r9, r10, r5, r8			\n\t"//25
			"umaal r9, r11, r4, r6			\n\t"//16
			"vmlal.u32 q12, d4, d6[0] 		\n\t"	//
			"umaal r9, r12, r3, r7			\n\t"//07
			"str r9, [r2, #4 * 31]			\n\t"//C7
			"vmlal.u32 q11, d2, d6[0] 		\n\t"	//
			"ldr r9, [r2, #4 * 32]			\n\t"//C8	
			"ldr r8, [r0, #4 * 20]			\n\t"//B8
			"vmlal.u32 q10, d5, d6[0] 		\n\t"	//
			"umaal r9, r10, r5, r6			\n\t"//26
			"umaal r9, r11, r4, r7			\n\t"//17
			"veor q4,q4,q4	 				\n\t"
			"umaal r9, r12, r3, r8			\n\t"//08
			"str r9, [r2, #4 * 32]			\n\t"//C8
			"veor q5,q5,q5	 				\n\t"
			"ldr r9, [r2, #4 * 33]			\n\t"//C9
			"ldr r6, [r0, #4 * 21]			\n\t"//B9
			"veor q6,q6,q6	 				\n\t"
			"umaal r9, r10, r5, r7			\n\t"//27
			"umaal r9, r11, r4, r8			\n\t"//18
			"veor q7,q7,q7	 				\n\t"
			"umaal r9, r12, r3, r6			\n\t"//09
			"str r9, [r2, #4 * 33]			\n\t"//C9
			"veor q8,q8,q8	 				\n\t"
			"ldr r9, [r2, #4 * 34]			\n\t"//C10
			"ldr r7, [r0, #4 * 22]			\n\t"//B10
			"veor q9,q9,q9	 				\n\t"
			"umaal r9, r10, r5, r8			\n\t"//28
			"umaal r9, r11, r4, r6			\n\t"//19
			"vtrn.32 q15, q4 				\n\t"
			"umaal r9, r12, r3, r7			\n\t"//010
			"str r9, [r2, #4 * 34]			\n\t"//C10
			"vtrn.32 q14, q5 				\n\t"
			"ldr r9, [r2, #4 * 35]			\n\t"//C11
			"ldr r8, [r0, #4 * 23]			\n\t"//B11
			"vtrn.32 q13, q6 				\n\t"
			"umaal r9, r10, r5, r6			\n\t"//29
			"umaal r9, r11, r4, r7			\n\t"//110
			"vtrn.32 q12, q7 				\n\t"
			"umaal r9, r12, r3, r8			\n\t"//011
			"str r9, [r2, #4 * 35]			\n\t"//C11
			"vtrn.32 q11, q8 				\n\t"
			"ldr r9, [r2, #4 * 36]			\n\t"//C12
			"ldr r3, [r1, #4 * 15]			\n\t"//A3
			"vtrn.32 q10, q9 				\n\t"
			"umaal r9, r10, r3, r6			\n\t"//39
			"umaal r9, r11, r5, r7			\n\t"//210
			"vadd.i64  q14, q14, q4 		\n\t"	//q
			"umaal r9, r12, r4, r8			\n\t"//111
			"str r9, [r2, #4 * 36]			\n\t"//C12
			"vadd.i64  q13, q13, q5 		\n\t"	//q
			"ldr r9, [r2, #4 * 37]			\n\t"//C13
			"ldr r4, [r1, #4 * 16]			\n\t"//A4
			"vadd.i64  q12, q12, q6 		\n\t"	//q
			"umaal r9, r10, r4, r6			\n\t"//49
			"umaal r9, r11, r3, r7			\n\t"//310
			"vadd.i64  q11, q11, q7 		\n\t"	//q
			"umaal r9, r12, r5, r8			\n\t"//211
			"str r9, [r2, #4 * 37]			\n\t"//C13
			"vadd.i64  q10, q10, q8 		\n\t"	//q
			"ldr r9, [r2, #4 * 38]			\n\t"//C14
			"ldr r5, [r1, #4 * 17]			\n\t"//A5
			"vqadd.u64 d18, d18, d31 		\n\t"	//d
			"umaal r9, r10, r5, r6			\n\t"//59
			"umaal r9, r11, r4, r7			\n\t"//410
			"vext.8 d6, d6, d30, #4 		\n\t"
			"umaal r9, r12, r3, r8			\n\t"//311
			"str r9, [r2, #4 * 38]			\n\t"//C14
			"veor d6, d6, d7				\n\t"
			"umaal r10, r11, r5, r7			\n\t"//510
			"umaal r10, r12, r4, r8			\n\t"//411
			"vstr.64 d6, [sp, #4*18]			\n\t"
			"str r10, [r2, #4 * 39]			\n\t"//C15
			"umaal r11, r12, r5, r8			\n\t"//511
			"vldr.64 d6, [sp, #4*8]			\n\t"
			"str r11, [r2, #4 * 40]			\n\t"//C16
			"str r12, [r2, #4 * 41]			\n\t"//C17
			"vmlal.u32 q14, d0, d6[0] 		\n\t"	//
			"ldr r6, [r0, #4 * 12]			\n\t"//B0
			"ldr r7, [r0, #4 * 13]			\n\t"//B1
			"vmlal.u32 q13, d3, d6[0] 		\n\t"	//
			"ldr r8, [r0, #4 * 14]			\n\t"//B2
			"ldr r9, [r2, #4 * 27]			\n\t"//C3
			"vmlal.u32 q12, d1, d6[0] 		\n\t"	//
			"mov r10, #0					\n\t"
			"umaal r9, r10, r3, r6			\n\t"//30
			"vmlal.u32 q11, d4, d6[0] 		\n\t"	//
			"str r9, [r2, #4 * 27]			\n\t"//C3
			"ldr r9, [r2, #4 * 28]			\n\t"//C4
			"vmlal.u32 q10, d2, d6[0] 		\n\t"	//
			"mov r11, #0					\n\t"
			"umaal r9, r10, r4, r6			\n\t"//40
			"vmlal.u32 q9 , d5, d6[0] 		\n\t"	//
			"umaal r9, r11, r3, r7			\n\t"//31
			"str r9, [r2, #4 * 28]			\n\t"//C4
			"veor q4,q4,q4	 				\n\t"
			"ldr r9, [r2, #4 * 29]			\n\t"//C5
			"mov r12, #0					\n\t"
			"veor q5,q5,q5	 				\n\t"
			"umaal r9, r10, r5, r6			\n\t"//50
			"umaal r9, r11, r4, r7			\n\t"//41
			"veor q6,q6,q6	 				\n\t"
			"umaal r9, r12, r3, r8			\n\t"//32
			"str r9, [r2, #4 * 29]			\n\t"//C5
			"veor q7,q7,q7	 				\n\t"
			"ldr r9, [r2, #4 * 30]			\n\t"//C6
			"ldr r6, [r0, #4 * 15]			\n\t"//B3
			"veor q8,q8,q8	 				\n\t"
			"umaal r9, r10, r5, r7			\n\t"//51
			"umaal r9, r11, r4, r8			\n\t"//42
			"veor q15,q15,q15 				\n\t"
			"umaal r9, r12, r3, r6			\n\t"//33
			"str r9, [r2, #4 * 30]			\n\t"//C6
			"vtrn.32 q14, q4 				\n\t"
			"ldr r9, [r2, #4 * 31]			\n\t"//C7
			"ldr r7, [r0, #4 * 16]			\n\t"//B4
			"vtrn.32 q13, q5 				\n\t"
			"umaal r9, r10, r5, r8			\n\t"//52
			"umaal r9, r11, r4, r6			\n\t"//43
			"vtrn.32 q12, q6 				\n\t"
			"umaal r9, r12, r3, r7			\n\t"//34
			"str r9, [r2, #4 * 31]			\n\t"//C7
			"vtrn.32 q11, q7 				\n\t"
			"ldr r9, [r2, #4 * 32]			\n\t"//C8
			"ldr r8, [r0, #4 * 17]			\n\t"//B5
			"vtrn.32 q10, q8 				\n\t"
			"umaal r9, r10, r5, r6			\n\t"//53
			"umaal r9, r11, r4, r7			\n\t"//44
			"vtrn.32 q9 , q15 				\n\t"
			"umaal r9, r12, r3, r8			\n\t"//35
			"str r9, [r2, #4 * 32]			\n\t"//C8
			"vadd.i64  q13, q13, q4 		\n\t"	//q
			"ldr r9, [r2, #4 * 33]			\n\t"//C9
			"ldr r6, [r0, #4 * 18]			\n\t"//B6
			"vadd.i64  q12, q12, q5 		\n\t"	//q
			"umaal r9, r10, r5, r7			\n\t"//54
			"umaal r9, r11, r4, r8			\n\t"//45
			"vadd.i64  q11, q11, q6 		\n\t"	//q
			"umaal r9, r12, r3, r6			\n\t"//36
			"str r9, [r2, #4 * 33]			\n\t"//C9
			"vadd.i64  q10, q10, q7 		\n\t"	//q
			"ldr r9, [r2, #4 * 34]			\n\t"//C10
			"ldr r7, [r0, #4 * 19]			\n\t"//B7
			"vadd.i64  q9 , q9 , q8 		\n\t"	//q
			"umaal r9, r10, r5, r8			\n\t"//55
			"umaal r9, r11, r4, r6			\n\t"//46
			"vqadd.u64 d30, d30, d29 		\n\t"	//d
			"umaal r9, r12, r3, r7			\n\t"//37
			"str r9, [r2, #4 * 34]			\n\t"//C10
			"vext.8 d6, d6, d28, #4 		\n\t"
			"ldr r9, [r2, #4 * 35]			\n\t"//C11
			"ldr r8, [r0, #4 * 20]			\n\t"//B8
			"vmlal.u32 q13, d0, d6[0] 		\n\t"	//
			"umaal r9, r10, r5, r6			\n\t"//56
			"umaal r9, r11, r4, r7			\n\t"//47
			"vmlal.u32 q12, d3, d6[0] 		\n\t"	//
			"umaal r9, r12, r3, r8			\n\t"//38
			"str r9, [r2, #4 * 35]			\n\t"//C11
			"vmlal.u32 q11, d1, d6[0] 		\n\t"	//
			"ldr r9, [r2, #4 * 36]			\n\t"//C12
			"ldr r3, [r1, #4 * 18]			\n\t"//A6
			"vmlal.u32 q10, d4, d6[0] 		\n\t"	//
			"umaal r9, r10, r3, r6			\n\t"//66
			"umaal r9, r11, r5, r7			\n\t"//57
			"vmlal.u32 q9 , d2, d6[0] 		\n\t"	//
			"umaal r9, r12, r4, r8			\n\t"//48
			"str r9, [r2, #4 * 36]			\n\t"//C12
			"vmlal.u32 q15, d5, d6[0] 		\n\t"	//
			"ldr r9, [r2, #4 * 37]			\n\t"//C13
			"ldr r4, [r1, #4 * 19]			\n\t"//A7
			"veor q4,q4,q4	 				\n\t"
			"umaal r9, r10, r4, r6			\n\t"//76
			"umaal r9, r11, r3, r7			\n\t"//67
			"veor q5,q5,q5	 				\n\t"
			"umaal r9, r12, r5, r8			\n\t"//58
			"str r9, [r2, #4 * 37]			\n\t"//C13
			"veor q6,q6,q6	 				\n\t"
			"ldr r9, [r2, #4 * 38]			\n\t"//C14
			"ldr r5, [r1, #4 * 20]			\n\t"//A8
			"veor q7,q7,q7	 				\n\t"
			"umaal r9, r10, r5, r6			\n\t"//86
			"umaal r9, r11, r4, r7			\n\t"//77
			"veor q8,q8,q8	 				\n\t"
			"umaal r9, r12, r3, r8			\n\t"//68
			"str r9, [r2, #4 * 38]			\n\t"//C14
			"veor q14,q14,q14 				\n\t"
			"ldr r9, [r2, #4 * 39]			\n\t"//C15
			"ldr r6, [r0, #4 * 21]			\n\t"//B9
			"vtrn.32 q13, q4 				\n\t"
			"umaal r9, r10, r5, r7			\n\t"//87
			"umaal r9, r11, r4, r8			\n\t"//78
			"vtrn.32 q12, q5 				\n\t"
			"umaal r9, r12, r3, r6			\n\t"//69
			"str r9, [r2, #4 * 39]			\n\t"//C15
			"vtrn.32 q11, q6 				\n\t"
			"ldr r9, [r2, #4 * 40]			\n\t"//C16
			"ldr r7, [r0, #4 * 22]			\n\t"//B10
			"vtrn.32 q10, q7 				\n\t"
			"umaal r9, r10, r5, r8			\n\t"//88
			"umaal r9, r11, r4, r6			\n\t"//79
			"vtrn.32 q9 , q8 				\n\t"
			"umaal r9, r12, r3, r7			\n\t"//610
			"str r9, [r2, #4 * 40]			\n\t"//C16
			"vtrn.32 q15, q14 				\n\t"
			"ldr r9, [r2, #4 * 41]			\n\t"//C17
			"ldr r8, [r0, #4 * 23]			\n\t"//B11
			"vadd.i64  q12, q12, q4 		\n\t"	//q
			"umaal r9, r10, r5, r6			\n\t"//89
			"umaal r9, r11, r4, r7			\n\t"//710
			"vadd.i64  q11, q11, q5 		\n\t"	//q
			"umaal r9, r12, r3, r8			\n\t"//611
			"str r9, [r2, #4 * 41]			\n\t"//C17
			"vadd.i64  q10, q10, q6 		\n\t"	//q
			"umaal r10, r11, r5, r7			\n\t"//810
			"umaal r10, r12, r4, r8			\n\t"//711
			"vadd.i64  q9 , q9 , q7 		\n\t"	//q
			"str r10, [r2, #4 * 42]			\n\t"//C18
			"umaal r11, r12, r5, r8			\n\t"//811
			"vadd.i64  q15, q15, q8 		\n\t"	//q
			"str r11, [r2, #4 * 43]			\n\t"//C19
			"str r12, [r2, #4 * 44]			\n\t"//C20
			"vqadd.u64 d28, d28, d27 		\n\t"	//d
			"ldr r6, [r0, #4 * 12]			\n\t"//B0
			"ldr r7, [r0, #4 * 13]			\n\t"//B1
			"vext.8 d6, d6, d26, #4 		\n\t"
			"ldr r8, [r0, #4 * 14]			\n\t"//B2
			"ldr r9, [r2, #4 * 30]			\n\t"//C6
			"veor d6, d6, d7				\n\t"			
			"mov r10, #0					\n\t"
			"umaal r9, r10, r3, r6			\n\t"//60
			"vstr.64 d6, [sp, #4*20]			\n\t"
			"str r9, [r2, #4 * 30]			\n\t"//C6
			"ldr r9, [r2, #4 * 31]			\n\t"//C7
			"vldr.64 d6, [sp, #4*10]		\n\t"			
			"mov r11, #0					\n\t"
			"umaal r9, r10, r4, r6			\n\t"//70
			"vmlal.u32 q12, d0, d6[0] 		\n\t"	//
			"umaal r9, r11, r3, r7			\n\t"//61
			"str r9, [r2, #4 * 31]			\n\t"//C7
			"vmlal.u32 q11, d3, d6[0] 		\n\t"	//
			"ldr r9, [r2, #4 * 32]			\n\t"//C8
			"mov r12, #0					\n\t"
			"vmlal.u32 q10, d1, d6[0] 		\n\t"	//
			"umaal r9, r10, r5, r6			\n\t"//80
			"umaal r9, r11, r4, r7			\n\t"//71
			"vmlal.u32 q9 , d4, d6[0] 		\n\t"	//
			"umaal r9, r12, r3, r8			\n\t"//62
			"str r9, [r2, #4 * 32]			\n\t"//C8
			"vmlal.u32 q15, d2, d6[0] 		\n\t"	//
			"ldr r9, [r2, #4 * 33]			\n\t"//C9
			"ldr r6, [r0, #4 * 15]			\n\t"//B3
			"vmlal.u32 q14, d5, d6[0] 		\n\t"	//
			"umaal r9, r10, r5, r7			\n\t"//81
			"umaal r9, r11, r4, r8			\n\t"//72
			"veor q4,q4,q4	 				\n\t"
			"umaal r9, r12, r3, r6			\n\t"//63
			"str r9, [r2, #4 * 33]			\n\t"//C9
			"veor q5,q5,q5	 				\n\t"
			"ldr r9, [r2, #4 * 34]			\n\t"//C10
			"ldr r7, [r0, #4 * 16]			\n\t"//B4
			"veor q6,q6,q6	 				\n\t"
			"umaal r9, r10, r5, r8			\n\t"//82
			"umaal r9, r11, r4, r6			\n\t"//73
			"veor q7,q7,q7	 				\n\t"
			"umaal r9, r12, r3, r7			\n\t"//64
			"str r9, [r2, #4 * 34]			\n\t"//C10
			"veor q8,q8,q8	 				\n\t"
			"ldr r9, [r2, #4 * 35]			\n\t"//C11
			"ldr r8, [r0, #4 * 17]			\n\t"//B5
			"veor q13,q13,q13 				\n\t"
			"umaal r9, r10, r5, r6			\n\t"//83
			"umaal r9, r11, r4, r7			\n\t"//74
			"vtrn.32 q12, q4 				\n\t"
			"umaal r9, r12, r3, r8			\n\t"//65
			"str r9, [r2, #4 * 35]			\n\t"//C11
			"vtrn.32 q11, q5 				\n\t"
			"ldr r9, [r2, #4 * 36]			\n\t"//C12
			"ldr r3, [r1, #4 * 21]			\n\t"//A9
			"vtrn.32 q10, q6 				\n\t"
			"umaal r9, r10, r3, r6			\n\t"//93
			"umaal r9, r11, r5, r7			\n\t"//84
			"vtrn.32 q9 , q7 				\n\t"
			"umaal r9, r12, r4, r8			\n\t"//75
			"str r9, [r2, #4 * 36]			\n\t"//C12
			"vtrn.32 q15, q8 				\n\t"
			"ldr r9, [r2, #4 * 37]			\n\t"//C13
			"ldr r4, [r1, #4 * 22]			\n\t"//A10
			"vtrn.32 q14, q13 				\n\t"	
			"umaal r9, r10, r4, r6			\n\t"//103
			"umaal r9, r11, r3, r7			\n\t"//94
			"vadd.i64  q11, q11, q4 		\n\t"	//q
			"umaal r9, r12, r5, r8			\n\t"//85
			"str r9, [r2, #4 * 37]			\n\t"//C13
			"vadd.i64  q10, q10, q5 		\n\t"	//q
			"ldr r9, [r2, #4 * 38]			\n\t"//C14
			"ldr r5, [r1, #4 * 23]			\n\t"//A11
			"vadd.i64  q9 , q9 , q6 		\n\t"	//q
			"umaal r9, r10, r5, r6			\n\t"//113
			"umaal r9, r11, r4, r7			\n\t"//104
			"vadd.i64  q15, q15, q7 		\n\t"	//q
			"umaal r9, r12, r3, r8			\n\t"//95
			"str r9, [r2, #4 * 38]			\n\t"//C14
			"vadd.i64  q14, q14, q8 		\n\t"	//q
			"ldr r9, [r2, #4 * 39]			\n\t"//C15
			"ldr r6, [r0, #4 * 18]			\n\t"//B6
			"vqadd.u64 d26, d26, d25 		\n\t"	//d
			"umaal r9, r10, r5, r7			\n\t"//114
			"umaal r9, r11, r4, r8			\n\t"//105
			"vext.8 d6, d6, d24, #4 		\n\t"
			"umaal r9, r12, r3, r6			\n\t"//96
			"str r9, [r2, #4 * 39]			\n\t"//C15
			"vmlal.u32 q11, d0, d6[0] 		\n\t"	//
			"ldr r9, [r2, #4 * 40]			\n\t"//C16
			"ldr r7, [r0, #4 * 19]			\n\t"//B7
			"vmlal.u32 q10, d3, d6[0] 		\n\t"	//
			"umaal r9, r10, r5, r8			\n\t"//115
			"umaal r9, r11, r4, r6			\n\t"//106
			"vmlal.u32 q9 , d1, d6[0] 		\n\t"	//
			"umaal r9, r12, r3, r7			\n\t"//97
			"str r9, [r2, #4 * 40]			\n\t"//C16
			"vmlal.u32 q15, d4, d6[0] 		\n\t"	//
			"ldr r9, [r2, #4 * 41]			\n\t"//C17
			"ldr r8, [r0, #4 * 20]			\n\t"//B8
			"vmlal.u32 q14, d2, d6[0] 		\n\t"	//
			"umaal r9, r10, r5, r6			\n\t"//116
			"umaal r9, r11, r4, r7			\n\t"//107
			"vmlal.u32 q13, d5, d6[0] 		\n\t"	//
			"umaal r9, r12, r3, r8			\n\t"//98
			"str r9, [r2, #4 * 41]			\n\t"//C17
			"veor q4,q4,q4	 				\n\t"
			"ldr r9, [r2, #4 * 42]			\n\t"//C18
			"ldr r6, [r0, #4 * 21]			\n\t"//B9
			"veor q5,q5,q5	 				\n\t"
			"umaal r9, r10, r5, r7			\n\t"//117
			"umaal r9, r11, r4, r8			\n\t"//108
			"veor q6,q6,q6	 				\n\t"
			"umaal r9, r12, r3, r6			\n\t"//99
			"str r9, [r2, #4 * 42]			\n\t"//C18
			"veor q7,q7,q7	 				\n\t"
			"ldr r9, [r2, #4 * 43]			\n\t"//C19
			"ldr r7, [r0, #4 * 22]			\n\t"//B10
			"veor q8,q8,q8	 				\n\t"
			"umaal r9, r10, r5, r8			\n\t"//118
			"umaal r9, r11, r4, r6			\n\t"//109
			"veor q12,q12,q12 				\n\t"
			"umaal r9, r12, r3, r7			\n\t"//910
			"str r9, [r2, #4 * 43]			\n\t"//C19
			"vtrn.32 q11, q4 				\n\t"
			"ldr r9, [r2, #4 * 44]			\n\t"//C20
			"ldr r8, [r0, #4 * 23]			\n\t"//B11
			"vtrn.32 q10, q5 				\n\t"
			"umaal r9, r10, r5, r6			\n\t"//119
			"umaal r9, r11, r4, r7			\n\t"//1010
			"vtrn.32 q9 , q6 				\n\t"
			"umaal r9, r12, r3, r8			\n\t"//9111
			"str r9, [r2, #4 * 44]			\n\t"//C20
			"vtrn.32 q15, q7 				\n\t"
			"umaal r10, r11, r5, r7			\n\t"//1110
			"umaal r10, r12, r4, r8			\n\t"//1011
			"vtrn.32 q14, q8 				\n\t"
			"str r10, [r2, #4 * 45]			\n\t"//C21
			"umaal r11, r12, r5, r8			\n\t"//1111
			"vtrn.32 q13, q12 				\n\t"	
			"str r11, [r2, #4 * 46]			\n\t"//C22
			"str r12, [r2, #4 * 47]			\n\t"//C23
			"vadd.i64  q10, q10, q4 		\n\t"	//q
			"vadd.i64  q9 , q9 , q5 		\n\t"	//q
			"vadd.i64  q15, q15, q6 		\n\t"	//q
			"vadd.i64  q14, q14, q7 		\n\t"	//q
			"vadd.i64  q13, q13, q8 		\n\t"	//q
			"vqadd.u64 d24, d24, d23 		\n\t"	//d
			"vext.8 d6, d6, d22, #4 		\n\t"
			"veor d6, d6, d7				\n\t"
			"vstr.64 d6, [sp, #4*22]		\n\t"


			//integration
			//q14 q15

//integration
			//q14 q15

			"ldr r0, [sp, #4 * 42]			\n\t"//loading result
			"mov r12, #0					\n\t"
			"vmov r14, d7[0]	 			\n\t"//carry			

			
			//1
			"ldr r1, [r0, #4 * 0]			\n\t"
			"ldr r2, [r0, #4 * 1]			\n\t"
			"ldr r3, [r0, #4 * 2]			\n\t"
			"ldr r4, [r0, #4 * 3]			\n\t"

			"ldr r5, [sp, #4 * 12]			\n\t"
			"ldr r6, [sp, #4 * 13]			\n\t"
			"ldr r7, [sp, #4 * 14]			\n\t"
			"ldr r8, [sp, #4 * 15]			\n\t"

			"adds r11, r14, r14				\n\t"//carry

			"adcs r1, r5, r1				\n\t"
			"adcs r2, r6, r2				\n\t"
			"adcs r3, r7, r3				\n\t"
			"adcs r4, r8, r4				\n\t"
			"adcs r11, r12, #0				\n\t"

			"ldr r5, [r0, #4 * 24]			\n\t"
			"ldr r6, [r0, #4 * 25]			\n\t"
			"ldr r7, [r0, #4 * 26]			\n\t"
			"ldr r8, [r0, #4 * 27]			\n\t"

			"adds r1, r5, r1				\n\t"
			"adcs r2, r6, r2				\n\t"
			"adcs r3, r7, r3				\n\t"
			"adcs r4, r8, r4				\n\t"
			"adcs r11, r11, #0				\n\t"

			"str r1, [r0, #4 * 12]			\n\t"
			"str r2, [r0, #4 * 13]			\n\t"
			"str r3, [r0, #4 * 14]			\n\t"
			"str r4, [r0, #4 * 15]			\n\t"

			//2
			"ldr r1, [r0, #4 * 4]			\n\t"
			"ldr r2, [r0, #4 * 5]			\n\t"
			"ldr r3, [r0, #4 * 6]			\n\t"
			"ldr r4, [r0, #4 * 7]			\n\t"

			"ldr r5, [sp, #4 * 16]			\n\t"
			"ldr r6, [sp, #4 * 17]			\n\t"
			"ldr r7, [sp, #4 * 18]			\n\t"
			"ldr r8, [sp, #4 * 19]			\n\t"

			"adds r1, r1, r11				\n\t"
			"adcs r2, r2, #0				\n\t"
			"adcs r3, r3, #0				\n\t"
			"adcs r4, r4, #0				\n\t"
			"adcs r11, r12, #0				\n\t"

			"adds r1, r5, r1				\n\t"
			"adcs r2, r6, r2				\n\t"
			"adcs r3, r7, r3				\n\t"
			"adcs r4, r8, r4				\n\t"
			"adcs r11, r11, #0				\n\t"

			"ldr r5, [r0, #4 * 28]			\n\t"
			"ldr r6, [r0, #4 * 29]			\n\t"
			"ldr r7, [r0, #4 * 30]			\n\t"
			"ldr r8, [r0, #4 * 31]			\n\t"

			"adds r1, r5, r1				\n\t"
			"adcs r2, r6, r2				\n\t"
			"adcs r3, r7, r3				\n\t"
			"adcs r4, r8, r4				\n\t"
			"adcs r11, r11, #0				\n\t"

			"str r1, [r0, #4 * 16]			\n\t"
			"str r2, [r0, #4 * 17]			\n\t"
			"str r3, [r0, #4 * 18]			\n\t"
			"str r4, [r0, #4 * 19]			\n\t"


			//3
			"ldr r1, [r0, #4 * 8]			\n\t"
			"ldr r2, [r0, #4 * 9]			\n\t"
			"ldr r3, [r0, #4 * 10]			\n\t"
			"ldr r4, [r0, #4 * 11]			\n\t"

			"ldr r5, [sp, #4 * 20]			\n\t"
			"ldr r6, [sp, #4 * 21]			\n\t"
			"ldr r7, [sp, #4 * 22]			\n\t"
			"ldr r8, [sp, #4 * 23]			\n\t"

			"adds r1, r1, r11				\n\t"
			"adcs r2, r2, #0				\n\t"
			"adcs r3, r3, #0				\n\t"
			"adcs r4, r4, #0				\n\t"
			"adcs r11, r12, #0				\n\t"

			"adds r1, r5, r1				\n\t"
			"adcs r2, r6, r2				\n\t"
			"adcs r3, r7, r3				\n\t"
			"adcs r4, r8, r4				\n\t"
			"adcs r11, r11, #0				\n\t"

			"ldr r5, [r0, #4 * 32]			\n\t"
			"ldr r6, [r0, #4 * 33]			\n\t"
			"ldr r7, [r0, #4 * 34]			\n\t"
			"ldr r8, [r0, #4 * 35]			\n\t"

			"adds r1, r5, r1				\n\t"
			"adcs r2, r6, r2				\n\t"
			"adcs r3, r7, r3				\n\t"
			"adcs r4, r8, r4				\n\t"
			"adcs r11, r11, #0				\n\t"

			"vmov d6[0], r11				\n\t"//carry

			"str r1, [r0, #4 * 20]			\n\t"
			"str r2, [r0, #4 * 21]			\n\t"
			"str r3, [r0, #4 * 22]			\n\t"
			"str r4, [r0, #4 * 23]			\n\t"


//			"mov r14, #0					\n\t"
//			"adcs r14, r14, #0				\n\t"//carry

			//second part finished

			"vmov r1, r2, d20 				\n\t"
			"vmov r3, r4, d18 				\n\t"
			"vmov r5, r6, d30 				\n\t"
			"vmov r7, r8, d28 				\n\t"
			"vmov r9, r10, d26 				\n\t"
			"vmov r11, r12, d24 			\n\t"

			//"adds r1, r1, r14				\n\t"
			"adds r2, r2, r3				\n\t"
			"adcs r3, r4, r5				\n\t"
			"adcs r4, r6, r7				\n\t"
			"adcs r5, r8, r9				\n\t"
			"adcs r6, r10, r11				\n\t"
			"adcs r12, r12, #0				\n\t"

			"eor r1, r1, r14				\n\t"
			"eor r2, r2, r14				\n\t"
			"eor r3, r3, r14				\n\t"
			"eor r4, r4, r14				\n\t"
			"eor r5, r5, r14				\n\t"
			"eor r6, r6, r14				\n\t"

			"str r1,  [sp, #4 * 12]			\n\t"//ST IN
			"str r2,  [sp, #4 * 13]			\n\t"//ST IN
			"str r3,  [sp, #4 * 14]			\n\t"//ST IN
			"str r4,  [sp, #4 * 15]			\n\t"//ST IN
			"str r5,  [sp, #4 * 16]			\n\t"//ST IN
			"str r6,  [sp, #4 * 17]			\n\t"//ST IN


			"vmov r1, r2, d21 				\n\t"
			"vmov r3, r4, d19 				\n\t"
			"vmov r5, r6, d31 				\n\t"
			"vmov r7, r8, d29 				\n\t"
			"vmov r9, r10, d27 				\n\t"
			"vmov r11,		d25[0] 			\n\t"

			"adds r1, r1, r12				\n\t"
			"adcs r2, r2, r3				\n\t"
			"adcs r3, r4, r5				\n\t"
			"adcs r4, r6, r7				\n\t"
			"adcs r5, r8, r9				\n\t"
			"adcs r6, r10, r11				\n\t"

			"eor r1, r1, r14				\n\t"
			"eor r2, r2, r14				\n\t"
			"eor r3, r3, r14				\n\t"
			"eor r4, r4, r14				\n\t"
			"eor r5, r5, r14				\n\t"
			"eor r6, r6, r14				\n\t"

			"str r1,  [sp, #4 * 18]			\n\t"//ST IN
			"str r2,  [sp, #4 * 19]			\n\t"//ST IN
			"str r3,  [sp, #4 * 20]			\n\t"//ST IN
			"str r4,  [sp, #4 * 21]			\n\t"//ST IN
			"str r5,  [sp, #4 * 22]			\n\t"//ST IN
			"str r6,  [sp, #4 * 23]			\n\t"//ST IN


////////////////////////////////////////////////////////
			"vmov r11, d6[0]				\n\t"//carry
			"mov r12, #0					\n\t"
			//1
			"ldr r1, [r0, #4 * 24]			\n\t"
			"ldr r2, [r0, #4 * 25]			\n\t"
			"ldr r3, [r0, #4 * 26]			\n\t"
			"ldr r4, [r0, #4 * 27]			\n\t"

			"ldr r5, [sp, #4 * 12]			\n\t"
			"ldr r6, [sp, #4 * 13]			\n\t"
			"ldr r7, [sp, #4 * 14]			\n\t"
			"ldr r8, [sp, #4 * 15]			\n\t"

			"adds r1, r1, r11				\n\t"
			"adcs r2, r2, #0				\n\t"
			"adcs r3, r3, #0				\n\t"
			"adcs r4, r4, #0				\n\t"
			"adcs r11, r12, #0				\n\t"

			"adds r1, r5, r1				\n\t"
			"adcs r2, r6, r2				\n\t"
			"adcs r3, r7, r3				\n\t"
			"adcs r4, r8, r4				\n\t"
			"adcs r11, r11, #0				\n\t"

			"ldr r5, [r0, #4 * 36]			\n\t"
			"ldr r6, [r0, #4 * 37]			\n\t"
			"ldr r7, [r0, #4 * 38]			\n\t"
			"ldr r8, [r0, #4 * 39]			\n\t"

			"adds r1, r5, r1				\n\t"
			"adcs r2, r6, r2				\n\t"
			"adcs r3, r7, r3				\n\t"
			"adcs r4, r8, r4				\n\t"
			"adcs r11, r11, #0				\n\t"

			"str r1, [r0, #4 * 24]			\n\t"
			"str r2, [r0, #4 * 25]			\n\t"
			"str r3, [r0, #4 * 26]			\n\t"
			"str r4, [r0, #4 * 27]			\n\t"


			//2
			"ldr r1, [r0, #4 * 28]			\n\t"
			"ldr r2, [r0, #4 * 29]			\n\t"
			"ldr r3, [r0, #4 * 30]			\n\t"
			"ldr r4, [r0, #4 * 31]			\n\t"

			"ldr r5, [sp, #4 * 16]			\n\t"
			"ldr r6, [sp, #4 * 17]			\n\t"
			"ldr r7, [sp, #4 * 18]			\n\t"
			"ldr r8, [sp, #4 * 19]			\n\t"

			"adds r1, r1, r11				\n\t"
			"adcs r2, r2, #0				\n\t"
			"adcs r3, r3, #0				\n\t"
			"adcs r4, r4, #0				\n\t"
			"adcs r11, r12, #0				\n\t"

			"adds r1, r5, r1				\n\t"
			"adcs r2, r6, r2				\n\t"
			"adcs r3, r7, r3				\n\t"
			"adcs r4, r8, r4				\n\t"
			"adcs r11, r11, #0				\n\t"

			"ldr r5, [r0, #4 * 40]			\n\t"
			"ldr r6, [r0, #4 * 41]			\n\t"
			"ldr r7, [r0, #4 * 42]			\n\t"
			"ldr r8, [r0, #4 * 43]			\n\t"

			"adds r1, r5, r1				\n\t"
			"adcs r2, r6, r2				\n\t"
			"adcs r3, r7, r3				\n\t"
			"adcs r4, r8, r4				\n\t"
			"adcs r11, r11, #0				\n\t"

			"str r1, [r0, #4 * 28]			\n\t"
			"str r2, [r0, #4 * 29]			\n\t"
			"str r3, [r0, #4 * 30]			\n\t"
			"str r4, [r0, #4 * 31]			\n\t"

			//3
			"ldr r1, [r0, #4 * 32]			\n\t"
			"ldr r2, [r0, #4 * 33]			\n\t"
			"ldr r3, [r0, #4 * 34]			\n\t"
			"ldr r4, [r0, #4 * 35]			\n\t"

			"ldr r5, [sp, #4 * 20]			\n\t"
			"ldr r6, [sp, #4 * 21]			\n\t"
			"ldr r7, [sp, #4 * 22]			\n\t"
			"ldr r8, [sp, #4 * 23]			\n\t"

			"adds r1, r1, r11				\n\t"
			"adcs r2, r2, #0				\n\t"
			"adcs r3, r3, #0				\n\t"
			"adcs r4, r4, #0				\n\t"
			"adcs r11, r12, #0				\n\t"

			"adds r1, r5, r1				\n\t"
			"adcs r2, r6, r2				\n\t"
			"adcs r3, r7, r3				\n\t"
			"adcs r4, r8, r4				\n\t"
			"adcs r11, r11, #0				\n\t"

			"ldr r5, [r0, #4 * 44]			\n\t"
			"ldr r6, [r0, #4 * 45]			\n\t"
			"ldr r7, [r0, #4 * 46]			\n\t"
			"ldr r8, [r0, #4 * 47]			\n\t"

			"adds r1, r5, r1				\n\t"
			"adcs r2, r6, r2				\n\t"
			"adcs r3, r7, r3				\n\t"
			"adcs r4, r8, r4				\n\t"
			"adcs r11, r11, #0				\n\t"

			"str r1, [r0, #4 * 32]			\n\t"
			"str r2, [r0, #4 * 33]			\n\t"
			"str r3, [r0, #4 * 34]			\n\t"
			"str r4, [r0, #4 * 35]			\n\t"
			//////////////////////////////////////





			"vmov r12, d7[0]	 			\n\t"//carry
			"adds r11, r11, r12				\n\t"//carry
			"asr  r12, r11, #10				\n\t"

			"ldr r1, [r0, #4 * 36]			\n\t"
			"ldr r2, [r0, #4 * 37]			\n\t"
			"ldr r3, [r0, #4 * 38]			\n\t"
			"ldr r4, [r0, #4 * 39]			\n\t"
			"ldr r5, [r0, #4 * 40]			\n\t"
			"ldr r6, [r0, #4 * 41]			\n\t"

			"adds r1, r1, r11				\n\t"
			"adcs r2, r2, r12				\n\t"
			"adcs r3, r3, r12				\n\t"
			"adcs r4, r4, r12				\n\t"
			"adcs r5, r5, r12				\n\t"
			"adcs r6, r6, r12				\n\t"

			"str r1, [r0, #4 * 36]			\n\t"
			"str r2, [r0, #4 * 37]			\n\t"
			"str r3, [r0, #4 * 38]			\n\t"
			"str r4, [r0, #4 * 39]			\n\t"
			"str r5, [r0, #4 * 40]			\n\t"
			"str r6, [r0, #4 * 41]			\n\t"
						
			"ldr r1, [r0, #4 * 42]			\n\t"
			"ldr r2, [r0, #4 * 43]			\n\t"
			"ldr r3, [r0, #4 * 44]			\n\t"
			"ldr r4, [r0, #4 * 45]			\n\t"
			"ldr r5, [r0, #4 * 46]			\n\t"
			"ldr r6, [r0, #4 * 47]			\n\t"

			"adcs r1, r1, r12				\n\t"
			"adcs r2, r2, r12				\n\t"
			"adcs r3, r3, r12				\n\t"
			"adcs r4, r4, r12				\n\t"
			"adcs r5, r5, r12				\n\t"
			"adcs r6, r6, r12				\n\t"

			"str r1, [r0, #4 * 42]			\n\t"
			"str r2, [r0, #4 * 43]			\n\t"
			"str r3, [r0, #4 * 44]			\n\t"
			"str r4, [r0, #4 * 45]			\n\t"
			"str r5, [r0, #4 * 46]			\n\t"
			"str r6, [r0, #4 * 47]			\n\t"			
	/**/		
			/////////////////////////////////////epilogue
			"add sp, #4 * 2 * 12		\n\t"	//384 * 3 (A384, B384, C384)

			"vpop {q4-q7}					\n\t"
			"pop  {r0-r11,pc}				\n\t"

			
	:
	:
	:
	);
}





void __attribute__ ((noinline, naked)) RED512(const felm_t ma, felm_t mb, felm_t mc){




	asm(
			"push  {r3-r11,lr}			\n\t"
			"vpush {q4-q7}				\n\t"
			
			//1 2 3 4 5 6 7 8
			//9 10 11 12

			"sub sp, #4 * 12				\n\t"	//384

			///////////////////////////////////////////////
			"ldr r3, [r1, #4 * 11]			\n\t"//M11			//ROUND1		
			"ldr r4, [r1, #4 * 12]			\n\t"//M12			//ROUND1		
			"ldr r5, [r1, #4 * 13]			\n\t"//M13			//ROUND1		


			"ldr r9,  [r0, #4 * 11]			\n\t"//LD IN11
			"ldr r10,  [r0, #4 * 12]		\n\t"//LD IN12

			"ldr r6, [r0, #4 * 0]			\n\t"//Q0
			"ldr r7, [r0, #4 * 1]			\n\t"//Q1
			"ldr r8, [r0, #4 * 2]			\n\t"//Q2

			"mov r11, #0					\n\t"	
			"mov r12, #0					\n\t"			
			"mov r14, #0					\n\t"	

			"umlal r9, r12, r3, r6			\n\t"//#1
			"umaal r10, r12, r4, r6			\n\t"//#2
			"umaal r10, r14, r3, r7			\n\t"//#3
			"str r9,  [r0, #4 * 11]			\n\t"//ST IN0
			"str r10, [r0, #4 * 12]			\n\t"//ST IN1





			"ldr r9, [r0, #4 * 13]			\n\t"//LD IN2

			"umaal r9, r11, r5, r6			\n\t"//#3
			"umaal r9, r12, r4, r7			\n\t"//#4
			"umaal r9, r14, r3, r8			\n\t"//#5

			"str r9,  [r0, #4 * 13]			\n\t"//LD IN0



//q0-q2 op1
//q3	op2
//q4-q9 rst
//q10-q15
			"vldr.64 d0, [r1, #4*12]		\n\t"
			"vldr.64 d1, [r1, #4*14]		\n\t"
			"vldr.64 d2, [r1, #4*16]		\n\t"
			"vldr.64 d3, [r1, #4*18]		\n\t"
			"vldr.64 d4, [r1, #4*20]		\n\t"
			"vldr.64 d5, [r1, #4*22]		\n\t"

			"vtrn.32 d0, d3 				\n\t"
			"vtrn.32 d1, d4 				\n\t"
			"vtrn.32 d2, d5 				\n\t"

			"vldr.64 d8, [r0, #4*24]		\n\t"
			"vldr.64 d9, [r0, #4*26]		\n\t"
			"vldr.64 d10, [r0, #4*28]		\n\t"		
			"vldr.64 d11, [r0, #4*30]		\n\t"
			"vldr.64 d12, [r0, #4*32]		\n\t"
			"vldr.64 d13, [r0, #4*34]		\n\t"

			"vtrn.32 d8, d11 				\n\t"
			"vtrn.32 d9, d12 				\n\t"
			"vtrn.32 d10, d13 				\n\t"

			"vmov.i64     q3, 0xFFFFFFFF    \n\t"// mast 1
			"vshr.u64     q3, q3, #31		\n\t"

			"vmull.u32 q15, d8 , d7[0] 		\n\t"	//A6A0 X B0B0			
			"vmull.u32 q14, d11, d7[0] 		\n\t"	//A7A1 X B0B0
			"vmull.u32 q13, d9 , d7[0] 		\n\t"	//A8A2 X B0B0			
			"vmull.u32 q12, d12, d7[0] 		\n\t"	//A9A3 X B0B0		
			"vmull.u32 q11, d10, d7[0] 		\n\t"	//A10A4 X B0B0		
			"vmull.u32 q10, d13, d7[0] 		\n\t"	//A11A5 X B0B0		
			
			"vmov d6[0], r10 				\n\t"
			"vmov d6[1], r9 				\n\t"//NOW NEON

/*
			//TEST
			//"str r9,  [r2, #4 * 0]			\n\t"//ST IN0
			"vstr.64 d30, [r2, #4*0]		\n\t"
			"vstr.64 d28, [r2, #4*2]		\n\t"
			"vstr.64 d26, [r2, #4*4]		\n\t"
			"vstr.64 d24, [r2, #4*6]		\n\t"
			"vstr.64 d22, [r2, #4*8]		\n\t"
			"vstr.64 d20, [r2, #4*10]		\n\t"

			"vstr.64 d31, [r2, #4*12]		\n\t"
			"vstr.64 d29, [r2, #4*14]		\n\t"
			"vstr.64 d27, [r2, #4*16]		\n\t"
			"vstr.64 d25, [r2, #4*18]		\n\t"
			"vstr.64 d23, [r2, #4*20]		\n\t"
			"vstr.64 d21, [r2, #4*22]		\n\t"
			//TEST
			*/
////////////////////////////////////////////////////////////////////////Section#1
"ldr r9, [r0, #4 * 14]			\n\t"//LD IN
			"ldr r3, [r1, #4 * 14]			\n\t"//LD M
"vmlal.u32 q15, d0, d6[0] 		\n\t"	//
			"umaal r9, r11, r3, r6			\n\t"//
			"umaal r9, r12, r5, r7			\n\t"//
			"vmlal.u32 q14, d3, d6[0] 		\n\t"	//
			"umaal r9, r14, r4, r8			\n\t"//
			"str r9, [r0, #4 * 14]			\n\t"//ST IN
			"vmlal.u32 q13, d1, d6[0] 		\n\t"	//
			"ldr r9, [r0, #4 * 15]			\n\t"//LD IN
			"ldr r4, [r1, #4 * 15]			\n\t"//LD M
			"vmlal.u32 q12, d4, d6[0] 		\n\t"	//
			"umaal r9, r11, r4, r6			\n\t"//
			"umaal r9, r12, r3, r7			\n\t"//
			"vmlal.u32 q11, d2, d6[0] 		\n\t"	//
			"umaal r9, r14, r5, r8			\n\t"//
			"str r9, [r0, #4 * 15]			\n\t"//ST IN			
			"vmlal.u32 q10, d5, d6[0] 		\n\t"	//
			"ldr r9, [r0, #4 * 16]			\n\t"//LD IN
			"ldr r5, [r1, #4 * 16]			\n\t"//LD M
			"veor q4,q4,q4	 				\n\t"
			"umaal r9, r11, r5, r6			\n\t"//
			"umaal r9, r12, r4, r7			\n\t"//
			"veor q5,q5,q5	 				\n\t"
			"umaal r9, r14, r3, r8			\n\t"//
			"str r9, [r0, #4 * 16]			\n\t"//ST IN
			"veor q6,q6,q6	 				\n\t"
			"ldr r9, [r0, #4 * 17]			\n\t"//LD IN
			"ldr r3, [r1, #4 * 17]			\n\t"//LD M
			"veor q7,q7,q7	 				\n\t"
			"umaal r9, r11, r3, r6			\n\t"//
			"umaal r9, r12, r5, r7			\n\t"//
			"veor q8,q8,q8	 				\n\t"
			"umaal r9, r14, r4, r8			\n\t"//
			"str r9, [r0, #4 * 17]			\n\t"//ST IN
			"veor q9,q9,q9	 				\n\t"
			"ldr r9, [r0, #4 * 18]			\n\t"//LD IN
			"ldr r4, [r1, #4 * 18]			\n\t"//LD M
			"vtrn.32 q15, q4 				\n\t"
			"umaal r9, r11, r4, r6			\n\t"//
			"umaal r9, r12, r3, r7			\n\t"//
			"vtrn.32 q14, q5 				\n\t"
			"umaal r9, r14, r5, r8			\n\t"//
			"str r9, [r0, #4 * 18]			\n\t"//ST IN
			"vtrn.32 q13, q6 				\n\t"
			"ldr r9, [r0, #4 * 19]			\n\t"//LD IN
			"ldr r5, [r1, #4 * 19]			\n\t"//LD M
			"vtrn.32 q12, q7 				\n\t"
			"umaal r9, r11, r5, r6			\n\t"//
			"umaal r9, r12, r4, r7			\n\t"//
			"vtrn.32 q11, q8 				\n\t"
			"umaal r9, r14, r3, r8			\n\t"//
			"str r9, [r0, #4 * 19]			\n\t"//ST IN
			"vtrn.32 q10, q9 				\n\t"
			"ldr r9, [r0, #4 * 20]			\n\t"//LD IN
			"ldr r3, [r1, #4 * 20]			\n\t"//LD M
			"vadd.i64  q14, q14, q4 		\n\t"	//q
			"umaal r9, r11, r3, r6			\n\t"//
			"umaal r9, r12, r5, r7			\n\t"//
			"vadd.i64  q13, q13, q5 		\n\t"	//q
			"umaal r9, r14, r4, r8			\n\t"//
			"str r9, [r0, #4 * 20]			\n\t"//ST IN
			"vadd.i64  q12, q12, q6 		\n\t"	//q
			"ldr r9, [r0, #4 * 21]			\n\t"//LD IN
			"ldr r4, [r1, #4 * 21]			\n\t"//LD M
			"vadd.i64  q11, q11, q7 		\n\t"	//q
			"umaal r9, r11, r4, r6			\n\t"//
			"umaal r9, r12, r3, r7			\n\t"//
			"vadd.i64  q10, q10, q8 		\n\t"	//q
			"umaal r9, r14, r5, r8			\n\t"//
			"str r9, [r0, #4 * 21]			\n\t"//ST IN
			"vqadd.u64 d18, d18, d31 		\n\t"	//d
			"ldr r9, [r0, #4 * 22]			\n\t"//LD IN
			"ldr r5, [r1, #4 * 22]			\n\t"//LD M
			"vext.8 d6, d6, d30, #4 		\n\t"
			"umaal r9, r11, r5, r6			\n\t"//
			"umaal r9, r12, r4, r7			\n\t"//
			"vmlal.u32 q14, d0, d6[0] 		\n\t"	//
			"umaal r9, r14, r3, r8			\n\t"//
			"str r9, [r0, #4 * 22]			\n\t"//ST IN
			"vmlal.u32 q13, d3, d6[0] 		\n\t"	//
			"ldr r9, [r0, #4 * 23]			\n\t"//LD IN
			"ldr r3, [r1, #4 * 23]			\n\t"//LD M
			"vmlal.u32 q12, d1, d6[0] 		\n\t"	//
			"umaal r9, r11, r3, r6			\n\t"//
			"umaal r9, r12, r5, r7			\n\t"//
			"vmlal.u32 q11, d4, d6[0] 		\n\t"	//
			"umaal r9, r14, r4, r8			\n\t"//
			"str r9, [r0, #4 * 23]			\n\t"//ST IN
			"vmlal.u32 q10, d2, d6[0] 		\n\t"	//
			"umaal r11, r12, r3, r7			\n\t"//
			"umaal r11, r14, r5, r8			\n\t"//
			"vmlal.u32 q9 , d5, d6[0] 		\n\t"	//
			"str r11, [r0, #4 * 24]			\n\t"//ST IN
			"umaal r12, r14, r3, r8			\n\t"//
			"veor q4,q4,q4	 				\n\t"
			"str r12, [r0, #4 * 25]			\n\t"//ST IN
			"str r14, [r0, #4 * 26]			\n\t"//ST IN
			"veor q5,q5,q5	 				\n\t"
			"ldr r3, [r1, #4 * 11]			\n\t"//M
			"ldr r4, [r1, #4 * 12]			\n\t"//M
			"veor q6,q6,q6	 				\n\t"
			"ldr r5, [r1, #4 * 13]			\n\t"//M
			"ldr r9,  [r0, #4 * 14]			\n\t"//LD IN
			"veor q7,q7,q7	 				\n\t"
			"ldr r10,  [r0, #4 * 15]		\n\t"//LD IN
			"ldr r6, [r0, #4 * 3]			\n\t"//Q
			"veor q8,q8,q8	 				\n\t"
			"ldr r7, [r0, #4 * 4]			\n\t"//Q
			"ldr r8, [r0, #4 * 5]			\n\t"//Q
			"veor q15,q15,q15				\n\t"
			"mov r11, #0					\n\t"	
			"mov r12, #0					\n\t"			
			"vtrn.32 q14, q4 				\n\t"
			"mov r14, #0					\n\t"	
			"umlal r9, r12, r3, r6			\n\t"//#1
			"vtrn.32 q13, q5 				\n\t"
			"umaal r10, r12, r4, r6			\n\t"//#2
			"umaal r10, r14, r3, r7			\n\t"//#3
			"vtrn.32 q12, q6 				\n\t"
			"str r9,  [r0, #4 * 14]			\n\t"//ST IN
			"str r10, [r0, #4 * 15]			\n\t"//ST IN
			"vtrn.32 q11, q7 				\n\t"
			
			"vtrn.32 q10, q8 				\n\t"
			"vtrn.32 q9 , q15 				\n\t"
			"vadd.i64  q13, q13, q4 		\n\t"	//q
			"vadd.i64  q12, q12, q5 		\n\t"	//q
			"vadd.i64  q11, q11, q6 		\n\t"	//q
			"vadd.i64  q10, q10, q7 		\n\t"	//q
			"vadd.i64  q9, q9, q8 			\n\t"	//q
			"vqadd.u64 d30, d30, d29 		\n\t"	//d
			"vext.8 d6, d6, d28, #4 		\n\t"
			"vstr.64 d6, [sp, #4*0]			\n\t"	//saving


//////////////////////////////////////////////////////////////////////




			"vmov d7[0], r9 				\n\t"
			"vmov d7[1], r10 				\n\t"//NOW NEON

			"ldr r9, [r0, #4 * 16]			\n\t"//LD IN2

			"umaal r9, r11, r5, r6			\n\t"//#3
			"umaal r9, r12, r4, r7			\n\t"//#4
			"umaal r9, r14, r3, r8			\n\t"//#5

			"str r9,  [r0, #4 * 16]			\n\t"//LD IN0

			"vmov d6[0], r9 				\n\t"



////////////////////////////////////////////////////////////////////////Section#2

			


			//NEON

			//
			"ldr r9, [r0, #4 * 17]			\n\t"//LD IN
			"vmlal.u32 q13, d0, d7[0] 		\n\t"	//
			"ldr r3, [r1, #4 * 14]			\n\t"//LD M
			"vmlal.u32 q12, d3, d7[0] 		\n\t"	//
			"umaal r9, r11, r3, r6			\n\t"//
			"vmlal.u32 q11, d1, d7[0] 		\n\t"	//
			"umaal r9, r12, r5, r7			\n\t"//
			"vmlal.u32 q10, d4, d7[0] 		\n\t"	//
			"umaal r9, r14, r4, r8			\n\t"//
			"vmlal.u32 q9 , d2, d7[0] 		\n\t"	//
			"str r9, [r0, #4 * 17]			\n\t"//ST IN
			"vmlal.u32 q15, d5, d7[0] 		\n\t"	//
			"ldr r9, [r0, #4 * 18]			\n\t"//LD IN
			"veor q4,q4,q4	 				\n\t"
			"ldr r4, [r1, #4 * 15]			\n\t"//LD M
			"veor q5,q5,q5	 				\n\t"
			"umaal r9, r11, r4, r6			\n\t"//
			"veor q6,q6,q6	 				\n\t"
			"umaal r9, r12, r3, r7			\n\t"//
			"veor q7,q7,q7	 				\n\t"
			"umaal r9, r14, r5, r8			\n\t"//
			"veor q8,q8,q8	 				\n\t"
			"str r9, [r0, #4 * 18]			\n\t"//ST IN			
			"veor q14,q14,q14				\n\t"
			"ldr r9, [r0, #4 * 19]			\n\t"//LD IN
			"vtrn.32 q13, q4 				\n\t"
			"ldr r5, [r1, #4 * 16]			\n\t"//LD M
			"vtrn.32 q12, q5 				\n\t"
			"umaal r9, r11, r5, r6			\n\t"//
			"vtrn.32 q11, q6 				\n\t"
			"umaal r9, r12, r4, r7			\n\t"//
			"vtrn.32 q10, q7 				\n\t"
			"umaal r9, r14, r3, r8			\n\t"//
			"vtrn.32 q9 , q8 				\n\t"
			"str r9, [r0, #4 * 19]			\n\t"//ST IN
			"vtrn.32 q15, q14 				\n\t"
			"ldr r9, [r0, #4 * 20]			\n\t"//LD IN
			"vadd.i64  q12, q12, q4 		\n\t"	//q
			"ldr r3, [r1, #4 * 17]			\n\t"//LD M
			"vadd.i64  q11, q11, q5 		\n\t"	//q
			"umaal r9, r11, r3, r6			\n\t"//
			"vadd.i64  q10, q10, q6 		\n\t"	//q
			"umaal r9, r12, r5, r7			\n\t"//
			"vadd.i64  q9 , q9 , q7 		\n\t"	//q
			"umaal r9, r14, r4, r8			\n\t"//
			"vadd.i64  q15, q15, q8 			\n\t"	//q
			"str r9, [r0, #4 * 20]			\n\t"//ST IN
			"vqadd.u64 d28, d28, d27 		\n\t"	//d
			"ldr r9, [r0, #4 * 21]			\n\t"//LD IN
			"vext.8 d7, d7, d26, #4 		\n\t"
			"ldr r4, [r1, #4 * 18]			\n\t"//LD M
			"vmlal.u32 q12, d0, d7[0] 		\n\t"	//
			"umaal r9, r11, r4, r6			\n\t"//
			"vmlal.u32 q11, d3, d7[0] 		\n\t"	//
			"umaal r9, r12, r3, r7			\n\t"//
			"vmlal.u32 q10, d1, d7[0] 		\n\t"	//
			"umaal r9, r14, r5, r8			\n\t"//
			"vmlal.u32 q9 , d4, d7[0] 		\n\t"	//
			"str r9, [r0, #4 * 21]			\n\t"//ST IN
			"vmlal.u32 q15, d2, d7[0] 		\n\t"	//
			"ldr r9, [r0, #4 * 22]			\n\t"//LD IN
			"vmlal.u32 q14, d5, d7[0] 		\n\t"	//
			"ldr r5, [r1, #4 * 19]			\n\t"//LD M
			"veor q4,q4,q4	 				\n\t"
			"umaal r9, r11, r5, r6			\n\t"//
			"veor q5,q5,q5	 				\n\t"
			"umaal r9, r12, r4, r7			\n\t"//
			"veor q6,q6,q6	 				\n\t"
			"umaal r9, r14, r3, r8			\n\t"//
			"veor q7,q7,q7	 				\n\t"
			"str r9, [r0, #4 * 22]			\n\t"//ST IN
			"veor q8,q8,q8	 				\n\t"
			"ldr r9, [r0, #4 * 23]			\n\t"//LD IN
			"veor q13,q13,q13				\n\t"
			"ldr r3, [r1, #4 * 20]			\n\t"//LD M
			"vtrn.32 q12, q4 				\n\t"
			"umaal r9, r11, r3, r6			\n\t"//
			"vtrn.32 q11, q5 				\n\t"
			"umaal r9, r12, r5, r7			\n\t"//
			"vtrn.32 q10, q6 				\n\t"
			"umaal r9, r14, r4, r8			\n\t"//
			"vtrn.32 q9 , q7 				\n\t"
			"str r9, [r0, #4 * 23]			\n\t"//ST IN
			"vtrn.32 q15, q8 				\n\t"
			"ldr r9, [r0, #4 * 24]			\n\t"//LD IN
			"vtrn.32 q14, q13 				\n\t"
			"ldr r4, [r1, #4 * 21]			\n\t"//LD M
			"vadd.i64  q11, q11, q4 		\n\t"	//q
			"umaal r9, r11, r4, r6			\n\t"//
			"vadd.i64  q10, q10, q5 		\n\t"	//q
			"umaal r9, r12, r3, r7			\n\t"//
			"vadd.i64  q9 , q9 , q6 		\n\t"	//q
			"umaal r9, r14, r5, r8			\n\t"//
			"vadd.i64  q15, q15, q7 		\n\t"	//q
			"str r9, [r0, #4 * 24]			\n\t"//ST IN
			"vadd.i64  q14, q14, q8 		\n\t"	//q
			"ldr r9, [r0, #4 * 25]			\n\t"//LD IN
			"vqadd.u64 d26, d26, d25 		\n\t"	//d
			"ldr r5, [r1, #4 * 22]			\n\t"//LD M
			"vext.8 d7, d7, d24, #4 		\n\t"
			"umaal r9, r11, r5, r6			\n\t"//
			"vstr.64 d7, [sp, #4*2]			\n\t"	//saving
			"umaal r9, r12, r4, r7			\n\t"//
			"vmlal.u32 q11, d0, d6[0] 		\n\t"	//
			"umaal r9, r14, r3, r8			\n\t"//
			"vmlal.u32 q10, d3, d6[0] 		\n\t"	//
			"str r9, [r0, #4 * 25]			\n\t"//ST IN
			"vmlal.u32 q9 , d1, d6[0] 		\n\t"	//
			"ldr r9, [r0, #4 * 26]			\n\t"//LD IN
			"vmlal.u32 q15, d4, d6[0] 		\n\t"	//
			"ldr r3, [r1, #4 * 23]			\n\t"//LD M
			"vmlal.u32 q14, d2, d6[0] 		\n\t"	//
			"umaal r9, r11, r3, r6			\n\t"//
			"vmlal.u32 q13, d5, d6[0] 		\n\t"	//
			"umaal r9, r12, r5, r7			\n\t"//
			"veor q4,q4,q4	 				\n\t"
			"umaal r9, r14, r4, r8			\n\t"//
			"veor q5,q5,q5	 				\n\t"
			"str r9, [r0, #4 * 26]			\n\t"//ST IN
			"veor q6,q6,q6	 				\n\t"
			"umaal r11, r12, r3, r7			\n\t"//
			"veor q7,q7,q7	 				\n\t"
			"umaal r11, r14, r5, r8			\n\t"//
			"veor q8,q8,q8	 				\n\t"
			"str r11, [r0, #4 * 27]			\n\t"//ST IN
			"veor q12,q12,q12				\n\t"
			"umaal r12, r14, r3, r8			\n\t"//
			"vtrn.32 q11, q4 				\n\t"
			"str r12, [r0, #4 * 28]			\n\t"//ST IN
			"vtrn.32 q10, q5 				\n\t"
			"str r14, [r0, #4 * 29]			\n\t"//ST IN
			"vtrn.32 q9 , q6 				\n\t"
			"ldr r3, [r1, #4 * 11]			\n\t"//M
			"vtrn.32 q15, q7 				\n\t"
			"ldr r4, [r1, #4 * 12]			\n\t"//M
			"vtrn.32 q14, q8 				\n\t"
			"ldr r5, [r1, #4 * 13]			\n\t"//M
			"vtrn.32 q13, q12 				\n\t"
			"ldr r9,  [r0, #4 * 17]			\n\t"//LD IN
			"vadd.i64  q10, q10, q4 		\n\t"	//q
			"ldr r10,  [r0, #4 * 18]		\n\t"//LD IN
			"vadd.i64  q9 , q9 , q5 		\n\t"	//q
			"ldr r6, [r0, #4 * 6]			\n\t"//Q
			"vadd.i64  q15, q15, q6 		\n\t"	//q
			"ldr r7, [r0, #4 * 7]			\n\t"//Q
			"vadd.i64  q14, q14, q7 		\n\t"	//q
			"ldr r8, [r0, #4 * 8]			\n\t"//Q
			"vadd.i64  q13, q13, q8 		\n\t"	//q
			"mov r11, #0					\n\t"	
			"vqadd.u64 d24, d24, d23 		\n\t"	//d
			"mov r12, #0					\n\t"			
			"vext.8 d6, d6, d22, #4 		\n\t"
			"mov r14, #0					\n\t"	
			"umlal r9, r12, r3, r6			\n\t"//#1
			"umaal r10, r12, r4, r6			\n\t"//#2
			"umaal r10, r14, r3, r7			\n\t"//#3
			"str r9,  [r0, #4 * 17]			\n\t"//ST IN
			"str r10, [r0, #4 * 18]			\n\t"//ST IN			
			




			//





			"vmov d6[0], r9 				\n\t"
			"vmov d7[0], r10 				\n\t"//NOW NEON

			"ldr r9, [r0, #4 * 19]			\n\t"//LD IN2

			"umaal r9, r11, r5, r6			\n\t"//#3
			"umaal r9, r12, r4, r7			\n\t"//#4
			"umaal r9, r14, r3, r8			\n\t"//#5

			"str r9,  [r0, #4 * 19]			\n\t"//LD IN0

			"vmov d7[1], r9 				\n\t"

			
////////////////////////////////////////////////////////////////////////Section#3


			//NEON
			//d6[0], d7[0], d7[1]
			//
			"ldr r9, [r0, #4 * 20]			\n\t"//LD IN
			"vmlal.u32 q10, d0, d6[0] 		\n\t"	//
			"ldr r3, [r1, #4 * 14]			\n\t"//LD M
			"vmlal.u32 q9 , d3, d6[0] 		\n\t"	//
			"umaal r9, r11, r3, r6			\n\t"//
			"vmlal.u32 q15, d1, d6[0] 		\n\t"	//
			"umaal r9, r12, r5, r7			\n\t"//
			"vmlal.u32 q14, d4, d6[0] 		\n\t"	//
			"umaal r9, r14, r4, r8			\n\t"//
			"vmlal.u32 q13, d2, d6[0] 		\n\t"	//
			"str r9, [r0, #4 * 20]			\n\t"//ST IN
			"vmlal.u32 q12, d5, d6[0] 		\n\t"	//
			"ldr r9, [r0, #4 * 21]			\n\t"//LD IN
			"veor q4,q4,q4	 				\n\t"
			"ldr r4, [r1, #4 * 15]			\n\t"//LD M
			"veor q5,q5,q5	 				\n\t"
			"umaal r9, r11, r4, r6			\n\t"//
			"veor q6,q6,q6	 				\n\t"
			"umaal r9, r12, r3, r7			\n\t"//
			"veor q7,q7,q7	 				\n\t"
			"umaal r9, r14, r5, r8			\n\t"//
			"veor q8,q8,q8	 				\n\t"
			"str r9, [r0, #4 * 21]			\n\t"//ST IN			
			"veor q11,q11,q11				\n\t"
			"ldr r9, [r0, #4 * 22]			\n\t"//LD IN
			"vtrn.32 q10, q4 				\n\t"
			"ldr r5, [r1, #4 * 16]			\n\t"//LD M
			"vtrn.32 q9, q5 				\n\t"
			"umaal r9, r11, r5, r6			\n\t"//
			"vtrn.32 q15, q6 				\n\t"
			"umaal r9, r12, r4, r7			\n\t"//
			"vtrn.32 q14, q7 				\n\t"
			"umaal r9, r14, r3, r8			\n\t"//
			"vtrn.32 q13, q8 				\n\t"
			"str r9, [r0, #4 * 22]			\n\t"//ST IN
			"vtrn.32 q12, q11 				\n\t"
			"ldr r9, [r0, #4 * 23]			\n\t"//LD IN
			"vadd.i64  q9 , q9 , q4 		\n\t"	//q
			"ldr r3, [r1, #4 * 17]			\n\t"//LD M
			"vadd.i64  q15, q15, q5 		\n\t"	//q
			"umaal r9, r11, r3, r6			\n\t"//
			"vadd.i64  q14, q14, q6 		\n\t"	//q
			"umaal r9, r12, r5, r7			\n\t"//
			"vadd.i64  q13, q13, q7 		\n\t"	//q
			"umaal r9, r14, r4, r8			\n\t"//
			"vadd.i64  q12, q12, q8 		\n\t"	//q
			"str r9, [r0, #4 * 23]			\n\t"//ST IN
			"vqadd.u64 d22, d22, d21 		\n\t"	//d
			"ldr r9, [r0, #4 * 24]			\n\t"//LD IN
			"vext.8 d6, d6, d20, #4 		\n\t"
			"ldr r4, [r1, #4 * 18]			\n\t"//LD M
			"vstr.64 d6, [sp, #4*4]			\n\t"	//saving
			"umaal r9, r11, r4, r6			\n\t"//
			"vmlal.u32 q9 , d0, d7[0] 		\n\t"	//
			"umaal r9, r12, r3, r7			\n\t"//
			"vmlal.u32 q15, d3, d7[0] 		\n\t"	//
			"umaal r9, r14, r5, r8			\n\t"//
			"vmlal.u32 q14, d1, d7[0] 		\n\t"	//
			"str r9, [r0, #4 * 24]			\n\t"//ST IN
			"vmlal.u32 q13, d4, d7[0] 		\n\t"	//
			"ldr r9, [r0, #4 * 25]			\n\t"//LD IN
			"vmlal.u32 q12, d2, d7[0] 		\n\t"	//
			"ldr r5, [r1, #4 * 19]			\n\t"//LD M
			"vmlal.u32 q11, d5, d7[0] 		\n\t"	//
			"umaal r9, r11, r5, r6			\n\t"//
			"veor q4,q4,q4	 				\n\t"
			"umaal r9, r12, r4, r7			\n\t"//
			"veor q5,q5,q5	 				\n\t"
			"umaal r9, r14, r3, r8			\n\t"//
			"veor q6,q6,q6	 				\n\t"
			"str r9, [r0, #4 * 25]			\n\t"//ST IN
			"veor q7,q7,q7	 				\n\t"
			"ldr r9, [r0, #4 * 26]			\n\t"//LD IN
			"veor q8,q8,q8	 				\n\t"
			"ldr r3, [r1, #4 * 20]			\n\t"//LD M
			"veor q10,q10,q10				\n\t"
			"umaal r9, r11, r3, r6			\n\t"//
			"vtrn.32 q9 , q4 				\n\t"
			"umaal r9, r12, r5, r7			\n\t"//
			"vtrn.32 q15, q5 				\n\t"
			"umaal r9, r14, r4, r8			\n\t"//
			"vtrn.32 q14, q6 				\n\t"
			"str r9, [r0, #4 * 26]			\n\t"//ST IN
			"vtrn.32 q13, q7 				\n\t"
			"ldr r9, [r0, #4 * 27]			\n\t"//LD IN
			"vtrn.32 q12, q8 				\n\t"
			"ldr r4, [r1, #4 * 21]			\n\t"//LD M
			"vtrn.32 q11, q10 				\n\t"
			"umaal r9, r11, r4, r6			\n\t"//
			"vadd.i64  q15, q15, q4 		\n\t"	//q
			"umaal r9, r12, r3, r7			\n\t"//
			"vadd.i64  q14, q14, q5 		\n\t"	//q
			"umaal r9, r14, r5, r8			\n\t"//
			"vadd.i64  q13, q13, q6 		\n\t"	//q
			"str r9, [r0, #4 * 27]			\n\t"//ST IN
			"vadd.i64  q12, q12, q7 		\n\t"	//q
			"ldr r9, [r0, #4 * 28]			\n\t"//LD IN
			"vadd.i64  q11, q11, q8 		\n\t"	//q
			"ldr r5, [r1, #4 * 22]			\n\t"//LD M
			"vqadd.u64 d20, d20, d19 		\n\t"	//d
			"umaal r9, r11, r5, r6			\n\t"//
			"vext.8 d7, d7, d18, #4 		\n\t"		
			"umaal r9, r12, r4, r7			\n\t"//
			"vmlal.u32 q15, d0, d7[0] 		\n\t"	//
			"umaal r9, r14, r3, r8			\n\t"//
			"vmlal.u32 q14, d3, d7[0] 		\n\t"	//
			"str r9, [r0, #4 * 28]			\n\t"//ST IN
			"vmlal.u32 q13, d1, d7[0] 		\n\t"	//
			"ldr r9, [r0, #4 * 29]			\n\t"//LD IN
			"vmlal.u32 q12, d4, d7[0] 		\n\t"	//
			"ldr r3, [r1, #4 * 23]			\n\t"//LD M
			"vmlal.u32 q11, d2, d7[0] 		\n\t"	//
			"umaal r9, r11, r3, r6			\n\t"//
			"vmlal.u32 q10, d5, d7[0] 		\n\t"	//
			"umaal r9, r12, r5, r7			\n\t"//
			"veor q4,q4,q4	 				\n\t"
			"umaal r9, r14, r4, r8			\n\t"//
			"veor q5,q5,q5	 				\n\t"
			"str r9, [r0, #4 * 29]			\n\t"//ST IN
			"veor q6,q6,q6	 				\n\t"
			"umaal r11, r12, r3, r7			\n\t"//
			"veor q7,q7,q7	 				\n\t"
			"umaal r11, r14, r5, r8			\n\t"//
			"veor q8,q8,q8	 				\n\t"
			"str r11, [r0, #4 * 30]			\n\t"//ST IN
			"veor q9,q9,q9					\n\t"
			"umaal r12, r14, r3, r8			\n\t"//
			"vtrn.32 q15, q4 				\n\t"
			"str r12, [r0, #4 * 31]			\n\t"//ST IN
			"vtrn.32 q14, q5 				\n\t"
			"str r14, [r0, #4 * 32]			\n\t"//ST IN
			"vtrn.32 q13, q6 				\n\t"
						"ldr r3, [r1, #4 * 11]			\n\t"//M
			"vtrn.32 q12, q7 				\n\t"
			"ldr r4, [r1, #4 * 12]			\n\t"//M
			"vtrn.32 q11, q8 				\n\t"
			"ldr r5, [r1, #4 * 13]			\n\t"//M
			"vtrn.32 q10, q9 				\n\t"
			"ldr r9,  [r0, #4 * 20]			\n\t"//LD IN
			"vadd.i64  q14, q14, q4 		\n\t"	//q
			"ldr r10,  [r0, #4 * 21]		\n\t"//LD IN
			"vadd.i64  q13, q13, q5 		\n\t"	//q
			"ldr r6, [r0, #4 * 9]			\n\t"//Q
			"vadd.i64  q12, q12, q6 		\n\t"	//q
			"ldr r7, [r0, #4 * 10]			\n\t"//Q
			"vadd.i64  q11, q11, q7 		\n\t"	//q
			"ldr r8, [r0, #4 * 11]			\n\t"//Q
			"vadd.i64  q10, q10, q8 		\n\t"	//q
			"mov r11, #0					\n\t"	
			"vqadd.u64 d18, d18, d31 		\n\t"	//d
			"mov r12, #0					\n\t"			
			"vext.8 d7, d7, d30, #4 		\n\t"
			"mov r14, #0					\n\t"	
			"vstr.64 d7, [sp, #4*6]			\n\t"	//saving
			"umlal r9, r12, r3, r6			\n\t"//#1
			"umaal r10, r12, r4, r6			\n\t"//#2
			"umaal r10, r14, r3, r7			\n\t"//#3
			"str r9,  [r0, #4 * 20]			\n\t"//ST IN
			"str r10, [r0, #4 * 21]			\n\t"//ST IN





////////////////////////TEST START
//			"str r9,  [r2, #4 * 0]			\n\t"//ST IN0
//			"str r10, [r2, #4 * 1]			\n\t"//ST IN1
////////////////////////TEST END


			"vmov d6[0], r9 				\n\t"
			"vmov d6[1], r10 				\n\t"//NOW NEON

			"ldr r9, [r0, #4 * 22]			\n\t"//LD IN2

			"umaal r9, r11, r5, r6			\n\t"//#3
			"umaal r9, r12, r4, r7			\n\t"//#4
			"umaal r9, r14, r3, r8			\n\t"//#5

			"str r9,  [r0, #4 * 22]			\n\t"//LD IN0
////////////////////////TEST START
//			"str r9,  [r2, #4 * 2]			\n\t"//ST IN0
////////////////////////TEST END
			"vmov d7[0], r9 				\n\t"
			
			//
////////////////////////////////////////////////////////////////////////Section#4			
				"ldr r9, [r0, #4 * 23]			\n\t"//LD IN
			"vmlal.u32 q14, d0, d6[0] 		\n\t"	//
			"vmlal.u32 q13, d3, d6[0] 		\n\t"	//
			"ldr r3, [r1, #4 * 14]			\n\t"//LD M
			"vmlal.u32 q12, d1, d6[0] 		\n\t"	//
			"vmlal.u32 q11, d4, d6[0] 		\n\t"	//
			"umaal r9, r11, r3, r6			\n\t"//
			"vmlal.u32 q10, d2, d6[0] 		\n\t"	//
			"vmlal.u32 q9 , d5, d6[0] 		\n\t"	//
			"umaal r9, r12, r5, r7			\n\t"//
			"veor q4,q4,q4	 				\n\t"
			"veor q5,q5,q5	 				\n\t"
			"umaal r9, r14, r4, r8			\n\t"//
			"veor q6,q6,q6	 				\n\t"
			"veor q7,q7,q7	 				\n\t"
			"str r9, [r0, #4 * 23]			\n\t"//ST IN
			"veor q8,q8,q8	 				\n\t"
			"veor q15,q15,q15				\n\t"
			"ldr r9, [r0, #4 * 24]			\n\t"//LD IN
			"vtrn.32 q14, q4 				\n\t"
			"vtrn.32 q13, q5 				\n\t"
			"ldr r4, [r1, #4 * 15]			\n\t"//LD M
			"vtrn.32 q12, q6 				\n\t"
			"vtrn.32 q11, q7 				\n\t"
			"umaal r9, r11, r4, r6			\n\t"//
			"vtrn.32 q10, q8 				\n\t"
			"vtrn.32 q9 , q15 				\n\t"
			"umaal r9, r12, r3, r7			\n\t"//
			"vadd.i64  q13, q13, q4 		\n\t"	//q
			"vadd.i64  q12, q12, q5 		\n\t"	//q
			"umaal r9, r14, r5, r8			\n\t"//
			"vadd.i64  q11, q11, q6 		\n\t"	//q
			"vadd.i64  q10, q10, q7 		\n\t"	//q
			"str r9, [r0, #4 * 24]			\n\t"//ST IN			
			"vadd.i64  q9, q9, q8 			\n\t"	//q
			"vqadd.u64 d30, d30, d29 		\n\t"	//d
			"ldr r9, [r0, #4 * 25]			\n\t"//LD IN
			"vext.8 d6, d6, d28, #4 		\n\t"
			"vmlal.u32 q13, d0, d6[0] 		\n\t"	//
			"ldr r5, [r1, #4 * 16]			\n\t"//LD M
			"vmlal.u32 q12, d3, d6[0] 		\n\t"	//
			"vmlal.u32 q11, d1, d6[0] 		\n\t"	//
			"umaal r9, r11, r5, r6			\n\t"//
			"vmlal.u32 q10, d4, d6[0] 		\n\t"	//
			"vmlal.u32 q9 , d2, d6[0] 		\n\t"	//
			"umaal r9, r12, r4, r7			\n\t"//
			"vmlal.u32 q15, d5, d6[0] 		\n\t"
			"veor q4,q4,q4	 				\n\t"
			"umaal r9, r14, r3, r8			\n\t"//
			"veor q5,q5,q5	 				\n\t"
			"veor q6,q6,q6	 				\n\t"
			"str r9, [r0, #4 * 25]			\n\t"//ST IN
			"veor q7,q7,q7	 				\n\t"
			"veor q8,q8,q8	 				\n\t"
			"ldr r9, [r0, #4 * 26]			\n\t"//LD IN
			"veor q14,q14,q14					\n\t"
			"vtrn.32 q13, q4 				\n\t"
			"ldr r3, [r1, #4 * 17]			\n\t"//LD M
			"vtrn.32 q12, q5 				\n\t"
			"vtrn.32 q11, q6 				\n\t"
			"umaal r9, r11, r3, r6			\n\t"//
			"vtrn.32 q10, q7 				\n\t"
			"vtrn.32 q9 , q8 				\n\t"
			"umaal r9, r12, r5, r7			\n\t"//
			"vtrn.32 q15, q14 				\n\t"
			"vadd.i64  q12, q12, q4 		\n\t"	//q
			"umaal r9, r14, r4, r8			\n\t"//
			"vadd.i64  q11, q11, q5 		\n\t"	//q
			"vadd.i64  q10, q10, q6 		\n\t"	//q
			"str r9, [r0, #4 * 26]			\n\t"//ST IN
			"vadd.i64  q9 , q9 , q7 		\n\t"	//q
			"vadd.i64  q15, q15, q8 		\n\t"
			"ldr r9, [r0, #4 * 27]			\n\t"//LD IN
			"vqadd.u64 d28, d28, d27 		\n\t"	//d
			"vext.8 d6, d6, d26, #4 		\n\t"
			"ldr r4, [r1, #4 * 18]			\n\t"//LD M
			"vstr.64 d6, [sp, #4*8]			\n\t"
			"vmlal.u32 q12, d0, d7[0] 		\n\t"	//
			"umaal r9, r11, r4, r6			\n\t"//
			"vmlal.u32 q11, d3, d7[0] 		\n\t"	//
			"vmlal.u32 q10, d1, d7[0] 		\n\t"	//
			"umaal r9, r12, r3, r7			\n\t"//
			"vmlal.u32 q9 , d4, d7[0] 		\n\t"	//
			"vmlal.u32 q15, d2, d7[0] 		\n\t"	//
			"umaal r9, r14, r5, r8			\n\t"//
			"vmlal.u32 q14, d5, d7[0] 		\n\t"
			"veor q4,q4,q4	 				\n\t"
			"str r9, [r0, #4 * 27]			\n\t"//ST IN
			"veor q5,q5,q5	 				\n\t"
			"veor q6,q6,q6	 				\n\t"
			"ldr r9, [r0, #4 * 28]			\n\t"//LD IN
			"veor q7,q7,q7	 				\n\t"
			"veor q8,q8,q8	 				\n\t"
			"ldr r5, [r1, #4 * 19]			\n\t"//LD M
			"veor q13,q13,q13				\n\t"
			"vtrn.32 q12, q4 				\n\t"
			"umaal r9, r11, r5, r6			\n\t"//
			"vtrn.32 q11, q5 				\n\t"
			"vtrn.32 q10, q6 				\n\t"
			"umaal r9, r12, r4, r7			\n\t"//
			"vtrn.32 q9 , q7 				\n\t"
			"vtrn.32 q15, q8 				\n\t"
			"umaal r9, r14, r3, r8			\n\t"//
			"vtrn.32 q14, q13 				\n\t"
			"vadd.i64  q11, q11, q4 		\n\t"	//q
			"str r9, [r0, #4 * 28]			\n\t"//ST IN
			"vadd.i64  q10, q10, q5 		\n\t"	//q
			"vadd.i64  q9 , q9 , q6 		\n\t"	//q
			"ldr r9, [r0, #4 * 29]			\n\t"//LD IN
			"vadd.i64  q15, q15, q7 		\n\t"	//q
			"vadd.i64  q14, q14, q8 		\n\t"
			"ldr r3, [r1, #4 * 20]			\n\t"//LD M
			"vqadd.u64 d26, d26, d25 		\n\t"	//d
			"vext.8 d7, d7, d24, #4 		\n\t"
			"umaal r9, r11, r3, r6			\n\t"//
			"umaal r9, r12, r5, r7			\n\t"//
			"umaal r9, r14, r4, r8			\n\t"//
			"str r9, [r0, #4 * 29]			\n\t"//ST IN
			"ldr r9, [r0, #4 * 30]			\n\t"//LD IN
			"ldr r4, [r1, #4 * 21]			\n\t"//LD M
			"umaal r9, r11, r4, r6			\n\t"//
			"umaal r9, r12, r3, r7			\n\t"//
			"umaal r9, r14, r5, r8			\n\t"//
			"str r9, [r0, #4 * 30]			\n\t"//ST IN
			"ldr r9, [r0, #4 * 31]			\n\t"//LD IN
			"ldr r5, [r1, #4 * 22]			\n\t"//LD M
			"umaal r9, r11, r5, r6			\n\t"//
			"umaal r9, r12, r4, r7			\n\t"//
			"umaal r9, r14, r3, r8			\n\t"//
			"str r9, [r0, #4 * 31]			\n\t"//ST IN
			"ldr r9, [r0, #4 * 32]			\n\t"//LD IN
			"ldr r3, [r1, #4 * 23]			\n\t"//LD M
			"umaal r9, r11, r3, r6			\n\t"//
			"umaal r9, r12, r5, r7			\n\t"//
			"umaal r9, r14, r4, r8			\n\t"//
			"str r9, [r0, #4 * 32]			\n\t"//ST IN		
			"umaal r11, r12, r3, r7			\n\t"//
			"umaal r11, r14, r5, r8			\n\t"//
			"str r11, [r0, #4 * 33]			\n\t"//ST IN
			"umaal r12, r14, r3, r8			\n\t"//
			"str r12, [r0, #4 * 34]			\n\t"//ST IN
			"str r14, [r0, #4 * 35]			\n\t"//ST IN

			
			//////////////////////////last
			"ldr r3, [r1, #4 * 11]			\n\t"//M constant


			"ldr r9,  [r0, #4 * 23]			\n\t"//LD IN
			"mov r10, #0					\n\t"


			"ldr r4, [r0, #4 * 12]			\n\t"//Q
			"umlal r9, r10, r3, r4			\n\t"//#1
			"vmov d7[0], r9 				\n\t"//operand

			"str r9,  [r0, #4 * 23]			\n\t"//LD IN

			////////////////////////TEST START
			//"str r9,  [r2, #4 * 0]			\n\t"//ST IN0
			//"ldr r9,  [sp, #4 * 0]		\n\t"//LD IN
			//"str r9,  [r2, #4 * 0]			\n\t"//ST IN0

			//"ldr r9,  [sp, #4 * 1]		\n\t"//LD IN
			//"str r9,  [r2, #4 * 1]			\n\t"//ST IN0

			//"ldr r9,  [sp, #4 * 2]		\n\t"//LD IN
			//"str r9,  [r2, #4 * 2]			\n\t"//ST IN0

			////////////////////////TEST END


			/////////////////////
			"vmlal.u32 q11, d0, d7[0] 		\n\t"	//
			"vmlal.u32 q10, d3, d7[0] 		\n\t"	//
			"vmlal.u32 q9 , d1, d7[0] 		\n\t"	//
			"vmlal.u32 q15, d4, d7[0] 		\n\t"	//
			"vmlal.u32 q14, d2, d7[0] 		\n\t"	//
			"vmlal.u32 q13, d5, d7[0] 		\n\t"	//

			"veor q4,q4,q4	 				\n\t"
			"veor q5,q5,q5	 				\n\t"
			"veor q6,q6,q6	 				\n\t"
			"veor q7,q7,q7	 				\n\t"
			"veor q8,q8,q8	 				\n\t"
			"veor q12,q12,q12				\n\t"

			"vtrn.32 q11, q4 				\n\t"
			"vtrn.32 q10, q5 				\n\t"
			"vtrn.32 q9 , q6 				\n\t"
			"vtrn.32 q15, q7 				\n\t"
			"vtrn.32 q14, q8 				\n\t"
			"vtrn.32 q13, q12 				\n\t"

			"vadd.i64  q10, q10, q4 		\n\t"	//q
			"vadd.i64  q9 , q9 , q5 		\n\t"	//q
			"vadd.i64  q15, q15, q6 		\n\t"	//q
			"vadd.i64  q14, q14, q7 		\n\t"	//q
			"vadd.i64  q13, q13, q8 		\n\t"	//q

			"vqadd.u64 d24, d24, d23 		\n\t"	//d
			"vext.8 d7, d7, d22, #4 		\n\t"

			"vstr.64 d7, [sp, #4*10]		\n\t"	//saving





			"vldr.64 d0, [r0, #4*36]		\n\t"
			"vldr.64 d1, [r0, #4*38]		\n\t"
			"vldr.64 d2, [r0, #4*40]		\n\t"		
			"vldr.64 d3, [r0, #4*42]		\n\t"
			"vldr.64 d4, [r0, #4*44]		\n\t"
			"vldr.64 d5, [r0, #4*46]		\n\t"

			"vtrn.32 d0, d3 				\n\t"
			"vtrn.32 d1, d4 				\n\t"
			"vtrn.32 d2, d5 				\n\t"

			"vmov.i64     q3, 0xFFFFFFFF    \n\t"// mast 1
			"vshr.u64     q3, q3, #31		\n\t"


			"vmlal.u32 q10, d0, d7[0] 		\n\t"	//
			"vmlal.u32 q9 , d3, d7[0] 		\n\t"	//
			"vmlal.u32 q15, d1, d7[0] 		\n\t"	//
			"vmlal.u32 q14, d4, d7[0] 		\n\t"	//
			"vmlal.u32 q13, d2, d7[0] 		\n\t"	//
			"vmlal.u32 q12, d5, d7[0] 		\n\t"	//






			//2
			"ldr r4, [r0, #4 * 13]			\n\t"//Q
			"ldr r9,  [r0, #4 * 24]			\n\t"//LD IN
			"ldr r11,  [sp, #4 * 0]			\n\t"//LD IN
			"mov r12, #0					\n\t"

			"umaal r9, r10, r3, r4			\n\t"//#2

			"adds r9, r9, r11				\n\t"
			"adcs r10, r10, #0				\n\t"
			"adcs r12, r12, #0				\n\t"

			"str r9,  [r2, #4 * 0]			\n\t"//ST IN


			//3: 10-12
			"ldr r4, [r0, #4 * 14]			\n\t"//Q
			"ldr r9,  [r0, #4 * 25]			\n\t"//LD IN
			"ldr r11,  [sp, #4 * 1]			\n\t"//LD IN
			"mov r14, #0					\n\t"

			"umaal r9, r10, r3, r4			\n\t"//#2

			"adds r9, r9, r11				\n\t"
			"adcs r10, r10, r12				\n\t"
			"adcs r12, r14, #0				\n\t"

			"str r9,  [r2, #4 * 1]			\n\t"//ST IN




			//4: 10-12


	

			"ldr r4, [r0, #4 * 15]			\n\t"//Q
			"ldr r9,  [r0, #4 * 26]			\n\t"//LD IN
			"ldr r11,  [sp, #4 * 2]			\n\t"//LD IN
			"mov r14, #0					\n\t"

			"umaal r9, r10, r3, r4			\n\t"//#2

			"adds r9, r9, r11				\n\t"
			"adcs r10, r10, r12				\n\t"
			"adcs r12, r14, #0				\n\t"

			"str r9,  [r2, #4 * 2]			\n\t"//ST IN





			//5: 10-12
			"ldr r4, [r0, #4 * 16]			\n\t"//Q
			"ldr r9,  [r0, #4 * 27]			\n\t"//LD IN
			"ldr r11,  [sp, #4 * 3]			\n\t"//LD IN
			"mov r14, #0					\n\t"

			"umaal r9, r10, r3, r4			\n\t"//#2

			"adds r9, r9, r11				\n\t"
			"adcs r10, r10, r12				\n\t"
			"adcs r12, r14, #0				\n\t"

			"str r9,  [r2, #4 * 3]			\n\t"//ST IN

			//6: 10-12
			"ldr r4, [r0, #4 * 17]			\n\t"//Q
			"ldr r9,  [r0, #4 * 28]			\n\t"//LD IN
			"ldr r11,  [sp, #4 * 4]			\n\t"//LD IN
			"mov r14, #0					\n\t"

			"umaal r9, r10, r3, r4			\n\t"//#2

			"adds r9, r9, r11				\n\t"
			"adcs r10, r10, r12				\n\t"
			"adcs r12, r14, #0				\n\t"

			"str r9,  [r2, #4 * 4]			\n\t"//ST IN			
			
			//7: 10-12
			"ldr r4, [r0, #4 * 18]			\n\t"//Q
			"ldr r9,  [r0, #4 * 29]			\n\t"//LD IN
			"ldr r11,  [sp, #4 * 5]			\n\t"//LD IN
			"mov r14, #0					\n\t"

			"umaal r9, r10, r3, r4			\n\t"//#2

			"adds r9, r9, r11				\n\t"
			"adcs r10, r10, r12				\n\t"
			"adcs r12, r14, #0				\n\t"

			"str r9,  [r2, #4 * 5]			\n\t"//ST IN			

			//8: 10-12
			"ldr r4, [r0, #4 * 19]			\n\t"//Q
			"ldr r9,  [r0, #4 * 30]			\n\t"//LD IN
			"ldr r11,  [sp, #4 * 6]			\n\t"//LD IN
			"mov r14, #0					\n\t"

			"umaal r9, r10, r3, r4			\n\t"//#2

			"adds r9, r9, r11				\n\t"
			"adcs r10, r10, r12				\n\t"
			"adcs r12, r14, #0				\n\t"

			"str r9,  [r2, #4 * 6]			\n\t"//ST IN
			




			//9: 10-12
			"ldr r4, [r0, #4 * 20]			\n\t"//Q
			"ldr r9,  [r0, #4 * 31]			\n\t"//LD IN
			"ldr r11,  [sp, #4 * 7]			\n\t"//LD IN
			"mov r14, #0					\n\t"

			"umaal r9, r10, r3, r4			\n\t"//#2

			"adds r9, r9, r11				\n\t"
			"adcs r10, r10, r12				\n\t"
			"adcs r12, r14, #0				\n\t"

			"str r9,  [r2, #4 * 7]			\n\t"//ST IN	


			//10: 10-12
			"ldr r4, [r0, #4 * 21]			\n\t"//Q
			"ldr r9,  [r0, #4 * 32]			\n\t"//LD IN
			"ldr r11,  [sp, #4 * 8]			\n\t"//LD IN
			"mov r14, #0					\n\t"

			"umaal r9, r10, r3, r4			\n\t"//#2

			"adds r9, r9, r11				\n\t"
			"adcs r10, r10, r12				\n\t"
			"adcs r12, r14, #0				\n\t"

			"str r9,  [r2, #4 * 8]			\n\t"//ST IN

			//11: 10-12
			"ldr r4, [r0, #4 * 22]			\n\t"//Q
			"ldr r9,  [r0, #4 * 33]			\n\t"//LD IN
			"ldr r11,  [sp, #4 * 9]			\n\t"//LD IN
			"mov r14, #0					\n\t"

			"umaal r9, r10, r3, r4			\n\t"//#2

			"adds r9, r9, r11				\n\t"
			"adcs r10, r10, r12				\n\t"
			"adcs r12, r14, #0				\n\t"

			"str r9,  [r2, #4 * 9]			\n\t"//ST IN



			////////////////////////TEST END
//			"ldr r9,  [sp, #4 * 11]		\n\t"//LD IN
//			"str r9,  [r2, #4 * 13]			\n\t"//ST IN0
//			"str r10,  [r2, #4 * 11]			\n\t"//ST IN0
//			"ldr r9,  [r0, #4 * 34]			\n\t"//LD IN
//			"str r9,  [r2, #4 * 12]			\n\t"//ST IN0
			//"ldr r9,  [r0, #4 * 35]			\n\t"//LD IN
			//"str r9,  [r2, #4 * 13]			\n\t"//ST IN0
//			"ldr r9,  [sp, #4 * 10]			\n\t"//LD IN
//			"str r9,  [r2, #4 * 14]			\n\t"//ST IN0

//			"ldr r9,  [r0, #4 * 35]		\n\t"//LD IN
//			"str r9,  [r2, #4 * 14]			\n\t"//ST IN0

			////////////////////////TEST END



			//12: 10-12
			"ldr r4, [r0, #4 * 23]			\n\t"//Q
			"ldr r9,  [r0, #4 * 34]			\n\t"//LD IN
			"ldr r11,  [sp, #4 * 10]			\n\t"//LD IN
			"mov r14, #0					\n\t"






			"umaal r9, r10, r3, r4			\n\t"//#2



			"adds r9, r9, r11				\n\t"
			"adcs r10, r10, r12				\n\t"
			"adcs r12, r14, #0				\n\t"

			"str r9,  [r2, #4 * 10]			\n\t"//ST IN



			"ldr r9,  [r0, #4 * 35]			\n\t"//LD IN
			"ldr r11,  [sp, #4 * 11]			\n\t"//LD IN
			
			"adds r9, r9, r10				\n\t"
			"adcs r14, r14, #0				\n\t"

			"adds r9, r9, r11				\n\t"
			"adcs r14, r14, #0				\n\t"

			"str r9,  [r2, #4 * 11]			\n\t"//ST IN




			//TEST
			//"vstr.64 d20, [r2, #4*14]		\n\t"	//saving
			//

			//r14
			"mov r0, r2						\n\t"

			"vmov r1, r2, d20 				\n\t"
			"vmov r3, r4, d18 				\n\t"
			"vmov r5, r6, d30 				\n\t"
			"vmov r7, r8, d28 				\n\t"
			"vmov r9, r10, d26 				\n\t"
			"vmov r11, r12, d24 			\n\t"

			"adds r1, r1, r14				\n\t"
			"adcs r2, r2, r3				\n\t"
			"adcs r3, r4, r5				\n\t"
			"adcs r4, r6, r7				\n\t"
			"adcs r5, r8, r9				\n\t"
			"adcs r6, r10, r11				\n\t"
			"adcs r12, r12, #0				\n\t"

			"str r1,  [r0, #4 * 12]			\n\t"//ST IN
			"str r2,  [r0, #4 * 13]			\n\t"//ST IN
			"str r3,  [r0, #4 * 14]			\n\t"//ST IN
			"str r4,  [r0, #4 * 15]			\n\t"//ST IN
			"str r5,  [r0, #4 * 16]			\n\t"//ST IN
			"str r6,  [r0, #4 * 17]			\n\t"//ST IN


			"vmov r1, r2, d21 				\n\t"
			"vmov r3, r4, d19 				\n\t"
			"vmov r5, r6, d31 				\n\t"
			"vmov r7, r8, d29 				\n\t"
			"vmov r9, r10, d27 				\n\t"
			"vmov r11,		d25[0] 			\n\t"

			"adds r1, r1, r12				\n\t"
			"adcs r2, r2, r3				\n\t"
			"adcs r3, r4, r5				\n\t"
			"adcs r4, r6, r7				\n\t"
			"adcs r5, r8, r9				\n\t"
			"adcs r6, r10, r11				\n\t"

			"str r1,  [r0, #4 * 18]			\n\t"//ST IN
			"str r2,  [r0, #4 * 19]			\n\t"//ST IN
			"str r3,  [r0, #4 * 20]			\n\t"//ST IN
			"str r4,  [r0, #4 * 21]			\n\t"//ST IN
			"str r5,  [r0, #4 * 22]			\n\t"//ST IN
			"str r6,  [r0, #4 * 23]			\n\t"//ST IN
/**/



			//////////////////////////////////////
			"add sp, #4 * 12				\n\t"


			"vpop {q4-q7}					\n\t"
			"pop  {r3-r11,pc}				\n\t"
			
	:
	:
	:
	);
}

unsigned int modulo_p503[]={0x0
,  0x0
,  0x0
,  0x0
,  0x0
,  0x0
,  0x0
,  0x0
,  0x0
,  0x0
,  0x0
,  0xEEB00000
,  0x49F878A8
,  0xE3EC9685
,  0x13F7CC76
,  0xDA959B1A
,  0xD6EBE876
,  0x84E9867
,  0x5CB25748
,  0x8562B504
,  0x97BADC66
,  0xE12909F
,  0xD541F71C
,  0x6FE5};



void fpmul_mont(const felm_t ma, const felm_t mb, felm_t mc)
{ // Multiprecision multiplication, c = a*b mod p.
    dfelm_t temp = {0};
    int k;
    //mp_mul(ma, mb, temp, NWORDS_FIELD);
    MUL512(ma, mb, temp);
	RED512(temp,modulo_p503,mc);

    //rdc_mont(temp, mc);
}


void __attribute__ ((noinline, naked)) SQR512(const felm_t ma, felm_t mc){
	asm(
			
			"push  {r4-r11,lr}			\n\t"
			"vpush {q4-q7}				\n\t"
			//256-bit load

			////////////////////////////////////////////////////////
			"vldmia r0,  {q0-q3} 		\n\t"
			"vtrn.32 q0, q1 			\n\t"
			"ldr r2, [r0, #4 * 8]		\n\t"//A0//higher part
			"vshr.u64 q14, q2, #63  	\n\t"
			"vshr.u64 q4 , q3, #63  	\n\t"
			"ldr r3, [r0, #4 * 9]		\n\t"//A1
			"vshl.u64 q2, q2, #1  		\n\t"
			"vshl.u64 q3, q3, #1  		\n\t"
			"ldr r4, [r0, #4 * 10]			\n\t"//A2
			"vqadd.u64 d5, d5, d28 		\n\t"//d
			"vqadd.u64 d6, d6, d29 		\n\t"//d
			"ldr r5, [r0, #4 * 13]			\n\t"//A5
			"vqadd.u64 d7, d7, d8 		\n\t"//d	//carry in d9
			"vtrn.32 q2, q3 			\n\t "			//shuffle 6 2 4 0
			"ldr r6, [r0, #4 * 14]			\n\t"//A6
			"vmull.u32 q13, d0, d0[0] 	\n\t"	//A4A0 X B0B0
			"vmull.u32 q12, d2, d0[0] 	\n\t"	//A5A1 X B0B0
			"ldr r7, [r0, #4 * 15]			\n\t"//A7
			"vmull.u32 q11, d1, d0[0] 	\n\t"	//A6A2 X B0B0
			"vmull.u32 q10, d3, d0[0] 	\n\t"	//A7A3 X B0B0
			"umull r8, r9, r2, r5			\n\t"//05
			"veor q6,q6,q6 				\n\t"
			"veor q7,q7,q7 				\n\t"
			"str r8, [r1, #4 * 21]			\n\t"//C5
			"veor q8,q8,q8 				\n\t"
			"veor q9,q9,q9 				\n\t"
			"mov r8, #0						\n\t"
			"vtrn.32 q13, q6 			\n\t"
			"vtrn.32 q12, q7 			\n\t"
			"umaal r8, r9, r2, r6			\n\t"//06
			"vtrn.32 q11, q8 			\n\t"
			"vtrn.32 q10, q9 			\n\t"
			"str r8, [r1, #4 * 22]			\n\t"//C6
			"vadd.i64  q12, q12, q6 	\n\t"	//q
			"vadd.i64  q11, q11, q7 	\n\t"	//q
			"mov r8, #0						\n\t"
			"vadd.i64  q10, q10, q8 	\n\t"	//q
			"vqadd.u64 d18, d18, d27 	\n\t"	//d
			"umaal r8, r9, r2, r7			\n\t"//07
			"vext.8 d28, d28, d26, #4 	\n\t"
			"vmlal.u32 q12, d0, d2[0] 	\n\t"	//A4A0 X B0B0
			"mov r10, #0					\n\t"
			"vmlal.u32 q11, d2, d2[0] 	\n\t"	//A5A1 X B0B0
			"vmlal.u32 q10, d1, d2[0] 	\n\t"	//A6A2 X B0B0
			"umaal r8, r10, r3, r6			\n\t"//16
			"vmlal.u32 q9 , d3, d2[0] 	\n\t"	//A7A3 X B0B0
			"veor q13,q13,q13			\n\t"
			"str r8, [r1, #4 * 23]			\n\t"//C7
			"veor q6 ,q6 ,q6 				\n\t"
			"veor q7 ,q7 ,q7 				\n\t"
			"umaal r9, r10, r3, r7			\n\t"//17
			"veor q8 ,q8 ,q8 				\n\t"
			"vtrn.32 q12, q13 			\n\t"
			"str r9, [r1, #4 * 24]			\n\t"//C8
			"vtrn.32 q11, q6 			\n\t"
			"vtrn.32 q10, q7 			\n\t"
			"mov r9, #0						\n\t"
			"vtrn.32 q9 , q8 			\n\t"
			"vadd.i64  q11, q11, q13 	\n\t"	//q
			"umaal r9, r10, r4, r7			\n\t"//27
			"vadd.i64  q10, q10, q6 	\n\t"	//q
			"vadd.i64  q9 , q9 , q7 	\n\t"	//q
			"str r9, [r1, #4 * 25]			\n\t"//C9
			"vqadd.u64 d16, d16, d25 	\n\t"	//d
			"vext.8 d28, d28, d24, #4 	\n\t"
			"str r10, [r1, #4 * 26]			\n\t"//C10
			"vmlal.u32 q11, d0, d1[0] 	\n\t"	//A4A0 X B0B0
			"vmlal.u32 q10, d2, d1[0] 	\n\t"	//A5A1 X B0B0
			"umull r5, r6, r2, r2			\n\t"//00
			"vmlal.u32 q9 , d1, d1[0] 	\n\t"	//A6A2 X B0B0
			"vmlal.u32 q8 , d3, d1[0] 	\n\t"	//A7A3 X B0B0
			"umull r7, r8, r2, r3			\n\t"//01
			"veor q12,q12,q12			\n\t"
			"veor q13,q13,q13			\n\t"
			"mov r9, #0						\n\t"
			"veor q6 ,q6 ,q6 			\n\t"
			"veor q7 ,q7 ,q7 			\n\t"
			"adds r7, r7, r7				\n\t"
			"vtrn.32 q11, q12 			\n\t"
			"vtrn.32 q10, q13			\n\t"
			"adcs r8, r8, r8				\n\t"
			"vtrn.32 q9 , q6			\n\t"
			"vtrn.32 q8 , q7 			\n\t"
			"adcs r9, r9, #0				\n\t"
			"vadd.i64  q10, q10, q12 	\n\t"	//q
			"vadd.i64  q9 , q9 , q13	\n\t"	//q
			"adds r7, r7, r6				\n\t"
			"vadd.i64  q8 , q8 , q6 	\n\t"	//q
			"vqadd.u64 d14, d14, d23 	\n\t"	//d
			"adcs r11, r8, #0				\n\t"
			"vext.8 d29, d29, d22, #4 	\n\t"
			"vmlal.u32 q10, d0, d3[0] 	\n\t"	//A4A0 X B0B0
			"adcs r12, r9, #0				\n\t"
			"vmlal.u32 q9 , d2, d3[0] 	\n\t"	//A5A1 X B0B0
			"vmlal.u32 q8 , d1, d3[0] 	\n\t"	//A6A2 X B0B0
			"str r5, [r1, #4 * 16]			\n\t"//C0
			"vmlal.u32 q7 , d3, d3[0] 	\n\t"	//A7A3 X B0B0
			"veor q11,q11,q11			\n\t"
			"str r7, [r1, #4 * 17]			\n\t"//C1
			"veor q12,q12,q12			\n\t"
			"veor q13,q13,q13			\n\t"
			"ldr r5, [r0, #4 * 11]			\n\t"//A3
			"veor q6 ,q6 ,q6 			\n\t"
			"vtrn.32 q10, q11			\n\t"
			"umull r6, r7, r2, r4			\n\t"//02
			"vtrn.32 q9 , q12			\n\t"
			"vtrn.32 q8 , q13			\n\t"
			"mov r8, #0						\n\t"
			"vtrn.32 q7 , q6 			\n\t"
			"vadd.i64  q9 , q9 , q11	\n\t"	//q
			"umaal r7, r8, r2, r5			\n\t"//03
			"vadd.i64  q8 , q8 , q12	\n\t"	//q
			"vadd.i64  q7 , q7 , q13 	\n\t"	//q
			"mov r9, #0						\n\t"
			"vqadd.u64 d12, d12, d21 	\n\t"	//d
			"vext.8 d29, d29, d20, #4 	\n\t"
			"umaal r7, r9, r3, r4			\n\t"//12
			"vmlal.u32 q9 , d0, d0[1] 	\n\t"	//A4A0 X B0B0
			"vmlal.u32 q8 , d2, d0[1] 	\n\t"	//A5A1 X B0B0
			"mov r10, #0					\n\t"
			"vmlal.u32 q7 , d1, d0[1] 	\n\t"	//A6A2 X B0B0
			"vmlal.u32 q6 , d3, d0[1] 	\n\t"	//A7A3 X B0B0
			"adds r6, r6, r6				\n\t"
			"veor q10,q10,q10			\n\t"
			"veor q11,q11,q11			\n\t"
			"adcs r7, r7, r7				\n\t"
			"veor q12,q12,q12			\n\t"
			"veor q13,q13,q13			\n\t"
			"adcs r10, r10, #0				\n\t"
			"vtrn.32 q9 , q10			\n\t"
			"vtrn.32 q8 , q11			\n\t"
			"adds r6, r6, r11				\n\t"
			"vtrn.32 q7 , q12			\n\t"
			"vtrn.32 q6 , q13			\n\t"
			"adcs r7, r7, r12				\n\t"
			"vadd.i64  q8 , q8 , q10	\n\t"	//q
			"vadd.i64  q7 , q7 , q11 	\n\t"	//q
			"adcs r11, r10, #0				\n\t"
			"vadd.i64  q6 , q6 , q12	\n\t"	//q
			"vqadd.u64 d26, d26, d19 	\n\t"	//d
			"mov r12, #0					\n\t"
			"vext.8 d30, d30, d18, #4 	\n\t"
			"vmlal.u32 q8 , d0, d2[1] 	\n\t"	//A4A0 X B0B0
			"umaal r6, r12, r3, r3			\n\t"//11
			"vmlal.u32 q7 , d2, d2[1] 	\n\t"	//A5A1 X B0B0
			"vmlal.u32 q6 , d1, d2[1] 	\n\t"	//A6A2 X B0B0
			"adds r7, r7, r12				\n\t"
			"vmlal.u32 q13, d3, d2[1] 	\n\t"	//A7A3 X B0B0
			"veor q9 ,q9 ,q9			\n\t"
			"adcs r11, r10, #0				\n\t"
			"veor q10,q10,q10			\n\t"
			"veor q11,q11,q11			\n\t"
			"str r6, [r1, #4 * 18]			\n\t"//C2
			"veor q12,q12,q12			\n\t"
			"vtrn.32 q8 , q9			\n\t"
			"str r7, [r1, #4 * 19]			\n\t"//C3
			"vtrn.32 q7 , q10			\n\t"
			"vtrn.32 q6 , q11			\n\t"
			"ldr r7, [r0, #4 * 12]			\n\t"//A4
			"vtrn.32 q13, q12			\n\t"
			"vadd.i64  q7 , q7 , q9 	\n\t"	//q
			"mov r6, #0						\n\t"
			"vadd.i64  q6 , q6 , q10	\n\t"	//q
			"vadd.i64  q13, q13, q11	\n\t"	//q
			"umaal r6, r8, r2, r7			\n\t"//04
			"vqadd.u64 d24, d24, d17 	\n\t"	//d
			"vext.8 d30, d30, d16, #4 	\n\t"
			"mov r2, r7						\n\t"
			"vmlal.u32 q7 , d0, d1[1] 	\n\t"	//A4A0 X B0B0
			"vmlal.u32 q6 , d2, d1[1] 	\n\t"	//A5A1 X B0B0
			"umaal r6, r9, r3, r5			\n\t"//13
			"vmlal.u32 q13, d1, d1[1] 	\n\t"	//A6A2 X B0B0
			"vmlal.u32 q12, d3, d1[1] 	\n\t"	//A7A3 X B0B0
			"ldr r7, [r1, #4 * 21]			\n\t"//C5
			"veor q8 ,q8 ,q8			\n\t"
			"veor q9 ,q9 ,q9			\n\t"
			"umaal r7, r8, r3, r2			\n\t"//14
			"veor q10,q10,q10			\n\t"
			"veor q11,q11,q11			\n\t"
			"umaal r7, r9, r4, r5			\n\t"//23
			"vtrn.32 q7 , q8			\n\t"
			"vtrn.32 q6 , q9			\n\t"
			"mov r10, #0					\n\t"
			"vtrn.32 q13, q10			\n\t"
			"vtrn.32 q12, q11			\n\t"
			"adds r6, r6, r6				\n\t"
			"vadd.i64  q6 , q6 , q8		\n\t"	//q
			"vadd.i64  q13, q13, q9		\n\t"	//q
			"adcs r7, r7, r7				\n\t"
			"vadd.i64  q12, q12, q10	\n\t"	//q
			"vqadd.u64 d22, d22, d15 	\n\t"	//d
			"adcs r10, r10, #0				\n\t"
			"vext.8 d31, d31, d14, #4 	\n\t"
			"vmlal.u32 q6 , d0, d3[1] 	\n\t"	//A4A0 X B0B0
			"adds r6, r6, r11				\n\t"
			"vmlal.u32 q13, d2, d3[1] 	\n\t"	//A5A1 X B0B0
			"vmlal.u32 q12, d1, d3[1] 	\n\t"	//A6A2 X B0B0
			"adcs r7, r7, #0				\n\t"
			"vmlal.u32 q11, d3, d3[1] 	\n\t"	//A7A3 X B0B0
			"veor q7, q7, q7			\n\t"
			"adcs r11, r10, #0				\n\t"
			"veor q8 ,q8 ,q8			\n\t"
			"veor q9 ,q9 ,q9			\n\t"
			"mov r12, #0					\n\t"
			"veor q10,q10,q10			\n\t"
			"vtrn.32 q6 , q7			\n\t"
			"umaal r6, r12, r4, r4			\n\t"//22
			"vtrn.32 q13, q8			\n\t"
			"vtrn.32 q12, q9			\n\t"
			"adds r7, r7, r12				\n\t"
			"vtrn.32 q11, q10			\n\t"
			"vadd.i64  q13, q13, q7		\n\t"	//q
			"adcs r11, r10, #0				\n\t"
			"vadd.i64  q12, q12, q8		\n\t"	//q
			"vadd.i64  q11, q11, q9		\n\t"	//q
			"str r6, [r1, #4 * 20]			\n\t"//C4
			"vqadd.u64 d20, d20, d13 	\n\t"	//d
			"vext.8 d31, d31, d12, #4 	\n\t"
			"str r7, [r1, #4 * 21]			\n\t"//C5
			"vstr.64 d28, [r1]			\n\t"
			"vstr.64 d29, [r1, #4*2]	\n\t"
			"ldr r7, [r0, #4 * 13]			\n\t"//A5
			"vstr.64 d30, [r1, #4*4]	\n\t"
			"vstr.64 d31, [r1, #4*6]	\n\t"
			"ldr r6, [r1, #4 * 22]			\n\t"//C6
			"vmlal.u32 q13, d4, d0[0] 	\n\t"	//A4A0 X B0B0
			"vmlal.u32 q12, d6, d0[0] 	\n\t"	//A5A1 X B0B0
			"umaal r6, r8, r3, r7			\n\t"//15
			"vmlal.u32 q11, d5, d0[0] 	\n\t"	//A6A2 X B0B0
			"vmlal.u32 q10, d7, d0[0] 	\n\t"	//A7A3 X B0B0
			"mov r3, r7						\n\t"
			"veor q6,q6,q6 				\n\t"
			"veor q7,q7,q7 				\n\t"
			"umaal r6, r9, r4, r2			\n\t"//24
			"veor q8,q8,q8 				\n\t"
			"veor q9,q9,q9 				\n\t"
			"ldr r7, [r1, #4 * 23]			\n\t"//C7
			"vtrn.32 q13, q6 			\n\t"
			"vtrn.32 q12, q7 			\n\t"
			"umaal r7, r8, r4, r3			\n\t"//25
			"vtrn.32 q11, q8 			\n\t"
			"vtrn.32 q10, q9 			\n\t"
			"umaal r7, r9, r5, r2			\n\t"//34
			"vadd.i64  q12, q12, q6 	\n\t"	//q
			"vadd.i64  q11, q11, q7 	\n\t"	//q
			"mov r10, #0					\n\t"
			"vadd.i64  q10, q10, q8 	\n\t"	//q
			"vqadd.u64 d18, d18, d27 	\n\t"	//d
			"adds r6, r6, r6				\n\t"
			"vext.8 d28, d28, d26, #4 	\n\t"
			"vmlal.u32 q12, d4, d2[0] 	\n\t"	//A4A0 X B0B0
			"adcs r7, r7, r7				\n\t"
			"vmlal.u32 q11, d6, d2[0] 	\n\t"	//A5A1 X B0B0
			"vmlal.u32 q10, d5, d2[0] 	\n\t"	//A6A2 X B0B0
			"adcs r10, r10, #0				\n\t"
			"vmlal.u32 q9 , d7, d2[0] 	\n\t"	//A7A3 X B0B0
			"veor q13,q13,q13			\n\t"
			"adds r6, r6, r11				\n\t"
			"veor q6 ,q6 ,q6 				\n\t"
			"veor q7 ,q7 ,q7 				\n\t"
			"adcs r7, r7, #0				\n\t"
			"veor q8 ,q8 ,q8 				\n\t"
			"vtrn.32 q12, q13 			\n\t"
			"adcs r11, r10, #0				\n\t"
			"vtrn.32 q11, q6 			\n\t"
			"vtrn.32 q10, q7 			\n\t"
			"mov r12, #0					\n\t"
			"vtrn.32 q9 , q8 			\n\t"
			"vadd.i64  q11, q11, q13 	\n\t"	//q
			"umaal r6, r12, r5, r5			\n\t"//33
			"vadd.i64  q10, q10, q6 	\n\t"	//q
			"vadd.i64  q9 , q9 , q7 	\n\t"	//q
			"adds r7, r7, r12				\n\t"
			"vqadd.u64 d16, d16, d25 	\n\t"	//d
			"vext.8 d28, d28, d24, #4 	\n\t"
			"adcs r11, r10, #0				\n\t"
			"vmlal.u32 q11, d4, d1[0] 	\n\t"	//A4A0 X B0B0
			"vmlal.u32 q10, d6, d1[0] 	\n\t"	//A5A1 X B0B0
			"str r6, [r1, #4 * 22]			\n\t"//C6
			"vmlal.u32 q9 , d5, d1[0] 	\n\t"	//A6A2 X B0B0
			"vmlal.u32 q8 , d7, d1[0] 	\n\t"	//A7A3 X B0B0
			"str r7, [r1, #4 * 23]			\n\t"//C7
			"veor q12,q12,q12			\n\t"
			"veor q13,q13,q13			\n\t"
			"ldr r7, [r0, #4 * 14]			\n\t"//A6
			"veor q6 ,q6 ,q6 			\n\t"
			"veor q7 ,q7 ,q7 			\n\t"
			"ldr r6, [r1, #4 * 24]			\n\t"//C8
			"vtrn.32 q11, q12 			\n\t"
			"vtrn.32 q10, q13			\n\t"
			"umaal r6, r8, r4, r7			\n\t"//26
			"vtrn.32 q9 , q6			\n\t"
			"vtrn.32 q8 , q7 			\n\t"
			"mov r4, r7						\n\t"
			"vadd.i64  q10, q10, q12 	\n\t"	//q
			"vadd.i64  q9 , q9 , q13	\n\t"	//q
			"umaal r6, r9, r5, r3			\n\t"//35
			"vadd.i64  q8 , q8 , q6 	\n\t"	//q
			"vqadd.u64 d14, d14, d23 	\n\t"	//d
			"ldr r7, [r1, #4 * 25]			\n\t"//C9
			"vext.8 d29, d29, d22, #4 	\n\t"
			"vmlal.u32 q10, d4, d3[0] 	\n\t"	//A4A0 X B0B0
			"umaal r7, r8, r5, r4			\n\t"//36
			"vmlal.u32 q9 , d6, d3[0] 	\n\t"	//A5A1 X B0B0
			"vmlal.u32 q8 , d5, d3[0] 	\n\t"	//A6A2 X B0B0
			"umaal r7, r9, r2, r3			\n\t"//45
			"vmlal.u32 q7 , d7, d3[0] 	\n\t"	//A7A3 X B0B0
			"veor q11,q11,q11			\n\t"
			"mov r10, #0					\n\t"
			"veor q12,q12,q12			\n\t"
			"veor q13,q13,q13			\n\t"
			"adds r6, r6, r6				\n\t"
			"veor q6 ,q6 ,q6 			\n\t"
			"vtrn.32 q10, q11			\n\t"
			"adcs r7, r7, r7				\n\t"
			"vtrn.32 q9 , q12			\n\t"
			"vtrn.32 q8 , q13			\n\t"
			"adcs r10, r10, #0				\n\t"
			"vtrn.32 q7 , q6 			\n\t"
			"vadd.i64  q9 , q9 , q11	\n\t"	//q
			"adds r6, r6, r11				\n\t"
			"vadd.i64  q8 , q8 , q12	\n\t"	//q
			"vadd.i64  q7 , q7 , q13 	\n\t"	//q
			"adcs r7, r7, #0				\n\t"
			"vqadd.u64 d12, d12, d21 	\n\t"	//d
			"vext.8 d29, d29, d20, #4 	\n\t"
			"adcs r11, r10, #0				\n\t"
			"vmlal.u32 q9 , d4, d0[1] 	\n\t"	//A4A0 X B0B0
			"vmlal.u32 q8 , d6, d0[1] 	\n\t"	//A5A1 X B0B0
			"mov r12, #0					\n\t"
			"vmlal.u32 q7 , d5, d0[1] 	\n\t"	//A6A2 X B0B0
			"vmlal.u32 q6 , d7, d0[1] 	\n\t"	//A7A3 X B0B0
			"umaal r6, r12, r2, r2			\n\t"//44
			"veor q10,q10,q10			\n\t"
			"veor q11,q11,q11			\n\t"
			"adds r7, r7, r12				\n\t"
			"veor q12,q12,q12			\n\t"
			"veor q13,q13,q13			\n\t"
			"adcs r11, r10, #0				\n\t"
			"vtrn.32 q9 , q10			\n\t"
			"vtrn.32 q8 , q11			\n\t"
			"str r6, [r1, #4 * 24]			\n\t"//C8
			"vtrn.32 q7 , q12			\n\t"
			"vtrn.32 q6 , q13			\n\t"
			"str r7, [r1, #4 * 25]			\n\t"//C9
			"vadd.i64  q8 , q8 , q10	\n\t"	//q
			"vadd.i64  q7 , q7 , q11 	\n\t"	//q
			"ldr r7, [r0, #4 * 15]			\n\t"//A7
			"vadd.i64  q6 , q6 , q12	\n\t"	//q
			"vqadd.u64 d26, d26, d19 	\n\t"	//d
			"ldr r6, [r1, #4 * 26]			\n\t"//C10
			"vext.8 d30, d30, d18, #4 	\n\t"
			"vmlal.u32 q8 , d4, d2[1] 	\n\t"	//A4A0 X B0B0
			"umaal r6, r8, r5, r7			\n\t"//37
			"vmlal.u32 q7 , d6, d2[1] 	\n\t"	//A5A1 X B0B0
			"vmlal.u32 q6 , d5, d2[1] 	\n\t"	//A6A2 X B0B0
			"mov r5, r7						\n\t"
			"vmlal.u32 q13, d7, d2[1] 	\n\t"	//A7A3 X B0B0
			"veor q9 ,q9 ,q9			\n\t"
			"umaal r6, r9, r2, r4			\n\t"//46
			"veor q10,q10,q10			\n\t"
			"veor q11,q11,q11			\n\t"
			"mov r7, #0						\n\t"
			"veor q12,q12,q12			\n\t"
			"vtrn.32 q8 , q9			\n\t"
			"umaal r7, r8, r2, r5			\n\t"//47
			"vtrn.32 q7 , q10			\n\t"
			"vtrn.32 q6 , q11			\n\t"
			"umaal r7, r9, r3, r4			\n\t"//56
			"vtrn.32 q13, q12			\n\t"
			"vadd.i64  q7 , q7 , q9 	\n\t"	//q
			"mov r10, #0					\n\t"
			"vadd.i64  q6 , q6 , q10	\n\t"	//q
			"vadd.i64  q13, q13, q11	\n\t"	//q
			"adds r6, r6, r6				\n\t"
			"vqadd.u64 d24, d24, d17 	\n\t"	//d
			"vext.8 d30, d30, d16, #4 	\n\t"
			"adcs r7, r7, r7				\n\t"
			"vmlal.u32 q7 , d4, d1[1] 	\n\t"	//A4A0 X B0B0
			"vmlal.u32 q6 , d6, d1[1] 	\n\t"	//A5A1 X B0B0
			"adcs r10, r10, #0				\n\t"
			"vmlal.u32 q13, d5, d1[1] 	\n\t"	//A6A2 X B0B0
			"vmlal.u32 q12, d7, d1[1] 	\n\t"	//A7A3 X B0B0
			"adds r6, r6, r11				\n\t"
			"veor q8 ,q8 ,q8			\n\t"
			"veor q9 ,q9 ,q9			\n\t"
			"adcs r7, r7, #0				\n\t"
			"veor q10,q10,q10			\n\t"
			"veor q11,q11,q11			\n\t"
			"adcs r11, r10, #0				\n\t"
			"vtrn.32 q7 , q8			\n\t"
			"vtrn.32 q6 , q9			\n\t"
			"mov r12, #0					\n\t"
			"vtrn.32 q13, q10			\n\t"
			"vtrn.32 q12, q11			\n\t"
			"umaal r6, r12, r3, r3			\n\t"//55
			"vadd.i64  q6 , q6 , q8		\n\t"	//q
			"vadd.i64  q13, q13, q9		\n\t"	//q
			"adds r7, r7, r12				\n\t"
			"vadd.i64  q12, q12, q10	\n\t"	//q
			"vqadd.u64 d22, d22, d15 	\n\t"	//d
			"adcs r11, r10, #0				\n\t"
			"vext.8 d31, d31, d14, #4 	\n\t"
			"vmlal.u32 q6 , d4, d3[1] 	\n\t"	//A4A0 X B0B0
			"str r6, [r1, #4 * 26]			\n\t"//C10
			"vmlal.u32 q13, d6, d3[1] 	\n\t"	//A5A1 X B0B0
			"vmlal.u32 q12, d5, d3[1] 	\n\t"	//A6A2 X B0B0
			"str r7, [r1, #4 * 27]			\n\t"//C11
			"vmlal.u32 q11, d7, d3[1] 	\n\t"	//A7A3 X B0B0
			"veor q7, q7, q7			\n\t"
			"umaal r8, r9, r3, r5			\n\t"//57
			"veor q8 ,q8 ,q8			\n\t"
			"veor q9 ,q9 ,q9			\n\t"
			"mov r10, #0					\n\t"
			"veor q10,q10,q10			\n\t"
			"vtrn.32 q6 , q7			\n\t"
			"umaal r9, r10, r4, r5			\n\t"//67
			"vtrn.32 q13, q8			\n\t"
			"vtrn.32 q12, q9			\n\t"
			"mov r2, #0						\n\t"
			"vtrn.32 q11, q10			\n\t"
			"vadd.i64  q13, q13, q7		\n\t"	//q
			"adds r8, r8, r8				\n\t"
			"vadd.i64  q12, q12, q8		\n\t"	//q
			"vadd.i64  q11, q11, q9		\n\t"	//q
			"adcs r9, r9, r9				\n\t"
			"vqadd.u64 d20, d20, d13 	\n\t"	//d
			"vext.8 d31, d31, d12, #4 	\n\t"
			"adcs r10, r10, r10				\n\t"
			"vstr.64 d28, [r1, #4*8]	\n\t"
			"vstr.64 d29, [r1, #4*10]	\n\t"
			"adcs r2, r2, #0				\n\t"
			"vstr.64 d30, [r1, #4*12]	\n\t"
			"vstr.64 d31, [r1, #4*14]	\n\t"
			"mov r12, #0					\n\t"
			"vmlal.u32 q13, d0, d9[0] 	\n\t"	//A4A0 X B0B0
			"vmlal.u32 q12, d2, d9[0] 	\n\t"	//A5A1 X B0B0
			"mov r14, #0					\n\t"
			"vmlal.u32 q11, d1, d9[0] 	\n\t"	//A6A2 X B0B0
			"vmlal.u32 q10, d3, d9[0] 	\n\t"	//A7A3 X B0B0
			"umaal r8, r12, r4, r4			\n\t"//66
			"umaal r10, r14, r5, r5			\n\t"//77
			"adds r8, r8, r11				\n\t"
			"adcs r9, r9, r12				\n\t"
			"adcs r10, r10, #0				\n\t"
			"adcs r2, r2, r14				\n\t"
			"str r8, [r1, #4 * 28]			\n\t"//C12
			"str r9, [r1, #4 * 29]			\n\t"//C13
			"str r10, [r1, #4 * 30]			\n\t"//C14
			"str r2, [r1, #4 * 31]			\n\t"//C15





			"add r0, r1, #4*16				\n\t"
			//2
			"vmov r1, r2, d26 				\n\t"
			"vmov r3, r4, d24 				\n\t"
			"vmov r5, r6, d22 				\n\t"
			"vmov r7, r8, d20 				\n\t"

			//"adds r1, r1, r14				\n\t"//carry bit
			"adds r2, r2, r3 				\n\t"
			"adcs r3, r4, r5 				\n\t"
			"adcs r4, r6, r7 				\n\t"

			"vmov r5, r6, d27 				\n\t"
			"vmov r9, r10, d25 				\n\t"

			"adcs r5, r8, r5 				\n\t"
			"adcs r6, r6, r9 				\n\t"

			"vmov r7, r8, d23 				\n\t"
			"vmov r9, r11, d21 				\n\t"	//r11 is not used

			"adcs r7, r7, r10 				\n\t"
			"adcs r8, r8, r9 				\n\t"
			"adcs r14, r11, #0				\n\t"

			"ldr r9, [r0, #4 * 0]			\n\t"
			"ldr r10, [r0, #4 * 1]			\n\t"
			"ldr r11, [r0, #4 * 2]			\n\t"
			"ldr r12, [r0, #4 * 3]			\n\t"
			
			"adds r1 , r9 , r1				\n\t"
			"adcs r2, r10, r2				\n\t"
			"adcs r3, r11, r3				\n\t"
			"adcs r4, r12, r4				\n\t"

			"ldr r9, [r0, #4 * 4]			\n\t"
			"ldr r10, [r0, #4 * 5]			\n\t"
			"ldr r11, [r0, #4 * 6]			\n\t"
			"ldr r12, [r0, #4 * 7]			\n\t"

			"adcs r5, r9, r5				\n\t"
			"adcs r6, r10, r6				\n\t"
			"adcs r7, r11, r7				\n\t"
			"adcs r8, r12, r8				\n\t"

			//"mov r14, #0					\n\t"
			"adcs r14, r14, #0				\n\t"//carry

			"stmia r0!, {r1-r8} 			\n\t"


			"ldmia r0, {r1-r8} 				\n\t"
			"adds r1, r1 , r14				\n\t"
			"adcs r2, r2 , #0				\n\t"
			"adcs r3, r3 , #0				\n\t"
			"adcs r4, r4 , #0				\n\t"

			"adcs r5, r5 , #0				\n\t"
			"adcs r6, r6 , #0				\n\t"
			"adcs r7, r7 , #0				\n\t"
			"adcs r8, r8 , #0				\n\t"

			"stmia r0!, {r1-r8} 			\n\t"

			

			"vpop {q4-q7}				\n\t"
			"pop  {r4-r11,pc}			\n\t"
			:
			:
			:
			);	
}


void fpsqr_mont(const felm_t ma, felm_t mc)
{ // Multiprecision squaring, c = a^2 mod p.
    dfelm_t temp = {0};

    //mp_mul(ma, ma, temp, NWORDS_FIELD);
    MUL512(ma,ma,temp);
    //rdc_mont(temp, mc);
	RED512(temp,modulo_p503,mc);
}


void fpinv_mont(felm_t a)
{ // Field inversion using Montgomery arithmetic, a = a^(-1)*R mod p.
    felm_t tt;

    fpcopy(a, tt);
    fpinv_chain_mont(tt);
    fpsqr_mont(tt, tt);
    fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, a);
}


void fp2copy(const f2elm_t a, f2elm_t c)
{ // Copy a GF(p^2) element, c = a.
    fpcopy(a[0], c[0]);
    fpcopy(a[1], c[1]);
}


void fp2zero(f2elm_t a)
{ // Zero a GF(p^2) element, a = 0.
    fpzero(a[0]);
    fpzero(a[1]);
}


void fp2neg(f2elm_t a)
{ // GF(p^2) negation, a = -a in GF(p^2).
    fpneg(a[0]);
    fpneg(a[1]);
}


__inline void fp2add(const f2elm_t a, const f2elm_t b, f2elm_t c)           
{ // GF(p^2) addition, c = a+b in GF(p^2).
    fpadd(a[0], b[0], c[0]);
    fpadd(a[1], b[1], c[1]);
}


__inline void fp2sub(const f2elm_t a, const f2elm_t b, f2elm_t c)          
{ // GF(p^2) subtraction, c = a-b in GF(p^2).
    fpsub(a[0], b[0], c[0]);
    fpsub(a[1], b[1], c[1]);
}


void fp2div2(const f2elm_t a, f2elm_t c)          
{ // GF(p^2) division by two, c = a/2  in GF(p^2).
    fpdiv2(a[0], c[0]);
    fpdiv2(a[1], c[1]);
}


void fp2correction(f2elm_t a)
{ // Modular correction, a = a in GF(p^2).
    fpcorrection(a[0]);
    fpcorrection(a[1]);
}


__inline static void mp_addfast(const digit_t* a, const digit_t* b, digit_t* c)
{ // Multiprecision addition, c = a+b.    
#if (OS_TARGET == OS_WIN) || defined(GENERIC_IMPLEMENTATION)

    mp_add(a, b, c, NWORDS_FIELD);
    
#elif (OS_TARGET == OS_LINUX)                 
    
    mp_add_asm(a, b, c);    

#endif
}


__inline static void mp_addfastx2(const digit_t* a, const digit_t* b, digit_t* c)
{ // Double-length multiprecision addition, c = a+b.    
#if (OS_TARGET == OS_WIN) || defined(GENERIC_IMPLEMENTATION)

    mp_add(a, b, c, 2*NWORDS_FIELD);
    
#elif (OS_TARGET == OS_LINUX)                 
    
    mp_addx2_asm(a, b, c);    

#endif
}


void fp2sqr_mont(const f2elm_t a, f2elm_t c)
{ // GF(p^2) squaring using Montgomery arithmetic, c = a^2 in GF(p^2).
  // Inputs: a = a0+a1*i, where a0, a1 are in [0, 2*p-1] 
  // Output: c = c0+c1*i, where c0, c1 are in [0, 2*p-1] 
    felm_t t1, t2, t3;
    
    mp_addfast(a[0], a[1], t1);                      // t1 = a0+a1 
    fpsub(a[0], a[1], t2);                           // t2 = a0-a1
    mp_addfast(a[0], a[0], t3);                      // t3 = 2a0
    fpmul_mont(t1, t2, c[0]);                        // c0 = (a0+a1)(a0-a1)
    fpmul_mont(t3, a[1], c[1]);                      // c1 = 2a0*a1
}


__inline unsigned int mp_sub(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords)
{ // Multiprecision subtraction, c = a-b, where lng(a) = lng(b) = nwords. Returns the borrow bit.
    unsigned int i, borrow = 0;

    for (i = 0; i < nwords; i++) {
        SUBC(borrow, a[i], b[i], borrow, c[i]);
    }

    return borrow;
}


__inline static digit_t mp_subfast(const digit_t* a, const digit_t* b, digit_t* c)
{ // Multiprecision subtraction, c = a-b, where lng(a) = lng(b) = 2*NWORDS_FIELD. 
  // If c < 0 then returns mask = 0xFF..F, else mask = 0x00..0   
#if (OS_TARGET == OS_WIN) || defined(GENERIC_IMPLEMENTATION)

	return (0 - (digit_t)mp_sub(a, b, c, 2*NWORDS_FIELD));

#elif (OS_TARGET == OS_LINUX)                 

	return mp_subx2_asm(a, b, c);

#endif
}


void fp2mul_mont(const f2elm_t a, const f2elm_t b, f2elm_t c)
{ // GF(p^2) multiplication using Montgomery arithmetic, c = a*b in GF(p^2).
  // Inputs: a = a0+a1*i and b = b0+b1*i, where a0, a1, b0, b1 are in [0, 2*p-1] 
  // Output: c = c0+c1*i, where c0, c1 are in [0, 2*p-1] 
    felm_t t1, t2;
    dfelm_t tt1, tt2, tt3; 
    digit_t mask;
    unsigned int i, borrow = 0;
    
	MUL512(a[0], b[0], tt1);
	MUL512(a[1], b[1], tt2);
    //mp_mul(a[0], b[0], tt1, NWORDS_FIELD);           // tt1 = a0*b0
    //mp_mul(a[1], b[1], tt2, NWORDS_FIELD);           // tt2 = a1*b1
    mp_addfast(a[0], a[1], t1);                      // t1 = a0+a1
    mp_addfast(b[0], b[1], t2);                      // t2 = b0+b1
    mask = mp_subfast(tt1, tt2, tt3);                // tt3 = a0*b0 - a1*b1. If tt3 < 0 then mask = 0xFF..F, else if tt3 >= 0 then mask = 0x00..0
    for (i = 0; i < NWORDS_FIELD; i++) {
        ADDC(borrow, tt3[NWORDS_FIELD+i], ((digit_t*)PRIME)[i] & mask, borrow, tt3[NWORDS_FIELD+i]);
    }
    //rdc_mont(tt3, c[0]);                             // c[0] = a0*b0 - a1*b1
	RED512(tt3,modulo_p503,c[0]);
    mp_addfastx2(tt1, tt2, tt1);                     // tt1 = a0*b0 + a1*b1
    MUL512(t1, t2, tt2);
	//mp_mul(t1, t2, tt2, NWORDS_FIELD);               // tt2 = (a0+a1)*(b0+b1)
	mp_subfast(tt2, tt1, tt2);                       // tt2 = (a0+a1)*(b0+b1) - a0*b0 - a1*b1 
    //rdc_mont(tt2, c[1]);                             // c[1] = (a0+a1)*(b0+b1) - a0*b0 - a1*b1 
	RED512(tt2,modulo_p503,c[1]);
}


void fpinv_chain_mont(felm_t a)
{ // Chain to compute a^(p-3)/4 using Montgomery arithmetic.
    unsigned int i, j;
    
#if (NBITS_FIELD == 503)
    felm_t t[15], tt;

    // Precomputed table
    fpsqr_mont(a, tt);
    fpmul_mont(a, tt, t[0]);
    for (i = 0; i <= 13; i++) fpmul_mont(t[i], tt, t[i+1]);

    fpcopy(a, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[3], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 12; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[14], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[14], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (j = 0; j < 49; j++) {
        for (i = 0; i < 5; i++) fpsqr_mont(tt, tt);
        fpmul_mont(t[14], tt, tt);
    }
    fpcopy(tt, a);  

#elif (NBITS_FIELD == 751)
    felm_t t[27], tt;
    
    // Precomputed table
    fpsqr_mont(a, tt);
    fpmul_mont(a, tt, t[0]);
    fpmul_mont(t[0], tt, t[1]);
    fpmul_mont(t[1], tt, t[2]);
    fpmul_mont(t[2], tt, t[3]); 
    fpmul_mont(t[3], tt, t[3]);
    for (i = 3; i <= 8; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[9], tt, t[9]);
    for (i = 9; i <= 20; i++) fpmul_mont(t[i], tt, t[i+1]);
    fpmul_mont(t[21], tt, t[21]); 
    for (i = 21; i <= 24; i++) fpmul_mont(t[i], tt, t[i+1]); 
    fpmul_mont(t[25], tt, t[25]);
    fpmul_mont(t[25], tt, t[26]);

    fpcopy(a, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[20], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[24], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[23], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 9; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[15], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[26], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[20], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[14], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i = 0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[18], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[1], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[22], tt, tt);
    for (i = 0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[6], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[24], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[18], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[17], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(a, tt, tt);
    for (i = 0; i < 10; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[16], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[7], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[0], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[19], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[22], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[25], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[10], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[22], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[18], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[4], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[14], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[23], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[21], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[23], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[12], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[9], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[3], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[13], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[17], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[26], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[5], tt, tt);
    for (i = 0; i < 8; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[8], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[2], tt, tt);
    for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[11], tt, tt);
    for (i = 0; i < 7; i++) fpsqr_mont(tt, tt);
    fpmul_mont(t[20], tt, tt);
    for (j = 0; j < 61; j++) {
        for (i = 0; i < 6; i++) fpsqr_mont(tt, tt);
        fpmul_mont(t[26], tt, tt);
    }
    fpcopy(tt, a);  
#endif
}


void fp2inv_mont(f2elm_t a)
{// GF(p^2) inversion using Montgomery arithmetic, a = (a0-i*a1)/(a0^2+a1^2).
    f2elm_t t1;

    fpsqr_mont(a[0], t1[0]);                         // t10 = a0^2
    fpsqr_mont(a[1], t1[1]);                         // t11 = a1^2
    fpadd(t1[0], t1[1], t1[0]);                      // t10 = a0^2+a1^2
    fpinv_mont(t1[0]);                               // t10 = (a0^2+a1^2)^-1
    fpneg(a[1]);                                     // a = a0-i*a1
    fpmul_mont(a[0], t1[0], a[0]);
    fpmul_mont(a[1], t1[0], a[1]);                   // a = (a0-i*a1)*(a0^2+a1^2)^-1
}


void to_fp2mont(const f2elm_t a, f2elm_t mc)
{ // Conversion of a GF(p^2) element to Montgomery representation,
  // mc_i = a_i*R^2*R^(-1) = a_i*R in GF(p^2). 

    to_mont(a[0], mc[0]);
    to_mont(a[1], mc[1]);
}


void from_fp2mont(const f2elm_t ma, f2elm_t c)
{ // Conversion of a GF(p^2) element from Montgomery representation to standard representation,
  // c_i = ma_i*R^(-1) = a_i in GF(p^2).

    from_mont(ma[0], c[0]);
    from_mont(ma[1], c[1]);
}


__inline unsigned int mp_add(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords)
{ // Multiprecision addition, c = a+b, where lng(a) = lng(b) = nwords. Returns the carry bit.
    unsigned int i, carry = 0;
        
    for (i = 0; i < nwords; i++) {                      
        ADDC(carry, a[i], b[i], carry, c[i]);
    }

    return carry;
}


void mp_shiftleft(digit_t* x, unsigned int shift, const unsigned int nwords)
{
    unsigned int i, j = 0;

    while (shift > RADIX) {
        j += 1;
        shift -= RADIX;
    }

    for (i = 0; i < nwords-j; i++) 
        x[nwords-1-i] = x[nwords-1-i-j];
    for (i = nwords-j; i < nwords; i++) 
        x[nwords-1-i] = 0;
    if (shift != 0) {
        for (j = nwords-1; j > 0; j--) 
            SHIFTL(x[j], x[j-1], shift, x[j], RADIX);
        x[0] <<= shift;
    }
}


void mp_shiftr1(digit_t* x, const unsigned int nwords)
{ // Multiprecision right shift by one.
    unsigned int i;

    for (i = 0; i < nwords-1; i++) {
        SHIFTR(x[i+1], x[i], 1, x[i], RADIX);
    }
    x[nwords-1] >>= 1;
}


void mp_shiftl1(digit_t* x, const unsigned int nwords)
{ // Multiprecision left shift by one.
    int i;

    for (i = nwords-1; i > 0; i--) {
        SHIFTL(x[i], x[i-1], 1, x[i], RADIX);
    }
    x[0] <<= 1;
}
