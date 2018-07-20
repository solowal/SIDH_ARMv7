/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: core functions over GF(p) and GF(p^2)
*********************************************************************************************/
#include <arm_neon.h>
#include <stdio.h>

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
void __attribute__ ((noinline, naked)) MUL512(const digit_t* a, const digit_t* b, digit_t* c){
	asm(
			"push  {r0-r11,lr}			\n\t"
			"vpush {q4-q7}				\n\t"
			
			//1 2 3 4 5 6 7 8
			//9 10 11 12

			"ldmia r0!, {r1-r12} 		\n\t"
			"subs r1, r1, r9			\n\t"
			"sbcs r2, r2, r10			\n\t"
			"sbcs r3, r3, r11			\n\t"
			"sbcs r4, r4, r12			\n\t"

			"ldmia r0!, {r9-r12} 		\n\t"
			"sbcs r5, r5, r9			\n\t"
			"sbcs r6, r6, r10			\n\t"
			"sbcs r7, r7, r11			\n\t"
			"sbcs r8, r8, r12			\n\t"
			"sbcs r14, r14, r14			\n\t"

			"eor r1, r1, r14			\n\t"
			"eor r2, r2, r14			\n\t"
			"eor r3, r3, r14			\n\t"
			"eor r4, r4, r14			\n\t"

			"eor r5, r5, r14			\n\t"
			"eor r6, r6, r14			\n\t"
			"eor r7, r7, r14			\n\t"
			"eor r8, r8, r14			\n\t"

			"adds r1, r1, r14, LSR #31	\n\t"
			"adcs r2, r2, #0			\n\t"
			"adcs r3, r3, #0			\n\t"
			"adcs r4, r4, #0			\n\t"

			"adcs r5, r5, #0			\n\t"
			"adcs r6, r6, #0			\n\t"
			"adcs r7, r7, #0			\n\t"
			"adcs r8, r8, #0			\n\t"

			"vmov d0, r1, r2 			\n\t"
			"vmov d1, r3, r4 			\n\t"
			"vmov d2, r5, r6 			\n\t"
			"vmov d3, r7, r8 			\n\t"
			//carry bit: r14
			
			//second round
			"ldr r0, [sp, #4 * 17]		\n\t"//loading r1

			"ldmia r0!, {r1-r12} 		\n\t"
			"subs r1, r1, r9			\n\t"
			"sbcs r2, r2, r10			\n\t"
			"sbcs r3, r3, r11			\n\t"
			"sbcs r4, r4, r12			\n\t"

			"ldmia r0!, {r9-r12} 		\n\t"
			"sbcs r5, r5, r9			\n\t"
			"sbcs r6, r6, r10			\n\t"
			"sbcs r7, r7, r11			\n\t"
			"sbcs r8, r8, r12			\n\t"
			"sbcs r12, r12, r12			\n\t"

			"eor r1, r1, r12			\n\t"
			"eor r2, r2, r12			\n\t"
			"eor r3, r3, r12			\n\t"
			"eor r4, r4, r12			\n\t"

			"eor r5, r5, r12			\n\t"
			"eor r6, r6, r12			\n\t"
			"eor r7, r7, r12			\n\t"
			"eor r8, r8, r12			\n\t"			

			"adds r1, r1, r12, LSR #31	\n\t"
			"adcs r2, r2, #0			\n\t"
			"adcs r3, r3, #0			\n\t"
			"adcs r4, r4, #0			\n\t"

			"adcs r5, r5, #0			\n\t"
			"adcs r6, r6, #0			\n\t"
			"adcs r7, r7, #0			\n\t"
			"adcs r8, r8, #0			\n\t"
			//carry bit: r12
			"eor r12, r12, r14			\n\t"//carry integration
			"mvn r12, r12				\n\t"

			"vmov d4, r1, r2 			\n\t"
			"vmov d5, r3, r4 			\n\t"
			"vmov d6, r5, r6 			\n\t"
			"vmov d7, r7, r8 			\n\t"
			"vdup.32 q4, r12	 		\n\t"//duplication

			////////////////////////////////
			"ldr r0, [sp, #4 * 18]			\n\t"//loading result
			"ldr r1, [sp, #4 * 16]			\n\t"//loading op1
			"ldr r2, [sp, #4 * 17]			\n\t"//loading op2
			////////////////////////////////


			"ldr r3, [r1, #4 * 6]			\n\t"//A6
			"ldr r4, [r1, #4 * 7]			\n\t"//A7
			"vtrn.32 q0, q1 			\n\t "			
			"ldr r6, [r2, #4 * 0]			\n\t"//B0
			"ldr r7, [r2, #4 * 1]			\n\t"//B1
			"vmull.u32 q13, d0, d4[0] 	\n\t"	//A4A0 X B0B0
			"ldr r8, [r2, #4 * 2]			\n\t"//B2
			"umull r9, r10, r3, r6			\n\t"//60
			"vmull.u32 q12, d2, d4[0] 	\n\t"	//A5A1 X B0B0
			"str r9, [r0, #4 * 6]			\n\t"//C6
			"mov r9, #0						\n\t"
			"vmull.u32 q11, d1, d4[0] 	\n\t"	//A6A2 X B0B0
			"umaal r9, r10, r4, r6			\n\t"//70
			"mov r11, #0					\n\t"
			"vmull.u32 q10, d3, d4[0] 	\n\t"	//A7A3 X B0B0
			"umaal r9, r11, r3, r7			\n\t"//61
			"str r9, [r0, #4 * 7]			\n\t"//C7
			"veor q6,q6,q6 				\n\t"
			"umaal r10, r11, r4, r7			\n\t"//71
			"str r10, [r0, #4 * 8]			\n\t"//C8
			"veor q7,q7,q7 				\n\t"
			"str r11, [r0, #4 * 9]			\n\t"//C9
			"ldr r3, [r1, #4 * 0]			\n\t"//A0
			"veor q8,q8,q8 				\n\t"
			"ldr r4, [r1, #4 * 1]			\n\t"//A1
			"ldr r5, [r1, #4 * 2]			\n\t"//A2
			"veor q9,q9,q9 				\n\t"
			"umull r9, r10, r3, r6			\n\t"//00
			"str r9, [r0, #4 * 0]			\n\t"//C0
			"vtrn.32 q13, q6 			\n\t"
			"mov r9, #0						\n\t"
			"umaal r9, r10, r4, r6			\n\t"//10
			"vtrn.32 q12, q7 			\n\t"
			"mov r11, #0					\n\t"
			"umaal r9, r11, r3, r7			\n\t"//01
			"vtrn.32 q11, q8 			\n\t"
			"str r9, [r0, #4 * 1]			\n\t"//C1
			"mov r12, #0					\n\t"
			"vtrn.32 q10, q9 			\n\t"
			"umaal r12, r10, r5, r6			\n\t"//20
			"mov r9, #0						\n\t"
			"vadd.i64  q12, q12, q6 	\n\t"	//q
			"umaal r9, r11, r4, r7			\n\t"//11
			"umaal r12, r9, r3, r8			\n\t"//02
			"vadd.i64  q11, q11, q7 	\n\t"	//q
			"str r12, [r0, #4 * 2]			\n\t"//C2
			"mov r12, #0					\n\t"
			"vadd.i64  q10, q10, q8 	\n\t"	//q
			"ldr r6, [r2, #4 * 3]			\n\t"//B3
			"umaal r9, r10, r5, r7			\n\t"//21
			"vqadd.u64 d18, d18, d27 	\n\t"	//d
			"umaal r9, r11, r4, r8			\n\t"//12
			"umaal r9, r12, r3, r6			\n\t"//03
			"vext.8 d28, d28, d26, #4 	\n\t"
			"str r9, [r0, #4 * 3]			\n\t"//C3
			"mov r9, #0						\n\t"
			"vmlal.u32 q12, d0, d4[1] 	\n\t"	//A4A0 X B0B0
			"ldr r7, [r2, #4 * 4]			\n\t"//B4
			"umaal r9, r10, r5, r8			\n\t"//22
			"vmlal.u32 q11, d2, d4[1] 	\n\t"	//A5A1 X B0B0
			"umaal r9, r11, r4, r6			\n\t"//13
			"umaal r9, r12, r3, r7			\n\t"//04
			"vmlal.u32 q10, d1, d4[1] 	\n\t"	//A6A2 X B0B0
			"str r9, [r0, #4 * 4]			\n\t"//C4
			"mov r9, #0						\n\t"			
			"vmlal.u32 q9 , d3, d4[1] 	\n\t"	//A7A3 X B0B0
			"ldr r8, [r2, #4 * 5]			\n\t"//B5
			"umaal r9, r10, r5, r6			\n\t"//23
			"veor q13,q13,q13			\n\t"
			"umaal r9, r11, r4, r7			\n\t"//14
			"umaal r9, r12, r3, r8			\n\t"//05
			"veor q6 ,q6 ,q6 				\n\t"
			"str r9, [r0, #4 * 5]			\n\t"//C5
			"ldr r9, [r0, #4 * 6]			\n\t"//C6
			"veor q7 ,q7 ,q7 				\n\t"
			"ldr r6, [r2, #4 * 6]			\n\t"//B6
			"umaal r9, r10, r5, r7			\n\t"//24
			"veor q8 ,q8 ,q8 				\n\t"
			"umaal r9, r11, r4, r8			\n\t"//15
			"umaal r9, r12, r3, r6			\n\t"//06
			"vtrn.32 q12, q13 			\n\t"
			"str r9, [r0, #4 * 6]			\n\t"//C6
			"ldr r9, [r0, #4 * 7]			\n\t"//C7
			"vtrn.32 q11, q6 			\n\t"
			"ldr r7, [r2, #4 * 7]			\n\t"//B7
			"umaal r9, r10, r5, r8			\n\t"//25
			"vtrn.32 q10, q7 			\n\t"
			"umaal r9, r11, r4, r6			\n\t"//16
			"umaal r9, r12, r3, r7			\n\t"//07
			"vtrn.32 q9 , q8 			\n\t"
			"str r9, [r0, #4 * 7]			\n\t"//C7
			"ldr r9, [r0, #4 * 8]			\n\t"//C8		
			"vadd.i64  q11, q11, q13 	\n\t"	//q
			"ldr r3, [r1, #4 * 3]			\n\t"//A3
			"umaal r9, r10, r3, r8			\n\t"//35
			"vadd.i64  q10, q10, q6 	\n\t"	//q
			"umaal r9, r11, r5, r6			\n\t"//26
			"umaal r9, r12, r4, r7			\n\t"//17
			"vadd.i64  q9 , q9 , q7 	\n\t"	//q
			"str r9, [r0, #4 * 8]			\n\t"//C8
			"ldr r9, [r0, #4 * 9]			\n\t"//C9		
			"vqadd.u64 d16, d16, d25 	\n\t"	//d
			"ldr r4, [r1, #4 * 4]			\n\t"//A4
			"umaal r9, r10, r4, r8			\n\t"//45
			"vext.8 d28, d28, d24, #4 	\n\t"
			"umaal r9, r11, r3, r6			\n\t"//36
			"umaal r9, r12, r5, r7			\n\t"//27
			"vmlal.u32 q11, d0, d5[0] 	\n\t"	//A4A0 X B0B0
			"str r9, [r0, #4 * 9]			\n\t"//C9
			"umaal r10, r11, r4, r6			\n\t"//46
			"vmlal.u32 q10, d2, d5[0] 	\n\t"	//A5A1 X B0B0
			"umaal r10, r12, r3, r7			\n\t"//37
			"str r10, [r0, #4 * 10]			\n\t"//C10
			"vmlal.u32 q9 , d1, d5[0] 	\n\t"	//A6A2 X B0B0
			"umaal r11, r12, r4, r7			\n\t"//47
			"str r11, [r0, #4 * 11]			\n\t"//C11
			"vmlal.u32 q8 , d3, d5[0] 	\n\t"	//A7A3 X B0B0
			"str r12, [r0, #4 * 12]			\n\t"//C12
			"ldr r5, [r1, #4 * 5]			\n\t"//A5
			"veor q12,q12,q12			\n\t"
			"ldr r6, [r2, #4 * 0]			\n\t"//B0
			"ldr r7, [r2, #4 * 1]			\n\t"//B1
			"veor q13,q13,q13			\n\t"
			"ldr r8, [r2, #4 * 2]			\n\t"//B2
			"ldr r9, [r0, #4 * 3]			\n\t"//C3
			"veor q6 ,q6 ,q6 			\n\t"
			"mov r10, #0					\n\t"
			"umaal r9, r10, r3, r6			\n\t"//30
			"veor q7 ,q7 ,q7 			\n\t"
			"str r9, [r0, #4 * 3]			\n\t"//C3
			"ldr r9, [r0, #4 * 4]			\n\t"//C4
			"vtrn.32 q11, q12 			\n\t"
			"mov r11, #0					\n\t"
			"umaal r9, r10, r4, r6			\n\t"//40
			"vtrn.32 q10, q13			\n\t"
			"umaal r9, r11, r3, r7			\n\t"//31
			"str r9, [r0, #4 * 4]			\n\t"//C4
			"vtrn.32 q9 , q6			\n\t"
			"ldr r9, [r0, #4 * 5]			\n\t"//C5
			"mov r12, #0					\n\t"
			"vtrn.32 q8 , q7 			\n\t"
			"umaal r9, r10, r5, r6			\n\t"//50
			"umaal r9, r11, r4, r7			\n\t"//41
			"vadd.i64  q10, q10, q12 	\n\t"	//q
			"umaal r9, r12, r3, r8			\n\t"//32
			"str r9, [r0, #4 * 5]			\n\t"//C5
			"vadd.i64  q9 , q9 , q13	\n\t"	//q
			"ldr r9, [r0, #4 * 6]			\n\t"//C6
			"ldr r6, [r2, #4 * 3]			\n\t"//B3
			"vadd.i64  q8 , q8 , q6 	\n\t"	//q
			"umaal r9, r10, r5, r7			\n\t"//51
			"umaal r9, r11, r4, r8			\n\t"//42
			"vqadd.u64 d14, d14, d23 	\n\t"	//d
			"umaal r9, r12, r3, r6			\n\t"//33
			"str r9, [r0, #4 * 6]			\n\t"//C6
			"vext.8 d29, d29, d22, #4 	\n\t"
			"ldr r9, [r0, #4 * 7]			\n\t"//C7
			"ldr r7, [r2, #4 * 4]			\n\t"//B4
			"vmlal.u32 q10, d0, d5[1] 	\n\t"	//A4A0 X B0B0
			"umaal r9, r10, r5, r8			\n\t"//52
			"umaal r9, r11, r4, r6			\n\t"//43
			"vmlal.u32 q9 , d2, d5[1] 	\n\t"	//A5A1 X B0B0
			"umaal r9, r12, r3, r7			\n\t"//34
			"str r9, [r0, #4 * 7]			\n\t"//C7
			"vmlal.u32 q8 , d1, d5[1] 	\n\t"	//A6A2 X B0B0
			"ldr r9, [r0, #4 * 8]			\n\t"//C8
			"ldr r3, [r1, #4 * 6]			\n\t"//A6
			"vmlal.u32 q7 , d3, d5[1] 	\n\t"	//A7A3 X B0B0
			"umaal r9, r10, r3, r8			\n\t"//62
			"umaal r9, r11, r5, r6			\n\t"//53
			"veor q11,q11,q11			\n\t"
			"umaal r9, r12, r4, r7			\n\t"//44
			"str r9, [r0, #4 * 8]			\n\t"//C8
			"veor q12,q12,q12			\n\t"
			"ldr r9, [r0, #4 * 9]			\n\t"//C9
			"ldr r4, [r1, #4 * 7]			\n\t"//A7
			"veor q13,q13,q13			\n\t"
			"umaal r9, r10, r4, r8			\n\t"//72
			"umaal r9, r11, r3, r6			\n\t"//63
			"veor q6 ,q6 ,q6 			\n\t"
			"umaal r9, r12, r5, r7			\n\t"//54
			"str r9, [r0, #4 * 9]			\n\t"//C9
			"vtrn.32 q10, q11			\n\t"
			"ldr r9, [r0, #4 * 10]			\n\t"//C10
			"ldr r8, [r2, #4 * 5]			\n\t"//B5
			"vtrn.32 q9 , q12			\n\t"
			"umaal r9, r10, r4, r6			\n\t"//73
			"umaal r9, r11, r3, r7			\n\t"//64
			"vtrn.32 q8 , q13			\n\t"
			"umaal r9, r12, r5, r8			\n\t"//55
			"str r9, [r0, #4 * 10]			\n\t"//C10
			"vtrn.32 q7 , q6 			\n\t"
			"ldr r9, [r0, #4 * 11]			\n\t"//C11
			"ldr r6, [r2, #4 * 6]			\n\t"//B6
			"vadd.i64  q9 , q9 , q11	\n\t"	//q
			"umaal r9, r10, r4, r7			\n\t"//74
			"umaal r9, r11, r3, r8			\n\t"//65
			"vadd.i64  q8 , q8 , q12	\n\t"	//q
			"umaal r9, r12, r5, r6			\n\t"//56
			"str r9, [r0, #4 * 11]			\n\t"//C11
			"vadd.i64  q7 , q7 , q13 	\n\t"	//q
			"ldr r9, [r0, #4 * 12]			\n\t"//C12
			"ldr r7, [r2, #4 * 7]			\n\t"//B7
			"vqadd.u64 d12, d12, d21 	\n\t"	//d
			"umaal r9, r10, r4, r8			\n\t"//75
			"umaal r9, r11, r3, r6			\n\t"//66
			"vext.8 d29, d29, d20, #4 	\n\t"
			"umaal r9, r12, r5, r7			\n\t"//57
			"str r9, [r0, #4 * 12]			\n\t"//C12
			"vmlal.u32 q9 , d0, d6[0] 	\n\t"	//A4A0 X B0B0
			"umaal r10, r11, r4, r6			\n\t"//76
			"umaal r10, r12, r3, r7			\n\t"//67
			"vmlal.u32 q8 , d2, d6[0] 	\n\t"	//A5A1 X B0B0
			"str r10, [r0, #4 * 13]			\n\t"//C13
			"umaal r11, r12, r4, r7			\n\t"//77
			"vmlal.u32 q7 , d1, d6[0] 	\n\t"	//A6A2 X B0B0
			"str r11, [r0, #4 * 14]			\n\t"//C14
			"str r12, [r0, #4 * 15]			\n\t"//C15
			"vmlal.u32 q6 , d3, d6[0] 	\n\t"	//A7A3 X B0B0
			"ldr r3, [r1, #4 * 14]			\n\t"//A6
			"ldr r4, [r1, #4 * 15]			\n\t"//A7
			"veor q10,q10,q10			\n\t"
			"ldr r6, [r2, #4 * 8]			\n\t"//B0
			"ldr r7, [r2, #4 * 9]			\n\t"//B1
			"veor q11,q11,q11			\n\t"
			"ldr r8, [r2, #4 * 10]			\n\t"//B2
			"ldr r9, [r0, #4 * 14]			\n\t"//C6
			"veor q12,q12,q12			\n\t"
			"mov r10, #0					\n\t"
			"umaal r9, r10, r3, r6			\n\t"//60
			"veor q13,q13,q13			\n\t"
			"str r9, [r0, #4 * 14]			\n\t"//C6
			"ldr r9, [r0, #4 * 15]			\n\t"//C7
			"vtrn.32 q9 , q10			\n\t"
			"umaal r9, r10, r4, r6			\n\t"//70
			"mov r11, #0					\n\t"
			"vtrn.32 q8 , q11			\n\t"
			"umaal r9, r11, r3, r7			\n\t"//61
			"str r9, [r0, #4 * 15]			\n\t"//C7
			"vtrn.32 q7 , q12			\n\t"
			"umaal r10, r11, r4, r7			\n\t"//71
			"str r10, [r0, #4 * 24]			\n\t"//C8
			"vtrn.32 q6 , q13			\n\t"
			"str r11, [r0, #4 * 25]			\n\t"//C9
			"ldr r3, [r1, #4 * 8]			\n\t"//A0
			"vadd.i64  q8 , q8 , q10	\n\t"	//q
			"ldr r4, [r1, #4 * 9]			\n\t"//A1
			"ldr r5, [r1, #4 * 10]			\n\t"//A2
			"vadd.i64  q7 , q7 , q11 	\n\t"	//q
			"ldr r9, [r0, #4 * 8]			\n\t"//C0
			"mov r10, #0					\n\t"
			"vadd.i64  q6 , q6 , q12	\n\t"	//q
			"umaal r9, r10, r3, r6			\n\t"//00
			"str r9, [r0, #4 * 16]			\n\t"//C0
			"vqadd.u64 d26, d26, d19 	\n\t"	//d
			"ldr r9, [r0, #4 * 9]			\n\t"//C1
			"umaal r9, r10, r4, r6			\n\t"//10
			"vext.8 d30, d30, d18, #4 	\n\t"
			"mov r11, #0					\n\t"
			"umaal r9, r11, r3, r7			\n\t"//01
			"vmlal.u32 q8 , d0, d6[1] 	\n\t"	//A4A0 X B0B0
			"str r9, [r0, #4 * 17]			\n\t"//C1
			"ldr r12, [r0, #4 * 10]			\n\t"//C2
			"vmlal.u32 q7 , d2, d6[1] 	\n\t"	//A5A1 X B0B0
			"umaal r12, r10, r5, r6			\n\t"//20
			"mov r9, #0						\n\t"
			"vmlal.u32 q6 , d1, d6[1] 	\n\t"	//A6A2 X B0B0
			"umaal r9, r11, r4, r7			\n\t"//11
			"umaal r12, r9, r3, r8			\n\t"//02
			"vmlal.u32 q13, d3, d6[1] 	\n\t"	//A7A3 X B0B0
			"str r12, [r0, #4 * 18]			\n\t"//C2
			"ldr r12, [r0, #4 * 11]			\n\t"//C3
			"veor q9 ,q9 ,q9			\n\t"
			"ldr r6, [r2, #4 * 11]			\n\t"//B3
			"umaal r9, r10, r5, r7			\n\t"//21
			"veor q10,q10,q10			\n\t"
			"umaal r9, r11, r4, r8			\n\t"//12
			"umaal r9, r12, r3, r6			\n\t"//03
			"veor q11,q11,q11			\n\t"
			"str r9, [r0, #4 * 11]			\n\t"//C3
			"ldr r9, [r0, #4 * 12]			\n\t"//C4
			"veor q12,q12,q12			\n\t"
			"ldr r7, [r2, #4 * 12]			\n\t"//B4
			"umaal r9, r10, r5, r8			\n\t"//22
			"vtrn.32 q8 , q9			\n\t"
			"umaal r9, r11, r4, r6			\n\t"//13
			"umaal r9, r12, r3, r7			\n\t"//04
			"vtrn.32 q7 , q10			\n\t"
			"str r9, [r0, #4 * 12]			\n\t"//C4
			"ldr r9, [r0, #4 * 13]			\n\t"//C5		
			"vtrn.32 q6 , q11			\n\t"
			"ldr r8, [r2, #4 * 13]			\n\t"//B5
			"umaal r9, r10, r5, r6			\n\t"//23
			"vtrn.32 q13, q12			\n\t"
			"umaal r9, r11, r4, r7			\n\t"//14
			"umaal r9, r12, r3, r8			\n\t"//05
			"vadd.i64  q7 , q7 , q9 	\n\t"	//q
			"str r9, [r0, #4 * 13]			\n\t"//C5
			"ldr r9, [r0, #4 * 14]			\n\t"//C6
			"vadd.i64  q6 , q6 , q10	\n\t"	//q
			"ldr r6, [r2, #4 * 14]			\n\t"//B6
			"umaal r9, r10, r5, r7			\n\t"//24
			"vadd.i64  q13, q13, q11	\n\t"	//q
			"umaal r9, r11, r4, r8			\n\t"//15
			"umaal r9, r12, r3, r6			\n\t"//06
			"vqadd.u64 d24, d24, d17 	\n\t"	//d
			"str r9, [r0, #4 * 14]			\n\t"//C6
			"ldr r9, [r0, #4 * 15]			\n\t"//C7
			"vext.8 d30, d30, d16, #4 	\n\t"
			"ldr r7, [r2, #4 * 15]			\n\t"//B7
			"umaal r9, r10, r5, r8			\n\t"//25
			"vmlal.u32 q7 , d0, d7[0] 	\n\t"	//A4A0 X B0B0
			"umaal r9, r11, r4, r6			\n\t"//16
			"umaal r9, r12, r3, r7			\n\t"//07
			"vmlal.u32 q6 , d2, d7[0] 	\n\t"	//A5A1 X B0B0
			"str r9, [r0, #4 * 15]			\n\t"//C7
			"ldr r9, [r0, #4 * 24]			\n\t"//C8		
			"vmlal.u32 q13, d1, d7[0] 	\n\t"	//A6A2 X B0B0
			"ldr r3, [r1, #4 * 11]			\n\t"//A3
			"umaal r9, r10, r3, r8			\n\t"//35
			"vmlal.u32 q12, d3, d7[0] 	\n\t"	//A7A3 X B0B0
			"umaal r9, r11, r5, r6			\n\t"//26
			"umaal r9, r12, r4, r7			\n\t"//17
			"veor q8 ,q8 ,q8			\n\t"
			"str r9, [r0, #4 * 24]			\n\t"//C8
			"ldr r9, [r0, #4 * 25]			\n\t"//C9		
			"veor q9 ,q9 ,q9			\n\t"
			"ldr r4, [r1, #4 * 12]			\n\t"//A4
			"umaal r9, r10, r4, r8			\n\t"//45
			"veor q10,q10,q10			\n\t"
			"umaal r9, r11, r3, r6			\n\t"//36
			"umaal r9, r12, r5, r7			\n\t"//27
			"veor q11,q11,q11			\n\t"
			"str r9, [r0, #4 * 25]			\n\t"//C9
			"umaal r10, r11, r4, r6			\n\t"//46
			"vtrn.32 q7 , q8			\n\t"
			"umaal r10, r12, r3, r7			\n\t"//37
			"str r10, [r0, #4 * 26]			\n\t"//C10
			"vtrn.32 q6 , q9			\n\t"
			"umaal r11, r12, r4, r7			\n\t"//47
			"str r11, [r0, #4 * 27]			\n\t"//C11
			"vtrn.32 q13, q10			\n\t"
			"str r12, [r0, #4 * 28]			\n\t"//C12
			"ldr r5, [r1, #4 * 13]			\n\t"//A5
			"vtrn.32 q12, q11			\n\t"
			"ldr r6, [r2, #4 * 8]			\n\t"//B0
			"ldr r7, [r2, #4 * 9]			\n\t"//B1
			"vadd.i64  q6 , q6 , q8		\n\t"	//q
			"ldr r8, [r2, #4 * 10]			\n\t"//B2
			"ldr r9, [r0, #4 * 11]			\n\t"//C3
			"vadd.i64  q13, q13, q9		\n\t"	//q
			"mov r10, #0					\n\t"
			"umaal r9, r10, r3, r6			\n\t"//30
			"vadd.i64  q12, q12, q10	\n\t"	//q
			"str r9, [r0, #4 * 19]			\n\t"//C3
			"ldr r9, [r0, #4 * 12]			\n\t"//C4
			"vqadd.u64 d22, d22, d15 	\n\t"	//d
			"mov r11, #0					\n\t"
			"umaal r9, r10, r4, r6			\n\t"//40
			"vext.8 d31, d31, d14, #4 	\n\t"
			"umaal r9, r11, r3, r7			\n\t"//31
			"str r9, [r0, #4 * 20]			\n\t"//C4
			"vmlal.u32 q6 , d0, d7[1] 	\n\t"	//A4A0 X B0B0
			"ldr r9, [r0, #4 * 13]			\n\t"//C5
			"mov r12, #0					\n\t"
			"vmlal.u32 q13, d2, d7[1] 	\n\t"	//A5A1 X B0B0
			"umaal r9, r10, r5, r6			\n\t"//50
			"umaal r9, r11, r4, r7			\n\t"//41
			"vmlal.u32 q12, d1, d7[1] 	\n\t"	//A6A2 X B0B0
			"umaal r9, r12, r3, r8			\n\t"//32
			"str r9, [r0, #4 * 21]			\n\t"//C5
			"vmlal.u32 q11, d3, d7[1] 	\n\t"	//A7A3 X B0B0
			"ldr r9, [r0, #4 * 14]			\n\t"//C6
			"ldr r6, [r2, #4 * 11]			\n\t"//B3
			"veor q7, q7, q7			\n\t"
			"umaal r9, r10, r5, r7			\n\t"//51
			"umaal r9, r11, r4, r8			\n\t"//42
			"veor q8 ,q8 ,q8			\n\t"
			"umaal r9, r12, r3, r6			\n\t"//33
			"str r9, [r0, #4 * 22]			\n\t"//C6
			"veor q9 ,q9 ,q9			\n\t"
			"ldr r9, [r0, #4 * 15]			\n\t"//C7
			"ldr r7, [r2, #4 * 12]			\n\t"//B4
			"veor q10,q10,q10			\n\t"
			"umaal r9, r10, r5, r8			\n\t"//52
			"umaal r9, r11, r4, r6			\n\t"//43
			"vtrn.32 q6 , q7			\n\t"
			"umaal r9, r12, r3, r7			\n\t"//34
			"str r9, [r0, #4 * 23]			\n\t"//C7
			"vtrn.32 q13, q8			\n\t"
			"ldr r9, [r0, #4 * 24]			\n\t"//C8
			"ldr r3, [r1, #4 * 14]			\n\t"//A6
			"vtrn.32 q12, q9			\n\t"
			"umaal r9, r10, r3, r8			\n\t"//62
			"umaal r9, r11, r5, r6			\n\t"//53
			"vtrn.32 q11, q10			\n\t"
			"umaal r9, r12, r4, r7			\n\t"//44
			"str r9, [r0, #4 * 24]			\n\t"//C8
			"vadd.i64  q13, q13, q7		\n\t"	//q
			"ldr r9, [r0, #4 * 25]			\n\t"//C9
			"ldr r4, [r1, #4 * 15]			\n\t"//A7
			"vadd.i64  q12, q12, q8		\n\t"	//q
			"umaal r9, r10, r4, r8			\n\t"//72
			"umaal r9, r11, r3, r6			\n\t"//63
			"vadd.i64  q11, q11, q9		\n\t"	//q
			"umaal r9, r12, r5, r7			\n\t"//54
			"str r9, [r0, #4 * 25]			\n\t"//C9
			"vqadd.u64 d20, d20, d13 	\n\t"	//d
			"ldr r9, [r0, #4 * 26]			\n\t"//C10
			"ldr r8, [r2, #4 * 13]			\n\t"//B5
			"vext.8 d31, d31, d12, #4 	\n\t"
			"umaal r9, r10, r4, r6			\n\t"//73
			"umaal r9, r11, r3, r7			\n\t"//64
			"umaal r9, r12, r5, r8			\n\t"//55
			"str r9, [r0, #4 * 26]			\n\t"//C10
			"ldr r9, [r0, #4 * 27]			\n\t"//C11
			"ldr r6, [r2, #4 * 14]			\n\t"//B6
			"umaal r9, r10, r4, r7			\n\t"//74
			"umaal r9, r11, r3, r8			\n\t"//65
			"umaal r9, r12, r5, r6			\n\t"//56
			"str r9, [r0, #4 * 27]			\n\t"//C11
			"ldr r9, [r0, #4 * 28]			\n\t"//C12
			"ldr r7, [r2, #4 * 15]			\n\t"//B7
			"umaal r9, r10, r4, r8			\n\t"//75
			"umaal r9, r11, r3, r6			\n\t"//66
			"umaal r9, r12, r5, r7			\n\t"//57
			"str r9, [r0, #4 * 28]			\n\t"//C12
			"umaal r10, r11, r4, r6			\n\t"//76
			"umaal r10, r12, r3, r7			\n\t"//67
			"str r10, [r0, #4 * 29]			\n\t"//C13
			"umaal r11, r12, r4, r7			\n\t"//77
			"str r11, [r0, #4 * 30]			\n\t"//C14
			"str r12, [r0, #4 * 31]			\n\t"//C15

			//integration
			//q14 q15
			
			"veor q14, q14, q4				\n\t"//could be move to up
			"veor q15, q15, q4				\n\t"

			"ldmia r0!, {r1-r8} 			\n\t"
			
			"vmov r9, r10, d28 				\n\t"
			"vmov r11, r12, d29 			\n\t"

			"vmov r14, d8[0]	 			\n\t"//carry
			
			"adds r14, r14, r14				\n\t"//carry

			"adcs r1 , r9 , r1				\n\t"
			"adcs r2, r10, r2				\n\t"
			"adcs r3, r11, r3				\n\t"
			"adcs r4, r12, r4				\n\t"

			"vmov r9, r10, d30 				\n\t"
			"vmov r11, r12, d31 			\n\t"

			"adcs r5, r9, r5				\n\t"
			"adcs r6, r10, r6				\n\t"
			"adcs r7, r11, r7				\n\t"
			"adcs r8, r12, r8				\n\t"

			"mov r14, #0					\n\t"
			"adcs r14, r14, #0				\n\t"//carry


			"ldr r9, [r0, #4 * 8]			\n\t"//C16
			"ldr r10, [r0, #4 * 9]			\n\t"//C17
			"ldr r11, [r0, #4 * 10]			\n\t"//C18
			"ldr r12, [r0, #4 * 11]			\n\t"//C19

			"adds r1 , r9 , r1				\n\t"
			"adcs r2, r10, r2				\n\t"
			"adcs r3, r11, r3				\n\t"
			"adcs r4, r12, r4				\n\t"

			"ldr r9, [r0, #4 * 12]			\n\t"//C20
			"ldr r10, [r0, #4 * 13]			\n\t"//C21
			"ldr r11, [r0, #4 * 14]			\n\t"//C22
			"ldr r12, [r0, #4 * 15]			\n\t"//C23

			"adcs r5, r9, r5				\n\t"
			"adcs r6, r10, r6				\n\t"
			"adcs r7, r11, r7				\n\t"
			"adcs r8, r12, r8				\n\t"

			//"mov r14, #0					\n\t"
			"adcs r14, r14, #0				\n\t"//carry

			"stmia r0!, {r1-r8} 			\n\t"
			//second part finished

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


			"vmov r12, d8[0]	 			\n\t"//carry

			"eor r1, r1, r12				\n\t"
			"eor r2, r2, r12				\n\t"
			"eor r3, r3, r12				\n\t"
			"eor r4, r4, r12				\n\t"

			"eor r5, r5, r12				\n\t"
			"eor r6, r6, r12				\n\t"
			"eor r7, r7, r12				\n\t"
			"eor r8, r8, r12				\n\t"


			//"adds r12, r12, r12				\n\t"//carry

			"adds r1, r1, r14				\n\t"
			"adcs r2, r2, #0				\n\t"
			"adcs r3, r3, #0				\n\t"
			"adcs r4, r4, #0				\n\t"
			"adcs r5, r5, #0				\n\t"
			"adcs r6, r6, #0				\n\t"
			"adcs r7, r7, #0				\n\t"
			"adcs r8, r8, #0				\n\t"

			"mov r14, #0					\n\t"
			"adcs r14, r14, #0				\n\t"//carry
			/////////////////



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

			"ldr r9, [r0, #4 * 8]			\n\t"
			"ldr r10, [r0, #4 * 9]			\n\t"
			"ldr r11, [r0, #4 * 10]			\n\t"
			"ldr r12, [r0, #4 * 11]			\n\t"

			"adds r1 , r9 , r1				\n\t"
			"adcs r2, r10, r2				\n\t"
			"adcs r3, r11, r3				\n\t"
			"adcs r4, r12, r4				\n\t"

			"ldr r9, [r0, #4 * 12]			\n\t"
			"ldr r10, [r0, #4 * 13]			\n\t"
			"ldr r11, [r0, #4 * 14]			\n\t"
			"ldr r12, [r0, #4 * 15]			\n\t"

			"adcs r5, r9, r5				\n\t"
			"adcs r6, r10, r6				\n\t"
			"adcs r7, r11, r7				\n\t"
			"adcs r8, r12, r8				\n\t"

			"vmov r12, d8[0]	 			\n\t"//carry
			"adcs r14, r14, r12				\n\t"//carry
			"asr r12, r14,#10					\n\t"

			"stmia r0!, {r1-r8} 			\n\t"

			"ldmia r0!, {r1-r8} 			\n\t"
			"sub r0, r0, #4*8				\n\t"

			"adds r1, r1, r14		\n\t"
			"adcs r2, r2, r12				\n\t"
			"adcs r3, r3, r12				\n\t"
			"adcs r4, r4, r12				\n\t"

			"adcs r5, r5, r12				\n\t"
			"adcs r6, r6, r12				\n\t"
			"adcs r7, r7, r12				\n\t"
			"adcs r8, r8, r12				\n\t"

			"stmia r0!, {r1-r8} 			\n\t"

			"vpop {q4-q7}					\n\t"
			"pop  {r0-r11,pc}				\n\t"

			
	:
	:
	:
	);
}





void __attribute__ ((noinline, naked)) RED512(const felm_t ma, felm_t mb, felm_t mc){
	asm(
			//general purpose registers (r0-r12, r14)
			"push  {r4-r11,lr}				\n\t"
			"vpush {q4-q7}					\n\t"

			///////////////////////////////////////////////
			"ldr r3, [r1, #4 * 7]			\n\t"//M0			//ROUND1		
			"vldr.64 d0, [r1, #4*8]		\n\t"
			"ldr r4, [r1, #4 * 8]			\n\t"//M1
			"vldr.64 d1, [r1, #4*10]	\n\t"
			"ldr r5, [r1, #4 * 9]			\n\t"//M1
			"vldr.64 d2, [r1, #4*12]	\n\t"
			"ldr r9,  [r0, #4 * 7]			\n\t"//LD IN0
			"vldr.64 d3, [r1, #4*14]	\n\t"
			"ldr r10, [r0, #4 * 8]			\n\t"//LD IN1
			"vtrn.32 q0, q1 			\n\t"
			"ldr r6, [r0, #4 * 0]			\n\t"//Q0
			"vldr.64 d8, [r0, #4*16]		\n\t"
			"ldr r7, [r0, #4 * 1]			\n\t"//Q1
			"vldr.64 d9, [r0, #4*18]	\n\t"
			"ldr r8, [r0, #4 * 2]			\n\t"//Q2
			"vldr.64 d10, [r0, #4*20]	\n\t"
			"mov r11, #0					\n\t"			
			"vldr.64 d11, [r0, #4*22]	\n\t"
			"mov r12, #0					\n\t"			
			"vtrn.32 q4, q5 			\n\t"
			"mov r14, #0					\n\t"	
			"vmov.i64     q6, 0xFFFFFFFF         	\n\t"// mast 1
			"umlal r9, r12, r3, r6			\n\t"//#1
			"vshr.u64     q6, q6, #31		        \n\t"
			"umaal r10, r12, r4, r6			\n\t"//#2
			"vmull.u32 q13, d8, d12[0] 	\n\t"	//A4A0 X B0B0
			"umaal r10, r14, r3, r7			\n\t"//#3
			"vmull.u32 q12, d10, d12[0] 	\n\t"	//A5A1 X B0B0
			"vmov d4[0], r10 				\n\t"
			"vmull.u32 q11, d9, d12[0] 	\n\t"	//A6A2 X B0B0
			"str r9,  [r0, #4 * 7]			\n\t"//LD IN0
			"vmull.u32 q10, d11, d12[0] 	\n\t"	//A7A3 X B0B0		
			"str r10, [r0, #4 * 8]			\n\t"//LD IN1
			"ldr r9, [r0, #4 * 9]			\n\t"//LD IN2
			"umaal r9, r11, r5, r6			\n\t"//#3
			"umaal r9, r12, r4, r7			\n\t"//#4
			"umaal r9, r14, r3, r8			\n\t"//#5
			"vmov d4[1], r9 				\n\t"//NOW NEON

//
			"str r9,  [r0, #4 * 9]			\n\t"//LD IN0			
			"vmlal.u32 q13, d0, d4[0] 	\n\t"	//A4A0 X B0B0
			"mov r10, r14					\n\t"
			"vmlal.u32 q12, d2, d4[0] 	\n\t"	//A5A1 X B0B0
			"ldr r3, [r1, #4 * 10]			\n\t"//M
			"vmlal.u32 q11, d1, d4[0] 	\n\t"	//A6A2 X B0B0
			"ldr r9, [r0, #4 * 10]			\n\t"//IN2
			"vmlal.u32 q10, d3, d4[0] 	\n\t"	//A7A3 X B0B0
			"umaal r9, r10, r3, r6			\n\t"//
			"veor q6,q6,q6 				\n\t"
			"umaal r9, r11, r5, r7			\n\t"//
			"veor q7,q7,q7 				\n\t"
			"umaal r9, r12, r4, r8			\n\t"//
			"veor q8,q8,q8 				\n\t"
			"str r9, [r0, #4 * 10]			\n\t"
			"veor q9,q9,q9 				\n\t"
			"ldr r4, [r1, #4 * 11]			\n\t"//M
			"vtrn.32 q13, q6 			\n\t"
			"ldr r9, [r0, #4 * 11]			\n\t"//IN2
			"vtrn.32 q12, q7 			\n\t"
			"umaal r9, r10, r4, r6			\n\t"//
			"vtrn.32 q11, q8 			\n\t"
			"umaal r9, r11, r3, r7			\n\t"//
			"vtrn.32 q10, q9 			\n\t"
			"umaal r9, r12, r5, r8			\n\t"//
			"vadd.i64  q12, q12, q6 	\n\t"	//q
			"str r9, [r0, #4 * 11]			\n\t"
			"vadd.i64  q11, q11, q7 	\n\t"	//q
			"ldr r5, [r1, #4 * 12]			\n\t"//M
			"vadd.i64  q10, q10, q8 	\n\t"	//q
			"ldr r9, [r0, #4 * 12]			\n\t"//IN
			"vqadd.u64 d18, d18, d27 	\n\t"	//d
			"umaal r9, r10, r5, r6			\n\t"//
			"vext.8 d28, d28, d26, #4 	\n\t"
			"umaal r9, r11, r4, r7			\n\t"//
			"vmlal.u32 q12, d0, d4[1] 	\n\t"	//A4A0 X B0B0
			"umaal r9, r12, r3, r8			\n\t"//			
			"vmlal.u32 q11, d2, d4[1] 	\n\t"	//A5A1 X B0B0
			"str r9, [r0, #4 * 12]			\n\t"
			"vmlal.u32 q10, d1, d4[1] 	\n\t"	//A6A2 X B0B0
			"ldr r3, [r1, #4 * 13]			\n\t"//M
			"vmlal.u32 q9 , d3, d4[1] 	\n\t"	//A7A3 X B0B0
			"ldr r9, [r0, #4 * 13]			\n\t"//IN
			"veor q13,q13,q13			\n\t"
			"umaal r9, r10, r3, r6			\n\t"//
			"veor q6 ,q6 ,q6 				\n\t"
			"umaal r9, r11, r5, r7			\n\t"//
			"veor q7 ,q7 ,q7 				\n\t"
			"umaal r9, r12, r4, r8			\n\t"//
			"veor q8 ,q8 ,q8 				\n\t"
			"str r9, [r0, #4 * 13]			\n\t"			
			"vtrn.32 q12, q13 			\n\t"
			"ldr r4, [r1, #4 * 14]			\n\t"//M
			"vtrn.32 q11, q6 			\n\t"
			"ldr r9, [r0, #4 * 14]			\n\t"//IN
			"vtrn.32 q10, q7 			\n\t"
			"umaal r9, r10, r4, r6			\n\t"//
			"vtrn.32 q9 , q8 			\n\t"
			"umaal r9, r11, r3, r7			\n\t"//
			"vadd.i64  q11, q11, q13 	\n\t"	//q
			"umaal r9, r12, r5, r8			\n\t"//
			"vadd.i64  q10, q10, q6 	\n\t"	//q
			"str r9, [r0, #4 * 14]			\n\t"
			"vadd.i64  q9 , q9 , q7 	\n\t"	//q
			"ldr r5, [r1, #4 * 15]			\n\t"//M
			"vqadd.u64 d16, d16, d25 	\n\t"	//d
			"ldr r9, [r0, #4 * 15]			\n\t"//IN
			"vext.8 d28, d28, d24, #4 	\n\t"
			"umaal r9, r10, r5, r6			\n\t"//
			"umaal r9, r11, r4, r7			\n\t"//
			"umaal r9, r12, r3, r8			\n\t"//
			"str r9, [r0, #4 * 15]			\n\t"			
			"umaal r10, r11, r5, r7			\n\t"//
			"umaal r10, r12, r4, r8			\n\t"//
			"str r10, [r1, #4 * 16]			\n\t"//middle
			"umaal r11, r12, r5, r8			\n\t"//
			"str r11, [r1, #4 * 17]			\n\t"//middle
			"str r12, [r1, #4 * 18]			\n\t"//middle
						"ldr r3, [r1, #4 * 7]			\n\t"//M0			//ROUND2
			"ldr r4, [r1, #4 * 8]			\n\t"//M1
			"ldr r5, [r1, #4 * 9]			\n\t"//M1
			"ldr r9 , [r0, #4 * 10]			\n\t"//LD IN0
			"ldr r10, [r0, #4 * 11]			\n\t"//LD IN1
			"ldr r6, [r0, #4 * 3]			\n\t"//Q0
			"ldr r7, [r0, #4 * 4]			\n\t"//Q1
			"ldr r8, [r0, #4 * 5]			\n\t"//Q2	
			"mov r12, #0					\n\t"			
			"mov r14, #0					\n\t"	
			"umlal r9, r12, r3, r6			\n\t"//#1
			"umaal r10, r12, r4, r6			\n\t"//#2
			"umaal r10, r14, r3, r7			\n\t"//#3
			"vmov d5, r9,r10 				\n\t"
			"str r9,  [r0, #4 * 10]			\n\t"//LD IN0
			"ldr r11, [r0, #4 * 12]			\n\t"//LD IN2
			"mov r9, #0						\n\t"			
			"umaal r9, r11, r5, r6			\n\t"//#3
			"umaal r9, r12, r4, r7			\n\t"//#4
			"umaal r9, r14, r3, r8			\n\t"//#5
			"vmov d6[0], r9 				\n\t"//NOW NEON
			
//
"str r10,  [r0, #4 * 11]			\n\t"//LD IN0
"vmlal.u32 q11, d0, d5[0] 	\n\t"	//A4A0 X B0B0
			"str r9,  [r0, #4 * 12]			\n\t"//LD IN0
			"vmlal.u32 q10, d2, d5[0] 	\n\t"	//A5A1 X B0B0
			"mov r10, r14					\n\t"
			"vmlal.u32 q9 , d1, d5[0] 	\n\t"	//A6A2 X B0B0
			"ldr r3, [r1, #4 * 10]			\n\t"//M
			"vmlal.u32 q8 , d3, d5[0] 	\n\t"	//A7A3 X B0B0
			"ldr r9, [r0, #4 * 13]			\n\t"//IN2
			"veor q12,q12,q12			\n\t"
			"umaal r9, r10, r3, r6			\n\t"//
			"veor q13,q13,q13			\n\t"
			"umaal r9, r11, r5, r7			\n\t"//
			"veor q6 ,q6 ,q6 			\n\t"
			"umaal r9, r12, r4, r8			\n\t"//
			"veor q7 ,q7 ,q7 			\n\t"
			"str r9, [r0, #4 * 13]			\n\t"
			"vtrn.32 q11, q12 			\n\t"
			"ldr r4, [r1, #4 * 11]			\n\t"//M
			"vtrn.32 q10, q13			\n\t"
			"ldr r9, [r0, #4 * 14]			\n\t"//IN2
			"vtrn.32 q9 , q6			\n\t"
			"umaal r9, r10, r4, r6			\n\t"//
			"vtrn.32 q8 , q7 			\n\t"
			"umaal r9, r11, r3, r7			\n\t"//
			"vadd.i64  q10, q10, q12 	\n\t"	//q
			"umaal r9, r12, r5, r8			\n\t"//
			"vadd.i64  q9 , q9 , q13	\n\t"	//q
			"str r9, [r0, #4 * 14]			\n\t"
			"vadd.i64  q8 , q8 , q6 	\n\t"	//q
			"ldr r5, [r1, #4 * 12]			\n\t"//M
			"vqadd.u64 d14, d14, d23 	\n\t"	//d
			"ldr r9, [r0, #4 * 15]			\n\t"//IN
			"vext.8 d29, d29, d22, #4 	\n\t"
			"umaal r9, r10, r5, r6			\n\t"//
			"vmlal.u32 q10, d0, d5[1] 	\n\t"	//A4A0 X B0B0
			"umaal r9, r11, r4, r7			\n\t"//
			"vmlal.u32 q9 , d2, d5[1] 	\n\t"	//A5A1 X B0B0
			"umaal r9, r12, r3, r8			\n\t"//
			"vmlal.u32 q8 , d1, d5[1] 	\n\t"	//A6A2 X B0B0
			"str r9, [r0, #4 * 15]			\n\t"
			"vmlal.u32 q7 , d3, d5[1] 	\n\t"	//A7A3 X B0B0
			"ldr r3, [r1, #4 * 13]			\n\t"//M
			"veor q11,q11,q11			\n\t"
			"ldr r9, [r1, #4 * 16]			\n\t"//IN
			"veor q12,q12,q12			\n\t"
			"umaal r9, r10, r3, r6			\n\t"//
			"veor q13,q13,q13			\n\t"
			"umaal r9, r11, r5, r7			\n\t"//
			"veor q6 ,q6 ,q6 			\n\t"
			"umaal r9, r12, r4, r8			\n\t"//
			"vtrn.32 q10, q11			\n\t"
			"str r9, [r1, #4 * 16]			\n\t"
			"vtrn.32 q9 , q12			\n\t"
			"ldr r4, [r1, #4 * 14]			\n\t"//M
			"vtrn.32 q8 , q13			\n\t"
			"ldr r9, [r1, #4 * 17]			\n\t"//IN
			"vtrn.32 q7 , q6 			\n\t"
			"umaal r9, r10, r4, r6			\n\t"//
			"vadd.i64  q9 , q9 , q11	\n\t"	//q
			"umaal r9, r11, r3, r7			\n\t"//
			"vadd.i64  q8 , q8 , q12	\n\t"	//q
			"umaal r9, r12, r5, r8			\n\t"//
			"vadd.i64  q7 , q7 , q13 	\n\t"	//q
			"str r9, [r1, #4 * 17]			\n\t"
			"vqadd.u64 d12, d12, d21 	\n\t"	//d
			"ldr r5, [r1, #4 * 15]			\n\t"//M
			"vext.8 d29, d29, d20, #4 	\n\t"
			"ldr r9, [r1, #4 * 18]			\n\t"//IN
			"vmlal.u32 q9 , d0, d6[0] 	\n\t"	//A4A0 X B0B0
			"umaal r9, r10, r5, r6			\n\t"//
			"vmlal.u32 q8 , d2, d6[0] 	\n\t"	//A5A1 X B0B0
			"umaal r9, r11, r4, r7			\n\t"//
			"vmlal.u32 q7 , d1, d6[0] 	\n\t"	//A6A2 X B0B0
			"umaal r9, r12, r3, r8			\n\t"//
			"vmlal.u32 q6 , d3, d6[0] 	\n\t"	//A7A3 X B0B0
			"str r9, [r1, #4 * 18]			\n\t"
			"veor q10,q10,q10			\n\t"
			"umaal r10, r11, r5, r7			\n\t"//
			"veor q11,q11,q11			\n\t"
			"umaal r10, r12, r4, r8			\n\t"//
			"veor q12,q12,q12			\n\t"
			"str r10, [r1, #4 * 19]			\n\t"//middle
			"veor q13,q13,q13			\n\t"
			"umaal r11, r12, r5, r8			\n\t"//
			"vtrn.32 q9 , q10			\n\t"
			"str r11, [r1, #4 * 20]			\n\t"//middle
			"vtrn.32 q8 , q11			\n\t"
			"str r12, [r1, #4 * 21]			\n\t"//middle
			"vtrn.32 q7 , q12			\n\t"
			"ldr r3, [r1, #4 * 7]			\n\t"//M0	//ROUND3		
			"vtrn.32 q6 , q13			\n\t"
			"ldr r4, [r1, #4 * 8]			\n\t"//M1
			"vadd.i64  q8 , q8 , q10	\n\t"	//q
			"ldr r5, [r0, #4 * 6]			\n\t"//Q0
			"vadd.i64  q7 , q7 , q11 	\n\t"	//q
			"ldr r6, [r0, #4 * 7]			\n\t"//Q1
			"vadd.i64  q6 , q6 , q12	\n\t"	//q
			"ldr r7, [r0, #4 * 13]			\n\t"//IN0
			"vqadd.u64 d26, d26, d19 	\n\t"	//d
			"ldr r8, [r0, #4 * 14]			\n\t"//IN1
			"vext.8 d30, d30, d18, #4 	\n\t"
			"mov r9, #0						\n\t"			
			"mov r10, #0					\n\t"	
			"umlal r7, r9, r3, r5			\n\t"//#1
			"umaal r8, r9, r4, r5			\n\t"//#2
			"umaal r8, r10, r3, r6			\n\t"//#3
			"vmov d6[1], r7 				\n\t"
			"vmov d7[0], r8 				\n\t"
			"str r7,  [r0, #4 * 13]			\n\t"//LD IN0
			"str r8,  [r0, #4 * 14]			\n\t"//LD IN0
			"ldr r8, [r0, #4 * 15]			\n\t"//IN_correct
			"ldr r14, [r0, #4 * 8]			\n\t"//Q2
			"mov r11, #0					\n\t"	
			"umaal r10, r11, r3, r14		\n\t"//#3
			"ldr r3, [r1, #4 * 9]			\n\t"//M0
			"umaal r10, r8, r4, r6			\n\t"//#2
			"umaal r10, r9, r3, r5			\n\t"//#3
			"vmov d7[1], r10 				\n\t"
			"str r10, [r0, #4 * 15]			\n\t"//IN
//
			"vmlal.u32 q8 , d0, d6[1] 	\n\t"	//A4A0 X B0B0
			"ldr r14, [r1, #4 * 7]			\n\t"//M7
			"vmlal.u32 q7 , d2, d6[1] 	\n\t"	//A5A1 X B0B0
			"ldr r12, [r0, #4 * 9]			\n\t"//Q2
			"vmlal.u32 q6 , d1, d6[1] 	\n\t"	//A6A2 X B0B0
			"ldr r7, [r1, #4 * 16]			\n\t"//IN
			"vmlal.u32 q13, d3, d6[1] 	\n\t"	//A7A3 X B0B0
			"vmov r10, d28[0] 				\n\t"//7 8 9 10 11	
			"veor q9 ,q9 ,q9			\n\t"
			"umaal r7, r8, r14, r12			\n\t"//#4
			"veor q10,q10,q10			\n\t"
			"ldr r4, [r1, #4 * 10]			\n\t"//M			
			"veor q11,q11,q11			\n\t"
			"umaal r7, r9, r4, r5			\n\t"//#4
			"veor q12,q12,q12			\n\t"
			"umaal r7, r10, r3, r6			\n\t"//#5
			"vtrn.32 q8 , q9			\n\t"
			"adds r7, r7, r11				\n\t"
			"vtrn.32 q7 , q10			\n\t"
			"adcs r8, r8, r10				\n\t"
			"vtrn.32 q6 , q11			\n\t"
			"mov r11, #0					\n\t"
			"vtrn.32 q13, q12			\n\t"
			"adcs r11, r11, #0				\n\t"
			"vadd.i64  q7 , q7 , q9 	\n\t"	//q
			"str r7, [r2, #4 * 0]			\n\t"//IN
			"vadd.i64  q6 , q6 , q10	\n\t"	//q
			"ldr r12, [r0, #4 * 10]			\n\t"//Q2
			"vadd.i64  q13, q13, q11	\n\t"	//q
			"ldr r7, [r1, #4 * 17]			\n\t"//IN
			"vqadd.u64 d24, d24, d17 	\n\t"	//d
			"vmov r10, d28[1] 				\n\t"//7 8 9 10 11	
			"vext.8 d30, d30, d16, #4 	\n\t"
			"umaal r7, r8, r14, r12			\n\t"//#4
			"vmlal.u32 q7 , d0, d7[0] 	\n\t"	//A4A0 X B0B0
			"ldr r3, [r1, #4 * 11]			\n\t"//M			
			"vmlal.u32 q6 , d2, d7[0] 	\n\t"	//A5A1 X B0B0
			"umaal r7, r9, r3, r5			\n\t"//#4
			"vmlal.u32 q13, d1, d7[0] 	\n\t"	//A6A2 X B0B0
			"umaal r7, r10, r4, r6			\n\t"//#5
			"vmlal.u32 q12, d3, d7[0] 	\n\t"	//A7A3 X B0B0
			"adds r8, r8, r11				\n\t"
			"veor q8 ,q8 ,q8			\n\t"
			"mov r11, #0					\n\t"
			"veor q9 ,q9 ,q9			\n\t"
			"adcs r11, r11, #0				\n\t"
			"veor q10,q10,q10			\n\t"
			"adds r8, r8, r10				\n\t"
			"veor q11,q11,q11			\n\t"
			"adcs r11, r11, #0				\n\t"
			"vtrn.32 q7 , q8			\n\t"
			"str r7, [r2, #4 * 1]			\n\t"//IN
			"vtrn.32 q6 , q9			\n\t"
			"ldr r12, [r0, #4 * 11]			\n\t"//Q2
			"vtrn.32 q13, q10			\n\t"
			"ldr r7, [r1, #4 * 18]			\n\t"//IN
			"vtrn.32 q12, q11			\n\t"
			"vmov r10, d29[0] 				\n\t"//7 8 9 10 11	
			"vadd.i64  q6 , q6 , q8		\n\t"	//q
			"umaal r7, r8, r14, r12			\n\t"//#4
			"vadd.i64  q13, q13, q9		\n\t"	//q
			"ldr r4, [r1, #4 * 12]			\n\t"//M			
			"vadd.i64  q12, q12, q10	\n\t"	//q
			"umaal r7, r9, r4, r5			\n\t"//#4
			"vqadd.u64 d22, d22, d15 	\n\t"	//d
			"umaal r7, r10, r3, r6			\n\t"//#5
			"vext.8 d31, d31, d14, #4 	\n\t"
			"adds r8, r8, r11				\n\t"
			"vmlal.u32 q6 , d0, d7[1] 	\n\t"	//A4A0 X B0B0
			"mov r11, #0					\n\t"
			"vmlal.u32 q13, d2, d7[1] 	\n\t"	//A5A1 X B0B0
			"adcs r11, r11, #0				\n\t"
			"vmlal.u32 q12, d1, d7[1] 	\n\t"	//A6A2 X B0B0
			"adds r8, r8, r10				\n\t"
			"vmlal.u32 q11, d3, d7[1] 	\n\t"	//A7A3 X B0B0
			"adcs r11, r11, #0				\n\t"
			"veor q7, q7, q7			\n\t"
			"str r7, [r2, #4 * 2]			\n\t"//IN
			"veor q8 ,q8 ,q8			\n\t"
			"ldr r12, [r0, #4 * 12]			\n\t"//Q2
			"veor q9 ,q9 ,q9			\n\t"
			"ldr r7, [r1, #4 * 19]			\n\t"//IN
			"veor q10,q10,q10			\n\t"
			"vmov r10, d29[1] 				\n\t"//7 8 9 10 11	
			"vtrn.32 q6 , q7			\n\t"
			"umaal r7, r8, r14, r12			\n\t"//#4
			"vtrn.32 q13, q8			\n\t"
			"ldr r3, [r1, #4 * 13]			\n\t"//M			
			"vtrn.32 q12, q9			\n\t"
			"umaal r7, r9, r3, r5			\n\t"//#4
			"vtrn.32 q11, q10			\n\t"
			"umaal r7, r10, r4, r6			\n\t"//#5
			"vadd.i64  q13, q13, q7		\n\t"	//q
			"adds r8, r8, r11				\n\t"
			"vadd.i64  q12, q12, q8		\n\t"	//q
			"mov r11, #0					\n\t"
			"vadd.i64  q11, q11, q9		\n\t"	//q
			"adcs r11, r11, #0				\n\t"
			"vqadd.u64 d20, d20, d13 	\n\t"	//d
			"adds r8, r8, r10				\n\t"
			"vext.8 d31, d31, d12, #4 	\n\t"
			"adcs r11, r11, #0				\n\t"
			"vldr.64 d4, [r0, #4*24]		\n\t"
			"str r7, [r2, #4 * 3]			\n\t"//IN
			"vldr.64 d5, [r0, #4*26]	\n\t"
			"ldr r12, [r0, #4 * 13]			\n\t"//Q2
			"vldr.64 d6, [r0, #4*28]	\n\t"
			"ldr r7, [r1, #4 * 20]			\n\t"//IN
			"vldr.64 d7, [r0, #4*30]	\n\t"
			"vmov r10, d30[0] 				\n\t"//7 8 9 10 11	
			"vtrn.32 q2, q3 			\n\t"
			"umaal r7, r8, r14, r12			\n\t"//#4
			"vmov.i64     q4, 0xFFFFFFFF      	\n\t"// mast 1
			"ldr r4, [r1, #4 * 14]			\n\t"//M			
			"vshr.u64     q4, q4, #31		    \n\t"
			"umaal r7, r9, r4, r5			\n\t"//#4
			"vmlal.u32 q13, d4, d8[0] 	\n\t"	//A4A0 X B0B0
			"umaal r7, r10, r3, r6			\n\t"//#5
			"vmlal.u32 q12, d6, d8[0] 	\n\t"	//A5A1 X B0B0
			"adds r8, r8, r11				\n\t"
			"vmlal.u32 q11, d5, d8[0] 	\n\t"	//A6A2 X B0B0
			"mov r11, #0					\n\t"
			"vmlal.u32 q10, d7, d8[0] 	\n\t"	//A7A3 X B0B0	
			"adcs r11, r11, #0				\n\t"
			"adds r8, r8, r10				\n\t"
			"adcs r11, r11, #0				\n\t"
			"str r7, [r2, #4 * 4]			\n\t"//IN
			"ldr r12, [r0, #4 * 14]			\n\t"//Q2
			"ldr r7, [r1, #4 * 21]			\n\t"//IN
			"vmov r10, d30[1] 				\n\t"//7 8 9 10 11	
			"umaal r7, r8, r14, r12			\n\t"//#4
			"ldr r3, [r1, #4 * 15]			\n\t"//M			
			"umaal r7, r9, r3, r5			\n\t"//#4
			"umaal r7, r10, r4, r6			\n\t"//#5
			"adds r8, r8, r11				\n\t"
			"mov r11, #0					\n\t"
			"adcs r11, r11, #0				\n\t"
			"adds r8, r8, r10				\n\t"
			"adcs r11, r11, #0				\n\t"

			"str r7, [r2, #4 * 5]			\n\t"//IN
			"ldr r12, [r0, #4 * 15]			\n\t"//Q2
			"vmov r7, r10, d31 				\n\t"//7 8 9 10 11	
			"umaal r7, r8, r14, r12			\n\t"//#4
			"umaal r7, r9, r3, r6			\n\t"//#4
			"adds r8, r8, r11				\n\t"
			"mov r11, #0					\n\t"
			"adcs r11, r11, #0				\n\t"
			"adds r8, r8, r9				\n\t"
			"adcs r11, r11, #0				\n\t"
			"adds r8, r8, r10				\n\t"
			"adcs r11, r11, #0				\n\t"
			"str r7, [r2, #4 * 6]			\n\t"//IN
			"str r8, [r2, #4 * 7]			\n\t"//IN
		

//8 9 11 have values
//////////////////////////////////////////////////////////////////////////////////

			

/////////////////////////////////////////////////////////////////////////////////////////
			"vmov r3, r4, d26 			\n\t"
			"vmov r5, r6, d24 			\n\t"
			"vmov r7, r8, d22 			\n\t"
			"vmov r9, r10, d20 			\n\t"

			"adcs r3, r3, r11 			\n\t"
			"adcs r4, r4, r5 			\n\t"
			"adcs r5, r6, r7 			\n\t"
			"adcs r6, r8, r9 			\n\t"

			"str r3, [r2, #4 * 8]			\n\t"//IN0
			"str r4, [r2, #4 * 9]			\n\t"//IN0
			"str r5, [r2, #4 * 10]			\n\t"//IN0
			"str r6, [r2, #4 * 11]			\n\t"//IN0

			"vmov r3, r4, d27 			\n\t"
			"vmov r5, r6, d25 			\n\t"
			"vmov r7, r8, d23 			\n\t"
			"vmov r9, r11, d21 			\n\t"

			"adcs r3, r3, r10 			\n\t"
			"adcs r4, r4, r5 			\n\t"
			"adcs r5, r6, r7 			\n\t"
			"adcs r6, r8, r9 			\n\t"

			"str r3, [r2, #4 * 12]			\n\t"//IN0
			"str r4, [r2, #4 * 13]			\n\t"//IN0
			"str r5, [r2, #4 * 14]			\n\t"//IN0
			"str r6, [r2, #4 * 15]			\n\t"//IN0

			"vpop {q4-q7}					\n\t"
			"pop  {r4-r11,pc}				\n\t"
	:
	:
	:
	);
}


static unsigned int modulo_p503[64]={0x0,  0x0,  0x0,  0x0,  0x0,  0x0,  0x0,  0xAC000000,  0x2211E7A0,  0x13085BDA,  0x7B7E7DAF,  0x1B9BF6C8,  0xDA77A4D0,  0x6045C6BD,  0x41811E1E,  0x4066F5};



void fpmul_mont(const felm_t ma, const felm_t mb, felm_t mc)
{ // Multiprecision multiplication, c = a*b mod p.
    dfelm_t temp = {0};
    int k;
    //mp_mul(ma, mb, temp, NWORDS_FIELD);
    MUL512(ma, mb, temp);
	RED512(temp,modulo_p503,mc);

    //rdc_mont(temp, mc);
}



void fpsqr_mont(const felm_t ma, felm_t mc)
{ // Multiprecision squaring, c = a^2 mod p.
    dfelm_t temp = {0};

    //mp_mul(ma, ma, temp, NWORDS_FIELD);
    MUL512(ma,ma,temp);
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
    
    //mp_mul(a[0], b[0], tt1, NWORDS_FIELD);           // tt1 = a0*b0
    //mp_mul(a[1], b[1], tt2, NWORDS_FIELD);           // tt2 = a1*b1
	MUL512(a[0], b[0], tt1);
	MUL512(a[1], b[1], tt2);

    mp_addfast(a[0], a[1], t1);                      // t1 = a0+a1
    mp_addfast(b[0], b[1], t2);                      // t2 = b0+b1
    mask = mp_subfast(tt1, tt2, tt3);                // tt3 = a0*b0 - a1*b1. If tt3 < 0 then mask = 0xFF..F, else if tt3 >= 0 then mask = 0x00..0
    for (i = 0; i < NWORDS_FIELD; i++) {
        ADDC(borrow, tt3[NWORDS_FIELD+i], ((digit_t*)PRIME)[i] & mask, borrow, tt3[NWORDS_FIELD+i]);
    }
    //rdc_mont(tt3, c[0]);                             // c[0] = a0*b0 - a1*b1
	RED512(tt3,modulo_p503,c[0]);
    mp_addfastx2(tt1, tt2, tt1);                     // tt1 = a0*b0 + a1*b1
    //mp_mul(t1, t2, tt2, NWORDS_FIELD);               // tt2 = (a0+a1)*(b0+b1)
	MUL512(t1, t2, tt2);
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