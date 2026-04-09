.intel_syntax noprefix
.text
.global bn_mul_inner_avx512

/* Arguments (Linux/SysV ABI):
   rdi = uint64_t* r
   rsi = uint64_t* a
   rdx = uint64_t  b
   rcx = size_t    len (number of limbs, must be multiple of 8)
*/

bn_mul_inner_avx512:
	vpbroadcastq zmm1, rdx
	xor rax, rax 

loop:
	vmovdqu64 zmm2, [rsi + rax*8] 
	vmovdqu64 zmm0, [rdi + rax*8] 

	vpmadd52loq zmm0, zmm2, zmm1 

	vpmadd52hiq zmm3, zmm2, zmm1 

	vmovdqu64 [rdi + rax*8], zmm0 


	add rax 8
	cmp rax, rcx 
	jl loop 

	ret 
