; uint64_t bn_add_inner(uint64_t* r, const uint64_t* a, uint64_t a_size, const uint64_t* b, uint64_t b_size); 
; calling convention: rdi = r, rsi = a, rdx = a_size, rcx = b; r8 = b_size 
; returns possible carry 

section .text 
global bn_add_inner 


bn_add_inner: 

	mov r10, r8 
	mov r11, rdx 
	sub r11, r8

	xor r9, r9	; r9 = i = 0
	clc

	test r10, r10 	; check if b_size == 0
	jz .start_prop

; add a and b
.loop_add: 
	mov rax, [rsi]	; rax = a[i] 
	adc rax, [rcx] 	; rax += b[i] + C 
	mov [rdi], rax 	; r[i] = rax 

	lea     rsi, [rsi + 8]
    lea     rcx, [rcx + 8]
    lea     rdi, [rdi + 8]

	dec r10 
	jnz .loop_add

.start_prop:
	setc al 

	test r11, r11 
	jz	.done 	; if i >= a_size, return 

	shr al, 1

.loop_prop: 
	mov r10, [rsi] 
	adc r10 , 0 	; add carry 
	mov [rdi], r10 	; r[i] = rax 

	lea     rsi, [rsi + 8]
    lea     rdi, [rdi + 8]

	dec r11
	jnz .loop_prop 

	setc al
.done: 
	movzx rax,al
	ret 


