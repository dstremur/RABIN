; uint64_t bn_add_inner(uint64_t* r, const uint64_t* a, uint64_t a_size, const uint64_t* b, uint64_t b_size); 
; calling convention: rdi = r, rsi = a, rdx = a_size, rcx = b; r8 = b_size 
; returns possible carry 

section .text 
global bn_add_inner 


bn_add_inner: 
	push rbx 

	mov r10, r8 
	mov r11, rdx 
	sub r11, r8

	xor r9, r9	; r9 = i = 0
	clc

	test r10, r10 	; check if b_size == 0
	jz .start_prop

; add a and b
.loop_add: 
	mov rax, [rsi + r9*8]	; rax = a[i] 
	adc rax, [rcx + r9*8] 	; rax += b[i] + C 
	mov [rdi + r9*8], rax 	; r[i] = rax 

	inc r9 
	dec r10 
	jnz .loop_add

.start_prop:
	setc bl 

	test r11, r11 
	jz	.done 	; if i >= a_size, return 

	shr bl, 1

.loop_prop: 
	mov rax, [rsi + r9*8] 
	adc rax , 0 	; add carry 
	mov [rdi + r9*8], rax 	; r[i] = rax 

	inc r9 
	dec r11
	jnz .loop_prop 

	setc bl
.done: 
	movzx rax, bl
	pop rbx 

	ret 


