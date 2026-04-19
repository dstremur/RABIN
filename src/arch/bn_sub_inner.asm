; uint64_t bn_sub_inner(uint64_t* r, const uint64_t* a, uint64_t a_size, const uint64_t* b, uint64_t b_size); 
; calling convention: rdi = r, rsi = a, rdx = a_size, rcx = b; r8 = b_size 
; returns final borrow 

section .text 
global bn_sub_inner 

bn_sub_inner: 
	push rbx 

	mov r10, r8
	mov r11, rdx 
	sub r11, r8		; r11 = a_size - b_size 

	xor r9, r9 
	clc 

	test r10, r10 
	jz .start_prop 

; subtract b from a 
.loop_sub: 
	mov rax, [rsi + r9*8] 
	sbb rax, [rcx + r9*8] 
	mov [rdi + r9*8], rax 

	inc r9 
	dec r10 
	jnz .loop_sub 

.start_prop: 
	setc bl 

	test r11, r11 
	jz .done 

	shr bl, 1 

.loop_prop:
	mov rax, [rsi + r9*8] 
	sbb rax, 0 
	mov [rdi + r9*8], rax 

	inc r9
	dec r11 
	jnz .loop_prop 

	setc bl 

.done: 
	movzx rax, bl 
	pop rbx 
	ret 
