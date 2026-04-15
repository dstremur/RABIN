; Target: x86_64, Intel Syntax
; ABI: System V -> rdi=r, rsi=b, rdx=a_limb, rcx=len

section .text
global bn_mul_add_inner

bn_mul_add_inner:
    test rcx, rcx
    jz .done_early

    ; No push/pop needed, we only use volatile registers
	push rbx
    clc             ; Clear Carry Flag (CF)
    xor rax, rax    ; Clears Overflow Flag (OF)
    xor r8, r8      ; r8 will hold our Previous High Limb

    align 16        ; Align loop for instruction cache
.loop:
    mov r9, [rsi]
    
    ; mulx: r11:r10 = a_limb * b[i]
    mulx r11, r10, r9

    ; Chain 1 (CF): Add previous high limb to current low limb
    adcx r10, r8
    
    ; Chain 2 (OF): Add destination limb to current low limb
    adox r10, [rdi]
    
    ; Store low limb and carry the new high limb to the next iteration
    mov [rdi], r10
    mov r8, r11
    
    ; Advance pointers
    lea rdi, [rdi + 8]
    lea rsi, [rsi + 8]

    ; THE MAGIC LOOP:
    ; 'lea' decrements rcx without touching CF or OF!
    ; 'jrcxz' jumps out if rcx == 0 without touching CF or OF!
	loop .loop

.done:
    ; Finalize the carries. 
    ; adcx and adox do NOT interfere with each other's flags!
    mov r10, 0
    adcx r8, r10    ; Fold in remaining CF
    adox r8, r10    ; Fold in remaining OF
    add [rdi], r8

	pop rbx

.done_early:
    ret

section .note.GNU-stack noalloc noexec nowrite progbits
