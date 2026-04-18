; Target: x86_64, Intel Syntax
; ABI: System V -> rdi = r, rsi = b, rdx = a_limb, rcx = len

section .text
global bn_mul_add_inner

bn_mul_add_inner:
    test    rcx, rcx
    jz      .done_early

    mov     rax, rcx
    and     rax, 3          ; remainder = len % 4
    shr     rcx, 2          ; blocks = len / 4

    xor     r8, r8          ; carry from previous mulx high part

    xor     r9, r9
    adcx    r9, r9          ; CF = 0
    adox    r9, r9          ; OF = 0

	test rcx, rcx
    jz   .tail_setup

align 16
.loop:
    ; limb 0
    mulx    r11, r10, [rsi]
    adcx    r10, r8
    adox    r10, [rdi]
    mov     [rdi], r10
    mov     r8, r11

    ; limb 1
    mulx    r11, r10, [rsi + 8]
    adcx    r10, r8
    adox    r10, [rdi + 8]
    mov     [rdi + 8], r10
    mov     r8, r11

    ; limb 2
    mulx    r11, r10, [rsi + 16]
    adcx    r10, r8
    adox    r10, [rdi + 16]
    mov     [rdi + 16], r10
    mov     r8, r11

    ; limb 3
    mulx    r11, r10, [rsi + 24]
    adcx    r10, r8
    adox    r10, [rdi + 24]
    mov     [rdi + 24], r10
    mov     r8, r11

    lea     rdi, [rdi + 32]
    lea     rsi, [rsi + 32]

    loop    .loop

.tail_setup:
    mov rcx, rax	; switch to rcx for loop to work 
    jrcxz  .done

.tail_loop:
    mulx    r11, r10, [rsi]
    adcx    r10, r8
    adox    r10, [rdi]
    mov     [rdi], r10
    mov     r8, r11

    lea     rdi, [rdi + 8]
    lea     rsi, [rsi + 8]

    loop    .tail_loop

.done:
    mov     r10, 0
    adcx    r8, r10
    adox    r8, r10
    add     [rdi], r8

.done_early:
    ret

section .note.GNU-stack noalloc noexec nowrite progbits
