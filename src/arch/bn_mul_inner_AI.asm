; Target: x86_64, Intel Syntax
; ABI: System V -> rdi = r, rsi = b, rdx = a_limb, rcx = len
;
; Computes: r[i] += a_limb * b[i] for i = 0..len-1
; Adds carry to r[len], returns carry in rax
;
; Optimizations:
;   - 8x unrolled main loop (was 4x)
;   - Alternating carry registers r8/r9 (eliminates mov per limb)
;   - ADX dual carry chains (CF via adcx, OF via adox)
;   - Carry consumption between blocks (fixes lost-CF bug)
;   - Flag-preserving loop control in tail (lea+jrcxz+jmp)
;   - Proper return value in rax

section .text
global bn_mul_add_inner

bn_mul_add_inner:
    test    rcx, rcx
    jz      .done

    mov     rax, rcx
    and     rax, 7                ; remainder = len % 8
    shr     rcx, 3                ; blocks = len / 8

    xor     r8, r8                ; carry = 0
    xor     r11, r11              ; permanent zero (for carry consumption + flag clear)

    test    rcx, rcx              ; also clears CF and OF
    jz      .tail_setup

align 16
.loop8:
    xor     r11d, r11d            ; clear CF, OF (may be corrupted by dec)

    ; ---- Even limbs: mulx high -> r9 ----
    ; ---- Odd  limbs: mulx high -> r8 ----

    ; Limb 0 (even)
    mulx    r9, r10, [rsi]
    adcx    r10, r8
    adox    r10, [rdi]
    mov     [rdi], r10

    ; Limb 1 (odd)
    mulx    r8, r10, [rsi + 8]
    adcx    r10, r9
    adox    r10, [rdi + 8]
    mov     [rdi + 8], r10

    ; Limb 2 (even)
    mulx    r9, r10, [rsi + 16]
    adcx    r10, r8
    adox    r10, [rdi + 16]
    mov     [rdi + 16], r10

    ; Limb 3 (odd)
    mulx    r8, r10, [rsi + 24]
    adcx    r10, r9
    adox    r10, [rdi + 24]
    mov     [rdi + 24], r10

    ; Limb 4 (even)
    mulx    r9, r10, [rsi + 32]
    adcx    r10, r8
    adox    r10, [rdi + 32]
    mov     [rdi + 32], r10

    ; Limb 5 (odd)
    mulx    r8, r10, [rsi + 40]
    adcx    r10, r9
    adox    r10, [rdi + 40]
    mov     [rdi + 40], r10

    ; Limb 6 (even)
    mulx    r9, r10, [rsi + 48]
    adcx    r10, r8
    adox    r10, [rdi + 48]
    mov     [rdi + 48], r10

    ; Limb 7 (odd)
    mulx    r8, r10, [rsi + 56]
    adcx    r10, r9
    adox    r10, [rdi + 56]
    mov     [rdi + 56], r10

    ; Combine carry chains: r8 = r8 + CF + OF
    ; Safe: total carry always fits in 64 bits (r8 + CF + OF < 2^64)
    adcx    r8, r11
    adox    r8, r11

    lea     rdi, [rdi + 64]
    lea     rsi, [rsi + 64]

    dec     rcx                   ; modifies OF, but CF/OF already consumed
    jnz     .loop8

.tail_setup:
    ; r8 = carry, CF = 0, OF = 0 (from carry consumption or test)
    mov     rcx, rax              ; remainder
    jrcxz   .final_carry

.tail_loop:
    mulx    r9, r10, [rsi]
    adcx    r10, r8
    adox    r10, [rdi]
    mov     [rdi], r10
    mov     r8, r9

    lea     rdi, [rdi + 8]
    lea     rsi, [rsi + 8]

    ; Flag-preserving loop control (lea+jrcxz don't touch CF/OF)
    lea     rcx, [rcx - 1]
    jrcxz   .tail_done
    jmp     .tail_loop

.tail_done:
    ; Consume final carries into r8
    adcx    r8, r11
    adox    r8, r11

.final_carry:
    add     [rdi], r8
    mov     rax, r8               ; return carry in rax

.done:
    ret

section .note.GNU-stack noalloc noexec nowrite progbits