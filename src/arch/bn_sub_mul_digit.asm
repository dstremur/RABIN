; Target: x86_64, System V ABI
; uint64_t bn_sub_mul_digit(uint64_t* r, const uint64_t* b, uint64_t a_limb, uint64_t len);

section .text
global bn_sub_mul_digit

bn_sub_mul_digit:
    ; rdi = C_limbs
    ; rsi = A_limbs
    ; rdx = b_k (a_limb)
    ; rcx = len

    mov r9, rdx         ; Move b_k to r9 so we can use rdx for the mul instruction
    xor r8, r8          ; r8 will be our unified carry/borrow accumulator

    test rcx, rcx
    jz .done

align 16
.loop:
    mov rax, [rsi]
    mul r9              ; rdx:rax = A[i] * b_k
    
    add rax, r8         ; Add previous carry to the low product
    adc rdx, 0          ; Add any carry-out to the high product (rdx)
    
    sub [rdi], rax      ; Subtract the low part from C[i]
    adc rdx, 0          ; WIZARDRY: If 'sub' caused a borrow, CF is set. 
                        ; This adds 1 to the high product for the next loop!
    
    mov r8, rdx         ; Save unified carry/borrow for next iteration

    lea rsi, [rsi + 8]  ; Advance A_limbs pointer
    lea rdi, [rdi + 8]  ; Advance C_limbs pointer
    dec rcx
    jnz .loop

.done:
    mov rax, r8         ; Return the final carry/borrow to C
    ret
