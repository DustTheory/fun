            global      naive
            extern      printf
            section     .text

naive:    
            push        rbp
            mov         rbp, rsp
            mov         QWORD [rbp-8], rdi      ; needle
            mov         DWORD [rbp-12], esi     ; needle length
            mov         QWORD [rbp-24], rdx     ; haystack
            mov         DWORD [rbp-16], ecx     ; haystack length
            dec         esi
            sub         ecx, esi
            jl          .afterSearch
            inc         esi

            xor         rax, rax                ; rax is our main loop index
.outerLoop:
            xor         r8, r8                  ; r8 is our inner loop index
.innerLoop:
            mov         r9, rdx
            add         r9, rax
            add         r9, r8                  ; r9 = &haystack[rax + r8] = rdx + rax + r8
            mov         r9b, BYTE [r9]
            mov         r10, rdi
            add         r10, r8                 ; r10 = &needle[r8] = rdi + r8
            mov         r10b, BYTE[r10]
            cmp         r9b, r10b                 ; if( needle[i] != haystack[r8]) break; (jump to outerLoopEnd)
            jne         .outerLoopEnd
            inc         r8
            cmp         esi, r8d                 ; else if(r8 == needle_length) rax = rdx + rax + r8 (haystack[rax + r8]), goto returnrax;
            dec         r8
            jne         .innerLoopEnd
            add         rax, rdx
            add         eax, r8d
            jmp         .returnrax
.innerLoopEnd:
            inc         r8
            cmp         r8d, esi
            jne         .innerLoop
.outerLoopEnd:
            inc         rax
            cmp         eax, ecx
            jne         .outerLoop
            
.afterSearch:            
            mov         rax, 0                  ; set rax to NULL
.returnrax:
            pop         rbp
            ret

            section     .rodata
format:
            db           "%s", 10, 0