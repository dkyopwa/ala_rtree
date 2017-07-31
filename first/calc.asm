;

.data
res12	dq	0
res12_	dq	0

.code
distance_sse_v4 proc
	; load parameters
	movaps      xmm0, xmmword ptr [rcx]	; xmm0 = vec1
	movaps		xmm1, xmmword ptr [rdx]	; xmm1 = vec2

	; start funstion and load params
	mov         qword ptr [rsp+10h],rdx  
	mov         qword ptr [rsp+8],rcx  
	push        rbp  
	push        rdi  
	sub         rsp,838h  
	lea         rbp,[rsp+30h]  
	mov         rdi,rsp  
	mov         ecx,20Eh  
	mov         eax,0CCCCCCCCh  
	rep stos    dword ptr [rdi]  
	mov         rcx,qword ptr [rsp+858h]  
;	mov         rax,qword ptr [__security_cookie (07FF62D7EB010h)]  
	xor         rax,rbp  
	mov         qword ptr [rbp+7F8h],rax

	; __m128 res2 = _mm_movehl_ps(*vec1, *vec1); // px, py, px, py
	movaps		xmm2,xmm0
	movhlps     xmm2,xmm0	; xmm2 = res2

	; 	__m128 res1 = _mm_sub_ps(*vec1, *vec2); // where 1 - vx, 2 - vy, 3 - wx(t1(c1)), 4 - wy(t2(c1))
	movaps		xmm3,xmm0
	subps       xmm3,xmm1	; xmm3 = res1

	; 	__m128 res3 = _mm_movelh_ps(res1, res1); // vx, vy, vx, vy
	movaps		xmm4,xmm3
	movlhps     xmm4,xmm4	; xmm4 = res3

	; __m128 res4 = _mm_sub_ps(res2, *vec1); // t1(c2), t2(c2), 0(unk), 0(unk)
	subps       xmm2,xmm1	; xmm2 (res2) = res4

	; __m128 res5 = _mm_mul_ps(res3, res1); // vx*vx, vy*vy, wx*vx, wy*vy
	mulps       xmm4,xmm3	; xmm4 (res3) = res5

	; __m128 res8 = _mm_shuffle_ps(res1, res4, 78); // for sqrt: t1(c1), t2(c1), t1(c2), t2(c2)
	shufps      xmm3,xmm2,4Eh	; xmm3 = res8

	; __m128 res6 = _mm_shuffle_ps(res5, res5, 245); // vy*vy, vy*vy, wy*wy, wy*wy
	movaps		xmm5,xmm4
	shufps      xmm5,xmm5,0F5h	; xmm5 = res6

	; __m128 res9 = _mm_mul_ps(res8, res8); // for sqrt: t1(c1)^2, t2(c1)^2, t1(c2)^2, t2(c2)^2
	mulps       xmm3,xmm3	; xmm3 (res8) = res9

	; __m128 res7 = _mm_add_ps(res5, res6); // c1, unk, c2, unk
	addps       xmm4,xmm5	; xmm4 (res5) = res7

	; __m128 res10 = _mm_shuffle_ps(res9, res9, 245); // t2(c1)^2, t2(c1)^2, t2(c2)^2, t2(c2)^2
	movaps		xmm6,xmm3
	shufps      xmm6,xmm6,0F5h	; xmm6 = res10

	; __m128 res11 = _mm_add_ps(res9, res10); // t1 * t1 + t2 * t2 (c1), unk, t1 * t1 + t2 * t2 (c2), unk
	addps       xmm3,xmm6	; xmm3 = res11

	; __m128 res12 = _mm_sqrt_ps(res11); // sqrt(c1), unk, sqrt(c2), unk
	sqrtps      xmm3,xmm3
	movaps      xmmword ptr [res12],xmm3

	mov			rax,0FFFFFFFFh ; result return

	; end function
	movdqu      xmmword ptr [rsp+20h],xmm0  
	lea         rcx,[rbp-30h]  
;	lea         rdx,[__real@447a0000+74Ch (07FF62D7E8040h)]  
;	call        _RTC_CheckStackVars (07FF62D7D118Bh)  
	movdqu      xmm0,xmmword ptr [rsp+20h]  
	mov         rcx,qword ptr [rbp+7F8h]  
	xor         rcx,rbp  
;	call        __security_check_cookie (07FF62D7D13CAh)  
	lea         rsp,[rbp+808h]  
	pop         rdi  
	pop         rbp  
	ret
distance_sse_v4 endp
end
