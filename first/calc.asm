;

.data
res12		dd		0
res12_1		dd		0
res12_2		dd		0
res12_3		dd		0

res7		dd		0
res7_1		dd		0
res7_2		dd		0
res7_3		dd		0

res2		dd		0
res2_1		dd		0
res2_2		dd		0
res2_3		dd		0

tmp			dd		0
tmp_1		dd		0
tmp_2		dd		0
tmp_3		dd		0

.code
; !!! --- current function is only for one thread (distance_sse_v3) --- !!!
distance_sse_v4 proc
	; load parameters
	movaps      xmm0, xmmword ptr [rcx]	; xmm0 = vec1
	movaps		xmm1, xmmword ptr [rdx]	; xmm1 = vec2

	; start funstion and load params
;	mov         qword ptr [rsp+10h],rdx  
;	mov         qword ptr [rsp+8],rcx  
;	push        rbp  
;	push        rdi  
;	sub         rsp,838h  
;	lea         rbp,[rsp+30h]  
;	mov         rdi,rsp  
;	mov         ecx,20Eh  
;	mov         eax,0CCCCCCCCh  
;	rep stos    dword ptr [rdi]  
;	mov         rcx,qword ptr [rsp+858h]  
;;	mov         rax,qword ptr [__security_cookie (07FF62D7EB010h)]  
;	xor         rax,rbp  
;	mov         qword ptr [rbp+7F8h],rax

	; __m128 res2 = _mm_movehl_ps(*vec1, *vec1); // px, py, px, py
	movaps		xmm2,xmm0
	movhlps     xmm2,xmm0	; xmm2 = res2
	movaps		xmmword ptr [res2],xmm2	; store res2 = res2

	; 	__m128 res1 = _mm_sub_ps(*vec1, *vec2); // where 1 - vx, 2 - vy, 3 - wx(t1(c1)), 4 - wy(t2(c1))
	movaps		xmm3,xmm0
	subps       xmm3,xmm1	; xmm3 = res1
	movaps		xmm7,xmm3	; store res1 => mxx7 = res1

	; 	__m128 res3 = _mm_movelh_ps(res1, res1); // vx, vy, vx, vy
	movaps		xmm4,xmm3	; store res1
	movlhps     xmm4,xmm4	; xmm4 = res3

	; __m128 res4 = _mm_sub_ps(res2, *vec1); // t1(c2), t2(c2), 0(unk), 0(unk)
	subps       xmm2,xmm0	; xmm2 (res2) = res4

	; __m128 res5 = _mm_mul_ps(res3, res1); // vx*vx, vy*vy, wx*vx, wy*vy
	mulps       xmm4,xmm3	; xmm4 (res3) = res5

	; __m128 res8 = _mm_shuffle_ps(res1, res4, 78); // for sqrt: t1(c1), t2(c1), t1(c2), t2(c2)
	shufps      xmm3,xmm2,4Eh	; xmm3 (res1) = res8

	; __m128 res6 = _mm_shuffle_ps(res5, res5, 245); // vy*vy, vy*vy, wy*wy, wy*wy
	movaps		xmm5,xmm4
	shufps      xmm5,xmm5,0F5h	; xmm5 = res6

	; __m128 res9 = _mm_mul_ps(res8, res8); // for sqrt: t1(c1)^2, t2(c1)^2, t1(c2)^2, t2(c2)^2
	mulps       xmm3,xmm3	; xmm3 (res8) = res9

	; __m128 res7 = _mm_add_ps(res5, res6); // c1, unk, c2, unk
	addps       xmm4,xmm5	; xmm4 (res5) = res7
	movaps      xmmword ptr [res7],xmm4

	; __m128 res10 = _mm_shuffle_ps(res9, res9, 245); // t2(c1)^2, t2(c1)^2, t2(c2)^2, t2(c2)^2
	movaps		xmm6,xmm3
	shufps      xmm6,xmm6,0F5h	; xmm6 = res10

	; __m128 res11 = _mm_add_ps(res9, res10); // t1 * t1 + t2 * t2 (c1), unk, t1 * t1 + t2 * t2 (c2), unk
	addps       xmm3,xmm6	; xmm3 = res11

	; __m128 res12 = _mm_sqrt_ps(res11); // sqrt(c1), unk, sqrt(c2), unk
	sqrtps      xmm3,xmm3			; xmm3 = res12
	movaps      xmmword ptr [res12],xmm3

	; if (*(float*)&res7 <= 0) {
	; if (0 > *(float*)&res7) than go to next step
	xorps       xmm0,xmm0  
	comiss      xmm0,dword ptr [res7]
	;jnb			next1
	jb			next1
	; return ((float*)&res12);
	;mov			eax, dword ptr [res12]
	;movss		xmm0,dword ptr [res12]
	movss		xmm0,xmm3
	jmp			end_func

next1:

	; if (((float*)&res7)[2] <= *(float*)&res7) {
	movss       xmm0,dword ptr [res7]  
	comiss      xmm0,dword ptr [res7_2]  
	;jnb          next2
	jb			next2
		; return ((float*)&res12)[2];
	;mov       eax,dword ptr [res12_2]  
	shufps		xmm0,xmm3,0E0h
	shufps		xmm0,xmm0,0E6h
	jmp         end_func

next2:
	; //coord b = c1 / c2;
	; coord b = *(float*)&res7 / ((float*)&res7)[2];
	divss       xmm4,dword ptr [res7_2]	; xmm4 (res7) = b

	;__m128 bb = _mm_set1_ps(b);
	shufps      xmm4,xmm4,0		; xmm4 (b) = xmm4 bb (b, b, b, b)

	; //pbx = line_p0->x + b * vx;
	; //pby = line_p0->y + b * vy;
	; res9 = _mm_mul_ps(bb, res1); // vx * b, vy * b, unk, unk
	mulps       xmm4,xmm7	; xmm4 = res9

	; res10 = _mm_add_ps(*vec2, res9); // pbx, pby, unk, unk
	movaps      xmm5,xmmword ptr [res2]	; for next step
	addps       xmm4,xmm1	; xmm4 = res10

	; 	//coord t1 = p->x - pbx;
	; //coord t2 = p->y - pby;
	; res11 = _mm_sub_ps(res2, res10); // t1, t2, unk, unk
	subps       xmm5,xmm4		; xmm5 = res11
	; res12 = _mm_mul_ps(res11, res11); // t1^2, t2^2, unk, unk
	mulps       xmm5,xmm5
	movaps      xmmword ptr [res12],xmm5	; xmm5 = res12

	; 	return sqrt(*((float*)&res12) + ((float*)&res12)[1]);
	addss       xmm5,dword ptr [res12_1]
	sqrtss		xmm5,xmm5

	; prepare for return
	;cvtss2sd	xmm5,xmm5
	;movaps		xmmword ptr [res12],xmm5

	movss		xmm0,xmm5
	;movss		dword ptr [res12],xmm5 ; result return
	;mov			eax,dword ptr [res12]

	
end_func:

	; end function
;	movdqu      xmmword ptr [rsp+20h],xmm0  
;	lea         rcx,[rbp-30h]  
;;	lea         rdx,[__real@447a0000+74Ch (07FF62D7E8040h)]  
;;	call        _RTC_CheckStackVars (07FF62D7D118Bh)  
;	movdqu      xmm0,xmmword ptr [rsp+20h]  
;	mov         rcx,qword ptr [rbp+7F8h]  
;	xor         rcx,rbp  
;;	call        __security_check_cookie (07FF62D7D13CAh)  
;	lea         rsp,[rbp+808h]  
;	pop         rdi  
;	pop         rbp
	
  
	ret
distance_sse_v4 endp


; --------------------------------------------------------------------------------------
; fuctions calculate distance for multithreads (distance_sse_v3)
distance_sse_v5 proc
	; load parameters
	movaps      xmm0, xmmword ptr [rcx]	; xmm0 = vec1
	movaps		xmm1, xmmword ptr [rdx]	; xmm1 = vec2
	mov			rdx,r8

	; __m128 res2 = _mm_movehl_ps(*vec1, *vec1); // px, py, px, py
	movaps		xmm2,xmm0
	movhlps     xmm2,xmm0	; xmm2 = res2
	movaps		xmmword ptr [rdx],xmm2	; store rdx = res2

	; 	__m128 res1 = _mm_sub_ps(*vec1, *vec2); // where 1 - vx, 2 - vy, 3 - wx(t1(c1)), 4 - wy(t2(c1))
	movaps		xmm3,xmm0
	subps       xmm3,xmm1	; xmm3 = res1
	movaps		xmm7,xmm3	; store res1 => mxx7 = res1

	; 	__m128 res3 = _mm_movelh_ps(res1, res1); // vx, vy, vx, vy
	movaps		xmm4,xmm3	; store res1
	movlhps     xmm4,xmm4	; xmm4 = res3

	; __m128 res4 = _mm_sub_ps(res2, *vec1); // t1(c2), t2(c2), 0(unk), 0(unk)
	subps       xmm2,xmm0	; xmm2 (res2) = res4

	; __m128 res5 = _mm_mul_ps(res3, res1); // vx*vx, vy*vy, wx*vx, wy*vy
	mulps       xmm4,xmm3	; xmm4 (res3) = res5

	; __m128 res8 = _mm_shuffle_ps(res1, res4, 78); // for sqrt: t1(c1), t2(c1), t1(c2), t2(c2)
	shufps      xmm3,xmm2,4Eh	; xmm3 (res1) = res8

	; __m128 res6 = _mm_shuffle_ps(res5, res5, 245); // vy*vy, vy*vy, wy*wy, wy*wyxmmword ptr
	movaps		xmm5,xmm4
	shufps      xmm5,xmm5,0F5h	; xmm5 = res6

	; __m128 res9 = _mm_mul_ps(res8, res8); // for sqrt: t1(c1)^2, t2(c1)^2, t1(c2)^2, t2(c2)^2
	mulps       xmm3,xmm3	; xmm3 (res8) = res9

	; __m128 res7 = _mm_add_ps(res5, res6); // c1, unk, c2, unk
	addps       xmm4,xmm5	; xmm4 (res5) = res7
	mov         eax,10h
	imul        rax,rax,1
	movaps      xmmword ptr [rdx[rax]],xmm4
	;movaps      xmmword ptr [res7],xmm4

	; __m128 res10 = _mm_shuffle_ps(res9, res9, 245); // t2(c1)^2, t2(c1)^2, t2(c2)^2, t2(c2)^2
	movaps		xmm6,xmm3
	shufps      xmm6,xmm6,0F5h	; xmm6 = res10

	; __m128 res11 = _mm_add_ps(res9, res10); // t1 * t1 + t2 * t2 (c1), unk, t1 * t1 + t2 * t2 (c2), unk
	addps       xmm3,xmm6	; xmm3 = res11

	; __m128 res12 = _mm_sqrt_ps(res11); // sqrt(c1), unk, sqrt(c2), unk
	sqrtps      xmm3,xmm3			; xmm3 = res12
	mov         eax,10h
	imul        rax,rax,2
	movaps      xmmword ptr [rdx[rax]],xmm3
	;movaps      xmmword ptr [res12],xmm3

	; if (*(float*)&res7 <= 0) {
	; if (0 > *(float*)&res7) than go to next step
	xorps       xmm0,xmm0  
	mov         eax,10h
	imul        rax,rax,1
	comiss      xmm0,dword ptr [rdx[rax]]
	;comiss      xmm0,dword ptr [res7]
	jb			next1
	; return ((float*)&res12);
	movss		xmm0,xmm3
	jmp			end_func

next1:

	; if (((float*)&res7)[2] <= *(float*)&res7) {
	;movss       xmm0,dword ptr [res7]
	movss       xmm0,dword ptr [rdx[rax]]
	;comiss      xmm0,dword ptr [res7_2]
	comiss      xmm0,dword ptr [rdx[rax] + 8]
	;jnb          next2
	jb			next2
		; return ((float*)&res12)[2];
	shufps		xmm0,xmm3,0E0h
	shufps		xmm0,xmm0,0E6h
	jmp         end_func

next2:
	; //coord b = c1 / c2;
	; coord b = *(float*)&res7 / ((float*)&res7)[2];
	;divss       xmm4,dword ptr [res7_2]	; xmm4 (res7) = b
	divss       xmm4,dword ptr [rdx[rax] + 8]	; xmm4 (res7) = b

	;__m128 bb = _mm_set1_ps(b);
	shufps      xmm4,xmm4,0		; xmm4 (b) = xmm4 bb (b, b, b, b)

	; //pbx = line_p0->x + b * vx;
	; //pby = line_p0->y + b * vy;
	; res9 = _mm_mul_ps(bb, res1); // vx * b, vy * b, unk, unk
	mulps       xmm4,xmm7	; xmm4 = res9

	; res10 = _mm_add_ps(*vec2, res9); // pbx, pby, unk, unk
	;movaps      xmm5,xmmword ptr [res2]	; for next step
	movaps      xmm5,xmmword ptr [rdx]	; for next step
	addps       xmm4,xmm1	; xmm4 = res10

	; 	//coord t1 = p->x - pbx;
	; //coord t2 = p->y - pby;
	; res11 = _mm_sub_ps(res2, res10); // t1, t2, unk, unk
	subps       xmm5,xmm4		; xmm5 = res11
	mov         eax,10h
	imul        rax,rax,2
	; res12 = _mm_mul_ps(res11, res11); // t1^2, t2^2, unk, unk
	mulps       xmm5,xmm5
	;movaps      xmmword ptr [res12],xmm5	; xmm5 = res12
	movaps      xmmword ptr [rdx[rax]],xmm5	; xmm5 = res12

	; 	return sqrt(*((float*)&res12) + ((float*)&res12)[1]);
	;addss       xmm5,dword ptr [res12_1]
	addss       xmm5,dword ptr [rdx[rax] + 4]
	sqrtss		xmm5,xmm5

	movss		xmm0,xmm5
	
end_func:

	; end function	
  
	ret

distance_sse_v5 endp


; --------------------------------------------------------------------------------------
; fuctions calculate distance  (distance_sse_v2)
distance_sse_v6 proc
	; load parameters
	movaps      xmm0, xmmword ptr [rcx]	; xmm0 (xmm2) = vec1
	movaps		xmm1, xmmword ptr [rdx]	; xmm1 = vec2

	movaps		xmm2, xmm0			; copy for other steps xmm2 = vec1
	; __m128 res1 = _mm_sub_ps(*vec1, *vec2); // where 1 - vx, 2 - vy, 3 - wx, 4 - wy
	subps		xmm0, xmm1			; xmm0 = res1

	; __m128 tt1 = _mm_movehl_ps(res1, res1); // wx, wy, wx, wy
	movaps		xmm3, xmm0
	movhlps		xmm3, xmm3			; xmm3 = tt1

	; __m128 res2 = _mm_mul_ps(res1, tt1); // c1.1, c1.2, ret_c1 t1*t1, ret_c1 t2*t2
	mulps		xmm3, xmm0			; xmm3 = res2

	; __m128 res3 = _mm_shuffle_ps(res2, res2, 177); // c1.2, c1.1, ret_c1 t2*t2, ret_c1 t1*t1
	movaps		xmm4, xmm3
	shufps      xmm4, xmm4, 0B1h	; xmm4 = res3

	;__m128 res4 = _mm_add_ps(res2, res3); // c1, c1, for sqrt, for sqrt
	addps		xmm4, xmm3			; free xmm3, xmm4 = res4

	; if ((*((float*)&res4)) <= 0) {
	xorps		xmm3, xmm3			; xmm3 = 0
	comiss		xmm3, xmm4
	jb			next1
	; return sqrt(((float*)&res4)[3]);
	shufps		xmm4, xmm4, 0E7h
	sqrtss		xmm0, xmm4
	; movaps		xmm0, xmm4
	jmp			end_func

next1:

	; __m128 res5 = _mm_mul_ps(res1, res1); // vx^2, vy^2, unk, unk
	movaps		xmm3, xmm0
	mulps		xmm3, xmm3			; xmm3 = res5

	; __m128 res51 = _mm_shuffle_ps(res5, res5, 225); // vy^2, vx^2, unk, unk
	movaps		xmm5, xmm3			; xmm5 = res5
	shufps		xmm3, xmm3, 0E1h	; xmm3 = res51

	; __m128 res52 = _mm_add_ss(res5, res51); // c2, unk, unk, unk
	addss		xmm3, xmm5			; free xmm5; xmm3 = res52

	; __m128 res6 = _mm_shuffle_ps(*vec1, *vec1, 78); // p.x, p.y, p1.x, p1.y
	shufps		xmm2, xmm2, 4Eh		; xmm2 = res6

	; if (*((float*)&res52) <= *((float*)&res4)) {
	comiss		xmm4, xmm3
	jb			next2
	; __m128 res7 = _mm_sub_ps(res6, res1); // t1, t2, unk, unk
	subps		xmm2, xmm0			; xmm2 = res7
	; __m128 res8 = _mm_mul_ps(res7, res7);
	mulps		xmm2, xmm2
	; return sqrt(*((float*)&res8) + ((float*)&res8)[1]);
	movaps		xmm0, xmm2
	shufps		xmm2, xmm2, 0E1h
	addss		xmm0, xmm2
	sqrtss		xmm0, xmm0
	jmp			end_func

next2:

	; coord b = *((float*)&res4) / *((float*)&res52);
	divss		xmm4, xmm3			; xmm4 = b
	; __m128 bb = _mm_set1_ps(b);
	shufps		xmm4, xmm4, 0		; xmm4 = bb

	; __m128 res9 = _mm_mul_ps(bb, res1); // vx * b, vy * b, unk, unk
	mulps		xmm0, xmm4			; xmm0 = res9
	; __m128 res10 = _mm_add_ps(*vec2, res9); // pbx, pby, unk, unk
	addps		xmm0, xmm1			; xmm0 = res10

	; __m128 res11 = _mm_sub_ps(res6, res10); // t1, t2, unk, unk
	subps		xmm2, xmm0			; xmm2 = res11

	; __m128 res12 = _mm_mul_ps(res11, res11); // t1^2, t2^2, unk, unk
	mulps		xmm2, xmm2			; xmm2 = res12

	; return sqrt(*((float*)&res12) + ((float*)&res12)[1]);
	movaps		xmm0, xmm2			; xmm0 = res12
	shufps		xmm2, xmm2, 0E1h
	addss		xmm0, xmm2
	sqrtss		xmm0, xmm0


end_func:

	; end function	
	ret
distance_sse_v6 endp

end
