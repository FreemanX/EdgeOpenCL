	.cpu generic+fp+simd
	.file	"util.c"
	.text
	.align	2
	.global	_Z18aligned_block_copyPlS_i
	.type	_Z18aligned_block_copyPlS_i, %function
_Z18aligned_block_copyPlS_i:
.LFB5:
	.cfi_startproc
	cmp	w2, #64
	bmi	.L1
	sub	w2, w2, #64
	and	x6, x2, 4294967232
	add	x6, x6, 64
	add	x6, x1, x6
.L3:
	ldr	x5, [x1]
	str	x5, [x0]
	ldr	x4, [x1,8]
	str	x4, [x0,8]
	ldr	x3, [x1,16]
	str	x3, [x0,16]
	ldr	x2, [x1,24]
	str	x2, [x0,24]
	ldr	x5, [x1,32]
	str	x5, [x0,32]
	ldr	x4, [x1,40]
	str	x4, [x0,40]
	ldr	x3, [x1,48]
	str	x3, [x0,48]
	add	x1, x1, 64
	ldr	x2, [x1,-8]
	str	x2, [x0,56]
	add	x0, x0, 64
	cmp	x6, x1
	bne	.L3
.L1:
	ret
	.cfi_endproc
.LFE5:
	.size	_Z18aligned_block_copyPlS_i, .-_Z18aligned_block_copyPlS_i
	.align	2
	.global	_Z28aligned_block_copy_backwardsPlS_i
	.type	_Z28aligned_block_copy_backwardsPlS_i, %function
_Z28aligned_block_copy_backwardsPlS_i:
.LFB6:
	.cfi_startproc
	add	w3, w2, 7
	cmp	w2, wzr
	csel	w3, w3, w2, lt
	asr	w3, w3, 3
	sbfiz	x3, x3, 3, 32
	sub	x3, x3, #8
	add	x1, x1, x3
	add	x0, x0, x3
	cmp	w2, #64
	bmi	.L5
	sub	w2, w2, #64
	and	x2, x2, 4294967232
	sub	x2, x1, x2
	sub	x2, x2, #64
.L7:
	ldr	x6, [x1]
	ldr	x5, [x1,-8]
	ldr	x4, [x1,-16]
	ldr	x3, [x1,-24]
	str	x6, [x0]
	str	x5, [x0,-8]
	str	x4, [x0,-16]
	str	x3, [x0,-24]
	ldr	x6, [x1,-32]
	ldr	x5, [x1,-40]
	ldr	x4, [x1,-48]
	sub	x1, x1, #64
	ldr	x3, [x1,8]
	str	x6, [x0,-32]
	str	x5, [x0,-40]
	str	x4, [x0,-48]
	str	x3, [x0,-56]
	sub	x0, x0, #64
	cmp	x2, x1
	bne	.L7
.L5:
	ret
	.cfi_endproc
.LFE6:
	.size	_Z28aligned_block_copy_backwardsPlS_i, .-_Z28aligned_block_copy_backwardsPlS_i
	.align	2
	.global	_Z33aligned_block_copy_backwards_bs32PlS_i
	.type	_Z33aligned_block_copy_backwards_bs32PlS_i, %function
_Z33aligned_block_copy_backwards_bs32PlS_i:
.LFB7:
	.cfi_startproc
	add	w3, w2, 7
	cmp	w2, wzr
	csel	w3, w3, w2, lt
	asr	w3, w3, 3
	sbfiz	x3, x3, 3, 32
	sub	x3, x3, #64
	add	x1, x1, x3
	add	x0, x0, x3
	cmp	w2, #64
	bmi	.L9
	sub	w2, w2, #64
	and	x2, x2, 4294967232
	sub	x2, x1, x2
	sub	x2, x2, #64
.L11:
	ldr	x6, [x1,32]
	ldr	x5, [x1,40]
	ldr	x4, [x1,48]
	ldr	x3, [x1,56]
	str	x6, [x0,32]
	str	x5, [x0,40]
	str	x4, [x0,48]
	str	x3, [x0,56]
	ldr	x6, [x1]
	ldr	x5, [x1,8]
	ldr	x4, [x1,16]
	ldr	x3, [x1,24]
	str	x6, [x0]
	str	x5, [x0,8]
	str	x4, [x0,16]
	str	x3, [x0,24]
	sub	x1, x1, #64
	sub	x0, x0, #64
	cmp	x1, x2
	bne	.L11
.L9:
	ret
	.cfi_endproc
.LFE7:
	.size	_Z33aligned_block_copy_backwards_bs32PlS_i, .-_Z33aligned_block_copy_backwards_bs32PlS_i
	.align	2
	.global	_Z33aligned_block_copy_backwards_bs64PlS_i
	.type	_Z33aligned_block_copy_backwards_bs64PlS_i, %function
_Z33aligned_block_copy_backwards_bs64PlS_i:
.LFB8:
	.cfi_startproc
	add	w3, w2, 7
	cmp	w2, wzr
	csel	w3, w3, w2, lt
	asr	w3, w3, 3
	sbfiz	x3, x3, 3, 32
	sub	x3, x3, #64
	add	x1, x1, x3
	add	x0, x0, x3
	cmp	w2, #64
	bmi	.L13
	sub	w2, w2, #64
	and	x2, x2, 4294967232
	sub	x2, x1, x2
	sub	x2, x2, #64
.L15:
	ldr	x6, [x1]
	ldr	x5, [x1,8]
	ldr	x4, [x1,16]
	ldr	x3, [x1,24]
	str	x6, [x0]
	str	x5, [x0,8]
	str	x4, [x0,16]
	str	x3, [x0,24]
	ldr	x6, [x1,32]
	ldr	x5, [x1,40]
	ldr	x4, [x1,48]
	ldr	x3, [x1,56]
	str	x6, [x0,32]
	str	x5, [x0,40]
	str	x4, [x0,48]
	str	x3, [x0,56]
	sub	x1, x1, #64
	sub	x0, x0, #64
	cmp	x1, x2
	bne	.L15
.L13:
	ret
	.cfi_endproc
.LFE8:
	.size	_Z33aligned_block_copy_backwards_bs64PlS_i, .-_Z33aligned_block_copy_backwards_bs64PlS_i
	.align	2
	.global	_Z23aligned_block_copy_pf32PlS_i
	.type	_Z23aligned_block_copy_pf32PlS_i, %function
_Z23aligned_block_copy_pf32PlS_i:
.LFB9:
	.cfi_startproc
	cmp	w2, #64
	bmi	.L17
	add	x3, x1, 256
	sub	w6, w2, #64
	and	x5, x6, 4294967232
	add	x4, x5, 320
	add	x1, x1, x4
.L19:
	ldr	x6, [x3,-256]
	ldr	x5, [x3,-248]
	ldr	x4, [x3,-240]
	ldr	x2, [x3,-232]
	str	x6, [x0]
	str	x5, [x0,8]
	str	x4, [x0,16]
	str	x2, [x0,24]
	ldr	x6, [x3,-224]
	ldr	x5, [x3,-216]
	ldr	x4, [x3,-208]
	ldr	x2, [x3,-200]
	str	x6, [x0,32]
	str	x5, [x0,40]
	str	x4, [x0,48]
	str	x2, [x0,56]
	add	x3, x3, 64
	add	x0, x0, 64
	cmp	x3, x1
	bne	.L19
.L17:
	ret
	.cfi_endproc
.LFE9:
	.size	_Z23aligned_block_copy_pf32PlS_i, .-_Z23aligned_block_copy_pf32PlS_i
	.align	2
	.global	_Z23aligned_block_copy_pf64PlS_i
	.type	_Z23aligned_block_copy_pf64PlS_i, %function
_Z23aligned_block_copy_pf64PlS_i:
.LFB10:
	.cfi_startproc
	cmp	w2, #64
	bmi	.L21
	add	x3, x1, 256
	sub	w6, w2, #64
	and	x5, x6, 4294967232
	add	x4, x5, 320
	add	x1, x1, x4
.L23:
	ldr	x6, [x3,-256]
	ldr	x5, [x3,-248]
	ldr	x4, [x3,-240]
	ldr	x2, [x3,-232]
	str	x6, [x0]
	str	x5, [x0,8]
	str	x4, [x0,16]
	str	x2, [x0,24]
	ldr	x6, [x3,-224]
	ldr	x5, [x3,-216]
	ldr	x4, [x3,-208]
	ldr	x2, [x3,-200]
	str	x6, [x0,32]
	str	x5, [x0,40]
	str	x4, [x0,48]
	str	x2, [x0,56]
	add	x3, x3, 64
	add	x0, x0, 64
	cmp	x3, x1
	bne	.L23
.L21:
	ret
	.cfi_endproc
.LFE10:
	.size	_Z23aligned_block_copy_pf64PlS_i, .-_Z23aligned_block_copy_pf64PlS_i
	.align	2
	.global	_Z18aligned_block_fillPlS_i
	.type	_Z18aligned_block_fillPlS_i, %function
_Z18aligned_block_fillPlS_i:
.LFB11:
	.cfi_startproc
	ldr	x1, [x1]
	cmp	w2, #64
	bmi	.L25
	sub	w2, w2, #64
	and	x3, x2, 4294967232
	add	x3, x3, 64
	add	x3, x0, x3
.L27:
	str	x1, [x0]
	str	x1, [x0,8]
	str	x1, [x0,16]
	str	x1, [x0,24]
	str	x1, [x0,32]
	str	x1, [x0,40]
	str	x1, [x0,48]
	add	x2, x0, 64
	str	x1, [x0,56]
	mov	x0, x2
	cmp	x3, x2
	bne	.L27
.L25:
	ret
	.cfi_endproc
.LFE11:
	.size	_Z18aligned_block_fillPlS_i, .-_Z18aligned_block_fillPlS_i
	.align	2
	.global	_Z28aligned_block_fill_shuffle16PlS_i
	.type	_Z28aligned_block_fill_shuffle16PlS_i, %function
_Z28aligned_block_fill_shuffle16PlS_i:
.LFB12:
	.cfi_startproc
	ldr	x1, [x1]
	cmp	w2, #64
	bmi	.L29
	sub	w2, w2, #64
	and	x2, x2, 4294967232
	add	x2, x2, 64
	add	x2, x0, x2
.L31:
	str	x1, [x0]
	str	x1, [x0,8]
	str	x1, [x0,24]
	str	x1, [x0,16]
	str	x1, [x0,40]
	str	x1, [x0,32]
	str	x1, [x0,48]
	str	x1, [x0,56]
	add	x0, x0, 64
	cmp	x0, x2
	bne	.L31
.L29:
	ret
	.cfi_endproc
.LFE12:
	.size	_Z28aligned_block_fill_shuffle16PlS_i, .-_Z28aligned_block_fill_shuffle16PlS_i
	.align	2
	.global	_Z28aligned_block_fill_shuffle32PlS_i
	.type	_Z28aligned_block_fill_shuffle32PlS_i, %function
_Z28aligned_block_fill_shuffle32PlS_i:
.LFB13:
	.cfi_startproc
	ldr	x1, [x1]
	cmp	w2, #64
	bmi	.L33
	sub	w2, w2, #64
	and	x2, x2, 4294967232
	add	x2, x2, 64
	add	x2, x0, x2
.L35:
	str	x1, [x0,24]
	str	x1, [x0]
	str	x1, [x0,16]
	str	x1, [x0,8]
	str	x1, [x0,56]
	str	x1, [x0,32]
	str	x1, [x0,48]
	str	x1, [x0,40]
	add	x0, x0, 64
	cmp	x0, x2
	bne	.L35
.L33:
	ret
	.cfi_endproc
.LFE13:
	.size	_Z28aligned_block_fill_shuffle32PlS_i, .-_Z28aligned_block_fill_shuffle32PlS_i
	.align	2
	.global	_Z28aligned_block_fill_shuffle64PlS_i
	.type	_Z28aligned_block_fill_shuffle64PlS_i, %function
_Z28aligned_block_fill_shuffle64PlS_i:
.LFB14:
	.cfi_startproc
	ldr	x1, [x1]
	cmp	w2, #64
	bmi	.L37
	sub	w2, w2, #64
	and	x2, x2, 4294967232
	add	x2, x2, 64
	add	x2, x0, x2
.L39:
	str	x1, [x0,40]
	str	x1, [x0,16]
	str	x1, [x0,56]
	str	x1, [x0,48]
	str	x1, [x0,8]
	str	x1, [x0,24]
	str	x1, [x0]
	str	x1, [x0,32]
	add	x0, x0, 64
	cmp	x0, x2
	bne	.L39
.L37:
	ret
	.cfi_endproc
.LFE14:
	.size	_Z28aligned_block_fill_shuffle64PlS_i, .-_Z28aligned_block_fill_shuffle64PlS_i
	.align	2
	.global	_Z7gettimev
	.type	_Z7gettimev, %function
_Z7gettimev:
.LFB15:
	.cfi_startproc
	stp	x29, x30, [sp, -32]!
	.cfi_def_cfa_offset 32
	.cfi_offset 29, -32
	.cfi_offset 30, -24
	add	x29, sp, 0
	.cfi_def_cfa_register 29
	add	x0, x29, 16
	mov	x1, 0
	bl	gettimeofday
	ldr	x0, [x29,16]
	lsl	x1, x0, 8
	sub	x1, x1, x0, lsl 3
	add	x0, x0, x1, lsl 6
	sub	x0, x0, x1
	ldr	x1, [x29,24]
	add	x0, x1, x0, lsl 6
	scvtf	d0, x0
	ldr	d1, .LC0
	fdiv	d0, d0, d1
	ldp	x29, x30, [sp], 32
	.cfi_restore 30
	.cfi_restore 29
	.cfi_def_cfa 31, 0
	ret
	.cfi_endproc
.LFE15:
	.size	_Z7gettimev, .-_Z7gettimev
	.align	3
.LC0:
	.word	0
	.word	1093567616
	.align	2
	.global	_Z4fmindd
	.type	_Z4fmindd, %function
_Z4fmindd:
.LFB16:
	.cfi_startproc
	fmov	d2, d0
	fcmpe	d0, d1
	bmi	.L43
	fmov	d2, d1
.L43:
	fmov	d0, d2
	ret
	.cfi_endproc
.LFE16:
	.size	_Z4fmindd, .-_Z4fmindd
	.align	2
	.global	_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i
	.type	_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i, %function
_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i:
.LFB18:
	.cfi_startproc
	stp	x29, x30, [sp, -96]!
	.cfi_def_cfa_offset 96
	.cfi_offset 29, -96
	.cfi_offset 30, -88
	add	x29, sp, 0
	.cfi_def_cfa_register 29
	stp	x19, x20, [sp,16]
	stp	x21, x22, [sp,32]
	stp	x23, x24, [sp,48]
	stp	x25, x26, [sp,64]
	str	x27, [sp,80]
	.cfi_offset 19, -80
	.cfi_offset 20, -72
	.cfi_offset 21, -64
	.cfi_offset 22, -56
	.cfi_offset 23, -48
	.cfi_offset 24, -40
	.cfi_offset 25, -32
	.cfi_offset 26, -24
	.cfi_offset 27, -16
	mov	x26, x0
	mov	w21, w1
	mov	x25, x2
	mov	w20, w3
	mov	x23, x4
	mov	w19, w5
	mov	x22, x6
	cbz	x0, .L58
	tbz	w1, #31, .L46
.L58:
	mov	w21, 0
.L46:
	cbz	x25, .L59
	tbz	w20, #31, .L48
.L59:
	mov	w20, 0
.L48:
	cbz	x23, .L60
	tbz	w19, #31, .L50
.L60:
	mov	w19, 0
.L50:
	cbz	x22, .L61
	tbz	w7, #31, .L52
.L61:
	mov	w7, 0
.L52:
	add	w0, w21, w20
	add	w0, w0, w19
	add	w7, w0, w7
	add	w27, w7, 9437184
	sxtw	x27, w27
	mov	x0, x27
	bl	malloc
	mov	x24, x0
	mov	w1, 204
	mov	x2, x27
	bl	memset
	add	x0, x24, 1044480
	add	x0, x0, 4095
	and	x0, x0, -1048576
	cbz	x26, .L54
	add	x0, x0, 696320
	add	x0, x0, 2688
	str	x0, [x26]
	add	x0, x0, x21, sxtw
	add	x0, x0, 1044480
	add	x0, x0, 4095
	and	x0, x0, -1048576
.L54:
	cbz	x25, .L55
	add	x0, x0, 348160
	add	x0, x0, 1280
	str	x0, [x25]
	add	x0, x0, x20, sxtw
	add	x0, x0, 1044480
	add	x0, x0, 4095
	and	x0, x0, -1048576
.L55:
	cbz	x23, .L56
	add	x0, x0, 835584
	add	x0, x0, 3200
	str	x0, [x23]
	add	x0, x0, x19, sxtw
	add	x0, x0, 1044480
	add	x0, x0, 4095
	and	x0, x0, -1048576
.L56:
	cbz	x22, .L57
	add	x0, x0, 208896
	add	x0, x0, 768
	str	x0, [x22]
.L57:
	mov	x0, x24
	ldp	x19, x20, [sp,16]
	.cfi_restore 20
	.cfi_restore 19
	ldp	x21, x22, [sp,32]
	.cfi_restore 22
	.cfi_restore 21
	ldp	x23, x24, [sp,48]
	.cfi_restore 24
	.cfi_restore 23
	ldp	x25, x26, [sp,64]
	.cfi_restore 26
	.cfi_restore 25
	ldr	x27, [sp,80]
	.cfi_restore 27
	ldp	x29, x30, [sp], 96
	.cfi_restore 30
	.cfi_restore 29
	.cfi_def_cfa 31, 0
	ret
	.cfi_endproc
.LFE18:
	.size	_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i, .-_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i
	.ident	"GCC: (GNU) 4.9.x 20150123 (prerelease)"
	.section	.note.GNU-stack,"",%progbits
