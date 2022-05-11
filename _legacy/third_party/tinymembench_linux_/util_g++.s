	.cpu generic+fp+simd
	.file	"util.c"
	.text
	.align	2
	.global	_Z18aligned_block_copyPlS_i
	.type	_Z18aligned_block_copyPlS_i, %function
_Z18aligned_block_copyPlS_i:
.LFB5:
	.cfi_startproc
	cmp	w2, 63
	ble	.L1
	sub	w5, w2, #64
	and	x5, x5, 4294967232
	add	x5, x5, 64
	add	x5, x1, x5
.L3:
	ldr	x2, [x1]
	add	x1, x1, 64
	str	x2, [x0]
	add	x0, x0, 64
	cmp	x5, x1
	ldr	x2, [x1,-56]
	str	x2, [x0,-56]
	ldr	x2, [x1,-48]
	str	x2, [x0,-48]
	ldr	x2, [x1,-40]
	str	x2, [x0,-40]
	ldr	x4, [x1,-32]
	str	x4, [x0,-32]
	ldr	x3, [x1,-24]
	ldr	x2, [x1,-16]
	str	x3, [x0,-24]
	str	x2, [x0,-16]
	ldr	x2, [x1,-8]
	str	x2, [x0,-8]
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
	cmp	w2, wzr
	add	w3, w2, 7
	csel	w3, w3, w2, lt
	cmp	w2, 63
	asr	w3, w3, 3
	mov	x4, -8
	add	x3, x4, x3, sxtw 3
	add	x1, x1, x3
	add	x0, x0, x3
	ble	.L7
	sub	w2, w2, #64
	mov	x3, -64
	lsr	w4, w2, 6
	mov	w2, 64
	umsubl	x2, w4, w2, x3
	add	x2, x1, x2
.L9:
	ldr	x3, [x1]
	sub	x1, x1, #64
	str	x3, [x0]
	sub	x0, x0, #64
	cmp	x2, x1
	ldr	x3, [x1,56]
	str	x3, [x0,56]
	ldr	x3, [x1,48]
	str	x3, [x0,48]
	ldr	x3, [x1,40]
	str	x3, [x0,40]
	ldr	x5, [x1,32]
	str	x5, [x0,32]
	ldr	x4, [x1,24]
	ldr	x3, [x1,16]
	str	x4, [x0,24]
	str	x3, [x0,16]
	ldr	x3, [x1,8]
	str	x3, [x0,8]
	bne	.L9
.L7:
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
	cmp	w2, wzr
	add	w3, w2, 7
	csel	w3, w3, w2, lt
	cmp	w2, 63
	mov	x4, -64
	asr	w3, w3, 3
	add	x3, x4, x3, sxtw 3
	add	x1, x1, x3
	add	x0, x0, x3
	ble	.L12
	sub	w2, w2, #64
	mov	w3, 64
	lsr	w2, w2, 6
	umsubl	x2, w2, w3, x4
	add	x2, x1, x2
.L14:
	ldr	x3, [x1,32]
	sub	x1, x1, #64
	str	x3, [x0,32]
	sub	x0, x0, #64
	cmp	x1, x2
	ldr	x3, [x1,104]
	str	x3, [x0,104]
	ldr	x3, [x1,112]
	str	x3, [x0,112]
	ldr	x3, [x1,120]
	str	x3, [x0,120]
	ldr	x3, [x1,64]
	str	x3, [x0,64]
	ldr	x3, [x1,72]
	str	x3, [x0,72]
	ldr	x4, [x1,80]
	ldr	x3, [x1,88]
	str	x4, [x0,80]
	str	x3, [x0,88]
	bne	.L14
.L12:
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
	cmp	w2, wzr
	add	w3, w2, 7
	csel	w3, w3, w2, lt
	cmp	w2, 63
	mov	x4, -64
	asr	w3, w3, 3
	add	x3, x4, x3, sxtw 3
	add	x1, x1, x3
	add	x0, x0, x3
	ble	.L17
	sub	w2, w2, #64
	mov	w3, 64
	lsr	w2, w2, 6
	umsubl	x2, w2, w3, x4
	add	x2, x1, x2
.L19:
	ldr	x3, [x1]
	sub	x1, x1, #64
	str	x3, [x0]
	sub	x0, x0, #64
	cmp	x1, x2
	ldr	x3, [x1,72]
	str	x3, [x0,72]
	ldr	x3, [x1,80]
	str	x3, [x0,80]
	ldr	x3, [x1,88]
	str	x3, [x0,88]
	ldr	x3, [x1,96]
	str	x3, [x0,96]
	ldr	x3, [x1,104]
	str	x3, [x0,104]
	ldr	x4, [x1,112]
	ldr	x3, [x1,120]
	str	x4, [x0,112]
	str	x3, [x0,120]
	bne	.L19
.L17:
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
	cmp	w2, 63
	ble	.L22
	sub	w4, w2, #64
	add	x3, x1, 256
	and	x4, x4, 4294967232
	add	x4, x4, 320
	add	x1, x1, x4
.L24:
	add	x0, x0, 64
	ldr	x2, [x3,-256]
	ldr	x4, [x3,-208]
	str	x2, [x0,-64]
	ldr	x2, [x3,-248]
	str	x2, [x0,-56]
	ldr	x2, [x3,-240]
	str	x2, [x0,-48]
	ldr	x2, [x3,-232]
	str	x2, [x0,-40]
	ldr	x2, [x3,-224]
	str	x2, [x0,-32]
	ldr	x2, [x3,-216]
	str	x2, [x0,-24]
	str	x4, [x0,-16]
	ldr	x2, [x3,-200]
	add	x3, x3, 64
	cmp	x3, x1
	str	x2, [x0,-8]
	bne	.L24
.L22:
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
	cmp	w2, 63
	ble	.L27
	sub	w4, w2, #64
	add	x3, x1, 256
	and	x4, x4, 4294967232
	add	x4, x4, 320
	add	x1, x1, x4
.L29:
	add	x0, x0, 64
	ldr	x2, [x3,-256]
	ldr	x4, [x3,-208]
	str	x2, [x0,-64]
	ldr	x2, [x3,-248]
	str	x2, [x0,-56]
	ldr	x2, [x3,-240]
	str	x2, [x0,-48]
	ldr	x2, [x3,-232]
	str	x2, [x0,-40]
	ldr	x2, [x3,-224]
	str	x2, [x0,-32]
	ldr	x2, [x3,-216]
	str	x2, [x0,-24]
	str	x4, [x0,-16]
	ldr	x2, [x3,-200]
	add	x3, x3, 64
	cmp	x3, x1
	str	x2, [x0,-8]
	bne	.L29
.L27:
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
	cmp	w2, 63
	ldr	x1, [x1]
	ble	.L32
	sub	w2, w2, #64
	and	x2, x2, 4294967232
	add	x2, x2, 64
	add	x2, x0, x2
.L34:
	str	x1, [x0]
	add	x0, x0, 64
	cmp	x2, x0
	str	x1, [x0,-56]
	str	x1, [x0,-48]
	str	x1, [x0,-40]
	str	x1, [x0,-32]
	str	x1, [x0,-24]
	str	x1, [x0,-16]
	str	x1, [x0,-8]
	bne	.L34
.L32:
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
	cmp	w2, 63
	ldr	x1, [x1]
	ble	.L37
	sub	w2, w2, #64
	and	x2, x2, 4294967232
	add	x2, x2, 64
	add	x2, x0, x2
.L39:
	str	x1, [x0]
	add	x0, x0, 64
	cmp	x0, x2
	str	x1, [x0,-56]
	str	x1, [x0,-40]
	str	x1, [x0,-48]
	str	x1, [x0,-24]
	str	x1, [x0,-32]
	str	x1, [x0,-16]
	str	x1, [x0,-8]
	bne	.L39
.L37:
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
	cmp	w2, 63
	ldr	x1, [x1]
	ble	.L42
	sub	w2, w2, #64
	and	x2, x2, 4294967232
	add	x2, x2, 64
	add	x2, x0, x2
.L44:
	str	x1, [x0,24]
	add	x0, x0, 64
	cmp	x0, x2
	str	x1, [x0,-64]
	str	x1, [x0,-48]
	str	x1, [x0,-56]
	str	x1, [x0,-8]
	str	x1, [x0,-32]
	str	x1, [x0,-16]
	str	x1, [x0,-24]
	bne	.L44
.L42:
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
	cmp	w2, 63
	ldr	x1, [x1]
	ble	.L47
	sub	w2, w2, #64
	and	x2, x2, 4294967232
	add	x2, x2, 64
	add	x2, x0, x2
.L49:
	str	x1, [x0,40]
	add	x0, x0, 64
	cmp	x0, x2
	str	x1, [x0,-48]
	str	x1, [x0,-8]
	str	x1, [x0,-16]
	str	x1, [x0,-56]
	str	x1, [x0,-40]
	str	x1, [x0,-64]
	str	x1, [x0,-32]
	bne	.L49
.L47:
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
	mov	x1, 0
	add	x29, sp, 0
	.cfi_def_cfa_register 29
	add	x0, x29, 16
	bl	gettimeofday
	ldr	x0, [x29,16]
	ldr	d1, .LC0
	lsl	x1, x0, 8
	sub	x1, x1, x0, lsl 3
	add	x0, x0, x1, lsl 6
	sub	x0, x0, x1
	ldr	x1, [x29,24]
	ldp	x29, x30, [sp], 32
	.cfi_restore 30
	.cfi_restore 29
	.cfi_def_cfa 31, 0
	add	x0, x1, x0, lsl 6
	scvtf	d0, x0
	fdiv	d0, d0, d1
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
	fcmpe	d0, d1
	bmi	.L54
	fmov	d0, d1
.L54:
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
	str	x21, [sp,32]
	.cfi_offset 19, -80
	.cfi_offset 20, -72
	.cfi_offset 21, -64
	mov	w19, w1
	mov	x21, x0
	mov	x20, x2
	cbz	x0, .L69
	tbnz	w1, #31, .L69
.L57:
	cbz	x20, .L70
	tbnz	w3, #31, .L70
.L59:
	cbz	x4, .L71
	tbnz	w5, #31, .L71
	cbz	x6, .L72
.L97:
	tbnz	w7, #31, .L72
.L63:
	add	w0, w19, w3
	str	x6, [x29,56]
	add	w0, w0, w5
	str	x4, [x29,64]
	add	w7, w0, w7
	str	x3, [x29,72]
	add	w2, w7, 9437184
	str	x5, [x29,80]
	sxtw	x2, w2
	str	x2, [x29,88]
	mov	x0, x2
	bl	malloc
	ldr	x2, [x29,88]
	mov	w1, 204
	bl	memset
	mov	x7, x0
	add	x0, x0, 1044480
	ldr	x5, [x29,80]
	add	x0, x0, 4095
	ldr	x3, [x29,72]
	and	x0, x0, -1048576
	ldr	x4, [x29,64]
	ldr	x6, [x29,56]
	cbz	x21, .L65
	add	x0, x0, 696320
	add	x0, x0, 2688
	add	x19, x0, x19, sxtw
	add	x19, x19, 1044480
	str	x0, [x21]
	add	x19, x19, 4095
	and	x0, x19, -1048576
.L65:
	cbz	x20, .L66
	add	x0, x0, 348160
	add	x0, x0, 1280
	add	x3, x0, x3, sxtw
	add	x3, x3, 1044480
	str	x0, [x20]
	add	x3, x3, 4095
	and	x0, x3, -1048576
.L66:
	cbz	x4, .L67
	add	x0, x0, 835584
	add	x0, x0, 3200
	add	x5, x0, x5, sxtw
	add	x5, x5, 1044480
	str	x0, [x4]
	add	x5, x5, 4095
	and	x0, x5, -1048576
.L67:
	cbz	x6, .L68
	add	x0, x0, 208896
	add	x0, x0, 768
	str	x0, [x6]
.L68:
	ldp	x19, x20, [sp,16]
	.cfi_remember_state
	.cfi_restore 20
	.cfi_restore 19
	ldr	x21, [sp,32]
	.cfi_restore 21
	mov	x0, x7
	ldp	x29, x30, [sp], 96
	.cfi_restore 30
	.cfi_restore 29
	.cfi_def_cfa 31, 0
	ret
.L71:
	.cfi_restore_state
	mov	w5, 0
	cbnz	x6, .L97
.L72:
	mov	w7, 0
	b	.L63
.L70:
	mov	w3, 0
	b	.L59
.L69:
	mov	w19, 0
	b	.L57
	.cfi_endproc
.LFE18:
	.size	_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i, .-_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i
	.ident	"GCC: (GNU) 4.9.x 20150123 (prerelease)"
	.section	.note.GNU-stack,"",%progbits
