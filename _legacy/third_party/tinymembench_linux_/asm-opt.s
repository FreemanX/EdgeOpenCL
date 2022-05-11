	.cpu generic+fp+simd
	.file	"asm-opt.c"
	.text
	.align	2
	.global	_Z17check_cpu_featurePKc
	.type	_Z17check_cpu_featurePKc, %function
_Z17check_cpu_featurePKc:
.LFB7:
	.cfi_startproc
	stp	x29, x30, [sp, -112]!
	.cfi_def_cfa_offset 112
	.cfi_offset 29, -112
	.cfi_offset 30, -104
	add	x29, sp, 0
	.cfi_def_cfa_register 29
	stp	x19, x20, [sp,16]
	stp	x21, x22, [sp,32]
	stp	x23, x24, [sp,48]
	stp	x25, x26, [sp,64]
	stp	x27, x28, [sp,80]
	.cfi_offset 19, -96
	.cfi_offset 20, -88
	.cfi_offset 21, -80
	.cfi_offset 22, -72
	.cfi_offset 23, -64
	.cfi_offset 24, -56
	.cfi_offset 25, -48
	.cfi_offset 26, -40
	.cfi_offset 27, -32
	.cfi_offset 28, -24
	mov	x22, x0
	mov	w24, 1024
	adrp	x0, .LC0
	add	x0, x0, :lo12:.LC0
	str	x0, [x29,104]
	adrp	x0, .LC1
	add	x0, x0, :lo12:.LC1
	str	x0, [x29,96]
	mov	w26, 10
.L7:
	sxtw	x0, w24
	bl	malloc
	mov	x23, x0
	cbz	x0, .L16
	ldr	x0, [x29,104]
	ldr	x1, [x29,96]
	bl	fopen
	mov	x25, x0
	cbnz	x0, .L17
	mov	w27, 0
	b	.L4
.L15:
	mov	x0, x23
	mov	w1, w26
	bl	strchr
	cbnz	x0, .L5
	mov	x0, x25
	bl	feof
	mov	w19, w0
	cbnz	w0, .L5
	mov	x0, x25
	bl	fclose
	mov	x0, x23
	bl	free
	lsl	w24, w24, 1
	cmp	w24, 1048576
	ble	.L7
	b	.L6
.L5:
	ldrb	w0, [x22]
	cbz	w0, .L3
	ldrb	w0, [x23]
	add	x19, x23, 1
	cmp	w0, 63
	beq	.L9
	b	.L3
.L10:
	add	x19, x19, 1
.L9:
	ldrb	w0, [x19]
	bl	isspace
	cbnz	w0, .L10
	add	x19, x19, 1
	b	.L11
.L14:
	cmp	x20, x21
	bls	.L12
	ldrb	w0, [x20,-1]
	bl	isspace
	cbz	w0, .L13
.L12:
	mov	x0, x22
	bl	strlen
	ldrb	w0, [x20,x0]
	cbz	w0, .L18
	bl	isspace
	cbnz	w0, .L19
.L13:
	add	x19, x19, 1
.L11:
	sub	x21, x19, #1
	mov	x0, x21
	mov	x1, x22
	bl	strstr
	mov	x20, x0
	cbnz	x0, .L14
	b	.L3
.L18:
	mov	w27, w28
	b	.L3
.L19:
	mov	w27, w28
	b	.L3
.L17:
	mov	w27, 0
	mov	w28, 1
.L3:
	mov	x0, x23
	mov	w1, w24
	mov	x2, x25
	bl	fgets
	cbnz	x0, .L15
	mov	x0, x25
	bl	fclose
.L4:
	mov	x0, x23
	bl	free
	b	.L2
.L16:
	mov	w27, 0
.L2:
	mov	w19, w27
.L6:
	mov	w0, w19
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
	ldp	x27, x28, [sp,80]
	.cfi_restore 28
	.cfi_restore 27
	ldp	x29, x30, [sp], 112
	.cfi_restore 30
	.cfi_restore 29
	.cfi_def_cfa 31, 0
	ret
	.cfi_endproc
.LFE7:
	.size	_Z17check_cpu_featurePKc, .-_Z17check_cpu_featurePKc
	.align	2
	.global	_Z18get_asm_benchmarksv
	.type	_Z18get_asm_benchmarksv, %function
_Z18get_asm_benchmarksv:
.LFB8:
	.cfi_startproc
	adrp	x0, .LANCHOR0
	add	x0, x0, :lo12:.LANCHOR0
	ret
	.cfi_endproc
.LFE8:
	.size	_Z18get_asm_benchmarksv, .-_Z18get_asm_benchmarksv
	.align	2
	.global	_Z30get_asm_framebuffer_benchmarksv
	.type	_Z30get_asm_framebuffer_benchmarksv, %function
_Z30get_asm_framebuffer_benchmarksv:
.LFB9:
	.cfi_startproc
	adrp	x0, .LANCHOR0
	add	x0, x0, :lo12:.LANCHOR0
	add	x0, x0, 432
	ret
	.cfi_endproc
.LFE9:
	.size	_Z30get_asm_framebuffer_benchmarksv, .-_Z30get_asm_framebuffer_benchmarksv
	.section	.rodata.str1.8,"aMS",%progbits,1
	.align	3
.LC0:
	.string	"/proc/cpuinfo"
	.zero	2
.LC1:
	.string	"r"
	.zero	6
.LC2:
	.string	"NEON LDP/STP copy (from framebuffer)"
	.zero	3
.LC3:
	.string	"NEON LDP/STP 2-pass copy (from framebuffer)"
	.zero	4
.LC4:
	.string	"NEON LD1/ST1 copy (from framebuffer)"
	.zero	3
.LC5:
	.string	"NEON LD1/ST1 2-pass copy (from framebuffer)"
	.zero	4
.LC6:
	.string	"ARM LDP/STP copy (from framebuffer)"
	.zero	4
.LC7:
	.string	"ARM LDP/STP 2-pass copy (from framebuffer)"
	.zero	5
.LC8:
	.string	"NEON LDP/STP copy"
	.zero	6
.LC9:
	.string	"NEON LDP/STP copy pldl2strm (32 bytes step)"
	.zero	4
.LC10:
	.string	"NEON LDP/STP copy pldl2strm (64 bytes step)"
	.zero	4
.LC11:
	.string	"NEON LDP/STP copy pldl1keep (32 bytes step)"
	.zero	4
.LC12:
	.string	"NEON LDP/STP copy pldl1keep (64 bytes step)"
	.zero	4
.LC13:
	.string	"NEON LD1/ST1 copy"
	.zero	6
.LC14:
	.string	"NEON STP fill"
	.zero	2
.LC15:
	.string	"NEON STNP fill"
	.zero	1
.LC16:
	.string	"ARM LDP/STP copy"
	.zero	7
.LC17:
	.string	"ARM STP fill"
	.zero	3
.LC18:
	.string	"ARM STNP fill"
	.section	.data.rel,"aw",%progbits
	.align	3
.LANCHOR0 = . + 0
	.type	_ZL12aarch64_neon, %object
	.size	_ZL12aarch64_neon, 432
_ZL12aarch64_neon:
	.xword	.LC8
	.word	0
	.zero	4
	.xword	aligned_block_copy_ldpstp_q_aarch64
	.xword	.LC9
	.word	0
	.zero	4
	.xword	aligned_block_copy_ldpstp_q_pf32_l2strm_aarch64
	.xword	.LC10
	.word	0
	.zero	4
	.xword	aligned_block_copy_ldpstp_q_pf64_l2strm_aarch64
	.xword	.LC11
	.word	0
	.zero	4
	.xword	aligned_block_copy_ldpstp_q_pf32_l1keep_aarch64
	.xword	.LC12
	.word	0
	.zero	4
	.xword	aligned_block_copy_ldpstp_q_pf64_l1keep_aarch64
	.xword	.LC13
	.word	0
	.zero	4
	.xword	aligned_block_copy_ld1st1_aarch64
	.xword	.LC14
	.word	0
	.zero	4
	.xword	aligned_block_fill_stp_q_aarch64
	.xword	.LC15
	.word	0
	.zero	4
	.xword	aligned_block_fill_stnp_q_aarch64
	.xword	.LC16
	.word	0
	.zero	4
	.xword	aligned_block_copy_ldpstp_x_aarch64
	.xword	.LC17
	.word	0
	.zero	4
	.xword	aligned_block_fill_stp_x_aarch64
	.xword	.LC18
	.word	0
	.zero	4
	.xword	aligned_block_fill_stnp_x_aarch64
	.xword	.LC2
	.word	0
	.zero	4
	.xword	aligned_block_copy_ldpstp_q_aarch64
	.xword	.LC3
	.word	1
	.zero	4
	.xword	aligned_block_copy_ldpstp_q_aarch64
	.xword	.LC4
	.word	0
	.zero	4
	.xword	aligned_block_copy_ld1st1_aarch64
	.xword	.LC5
	.word	1
	.zero	4
	.xword	aligned_block_copy_ld1st1_aarch64
	.xword	.LC6
	.word	0
	.zero	4
	.xword	aligned_block_copy_ldpstp_x_aarch64
	.xword	.LC7
	.word	1
	.zero	4
	.xword	aligned_block_copy_ldpstp_x_aarch64
	.xword	0
	.word	0
	.zero	4
	.xword	0
	.type	_ZL15aarch64_neon_fb, %object
	.size	_ZL15aarch64_neon_fb, 168
_ZL15aarch64_neon_fb:
	.xword	.LC2
	.word	0
	.zero	4
	.xword	aligned_block_copy_ldpstp_q_aarch64
	.xword	.LC3
	.word	1
	.zero	4
	.xword	aligned_block_copy_ldpstp_q_aarch64
	.xword	.LC4
	.word	0
	.zero	4
	.xword	aligned_block_copy_ld1st1_aarch64
	.xword	.LC5
	.word	1
	.zero	4
	.xword	aligned_block_copy_ld1st1_aarch64
	.xword	.LC6
	.word	0
	.zero	4
	.xword	aligned_block_copy_ldpstp_x_aarch64
	.xword	.LC7
	.word	1
	.zero	4
	.xword	aligned_block_copy_ldpstp_x_aarch64
	.xword	0
	.word	0
	.zero	4
	.xword	0
	.ident	"GCC: (GNU) 4.9.x 20150123 (prerelease)"
	.section	.note.GNU-stack,"",%progbits
