	.text
	.file	"util.c"
	.globl	_Z18aligned_block_copyPlS_i // -- Begin function _Z18aligned_block_copyPlS_i
	.p2align	2
	.type	_Z18aligned_block_copyPlS_i,@function
_Z18aligned_block_copyPlS_i:            // @_Z18aligned_block_copyPlS_i
// BB#0:
	cmp		w2, #64         // =64
	b.lt	.LBB0_3
// BB#1:
	add	w8, w2, #64             // =64
.LBB0_2:                                // =>This Inner Loop Header: Depth=1
	ldp		x9, x10, [x1]
	ldp	x11, x12, [x1, #16]
	sub	w8, w8, #64             // =64
	cmp		w8, #127        // =127
	str		x9, [x0]
	str	x10, [x0, #8]
	str	x11, [x0, #16]
	str	x12, [x0, #24]
	ldp	x9, x10, [x1, #32]
	ldp	x11, x12, [x1, #48]
	add	x1, x1, #64             // =64
	str	x9, [x0, #32]
	str	x10, [x0, #40]
	str	x11, [x0, #48]
	str	x12, [x0, #56]
	add	x0, x0, #64             // =64
	b.gt	.LBB0_2
.LBB0_3:
	ret
.Lfunc_end0:
	.size	_Z18aligned_block_copyPlS_i, .Lfunc_end0-_Z18aligned_block_copyPlS_i
                                        // -- End function
	.globl	_Z28aligned_block_copy_backwardsPlS_i // -- Begin function _Z28aligned_block_copy_backwardsPlS_i
	.p2align	2
	.type	_Z28aligned_block_copy_backwardsPlS_i,@function
_Z28aligned_block_copy_backwardsPlS_i:  // @_Z28aligned_block_copy_backwardsPlS_i
// BB#0:
	add	w8, w2, #7              // =7
	cmp		w2, #0          // =0
	csel	w8, w8, w2, lt
	cmp		w2, #64         // =64
	b.lt	.LBB1_3
// BB#1:
	asr	w8, w8, #3
	sub	w8, w8, #1              // =1
	sxtw	x8, w8
	lsl	x9, x8, #3
	add		x8, x0, x9
	add		x9, x1, x9
	add	w10, w2, #64            // =64
.LBB1_2:                                // =>This Inner Loop Header: Depth=1
	ldp	x12, x11, [x9, #-8]
	ldp	x13, x14, [x9, #-24]
	sub	w10, w10, #64           // =64
	cmp		w10, #127       // =127
	str		x11, [x8]
	stur	x12, [x8, #-8]
	stur	x14, [x8, #-16]
	stur	x13, [x8, #-24]
	ldp	x11, x12, [x9, #-40]
	ldp	x13, x14, [x9, #-56]
	sub	x9, x9, #64             // =64
	stur	x12, [x8, #-32]
	stur	x11, [x8, #-40]
	stur	x14, [x8, #-48]
	stur	x13, [x8, #-56]
	sub	x8, x8, #64             // =64
	b.gt	.LBB1_2
.LBB1_3:
	ret
.Lfunc_end1:
	.size	_Z28aligned_block_copy_backwardsPlS_i, .Lfunc_end1-_Z28aligned_block_copy_backwardsPlS_i
                                        // -- End function
	.globl	_Z33aligned_block_copy_backwards_bs32PlS_i // -- Begin function _Z33aligned_block_copy_backwards_bs32PlS_i
	.p2align	2
	.type	_Z33aligned_block_copy_backwards_bs32PlS_i,@function
_Z33aligned_block_copy_backwards_bs32PlS_i: // @_Z33aligned_block_copy_backwards_bs32PlS_i
// BB#0:
	add	w8, w2, #7              // =7
	cmp		w2, #0          // =0
	csel	w8, w8, w2, lt
	cmp		w2, #64         // =64
	b.lt	.LBB2_3
// BB#1:
	asr	w9, w8, #3
	sub	w9, w9, #8              // =8
	sxtw	x9, w9
	lsl	x9, x9, #3
	add	x10, x9, #32            // =32
	add	w8, w2, #64             // =64
	add		x9, x0, x10
	add		x10, x1, x10
.LBB2_2:                                // =>This Inner Loop Header: Depth=1
	ldp		x11, x12, [x10]
	ldp	x13, x14, [x10, #16]
	sub	w8, w8, #64             // =64
	cmp		w8, #127        // =127
	str		x11, [x9]
	str	x12, [x9, #8]
	str	x13, [x9, #16]
	str	x14, [x9, #24]
	ldp	x11, x12, [x10, #-32]
	ldp	x13, x14, [x10, #-16]
	sub	x10, x10, #64           // =64
	stur	x11, [x9, #-32]
	stur	x12, [x9, #-24]
	stur	x13, [x9, #-16]
	stur	x14, [x9, #-8]
	sub	x9, x9, #64             // =64
	b.gt	.LBB2_2
.LBB2_3:
	ret
.Lfunc_end2:
	.size	_Z33aligned_block_copy_backwards_bs32PlS_i, .Lfunc_end2-_Z33aligned_block_copy_backwards_bs32PlS_i
                                        // -- End function
	.globl	_Z33aligned_block_copy_backwards_bs64PlS_i // -- Begin function _Z33aligned_block_copy_backwards_bs64PlS_i
	.p2align	2
	.type	_Z33aligned_block_copy_backwards_bs64PlS_i,@function
_Z33aligned_block_copy_backwards_bs64PlS_i: // @_Z33aligned_block_copy_backwards_bs64PlS_i
// BB#0:
	add	w8, w2, #7              // =7
	cmp		w2, #0          // =0
	csel	w8, w8, w2, lt
	cmp		w2, #64         // =64
	b.lt	.LBB3_3
// BB#1:
	asr	w8, w8, #3
	sub	w8, w8, #8              // =8
	sxtw	x8, w8
	lsl	x9, x8, #3
	add		x8, x0, x9
	add		x9, x1, x9
	add	w10, w2, #64            // =64
.LBB3_2:                                // =>This Inner Loop Header: Depth=1
	ldp		x11, x12, [x9]
	ldp	x13, x14, [x9, #16]
	sub	w10, w10, #64           // =64
	cmp		w10, #127       // =127
	str		x11, [x8]
	str	x12, [x8, #8]
	str	x13, [x8, #16]
	str	x14, [x8, #24]
	ldp	x11, x12, [x9, #32]
	ldp	x13, x14, [x9, #48]
	sub	x9, x9, #64             // =64
	str	x11, [x8, #32]
	str	x12, [x8, #40]
	str	x13, [x8, #48]
	str	x14, [x8, #56]
	sub	x8, x8, #64             // =64
	b.gt	.LBB3_2
.LBB3_3:
	ret
.Lfunc_end3:
	.size	_Z33aligned_block_copy_backwards_bs64PlS_i, .Lfunc_end3-_Z33aligned_block_copy_backwards_bs64PlS_i
                                        // -- End function
	.globl	_Z23aligned_block_copy_pf32PlS_i // -- Begin function _Z23aligned_block_copy_pf32PlS_i
	.p2align	2
	.type	_Z23aligned_block_copy_pf32PlS_i,@function
_Z23aligned_block_copy_pf32PlS_i:       // @_Z23aligned_block_copy_pf32PlS_i
// BB#0:
	cmp		w2, #64         // =64
	b.lt	.LBB4_3
// BB#1:
	add	w8, w2, #64             // =64
	add	x9, x1, #256            // =256
.LBB4_2:                                // =>This Inner Loop Header: Depth=1
	prfm	 pldl1strm, [x9]
	ldp	x10, x11, [x9, #-256]
	ldp	x12, x13, [x9, #-240]
	sub	w8, w8, #64             // =64
	cmp		w8, #127        // =127
	str		x10, [x0]
	str	x11, [x0, #8]
	str	x12, [x0, #16]
	str	x13, [x0, #24]
	prfm	pldl1strm, [x9, #32]
	ldp	x10, x11, [x9, #-224]
	ldp	x12, x13, [x9, #-208]
	add	x9, x9, #64             // =64
	str	x10, [x0, #32]
	str	x11, [x0, #40]
	str	x12, [x0, #48]
	str	x13, [x0, #56]
	add	x0, x0, #64             // =64
	b.gt	.LBB4_2
.LBB4_3:
	ret
.Lfunc_end4:
	.size	_Z23aligned_block_copy_pf32PlS_i, .Lfunc_end4-_Z23aligned_block_copy_pf32PlS_i
                                        // -- End function
	.globl	_Z23aligned_block_copy_pf64PlS_i // -- Begin function _Z23aligned_block_copy_pf64PlS_i
	.p2align	2
	.type	_Z23aligned_block_copy_pf64PlS_i,@function
_Z23aligned_block_copy_pf64PlS_i:       // @_Z23aligned_block_copy_pf64PlS_i
// BB#0:
	cmp		w2, #64         // =64
	b.lt	.LBB5_3
// BB#1:
	add	x8, x1, #256            // =256
	add	w9, w2, #64             // =64
.LBB5_2:                                // =>This Inner Loop Header: Depth=1
	prfm	 pldl1strm, [x8]
	ldp	x10, x11, [x8, #-256]
	ldp	x12, x13, [x8, #-240]
	sub	w9, w9, #64             // =64
	cmp		w9, #127        // =127
	str		x10, [x0]
	str	x11, [x0, #8]
	str	x12, [x0, #16]
	str	x13, [x0, #24]
	ldp	x10, x11, [x8, #-224]
	ldp	x12, x13, [x8, #-208]
	add	x8, x8, #64             // =64
	str	x10, [x0, #32]
	str	x11, [x0, #40]
	str	x12, [x0, #48]
	str	x13, [x0, #56]
	add	x0, x0, #64             // =64
	b.gt	.LBB5_2
.LBB5_3:
	ret
.Lfunc_end5:
	.size	_Z23aligned_block_copy_pf64PlS_i, .Lfunc_end5-_Z23aligned_block_copy_pf64PlS_i
                                        // -- End function
	.globl	_Z18aligned_block_fillPlS_i // -- Begin function _Z18aligned_block_fillPlS_i
	.p2align	2
	.type	_Z18aligned_block_fillPlS_i,@function
_Z18aligned_block_fillPlS_i:            // @_Z18aligned_block_fillPlS_i
// BB#0:
	cmp		w2, #64         // =64
	b.lt	.LBB6_3
// BB#1:
	ldr		x8, [x1]
	add	w9, w2, #64             // =64
.LBB6_2:                                // =>This Inner Loop Header: Depth=1
	sub	w9, w9, #64             // =64
	str		x8, [x0]
	str	x8, [x0, #8]
	str	x8, [x0, #16]
	str	x8, [x0, #24]
	str	x8, [x0, #32]
	str	x8, [x0, #40]
	str	x8, [x0, #48]
	str	x8, [x0, #56]
	cmp		w9, #127        // =127
	add	x0, x0, #64             // =64
	b.gt	.LBB6_2
.LBB6_3:
	ret
.Lfunc_end6:
	.size	_Z18aligned_block_fillPlS_i, .Lfunc_end6-_Z18aligned_block_fillPlS_i
                                        // -- End function
	.globl	_Z28aligned_block_fill_shuffle16PlS_i // -- Begin function _Z28aligned_block_fill_shuffle16PlS_i
	.p2align	2
	.type	_Z28aligned_block_fill_shuffle16PlS_i,@function
_Z28aligned_block_fill_shuffle16PlS_i:  // @_Z28aligned_block_fill_shuffle16PlS_i
// BB#0:
	cmp		w2, #64         // =64
	b.lt	.LBB7_3
// BB#1:
	ldr		x8, [x1]
	add	w9, w2, #64             // =64
.LBB7_2:                                // =>This Inner Loop Header: Depth=1
	sub	w9, w9, #64             // =64
	str		x8, [x0]
	str	x8, [x0, #8]
	str	x8, [x0, #24]
	str	x8, [x0, #16]
	str	x8, [x0, #40]
	str	x8, [x0, #32]
	str	x8, [x0, #48]
	str	x8, [x0, #56]
	cmp		w9, #127        // =127
	add	x0, x0, #64             // =64
	b.gt	.LBB7_2
.LBB7_3:
	ret
.Lfunc_end7:
	.size	_Z28aligned_block_fill_shuffle16PlS_i, .Lfunc_end7-_Z28aligned_block_fill_shuffle16PlS_i
                                        // -- End function
	.globl	_Z28aligned_block_fill_shuffle32PlS_i // -- Begin function _Z28aligned_block_fill_shuffle32PlS_i
	.p2align	2
	.type	_Z28aligned_block_fill_shuffle32PlS_i,@function
_Z28aligned_block_fill_shuffle32PlS_i:  // @_Z28aligned_block_fill_shuffle32PlS_i
// BB#0:
	cmp		w2, #64         // =64
	b.lt	.LBB8_3
// BB#1:
	ldr		x8, [x1]
	add	w9, w2, #64             // =64
	add	x10, x0, #32            // =32
.LBB8_2:                                // =>This Inner Loop Header: Depth=1
	sub	w9, w9, #64             // =64
	stur	x8, [x10, #-8]
	stur	x8, [x10, #-32]
	stur	x8, [x10, #-16]
	stur	x8, [x10, #-24]
	str	x8, [x10, #24]
	str		x8, [x10]
	str	x8, [x10, #16]
	str	x8, [x10, #8]
	cmp		w9, #127        // =127
	add	x10, x10, #64           // =64
	b.gt	.LBB8_2
.LBB8_3:
	ret
.Lfunc_end8:
	.size	_Z28aligned_block_fill_shuffle32PlS_i, .Lfunc_end8-_Z28aligned_block_fill_shuffle32PlS_i
                                        // -- End function
	.globl	_Z28aligned_block_fill_shuffle64PlS_i // -- Begin function _Z28aligned_block_fill_shuffle64PlS_i
	.p2align	2
	.type	_Z28aligned_block_fill_shuffle64PlS_i,@function
_Z28aligned_block_fill_shuffle64PlS_i:  // @_Z28aligned_block_fill_shuffle64PlS_i
// BB#0:
	cmp		w2, #64         // =64
	b.lt	.LBB9_3
// BB#1:
	ldr		x8, [x1]
	add	w9, w2, #64             // =64
	add	x10, x0, #32            // =32
.LBB9_2:                                // =>This Inner Loop Header: Depth=1
	sub	w9, w9, #64             // =64
	str	x8, [x10, #8]
	stur	x8, [x10, #-16]
	str	x8, [x10, #24]
	str	x8, [x10, #16]
	stur	x8, [x10, #-24]
	stur	x8, [x10, #-8]
	stur	x8, [x10, #-32]
	str	x8, [x10], #64
	cmp		w9, #127        // =127
	b.gt	.LBB9_2
.LBB9_3:
	ret
.Lfunc_end9:
	.size	_Z28aligned_block_fill_shuffle64PlS_i, .Lfunc_end9-_Z28aligned_block_fill_shuffle64PlS_i
                                        // -- End function
	.section	.rodata.cst8,"aM",@progbits,8
	.p2align	3               // -- Begin function _Z7gettimev
.LCPI10_0:
	.xword	4696837146684686336     // double 1.0E+6
	.text
	.globl	_Z7gettimev
	.p2align	2
	.type	_Z7gettimev,@function
_Z7gettimev:                            // @_Z7gettimev
// BB#0:
	sub	sp, sp, #32             // =32
	mov	 x0, sp
	mov	 x1, xzr
	stp	x29, x30, [sp, #16]     // 8-byte Folded Spill
	add	x29, sp, #16            // =16
	bl	gettimeofday
	ldp		x8, x9, [sp]
	adrp	x10, .LCPI10_0
	ldr	d0, [x10, :lo12:.LCPI10_0]
	mov	w10, #16960
	movk	w10, #15, lsl #16
	ldp	x29, x30, [sp, #16]     // 8-byte Folded Reload
	nop
	madd	x8, x8, x10, x9
	scvtf	d1, x8
	fdiv	d0, d1, d0
	add	sp, sp, #32             // =32
	ret
.Lfunc_end10:
	.size	_Z7gettimev, .Lfunc_end10-_Z7gettimev
                                        // -- End function
	.globl	_Z4fmindd               // -- Begin function _Z4fmindd
	.p2align	2
	.type	_Z4fmindd,@function
_Z4fmindd:                              // @_Z4fmindd
// BB#0:
	fcmp	d0, d1
	fcsel	d0, d0, d1, mi
	ret
.Lfunc_end11:
	.size	_Z4fmindd, .Lfunc_end11-_Z4fmindd
                                        // -- End function
	.globl	_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i // -- Begin function _Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i
	.p2align	2
	.type	_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i,@function
_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i: // @_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i
// BB#0:
	str	x27, [sp, #-96]!        // 8-byte Folded Spill
	stp	x24, x23, [sp, #32]     // 8-byte Folded Spill
	mov	 x23, x0
	cmp		w1, #0          // =0
	ccmp	x23, #0, #4, ge
	stp	x22, x21, [sp, #48]     // 8-byte Folded Spill
	mov	 x22, x2
	csel	w27, wzr, w1, eq
	cmp		w3, #0          // =0
	ccmp	x22, #0, #4, ge
	stp	x26, x25, [sp, #16]     // 8-byte Folded Spill
	stp	x20, x19, [sp, #64]     // 8-byte Folded Spill
	mov	 x20, x4
	csel	w26, wzr, w3, eq
	cmp		w5, #0          // =0
	ccmp	x20, #0, #4, ge
	mov	 x19, x6
	csel	w25, wzr, w5, eq
	cmp		w7, #0          // =0
	add		w8, w27, w26
	ccmp	x19, #0, #4, ge
	csel	w9, wzr, w7, eq
	add		w8, w8, w25
	add		w8, w8, w9
	add	w8, w8, #2304, lsl #12  // =9437184
	sxtw	x24, w8
	mov	 x0, x24
	stp	x29, x30, [sp, #80]     // 8-byte Folded Spill
	add	x29, sp, #80            // =80
	bl	malloc
	mov	w1, #204
	mov	 x2, x24
	mov	 x21, x0
	bl	memset
	orr	w8, wzr, #0xfffff
	add		x8, x21, x8
	and	x8, x8, #0xfffffffffff00000
	cbz	x23, .LBB12_2
// BB#1:
	mov	w9, #43648
	movk	w9, #10, lsl #16
	add		x8, x8, x9
	str		x8, [x23]
	add	x8, x8, w27, sxtw
	orr	w9, wzr, #0xfffff
	add		x8, x8, x9
	and	x8, x8, #0xfffffffffff00000
.LBB12_2:
	cbz	x22, .LBB12_4
// BB#3:
	mov	w9, #21760
	movk	w9, #5, lsl #16
	add		x8, x8, x9
	str		x8, [x22]
	add	x8, x8, w26, sxtw
	orr	w9, wzr, #0xfffff
	add		x8, x8, x9
	and	x8, x8, #0xfffffffffff00000
.LBB12_4:
	cbz	x20, .LBB12_6
// BB#5:
	mov	w9, #52352
	movk	w9, #12, lsl #16
	add		x8, x8, x9
	str		x8, [x20]
	add	x8, x8, w25, sxtw
	orr	w9, wzr, #0xfffff
	add		x8, x8, x9
	and	x8, x8, #0xfffffffffff00000
.LBB12_6:
	cbz	x19, .LBB12_8
// BB#7:
	mov	w9, #13056
	movk	w9, #3, lsl #16
	add		x8, x8, x9
	str		x8, [x19]
.LBB12_8:
	mov	 x0, x21
	ldp	x29, x30, [sp, #80]     // 8-byte Folded Reload
	ldp	x20, x19, [sp, #64]     // 8-byte Folded Reload
	ldp	x22, x21, [sp, #48]     // 8-byte Folded Reload
	ldp	x24, x23, [sp, #32]     // 8-byte Folded Reload
	ldp	x26, x25, [sp, #16]     // 8-byte Folded Reload
	ldr	x27, [sp], #96          // 8-byte Folded Reload
	ret
.Lfunc_end12:
	.size	_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i, .Lfunc_end12-_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i
                                        // -- End function

	.ident	"Android (4691093 based on r316199) clang version 6.0.2 (https://android.googlesource.com/toolchain/clang 183abd29fc496f55536e7d904e0abae47888fc7f) (https://android.googlesource.com/toolchain/llvm 34361f192e41ed6e4e8f9aca80a4ea7e9856f327) (based on LLVM 6.0.2svn)"
	.section	".note.GNU-stack","",@progbits
