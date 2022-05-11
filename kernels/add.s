	.section	__TEXT,__text,regular,pure_instructions
	.build_version macos, 10, 15	sdk_version 10, 15, 6
	.globl	_vecAdd                 ## -- Begin function vecAdd
	.p2align	4, 0x90
_vecAdd:                                ## @vecAdd
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%r15
	pushq	%r14
	pushq	%r12
	pushq	%rbx
	.cfi_offset %rbx, -48
	.cfi_offset %r12, -40
	.cfi_offset %r14, -32
	.cfi_offset %r15, -24
	movl	%ecx, %ebx
	movq	%rdx, %r14
	movq	%rsi, %r15
	movq	%rdi, %r12
	xorl	%edi, %edi
	callq	__Z13get_global_idj
	cmpl	%ebx, %eax
	jae	LBB0_2
## %bb.1:
	cltq
	movss	(%r12,%rax,4), %xmm0    ## xmm0 = mem[0],zero,zero,zero
	addss	(%r15,%rax,4), %xmm0
	movss	%xmm0, (%r14,%rax,4)
LBB0_2:
	popq	%rbx
	popq	%r12
	popq	%r14
	popq	%r15
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_matAdd                 ## -- Begin function matAdd
	.p2align	4, 0x90
_matAdd:                                ## @matAdd
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%r15
	pushq	%r14
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	pushq	%rax
	.cfi_offset %rbx, -56
	.cfi_offset %r12, -48
	.cfi_offset %r13, -40
	.cfi_offset %r14, -32
	.cfi_offset %r15, -24
	movl	%r8d, %r14d
	movq	%rdx, %rbx
	movq	%rsi, %r15
	movq	%rdi, %r12
	xorl	%edi, %edi
	callq	__Z13get_global_idj
	movq	%rax, %r13
	movl	$1, %edi
	callq	__Z13get_global_idj
	movl	%eax, %edi
	movl	%r14d, %esi
	callq	__Z5mul24ii
	movl	%r13d, %edi
	movl	$4, %esi
	movl	%eax, %edx
	callq	__Z5mad24iii
	cltq
	movaps	(%r12,%rax,4), %xmm0
	addps	(%r15,%rax,4), %xmm0
	addps	(%rbx,%rax,4), %xmm0
	movaps	%xmm0, (%rbx,%rax,4)
	addq	$8, %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function

.subsections_via_symbols
