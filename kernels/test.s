	.section	__TEXT,__text,regular,pure_instructions
	.build_version macos, 10, 15	sdk_version 10, 15, 6
	.section	__TEXT,__literal4,4byte_literals
	.p2align	2               ## -- Begin function kernel1D
LCPI0_0:
	.long	1120403456              ## float 100
	.section	__TEXT,__text,regular,pure_instructions
	.globl	_kernel1D
	.p2align	4, 0x90
_kernel1D:                              ## @kernel1D
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%rbx
	pushq	%rax
	.cfi_offset %rbx, -24
	movq	%rdi, %rbx
	xorl	%edi, %edi
	callq	__Z13get_global_idj
	cvtsi2ss	%eax, %xmm0
	divss	LCPI0_0(%rip), %xmm0
	cltq
	movss	%xmm0, (%rbx,%rax,4)
	addq	$8, %rsp
	popq	%rbx
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.section	__TEXT,__literal4,4byte_literals
	.p2align	2               ## -- Begin function kernel2D
LCPI1_0:
	.long	1120403456              ## float 100
	.section	__TEXT,__text,regular,pure_instructions
	.globl	_kernel2D
	.p2align	4, 0x90
_kernel2D:                              ## @kernel2D
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%r15
	pushq	%r14
	pushq	%rbx
	pushq	%rax
	.cfi_offset %rbx, -40
	.cfi_offset %r14, -32
	.cfi_offset %r15, -24
	movq	%rdi, %r14
	xorl	%edi, %edi
	callq	__Z13get_global_idj
	movq	%rax, %r15
	movl	$1, %edi
	callq	__Z13get_global_idj
	movq	%rax, %rbx
	xorl	%edi, %edi
	callq	__Z15get_global_sizej
	imull	%ebx, %eax
	addl	%r15d, %eax
	cvtsi2ss	%eax, %xmm0
	divss	LCPI1_0(%rip), %xmm0
	cltq
	movss	%xmm0, (%r14,%rax,4)
	addq	$8, %rsp
	popq	%rbx
	popq	%r14
	popq	%r15
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_kernel3D               ## -- Begin function kernel3D
	.p2align	4, 0x90
_kernel3D:                              ## @kernel3D
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
	movq	%rdi, %r14
	xorl	%edi, %edi
	callq	__Z13get_global_idj
	movq	%rax, %r15
	movl	$1, %edi
	callq	__Z13get_global_idj
	movq	%rax, %r12
	movl	$2, %edi
	callq	__Z13get_global_idj
	movq	%rax, %r13
	xorl	%edi, %edi
	callq	__Z15get_global_sizej
	movq	%rax, %rbx
	movl	$1, %edi
	callq	__Z15get_global_sizej
	imull	%r12d, %ebx
	addl	%r15d, %ebx
	imull	%r13d, %eax
	addl	%ebx, %eax
	cltq
	imulq	$274877907, %rax, %rcx  ## imm = 0x10624DD3
	movq	%rcx, %rdx
	shrq	$63, %rdx
	sarq	$38, %rcx
	addl	%edx, %ecx
	cvtsi2ss	%ecx, %xmm0
	movss	%xmm0, (%r14,%rax,4)
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
