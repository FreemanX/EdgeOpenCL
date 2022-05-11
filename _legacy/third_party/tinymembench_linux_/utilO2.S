
util.o:     file format elf64-littleaarch64


Disassembly of section .text:

0000000000000000 <_Z18aligned_block_copyPlS_i>:
   0:	7100fc5f 	cmp	w2, #0x3f
   4:	5400032d 	b.le	68 <_Z18aligned_block_copyPlS_i+0x68>
   8:	51010045 	sub	w5, w2, #0x40
   c:	927a64a5 	and	x5, x5, #0xffffffc0
  10:	910100a5 	add	x5, x5, #0x40
  14:	8b050025 	add	x5, x1, x5
  18:	f9400022 	ldr	x2, [x1]
  1c:	91010021 	add	x1, x1, #0x40
  20:	f9000002 	str	x2, [x0]
  24:	91010000 	add	x0, x0, #0x40
  28:	eb0100bf 	cmp	x5, x1
  2c:	f85c8022 	ldur	x2, [x1,#-56]
  30:	f81c8002 	stur	x2, [x0,#-56]
  34:	f85d0022 	ldur	x2, [x1,#-48]
  38:	f81d0002 	stur	x2, [x0,#-48]
  3c:	f85d8022 	ldur	x2, [x1,#-40]
  40:	f81d8002 	stur	x2, [x0,#-40]
  44:	f85e0024 	ldur	x4, [x1,#-32]
  48:	f81e0004 	stur	x4, [x0,#-32]
  4c:	f85e8023 	ldur	x3, [x1,#-24]
  50:	f85f0022 	ldur	x2, [x1,#-16]
  54:	f81e8003 	stur	x3, [x0,#-24]
  58:	f81f0002 	stur	x2, [x0,#-16]
  5c:	f85f8022 	ldur	x2, [x1,#-8]
  60:	f81f8002 	stur	x2, [x0,#-8]
  64:	54fffda1 	b.ne	18 <_Z18aligned_block_copyPlS_i+0x18>
  68:	d65f03c0 	ret

000000000000006c <_Z28aligned_block_copy_backwardsPlS_i>:
  6c:	6b1f005f 	cmp	w2, wzr
  70:	11001c43 	add	w3, w2, #0x7
  74:	1a82b063 	csel	w3, w3, w2, lt
  78:	7100fc5f 	cmp	w2, #0x3f
  7c:	13037c63 	asr	w3, w3, #3
  80:	928000e4 	mov	x4, #0xfffffffffffffff8    	// #-8
  84:	8b23cc83 	add	x3, x4, w3, sxtw #3
  88:	8b030021 	add	x1, x1, x3
  8c:	8b030000 	add	x0, x0, x3
  90:	5400036d 	b.le	fc <_Z28aligned_block_copy_backwardsPlS_i+0x90>
  94:	51010042 	sub	w2, w2, #0x40
  98:	928007e3 	mov	x3, #0xffffffffffffffc0    	// #-64
  9c:	53067c44 	lsr	w4, w2, #6
  a0:	52800802 	mov	w2, #0x40                  	// #64
  a4:	9ba28c82 	umsubl	x2, w4, w2, x3
  a8:	8b020022 	add	x2, x1, x2
  ac:	f9400023 	ldr	x3, [x1]
  b0:	d1010021 	sub	x1, x1, #0x40
  b4:	f9000003 	str	x3, [x0]
  b8:	d1010000 	sub	x0, x0, #0x40
  bc:	eb01005f 	cmp	x2, x1
  c0:	f9401c23 	ldr	x3, [x1,#56]
  c4:	f9001c03 	str	x3, [x0,#56]
  c8:	f9401823 	ldr	x3, [x1,#48]
  cc:	f9001803 	str	x3, [x0,#48]
  d0:	f9401423 	ldr	x3, [x1,#40]
  d4:	f9001403 	str	x3, [x0,#40]
  d8:	f9401025 	ldr	x5, [x1,#32]
  dc:	f9001005 	str	x5, [x0,#32]
  e0:	f9400c24 	ldr	x4, [x1,#24]
  e4:	f9400823 	ldr	x3, [x1,#16]
  e8:	f9000c04 	str	x4, [x0,#24]
  ec:	f9000803 	str	x3, [x0,#16]
  f0:	f9400423 	ldr	x3, [x1,#8]
  f4:	f9000403 	str	x3, [x0,#8]
  f8:	54fffda1 	b.ne	ac <_Z28aligned_block_copy_backwardsPlS_i+0x40>
  fc:	d65f03c0 	ret

0000000000000100 <_Z33aligned_block_copy_backwards_bs32PlS_i>:
 100:	6b1f005f 	cmp	w2, wzr
 104:	11001c43 	add	w3, w2, #0x7
 108:	1a82b063 	csel	w3, w3, w2, lt
 10c:	7100fc5f 	cmp	w2, #0x3f
 110:	928007e4 	mov	x4, #0xffffffffffffffc0    	// #-64
 114:	13037c63 	asr	w3, w3, #3
 118:	8b23cc83 	add	x3, x4, w3, sxtw #3
 11c:	8b030021 	add	x1, x1, x3
 120:	8b030000 	add	x0, x0, x3
 124:	5400034d 	b.le	18c <_Z33aligned_block_copy_backwards_bs32PlS_i+0x8c>
 128:	51010042 	sub	w2, w2, #0x40
 12c:	52800803 	mov	w3, #0x40                  	// #64
 130:	53067c42 	lsr	w2, w2, #6
 134:	9ba39042 	umsubl	x2, w2, w3, x4
 138:	8b020022 	add	x2, x1, x2
 13c:	f9401023 	ldr	x3, [x1,#32]
 140:	d1010021 	sub	x1, x1, #0x40
 144:	f9001003 	str	x3, [x0,#32]
 148:	d1010000 	sub	x0, x0, #0x40
 14c:	eb02003f 	cmp	x1, x2
 150:	f9403423 	ldr	x3, [x1,#104]
 154:	f9003403 	str	x3, [x0,#104]
 158:	f9403823 	ldr	x3, [x1,#112]
 15c:	f9003803 	str	x3, [x0,#112]
 160:	f9403c23 	ldr	x3, [x1,#120]
 164:	f9003c03 	str	x3, [x0,#120]
 168:	f9402023 	ldr	x3, [x1,#64]
 16c:	f9002003 	str	x3, [x0,#64]
 170:	f9402423 	ldr	x3, [x1,#72]
 174:	f9002403 	str	x3, [x0,#72]
 178:	f9402824 	ldr	x4, [x1,#80]
 17c:	f9402c23 	ldr	x3, [x1,#88]
 180:	f9002804 	str	x4, [x0,#80]
 184:	f9002c03 	str	x3, [x0,#88]
 188:	54fffda1 	b.ne	13c <_Z33aligned_block_copy_backwards_bs32PlS_i+0x3c>
 18c:	d65f03c0 	ret

0000000000000190 <_Z33aligned_block_copy_backwards_bs64PlS_i>:
 190:	6b1f005f 	cmp	w2, wzr
 194:	11001c43 	add	w3, w2, #0x7
 198:	1a82b063 	csel	w3, w3, w2, lt
 19c:	7100fc5f 	cmp	w2, #0x3f
 1a0:	928007e4 	mov	x4, #0xffffffffffffffc0    	// #-64
 1a4:	13037c63 	asr	w3, w3, #3
 1a8:	8b23cc83 	add	x3, x4, w3, sxtw #3
 1ac:	8b030021 	add	x1, x1, x3
 1b0:	8b030000 	add	x0, x0, x3
 1b4:	5400034d 	b.le	21c <_Z33aligned_block_copy_backwards_bs64PlS_i+0x8c>
 1b8:	51010042 	sub	w2, w2, #0x40
 1bc:	52800803 	mov	w3, #0x40                  	// #64
 1c0:	53067c42 	lsr	w2, w2, #6
 1c4:	9ba39042 	umsubl	x2, w2, w3, x4
 1c8:	8b020022 	add	x2, x1, x2
 1cc:	f9400023 	ldr	x3, [x1]
 1d0:	d1010021 	sub	x1, x1, #0x40
 1d4:	f9000003 	str	x3, [x0]
 1d8:	d1010000 	sub	x0, x0, #0x40
 1dc:	eb02003f 	cmp	x1, x2
 1e0:	f9402423 	ldr	x3, [x1,#72]
 1e4:	f9002403 	str	x3, [x0,#72]
 1e8:	f9402823 	ldr	x3, [x1,#80]
 1ec:	f9002803 	str	x3, [x0,#80]
 1f0:	f9402c23 	ldr	x3, [x1,#88]
 1f4:	f9002c03 	str	x3, [x0,#88]
 1f8:	f9403023 	ldr	x3, [x1,#96]
 1fc:	f9003003 	str	x3, [x0,#96]
 200:	f9403423 	ldr	x3, [x1,#104]
 204:	f9003403 	str	x3, [x0,#104]
 208:	f9403824 	ldr	x4, [x1,#112]
 20c:	f9403c23 	ldr	x3, [x1,#120]
 210:	f9003804 	str	x4, [x0,#112]
 214:	f9003c03 	str	x3, [x0,#120]
 218:	54fffda1 	b.ne	1cc <_Z33aligned_block_copy_backwards_bs64PlS_i+0x3c>
 21c:	d65f03c0 	ret

0000000000000220 <_Z23aligned_block_copy_pf32PlS_i>:
 220:	7100fc5f 	cmp	w2, #0x3f
 224:	5400034d 	b.le	28c <_Z23aligned_block_copy_pf32PlS_i+0x6c>
 228:	51010044 	sub	w4, w2, #0x40
 22c:	91040023 	add	x3, x1, #0x100
 230:	927a6484 	and	x4, x4, #0xffffffc0
 234:	91050084 	add	x4, x4, #0x140
 238:	8b040021 	add	x1, x1, x4
 23c:	91010000 	add	x0, x0, #0x40
 240:	f8500062 	ldur	x2, [x3,#-256]
 244:	f8530064 	ldur	x4, [x3,#-208]
 248:	f81c0002 	stur	x2, [x0,#-64]
 24c:	f8508062 	ldur	x2, [x3,#-248]
 250:	f81c8002 	stur	x2, [x0,#-56]
 254:	f8510062 	ldur	x2, [x3,#-240]
 258:	f81d0002 	stur	x2, [x0,#-48]
 25c:	f8518062 	ldur	x2, [x3,#-232]
 260:	f81d8002 	stur	x2, [x0,#-40]
 264:	f8520062 	ldur	x2, [x3,#-224]
 268:	f81e0002 	stur	x2, [x0,#-32]
 26c:	f8528062 	ldur	x2, [x3,#-216]
 270:	f81e8002 	stur	x2, [x0,#-24]
 274:	f81f0004 	stur	x4, [x0,#-16]
 278:	f8538062 	ldur	x2, [x3,#-200]
 27c:	91010063 	add	x3, x3, #0x40
 280:	eb01007f 	cmp	x3, x1
 284:	f81f8002 	stur	x2, [x0,#-8]
 288:	54fffda1 	b.ne	23c <_Z23aligned_block_copy_pf32PlS_i+0x1c>
 28c:	d65f03c0 	ret

0000000000000290 <_Z23aligned_block_copy_pf64PlS_i>:
 290:	7100fc5f 	cmp	w2, #0x3f
 294:	5400034d 	b.le	2fc <_Z23aligned_block_copy_pf64PlS_i+0x6c>
 298:	51010044 	sub	w4, w2, #0x40
 29c:	91040023 	add	x3, x1, #0x100
 2a0:	927a6484 	and	x4, x4, #0xffffffc0
 2a4:	91050084 	add	x4, x4, #0x140
 2a8:	8b040021 	add	x1, x1, x4
 2ac:	91010000 	add	x0, x0, #0x40
 2b0:	f8500062 	ldur	x2, [x3,#-256]
 2b4:	f8530064 	ldur	x4, [x3,#-208]
 2b8:	f81c0002 	stur	x2, [x0,#-64]
 2bc:	f8508062 	ldur	x2, [x3,#-248]
 2c0:	f81c8002 	stur	x2, [x0,#-56]
 2c4:	f8510062 	ldur	x2, [x3,#-240]
 2c8:	f81d0002 	stur	x2, [x0,#-48]
 2cc:	f8518062 	ldur	x2, [x3,#-232]
 2d0:	f81d8002 	stur	x2, [x0,#-40]
 2d4:	f8520062 	ldur	x2, [x3,#-224]
 2d8:	f81e0002 	stur	x2, [x0,#-32]
 2dc:	f8528062 	ldur	x2, [x3,#-216]
 2e0:	f81e8002 	stur	x2, [x0,#-24]
 2e4:	f81f0004 	stur	x4, [x0,#-16]
 2e8:	f8538062 	ldur	x2, [x3,#-200]
 2ec:	91010063 	add	x3, x3, #0x40
 2f0:	eb01007f 	cmp	x3, x1
 2f4:	f81f8002 	stur	x2, [x0,#-8]
 2f8:	54fffda1 	b.ne	2ac <_Z23aligned_block_copy_pf64PlS_i+0x1c>
 2fc:	d65f03c0 	ret

0000000000000300 <_Z18aligned_block_fillPlS_i>:
 300:	7100fc5f 	cmp	w2, #0x3f
 304:	f9400021 	ldr	x1, [x1]
 308:	5400020d 	b.le	348 <_Z18aligned_block_fillPlS_i+0x48>
 30c:	51010042 	sub	w2, w2, #0x40
 310:	927a6442 	and	x2, x2, #0xffffffc0
 314:	91010042 	add	x2, x2, #0x40
 318:	8b020002 	add	x2, x0, x2
 31c:	f9000001 	str	x1, [x0]
 320:	91010000 	add	x0, x0, #0x40
 324:	eb00005f 	cmp	x2, x0
 328:	f81c8001 	stur	x1, [x0,#-56]
 32c:	f81d0001 	stur	x1, [x0,#-48]
 330:	f81d8001 	stur	x1, [x0,#-40]
 334:	f81e0001 	stur	x1, [x0,#-32]
 338:	f81e8001 	stur	x1, [x0,#-24]
 33c:	f81f0001 	stur	x1, [x0,#-16]
 340:	f81f8001 	stur	x1, [x0,#-8]
 344:	54fffec1 	b.ne	31c <_Z18aligned_block_fillPlS_i+0x1c>
 348:	d65f03c0 	ret

000000000000034c <_Z28aligned_block_fill_shuffle16PlS_i>:
 34c:	7100fc5f 	cmp	w2, #0x3f
 350:	f9400021 	ldr	x1, [x1]
 354:	5400020d 	b.le	394 <_Z28aligned_block_fill_shuffle16PlS_i+0x48>
 358:	51010042 	sub	w2, w2, #0x40
 35c:	927a6442 	and	x2, x2, #0xffffffc0
 360:	91010042 	add	x2, x2, #0x40
 364:	8b020002 	add	x2, x0, x2
 368:	f9000001 	str	x1, [x0]
 36c:	91010000 	add	x0, x0, #0x40
 370:	eb02001f 	cmp	x0, x2
 374:	f81c8001 	stur	x1, [x0,#-56]
 378:	f81d8001 	stur	x1, [x0,#-40]
 37c:	f81d0001 	stur	x1, [x0,#-48]
 380:	f81e8001 	stur	x1, [x0,#-24]
 384:	f81e0001 	stur	x1, [x0,#-32]
 388:	f81f0001 	stur	x1, [x0,#-16]
 38c:	f81f8001 	stur	x1, [x0,#-8]
 390:	54fffec1 	b.ne	368 <_Z28aligned_block_fill_shuffle16PlS_i+0x1c>
 394:	d65f03c0 	ret

0000000000000398 <_Z28aligned_block_fill_shuffle32PlS_i>:
 398:	7100fc5f 	cmp	w2, #0x3f
 39c:	f9400021 	ldr	x1, [x1]
 3a0:	5400020d 	b.le	3e0 <_Z28aligned_block_fill_shuffle32PlS_i+0x48>
 3a4:	51010042 	sub	w2, w2, #0x40
 3a8:	927a6442 	and	x2, x2, #0xffffffc0
 3ac:	91010042 	add	x2, x2, #0x40
 3b0:	8b020002 	add	x2, x0, x2
 3b4:	f9000c01 	str	x1, [x0,#24]
 3b8:	91010000 	add	x0, x0, #0x40
 3bc:	eb02001f 	cmp	x0, x2
 3c0:	f81c0001 	stur	x1, [x0,#-64]
 3c4:	f81d0001 	stur	x1, [x0,#-48]
 3c8:	f81c8001 	stur	x1, [x0,#-56]
 3cc:	f81f8001 	stur	x1, [x0,#-8]
 3d0:	f81e0001 	stur	x1, [x0,#-32]
 3d4:	f81f0001 	stur	x1, [x0,#-16]
 3d8:	f81e8001 	stur	x1, [x0,#-24]
 3dc:	54fffec1 	b.ne	3b4 <_Z28aligned_block_fill_shuffle32PlS_i+0x1c>
 3e0:	d65f03c0 	ret

00000000000003e4 <_Z28aligned_block_fill_shuffle64PlS_i>:
 3e4:	7100fc5f 	cmp	w2, #0x3f
 3e8:	f9400021 	ldr	x1, [x1]
 3ec:	5400020d 	b.le	42c <_Z28aligned_block_fill_shuffle64PlS_i+0x48>
 3f0:	51010042 	sub	w2, w2, #0x40
 3f4:	927a6442 	and	x2, x2, #0xffffffc0
 3f8:	91010042 	add	x2, x2, #0x40
 3fc:	8b020002 	add	x2, x0, x2
 400:	f9001401 	str	x1, [x0,#40]
 404:	91010000 	add	x0, x0, #0x40
 408:	eb02001f 	cmp	x0, x2
 40c:	f81d0001 	stur	x1, [x0,#-48]
 410:	f81f8001 	stur	x1, [x0,#-8]
 414:	f81f0001 	stur	x1, [x0,#-16]
 418:	f81c8001 	stur	x1, [x0,#-56]
 41c:	f81d8001 	stur	x1, [x0,#-40]
 420:	f81c0001 	stur	x1, [x0,#-64]
 424:	f81e0001 	stur	x1, [x0,#-32]
 428:	54fffec1 	b.ne	400 <_Z28aligned_block_fill_shuffle64PlS_i+0x1c>
 42c:	d65f03c0 	ret

0000000000000430 <_Z7gettimev>:
 430:	a9be7bfd 	stp	x29, x30, [sp,#-32]!
 434:	d2800001 	mov	x1, #0x0                   	// #0
 438:	910003fd 	mov	x29, sp
 43c:	910043a0 	add	x0, x29, #0x10
 440:	94000000 	bl	0 <gettimeofday>
 444:	f9400ba0 	ldr	x0, [x29,#16]
 448:	5c000181 	ldr	d1, 478 <_Z7gettimev+0x48>
 44c:	d378dc01 	lsl	x1, x0, #8
 450:	cb000c21 	sub	x1, x1, x0, lsl #3
 454:	8b011800 	add	x0, x0, x1, lsl #6
 458:	cb010000 	sub	x0, x0, x1
 45c:	f9400fa1 	ldr	x1, [x29,#24]
 460:	a8c27bfd 	ldp	x29, x30, [sp],#32
 464:	8b001820 	add	x0, x1, x0, lsl #6
 468:	9e620000 	scvtf	d0, x0
 46c:	1e611800 	fdiv	d0, d0, d1
 470:	d65f03c0 	ret
 474:	d503201f 	nop
 478:	00000000 	.word	0x00000000
 47c:	412e8480 	.word	0x412e8480

0000000000000480 <_Z4fmindd>:
 480:	1e612010 	fcmpe	d0, d1
 484:	54000044 	b.mi	48c <_Z4fmindd+0xc>
 488:	1e604020 	fmov	d0, d1
 48c:	d65f03c0 	ret

0000000000000490 <_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i>:
 490:	a9ba7bfd 	stp	x29, x30, [sp,#-96]!
 494:	910003fd 	mov	x29, sp
 498:	a90153f3 	stp	x19, x20, [sp,#16]
 49c:	f90013f5 	str	x21, [sp,#32]
 4a0:	2a0103f3 	mov	w19, w1
 4a4:	aa0003f5 	mov	x21, x0
 4a8:	aa0203f4 	mov	x20, x2
 4ac:	b40008c0 	cbz	x0, 5c4 <_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i+0x134>
 4b0:	37f808a1 	tbnz	w1, #31, 5c4 <_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i+0x134>
 4b4:	b4000854 	cbz	x20, 5bc <_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i+0x12c>
 4b8:	37f80823 	tbnz	w3, #31, 5bc <_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i+0x12c>
 4bc:	b4000784 	cbz	x4, 5ac <_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i+0x11c>
 4c0:	37f80765 	tbnz	w5, #31, 5ac <_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i+0x11c>
 4c4:	b4000786 	cbz	x6, 5b4 <_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i+0x124>
 4c8:	37f80767 	tbnz	w7, #31, 5b4 <_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i+0x124>
 4cc:	0b030260 	add	w0, w19, w3
 4d0:	f9001fa6 	str	x6, [x29,#56]
 4d4:	0b050000 	add	w0, w0, w5
 4d8:	f90023a4 	str	x4, [x29,#64]
 4dc:	0b070007 	add	w7, w0, w7
 4e0:	f90027a3 	str	x3, [x29,#72]
 4e4:	116400e2 	add	w2, w7, #0x900, lsl #12
 4e8:	f9002ba5 	str	x5, [x29,#80]
 4ec:	93407c42 	sxtw	x2, w2
 4f0:	f9002fa2 	str	x2, [x29,#88]
 4f4:	aa0203e0 	mov	x0, x2
 4f8:	94000000 	bl	0 <malloc>
 4fc:	f9402fa2 	ldr	x2, [x29,#88]
 500:	52801981 	mov	w1, #0xcc                  	// #204
 504:	94000000 	bl	0 <memset>
 508:	aa0003e7 	mov	x7, x0
 50c:	9143fc00 	add	x0, x0, #0xff, lsl #12
 510:	f9402ba5 	ldr	x5, [x29,#80]
 514:	913ffc00 	add	x0, x0, #0xfff
 518:	f94027a3 	ldr	x3, [x29,#72]
 51c:	926cac00 	and	x0, x0, #0xfffffffffff00000
 520:	f94023a4 	ldr	x4, [x29,#64]
 524:	f9401fa6 	ldr	x6, [x29,#56]
 528:	b4000115 	cbz	x21, 548 <_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i+0xb8>
 52c:	9142a800 	add	x0, x0, #0xaa, lsl #12
 530:	912a0000 	add	x0, x0, #0xa80
 534:	8b33c013 	add	x19, x0, w19, sxtw
 538:	9143fe73 	add	x19, x19, #0xff, lsl #12
 53c:	f90002a0 	str	x0, [x21]
 540:	913ffe73 	add	x19, x19, #0xfff
 544:	926cae60 	and	x0, x19, #0xfffffffffff00000
 548:	b4000114 	cbz	x20, 568 <_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i+0xd8>
 54c:	91415400 	add	x0, x0, #0x55, lsl #12
 550:	91140000 	add	x0, x0, #0x500
 554:	8b23c003 	add	x3, x0, w3, sxtw
 558:	9143fc63 	add	x3, x3, #0xff, lsl #12
 55c:	f9000280 	str	x0, [x20]
 560:	913ffc63 	add	x3, x3, #0xfff
 564:	926cac60 	and	x0, x3, #0xfffffffffff00000
 568:	b4000104 	cbz	x4, 588 <_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i+0xf8>
 56c:	91433000 	add	x0, x0, #0xcc, lsl #12
 570:	91320000 	add	x0, x0, #0xc80
 574:	8b25c005 	add	x5, x0, w5, sxtw
 578:	9143fca5 	add	x5, x5, #0xff, lsl #12
 57c:	f9000080 	str	x0, [x4]
 580:	913ffca5 	add	x5, x5, #0xfff
 584:	926caca0 	and	x0, x5, #0xfffffffffff00000
 588:	b4000086 	cbz	x6, 598 <_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i+0x108>
 58c:	9140cc00 	add	x0, x0, #0x33, lsl #12
 590:	910c0000 	add	x0, x0, #0x300
 594:	f90000c0 	str	x0, [x6]
 598:	a94153f3 	ldp	x19, x20, [sp,#16]
 59c:	f94013f5 	ldr	x21, [sp,#32]
 5a0:	aa0703e0 	mov	x0, x7
 5a4:	a8c67bfd 	ldp	x29, x30, [sp],#96
 5a8:	d65f03c0 	ret
 5ac:	52800005 	mov	w5, #0x0                   	// #0
 5b0:	b5fff8c6 	cbnz	x6, 4c8 <_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i+0x38>
 5b4:	52800007 	mov	w7, #0x0                   	// #0
 5b8:	17ffffc5 	b	4cc <_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i+0x3c>
 5bc:	52800003 	mov	w3, #0x0                   	// #0
 5c0:	17ffffbf 	b	4bc <_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i+0x2c>
 5c4:	52800013 	mov	w19, #0x0                   	// #0
 5c8:	17ffffbb 	b	4b4 <_Z29alloc_four_nonaliased_buffersPPviS0_iS0_iS0_i+0x24>
