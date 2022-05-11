; ModuleID = 'add.cl'
source_filename = "add.cl"
target datalayout = "e-m:o-i64:64-f80:128-n8:16:32:64-S128"
target triple = "x86_64-apple-macosx10.15.0"

; Function Attrs: convergent nounwind ssp uwtable
define spir_kernel void @vecAdd(float* nocapture readonly, float* nocapture readonly, float* nocapture, i32) local_unnamed_addr #0 !kernel_arg_addr_space !5 !kernel_arg_access_qual !6 !kernel_arg_type !7 !kernel_arg_base_type !7 !kernel_arg_type_qual !8 {
  %5 = tail call i64 @_Z13get_global_idj(i32 0) #2
  %6 = trunc i64 %5 to i32
  %7 = icmp ult i32 %6, %3
  br i1 %7, label %8, label %17

8:                                                ; preds = %4
  %9 = shl i64 %5, 32
  %10 = ashr exact i64 %9, 32
  %11 = getelementptr inbounds float, float* %0, i64 %10
  %12 = load float, float* %11, align 4, !tbaa !9
  %13 = getelementptr inbounds float, float* %1, i64 %10
  %14 = load float, float* %13, align 4, !tbaa !9
  %15 = fadd float %12, %14
  %16 = getelementptr inbounds float, float* %2, i64 %10
  store float %15, float* %16, align 4, !tbaa !9
  br label %17

17:                                               ; preds = %8, %4
  ret void
}

; Function Attrs: convergent nounwind readnone
declare i64 @_Z13get_global_idj(i32) local_unnamed_addr #1

; Function Attrs: convergent nounwind ssp uwtable
define spir_kernel void @matAdd(float* nocapture readonly, float* nocapture readonly, float* nocapture, i32, i32) local_unnamed_addr #0 !kernel_arg_addr_space !13 !kernel_arg_access_qual !14 !kernel_arg_type !15 !kernel_arg_base_type !15 !kernel_arg_type_qual !16 {
  %6 = tail call i64 @_Z13get_global_idj(i32 0) #2
  %7 = trunc i64 %6 to i32
  %8 = tail call i64 @_Z13get_global_idj(i32 1) #2
  %9 = trunc i64 %8 to i32
  %10 = tail call i32 @_Z5mul24ii(i32 %9, i32 %4) #2
  %11 = tail call i32 @_Z5mad24iii(i32 %7, i32 4, i32 %10) #2
  %12 = sext i32 %11 to i64
  %13 = getelementptr inbounds float, float* %0, i64 %12
  %14 = bitcast float* %13 to <4 x float>*
  %15 = load <4 x float>, <4 x float>* %14, align 16, !tbaa !17
  %16 = getelementptr inbounds float, float* %1, i64 %12
  %17 = bitcast float* %16 to <4 x float>*
  %18 = load <4 x float>, <4 x float>* %17, align 16, !tbaa !17
  %19 = fadd <4 x float> %15, %18
  %20 = getelementptr inbounds float, float* %2, i64 %12
  %21 = bitcast float* %20 to <4 x float>*
  %22 = load <4 x float>, <4 x float>* %21, align 16, !tbaa !17
  %23 = fadd <4 x float> %22, %19
  store <4 x float> %23, <4 x float>* %21, align 16, !tbaa !17
  ret void
}

; Function Attrs: convergent nounwind readnone
declare i32 @_Z5mul24ii(i32, i32) local_unnamed_addr #1

; Function Attrs: convergent nounwind readnone
declare i32 @_Z5mad24iii(i32, i32, i32) local_unnamed_addr #1

attributes #0 = { convergent nounwind ssp uwtable "correctly-rounded-divide-sqrt-fp-math"="false" "darwin-stkchk-strong-link" "denorms-are-zero"="false" "disable-tail-calls"="false" "frame-pointer"="all" "less-precise-fpmad"="false" "min-legal-vector-width"="0" "no-infs-fp-math"="false" "no-jump-tables"="false" "no-nans-fp-math"="false" "no-signed-zeros-fp-math"="false" "no-trapping-math"="false" "probe-stack"="___chkstk_darwin" "stack-protector-buffer-size"="8" "target-cpu"="penryn" "target-features"="+cx16,+cx8,+fxsr,+mmx,+sahf,+sse,+sse2,+sse3,+sse4.1,+ssse3,+x87" "uniform-work-group-size"="true" "unsafe-fp-math"="false" "use-soft-float"="false" }
attributes #1 = { convergent nounwind readnone "correctly-rounded-divide-sqrt-fp-math"="false" "darwin-stkchk-strong-link" "denorms-are-zero"="false" "disable-tail-calls"="false" "frame-pointer"="all" "less-precise-fpmad"="false" "no-infs-fp-math"="false" "no-nans-fp-math"="false" "no-signed-zeros-fp-math"="false" "no-trapping-math"="false" "probe-stack"="___chkstk_darwin" "stack-protector-buffer-size"="8" "target-cpu"="penryn" "target-features"="+cx16,+cx8,+fxsr,+mmx,+sahf,+sse,+sse2,+sse3,+sse4.1,+ssse3,+x87" "unsafe-fp-math"="false" "use-soft-float"="false" }
attributes #2 = { convergent nounwind readnone }

!llvm.module.flags = !{!0, !1, !2}
!opencl.ocl.version = !{!3}
!llvm.ident = !{!4}

!0 = !{i32 2, !"SDK Version", [3 x i32] [i32 10, i32 15, i32 6]}
!1 = !{i32 1, !"wchar_size", i32 4}
!2 = !{i32 7, !"PIC Level", i32 2}
!3 = !{i32 1, i32 0}
!4 = !{!"Apple clang version 11.0.3 (clang-1103.0.32.62)"}
!5 = !{i32 1, i32 1, i32 1, i32 0}
!6 = !{!"none", !"none", !"none", !"none"}
!7 = !{!"float*", !"float*", !"float*", !"uint"}
!8 = !{!"", !"", !"", !""}
!9 = !{!10, !10, i64 0}
!10 = !{!"float", !11, i64 0}
!11 = !{!"omnipotent char", !12, i64 0}
!12 = !{!"Simple C/C++ TBAA"}
!13 = !{i32 1, i32 1, i32 1, i32 0, i32 0}
!14 = !{!"none", !"none", !"none", !"none", !"none"}
!15 = !{!"float*", !"float*", !"float*", !"int", !"int"}
!16 = !{!"", !"", !"", !"", !""}
!17 = !{!11, !11, i64 0}
