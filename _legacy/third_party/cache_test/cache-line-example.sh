export compile="aarch64-linux-android-gcc -pie -fPIC -fPIE -O3"
export upload="adb push ./cache-line-example /data/local/tmp"
export run="adb shell /data/local/tmp/cache-line-example"
$compile -DBATCH_SIZE=1                         -DCACHE_LINE_FRIENDLY     -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run
$compile -DBATCH_SIZE=2                         -DCACHE_LINE_FRIENDLY     -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=4                         -DCACHE_LINE_FRIENDLY     -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=8                         -DCACHE_LINE_FRIENDLY     -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=16                        -DCACHE_LINE_FRIENDLY     -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=32                        -DCACHE_LINE_FRIENDLY     -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=64                        -DCACHE_LINE_FRIENDLY     -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=128                       -DCACHE_LINE_FRIENDLY     -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=256                       -DCACHE_LINE_FRIENDLY     -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=1   -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLY     -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=2   -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLY     -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=4   -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLY     -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=8   -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLY     -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=16  -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLY     -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=32  -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLY     -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=64  -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLY     -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=128 -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLY     -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=256 -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLY     -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
                 
$compile -DBATCH_SIZE=1                         -DCACHE_LINE_FRIENDLY_ISH -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=2                         -DCACHE_LINE_FRIENDLY_ISH -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=4                         -DCACHE_LINE_FRIENDLY_ISH -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=8                         -DCACHE_LINE_FRIENDLY_ISH -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=16                        -DCACHE_LINE_FRIENDLY_ISH -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=32                        -DCACHE_LINE_FRIENDLY_ISH -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=64                        -DCACHE_LINE_FRIENDLY_ISH -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=128                       -DCACHE_LINE_FRIENDLY_ISH -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=256                       -DCACHE_LINE_FRIENDLY_ISH -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=1   -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLY_ISH -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=2   -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLY_ISH -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=4   -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLY_ISH -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=8   -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLY_ISH -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=16  -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLY_ISH -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=32  -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLY_ISH -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=64  -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLY_ISH -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=128 -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLY_ISH -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=256 -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLY_ISH -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 

$compile -DBATCH_SIZE=1                         -DCACHE_LINE_FRIENDLIER   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=2                         -DCACHE_LINE_FRIENDLIER   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=4                         -DCACHE_LINE_FRIENDLIER   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=8                         -DCACHE_LINE_FRIENDLIER   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=16                        -DCACHE_LINE_FRIENDLIER   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=32                        -DCACHE_LINE_FRIENDLIER   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=64                        -DCACHE_LINE_FRIENDLIER   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=128                       -DCACHE_LINE_FRIENDLIER   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=256                       -DCACHE_LINE_FRIENDLIER   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=1   -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLIER   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=2   -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLIER   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=4   -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLIER   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=8   -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLIER   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=16  -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLIER   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=32  -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLIER   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=64  -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLIER   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=128 -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLIER   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=256 -DCACHE_LINE_PREFETCH -DCACHE_LINE_FRIENDLIER   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
                 
$compile -DBATCH_SIZE=1                         -DCACHE_LINE_UNFRIENDLY   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=2                         -DCACHE_LINE_UNFRIENDLY   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=4                         -DCACHE_LINE_UNFRIENDLY   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=8                         -DCACHE_LINE_UNFRIENDLY   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=16                        -DCACHE_LINE_UNFRIENDLY   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=32                        -DCACHE_LINE_UNFRIENDLY   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=64                        -DCACHE_LINE_UNFRIENDLY   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=128                       -DCACHE_LINE_UNFRIENDLY   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=256                       -DCACHE_LINE_UNFRIENDLY   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=1   -DCACHE_LINE_PREFETCH -DCACHE_LINE_UNFRIENDLY   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=2   -DCACHE_LINE_PREFETCH -DCACHE_LINE_UNFRIENDLY   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=4   -DCACHE_LINE_PREFETCH -DCACHE_LINE_UNFRIENDLY   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=8   -DCACHE_LINE_PREFETCH -DCACHE_LINE_UNFRIENDLY   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=16  -DCACHE_LINE_PREFETCH -DCACHE_LINE_UNFRIENDLY   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=32  -DCACHE_LINE_PREFETCH -DCACHE_LINE_UNFRIENDLY   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=64  -DCACHE_LINE_PREFETCH -DCACHE_LINE_UNFRIENDLY   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=128 -DCACHE_LINE_PREFETCH -DCACHE_LINE_UNFRIENDLY   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
$compile -DBATCH_SIZE=256 -DCACHE_LINE_PREFETCH -DCACHE_LINE_UNFRIENDLY   -DCACHE_LINE_SIZE=64 -o cache-line-example cache-line-example.c && $upload &> push.log && $run 
