**Title:** `-Waggressive-loop-optimizations` false positive in `internal_multialloc_arrays` (GCC 15, `-O3`/`-fsplit-loops`)

---

Building `libs/container/src/alloc_lib.c` with GCC 15 at `-O3` produces a `-Waggressive-loop-optimizations` warning. I believe it is a false positive — the code looks correct to me — but the loop is written in a way GCC cannot bound, and a one-token change silences it without altering semantics.

### Reproducing

No build system needed — four lines from an empty directory:

```sh
curl -LO https://github.com/boostorg/boost/releases/download/boost-1.88.0/boost-1.88.0-cmake.tar.xz
tar xf boost-1.88.0-cmake.tar.xz
cd boost-1.88.0
gcc -O3 -c libs/container/src/alloc_lib.c -o /dev/null \
  -Ilibs/container/include -Ilibs/assert/include -Ilibs/config/include \
  -Ilibs/intrusive/include -Ilibs/move/include
```

### Warning

```
In file included from libs/container/src/dlmalloc_ext_2_8_6.c:52,
                 from libs/container/src/alloc_lib.c:24:
In function ‘internal_multialloc_arrays’,
    inlined from ‘boost_cont_multialloc_arrays’ at libs/container/src/dlmalloc_ext_2_8_6.c:1112:13:
libs/container/src/dlmalloc_ext_2_8_6.c:1085:41: warning: iteration 2305843009213693951 invokes undefined behavior [-Waggressive-loop-optimizations]
 1085 |                size = request2size(sizes[i]*element_size);
      |                                         ^
libs/container/src/dlmalloc_2_8_6.c:2231:6: note: in definition of macro ‘request2size’
 2231 |   (((req) < MIN_REQUEST)? MIN_CHUNK_SIZE : pad_request(req))
      |      ^~~
libs/container/src/dlmalloc_ext_2_8_6.c:1083:24: note: within this loop
 1083 |             for(++i; i != next_i; ++i) {
      |                      ~~^~~~~~~~~
```

### Environment

- Boost 1.88.0 (`boost-1.88.0-cmake.tar.xz` release archive)
- GCC 15.2.0 (Ubuntu 15.2.0-14ubuntu1~24~ppa1), x86-64, Ubuntu 24.04
- No `-Wall` — `-Waggressive-loop-optimizations` is on by default

### Which optimization levels

I compiled that one translation unit at each level. It is `-O3` only:

| Flags | Result |
|---|---|
| `-O0`, `-O1`, `-O2`, `-Os` | clean |
| `-O3` | warns |
| `-O2 -fsplit-loops` | warns |
| `-O3 -fno-split-loops` | clean |

So the trigger is specifically `-fsplit-loops`, which `-O3` enables and `-O2` does not. That also explains why this has stayed relatively quiet: distributions that build Boost at `-O2` never see it.

### Why I think it is a false positive

The number is the tell: 2305843009213693951 is 2^61 - 1, which is how many 8-byte `size_t`s could be indexed before running off any conceivable object. So GCC is not reporting a wrong computation — it is saying *"I cannot prove this loop exits before `sizes[i]` leaves the array, and if it did not, that would be UB."*

It cannot prove it because the condition is `i != next_i` on unsigned operands, so GCC must consider `i` incrementing past `next_i` and wrapping. The invariant that rules that out is established in the enclosing loop rather than locally:

```c
for(i = 0, next_i = 0; i != n_elements; i = next_i)
{
   int error = 0;
   size_t accum_size;
   for(accum_size = 0; next_i != n_elements; ++next_i){
      size_t cur_array_size = sizes[next_i];
      if(max_size < cur_array_size){
         error = 1;
         break;                       /* next_i unchanged, so next_i == i */
      }
      else{
         size_t reqsize = request2size(cur_array_size*element_size);
         if(((boost_cont_multialloc_segmented_malloc_size - CHUNK_OVERHEAD) - accum_size) < reqsize){
            if(!accum_size){
               accum_size += reqsize;
               ++next_i;              /* next_i == i + 1 */
            }
            break;                    /* else accum_size > 0, so next_i > i already */
         }
         accum_size += reqsize;
      }
   }

   mem = error ? 0 : mspace_malloc_lockless(m, accum_size - CHUNK_OVERHEAD);
   if (mem == 0){
      /* ... early return 0 ... */
   }
```

Each pass of the outer loop begins with `i == next_i`. The inner loop then either advances `next_i` at least once, or sets `error = 1` and leaves it equal to `i`. The `error` case takes `mem = 0` and returns before reaching the split-out loop. So by the time line 1083 runs, `next_i >= i + 1` always holds, and the trip count is exactly `next_i - i - 1`. No wraparound and no out-of-bounds access — but none of that is visible to the optimizer at the loop itself.

### Suggested fix

Changing the comparison gives GCC a bound it can prove, and is equivalent given `next_i >= i + 1`:

```diff
-            for(++i; i != next_i; ++i) {
+            for(++i; i < next_i; ++i) {
```

That is the only change needed to silence the diagnostic here. Happy to open a PR if you would like it in that form, or to test a different approach if you would rather keep `!=` and hint the bound another way.

### Workaround for anyone hitting this

The warning is confined to this one target, so it can be scoped rather than disabled globally. With the CMake build:

```cmake
if (TARGET boost_container)
    target_compile_options(boost_container PRIVATE
            $<$<AND:$<COMPILE_LANGUAGE:C>,$<C_COMPILER_ID:GNU>>:-Wno-aggressive-loop-optimizations>)
endif ()
```
