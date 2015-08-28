V24 utils
9 utils.f90 S582 0
08/27/2015  13:17:39
use mpi_constants public 0 indirect
use mpi_sizeofs public 0 indirect
use mpi_base public 0 indirect
use mpi public 0 indirect
use data_type public 0 direct
use rjmcmc_com public 0 direct
enduse
D 346 24 2304 408 2302 7
D 376 20 7
D 378 20 7
D 380 20 7
D 382 20 7
D 384 24 2338 808 2335 7
D 438 20 7
D 440 20 7
D 442 20 7
D 444 20 7
D 446 20 7
D 448 20 7
D 450 20 7
D 452 20 7
D 532 21 346 1 3 869 0 0 1 0 0
 0 868 3 3 869 869
D 535 18 23
D 537 21 346 1 3 871 0 0 1 0 0
 0 870 3 3 871 871
D 540 21 9 1 872 875 1 1 0 0 1
 3 873 3 3 873 874
D 543 21 9 1 3 879 0 0 1 0 0
 0 878 3 3 879 879
D 546 21 9 1 881 884 1 1 0 0 1
 3 882 3 3 882 883
D 549 21 9 1 3 887 0 0 1 0 0
 0 886 3 3 887 887
S 582 24 0 0 0 8 1 0 4667 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 utils
S 595 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 2302 25 180 rjmcmc_com objstruc
R 2304 5 182 rjmcmc_com par objstruc
R 2305 5 183 rjmcmc_com par$sd objstruc
R 2306 5 184 rjmcmc_com par$p objstruc
R 2307 5 185 rjmcmc_com par$o objstruc
R 2310 5 188 rjmcmc_com sdpar objstruc
R 2311 5 189 rjmcmc_com sdpar$sd objstruc
R 2312 5 190 rjmcmc_com sdpar$p objstruc
R 2313 5 191 rjmcmc_com sdpar$o objstruc
R 2316 5 194 rjmcmc_com z objstruc
R 2317 5 195 rjmcmc_com z$sd objstruc
R 2318 5 196 rjmcmc_com z$p objstruc
R 2319 5 197 rjmcmc_com z$o objstruc
R 2322 5 200 rjmcmc_com h objstruc
R 2323 5 201 rjmcmc_com h$sd objstruc
R 2324 5 202 rjmcmc_com h$p objstruc
R 2325 5 203 rjmcmc_com h$o objstruc
R 2327 5 205 rjmcmc_com k objstruc
R 2328 5 206 rjmcmc_com nfp objstruc
R 2329 5 207 rjmcmc_com nfail objstruc
R 2330 5 208 rjmcmc_com logl objstruc
R 2331 5 209 rjmcmc_com logp objstruc
R 2332 5 210 rjmcmc_com logwt objstruc
R 2333 5 211 rjmcmc_com logpr objstruc
R 2334 5 212 rjmcmc_com lognorm objstruc
R 2335 25 213 rjmcmc_com datastruc
R 2338 5 216 rjmcmc_com robs datastruc
R 2339 5 217 rjmcmc_com robs$sd datastruc
R 2340 5 218 rjmcmc_com robs$p datastruc
R 2341 5 219 rjmcmc_com robs$o datastruc
R 2344 5 222 rjmcmc_com angobs datastruc
R 2345 5 223 rjmcmc_com angobs$sd datastruc
R 2346 5 224 rjmcmc_com angobs$p datastruc
R 2347 5 225 rjmcmc_com angobs$o datastruc
R 2351 5 229 rjmcmc_com rrep datastruc
R 2352 5 230 rjmcmc_com rrep$sd datastruc
R 2353 5 231 rjmcmc_com rrep$p datastruc
R 2354 5 232 rjmcmc_com rrep$o datastruc
R 2358 5 236 rjmcmc_com res datastruc
R 2359 5 237 rjmcmc_com res$sd datastruc
R 2360 5 238 rjmcmc_com res$p datastruc
R 2361 5 239 rjmcmc_com res$o datastruc
R 2365 5 243 rjmcmc_com sigma datastruc
R 2366 5 244 rjmcmc_com sigma$sd datastruc
R 2367 5 245 rjmcmc_com sigma$p datastruc
R 2368 5 246 rjmcmc_com sigma$o datastruc
R 2371 5 249 rjmcmc_com lognorm datastruc
R 2372 5 250 rjmcmc_com lognorm$sd datastruc
R 2373 5 251 rjmcmc_com lognorm$p datastruc
R 2374 5 252 rjmcmc_com lognorm$o datastruc
R 2377 5 255 rjmcmc_com logdet datastruc
R 2378 5 256 rjmcmc_com logdet$sd datastruc
R 2379 5 257 rjmcmc_com logdet$p datastruc
R 2380 5 258 rjmcmc_com logdet$o datastruc
R 2383 5 261 rjmcmc_com ndpf datastruc
R 2384 5 262 rjmcmc_com ndpf$sd datastruc
R 2385 5 263 rjmcmc_com ndpf$p datastruc
R 2386 5 264 rjmcmc_com ndpf$o datastruc
R 2388 5 266 rjmcmc_com nang datastruc
R 2389 5 267 rjmcmc_com nband datastruc
S 2492 23 5 0 0 8 2499 582 15496 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 resample
S 2493 7 3 1 0 532 1 2492 15505 800204 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 oldsample
S 2494 6 3 1 0 6 1 2492 15515 800004 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 oldsize
S 2495 6 3 1 0 6 1 2492 15523 800004 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 newsize
S 2496 1 3 0 0 20 1 2492 15531 4 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 mode
S 2497 1 3 0 0 535 1 2492 15536 4 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 k_mode
S 2498 7 3 0 0 537 1 2492 15543 800204 1003000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 newsample
S 2499 14 5 0 0 537 1 2492 15496 204 1400000 0 0 0 1016 5 0 0 2498 0 0 0 0 0 0 0 0 0 6 0 582 0 0 0 0 resample
F 2499 5 2493 2494 2495 2496 2497
S 2500 6 1 0 0 6 1 2492 15553 40800006 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_868
S 2501 6 1 0 0 6 1 2492 15561 40800006 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_870
S 2502 23 5 0 0 8 2505 582 15569 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 cumsum
S 2503 7 3 1 0 540 1 2502 15576 20400004 10003000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array
S 2504 1 3 1 0 9 1 2502 2777 80000004 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 shift
S 2505 14 5 0 0 543 1 2502 15569 20000204 400000 0 0 0 1022 2 0 0 2506 0 0 0 0 0 0 0 0 0 162 0 582 0 0 0 0 cumsum
F 2505 2 2503 2504
S 2506 7 3 0 0 543 1 2502 15569 800204 1003000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 cumsum
S 2507 6 1 0 0 6 1 2502 15582 40800006 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_0_1
S 2508 6 1 0 0 6 1 2502 15590 40800006 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_2_1
S 2509 6 1 0 0 6 1 2502 15598 40800006 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_3_1
S 2510 6 1 0 0 6 1 2502 15606 40800006 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_877
S 2511 6 1 0 0 6 1 2502 15614 40800006 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_881
S 2513 23 5 0 0 6 2517 582 15633 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 logsum
S 2514 1 3 1 0 9 1 2513 12574 4 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 a
S 2515 1 3 1 0 9 1 2513 12576 4 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 b
S 2516 1 3 0 0 9 1 2513 12578 4 1003000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 c
S 2517 14 5 0 0 9 1 2513 15633 4 1400000 0 0 0 1025 2 0 0 2516 0 0 0 0 0 0 0 0 0 187 0 582 0 0 0 0 logsum
F 2517 2 2514 2515
S 2518 23 5 0 0 6 2520 582 15640 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 logcumsum
S 2519 7 3 1 0 546 1 2518 15576 20400004 10003000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array
S 2520 14 5 0 0 549 1 2518 15640 20000204 400000 0 0 0 1028 1 0 0 2521 0 0 0 0 0 0 0 0 0 201 0 582 0 0 0 0 logcumsum
F 2520 1 2519
S 2521 7 3 0 0 549 1 2518 15640 800204 1003000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 logcumsum
S 2522 6 1 0 0 6 1 2518 15582 40800006 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_0_1
S 2523 6 1 0 0 6 1 2518 15590 40800006 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_2_1
S 2524 6 1 0 0 6 1 2518 15598 40800006 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_3_1
S 2525 6 1 0 0 6 1 2518 15650 40800006 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_886
S 2526 6 1 0 0 6 1 2518 15658 40800006 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_890
A 23 2 0 0 0 6 595 0 0 0 23 0 0 0 0 0 0 0 0 0
A 868 1 0 0 387 6 2494 0 0 0 0 0 0 0 0 0 0 0 0 0
A 869 1 0 0 397 6 2500 0 0 0 0 0 0 0 0 0 0 0 0 0
A 870 1 0 0 388 6 2495 0 0 0 0 0 0 0 0 0 0 0 0 0
A 871 1 0 0 399 6 2501 0 0 0 0 0 0 0 0 0 0 0 0 0
A 872 1 0 0 408 6 2509 0 0 0 0 0 0 0 0 0 0 0 0 0
A 873 1 0 0 406 6 2507 0 0 0 0 0 0 0 0 0 0 0 0 0
A 874 1 0 0 751 6 2510 0 0 0 0 0 0 0 0 0 0 0 0 0
A 875 1 0 0 407 6 2508 0 0 0 0 0 0 0 0 0 0 0 0 0
A 876 1 0 0 297 0 426 0 0 0 0 0 0 0 0 0 0 0 0 0
A 877 1 0 7 328 540 2503 0 0 0 0 0 0 0 0 0 0 0 0 0
A 878 14 0 0 0 6 876 0 0 0 0 0 0 243 2 1 0 0 0 0
W 2 877 5
A 879 1 0 0 410 6 2511 0 0 0 0 0 0 0 0 0 0 0 0 0
A 881 1 0 0 433 6 2524 0 0 0 0 0 0 0 0 0 0 0 0 0
A 882 1 0 0 429 6 2522 0 0 0 0 0 0 0 0 0 0 0 0 0
A 883 1 0 0 435 6 2525 0 0 0 0 0 0 0 0 0 0 0 0 0
A 884 1 0 0 431 6 2523 0 0 0 0 0 0 0 0 0 0 0 0 0
A 885 1 0 9 209 546 2519 0 0 0 0 0 0 0 0 0 0 0 0 0
A 886 14 0 0 0 6 876 0 0 0 0 0 0 243 2 4 0 0 0 0
W 2 885 5
A 887 1 0 0 437 6 2526 0 0 0 0 0 0 0 0 0 0 0 0 0
Z
T 2302 346 0 0 0 0
A 2306 7 376 0 1 2 1
A 2312 7 378 0 1 2 1
A 2318 7 380 0 1 2 1
A 2324 7 382 0 1 2 0
T 2335 384 0 0 0 0
A 2340 7 438 0 1 2 1
A 2346 7 440 0 1 2 1
A 2353 7 442 0 1 2 1
A 2360 7 444 0 1 2 1
A 2367 7 446 0 1 2 1
A 2373 7 448 0 1 2 1
A 2379 7 450 0 1 2 1
A 2385 7 452 0 1 2 0
Z
