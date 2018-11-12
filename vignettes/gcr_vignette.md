# Analysis

We will have a data namely “Yan” for the example of LTMG - GCR
(biclustering) pipeline.

Basically, we may need the following steps for this
analysis

## (i) a standard data loading function

``` r
data0 <- log(as.matrix(read.delim("Yan_expression_RPKM.txt", row.names = 1)))
```

## (ii) running LTMG -\> a standard output of LTMG parameters

### Select genes

Genes with non-zero expression in more than 5 samples in `data0`

``` r
selected.genes <- which(rowSums(data0 > 0) > 5)
print(head(selected.genes))
```

    ##     RPS11     ELMO2   CREB3L1     PNMA1   TMEM216 LOC653712 
    ##         2         3         4         5         7         8

### Run LTMG for the selected genes

LTMG -\> output list(\(N\): number of peaks; a \(3N\) matrix: \(A\),
\(U\), \(S\); If Zcut (If 0 expression is more than 5, Zcut); Iteration
number, upper limit 1000)

``` r
library(LTMGSCA)
for (gene in head(selected.genes, 3)) {
  for (k in 1:5) {
    print(SeparateKRpkmNew(x = data0[gene, ], n = 100, q = 0, k = k, err = 1e-10))
  }
}
```

    ##      p     mean         sd
    ## [1,] 1 7.802228 0.04416125
    ##                                  p     mean        sd
    ## Late_blastocyst.2_Cell.7 0.5189909 7.632955 0.5688547
    ## X4.cell_embryo.1_Cell.4  0.4810091 7.872066 0.5948461
    ##                                  p     mean        sd
    ## X2.cell_embryo.1_Cell.2  0.1210234 7.083230 0.1913294
    ## Late_blastocyst.3_Cell.5 0.7937210 7.768004 0.5549263
    ## X4.cell_embryo.3_Cell.2  0.0852556 8.505077 0.1026225
    ##                                    p     mean         sd
    ## Late_blastocyst.2_Cell.4 0.006933978 7.788798 0.36486339
    ## Morulae.1_Cell.8         0.448993364 7.267528 0.35678557
    ## Morulae.2_Cell.3         0.480661380 8.095959 0.44529688
    ## X4.cell_embryo.3_Cell.1  0.063411278 8.507573 0.09591338
    ##                                    p     mean        sd
    ## Oocyte.1                 0.005328433 7.727174 0.3884293
    ## Late_blastocyst.2_Cell.7 0.356379577 7.199494 0.3289284
    ## Late_blastocyst.3_Cell.5 0.309216553 7.755904 0.3974516
    ## X4.cell_embryo.1_Cell.4  0.322192465 8.343583 0.3577630
    ## X8.cell_embryo.1_Cell.1  0.006882971 7.925339 0.4211202
    ##        p mean  sd
    ## [1,] NaN  NaN NaN
    ##                                 p      mean       sd
    ## X8.cell_embryo.2_Cell.4 0.3138683 -2.385534 2.460532
    ## X8.cell_embryo.1_Cell.1 0.6861317  2.645441 1.280540
    ##                                  p      mean        sd
    ## Late_blastocyst.2_Cell.1 0.2726636 -2.943190 1.0155455
    ## X8.cell_embryo.3_Cell.1  0.4874510  1.851239 0.8095599
    ## X4.cell_embryo.1_Cell.2  0.2398854  4.123094 0.1710736
    ##                                 p      mean        sd
    ## X8.cell_embryo.2_Cell.6 0.2722165 -2.584755 0.8242646
    ## X8.cell_embryo.2_Cell.7 0.2930781  1.687644 0.7978886
    ## X8.cell_embryo.2_Cell.3 0.1945976  2.090904 0.7665579
    ## X2.cell_embryo.3_Cell.1 0.2401079  4.123058 0.1711315
    ##                                  p      mean        sd
    ## Morulae.1_Cell.2        0.27343650 -2.346844 0.7156853
    ## X8.cell_embryo.2_Cell.4 0.22096460  1.561465 0.7070268
    ## X8.cell_embryo.3_Cell.1 0.16654932  1.775654 0.7409112
    ## X8.cell_embryo.1_Cell.1 0.09679837  2.609410 0.5179847
    ## X4.cell_embryo.2_Cell.4 0.24225121  4.122422 0.1716394
    ##      p           mean        sd
    ## [1,] 1 -1.154028e+138 0.6502746
    ##                                 p       mean        sd
    ## Zygote.1                0.8524664 -1.4311960 1.4985531
    ## X8.cell_embryo.1_Cell.2 0.1475336  0.9262817 0.9096018
    ##                                   p       mean           sd
    ## X4.cell_embryo.2_Cell.3  0.70557942 -1.2865546 0.4785381776
    ## Late_blastocyst.2_Cell.3 0.27219997  0.7239671 0.5938110246
    ## X4.cell_embryo.3_Cell.3  0.02222061  2.7078834 0.0005667612
    ##                                  p       mean           sd
    ## Zygote.2                0.72496505 -1.2936882 0.4602185469
    ## X2.cell_embryo.2_Cell.1 0.11346035  0.4516975 0.2993225110
    ## X8.cell_embryo.1_Cell.4 0.13935362  1.0987104 0.4882202184
    ## X4.cell_embryo.1_Cell.2 0.02222099  2.7078834 0.0005667612
    ##                                   p       mean           sd
    ## Zygote.2                 0.69349422 -1.3076762 0.4388586900
    ## Zygote.1                 0.06286983  0.1033333 0.7265148824
    ## Late_blastocyst.2_Cell.3 0.08060027  0.4484580 0.3466538911
    ## X8.cell_embryo.1_Cell.2  0.14081500  1.0348833 0.5179808674
    ## X4.cell_embryo.3_Cell.4  0.02222069  2.7078834 0.0005667612

### Here we have the BIC functions:

``` r
BIC_f_zcut <- function(y, rrr, Zcut) {
  n <- length(y)
  nparams <- nrow(rrr) * 3
  w <- rrr[, 1]
  u <- rrr[, 2]
  sig <- rrr[, 3]
  cc <- c()
  y0 <- y[which(y >= Zcut)]
  y1 <- y[which(y < Zcut)]
  y1 <- y1 * 0 + Zcut
  for (i in 1:nrow(rrr)) {
    c0 <- dnorm(y0, u[i], sig[i]) * w[i]
    c1 <- (1 - pnorm(y1, u[i], sig[i])) * w[i]
    c <- c(c0, c1)
    cc <- rbind(cc, c)
  }
  d <- apply(cc, 2, sum)
  e <- sum(log(d))
  f <- e * 2 - nparams * log(n)
  return (f)
}

BIC_f_zcut2 <- function(y, rrr, Zcut) {
  n <- length(y)
  nparams <- nrow(rrr) * 3
  w <- rrr[, 1]
  u <- rrr[, 2]
  sig <- rrr[, 3]
  y0 <- y[which(y >= Zcut)]
  cc <- c()
  for (i in 1:nrow(rrr)) {
    c <- dnorm(y0, u[i], sig[i]) * w[i]
    cc <- rbind(cc, c)
  }
  d <- apply(cc, 2, sum)
  e <- sum(log(d))
  f <- e * 2 - nparams * log(n)
  return (f)
}
```

We can now get `f` value using `BIC_f_zcut2()`.

``` r
for (k in 1:5) {
  rrr <- SeparateKRpkmNew(x = data0[selected.genes[1], ], n = 100, q = 0, k = k, err = 1e-10)
  print(BIC_f_zcut2(y = data0[selected.genes[1], ], rrr, 0))
}
```

    ## [1] -16016.58
    ## [1] -188.4811
    ## [1] -195.0944
    ## [1] -208.8539
    ## [1] -223.8346

We are only print while `k != 2`.

``` r
GetBestK <- function(x, n, q, err = 1e-10){
  best.bic <- -Inf
  best.k <- 0
  best.result <- c(0, 0, 0)
  for (k in 1:7) {
    rrr <- SeparateKRpkmNew(x = x, n = n, q = q, k = k, err = err)
    bic <- BIC_f_zcut2(y = x, rrr, q)
    if(is.nan(bic)) {
      bic <- -Inf
    }
    if (bic >= best.bic) {
      best.bic <- bic
      best.k <- k
      best.result <- rrr
    } else {
      return(list(k = best.k, bic = best.bic, result = best.result))
    }
  }
  return(list(k = 0, bic = 0, result = c(0, 0, 0)))
}
```

``` r
for (gene in head(selected.genes, 30)) {
  best <- GetBestK(x = data0[gene, ], n = 100, q = 0, err = 1e-10)
  if (best[1] != 2) {
    print(gene)
  }
}
```

    ## [1] 3
    ## [1] 4
    ## [1] 42

This is the 3rd one:

``` r
best <- GetBestK(x = data0[3, ], n = 100, q = 0, err = 1e-10)
print(best)
```

    ## $k
    ## [1] 3
    ## 
    ## $bic
    ## [1] -247.7207
    ## 
    ## $result
    ##                                  p      mean        sd
    ## Late_blastocyst.2_Cell.1 0.2726636 -2.943190 1.0155455
    ## X8.cell_embryo.3_Cell.1  0.4874510  1.851239 0.8095599
    ## X4.cell_embryo.1_Cell.2  0.2398854  4.123094 0.1710736

``` r
hist(data0[3,], breaks = 60)
```

![](C:/Users/yz116/AppData/Local/Temp/RtmpKqIkat/preview-429c2ea323b2.dir/gcr_vignette_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

## (iii) LTMG -\> discretization

We have the following functions ready for this step:
`calculate_prob_sep_Zcut`, `discretization_method_1_LLR_mean`, and
`Build_R_matrix`.

``` r
calculate_prob_sep_Zcut <- function(data1, Zcut, a, u, sig) {
  cc <- matrix(0, length(a), length(data1))
  colnames(cc) <- names(data1)
  for (i in 1:length(a)) {
    c <- a[i] / sig[i] * exp(-(data1 - u[i]) ^ 2 / (2 * sig[i] ^ 2))
    cc[i, ] <- c
  }
  cut_p <- rep(0, length(a))
  for (i in 1:length(a)) {
    cut_p[i] <- a[i] * pnorm(Zcut, u[i], sig[i])
  }
  for (i in 1:ncol(cc)) {
    if (data1[i] < Zcut) {
      cc[, i] <- cut_p
    }
  }
  cc[which(is.na(cc) == 1)] <- 0
  return(cc)
}
```

``` r
discretization_method_1_LLR_mean <- function(y, aaa, ccc, LLR_cut = 2) {
  K <- 1 / LLR_cut + 1
  if (nrow(aaa) == 1) {
    print("Only one class")
    return(y)
  } else {
    discretized_y <- rep(0, length(y))
    for (i in 1:ncol(ccc)) {
      ll <- which(ccc[, i] == max(ccc[, i]))[1]
      if ((max(ccc[, i])/sum(ccc[, i])) > (1/K)) {
        discretized_y[i] <- ll
      }
    }
    blocks <- c()
    st_c <- 1
    end_c <- 1
    st_c_v <- y[order(y)[1]]
    end_c_v <- y[order(y)[1]]
    label_c <- discretized_y[order(y)[1]]
    for (i in 2:length(order(y))) {
      if (discretized_y[order(y)[i]] == discretized_y[order(y)[i - 1]]) {
        end_c <- i
        end_c_v <- y[order(y)[i]]
        if (i == length(order(y))) {
          end_c <- i
          end_c_v <- y[order(y)[i]]
          blocks <- rbind(blocks, c(st_c, end_c, st_c_v, end_c_v, label_c))
        }
      } else {
        blocks <- rbind(blocks, c(st_c, end_c, st_c_v, end_c_v, 
          label_c))
        label_c <- discretized_y[order(y)[i]]
        st_c <- i
        end_c <- i
        st_c_v <- y[order(y)[i]]
        end_c_v <- y[order(y)[i]]
        if (i == length(order(y))) {
          end_c <- i
          end_c_v <- y[order(y)[i]]
          blocks <- rbind(blocks, c(st_c, end_c, st_c_v, end_c_v, label_c))
        }
      }
    }
    if (nrow(blocks) > 1) {
      for (i in 1:nrow(blocks)) {
        if (blocks[i, 5] != 0) {
          tg_i <- blocks[i, 5]
          if (!((blocks[i, 3] <= aaa[tg_i, 2]) & (blocks[i, 4] >= aaa[tg_i, 2]))) {
          blocks[i, 5] <- 0
          }
        }
      }
      for (i in 1:nrow(blocks)) {
        discretized_y[order(y)[blocks[i, 1]:blocks[i, 2]]] <- blocks[i, 5]
      }
    }
    return(discretized_y)
  }
}
```

``` r
Build_R_matrix <- function(cc, Zcut0, U, Gname) {
  tg_s <- intersect(which(U > Zcut0), unique(cc))
  dd <- c()
  nc <- c()
  if (length(tg_s) > 0) {
    for (i in 1:length(tg_s)) {
      nc <- c(nc, paste(Gname, tg_s[i], sep = "__"))
      ccc <- (cc == tg_s[i]) * 1
      dd <- rbind(dd, ccc)
    }
  }
  rownames(dd) <- nc
  return(dd)
}
```

best$result is a K\*3 matrix with 1st, 2nd and 3rd columns are the A, U,
S of the gene x is the normalized expression level

``` r
i <- 4
x <- data0[i, ]
Zcut0 <- 0
best <- GetBestK(x = x, n = 1000, q = Zcut0, err = 1e-10)

pp <- calculate_prob_sep_Zcut(x, Zcut0, best$result[, 1], best$result[, 2], best$result[, 3])
cc <- discretization_method_1_LLR_mean(x, best$result, pp, LLR_cut = 0.1)
dd <- Build_R_matrix(cc, Zcut0, best$result[, 2], rownames(data0)[i])

print(x)
```

    ##                  Oocyte.1                  Oocyte.2 
    ##                 0.0000000                 0.3534698 
    ##                  Oocyte.3                  Zygote.1 
    ##                 0.6657760                 0.5905606 
    ##                  Zygote.2                  Zygote.3 
    ##                 0.3611648                -0.3133418 
    ##   X2.cell_embryo.1_Cell.1   X2.cell_embryo.1_Cell.2 
    ##                 0.4643627                 0.2342813 
    ##   X2.cell_embryo.2_Cell.1   X2.cell_embryo.2_Cell.2 
    ##                 0.6339278                 0.3708737 
    ##   X2.cell_embryo.3_Cell.1   X2.cell_embryo.3_Cell.2 
    ##                -0.3523984                -0.7052198 
    ##   X4.cell_embryo.1_Cell.1   X4.cell_embryo.1_Cell.2 
    ##                      -Inf                 1.5526559 
    ##   X4.cell_embryo.1_Cell.3   X4.cell_embryo.1_Cell.4 
    ##                      -Inf                      -Inf 
    ##   X4.cell_embryo.2_Cell.1   X4.cell_embryo.2_Cell.2 
    ##                -1.3586792                -0.2943711 
    ##   X4.cell_embryo.2_Cell.3   X4.cell_embryo.2_Cell.4 
    ##                 0.4491630                 1.0217312 
    ##   X4.cell_embryo.3_Cell.1   X4.cell_embryo.3_Cell.2 
    ##                 1.1177611                 0.9250522 
    ##   X4.cell_embryo.3_Cell.3   X4.cell_embryo.3_Cell.4 
    ##                 1.4548872                 1.6122340 
    ##   X8.cell_embryo.1_Cell.1   X8.cell_embryo.1_Cell.2 
    ##                -0.7635696                 1.0328285 
    ##   X8.cell_embryo.1_Cell.3   X8.cell_embryo.1_Cell.4 
    ##                 0.6559644                 0.9266370 
    ##   X8.cell_embryo.2_Cell.1   X8.cell_embryo.2_Cell.2 
    ##                      -Inf                      -Inf 
    ##   X8.cell_embryo.2_Cell.3   X8.cell_embryo.2_Cell.4 
    ##                      -Inf                      -Inf 
    ##   X8.cell_embryo.2_Cell.5   X8.cell_embryo.2_Cell.6 
    ##                      -Inf                      -Inf 
    ##   X8.cell_embryo.2_Cell.7   X8.cell_embryo.2_Cell.8 
    ##                      -Inf                      -Inf 
    ##   X8.cell_embryo.3_Cell.1   X8.cell_embryo.3_Cell.2 
    ##                      -Inf                      -Inf 
    ##   X8.cell_embryo.3_Cell.3   X8.cell_embryo.3_Cell.4 
    ##                      -Inf                      -Inf 
    ##   X8.cell_embryo.3_Cell.5   X8.cell_embryo.3_Cell.6 
    ##                      -Inf                      -Inf 
    ##   X8.cell_embryo.3_Cell.7   X8.cell_embryo.3_Cell.8 
    ##                      -Inf                 1.7516317 
    ##          Morulae.1_Cell.1          Morulae.1_Cell.2 
    ##                      -Inf                      -Inf 
    ##          Morulae.1_Cell.3          Morulae.1_Cell.4 
    ##                      -Inf                -0.2930297 
    ##          Morulae.1_Cell.5          Morulae.1_Cell.6 
    ##                      -Inf                      -Inf 
    ##          Morulae.1_Cell.7          Morulae.1_Cell.8 
    ##                      -Inf                      -Inf 
    ##          Morulae.2_Cell.1          Morulae.2_Cell.2 
    ##                      -Inf                      -Inf 
    ##          Morulae.2_Cell.3          Morulae.2_Cell.4 
    ##                      -Inf                      -Inf 
    ##          Morulae.2_Cell.5          Morulae.2_Cell.6 
    ##                      -Inf                -0.9702191 
    ##          Morulae.2_Cell.7          Morulae.2_Cell.8 
    ##                -0.9597203                      -Inf 
    ##  Late_blastocyst.1_Cell.1  Late_blastocyst.1_Cell.2 
    ##                      -Inf                      -Inf 
    ##  Late_blastocyst.1_Cell.3  Late_blastocyst.1_Cell.4 
    ##                      -Inf                      -Inf 
    ##  Late_blastocyst.1_Cell.5  Late_blastocyst.1_Cell.6 
    ##                      -Inf                      -Inf 
    ##  Late_blastocyst.1_Cell.7  Late_blastocyst.1_Cell.8 
    ##                      -Inf                      -Inf 
    ##  Late_blastocyst.1_Cell.9 Late_blastocyst.1_Cell.10 
    ##                      -Inf                      -Inf 
    ## Late_blastocyst.1_Cell.11 Late_blastocyst.1_Cell.12 
    ##                      -Inf                      -Inf 
    ##  Late_blastocyst.2_Cell.1  Late_blastocyst.2_Cell.2 
    ##                -0.3495575                 1.6981811 
    ##  Late_blastocyst.2_Cell.3  Late_blastocyst.2_Cell.4 
    ##                 0.6714127                      -Inf 
    ##  Late_blastocyst.2_Cell.5  Late_blastocyst.2_Cell.6 
    ##                      -Inf                      -Inf 
    ##  Late_blastocyst.2_Cell.7  Late_blastocyst.2_Cell.8 
    ##                      -Inf                      -Inf 
    ##  Late_blastocyst.2_Cell.9 Late_blastocyst.2_Cell.10 
    ##                      -Inf                      -Inf 
    ##  Late_blastocyst.3_Cell.1  Late_blastocyst.3_Cell.2 
    ##                      -Inf                      -Inf 
    ##  Late_blastocyst.3_Cell.3  Late_blastocyst.3_Cell.4 
    ##                      -Inf                      -Inf 
    ##  Late_blastocyst.3_Cell.5  Late_blastocyst.3_Cell.6 
    ##                 2.7084501                      -Inf 
    ##  Late_blastocyst.3_Cell.7  Late_blastocyst.3_Cell.8 
    ##                      -Inf                 2.7073166

``` r
print(pp)
```

    ##        Oocyte.1    Oocyte.2     Oocyte.3     Zygote.1
    ## [1,] 0.02453824 0.001542643 0.0000750863 0.0001633965
    ## [2,] 0.22504345 0.382762784 0.4578324387 0.4495556219
    ## [3,] 0.00000000 0.000000000 0.0000000000 0.0000000000
    ##        Zygote.2   Zygote.3 X2.cell_embryo.1_Cell.1
    ## [1,] 0.00144129 0.70151041            0.0005612054
    ## [2,] 0.38571522 0.03192687            0.4208080150
    ## [3,] 0.00000000 0.00000000            0.0000000000
    ##      X2.cell_embryo.1_Cell.2 X2.cell_embryo.2_Cell.1
    ## [1,]             0.004237879            0.0001047639
    ## [2,]             0.332717460            0.4551859299
    ## [3,]             0.000000000            0.0000000000
    ##      X2.cell_embryo.2_Cell.2 X2.cell_embryo.3_Cell.1
    ## [1,]             0.001322238              0.70151041
    ## [2,]             0.389380992              0.03192687
    ## [3,]             0.000000000              0.00000000
    ##      X2.cell_embryo.3_Cell.2 X4.cell_embryo.1_Cell.1
    ## [1,]              0.70151041              0.70151041
    ## [2,]              0.03192687              0.03192687
    ## [3,]              0.00000000              0.00000000
    ##      X4.cell_embryo.1_Cell.2 X4.cell_embryo.1_Cell.3
    ## [1,]            7.305006e-10              0.70151041
    ## [2,]            1.725700e-01              0.03192687
    ## [3,]            0.000000e+00              0.00000000
    ##      X4.cell_embryo.1_Cell.4 X4.cell_embryo.2_Cell.1
    ## [1,]              0.70151041              0.70151041
    ## [2,]              0.03192687              0.03192687
    ## [3,]              0.00000000              0.00000000
    ##      X4.cell_embryo.2_Cell.2 X4.cell_embryo.2_Cell.3
    ## [1,]              0.70151041            0.0006472473
    ## [2,]              0.03192687            0.4162216657
    ## [3,]              0.00000000            0.0000000000
    ##      X4.cell_embryo.2_Cell.4 X4.cell_embryo.3_Cell.1
    ## [1,]            1.236439e-06            3.619440e-07
    ## [2,]            4.028836e-01            3.663466e-01
    ## [3,]            0.000000e+00            0.000000e+00
    ##      X4.cell_embryo.3_Cell.2 X4.cell_embryo.3_Cell.3
    ## [1,]            4.044036e-06            3.230622e-09
    ## [2,]            4.319703e-01            2.140125e-01
    ## [3,]            0.000000e+00            0.000000e+00
    ##      X4.cell_embryo.3_Cell.4 X8.cell_embryo.1_Cell.1
    ## [1,]            2.876387e-10              0.70151041
    ## [2,]            1.493903e-01              0.03192687
    ## [3,]            0.000000e+00              0.00000000
    ##      X8.cell_embryo.1_Cell.2 X8.cell_embryo.1_Cell.3
    ## [1,]            1.075611e-06            8.325015e-05
    ## [2,]            3.990062e-01            4.571534e-01
    ## [3,]            0.000000e+00            0.000000e+00
    ##      X8.cell_embryo.1_Cell.4 X8.cell_embryo.2_Cell.1
    ## [1,]            3.967901e-06              0.70151041
    ## [2,]            4.315677e-01              0.03192687
    ## [3,]            0.000000e+00              0.00000000
    ##      X8.cell_embryo.2_Cell.2 X8.cell_embryo.2_Cell.3
    ## [1,]              0.70151041              0.70151041
    ## [2,]              0.03192687              0.03192687
    ## [3,]              0.00000000              0.00000000
    ##      X8.cell_embryo.2_Cell.4 X8.cell_embryo.2_Cell.5
    ## [1,]              0.70151041              0.70151041
    ## [2,]              0.03192687              0.03192687
    ## [3,]              0.00000000              0.00000000
    ##      X8.cell_embryo.2_Cell.6 X8.cell_embryo.2_Cell.7
    ## [1,]              0.70151041              0.70151041
    ## [2,]              0.03192687              0.03192687
    ## [3,]              0.00000000              0.00000000
    ##      X8.cell_embryo.2_Cell.8 X8.cell_embryo.3_Cell.1
    ## [1,]              0.70151041              0.70151041
    ## [2,]              0.03192687              0.03192687
    ## [3,]              0.00000000              0.00000000
    ##      X8.cell_embryo.3_Cell.2 X8.cell_embryo.3_Cell.3
    ## [1,]              0.70151041              0.70151041
    ## [2,]              0.03192687              0.03192687
    ## [3,]              0.00000000              0.00000000
    ##      X8.cell_embryo.3_Cell.4 X8.cell_embryo.3_Cell.5
    ## [1,]              0.70151041              0.70151041
    ## [2,]              0.03192687              0.03192687
    ## [3,]              0.00000000              0.00000000
    ##      X8.cell_embryo.3_Cell.6 X8.cell_embryo.3_Cell.7
    ## [1,]              0.70151041              0.70151041
    ## [2,]              0.03192687              0.03192687
    ## [3,]              0.00000000              0.00000000
    ##      X8.cell_embryo.3_Cell.8 Morulae.1_Cell.1
    ## [1,]            3.008084e-11       0.70151041
    ## [2,]            1.025517e-01       0.03192687
    ## [3,]            0.000000e+00       0.00000000
    ##      Morulae.1_Cell.2 Morulae.1_Cell.3 Morulae.1_Cell.4
    ## [1,]       0.70151041       0.70151041       0.70151041
    ## [2,]       0.03192687       0.03192687       0.03192687
    ## [3,]       0.00000000       0.00000000       0.00000000
    ##      Morulae.1_Cell.5 Morulae.1_Cell.6 Morulae.1_Cell.7
    ## [1,]       0.70151041       0.70151041       0.70151041
    ## [2,]       0.03192687       0.03192687       0.03192687
    ## [3,]       0.00000000       0.00000000       0.00000000
    ##      Morulae.1_Cell.8 Morulae.2_Cell.1 Morulae.2_Cell.2
    ## [1,]       0.70151041       0.70151041       0.70151041
    ## [2,]       0.03192687       0.03192687       0.03192687
    ## [3,]       0.00000000       0.00000000       0.00000000
    ##      Morulae.2_Cell.3 Morulae.2_Cell.4 Morulae.2_Cell.5
    ## [1,]       0.70151041       0.70151041       0.70151041
    ## [2,]       0.03192687       0.03192687       0.03192687
    ## [3,]       0.00000000       0.00000000       0.00000000
    ##      Morulae.2_Cell.6 Morulae.2_Cell.7 Morulae.2_Cell.8
    ## [1,]       0.70151041       0.70151041       0.70151041
    ## [2,]       0.03192687       0.03192687       0.03192687
    ## [3,]       0.00000000       0.00000000       0.00000000
    ##      Late_blastocyst.1_Cell.1 Late_blastocyst.1_Cell.2
    ## [1,]               0.70151041               0.70151041
    ## [2,]               0.03192687               0.03192687
    ## [3,]               0.00000000               0.00000000
    ##      Late_blastocyst.1_Cell.3 Late_blastocyst.1_Cell.4
    ## [1,]               0.70151041               0.70151041
    ## [2,]               0.03192687               0.03192687
    ## [3,]               0.00000000               0.00000000
    ##      Late_blastocyst.1_Cell.5 Late_blastocyst.1_Cell.6
    ## [1,]               0.70151041               0.70151041
    ## [2,]               0.03192687               0.03192687
    ## [3,]               0.00000000               0.00000000
    ##      Late_blastocyst.1_Cell.7 Late_blastocyst.1_Cell.8
    ## [1,]               0.70151041               0.70151041
    ## [2,]               0.03192687               0.03192687
    ## [3,]               0.00000000               0.00000000
    ##      Late_blastocyst.1_Cell.9 Late_blastocyst.1_Cell.10
    ## [1,]               0.70151041                0.70151041
    ## [2,]               0.03192687                0.03192687
    ## [3,]               0.00000000                0.00000000
    ##      Late_blastocyst.1_Cell.11 Late_blastocyst.1_Cell.12
    ## [1,]                0.70151041                0.70151041
    ## [2,]                0.03192687                0.03192687
    ## [3,]                0.00000000                0.00000000
    ##      Late_blastocyst.2_Cell.1 Late_blastocyst.2_Cell.2
    ## [1,]               0.70151041             7.241376e-11
    ## [2,]               0.03192687             1.192267e-01
    ## [3,]               0.00000000             0.000000e+00
    ##      Late_blastocyst.2_Cell.3 Late_blastocyst.2_Cell.4
    ## [1,]             7.074638e-05               0.70151041
    ## [2,]             4.581673e-01               0.03192687
    ## [3,]             0.000000e+00               0.00000000
    ##      Late_blastocyst.2_Cell.5 Late_blastocyst.2_Cell.6
    ## [1,]               0.70151041               0.70151041
    ## [2,]               0.03192687               0.03192687
    ## [3,]               0.00000000               0.00000000
    ##      Late_blastocyst.2_Cell.7 Late_blastocyst.2_Cell.8
    ## [1,]               0.70151041               0.70151041
    ## [2,]               0.03192687               0.03192687
    ## [3,]               0.00000000               0.00000000
    ##      Late_blastocyst.2_Cell.9 Late_blastocyst.2_Cell.10
    ## [1,]               0.70151041                0.70151041
    ## [2,]               0.03192687                0.03192687
    ## [3,]               0.00000000                0.00000000
    ##      Late_blastocyst.3_Cell.1 Late_blastocyst.3_Cell.2
    ## [1,]               0.70151041               0.70151041
    ## [2,]               0.03192687               0.03192687
    ## [3,]               0.00000000               0.00000000
    ##      Late_blastocyst.3_Cell.3 Late_blastocyst.3_Cell.4
    ## [1,]               0.70151041               0.70151041
    ## [2,]               0.03192687               0.03192687
    ## [3,]               0.00000000               0.00000000
    ##      Late_blastocyst.3_Cell.5 Late_blastocyst.3_Cell.6
    ## [1,]             3.029268e-19               0.70151041
    ## [2,]             1.793851e-03               0.03192687
    ## [3,]             2.377976e+01               0.00000000
    ##      Late_blastocyst.3_Cell.7 Late_blastocyst.3_Cell.8
    ## [1,]               0.70151041             3.105414e-19
    ## [2,]               0.03192687             1.805197e-03
    ## [3,]               0.00000000             2.377974e+01

``` r
print(cc)
```

    ##  [1] 2 2 2 2 2 1 2 2 2 2 1 1 1 2 1 1 1 1 2 2 2 2 2 2 1 2 2 2
    ## [29] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1
    ## [57] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 1 1 1
    ## [85] 1 1 3 1 1 3

``` r
print(dd)
```

    ##            [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
    ## CREB3L1__2    1    1    1    1    1    0    1    1    1
    ## CREB3L1__3    0    0    0    0    0    0    0    0    0
    ##            [,10] [,11] [,12] [,13] [,14] [,15] [,16] [,17]
    ## CREB3L1__2     1     0     0     0     1     0     0     0
    ## CREB3L1__3     0     0     0     0     0     0     0     0
    ##            [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25]
    ## CREB3L1__2     0     1     1     1     1     1     1     0
    ## CREB3L1__3     0     0     0     0     0     0     0     0
    ##            [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33]
    ## CREB3L1__2     1     1     1     0     0     0     0     0
    ## CREB3L1__3     0     0     0     0     0     0     0     0
    ##            [,34] [,35] [,36] [,37] [,38] [,39] [,40] [,41]
    ## CREB3L1__2     0     0     0     0     0     0     0     0
    ## CREB3L1__3     0     0     0     0     0     0     0     0
    ##            [,42] [,43] [,44] [,45] [,46] [,47] [,48] [,49]
    ## CREB3L1__2     0     0     1     0     0     0     0     0
    ## CREB3L1__3     0     0     0     0     0     0     0     0
    ##            [,50] [,51] [,52] [,53] [,54] [,55] [,56] [,57]
    ## CREB3L1__2     0     0     0     0     0     0     0     0
    ## CREB3L1__3     0     0     0     0     0     0     0     0
    ##            [,58] [,59] [,60] [,61] [,62] [,63] [,64] [,65]
    ## CREB3L1__2     0     0     0     0     0     0     0     0
    ## CREB3L1__3     0     0     0     0     0     0     0     0
    ##            [,66] [,67] [,68] [,69] [,70] [,71] [,72] [,73]
    ## CREB3L1__2     0     0     0     0     0     0     0     0
    ## CREB3L1__3     0     0     0     0     0     0     0     0
    ##            [,74] [,75] [,76] [,77] [,78] [,79] [,80] [,81]
    ## CREB3L1__2     1     1     0     0     0     0     0     0
    ## CREB3L1__3     0     0     0     0     0     0     0     0
    ##            [,82] [,83] [,84] [,85] [,86] [,87] [,88] [,89]
    ## CREB3L1__2     0     0     0     0     0     0     0     0
    ## CREB3L1__3     0     0     0     0     0     1     0     0
    ##            [,90]
    ## CREB3L1__2     0
    ## CREB3L1__3     1

``` r
i <- 5
x <- data0[i, ]
Zcut0 <- 0
best <- GetBestK(x = x, n = 1000, q = Zcut0, err = 1e-10)

pp <- calculate_prob_sep_Zcut(x, Zcut0, best$result[, 1], best$result[, 2], best$result[, 3])
cc <- discretization_method_1_LLR_mean(x, best$result, pp, LLR_cut = 0.1)
dd <- Build_R_matrix(cc, Zcut0, best$result[, 2], rownames(data0)[i])

print(x)
```

    ##                  Oocyte.1                  Oocyte.2 
    ##                -0.3871342                 0.2949059 
    ##                  Oocyte.3                  Zygote.1 
    ##                 0.7537718                 1.7446679 
    ##                  Zygote.2                  Zygote.3 
    ##                 1.5871923                 1.7313015 
    ##   X2.cell_embryo.1_Cell.1   X2.cell_embryo.1_Cell.2 
    ##                 1.4611699                 1.4011830 
    ##   X2.cell_embryo.2_Cell.1   X2.cell_embryo.2_Cell.2 
    ##                 1.4731599                 1.5475625 
    ##   X2.cell_embryo.3_Cell.1   X2.cell_embryo.3_Cell.2 
    ##                 1.0501221                 0.8742180 
    ##   X4.cell_embryo.1_Cell.1   X4.cell_embryo.1_Cell.2 
    ##                 1.4060970                -1.8325815 
    ##   X4.cell_embryo.1_Cell.3   X4.cell_embryo.1_Cell.4 
    ##                 0.1856493                      -Inf 
    ##   X4.cell_embryo.2_Cell.1   X4.cell_embryo.2_Cell.2 
    ##                 1.1413524                      -Inf 
    ##   X4.cell_embryo.2_Cell.3   X4.cell_embryo.2_Cell.4 
    ##                 2.3799164                -0.9288695 
    ##   X4.cell_embryo.3_Cell.1   X4.cell_embryo.3_Cell.2 
    ##                 1.5411591                 1.7523254 
    ##   X4.cell_embryo.3_Cell.3   X4.cell_embryo.3_Cell.4 
    ##                 1.5091755                 1.7288197 
    ##   X8.cell_embryo.1_Cell.1   X8.cell_embryo.1_Cell.2 
    ##                 3.5036336                 3.6047095 
    ##   X8.cell_embryo.1_Cell.3   X8.cell_embryo.1_Cell.4 
    ##                 2.8523237                 3.7815954 
    ##   X8.cell_embryo.2_Cell.1   X8.cell_embryo.2_Cell.2 
    ##                 2.7822296                 1.7838953 
    ##   X8.cell_embryo.2_Cell.3   X8.cell_embryo.2_Cell.4 
    ##                 2.0550208                 3.6057967 
    ##   X8.cell_embryo.2_Cell.5   X8.cell_embryo.2_Cell.6 
    ##                 2.8845213                 2.2490788 
    ##   X8.cell_embryo.2_Cell.7   X8.cell_embryo.2_Cell.8 
    ##                 2.2758304                 2.9764475 
    ##   X8.cell_embryo.3_Cell.1   X8.cell_embryo.3_Cell.2 
    ##                 3.2518463                 3.4940196 
    ##   X8.cell_embryo.3_Cell.3   X8.cell_embryo.3_Cell.4 
    ##                 4.0421210                 3.8096580 
    ##   X8.cell_embryo.3_Cell.5   X8.cell_embryo.3_Cell.6 
    ##                 3.6084554                 2.3880286 
    ##   X8.cell_embryo.3_Cell.7   X8.cell_embryo.3_Cell.8 
    ##                 2.9516759                 4.0716217 
    ##          Morulae.1_Cell.1          Morulae.1_Cell.2 
    ##                 3.9599366                 3.4106195 
    ##          Morulae.1_Cell.3          Morulae.1_Cell.4 
    ##                -0.1086994                 3.2040464 
    ##          Morulae.1_Cell.5          Morulae.1_Cell.6 
    ##                 3.0037004                 4.1114475 
    ##          Morulae.1_Cell.7          Morulae.1_Cell.8 
    ##                 3.7993021                -0.8141855 
    ##          Morulae.2_Cell.1          Morulae.2_Cell.2 
    ##                 3.0061775                 4.1204670 
    ##          Morulae.2_Cell.3          Morulae.2_Cell.4 
    ##                 3.9427656                 3.4452781 
    ##          Morulae.2_Cell.5          Morulae.2_Cell.6 
    ##                 3.0906333                 3.1098642 
    ##          Morulae.2_Cell.7          Morulae.2_Cell.8 
    ##                 3.8282935                 4.1413236 
    ##  Late_blastocyst.1_Cell.1  Late_blastocyst.1_Cell.2 
    ##                 0.3133498                 2.3487050 
    ##  Late_blastocyst.1_Cell.3  Late_blastocyst.1_Cell.4 
    ##                 3.4520489                 2.5867861 
    ##  Late_blastocyst.1_Cell.5  Late_blastocyst.1_Cell.6 
    ##                 4.2983325                 3.2938349 
    ##  Late_blastocyst.1_Cell.7  Late_blastocyst.1_Cell.8 
    ##                 2.2683042                 1.4092782 
    ##  Late_blastocyst.1_Cell.9 Late_blastocyst.1_Cell.10 
    ##                      -Inf                 3.0728322 
    ## Late_blastocyst.1_Cell.11 Late_blastocyst.1_Cell.12 
    ##                 1.8878269                 2.3270826 
    ##  Late_blastocyst.2_Cell.1  Late_blastocyst.2_Cell.2 
    ##                 2.2086040                 2.9670756 
    ##  Late_blastocyst.2_Cell.3  Late_blastocyst.2_Cell.4 
    ##                 3.2008340                -0.3768777 
    ##  Late_blastocyst.2_Cell.5  Late_blastocyst.2_Cell.6 
    ##                -0.2256467                 2.8136107 
    ##  Late_blastocyst.2_Cell.7  Late_blastocyst.2_Cell.8 
    ##                -0.8462984                 2.6758713 
    ##  Late_blastocyst.2_Cell.9 Late_blastocyst.2_Cell.10 
    ##                 1.4768204                 2.1813210 
    ##  Late_blastocyst.3_Cell.1  Late_blastocyst.3_Cell.2 
    ##                 2.0725428                 1.3790180 
    ##  Late_blastocyst.3_Cell.3  Late_blastocyst.3_Cell.4 
    ##                 2.3051817                 2.4216118 
    ##  Late_blastocyst.3_Cell.5  Late_blastocyst.3_Cell.6 
    ##                 3.4375293                -1.9661129 
    ##  Late_blastocyst.3_Cell.7  Late_blastocyst.3_Cell.8 
    ##                -0.9314044                 1.9226415

``` r
print(pp)
```

    ##         Oocyte.1   Oocyte.2   Oocyte.3   Zygote.1
    ## [1,] 0.142058828 0.10099843 0.08286467 0.04255067
    ## [2,] 0.002148516 0.03887417 0.11327751 0.52103074
    ##        Zygote.2   Zygote.3 X2.cell_embryo.1_Cell.1
    ## [1,] 0.04835016 0.04302846              0.05323813
    ## [2,] 0.43916973 0.51406614              0.37563368
    ##      X2.cell_embryo.1_Cell.2 X2.cell_embryo.2_Cell.1
    ## [1,]                0.055632              0.05276456
    ## [2,]                0.346590              0.38154475
    ##      X2.cell_embryo.2_Cell.2 X2.cell_embryo.3_Cell.1
    ## [1,]              0.04986525              0.07025561
    ## [2,]              0.41889125              0.20003396
    ##      X2.cell_embryo.3_Cell.2 X4.cell_embryo.1_Cell.1
    ## [1,]              0.07776154              0.05543441
    ## [2,]              0.14438931              0.34893383
    ##      X4.cell_embryo.1_Cell.2 X4.cell_embryo.1_Cell.3
    ## [1,]             0.142058828              0.10478288
    ## [2,]             0.002148516              0.02913153
    ##      X4.cell_embryo.1_Cell.4 X4.cell_embryo.2_Cell.1
    ## [1,]             0.142058828              0.06638278
    ## [2,]             0.002148516              0.23375092
    ##      X4.cell_embryo.2_Cell.2 X4.cell_embryo.2_Cell.3
    ## [1,]             0.142058828              0.02336931
    ## [2,]             0.002148516              0.78888868
    ##      X4.cell_embryo.2_Cell.4 X4.cell_embryo.3_Cell.1
    ## [1,]             0.142058828              0.05011202
    ## [2,]             0.002148516              0.41563682
    ##      X4.cell_embryo.3_Cell.2 X4.cell_embryo.3_Cell.3
    ## [1,]              0.04227821              0.05135246
    ## [2,]              0.52501705              0.39948877
    ##      X4.cell_embryo.3_Cell.4 X8.cell_embryo.1_Cell.1
    ## [1,]              0.04311748             0.005825674
    ## [2,]              0.51277225             0.559069231
    ##      X8.cell_embryo.1_Cell.2 X8.cell_embryo.1_Cell.3
    ## [1,]              0.00503651              0.01371725
    ## [2,]              0.50662638              0.80729450
    ##      X8.cell_embryo.1_Cell.4 X8.cell_embryo.2_Cell.1
    ## [1,]             0.003872077              0.01491548
    ## [2,]             0.415124284              0.81700745
    ##      X8.cell_embryo.2_Cell.2 X8.cell_embryo.2_Cell.3
    ## [1,]              0.04116475              0.03228834
    ## [2,]              0.54141032              0.67416070
    ##      X8.cell_embryo.2_Cell.4 X8.cell_embryo.2_Cell.5
    ## [1,]             0.005028538              0.01319236
    ## [2,]             0.506059230              0.80143059
    ##      X8.cell_embryo.2_Cell.6 X8.cell_embryo.2_Cell.7
    ## [1,]              0.02673139               0.0260188
    ## [2,]              0.75083533               0.7596164
    ##      X8.cell_embryo.2_Cell.8 X8.cell_embryo.3_Cell.1
    ## [1,]               0.0117793             0.008248957
    ## [2,]               0.7800521             0.680740316
    ##      X8.cell_embryo.3_Cell.2 X8.cell_embryo.3_Cell.3
    ## [1,]             0.005905847             0.002579446
    ## [2,]             0.564004259             0.290916028
    ##      X8.cell_embryo.3_Cell.4 X8.cell_embryo.3_Cell.5
    ## [1,]             0.003710332             0.005009086
    ## [2,]             0.400950808             0.504672239
    ##      X8.cell_embryo.3_Cell.6 X8.cell_embryo.3_Cell.7
    ## [1,]              0.02317101              0.01214783
    ## [2,]              0.79082419              0.78646959
    ##      X8.cell_embryo.3_Cell.8 Morulae.1_Cell.1
    ## [1,]             0.002459978      0.002939261
    ## [2,]             0.278135320      0.328057092
    ##      Morulae.1_Cell.2 Morulae.1_Cell.3 Morulae.1_Cell.4
    ## [1,]      0.006640721      0.142058828      0.008791032
    ## [2,]      0.606114778      0.002148516      0.701170601
    ##      Morulae.1_Cell.5 Morulae.1_Cell.6 Morulae.1_Cell.7
    ## [1,]       0.01138407      0.002306365      0.003769333
    ## [2,]       0.77245475      0.261373546      0.406164587
    ##      Morulae.1_Cell.8 Morulae.2_Cell.1 Morulae.2_Cell.2
    ## [1,]      0.142058828       0.01134867      0.002272763
    ## [2,]      0.002148516       0.77173688      0.257658019
    ##      Morulae.2_Cell.3 Morulae.2_Cell.4 Morulae.2_Cell.5
    ## [1,]        0.0030197      0.006326611       0.01019296
    ## [2,]        0.3360840      0.588788281       0.74466823
    ##      Morulae.2_Cell.6 Morulae.2_Cell.7 Morulae.2_Cell.8
    ## [1,]        0.0099434      0.003606151      0.002196698
    ## [2,]        0.7378356      0.391621125      0.249182691
    ##      Late_blastocyst.1_Cell.1 Late_blastocyst.1_Cell.2
    ## [1,]               0.10033384               0.02414328
    ## [2,]               0.04076196               0.78096288
    ##      Late_blastocyst.1_Cell.3 Late_blastocyst.1_Cell.4
    ## [1,]              0.006266712               0.01867651
    ## [2,]              0.585372106               0.82113503
    ##      Late_blastocyst.1_Cell.5 Late_blastocyst.1_Cell.6
    ## [1,]              0.001692224              0.007795538
    ## [2,]              0.190788203              0.661922969
    ##      Late_blastocyst.1_Cell.7 Late_blastocyst.1_Cell.8
    ## [1,]               0.02621797               0.05530664
    ## [2,]               0.75719541               0.35045467
    ##      Late_blastocyst.1_Cell.9 Late_blastocyst.1_Cell.10
    ## [1,]              0.142058828                 0.0104284
    ## [2,]              0.002148516                 0.7507792
    ##      Late_blastocyst.1_Cell.11 Late_blastocyst.1_Cell.12
    ## [1,]                0.03761386                0.02468975
    ## [2,]                0.59449805                0.77503573
    ##      Late_blastocyst.2_Cell.1 Late_blastocyst.2_Cell.2
    ## [1,]               0.02783421               0.01191768
    ## [2,]               0.73664808               0.78253551
    ##      Late_blastocyst.2_Cell.3 Late_blastocyst.2_Cell.4
    ## [1,]              0.008828473              0.142058828
    ## [2,]              0.702502598              0.002148516
    ##      Late_blastocyst.2_Cell.5 Late_blastocyst.2_Cell.6
    ## [1,]              0.142058828               0.01436952
    ## [2,]              0.002148516               0.81318328
    ##      Late_blastocyst.2_Cell.7 Late_blastocyst.2_Cell.8
    ## [1,]              0.142058828               0.01688368
    ## [2,]              0.002148516               0.82349604
    ##      Late_blastocyst.2_Cell.9 Late_blastocyst.2_Cell.10
    ## [1,]               0.05262032                0.02859434
    ## [2,]               0.38335586                0.72650345
    ##      Late_blastocyst.3_Cell.1 Late_blastocyst.3_Cell.2
    ## [1,]               0.03175874               0.05652637
    ## [2,]               0.68189923               0.33610218
    ##      Late_blastocyst.3_Cell.3 Late_blastocyst.3_Cell.4
    ## [1,]               0.02525186               0.02236261
    ## [2,]               0.76867837               0.79827767
    ##      Late_blastocyst.3_Cell.5 Late_blastocyst.3_Cell.6
    ## [1,]              0.006395745              0.142058828
    ## [2,]              0.592685944              0.002148516
    ##      Late_blastocyst.3_Cell.7 Late_blastocyst.3_Cell.8
    ## [1,]              0.142058828               0.03646491
    ## [2,]              0.002148516               0.61180664

``` r
print(cc)
```

    ##  [1] 1 1 2 2 2 2 2 2 2 2 2 2 2 1 1 1 2 1 2 1 2 2 2 2 2 2 2 2
    ## [29] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 1 2 2 2 2
    ## [57] 2 2 2 2 1 2 2 2 2 2 2 2 1 2 2 2 2 2 2 1 1 2 1 2 2 2 2 2
    ## [85] 2 2 2 1 1 2

``` r
print(dd)
```

    ##          [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
    ## PNMA1__2    0    0    1    1    1    1    1    1    1     1
    ##          [,11] [,12] [,13] [,14] [,15] [,16] [,17] [,18]
    ## PNMA1__2     1     1     1     0     0     0     1     0
    ##          [,19] [,20] [,21] [,22] [,23] [,24] [,25] [,26]
    ## PNMA1__2     1     0     1     1     1     1     1     1
    ##          [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34]
    ## PNMA1__2     1     1     1     1     1     1     1     1
    ##          [,35] [,36] [,37] [,38] [,39] [,40] [,41] [,42]
    ## PNMA1__2     1     1     1     1     1     1     1     1
    ##          [,43] [,44] [,45] [,46] [,47] [,48] [,49] [,50]
    ## PNMA1__2     1     1     1     1     0     1     1     1
    ##          [,51] [,52] [,53] [,54] [,55] [,56] [,57] [,58]
    ## PNMA1__2     1     0     1     1     1     1     1     1
    ##          [,59] [,60] [,61] [,62] [,63] [,64] [,65] [,66]
    ## PNMA1__2     1     1     0     1     1     1     1     1
    ##          [,67] [,68] [,69] [,70] [,71] [,72] [,73] [,74]
    ## PNMA1__2     1     1     0     1     1     1     1     1
    ##          [,75] [,76] [,77] [,78] [,79] [,80] [,81] [,82]
    ## PNMA1__2     1     0     0     1     0     1     1     1
    ##          [,83] [,84] [,85] [,86] [,87] [,88] [,89] [,90]
    ## PNMA1__2     1     1     1     1     1     0     0     1

## (iv) directly apply qubic from the QUBIC package

We are trying save the result to a file
first.

``` r
WriteQubicInput <- function(file.name, data0, genes, q = 0, err = 1e-10) {
  cat("o", colnames(data0), "\n", file = file.name)
  for (i in genes) {
    cat(i, colnames(data0), "\n", file = "progress")
    x <- data0[i, ]
    Zcut0 <- q
    best <- GetBestK(x = x, n = 1000, q = Zcut0, err = 1e-10)
    if (best$k == 0) {
      next
    }
    pp <- calculate_prob_sep_Zcut(x, Zcut0, best$result[, 1], best$result[, 2], best$result[, 3])
    cc <- discretization_method_1_LLR_mean(x, best$result, pp, LLR_cut = 0.1)
    dd <- Build_R_matrix(cc, Zcut0, best$result[, 2], rownames(data0)[i])
    write.table(dd, file = file.name, col.names = FALSE, append = TRUE, quote = FALSE)
  }
}
```

``` r
system.time(WriteQubicInput("qubic_input_head", data0, head(selected.genes)))
```

    ##    user  system elapsed 
    ##    2.60    0.02    2.72

``` r
system.time(WriteQubicInput("qubic_input_head30", data0, head(selected.genes, 30)))
```

    ##    user  system elapsed 
    ##    8.94    0.06    9.22

This may be slow…

``` r
print(length(selected.genes))
```

    ## [1] 14542

``` r
qubic.file = "qubic_input"
if (!file.exists(qubic.file)) {
  WriteQubicInput(qubic.file, data0, selected.genes)
}
```

It is the time to read all the data
back.

``` r
qubic.input <- as.matrix(read.table(qubic.file, row.names = 1, header = TRUE))
```

Run QUBIC(Zhang et al. 2017), need several minute.

``` r
library(QUBIC)
if (!file.exists("res.RData")) {
  res <- qubiclust_d(qubic.input)
  save(res, file="res.RData")
} else {
  load("res.RData")
}
```

``` r
save(res, file="res.RData")
```

## (v) results summary

``` r
res
```

    ## 
    ## An object of class Biclust 
    ## 
    ## call:
    ##  NULL
    ## 
    ## Number of Clusters found:  100 
    ## 
    ## First  5  Cluster sizes:
    ##                    BC 1 BC 2 BC 3 BC 4 BC 5
    ## Number of Rows:     934  931  921  959  941
    ## Number of Columns:   73   73   73   70   71

``` r
biclust::summary(res)
```

    ## 
    ## An object of class Biclust 
    ## 
    ## call:
    ##  NULL
    ## 
    ## Number of Clusters found:  100 
    ## 
    ## Cluster sizes:
    ##                    BC 1 BC 2 BC 3 BC 4 BC 5 BC 6 BC 7 BC 8
    ## Number of Rows:     934  931  921  959  941  911  909  884
    ## Number of Columns:   73   73   73   70   71   73   73   75
    ##                    BC 9 BC 10 BC 11 BC 12 BC 13 BC 14 BC 15
    ## Number of Rows:     882   917   902   966   885   872   872
    ## Number of Columns:   75    72    73    68    74    75    75
    ##                    BC 16 BC 17 BC 18 BC 19 BC 20 BC 21
    ## Number of Rows:      894   892   891   910   897   907
    ## Number of Columns:    73    73    73    71    72    71
    ##                    BC 22 BC 23 BC 24 BC 25 BC 26 BC 27
    ## Number of Rows:      870   870   905   866   850   849
    ## Number of Columns:    74    74    71    74    75    75
    ##                    BC 28 BC 29 BC 30 BC 31 BC 32 BC 33
    ## Number of Rows:      846   845   832   832   832   832
    ## Number of Columns:    75    74    75    75    75    75
    ##                    BC 34 BC 35 BC 36 BC 37 BC 38 BC 39
    ## Number of Rows:      831   831   831   831   831   831
    ## Number of Columns:    75    75    75    75    75    75
    ##                    BC 40 BC 41 BC 42 BC 43 BC 44 BC 45
    ## Number of Rows:      831   838   826   826   826   826
    ## Number of Columns:    75    74    75    75    75    75
    ##                    BC 46 BC 47 BC 48 BC 49 BC 50 BC 51
    ## Number of Rows:      832   831   863   808   814   792
    ## Number of Columns:    74    74    71    75    73    75
    ##                    BC 52 BC 53 BC 54 BC 55 BC 56 BC 57
    ## Number of Rows:      792   792   792   792   792   792
    ## Number of Columns:    75    75    75    75    75    75
    ##                    BC 58 BC 59 BC 60 BC 61 BC 62 BC 63
    ## Number of Rows:      792   792   792   792   792   792
    ## Number of Columns:    75    75    75    75    75    75
    ##                    BC 64 BC 65 BC 66 BC 67 BC 68 BC 69
    ## Number of Rows:      792   792   792   792   792   792
    ## Number of Columns:    75    75    75    75    75    75
    ##                    BC 70 BC 71 BC 72 BC 73 BC 74 BC 75
    ## Number of Rows:      792   792   792   792   792   792
    ## Number of Columns:    75    75    75    75    75    75
    ##                    BC 76 BC 77 BC 78 BC 79 BC 80 BC 81
    ## Number of Rows:      792   792   792   792   792   792
    ## Number of Columns:    75    75    75    75    75    75
    ##                    BC 82 BC 83 BC 84 BC 85 BC 86 BC 87
    ## Number of Rows:      792   792   792   792   792   792
    ## Number of Columns:    75    75    75    75    75    75
    ##                    BC 88 BC 89 BC 90 BC 91 BC 92 BC 93
    ## Number of Rows:      792   792   792   791   792   792
    ## Number of Columns:    75    75    75    75    74    74
    ##                    BC 94 BC 95 BC 96 BC 97 BC 98 BC 99
    ## Number of Rows:      792   792   739   725   699   683
    ## Number of Columns:    74    74    77    76    78    77
    ##                    BC 100
    ## Number of Rows:       668
    ## Number of Columns:     78

# References

<div id="refs" class="references">

<div id="ref-zhang16">

Zhang, Yu, Juan Xie, Jinyu Yang, Anne Fennell, Chi Zhang, and Qin Ma.
2017. “QUBIC: A Bioconductor Package for Qualitative Biclustering
Analysis of Gene Co- Expression Data.” *Bioinformatics* 33 (3): 450–52.
<https://doi.org/10.1093/bioinformatics/btw635>.

</div>

</div>
