#library(testthat)
#library(mvMORPH)
# Comment test_check in order to not run tests with R CMD check
# test_check("mvMORPH")

# test-suite for the functions mvBM/mvOU/mvEB/mvSHIFT/mvLL/mvSIM/mvOUTS/mvRWTS are under preparation

tests <- FALSE

if(tests==TRUE){
    require(mvMORPH)
    # General tree to use in fit and simulations
    set.seed(1)
    
    # Generating a random tree
    #tree<-phylo<-pbtree(n=50)
    # Setting the regime states of tip species
    #sta<-as.vector(c(rep("Forest",20),rep("Savannah",30))); names(sta)<-tree$tip.label
    # Making the simmap tree with mapped states
    #tree<-make.simmap(tree,sta , model="ER", nsim=1)
    
    # for reproducibility
    tree <- read.simmap(text="((t37:{Forest,0.17911043},t38:{Forest,0.17911043}):{Forest,4.19021273},(((t39:{Forest,0.1406569},t40:{Forest,0.1406569}):{Forest,0.07210538},t36:{Forest,0.21276228}):{Forest,3.77896996},((((((t21:{Forest,0.474973},t22:{Forest,0.474973}):{Forest,0.11994112},(t28:{Forest,0.3176896},(t41:{Forest,0.11335236},t42:{Forest,0.11335236}):{Forest,0.20433723}):{Forest,0.27722453}):{Forest,0.34242237},t8:{Forest,0.9373365}):{Forest,0.87640728},(((t12:{Forest,0.73431592},t13:{Forest,0.73431592}):{Forest,0.69159925},(t2:{Forest,1.31658647},(t18:{Forest,0.57994406},t19:{Forest,0.57994406}):{Forest,0.73664241}):{Forest,0.1093287}):{Forest,0.23284738},(t1:{Forest,1.54484389},((t16:{Forest,0.59318045},t17:{Forest,0.59318045}):{Forest,0.33458062},(t34:{Forest,0.25021113},t35:{Savannah,0.17781996:Forest,0.07239117}):{Forest,0.67754994}):{Forest,0.61708282}):{Forest,0.11391867}):{Forest,0.15498123}):{Forest,0.33873822},(t4:{Savannah,1.234117},t5:{Savannah,1.234117}):{Savannah,0.0748023:Forest,0.8435627}):{Forest,1.17266697},(((t3:{Savannah,1.25848874},((t14:{Savannah,0.60935948},t15:{Savannah,0.60935948}):{Savannah,0.08390608},(t20:{Savannah,0.57975625},((t24:{Savannah,0.32641206},(t47:{Savannah,0.01395592},t48:{Savannah,0.01395592}):{Savannah,0.31245613}):{Savannah,0.19000652},(t29:{Savannah,0.30478906},t30:{Savannah,0.30478906}):{Savannah,0.21162951}):{Savannah,0.06333768}):{Savannah,0.1135093}):{Savannah,0.56522318}):{Savannah,0.47279525},((((t31:{Savannah,0.28631849},(t45:{Savannah,0.01659224},t46:{Savannah,0.01659224}):{Savannah,0.26972625}):{Savannah,0.15599405},t23:{Savannah,0.44231253}):{Savannah,0.33983399},t11:{Savannah,0.78214652}):{Savannah,0.11039434},(t25:{Savannah,0.32607817},(t49:{Savannah,0.00334081},t50:{Savannah,0.00334081}):{Savannah,0.32273736}):{Savannah,0.5664627}):{Savannah,0.83874312}):{Savannah,0.87012285},(((t26:{Savannah,0.32102325},t27:{Savannah,0.32102325}):{Savannah,0.59890242},t9:{Savannah,0.91992568}):{Savannah,1.43882358},(((t6:{Savannah,1.06510523},t7:{Savannah,1.06510523}):{Savannah,0.02627092},(t43:{Savannah,0.02516062},t44:{Savannah,0.02516062}):{Savannah,1.06621553}):{Savannah,0.79310082},((t32:{Savannah,0.27743011},t33:{Savannah,0.27743011}):{Savannah,0.5706055},t10:{Savannah,0.84803561}):{Savannah,1.03644137}):{Savannah,0.47427228}):{Savannah,0.24265758}):{Savannah,0.19050932:Forest,0.53323281}):{Forest,0.66658327}):{Forest,0.37759092});")
    
    col<-c("blue","orange"); names(col)<-c("Forest","Savannah")
    # Plot of the phylogeny for illustration
    plotSimmap(tree,col,fsize=0.6,node.numbers=FALSE,lwd=3, pts=FALSE)
    
    # for reproducibility without seed
    #trait <-
   
    # General time-series
    timeseries <- 0:49
    ## tests on mvSIM
    
    ## Bivariate case
    # generate sigma
    sigma <- matrix(m <- runif(4), 2, 2) %*% t(matrix(m, 2, 2))
    sigbmm <- list(matrix(m <- runif(4), 2, 2) %*% t(matrix(m, 2, 2)),
                    matrix(m <- runif(4), 2, 2) %*% t(matrix(m, 2, 2)))
    # generate alpha
    alpha <- matrix(m <- runif(4), 2, 2) %*% t(matrix(m, 2, 2))
    # generate theta
    ou1 <- bm <- c(0, 2)
    oum <- bmD <- matrix(c(0,1,0.5,2),2,2) # traits are on columns
    # generate beta
    beta <- -0.05
    # generate trend
    trend <- c(0.02, 0.02)
    
    # param BM1
    parbm1 <- list(sigma=sigma, theta=bm)
    # param BM1 trend
    parbm1t <- list(sigma=sigma, theta=bm, trend=trend)
    # param BM1 multiple means
    parbm1D <- list(sigma=sigma, theta=bmD, smean=FALSE)
    # param BMM
    parbmm <- list(sigma=sigbmm, theta=bm)
    # param BMM trend
    parbmmt <- list(sigma=sigbmm, theta=bm, trend=trend)
    # param BMM multiple means
    parbmmD <- list(sigma=sigbmm, theta=bmD, smean=FALSE)
    
    # param OU1
    parou1 <- list(sigma=sigma, alpha=alpha, theta=ou1)
    # param OUM
    paroum <- list(sigma=sigma, alpha=alpha, theta=oum)
    # param OUM
    paroumR <- list(sigma=sigma, alpha=alpha, theta=ou1, root=FALSE)
    
    # param EB
    pareb <- list(sigma=sigma, theta=bm, beta=beta)
    
    # Simulate datasets
    dat1 <- mvSIM(tree, model="BM1", nsim=1, param=parbm1)
     dat2 <- mvSIM(tree, model="BM1", nsim=1, param=parbm1t)
      dat3 <- mvSIM(tree, model="BM1", nsim=1, param=parbm1D)
       dat4 <- mvSIM(tree, model="BMM", nsim=1, param=parbmm)
        dat5 <- mvSIM(tree, model="BMM", nsim=1, param=parbmmt)
         dat6 <- mvSIM(tree, model="BMM", nsim=1, param=parbmmD)
          dat7 <- mvSIM(tree, model="OU1", nsim=1, param=parou1)
           dat8 <- mvSIM(tree, model="OUM", nsim=1, param=paroum)
            dat9 <- mvSIM(timeseries, model="OUTS", nsim=1, param=paroum)
             dat10 <- mvSIM(timeseries, model="OUTS", nsim=1, param=paroumR)
              dat11 <- mvSIM(tree, model="EB", nsim=1, param=pareb)
             
    # Intentional errors
    test1 <- mvSIM(tree, model="OUTS", nsim=1, param=paroum) # tree instead of ts
     test2 <- mvSIM(timeseries, model="OU1", nsim=1, param=parou1) # ts instead of tree
      test3 <- mvSIM(timeseries, model="OUTS", nsim=1, param=parou1) # wrong number of theta
       test4 <- mvSIM(model="OUTS", nsim=1, param=paroum) # missing ts/tree
       
       
       ## Model fit
       
         fit1 <- mvBM(tree, dat1, model="BM1")
           fit2 <- mvBM(tree, dat2, model="BM1", param=list(trend=T)) # must be fitted on a n-ultrametric tree
             fit3 <- mvBM(tree, dat3, model="BM1", param=list(smean=F))
               fit4 <- mvBM(tree, dat4, model="BMM")
                 fit5 <- mvBM(tree, dat5, model="BMM", param=list(trend=T)) # must be fitted on a n-ultrametric tree
                   fit6 <- mvBM(tree, dat6, model="BMM", param=list(smean=F))
                     fit7 <- mvOU(tree, dat7, model="OU1")
                     fit8 <- mvOU(tree, dat8, model="OUM", param=list(vcv="randomRoot")) # avoid large values sometimes; known problem to fix?
                         fit9 <- mvOUTS(timeseries, dat9)
                           fit10 <- mvOUTS(timeseries, dat10, param=list(root=F))
                             fit11 <- mvRWTS(timeseries, dat10)
                               fit11 <- mvEB(tree, dat11)
                               
        # Tests for errors
        test5 <- mvBM(tree, dat2, model="BM1", param=list(trend=c(1,1,1,1), smean=F))
        
        ## Model fit with NA values and defaults methods
        
        dat1[8,2]<- dat1[25,1] <- NA
          dat2[8,2]<- dat2[25,1] <- NA
            dat3[8,2]<- dat3[25,1] <- NA
              dat4[8,2]<- dat4[25,1] <- NA
                dat5[8,2]<- dat5[25,1] <- NA
                  dat6[8,2]<- dat6[25,1] <- NA
                    dat7[8,2]<- dat7[25,1] <- NA
                      dat8[8,2]<- dat8[25,1] <- NA
                        dat9[8,2]<- dat9[25,1] <- NA
                          dat10[8,2]<- dat10[25,1] <- NA
                          # then fit
                          fit1 <- mvBM(tree, dat1, model="BM1")
                          fit2 <- mvBM(tree, dat2, model="BM1", param=list(trend=T))
                          fit3 <- mvBM(tree, dat3, model="BM1", param=list(smean=F))
                          fit4 <- mvBM(tree, dat4, model="BMM")
                          fit5 <- mvBM(tree, dat5, model="BMM", param=list(trend=T))
                          fit6 <- mvBM(tree, dat6, model="BMM", param=list(smean=F))
                          fit7 <- mvOU(tree, dat7, model="OU1")
                          fit8 <- mvOU(tree, dat8, model="OUM")
                          fit9 <- mvOUTS(timeseries, dat9)
                          fit10 <- mvOUTS(timeseries, dat10, param=list(root=F))
                          fit11 <- mvRWTS(timeseries, dat10)
                          fit11 <- mvEB(tree, dat11)


# Expected likelihood

# for reproducibility (dat1 converted to vector)
traits <- matrix(c(0.759912758269539, 0.614545055791663, 1.99523967616929, 2.06050645201129, 1.72847166301803, -2.74185937878331, -1.74984543950876, -2.06239818030884, -2.06041993850937, -1.79080629892104, -1.20081330679844, -0.54349593713857, -0.0607971864579239, -0.00233277135380217, -0.400693887343343, -1.65772161485507, -0.182544756090116, -0.654109529770913, -0.737499486430531, -1.72515976761732, -1.72562879593648, -0.25330550757395, 0.729897907633583, -0.246081522346491, 0.151078889032263, 0.0138479096177918, -0.848009457276881, -0.581243223593088, -0.667057834272238, -0.66847007470311, 0.0981507462333929, 0.172861363132646, -0.353171766681764, -0.451300105973155, -0.368139106479026, -0.0627361956502809, -0.709892073985014, -0.875683848248537, -0.570848855977929, -0.532265981711379, -0.369679020727696, 0.0776255560504741, 0.138142231533873, -0.833365792498594, -0.309928572800132, -1.35885805922583, -1.14218922343049, 1.48087206844568, 1.09136693791236, -0.0180692260844691, 3.23281192282276, 3.00136528638948, 5.32074612474669, 5.41402685056098, 4.913298892713, -2.23542682300963, -0.726407146095314, -1.1868584417652, -1.24452289327373, -0.784150890488637, 0.146694082761811, 1.32180581471685, 2.05153383820109, 2.03959665464931, 1.48623453081428, -0.518344824420335, 1.70143196101035, 1.04697657897734, 0.892292938362765, -0.630853524174634, -0.629376986734026, 1.6045525143076, 3.11635168256392, 1.62949197480133, 2.32593600511521, 2.02205029936445, 0.748804111684117, 1.16117939166309, 1.05667448503384, 1.05133732067517, 2.22897875073374, 2.34707599071464, 1.45209863662888, 1.34341733329564, 1.48063807048236, 1.95499546327074, 1.0100600946259, 0.727548824795377, 1.15112488908316, 1.20921909380616, 1.36592576262338, 2.05376926747503, 2.14226651915817, 0.712214997646626, 1.47055230647709, -0.116941097964031, 0.212904236741522, 4.41536257330546, 3.8195839768853, 2.09983331019748), ncol=2)

# univariate:
test6 <- mvLL(tree, traits[,1], method="pic")

    #$logl
    #[1] -37.28764
    #
    #$theta
    #[1] 0.3443453
    #
    #$sigma
    #[,1]
    #[1,] 0.3028476
    #
    #attr(,"class")
    #[1] "mvmorph" "loglik"

 test7 <- mvLL(tree, traits[,1], method="rpf")
  test8 <- mvLL(tree, traits[,1], method="sparse")
   test9 <- mvLL(tree, traits[,1], method="inverse")
    test10 <- mvLL(tree, traits[,1], method="pseudoinverse")
    
    all.equal(test6$loglik, test7$loglik)
    all.equal(test6$loglik, test8$loglik)
    all.equal(test6$loglik, test9$loglik)
    all.equal(test6$loglik, test10$loglik)
    
# multivariate:
    test11 <- mvLL(tree, traits, method="pic")
    
    #$logl
    #[1] 56.40704
    #
    #$theta
    #[1] 0.3443453 2.6059401
    #
    #$sigma
    #[,1]      [,2]
    #[1,] 0.3028476 0.4721665
    #[2,] 0.4721665 0.7377560
    #
    #attr(,"class")
    #[1] "mvmorph" "loglik"
    
    test12 <- mvLL(tree, traits, method="pic", param=list(estim=FALSE, mu=c(1,1)))
    
    #$logl
    #[1] -1983.862
    #
    #$theta
    #[1] 1 1
    #
    #attr(,"class")
    #[1] "mvmorph" "loglik"
    
    test13 <- mvLL(tree, traits, method="pic", param=list(estim=FALSE, mu=c(1,1), sigma=matrix(c(2,1,1,1.5),2,2)))
    
    #$logl
    #[1] -115.824
    #
    #$theta
    #[1] 1 1
    #
    #attr(,"class")
    #[1] "mvmorph" "loglik"

## ----------------- 3 traits ------------------- ##

## Bivariate case
# generate sigma
sigma <- matrix(m <- runif(9), 3, 3) %*% t(matrix(m, 3, 3))
sigbmm <- list(matrix(m <- runif(9), 3, 3) %*% t(matrix(m, 3, 3)),
matrix(m <- runif(9), 3, 3) %*% t(matrix(m, 3, 3)))
# generate alpha
alpha <- matrix(m <- runif(9), 3, 3) %*% t(matrix(m, 3, 3))
# generate theta
ou1 <- bm <- c(0, 2, 1)
oum <- bmD <- matrix(c(0,1,0,1.5,0.5,2),ncol=3, nrow=2) # traits are on columns
# generate beta
beta <- -0.05
# generate trend
trend <- c(0.02, 0.02, 0.02)

# param BM1
parbm1 <- list(sigma=sigma, theta=bm)
# param BM1 trend
parbm1t <- list(sigma=sigma, theta=bm, trend=trend)
# param BM1 multiple means
parbm1D <- list(sigma=sigma, theta=bmD, smean=FALSE)
# param BMM
parbmm <- list(sigma=sigbmm, theta=bm)
# param BMM trend
parbmmt <- list(sigma=sigbmm, theta=bm, trend=trend)
# param BMM multiple means
parbmmD <- list(sigma=sigbmm, theta=bmD, smean=FALSE)

# param OU1
parou1 <- list(sigma=sigma, alpha=alpha, theta=ou1)
# param OUM
paroum <- list(sigma=sigma, alpha=alpha, theta=oum)
# param OUM
paroumR <- list(sigma=sigma, alpha=alpha, theta=ou1, root=FALSE)

# param EB
pareb <- list(sigma=sigma, theta=bm, beta=beta)

# Simulate datasets
dat1 <- mvSIM(tree, model="BM1", nsim=1, param=parbm1)
 dat2 <- mvSIM(tree, model="BM1", nsim=1, param=parbm1t)
  dat3 <- mvSIM(tree, model="BM1", nsim=1, param=parbm1D)
   dat4 <- mvSIM(tree, model="BMM", nsim=1, param=parbmm)
    dat5 <- mvSIM(tree, model="BMM", nsim=1, param=parbmmt)
     dat6 <- mvSIM(tree, model="BMM", nsim=1, param=parbmmD)
      dat7 <- mvSIM(tree, model="OU1", nsim=1, param=parou1)
       dat8 <- mvSIM(tree, model="OUM", nsim=1, param=paroum)
        dat9 <- mvSIM(timeseries, model="OUTS", nsim=1, param=paroum)
         dat10 <- mvSIM(timeseries, model="OUTS", nsim=1, param=paroumR)
          dat11 <- mvSIM(tree, model="EB", nsim=1, param=pareb)

# Intentional errors
test1 <- mvSIM(tree, model="OUTS", nsim=1, param=paroum) # tree instead of ts
    test2 <- mvSIM(timeseries, model="OU1", nsim=1, param=parou1) # ts instead of tree
        test3 <- mvSIM(timeseries, model="OUTS", nsim=1, param=parou1) # wrong number of theta
            test4 <- mvSIM(model="OUTS", nsim=1, param=paroum) # missing ts/tree


## Model fit

fit1 <- mvBM(tree, dat1, model="BM1")
 fit2 <- mvBM(tree, dat2, model="BM1", param=list(trend=T))
  fit3 <- mvBM(tree, dat3, model="BM1", param=list(smean=F))
   fit4 <- mvBM(tree, dat4, model="BMM")
    fit5 <- mvBM(tree, dat5, model="BMM", param=list(trend=T))
     fit6 <- mvBM(tree, dat6, model="BMM", param=list(smean=F))
      fit7 <- mvOU(tree, dat7, model="OU1")
       fit8 <- mvOU(tree, dat8, model="OUM")
        fit9 <- mvOUTS(timeseries, dat9)
         fit10 <- mvOUTS(timeseries, dat10, param=list(root=F))
          fit11 <- mvRWTS(timeseries, dat10)
           fit11 <- mvEB(tree, dat11)
           
           fit12 <- mvOU(tree, dat8, model="OUM", param=list(vcv="randomRoot"))
           fit13 <- mvOU(tree, dat8, model="OUM", param=list(vcv="randomRoot", root=TRUE))
           fit14 <- mvOU(tree, dat8, model="OUM", param=list(vcv="randomRoot", alpha=alpha)) # error
           fit15 <- mvOU(tree, dat8, model="OUM", param=list(vcv="randomRoot", alpha=alpha, decomp="qr")) # error
           fit16 <- mvOU(tree, dat8, model="OUM", param=list(vcv="randomRoot", sigma=sigma)) # error
           fit17 <- mvBM(tree, dat1, model="BM1", param=list(sigma=sigma))
           fit18 <- mvBM(tree, dat1, model="BMM", param=list(sigma=sigbmm))
           fit19 <- mvOU(tree, dat8, model="OUM", param=list(vcv="randomRoot", alpha=runif(6), decomp="choles")) # error for the name
           fit20 <- mvOUTS(timeseries, dat9, param=list(vcv="randomRoot", alpha=alpha, decomp="qr"))
           fit21 <- mvOUTS(timeseries, dat9, param=list(vcv="randomRoot", decomp="qr", alpha=alpha, decompSigma="eigen+"))

## Model fit with NA values and defaults methods

dat1[8,2]<- dat1[25,1] <- NA
dat2[8,2]<- dat2[25,1] <- NA
dat3[8,2]<- dat3[25,1] <- NA
dat4[8,2]<- dat4[25,1] <- NA
dat5[8,2]<- dat5[25,1] <- NA
dat6[8,2]<- dat6[25,1] <- NA
dat7[8,2]<- dat7[25,1] <- NA
dat8[8,2]<- dat8[25,1] <- NA
dat9[8,2]<- dat9[25,1] <- NA
dat10[8,2]<- dat10[25,1] <- NA
# then fit
fit1 <- mvBM(tree, dat1, model="BM1")
    fit2 <- mvBM(tree, dat2, model="BM1", param=list(trend=T))
        fit3 <- mvBM(tree, dat3, model="BM1", param=list(smean=F))
            fit4 <- mvBM(tree, dat4, model="BMM")
                fit5 <- mvBM(tree, dat5, model="BMM", param=list(trend=T))
                    fit6 <- mvBM(tree, dat6, model="BMM", param=list(smean=F))
                        fit7 <- mvOU(tree, dat7, model="OU1")
                            fit8 <- mvOU(tree, dat8, model="OUM")
                                fit9 <- mvOUTS(timeseries, dat9)
                                    fit10 <- mvOUTS(timeseries, dat10, param=list(root=F))
                                        fit11 <- mvRWTS(timeseries, dat10)
                                            fit12 <- mvEB(tree, dat11)
                                              fit13 <- mvRWTS(timeseries, dat10, param=list(trend=TRUE))


# simulate generic
d1 <- simulate(fit1,nsim=1, tree=tree)
    d2 <- simulate(fit2,nsim=1, tree=tree)
        d3 <- simulate(fit3,nsim=1, tree=tree)
            d4 <- simulate(fit4,nsim=1, tree=tree)
                d5 <- simulate(fit5,nsim=1, tree=tree)
                    d6 <- simulate(fit6,nsim=1, tree=tree)
                        d7 <- simulate(fit7,nsim=1, tree=tree)
                            d8 <- simulate(fit8,nsim=1, tree=tree)
                                d9 <- simulate(fit9,nsim=1, tree=timeseries)
                                    d10 <- simulate(fit10,nsim=1, tree=timeseries)
                                        d11 <- simulate(fit11,nsim=1, tree=timeseries)
                                            d12 <- simulate(fit12,nsim=1, tree=tree)
                                                d13 <- simulate(fit13,nsim=1, tree=timeseries)

} # end of tests

