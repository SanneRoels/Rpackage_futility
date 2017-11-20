# remarks Sanne Roels

## for the implementation of the mixture part
*L 660*: the inclusion of 3 extra parameters (`MIXTURE`, `mix.weights` and `eventPriorWeightRobust` ) in the function `completeTrial.pooledArms`.

*L 720*: in our simulations i kept the drop-out rate fixed at 0.1/annually and not at the observed drop-out rate. In comment, this is the only difference with

*L 730*: the core of the mixture part: the `alpha` and `beta` become vectors in the case of MIXTURE = TRUE and remain scalars 

*L 750*: core of the implementation: the updating of the mixture weights. 

*L 780*: guarantee that the two eventrates get sampled with a appropriate weights in the nTrials simulations.  

*L 826*: save the mixture posterior weights (for QC); I added this to your output as well (remains NA if MIXTURE=FALSE).


## small generic remarks to further improve/clean up the code
*L 48*: in `FillinInterimdata.Pooled` suggestion on code cleaning for final edition -- crucial part is around line 168 
(I did not check whether the error also holds for the per-arm function `FillinInterimdata.byArm`)


```
geneRate <- function(N, rate,seed=NULL){
    ## changing the 'rate' parameter to a form that gives the
    ## prob(event in interval (0,1)) = rate, i.e., P(dropout time <=1) = 
    ## 1 - exp(-lambda) = rate
    ## via cumulative hazard.
    set.seed(seed)
    rexp(N, rate= -log( 1 - rate) )
}
VE=round(0.5*17/24+0.25*7/24,2)
# chunkc to test the code
    INC=0.02
    
    cat("\n", INC, "\n ")
    SEED = 1402
    
    #annual to weekly
    set.WINC <- INC/52
    
    # drop out weekly 
    do.WINC <- 0.1/52
    
    # n
    n.PLAC <- n.VAX <-1300
    
    # follow - up time   
    fu.W <- 104 
    
    # sample size re-est time point
    CHECK.W <- 52
    PP.W <- 28
    
    # Protocol accrual settings
    totAC.W <- 61
    # halved during first 13 weeks
    Half <- 13
    frac <- 1/(Half+(totAC.W-Half)*2)
    w.frac <- c(rep(frac,Half), rep(frac*2, totAC.W-Half))
    W.in <- 1300*w.frac
    
    # for now, no rounding
    Nhalf <- W.in[1:13]
    Nfull <- W.in[14:61]
    
    # uniform at half rate over first 13 weeks for vaccine and control
    set.seed(SEED)
    Thalf <- sort(runif(round(sum(Nhalf)), min=1, max=13)) 
    Tfull <- sort(runif(round(sum(Nfull)), min=14, max=61))
    IN.tC= c(Thalf, Tfull)
    rm(Thalf, Tfull)
    set.seed(SEED+3501)
    Thalf <- sort(runif(round(sum(Nhalf)), min=1, max=13)) 
    Tfull <- sort(runif(round(sum(Nfull)), min=14, max=61))
    IN.tV= c(Thalf, Tfull)
    rm(Thalf, Tfull)
    
    
    ## ---- Generate dropout/infection times ----
    inf.tC=geneRate(1300, set.WINC, seed=SEED*3)
    miss.tC=geneRate(1300, do.WINC, seed=SEED*5)
    inf.tV=geneRate(1300, (set.WINC*(1-VE)), seed=SEED*7)
    miss.tV=geneRate(1300, do.WINC,             seed=SEED*9)
    
    ## add entry time and infection/drop-out time ----
    INF.t=round(c(inf.tC, inf.tV) + c(IN.tC, IN.tV))
    MISS.t= round(c(miss.tC, miss.tV) + c(IN.tC, IN.tV))
    
    
    ## check for infection at cHECK.W
    inf.cw= ifelse((INF.t< MISS.t) & (INF.t < CHECK.W), 1, 0) 
    
    ## Determine the last observed time at Check.W
    last.ob = pmin(INF.t, MISS.t)
    last.ob = pmin(last.ob, rep(CHECK.W, times=length(last.ob)))
    
    ## identifier for `enrolled at check.w`
    IDX <- c(IN.tC, IN.tV) < CHECK.W    
    
    # only condider the the subjects enrolled at check.w
    inf.cw<- inf.cw[IDX]
    last.ob <- last.ob[IDX]
    
    # ---- create the interimDdata data.frame ----
    
    # assign treatment arm (2 datasets were inially concatenated + IDX -> correct subset)
    arm <- rep(c("C", "V"), each=1300)[IDX]
    
    # quasi eliminate the re-scheduling
    schedule <- rbinom(length(arm), 1, 0.00001)
    
    # entry time -> reference date is t0=0 
    entry <- c(IN.tC, IN.tV)[IDX]
    
    # event status
    event <- inf.cw
    
    # dropout status
    dropout <- ifelse(inf.cw==0 & MISS.t[IDX]< CHECK.W, 1, 0)
    
    # exit -> NA if in follow-up, last.ob else
    exit <- rep(NA, length(arm))
    exit[event==1] <- last.ob[event==1] 
    exit[dropout==1] <- last.ob[dropout==1] 
    
    # followup 
    followup <- ifelse(event==1 | dropout==1, 0, 1)
    
    # fill the interimData data.frame
    interimData <- data.frame(arm=arm, schedule2=schedule, entry=entry, exit=exit,
                              last_visit_dt=last.ob, event=event, dropout=dropout, complete=0,
                              followup=followup)

# robust mixture prior with uninformative part with w=1/1000 and informative part at w=0.5	
trialObj <- completeTrial.pooledArms(interimData=interimData, nTrials=2000, N=2600, enrollRate=NULL, enrollRatePeriod =52,
                                 eventPriorWeight=0.5, eventPriorRate=(1+0.57)/2*0.042, MIXTURE=TRUE, eventPriorWeightRobust = 1/1000,
                                 fuTime=104, visitSchedule=c(0,12,24,28,39,52,56,65,78,91,104),
                                 visitSchedule2=NULL, saveDir=NULL, randomSeed=SEED)

# calculte propbability to obtain more than e.g. 80 events
mean(unlist(lapply(trialObj$trialData, function(x){ sum(x$event)})) > 80)

# prior gamma with weight =0.5
trialObj <- completeTrial.pooledArms(interimData=interimData, nTrials=2000, N=2600, enrollRate=NULL, enrollRatePeriod =52,                                 eventPriorWeight=0.5, eventPriorRate=(1+0.57)/2*0.042, MIXTURE=FALSE, eventPriorWeightRobust = NULL,
                                 fuTime=104, visitSchedule=c(0,12,24,28,39,52,56,65,78,91,104),
                                 visitSchedule2=NULL, saveDir=NULL, randomSeed=SEED)
# calculte propbability to obtain more than e.g. 80 events
mean(unlist(lapply(trialObj$trialData, function(x){ sum(x$event)})) > 80)

``` 



