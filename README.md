# hartigan_Capstone
---
title: "Hartigan_Chrissy_Capstone"
output:
  html_document: default
  html_notebook: default
---

```{r}
library(tidyverse)
library(datapasta)
library(ez)
library(viridis)
library(Hmisc)
library(lubridate)
library(survival)
library(survminer)
library(knitr)
library(kableExtra)
library(coin)
```

Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Cmd+Option+I*.

When you save the notebook, an HTML file containing the code and output will be saved alongside it (click the *Preview* button or press *Cmd+Shift+K* to preview the HTML file).


## Background

Transplant is a life-saving cure for end stage organ failure. Patients who recieve organ trasnplants have to remain on immunosuppressive drugs for the remainder of their lives in order to prevent T cell mediated rejection of the graft. Small molecule immunosuppressants such as calcineurin inhibitors are extremely toxic due to their non-specific nature. Belatacept is a CTLA-4Ig fusion protein that specifically inhibits T cell costimulatory signals to inhibit T cell mediated graft rejection. CTLA-4Ig treatment is associated with increased instances of acute rejection in the weeks following transplant, so we are looking for other coinhibitory or costimulatory pathways to synergize with CTLA-4Ig to improve immunosuppression for transplant patients. 

TIGIT is an inhibitory molecule that is expressed on activated T cells as well as a subset of T regulatory cells. Treatment of cells with an agonistic antibody for TIGIT has been shown to enhance the ability of Tregs to suppress Th1 and Th17 responses (some of the cells involved in acute rejection after CTLA-4Ig treatment). **I hypothesize that TIGIT agonism can synergize with CTLA-4 to be an effective immunosuppressive regimen in transplant.**

To study this, we use a mouse model of skin graft. We use a minor antigenic mismatch modelto study survival of the graft over time to see if combination treatment with the TIGIT agonist and CTLA-4Ig improves graft survival compared to CTLA-4Ig treatment alone. In another iteration of the experiment, we sacrifice the mice 10 days after grafting in order to study the cellular response to treatment by analyzing phenotype and cytokine response by flow cytometry. For the purposes of this Capstone, I will describe the statistical evaluation of active-caspase 3 expression, just one out of many proteins probed in my flow cytometry pannel.

## The question:

It is unknown how TIGIT agonism enhances the ability of Tregs to suppress T cell responses. Additionally, it is unknown whether TIGIT agonism will be affected by treatment CTLA-4Ig. Further, this pathway has not been studied in the context of transplant/alloimmunity. We want to understnad if TIGIT agonism will improve graft survival, minimize allograft specific T cell responses, and if  alloreactive T cells are reduced, is the mechanism induction of apoptosis?

## Anticipated results:

If TIGIT agonism is an effective immunosuppressant, we should see an improvement in graft survival. If TIGIT agonism reduces alloreactive T cells in the draining lymph nodes or graft tissue, then we will be able to use flow cytometry to assess whether this reduction is due to cell death by detection of the active form of the enzyme caspase-3.

## The variables:

**Survival:** Outcome variable is a discrete measurement of time to rejection measured in days. We do bilateral grafts and count rejection as when both grafts have fallen off or are 10% of their original size. 

**Active Caspase-3 expression** Cell numbers and frequencies identified in flow cytometry experiments are discrete measured outcome variables. In the example below, the frequency of alloreactive cells that express Caspase-3 will be our measured outcome variable. The predictor variables will be the treatment groups: noRX, aTIGIT, CTLA-4Ig, and combo (CTLA-4Ig + aTIGIT). 

## Statistical Hypotheses

**Survival**
Null hypothesis: Treatment with CTLA-4Ig will have the same median survival time as TIGIT agonist + CTLA-4Ig: there will be no difference (MST)

$H_0 = \mu_{CTLA-4Ig}$ **=** $\mu_{combo}$

ALternative hypothesis: The treatments improve graft survival, independent use of TIGIT agonism will be less effective than CTLA-4Ig because the CD28 pathway is more ubiquitous and necessary for T cell activation than the TIGIT pathway is for inhibition. Combination therapy will be more effective than either of the treatments in extending graft survival.

$H_1 = \mu_{CTLA-4Ig}$ **<** $\mu_{combo}$

**Mechanism of suppression, cell death**
Null hypothesis: Immunosuppressive therapy will have no change on the expression of active caspase 3 expression in alloreactive T cells (therapy will not have an effect on cell death).

$H_0 = \mu_{noRX}$ **=** $\mu_{CTLA-4Ig}$ **=** $\mu_{aTIGIT}$ **=** $\mu_{combo}$

Alternative hypothesis: combination theraly will have increased active caspase 3 expression in alloreactive T cells compared to CTLA-4Ig treatment or TIGIT agonism alone, but all will have increased caspase activity compared to no treatment.

$H_1 = \mu_{noRX}$  **<** $\mu_{aTIGIT}$ ***<** $\mu_{CTLA-4Ig}$ **<** $\mu_{combo}$


## Statistical Methods

**Survival**
I will assess the survival of grafts using Kaplan-Meier analysis to look at differences in median survival time. I will compare hazard ratios between groups to determine the effect of treatment on each of the groups. Although differential graft survival will not tell us if there is or is not a cellular mechanism to TIGIT agonism, it will guide our hypotheses in the effect of TIGIT agonism (as in, does it work on its own or does it rely on CTLA-4Ig treatment).

**Mechanisms**
I will use ANOVA to anaylze difference across groups to compare active caspase 3 expression as an indication of cell death. I will use the ANOVA results to guide post-hoc analysis, likely using TUKEY to compare across all groups. Alternatively, I can do pairwise analyses to determine differences from the no-treatment group. However, the most important comparison for my question of whether or not TIGIT agonism will synergize with CTLA-4Ig is to compate the CTLA-4Ig treatment alone to the combination therapy group using a student's t test. We are only testing one dose so we do not need to do any sort of regression alaysis. If the data suggests differences between treatments but we cannot reject the null hypothesis, then I would consider doing dose respnse experiments in which case I would use mixed linear regression analyses to determine differences in treatment over different doses BUT I do not anticipate that at this stage.


## Experimental Details

Wildtype C67/Bl6 mice will receive antigen-specific T cells 1 day before graft. Day 0 mice will receive two skin grafts from mOVA (specific antigen) expressing mice (1 each ear and tail skin). Bandaids are cut at day 7 and mice will be monitored daily for graft survival. Technical failures will be considered grafts that did not adhere by day 7. In another iteration of the experiment, mice will be sacrificed at day 10. Each mouse will count as a biological replicate. I will do a power analysis to determine the number of mice necessary per group to be able to make a sound, rigerous conclusion. For the analysis of cells by flow cytometry, I can do technical replicates for each biological replicate by taking two aliquots of cells from the organ of interest (graft or draining lymph nodes) and plating them for in vitro stimulation. 

**Survival** 
We define rejection by both grafts being lost (fall off, necrotized), or less than 10% of their original size.

**Cellular Markers by flow cytometry**
The measurement outputs will be frequency and numbers of cells. Cell numbers will be calculated using CountBrite beads. I will set a type 1 error of 5%, meaning that differences between groups must have a p value <0.05 in order for us to reject the null hypothesis. For type 2 error, I want to set my experiment at 80% power. I will base my power calculation on the values for mean and standard deviation I got in a preliminary set up of the experiment. **The power calculation below is based on an experiment with a type-1 error of 95%, based on values from an experiment of n=5/group. I plugged in numbers to estimate how to acheive 80% power, or 20% type-2 error, of t tests comparing each drug to no treatment, and importantly CTLA-4Ig to combination therapy. The results ranged from 5-20 mice per group. I will therefore perform experimetns with 20 mice per group (n=5/group, repeat experiment 4 independent times).**


```{r}
## m1=3.528; sd1=1.090; m2=7.900; sd2=2.430; m3=4.834; sd3=1.687; m4=8.656; sd4=2.989

## power calculation for measuring differences in caspase expression.
t.pwr <- function(n){
  
  m3=4.834; sd3=1.687; m4=8.656; sd4=2.989
  						
  ssims=1000
  p.values <- c()
  i <- 1
  
  repeat{
    x=rnorm(n, m3, sd3); 
    y=rnorm(n, m4, sd4);
   
    p <- t.test(x, y, 
                paired=F, 
                alternative="two.sided", 
                var.equal=F,
                conf.level=0.95)$p.value
    p.values[i] <- p
    if (i==ssims) break
    i = i+1
    pwr <- length(which(p.values<0.05))/ssims
  }
  return(pwr)
}

frame <- data.frame(n=2:30)
data <- bind_cols(frame, 
                  power=apply(frame, 1, t.pwr))

ggplot(data, aes(n, power))+
  geom_point() +
  scale_y_continuous(breaks=c(seq(0, 1, 0.1)))+
  scale_x_continuous(breaks=c(seq(0,30,2)))+
  labs(x="n per group")

## 80% power- noRX, aTIGIT = 7 mice/group
## 80% power- noRX, CTLA-4Ig = 20 mice/group
## 80% power- noRX, CTLA-4Ig = 5 mice
## #80% power- CTLA4Ig, combo = 8 mice
```


## Survival Experiment
**Survival Data**

```{r}
## I am using real data because I am unsure how to simulate data of this kind, we usually use Kaplan-Meier analysis in GraphPad Prism to do these analyses. 

Graft <-read.csv("/Users/christinahartigan/IBS538/Survival2.csv"); Graft

##Summarize it
GraftSurv <- survfit(Surv(Day, Event) ~ Treatment, data = Graft)
summary(GraftSurv)

GraftSurv

##plot it
ggsurvplot(GraftSurv,
          conf.int = TRUE,
          surv.median.line = "hv", # Specify median survival
          ggtheme = theme_bw(), # Change ggplot2 theme
          palette = c("#002878", "#d28e00"))

##quantify it 
survdiff(Surv(Day, Event) ~ Treatment, data = Graft)
pchisq(12.28, df=1, lower.tail=F)

```

**Survival Write-Up**
In this model, we treat on days 0, 2, 4, and 6 and then stop treatment, so we anticipate rejection but look for differences in median survivial of the graft between treatment groups. Graft survival with CTLA-4Ig treatment (median survival = 26 days, 95%CI = 20 to undetermined) differs from survival with combination CTLA-4Ig + aTIGIT (median survival = 45 days, 95%CI = 31 to 3undetermind; Mantel-Haenszel log rank test, chisq= 11.2, p=0.0004578384). The elongation of median survival time indicated that treatment with TIGIT agonist is effective at mitigating alloreactivity. 

## Detection of Active Caspase-3 by Flow Cytometry
**The data: Active Caspase 3 expression - frequency active caspase 3 expressing cells of the total antigen specific CD8+ T cell population**

```{r}
##Simulating the data
set.seed(1234)
Treatment <- as.factor(c(rep("noRX", 20), rep("aTIGIT", 20), rep("CTLA-4Ig", 20), rep("Combo", 20)))
m1=3.528; sd1=1.090; m2=7.900; sd2=2.430; m3=4.834; sd3=1.687; m4=8.656; sd4=2.989
noRX_Casp3 <- rnorm(20, m1, sd1)
aTIGIT_Casp3 <- rnorm(20, m2, sd2)
CTLA4Ig_Casp3 <-rnorm(20, m3, sd3)
Combo_Casp3 <- rnorm(20, m4, sd4)
Freq <- c(noRX_Casp3, aTIGIT_Casp3, CTLA4Ig_Casp3, Combo_Casp3)

CaspData <- data.frame(Treatment, Freq)

##Visualizing the results:
ggplot(CaspData, aes(x=Treatment, y=Freq, group=Treatment)) +
  geom_point(size=1) +
  geom_boxplot(width = 0.5) +
  xlab("Treatment") +
  ylab("Frequency of alloreactive cells (% of CD8+)") 
 
  
##Run the ANOVA
ID <- as.factor(seq(1:80))

Casp2 <- CaspData %>%  add_column(
    ID=rep(ID, each=1), 
    .before=T)
results <- anova(lm(Freq ~ Treatment, Casp2)); results


##Post-hoc analysis: The F test shows that there is a significant difference between the groups. I am going to perform paitwise analyses across all groups to compaer the groups after this one-way completely randomized ANOVA.I chose bonferroni p value adjustment to be stringent.

post_hoc <- pairwise.t.test(Casp2$Freq, Casp2$Treatment, paired=FALSE, alternative="two.sided", pooled.sd=TRUE, p.adjust= "bonf"); post_hoc

```


**Active Caspase-3 Expression - Frequency of ALloreactive T cells**
Alloreactive T cells in the graft draining lymph nodes of mice will express active caspase 3 if they are undergoing apoptosis. There is no difference in levels of active caspase 3 between noRX and CTLA-4Ig teated mice (n=20 per group, p value = 1 by students T test with bonferroni p value adjustment). TIGIT agonism induces active caspase 3 in alloreactive cells compared to noRx (n=20, p value = 3.9e-5 by students T test with bonferroni p value adjustment). Interesting, the combination of CTLA-4Ig treatment with TIGIT agonism induces more active caspase 3 expression than either the TIGIT agonist (n=20, p value = 0.00014 with bonferroni p value adjustment) or CTLA-4Ig treatment (n=20, p value = 4e-11, with bonferroni o value adjustment). Together, these data show that TIGIT agonism synergizes with CTLA-4Ig to mitigate alloreactive T cell responses through the induction of apoptosis. 

