##################
# case-crossover design for the entire state in pediatric population
# for synergistic effects of wildfire and extreme heat on hospitalization in California
# including codes to create a mock health dataset
# written by Chen Chen on 9/9/24
##################

## set the stage
outdir1 <- indir1 <- "" ## main working directory
if (!file.exists(file.path(indir1, "data"))) dir.create(file.path(indir1, "data"))
if (!file.exists(file.path(outdir1, "results"))) dir.create(file.path(outdir1, "results"))
library(data.table)

### code to create a mock health dataset--no association between exposure and outcome generated
#################
## read in a list of zctas with health data
zctas <- fread(file.path(indir1, "data", "zcta_list_1772zcta.csv"))
dd <- data.frame(dd=seq(from = as.Date("2006-01-01"), to = as.Date("2019-12-31"), by = 1), month=NA)
dd$month <- month(dd$dd)
dd <- dd$dd[dd$month >= 5 & dd$month <= 9] ## restricted to warm season
nn <- length(dd)
total.pop <- sum(zctas$pop, na.rm = TRUE)

set.seed(1234)
ha <- data.frame(
  zcta = rep(zctas$zcta, each=nn),
  pop = rep(zctas$pop, each=nn),
  date = rep(dd, time = nrow(zctas))
  )
setDT(ha)
ha[, c("Infectious_Parasitic_Disease",
       "enteritis", ##"Infectious_enteritis",
       "blood_immune", ##"Disorders of Blood and The Immune System",
       "endocrine", ##"Endocrine Nutritional and Metabolic Disease",
       "mental", ##"Mental Behavioral and Neurodevelopmental Disorders",
       "nervous", ##"Nervous System Conditions",
       "otitis", ##"Otitis_Media_Externa",
       "circulatory",
       "respiratory",
       "asthma", ##"Asthma",
       "digestive", ##"Digestive_Disease",
       "skin", ##"Skin Infections",
       "msk", ##"MSK Conditions",
       "gu", ##"GU Conditions",
       "renal", ##"Renal Conditions",
       "other", ##"Other SIgns and symptoms",
       "heat", ##"Heat Related Illness",
       "suicide_dep"##"Suicidality_Depression"
) := list(rpois(nn, lambda = pop * (821074/total.pop/nn)),
          rpois(nn, lambda = pop * (82859/total.pop/nn)),
          rpois(nn, lambda = pop * (22824/total.pop/nn)),
          rpois(nn, lambda = pop * (115464/total.pop/nn)),
          rpois(nn, lambda = pop * (484859/total.pop/nn)),
          rpois(nn, lambda = pop * (170295/total.pop/nn)),
          rpois(nn, lambda = pop * (530941/total.pop/nn)),
          rpois(nn, lambda = pop * (51312/total.pop/nn)),
          rpois(nn, lambda = pop * (2369108/total.pop/nn)),
          rpois(nn, lambda = pop * (359107/total.pop/nn)),
          rpois(nn, lambda = pop * (178055/total.pop/nn)),
          rpois(nn, lambda = pop * (116982/total.pop/nn)),
          rpois(nn, lambda = pop * (63387/total.pop/nn)),
          rpois(nn, lambda = pop * (221829/total.pop/nn)),
          rpois(nn, lambda = pop * (158188/total.pop/nn)),
          rpois(nn, lambda = pop * (833000/total.pop/nn)),
          rpois(nn, lambda = pop * (87183/total.pop/nn)),
          rpois(nn, lambda = pop * (142832/total.pop/nn))),
by = zcta]
write.csv(ha, file.path(outdir1, "data", "PEDS_combo_05_19_patzip_all_ICD_EM.csv"), row.names = FALSE)
#################

## run case-crossover at state-level with different definitions
## of climate change hazards (eh: 85th and 95th threshold; wf: 0, 5, 15, and 35 threshold) without restriction on zcta numbers (all 1772)
library(survival)
library(msm)
##################
wf.ths <- c(0, 5, 15, 35)
eh.ths <- c(85, 95)
types <- "binary"
lags <- c("same_day", "lag1")
outcomes <- c(
  "Infectious_Parasitic_Disease",
  "enteritis", ##"Infectious_enteritis",
  "blood_immune", ##"Disorders of Blood and The Immune System",
  "endocrine", ##"Endocrine Nutritional and Metabolic Disease",
  "mental", ##"Mental Behavioral and Neurodevelopmental Disorders",
  "nervous", ##"Nervous System Conditions",
  "otitis", ##"Otitis_Media_Externa",
  "circulatory",
  "respiratory",
  "asthma", ##"Asthma",
  "digestive", ##"Digestive_Disease",
  "skin", ##"Skin Infections",
  "msk", ##"MSK Conditions",
  "gu", ##"GU Conditions",
  "renal", ##"Renal Conditions",
  "other", ##"Other SIgns and symptoms",
  "heat", ##"Heat Related Illness",
  "suicide_dep", ##"Suicidality_Depression",
  "all_14" ## no infectious enteritis, asthma, renal disease or suicidal depression due to redundant codes
)

all.zcta <- fread(file.path(indir1, "data", "zcta_list_1772zcta.csv"))[, .(zcta, pop)] ## starting point of 1772 zctas with exposure and health data and starting with 9 in zcta


## clean health data to create dataset for case-crossover
ha <- fread(file.path(indir1, "data", "PEDS_combo_05_19_patzip_all_ICD_EM.csv")) ## mocked dataset created from above
ha <- ha[, all_14 := rowSums(.SD), .SDcols = c(4, 6:12, 14:17, 19, 20)] ## sum of 14 diagnostic codes
names(ha)[c(1, 3)] <- c("zcta", "case_date")

ha <- ha[ha$case_date > as.Date("2005-12-31") & ha$case_date < as.Date("2020-01-01"), ]
ha <- ha[ha$zcta %in% all.zcta$zcta, ]
length(unique(ha$zcta)) ## 1772 zipcodes
## create variables for merging
ha$wday <- wday(ha$case_date)
ha$month <- month(ha$case_date)
ha$year <- year(ha$case_date)

## calculation of reri--based on SAS codes from VanderWeele and Knol 2014 Appendix
## variance for joint effect based on deltamethod (car package)
## Note "se" in print-outs are se of transformed coefficients (delta method for joint effect, VanderWeele & Knol 2014 method for RERI). 
## The rest are se for coefficients of the model.
reri_interaction <- function(b, v, exposure1, exposure2) {
  reri <- exp(b[1] + b[2] + b[3]) - exp(b[1]) - exp(b[2]) + 1
  k1 <- exp(b[1] + b[2] + b[3]) - exp(b[1])
  k2 <- exp(b[1] + b[2] + b[3]) - exp(b[2])
  k3 <- exp(b[1] + b[2] + b[3])
  vreri <- v[1, 1] * k1 * k1 + 
    v[2, 2] * k2 * k2 + 
    v[3, 3] * k3 * k3 + 
    2 * v[1, 2] * k1 * k2 +
    2 * v[1, 3] * k1 * k3 +
    2 * v[2, 3] * k2 * k3
  se_reri <- sqrt(vreri)
  reri_ci95_l <- reri - 1.96 * se_reri
  reri_ci95_u <- reri + 1.96 * se_reri
  est <- exp(b)
  se <- sqrt(c(v[1, 1], v[2, 2], v[3, 3]))
  est_l <- exp(b - 1.96 * se)
  est_u <- exp(b + 1.96 * se)
  k3_se <- deltamethod(~exp(x1+x2+x3), b, v)
  k3_l <- k3 - 1.96 * k3_se
  k3_u <- k3 + 1.96 * k3_se
  temp <- data.frame( 
    est=c(est, reri, k3),
    ll=c(est_l, reri_ci95_l, k3_l),
    ul=c(est_u, reri_ci95_u, k3_u),
    se=c(se, se_reri, k3_se))
  rownames(temp) <- c(exposure1, exposure2, "interaction_multiplicative", "interaction_reri", "joint_effet")
  return(temp)
}

## analysis including all zctas regardless of their population and exposure types
for (wf.th in wf.ths) {
  for (eh.th in eh.ths) {
    for (tp in types) {
      cat("\n\n", "eh", eh.th, "wf", wf.th, tp, "\n")
      dataset <- paste0("eh", eh.th, "_wf", wf.th, "_", tp, "")
      exposure <- readRDS(file.path(indir1, "data", paste0(dataset,"_1772zcta_0619.rds")))
      
      ## calculation of state-wise rrs and reri
      out <- numeric()
      for (lag in lags) {
        if (lag!="same_day") { ## create lags in exposure
          nn <- as.numeric(gsub("lag", "", lag))
          new <- copy(exposure)
          new <- new[order(new$date), ]
          new[, eh:=shift(eh, n = nn, fill = NA, type = "lag"), by=.(zcta)]
          new[, wf:=shift(wf, n = nn, fill = NA, type = "lag"), by=.(zcta)]
        } else {
          new <- copy(exposure)
          new <- new[order(new$date), ]
        }
        new$wday <- wday(new$date)
        for (outcome in outcomes) {
          
          ha_ <- ha[eval(as.name(outcome))>0 & !is.na(eval(as.name(outcome))), ]
          ha_$id <- 1:nrow(ha_)
          dt <- merge(ha_, new, by = c("zcta", "wday", "month", "year"), all.x = TRUE, allow.cartesian=TRUE)
          
          dt$case <- ifelse(dt$case_date==dt$date, 1, 0)
          loc <- grep(outcome, names(dt))
          names(dt)[loc] <- "outcome"
          
          m <- clogit(case ~ eh*wf + strata(id), data=dt, weights=outcome, method="approximate")

          dt <- dt[!is.na(dt$wf & !is.na(dt$eh)), ] ## remove missing values
          m.eh <- clogit(case ~ eh + strata(id), data=dt, weights=outcome, method="approximate")

          m.wf <- clogit(case ~ wf + strata(id), data=dt, weights=outcome, method="approximate")

          ## create results for plotting
          temp <- data.frame(rbind(summary(m.eh)$conf.int[, c(1, 3, 4)], 
                                   summary(m.wf)$conf.int[, c(1, 3, 4)]))
          row.names(temp) <- c("extreme_heat", "wildfire")
          names(temp) <- c("est", "ll", "ul")
          temp$se <- c(sqrt(m.eh$var), sqrt(m.wf$var))
          temp <- rbind(reri_interaction(b=m$coefficients, v=m$var, "extreme_heat only", "wildfire only"),
                        temp)
          temp$estimate <- row.names(temp)
          out <- rbind(out, cbind(outcome = outcome, lag = lag, temp))
          
          m <- m.wf <- m.eh <- ha_ <- dt <- temp <- NULL
        }
        
        write.csv(out, file.path(outdir1, "results", paste0(dataset, "_state_model_summary_cause_specific.csv")), row.names = FALSE)
      } 
    }
  }
}

##################
