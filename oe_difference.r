logrank_OE_per_time <- function(time, status, group) {
  # ensure inputs are vectors
  time <- as.numeric(time)
  status <- as.integer(status)  # 1=event, 0=censor
  group <- as.factor(group)
  groups <- levels(group)
  
  # prepare data frame and sort by time
  d <- data.frame(time = time, status = status, group = group, stringsAsFactors = FALSE)
  d <- d[order(d$time), ]
  
  # unique event times (times with at least one event)
  ev_times <- sort(unique(d$time[d$status == 1]))
  
  # init results
  res_list <- vector("list", length(ev_times))
  
  for (i in seq_along(ev_times)) {
    t <- ev_times[i]
    # risk set: those with time >= t (still at risk just before t)
    risk_idx <- which(d$time >= t)
    # events at time t
    event_idx <- which(d$time == t & d$status == 1)
    
    n_j <- length(risk_idx)
    d_j <- length(event_idx)
    
    # per-group counts
    n_g <- as.numeric(tapply(seq_along(d$time)[risk_idx], d$group[risk_idx], length))
    O_g <- as.numeric(tapply(d$status[event_idx], d$group[event_idx], sum))
    # tapply returns NA for groups with zero; replace with 0
    n_g[is.na(n_g)] <- 0
    O_g[is.na(O_g)] <- 0
    
    # expected under log-rank for each group
    E_g <- (n_g / sum(n_g)) * d_j
    
    res_list[[i]] <- data.frame(
      time = t,
      n_total = n_j,
      d_total = d_j,
      group = groups,
      n_at_risk = n_g,
      observed = O_g,
      expected = E_g,
      O_minus_E = O_g - E_g,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }
  
  # bind and return
  res <- do.call(rbind, res_list)
  # order columns nicely
  res <- res[, c("time","group","n_total","d_total","n_at_risk","observed","expected","O_minus_E")]
  return(res)
}


survdiff(formula = Surv(futime, fustat) ~ rx, data = ovarian)
oe1 <- subset( with(ovarian, logrank_OE_per_time(futime, fustat, rx)), group==1)


o1 <- ovarian[with(ovarian, !(rx==1 & futime > 400)), ]
 with(o1, logrank_OE_per_time(futime, fustat, rx)) 
 
 
o2 <-  ovarian
o2[o2$rx==1 & o2$futime > 400,'fustat'] <- 1

 with(o2, logrank_OE_per_time(futime, fustat, rx)) 
 
plot(survfit(Surv(futime, fustat) ~ rx, data=o2))
 
 
 o3 <-  ovarian
o3[o3$rx==1 & o3$futime > 400,'fustat'] <- 0

 with(o3, logrank_OE_per_time(futime, fustat, rx)) 
