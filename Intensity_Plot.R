#### Intensity plots for Paper ####

Hom_Pois_fit <- PR_data
Hom_Pois_fit %>% head()
saveRDS(Hom_Pois_fit,file = "Simulations/Hom_Pois_fit_college.RDS")
# save this and repeat for Hawkes
Hom_Pois_fit <- readRDS("Simulations/Hom_Pois_fit_college.RDS")
Hom_Pois_fit[454,]

Hom_Hawk_fit <- PR_data
Hom_Hawk_fit %>% head()
saveRDS(Hom_Hawk_fit,file = "Simulations/Hom_Hawk_fit_college.RDS")


# then plot the estimated intensities from the online and batch methods 
# and compare them....
uniHawkesIntensityNumeric<-function(lambda0,alpha,beta, events, time.vec=NULL){
  if(is.null(time.vec)){
    delta.t <- tail(events,1)/2000
    time.vec <- seq(delta.t,tail(events,1),delta.t)
  }
  lambda.t <- rep(lambda0,length(time.vec))
  event.idx <- 2
  start.idx <- sum(time.vec<=events[2])+1
  r <- 0
  for(i in c(start.idx:length(time.vec))){
    current.t <- time.vec[i]
    if(current.t>events[event.idx+1]){
      event.idx <- event.idx + 1
      r <- exp(-beta*(events[event.idx]-events[event.idx-1]))*(1+r)
    }
    lambda.t[i]<-lambda0+alpha*exp(-beta*(current.t-events[event.idx]))*(1+r)
  }
  if(is.null(time.vec)){
    return(list(lambda.t=lambda.t,delta.t=delta.t,time.vec=time.vec))
  }else{
    return(list(lambda.t=lambda.t,time.vec=time.vec))
  }
}

i <- 454 # 225 is ok but not brilliant
# 454 could work.... if normalized...
# 4->3
# 136, 478
Hom_Hawk_fit[i,]
est_intensity_online <- uniHawkesIntensityNumeric(lambda0 = Hom_Hawk_fit$Baseline_online[i],
                                                  alpha = Hom_Hawk_fit$alpha_online[i],
                                                  beta = Hom_Hawk_fit$beta_online[i],
                                                  events = c(unlist(Hom_Hawk_fit$train_times[i]),
                                                                    unlist(Hom_Hawk_fit$test_times[i])))
est_intensity_full <- uniHawkesIntensityNumeric(lambda0 = Hom_Hawk_fit$Baseline_full[i],
                                                alpha = Hom_Hawk_fit$alpha_full[i],
                                                beta = Hom_Hawk_fit$beta_full[i],
                                                events = c(unlist(Hom_Hawk_fit$train_times[i]),
                                                           unlist(Hom_Hawk_fit$test_times[i])))

max(est_intensity_online$lambda.t)
max(est_intensity_full$lambda.t)
plot(est_intensity_online$time.vec,est_intensity_online$lambda.t,typ='l')
events = c(unlist(Hom_Hawk_fit$train_times[i]),
         unlist(Hom_Hawk_fit$test_times[i]))
points(events,rep(0.1,length(events)),col="red")

plot(est_intensity_full$time.vec,est_intensity_full$lambda.t,typ='l')


online_intens =as_tibble(est_intensity_online) %>% mutate(int = lambda.t/max(lambda.t),
                                                          model= "online") %>%
  select(time.vec,int,model)
full_intens = as_tibble(est_intensity_full) %>% mutate(int = lambda.t/max(lambda.t),
                                                       model = "full") %>% 
  select(time.vec,int,model)

total_intens = rbind(online_intens,full_intens)

as_tibble(events,height=1)
events_df = tibble::enframe(events)
events_df$y = 0

total_intens %>% ggplot(aes(time.vec,int)) + geom_line(aes(colour=model)) +
  scale_color_manual(name = "Procedure",values=c("red","black")) +
  xlab("Time") + ylab("Intensity") +
  theme(axis.title=element_text(size=14),axis.text.y=element_blank(),
        axis.ticks.y = element_blank() ) + 
  geom_point(data = events_df,aes(value,y)) #+
  #geom_hline(yintercept = Hom_Pois_fit$Baseline_online[i],colour="blue") +
  #geom_hline(yintercept = Hom_Pois_fit$Baseline_batch[i],colour="green")

axis.ticks.x=element_blank()

### then repeat this for Inhom Hawkes
