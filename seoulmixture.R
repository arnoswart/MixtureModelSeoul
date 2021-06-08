library( tidyverse )
library( magrittr )
library( patchwork )
library( here )
library( readxl )
library( openxlsx )
library( rstan )
library( tidybayes )
library( distributional )
library( broom)
library( knitr )
library( ggrepel )
library( caret )
library( eulerr )
library( gt )

setwd( here() )
source( "functions.R" )

df_controls  <- read_controls()
df_OD        <- read_data_OD()
df_p2p       <- calc_p2p( df_controls )

# OD_avg = slope * OD_ctrl
df_controls %>% left_join( df_p2p, by="plate_id" ) %>%
  mutate( OD_corr= OD_ctrl*slope  ) %>%
  pivot_longer( c(OD_corr, OD_ctrl), names_to="is_corr", values_to="OD" ) %>%
  ggplot( ) +
  geom_line( aes( x=plate_id, y=OD, 
                  color=control_id, group=interaction( control_id, is_corr ),
                  linetype=is_corr) )+
  theme( axis.text.x = element_text( angle=45)) +
  scale_x_discrete( "Plate ID", labels=str_c("Plate ", 1:11) ) +
  scale_y_continuous( "log-OD", breaks=seq( -1.5, 0.5, by=0.5), limits=c(-1.5, 0.5),
                      sec.axis=sec_axis( trans = ~., name="OD", 
                                         breaks=seq( -1.5, 0.5, by=0.5), 
                                         labels=round(10^seq( -1.5, 0.5, by=0.5),2))) +
  scale_color_discrete( guide=FALSE ) +
  scale_linetype_discrete("Corrected", labels=c("Yes", "No"))  +
  my_theme

ggsave( "controls.png", width= 15, height=10, units="cm", dpi=600 )

df_controls %>% left_join( df_p2p, by="plate_id" ) %>%
           mutate( OD_corr= OD_ctrl*slope  ) %>%
           pivot_longer( c(OD_corr, OD_ctrl), names_to="is_corr", values_to="OD" ) %>% 
  group_by( control_id, is_corr ) %>% 
  summarize( sd = sd( OD )) %>% 
  pivot_wider( names_from=is_corr, values_from=sd ) %>% 
  select( control_id, OD_ctrl, OD_corr ) %>% 
  mutate( perc_decrease = 100*(OD_corr-OD_ctrl)/OD_ctrl ) 


#
# Generate Table 1, study design
#

df_OD %>% group_by( study, sample_matrix ) %>%
  summarize( n_rats = n_distinct(sample), .groups="drop" ) %>%
  pivot_wider( names_from="study", values_from="n_rats", values_fill=0) %>% 
  gt(rowname_col = "sample_matrix") %>%
  tab_stubhead(label = "Sample matrix") %>%
  cols_label(
    captive = md("Captive rat study"),
    feeder = md("Feeder rat study"),
    wildrats = md("Wild rat study")
  ) %>% 
  tab_spanner(label="2018", columns=vars(captive, wildrats)) %>% 
  tab_spanner(label="2016", columns=vars(feeder)) %>%
  as_latex() %>% 
  as.character() %>% 
  str_replace_all( "longtable", "tabular") %>% 
  cat( file="./revision2/table1.tex")
  
#
# Number of censored observations, needed in model
#
n_censored <- df_OD %>% 
  filter( !is.na( OD)) %>% 
  mutate( is_censored = as.numeric(OD <=0) ) %>% 
  summarize( n_censored = sum(is_censored) ) %>% 
  pull( n_censored )

#
# Apply plate to plate correction
# 
df_OD <- df_OD %>% 
  apply_p2p( df_p2p ) %>%
  filter( !(study=="feeder" & sample_matrix=="heart fluid"))

df_OD %>% group_by( study ) %>% 
  summarize( mean=mean(OD), sd=sd(OD), cutoff2=mean+2*sd,
             mean_corr=mean( 10^logOD_corr),
             sd_corr=sd( 10^logOD_corr), 
             cutoff2_corr=mean_corr+2*sd_corr,
             cutoff2_log =log10(cutoff2))

#
# Classical approach
#

# Take mean + X*SD of wild rat population
df_cutoffs <- df_OD %>% 
  filter( study=="wildrats" ) %>% 
  summarize( mean = mean(logOD_corr, na.rm=T), 
             sd   = sd(logOD_corr, na.rm=T),
             mean_linear = mean(10^logOD_corr, na.rm=T), 
             sd_linear   = sd(10^logOD_corr, na.rm=T )) %>%
  mutate( cutoff_2sd_linear = log10(mean_linear+2*sd_linear),
          cutoff_3sd_linear = log10(mean_linear+3*sd_linear),
          cutoff_2sd = mean+2*sd, 
          cutoff_3sd = mean+3*sd ) %>%
  pivot_longer( starts_with( "cutoff" ), names_to="cutoff" ) %>%
  select( -mean, -sd, -mean_linear, -sd_linear ) %>% 
  left_join(
    tribble( ~cutoff,              ~cutoff_label,
             "cutoff_2sd_linear",  "Cut-off 2SD, linear",
             "cutoff_3sd_linear",  "Cut-off 3SD, linear",
             "cutoff_2sd",         "Cut-off 2SD, log",
             "cutoff_3sd",        "Cut-off 3SD, log" ) )

df_OD %>% 
  pivot_longer( cols=c(logOD, logOD_corr), names_to="corr", values_to="ODval") %>%
  mutate( corr=fct_recode( corr, "No"="logOD", "Yes"="logOD_corr")) %>% 
  ggplot( ) + 
  geom_histogram( aes( ODval, fill=corr ), position="dodge", bins=30 ) +
  geom_vline( data=df_cutoffs, aes( xintercept=value) ) +
  geom_label( data=df_cutoffs %>% 
                mutate( y=60- 2.5*(1:n())), aes( x=value, y=y, label=cutoff_label), 
                  size=2, label.padding = unit(0.1, "lines"), nudge_x=-0.25) +
  scale_x_continuous( "Log10( ODvalue )" , limits=c(-3.2, 1.2) ) +
  scale_y_continuous( "Count" ) +
  scale_fill_discrete( "Correction") +
  facet_wrap( vars(study), labeller = as_labeller(study_labels) )   +
  my_theme

ggsave( "studies.png", width= 15, height=10, units="cm", dpi=600 )

#
# Bayesian Binary mixture model
#
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


m_seoul <- stan(file="seoulmixture.stan", 
                init=initfn, iter=3000, 
                data = compose_data(df_OD %>% filter(!is.na(logOD_corr)), 
                                    n_components=2, n_censored=n_censored ) )

stan_trace( m_seoul, pars =c("mu_raw[1]", "mu_serum", "p_study[1]" ) )+
  theme( text=element_text(size=12))
ggsave( "traceplot.png",  width= 15, height=10, units="cm", dpi=600 )

df_seoul <- m_seoul %>%
  recover_types(df_OD) %>%
  spread_draws(  mu_raw[component], sigma[component], mu_serum,  p_study[study] ) %>%
  ungroup()


#
# Prior influence
#

m_seoul_alternative <- stan(file="seoulmixture_altenative.stan", 
                init=initfn, iter=3000, 
                data = compose_data(df_OD %>% filter(!is.na(logOD_corr)), 
                                    n_components=2, n_censored=n_censored ) )

df_seoul_compare <- rbind(
  m_seoul_alternative %>%
    recover_types(df_OD) %>%
    spread_draws(  mu_raw[component], sigma[component], mu_serum,  p_study[study] ) %>%
    ungroup() %>% 
    mutate( scenario="alternative" ),
  df_seoul %>% 
    mutate( scenario="baseline" ) )

df_mu_raw <- tribble( 
  ~mu, ~sigma, ~component, ~scenario,
  0,   5,      1,          "baseline",
  0,   5,      2,          "baseline",
  -2,   3,     1,          "alternative",
  2,   3,      2,          "alternative"  )

ggplot(  ) +
  stat_halfeye(aes(mu_raw, fill=as.factor(component), y=scenario ), data=df_seoul_compare ) +
  stat_dist_halfeye( aes( y=scenario, dist = dist_normal(mu,sigma) ), alpha=0.2, data=df_mu_raw ) +
  scale_x_continuous("mu_baseline") +
  scale_fill_discrete( "Component" ) +
  my_theme2
ggsave("scenario_mu_raw.png",   width= 15, height=10, units="cm", dpi=600 )

df_sigma_raw <- tribble( 
  ~mu, ~sig, ~scenario,
  0,   5,     "baseline",
  1,   2,     "alternative"  )

ggplot(  ) +
  stat_halfeye(aes(sigma, fill=as.factor(component), y=scenario ), data=df_seoul_compare ) +
  stat_dist_halfeye( aes( y=scenario, dist = dist_normal(mu, sig) ), alpha=0.2, data=df_sigma_raw )+
  scale_x_continuous("sigma+ and sigma-")+
  scale_fill_discrete( "Component" )+
  my_theme2
ggsave("scenario_sigma.png",   width= 15, height=10, units="cm", dpi=600 )

df_mu_serum <- tribble( 
  ~mu, ~sig, ~scenario,
  0,   5,     "baseline",
  1,   2,     "alternative"  )

ggplot(  ) +
  stat_halfeye(aes( mu_serum, y=scenario ), data=df_seoul_compare ) +
  stat_dist_halfeye( aes( y=scenario, dist = dist_normal(mu, sig) ), alpha=0.2, data=df_mu_serum )+
  scale_x_continuous("mu_samplematrix")+
  my_theme2
ggsave("scenario_matrix.png", width= 15, height=10, units="cm", dpi=600 )

df_p_study <- tribble( 
  ~a, ~b, ~scenario,
  1,   1,     "baseline",
  1/2,   1/2, "alternative"  )

ggplot(  ) +
  stat_halfeye(aes( p_study, fill=study, y=scenario ), data=df_seoul_compare ) +
  stat_dist_halfeye( aes( y=scenario, dist = dist_beta(a,b) ), alpha=0.2, data=df_p_study )+
  scale_x_continuous("p_study")+
  my_theme2
ggsave("scenario_p_study.png", width= 15, height=10, units="cm", dpi=600 )


#
# Parameter Table
#
df_seoul %>% 
  mutate( component=factor(component, levels=c(1,2), labels=c("Negative", "Positive"))) %>%
  group_by( component, study ) %>%
  mean_qi( mu_raw, mu_serum, sigma, p_study ) %>%
  mutate( mu   = sprintf( "%1.2f (%1.2f,%1.2f)", mu_raw, mu_raw.lower, mu_raw.upper),
          mu_serum = sprintf( "%1.2f (%1.2f,%1.2f)", mu_serum, mu_serum.lower, mu_serum.upper),
          sigma = sprintf( "%1.2f (%1.2f,%1.2f)", sigma, sigma.lower, sigma.upper),
          p_study = sprintf( "%1.2f (%1.2f,%1.2f)", p_study, p_study.lower, p_study.upper)
         ) %>%
  ungroup() %>% 
  select( study, component, mu, mu_serum, sigma, p_study ) %>%
  pivot_longer( -c(study, component), names_to="parameter", values_to="estimate" ) %>%
  mutate( component = ifelse( parameter %in% c("p_study", "mu_serum"), "Both", as.character(component) ),
          study = ifelse( parameter %in% c( "p_study"), as.character(study), "All studies" )) %>% 
  unique() %>% 
  pivot_wider( names_from=component, values_from=estimate) %>% 
  select( study, parameter, Negative, Positive, Both ) %>% 
  mutate( parameter = case_when(
    parameter=="mu"       ~ '$\\mu^\\text{serum}$',
    parameter=="mu_serum" ~ "$\\mu^\\text{shift}$",
    parameter=="sigma"   ~ "$\\sigma$",
    parameter=="p_study" & study=="captive"  ~ "$p^\\text{captive}$",
    parameter=="p_study" & study=="wildrats"  ~ "$p^\\text{wild rats}$",
    parameter=="p_study" & study=="feeder"  ~ "$p^\\text{feeder}$",
    TRUE ~ parameter )) %>%
  select(-study) %>% 
  gt() %>%
  fmt_missing( columns=everything(), missing_text = "-") %>% 
  as_latex() %>% 
  as.character() %>% 
  remove_escape_latex() %>% 
  str_replace_all( "longtable", "tabular") %>% 
  cat( file="./paper/table2.tex")

#
# Figure parameter estimates
#
p1 <- ggplot( df_seoul ) + stat_halfeye( aes( x=mu_serum )) +
  scale_x_continuous( "mu_matrixtype") +
  scale_y_continuous( "density") +
  theme( text=element_text(size=15))
p2 <- ggplot( df_seoul ) +  stat_halfeye( aes( x=mu_raw, fill=as.factor(component) ))+
  scale_x_continuous( "mu_heartfluid") +
  scale_y_continuous( "density") +
  scale_fill_discrete("component")+
  theme( text=element_text(size=15))
p3 <- ggplot( df_seoul) + stat_halfeye( aes( x=p_study, fill=study)) + 
  scale_x_continuous( "prevalence") +
  scale_y_continuous( "density") +
  theme( text=element_text(size=15))

((p1 / p2 ) | p3) + 
  plot_annotation(tag_levels = 'A')
#ggsave("parameters.png", width= 15, height=10, units="cm", dpi=600 )
                  
# Pr( pos | OD ) = P( OD | pos ) P( pos ) / sum( P(OD) ) 
# proportional to P( OD | pos ) P(pos) + P( OD | pos ) P(pos) 
param_mean <- df_seoul %>%
  group_by(component, study ) %>%
  select( component, mu_raw,mu_serum,sigma, p_study, study ) %>%
  summarize( across( everything() , mean), .groups="drop" ) %>%
  pivot_wider( names_from=c(component), values_from=c(mu_raw, sigma) ) %>%
  mutate( mu_1 = ifelse( study=="feeder", mu_raw_1 + mu_serum, mu_raw_1 ),
          mu_2 = ifelse( study=="feeder", mu_raw_2 + mu_serum, mu_raw_2 ))

df_OD %>%
  left_join( param_mean ) %>%
  mutate( p.pos = p_study * dnorm( logOD_corr, mu_2, sigma_2) /
            mix( logOD_corr, mu_1, mu_2, sigma_1, sigma_2, p_study)) %>%
  group_by( study) %>%
  mean_qi( p.pos, .width=.99 )

test_pos <- function( logod, cutoff ){
  ifelse((logod > cutoff) | is.na(logod), "pos", "neg" ) %>% 
    factor( levels=c("pos", "neg"), labels=c("pos", "neg"))
}

df_seoul_pred <- df_OD %>% 
  left_join( param_mean ) %>%
  mutate( p.pos = p_study * dnorm( logOD_corr, mu_2, sigma_2) / 
                            mix( logOD_corr, mu_1, mu_2, sigma_1, sigma_2, p_study)) %>%
  group_by( sample, study, p_study ) %>%
  summarize( p.pos = mean(p.pos),
             logOD_corr=mean(logOD_corr),
             across( starts_with("is_"), first)) %>%
  mutate( is_PCR_pos = factor( is_PCR_pos, levels=c("pos", "neg"), labels=c("pos", "neg")),
          is_mixture_pos           = test_pos( p.pos, 0.5),
          is_cutoff2_linear_pos    = test_pos(logOD_corr, df_cutoffs[1, 2]),
          is_cutoff3_linear_pos    = test_pos(logOD_corr, df_cutoffs[2, 2]),
          is_cutoff2_pos           = test_pos(logOD_corr, df_cutoffs[3, 2]),
          is_cutoff3_pos           = test_pos(logOD_corr, df_cutoffs[4, 2] ))  %>%
  mutate( is_vnt_pos = ifelse( is.na(is_vnt_pos), "-", as.character(is_vnt_pos) ) ) # for plotting
write.xlsx( df_seoul_pred, file = "seoul_predictions.xlsx" )


#
# PCR-mixture boxplots
#
p1 <- ggplot( df_seoul_pred ) +
  geom_boxplot( aes( x=is_PCR_pos, y=p.pos ), outlier.color=NA ) +
  geom_jitter( aes( x=is_PCR_pos, y=p.pos, color=study ), width=0.1, height=0.01, size=0.5  )+
  scale_y_continuous( "Probability of seropositivity", limits=c(0,1)) +
  scale_x_discrete( "PCR outcome" ) +
  facet_wrap( vars(study), labeller=as_labeller( study_labels) ) +
  my_theme + 
  theme( legend.position = "none")
p2 <- ggplot( df_seoul_pred ) +
  geom_boxplot( aes( x=is_PCR_pos, y=logOD_corr ), outlier.color=NA ) +
  geom_jitter( aes( x=is_PCR_pos, y=logOD_corr, color=study ), width=0.1, height=0.01, size=0.5  )+
  scale_y_continuous( "log10( OD ) corrected", limits=c(-2,0.5) ) +
  scale_x_discrete( "PCR outcome", breaks=c("", "pos", "neg"), limits=c("", "pos", "neg")) + 
  geom_hline( data=df_cutoffs, aes( yintercept=value) ) +
  geom_label( data=df_cutoffs, 
              aes( y=value, x=1, label=cutoff_label), 
              label.padding = unit(0.1, "lines"), size=1.9, nudge_x = 0.1 ) +
  facet_wrap( vars(study), labeller=as_labeller( study_labels )) +
  theme( legend.position = "none") +
  my_theme
p2 / p1 + plot_annotation(tag_levels = 'a')
ggsave( "pcr_vs_serology.png",width= 15, height=20, units="cm", dpi=600 )


plotfit(df_OD, df_seoul, df_cutoffs )
ggsave( "seoulfit.png" ,width= 15, height=10, units="cm", dpi=600 )


df_seoul_pred %>% 
  filter( is_mixture_pos != is_PCR_pos ) %>% 
  group_by( is_mixture_pos, is_PCR_pos, study ) %>% 
  summarize( n=n(), minimum=min(p.pos), maximum=max(p.pos) ) %>% 
  arrange( study )

f <- function( x) {
  abs(dnorm( x, mu_1, sigma_1)/dnorm( x, mu_2, sigma_2) - p_study/(1-p_study)) }

df_seoul_pred %>% 
  droplevels() %>% 
  group_by( study ) %>%  
  group_modify(
    ~with( .,
      rbind(
        confusionMatrix( table( is_mixture_pos, is_PCR_pos ) ) %>% 
          tidy() %>%
          select(term, estimate) %>%
          pivot_wider(names_from=term, values_from=estimate ) %>%
          mutate( test="Binary mixture P(pos)>0.5" ) %>%
          select( test, accuracy, sensitivity, specificity, PPV=pos_pred_value, NPV=neg_pred_value ),
        confusionMatrix( table( is_cutoff2_linear_pos , is_PCR_pos ) ) %>% 
          tidy() %>%
          select(term, estimate) %>%
          pivot_wider(names_from=term, values_from=estimate ) %>%
          mutate( test="Cutoff2 linear") %>%
          select( test, accuracy, sensitivity, specificity, PPV=pos_pred_value, NPV=neg_pred_value ),
        confusionMatrix( table( is_cutoff3_linear_pos , is_PCR_pos ) ) %>% 
          tidy() %>%
          select(term, estimate) %>%
          pivot_wider(names_from=term, values_from=estimate ) %>%
          mutate( test="Cutoff3 linear") %>%
          select( test, accuracy, sensitivity, specificity, PPV=pos_pred_value, NPV=neg_pred_value ),
        confusionMatrix( table( is_cutoff2_pos , is_PCR_pos ) ) %>% 
          tidy() %>%
          select(term, estimate) %>%
          pivot_wider(names_from=term, values_from=estimate ) %>%
          mutate( test="Cutoff2") %>%
          select( test, accuracy, sensitivity, specificity, PPV=pos_pred_value, NPV=neg_pred_value ),
        confusionMatrix( table( is_cutoff3_pos, is_PCR_pos ) ) %>% 
          tidy() %>%
          select(term, estimate) %>%
          pivot_wider(names_from=term, values_from=estimate ) %>%
          mutate( test="Cutoff3", accuracy=NA, PPV=NA, NPV=NA ) %>%
          select( test, accuracy, sensitivity, specificity, PPV=pos_pred_value, NPV=neg_pred_value )
        ))) %>% 
        rbind(
          param_mean %>% 
            group_by( study ) %>% 
            mutate( 
              cutoff_mixture = optimise( 
                function( x) {
                  abs(dnorm( x, mu_1, sigma_1)/dnorm( x, mu_2, sigma_2) - p_study/(1-p_study)) }, 
                interval=c(-1,2) )$minimum,
              sensitivity = pnorm( cutoff_mixture, mu_2, sigma_2, lower.tail = FALSE ),
              specificity = pnorm( cutoff_mixture, mu_1, sigma_1, lower.tail = TRUE ),
              NPV = 1/(1+((p_study)/(1-p_study))*
                         (pnorm( cutoff_mixture, mu_2, sigma_2, lower.tail = TRUE )/
                            pnorm( cutoff_mixture, mu_1, sigma_1, lower.tail = TRUE ))),
              PPV = 1/(1+((1-p_study)/(p_study))*
                         (pnorm( cutoff_mixture, mu_1, sigma_1, lower.tail = FALSE )/
                            pnorm( cutoff_mixture, mu_2, sigma_2, lower.tail = FALSE ))),
              accuracy = p_study * pnorm( cutoff_mixture, mu_2, sigma_2, lower.tail = FALSE ) +
                (1-p_study) * pnorm( cutoff_mixture, mu_1, sigma_1, lower.tail = TRUE ),
              test = "Binary mixture theoretical" ) %>% 
            select( test, accuracy, NPV, PPV, sensitivity, specificity, study )) %>% 
  as_tibble() %>% 
  mutate( study=factor( study, levels=c("captive", "feeder", "wildrats"), 
                        labels=c("Captive study", "Feeder study", "Wild rats"))) %>% 
  filter( !(study=="Wild rats" & test=="Binary mixture theoretical" )) %>% 
  group_by(study) %>%
  arrange( test ) %>% 
  gt( rowname_col="test") %>% 
  fmt_number(
    columns = 3:7,
    decimals = 2
  ) %>% 
  fmt_missing( columns=everything(), missing_text = "-") %>% 
  as_latex() %>% 
  as.character() %>% 
  str_replace_all( "longtable", "tabular") %>% 
  str_replace_all( "_", " ") %>% 
  cat( file="./paper/table3.tex")

df_seoul_pred %>% filter( study=="feeder", is_PCR_pos=="neg" ) %>% summary()

df_seoul_pred %>% 
  select( study, is_mixture_pos, starts_with( "is_cutoff")) %>% 
  group_by( study ) %>% 
  summarize( across( .cols=starts_with("is_"), .fns=function(x) sum(x=="pos")), n=n() ) %>% 
  select( mixture=is_mixture_pos, cutoff_2SD_linear=is_cutoff2_linear_pos,
          cutoff_3SD_linear=is_cutoff3_linear_pos, cutoff_2SD=is_cutoff2_pos,
          cutoff_3SD=is_cutoff3_pos, n) %>% 
  gt( rowname_col="test") %>% 
    fmt_number(
      decimals = 0,
      columns = 2:6
  ) %>% 
  fmt_missing( columns=everything(), missing_text = "-") %>% 
  as_latex() %>% 
  as.character() %>% 
  str_replace_all( "longtable", "tabular") %>% 
  str_replace_all( "_", " ") %>% 
  cat( file="./paper/tableS1.tex")

#
# Venn diagram
#
library( ggvenn )
ggvenn(
  df_seoul_pred %>% filter( study=="captive") %>% 
    mutate( across( starts_with("is_"),  function(x) x=="pos" )) %>% 
    select( 'RT-qPCR'=is_PCR_pos, 'Mixture model'=is_mixture_pos, 
            'Cut-off'=is_cutoff3_linear_pos),
  columns=c( "RT-qPCR", "Mixture model", "Cut-off" ),
  show_percentage = FALSE
)
ggsave("venn.png")
