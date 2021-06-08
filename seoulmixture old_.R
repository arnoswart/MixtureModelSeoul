library( tidyverse )
library( patchwork )
library( here )
library( readxl )
library( openxlsx )
library( rstan )
library( tidybayes )
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

ggplot( df_controls, aes(x=OD_ctrl, y=OD_avg) ) +
  geom_line( ) +
  geom_point( aes( color=control_id ) ) +
  geom_abline( data=df_p2p, aes(slope=slope, intercept=0 ), alpha=0.5) +
  facet_wrap( vars(plate_id) )
ggsave( "linearmodel.png" )

# OD_avg = slope * OD_ctrl
df_controls %>% left_join( df_p2p, by="plate_id" ) %>%
  mutate( OD_corr= OD_ctrl*slope  ) %>%
  pivot_longer( c(OD_corr, OD_ctrl), names_to="is_corr", values_to="OD" ) %>%
  ggplot( ) +
  geom_line( aes( x=plate_id, y=OD, 
                  color=control_id, group=interaction( control_id, is_corr ),
                  linetype=is_corr) )+
  #geom_text_repel( data = . %>% group_by(control_id, .groups="drop") %>%
                    #   summarize( avg = mean(OD) ),
                    # aes(y=avg,label = control_id), x=14,
                    # size = 4,
                    # nudge_x = 1,
                    # segment.color = NA ) +
  theme( axis.text.x = element_text( angle=45)) +
  scale_x_discrete( "Plate ID", labels= str_c( "plate ", 1:11)) +
  scale_y_continuous( "log-OD" ) +
  scale_color_discrete( guide=FALSE ) +
  scale_linetype_discrete( "Corrected", labels=c("Yes", "No"))

ggsave( "controls.png" )

#
# Generate Table 1, study design
#

df_OD %>% group_by( study, sample_matrix ) %>%
  summarize( n_rats = n(), .groups="drop" ) %>%
  pivot_wider( names_from="study", values_from="n_rats", values_fill=0) %>% 
  gt(rowname_col = "sample_matrix") %>%
  tab_stubhead(label = "Sample matrix") %>%
  cols_label(
    captive = md("Captive 2018"),
    feeder = md("Feeder 2016"),
    wildrats = md("Wildrats 2018")
  ) %>% 
  as_latex() %>% 
  as.character() %>% 
  str_replace_all( "longtable", "tabular") %>% 
  cat( file="./paper/table1.tex")
  
#
# Number of censored observations, needed in model
#
n_censored <- df_OD %>% 
  filter( OD <=0 ) %>% 
  group_by( sample_matrix, study ) %>%
  summarize( n_censored = n(), .groups="drop" ) %>% 
  pull(n_censored)


#
# Apply plate to plate correction
# 
df_OD <- df_OD %>% 
  apply_p2p( df_p2p ) %>%
  filter( !(study=="feeder" & sample_matrix=="heart fluid"))

df_OD %>% 
  pivot_longer( cols=c(logOD, logOD_corr), names_to="corr", values_to="ODval") %>%
  mutate( corr=fct_recode( corr, "No"="logOD", "Yes"="logOD_corr")) %>% 
  ggplot( ) + 
    geom_histogram( aes( ODval, fill=corr ), position="dodge" ) +
    scale_x_continuous( "Log10( ODvalue)" ) +
    scale_fill_discrete( "Correction") +
    facet_wrap( vars(study))
ggsave( "studies.png" )

#
# Classical approach
#

# Take mean + 2xSD of wild rat population
df_cutoffs <- df_OD %>% 
  filter( study=="wildrats" ) %>% 
  summarize( mean = mean(logOD_corr, na.rm=T), 
             sd   = sd(logOD_corr, na.rm=T),
             mean_linear = mean(10^logOD_corr), 
             sd_linear   = sd(10^logOD_corr )) %>%
  mutate( cutoff_2sd_linear = log10(mean_linear+2*sd_linear),
          cutoff_3sd_linear = log10(mean_linear+3*sd_linear),
          cutoff_2sd = mean+2*sd, 
          cutoff_3sd = mean+3*sd ) %>%
  pivot_longer( starts_with( "cutoff" ), names_to="cutoff" ) %>%
  select( -mean, -sd, -mean_linear, -sd_linear )

#
# Bayesian Binary mixture model
#
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

initfn <- function( ){ 
  return( list( p_study=c(0.05, 0.7,0.01), mu_raw=c(-1.5,0), mu_serum=0.2 ))}

m_seoul <- stan(file="seoulmixture.stan", 
                init=initfn, iter=3000, 
                data = compose_data(df_OD, n_components=2, n_censored=n_censored ) )

stan_trace( m_seoul, pars =c("mu_raw[1]", "mu_serum", "p_study[1]" ) )
ggsave( "traceplot.png")

df_seoul <- m_seoul %>%
  recover_types(df_OD) %>%
  spread_draws(  mu_raw[component], sigma[component], mu_serum,  p_study[study] ) %>%
  ungroup()

df_seoul %>% 
  mutate( component=factor(component, levels=c(1,2), labels=c("Negative", "Positive"))) %>% 
  group_by( component, study ) %>%
  mean_qi( mu_raw, mu_serum, sigma, p_study ) %>%
  mutate( mu   = sprintf( "%1.2f (%1.2f,%1.2f)", mu_raw, mu_raw.lower, mu_raw.upper),
          mu_serum = sprintf( "%1.2f (%1.2f,%1.2f)", mu_serum, mu_serum.lower, mu_serum.upper),
          sigma = sprintf( "%1.2f (%1.2f,%1.2f)", sigma, sigma.lower, sigma.upper),
          p_study = sprintf( "%1.2f (%1.2f,%1.2f)", p_study, p_study.lower, p_study.upper)
         ) %>%
  select( study, mu, mu_serum, sigma, p_study ) %>% 
  gt() %>%
  cols_label(
    study = md("Study"),
    mu = md("mu"),
    mu_serum = md("muserum"),
    sigma = md("sigma"),
    p_study = md("pstudy")) %>% 
  as_latex() %>% 
  as.character() %>% 
  str_replace_all( "longtable", "tabular") %>% 
  cat( file="./paper/table2.tex")

p1 <- ggplot( df_seoul ) + stat_halfeye( aes( x=mu_serum )) +
  scale_x_continuous( "mu_matrixtype") +
  scale_y_continuous( "density")
p2 <- ggplot( df_seoul ) +  stat_halfeye( aes( x=mu_raw, fill=as.factor(component) ))+
  scale_x_continuous( "mu_heartfluid") +
  scale_y_continuous( "density")
p3 <- ggplot( df_seoul) + stat_halfeye( aes( x=p_study, fill=study)) + 
  scale_x_continuous( "prevalence") +
  scale_y_continuous( "density")

((p1 / p2 ) | p3) + 
  plot_annotation(tag_levels = 'A')
ggsave("parameters.png")

plotfit(df_OD, df_seoul, df_cutoffs )
ggsave( "seoulfit.png" )                   
          
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
  filter( study=="captive") %>%
  group_by( sample) %>%
  summarize( logOD_corr=mean(logOD_corr) ) %>%
  expand_grid( df_cutoffs ) %>%
  group_by( cutoff ) %>%
  summarize(   total = n(),
               pos   = sum(logOD_corr > value )) %>%
  mutate( prev  = ifelse( pos==0, 0,  pos / total  ) )

df_OD %>%
  left_join( param_mean ) %>%
  mutate( p.pos = p_study * dnorm( logOD_corr, mu_2, sigma_2) /
            mix( logOD_corr, mu_1, mu_2, sigma_1, sigma_2, p_study)) %>%
  group_by( study) %>%
  mean_qi( p.pos, .width=.99 )

df_seoul_pred <- df_OD %>% 
  left_join( param_mean ) %>%
  mutate( p.pos = p_study * dnorm( logOD_corr, mu_2, sigma_2) / 
                            mix( logOD_corr, mu_1, mu_2, sigma_1, sigma_2, p_study)) %>%
  group_by( sample, study, p_study ) %>%
  summarize( p.pos = mean(p.pos),
             logOD_corr=mean(logOD_corr),
             across( starts_with("is_"), first)) %>%
  mutate( is_PCR_pos = factor( is_PCR_pos, levels=c("pos", "neg"), labels=c("pos", "neg")),
          is_mixture_pos   = ifelse(p.pos>0.5, "pos", "neg" ) %>% 
            factor( levels=c("pos", "neg"), labels=c("pos", "neg")),
          is_cutoff2_linear_pos    = ifelse(logOD_corr > df_cutoffs[1, 2], "pos", "neg" ) %>% 
            factor( levels=c("pos", "neg"), labels=c("pos", "neg")),
          is_cutoff3_linear_pos    = ifelse(logOD_corr > df_cutoffs[2, 2], "pos", "neg" ) %>% 
            factor( levels=c("pos", "neg"), labels=c("pos", "neg")),
          is_cutoff2_pos    = ifelse(logOD_corr > df_cutoffs[3, 2], "pos", "neg" ) %>% 
            factor( levels=c("pos", "neg"), labels=c("pos", "neg")),
          is_cutoff3_pos    = ifelse(logOD_corr > df_cutoffs[4, 2], "pos", "neg" ) %>% 
            factor( levels=c("pos", "neg"), labels=c("pos", "neg"))
          ) %>%
  mutate( is_vnt_pos = ifelse( is.na(is_vnt_pos), "-", as.character(is_vnt_pos) ) ) # for plotting
write.xlsx( df_seoul_pred, file = "seoul_predictions.xlsx" )

ggplot( df_seoul_pred %>% filter( !is.na(is_PCR_pos))  ) +
  geom_boxplot( aes( x=is_PCR_pos, y=p.pos ), outlier.color=NA ) +
  geom_jitter( aes( x=is_PCR_pos, y=p.pos, color=study ), width=0.1, height=0.01, size=2  )+
  scale_y_continuous( "Probability of seropositivity", limits=c(0,1)) +
  scale_x_discrete( "PCR outcome" ) +
  facet_wrap( vars(study))
ggsave( "pcr_vs_serology.png")




df_OD %>% 
  left_join( param_mean ) %>%
  mutate( p.pos = p_study * dnorm( logOD_corr, mu_2, sigma_2) / 
            mix( logOD_corr, mu_1, mu_2, sigma_1, sigma_2, p_study)) %>%
  group_by( sample, study ) %>%
  summarize( p.pos = mean(p.pos),
             logOD_corr=mean(logOD_corr)) %>%
  mutate( positive_at = case_when(
    logOD_corr > df_cutoffs[4,2] ~ "Cutoff 3SD",
    logOD_corr > df_cutoffs[3,2] ~ "Cutoff 2SD",
    logOD_corr > df_cutoffs[2,2] ~ "Cutoff 3SD (linear)",
    logOD_corr > df_cutoffs[1,2] ~ "Cutoff 2SD (linear)",
    TRUE ~ "Negative"
  ) %>% factor( levels=c("Negative", "Cutoff 2SD (linear)","Cutoff 3SD (linear)","Cutoff 2SD", "Cutoff 3SD")) ) %>% 
  ggplot() +
  #geom_histogram( aes( x=p.pos, fill=positive_at), position="stack", bins=50) +
  geom_jitter(aes( x=positive_at, color=positive_at, y=p.pos), width = 0.2, height = 0 ) +
  #geom_boxplot( aes( x=positive_at, fill=positive_at, y=p.pos) ) +
  geom_hline( yintercept=0.5 ) +
  facet_grid( rows=vars(study), scales = "free_y" ) +
  scale_fill_viridis_d()
ggsave( "cutoff_vs_mixture.png" )

# Statistics for paper
#xtabs( ~ is_mixture_pos + is_PCR_pos + study, data=df_seoul_pred ) %>% as.data.frame() %>%
#  kable( format="latex" ) %>% cat( file="./paper/table2.tex")

##FIX!
# xtabs( ~ is_cutoff1_pos + is_PCR_pos +study, data=df_seoul_pred ) %>%
#   kable( format="latex") %>% cat( file="./paper/table2.5.tex")
# 
# xtabs( ~ is_cutoff2_pos + is_PCR_pos +study, data=df_seoul_pred ) %>%
#   kable( format="latex") %>% cat( file="./paper/table3.tex")
# 
# xtabs( ~ is_cutoff3_pos + is_PCR_pos +study, data=df_seoul_pred ) %>%
#   kable( format="latex") %>% cat( file="./paper/table4.tex")

df_seoul_pred %>% 
  filter( study!="wildrats") %>% # No PCR
  droplevels() %>% 
  group_by( study ) %>%  
  group_modify(
    ~with( .,
      rbind(
        confusionMatrix( table( is_mixture_pos, is_PCR_pos ) ) %>% 
          tidy() %>%
          select(term, estimate) %>%
          pivot_wider(names_from=term, values_from=estimate ) %>%
          mutate( test="Binary Mixture") %>%
          select( test, accuracy, sensitivity, specificity, PPV=pos_pred_value, NPV=neg_pred_value ),
        confusionMatrix( table( is_cutoff2_linear_pos , is_PCR_pos ) ) %>% 
          tidy() %>%
          select(term, estimate) %>%
          pivot_wider(names_from=term, values_from=estimate ) %>%
          mutate( test="Cutoff2_linear") %>%
          select( test, accuracy, sensitivity, specificity, PPV=pos_pred_value, NPV=neg_pred_value ),
        confusionMatrix( table( is_cutoff3_linear_pos , is_PCR_pos ) ) %>% 
          tidy() %>%
          select(term, estimate) %>%
          pivot_wider(names_from=term, values_from=estimate ) %>%
          mutate( test="Cutoff3_linear") %>%
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
          mutate( test="Cutoff3") %>%
          select( test, accuracy, sensitivity, specificity, PPV=pos_pred_value, NPV=neg_pred_value )
        ))) %>%
  as_tibble() %>% 
  mutate( study=factor( study, levels=c("captive", "feeder"), labels=c("Captive study", "Feeder study"))) %>% 
  group_by(study) %>%
  gt( rowname_col="test") %>% 
  fmt_number(
    columns = 3:7,
    decimals = 2
  ) %>% 
  as_latex() %>% 
  as.character() %>% 
  str_replace_all( "longtable", "tabular") %>% 
  str_replace_all( "_", " ") %>% 
  cat( file="./paper/table3.tex")

 
df_seoul_pred %>% 
  filter( study=="Captive" ) %>% with(.,
    confusionMatrix( table( is_PCR_pos, is_cutoff2_pos ), positive="pos" ) %>%
    tidy() %>%
    select(term, estimate) %>%
    pivot_wider(names_from=term, values_from=estimate ) %>%
    mutate( test="Cutoff2") %>%
    select( test, accuracy, sensitivity, specificity, pos_pred_value, neg_pred_value, prevalence ))

df_seoul_pred %>% filter( is_mixture_pos != is_PCR_pos, study=="captive" ) %>% ungroup() %>% 
  select( starts_with("is_"))

df_seoul_pred %>% filter( is_vnt_pos != "-") %>% 
  with(.,
                       xtabs( ~is_PCR_pos+is_vnt_pos+study))

df_seoul_pred %>% filter( study=="wildrats" )



df_seoul_pred %>% filter( study=="captive") %>% 
  mutate( across( starts_with("is_"),  function(x) x=="pos" )) %>% 
  filter( !(is_PCR_pos | is_mixture_pos | is_cutoff3_linear_pos) ) %>% 
  nrow()

library( ggvenn )
ggvenn(
  df_seoul_pred %>% filter( study=="captive") %>% 
    mutate( across( starts_with("is_"),  function(x) x=="pos" )) %>% 
    select( PCR=is_PCR_pos, mixture=is_mixture_pos, `linear cutoff 3SD`=is_cutoff3_linear_pos),
  columns=c( "PCR", "mixture", "linear cutoff 3SD" ),
  show_percentage = FALSE
)
ggsave("venn.png")
