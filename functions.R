
study_labels <- c('captive'="Captive",'feeder'="Feeder",'wildrats'="Wild rats")

my_theme <-   theme( axis.text = element_text(size=10),
                     axis.title = element_text(size=11),
                     legend.title = element_text(size=10),
                     legend.text = element_text(size=8))
my_theme2 <-   theme( axis.text = element_text(size=15),
                     axis.title = element_text(size=16),
                     legend.title = element_text(size=15),
                     legend.text = element_text(size=13))

inv_logit <- function( x ){
  1/(1+exp(-x))
}

# De 2018 uitslagen (behalve wilde ratten) zijn de huidige studie.
# De 2016 uitslagen hadden we later toegevoegd, omdat die positieven bevatten die we nodig hebben.
# 2016 controls: 16-2121 16-2123 16-2124 16-2129 16-2122
# 2016 let op: OD's zijn gemiddelden
c_controls <- c("PG", "RG", "16-2126", "16-2127", "16-2129",
                "16-2121", "16-2123", "16-2124", "16-2129", "16-2122" )

read_data <- function(){
  data <- rbind(
    read_xlsx( "./data/ELISAoverzichtjuni2019_v2.xlsx", range=cell_cols("A:D"), sheet="data_compleet" ),
    read_xlsx( "./data/ELISAoverzichtjuni2019_v2.xlsx", range=cell_cols("A:C"), sheet="2016 Arno") %>% 
      mutate( OD2=NA_real_) ) %>%
    filter( !(sample %in% c("18-2648", "18-2649"))) # Email 22-12-2020
}

read_vnt <- function(){
  df_vnt <- read_xlsx("./data/Copy of seoul_predictions_with_cutoff_wild_rats_17-12-2020_MM22_12.xlsx") %>% 
    select( sample, is_vnt_pos ) %>%
    filter( is_vnt_pos !="-" ) %>% 
    rbind( tribble( ~sample, ~is_vnt_pos,
                    "16-2127", "pos",
                    "16-2128", "pos",
                    "18-2355", "pos")) %>%
    mutate( is_vnt_pos = ifelse( sample=="16-2179", "pos", is_vnt_pos) ) %>% 
    mutate( is_vnt_pos=as.factor(is_vnt_pos))

}

read_controls <- function(){
  
  data <- read_data()
  
  data_controls <- data %>%
    filter( sample %in% c_controls, sample !="N" ) %>%
    select( plate_id=plaat, control_id=sample, OD1, OD2, sample ) %>%
    pivot_longer( starts_with("OD"), names_to="replicate", values_to="OD_ctrl" ) %>%
    select( -replicate ) %>%
    mutate( plate_id=as.factor(plate_id), control_id=as.factor(control_id) ) %>%
    filter( !is.na( OD_ctrl )) %>% 
    filter( OD_ctrl > 0 ) %>%
    mutate( OD_ctrl = log10(OD_ctrl)) %>%
    group_by( control_id) %>%
    mutate(OD_avg = mean( OD_ctrl, na.rm=T) ) %>%
    ungroup()
}

read_data_OD <- function(){
  data <- read_data()
  
  data_serumtype <- rbind(
    read_xlsx( "./data/ELISAoverzichtjuni2019_v2.xlsx", sheet="Arno" ) %>%
      select( sample, serumtype=starts_with("serum") ),
    # 29-1-2021 Miriam finds new wild rats, but in the end we drop them, since
    # plate information is lost
    read_xlsx( "./data/Copy of ELISAoverzichtjuni2019_v2_MM 2901.xlsx", sheet=1 ) %>%
      select( sample, serumtype=starts_with("serum") ) ) %>%
    filter( serumtype!="N" ) %>%
    mutate( serumtype=factor( serumtype, labels=c("heart fluid", "serum") ) )
  
  data_PCR <- read_xlsx( "./data/ELISAoverzichtjuni2019_v2.xlsx", range=cell_cols("A:D"), sheet=1 ) %>%
    select( sample, is_PCR_pos=starts_with("PCR") ) %>%
    filter( is_PCR_pos!="ND") 
  
  data_vnt <- read_vnt()
  
  data_OD <- data %>%
    filter( !(sample %in% c_controls ), !sample=="N" ) %>%
    left_join( data_serumtype, by="sample" ) %>%

    select( sample,plate=plaat, OD1, OD2, serumtype ) %>%
    pivot_longer( starts_with("OD"), names_to="replicate", values_to="OD" ) %>%
    select( -replicate )  %>%
    mutate( plate=as.factor(plate)) %>%
    mutate( study = case_when(
      str_detect( sample, "18-22.*") ~ "wildrats",
      str_detect( sample, "16-.*")   ~ "feeder",
      str_detect( sample, "17-.*")   ~ "feeder",
      T                              ~ "captive" ) %>%
        as.factor() )  %>%   
    left_join( data_PCR, by="sample" ) %>%
    left_join( data_vnt, by="sample" ) %>%
    rename( sample_matrix=serumtype ) %>%
    filter( !(sample %in% c("18-2648", "18-2649"))) %>%  # Email 22-12-2020
    mutate( is_PCR_pos = ifelse( study=="wildrats", "neg", is_PCR_pos )) %>% 
    mutate( is_PCR_pos = as.factor( is_PCR_pos ))
}

calc_p2p <- function( data_controls ){
  data_controls  %>%
    group_by( plate_id ) %>%
    group_modify(
      ~lm( OD_avg ~ OD_ctrl -1 , data=. ) %>% tidy() ) %>%
    select( plate_id, term, estimate ) %>%
    pivot_wider( names_from="term", values_from=estimate ) %>%
    select( plate_id, slope=`OD_ctrl` ) %>%
    ungroup()
}

apply_p2p <- function( df_OD, df_p2p ){
  df_OD %>% 
    left_join( df_p2p %>% select( plate=plate_id, slope ), by="plate") %>%
    filter( OD > 0 ) %>% 
    mutate( logOD_corr = log10(OD)*slope,
            logOD = log10(OD) ) %>% 
    mutate( logOD = ifelse( is.infinite(logOD), NA, logOD),
            logOD_corr = ifelse( is.infinite(logOD_corr), NA, logOD_corr))
}

mix <- function( x,mu1,mu2,sigma1,sigma2, p, delta=0 ){
  return( (1-p)*dnorm(x,mu1+delta,sigma1) + p*dnorm(x,mu2+delta,sigma2) )
}

initfn <- function( ){ 
  return( list( p_study=c(0.05, 0.7,0.01), 
                mu_raw=c(-1.5,0), mu_serum=0.2 ))}

plotfit <- function(data_OD, df_seoul, df_cutoffs ){
  
  df_cutoffs <- df_cutoffs %>% 
    mutate( y=3- (1:n())/(1.7*n()))
  
  p <- ggplot( data_OD )
  p <- p + geom_vline( data=df_cutoffs, aes( xintercept=value) )
  p <- p + geom_histogram( aes( logOD_corr, y=..density.., fill=study ), 
                           bins=30, position="dodge" )
 
  tmp <- df_seoul %>% 
    pivot_wider( names_from=component, values_from=c(mu_raw, sigma) ) 
  
  for( i in 1:50){
    p <- p + stat_function( data=tibble( x=c(-2,1)), aes(x), fun=mix, 
                            args=as.list(tmp %>% filter( study=="captive") %>% 
                                           select( mu1=`mu_raw_1`, mu2=`mu_raw_2`, sigma1=`sigma_1`, sigma2=`sigma_2`,
                                                   p=`p_study` ) %>% 
                                           sample_n(1)),
                            alpha=0.1, color="red" )
    
    p <- p +  stat_function( data=tibble( x=c(-2,1)), aes(x), fun=mix, 
                             args=as.list(tmp %>% filter( study=="feeder") %>%
                                          select( mu1=`mu_raw_1`, mu2=`mu_raw_2`, sigma1=`sigma_1`, sigma2=`sigma_2`,
                                                  p=`p_study`, delta=`mu_serum`  ) %>% 
                                            sample_n(1)),
                             alpha=0.1, color="green" )
    
    p <- p +  stat_function( data=tibble( x=c(-2,1)), aes(x), fun=mix, 
                             args=as.list(tmp %>% filter( study=="wildrats") %>%
                                            select( mu1=`mu_raw_1`, mu2=`mu_raw_2`, sigma1=`sigma_1`, sigma2=`sigma_2`,
                                                    p=`p_study` ) %>% 
                                            sample_n(1)),
                             alpha=0.1, color="blue" )
  }
  p <- p + geom_label( data=df_cutoffs, aes( x=value, y=y, label=cutoff_label),
                       label.padding = unit(0.1, "lines"), size=3 )
  p <- p + scale_y_continuous( "Density", limits=c(0,3) )
  p <- p + scale_x_continuous( "LogOD-value" )
  p <- p + scale_fill_discrete( "Study", labels=study_labels )
  p+my_theme
}

# See https://github.com/rstudio/gt/issues/619
remove_escape_latex <- function(x) {
  enable_special_characters = function(x) {
    gsub("\\\\([&%$#_{}])", "\\1", x, fixed = FALSE, ignore.case = TRUE) 
  }
  enable_backslash = function(x) {
    gsub("\\\\textbackslash([[:space:]])?", "\\\\", x, fixed = FALSE, ignore.case = TRUE)
  }
  enable_tilde = function(x) {
    gsub("\\\\textasciitilde([[:space:]])?", "~", x, fixed = FALSE, ignore.case = TRUE)
  }
  enable_exponents = function(x) {
    gsub("\\\\textasciicircum ", "\\^", x, fixed = FALSE, ignore.case = TRUE)
  }
  
  enable_backslash(enable_special_characters(enable_tilde(enable_exponents(x))))
}