library(stringr)

# Cumbersome
rangey <- function(x){(x-min(x))/(max(x)-min(x))}
load_phenotype_profile_file <- function(filepath) {
  read_tsv(filepath, col_names = FALSE) 
}
load_sim_matrix_file <- function(filepath) {
  mtx <- read_csv(filepath, col_names = FALSE)
  mtx$X3 <- rangey(mtx$X3)
  mtx$X4 <- as.factor(mtx$X4)
  return(mtx)
}
get_roc <- function(mtx) {
  return(pROC::roc(as.factor(mtx$X4), mtx$X3, ci=TRUE))
}
create_profile_counts_hash <- function(phens) {
  counts <- new.env(hash = TRUE) 
  total <- 0
  for(a in phens$X2) {
    for(i in as.list(strsplit(a, ",")[[1]])) {
      if(is.null(counts[[i]])) { 
        counts[[i]] <- 0 
      }
      total <- total + 1
      counts[[i]] = counts[[i]] + 1
    }
  }
  return(counts)
}
get_simpsons_index <- function(counts) {
 weights = data.frame(iri = character(), count = integer(), stringsAsFactors = FALSE)

  c = 0
  for(x in ls(counts)) {
    n <- ((counts[[x]] / total) ^ 2)
    c <- c + n
    weights <- rbind(weights, data.frame(iri=x, count=n))
  }

  return(1-c)
}
add_full_iris <- function(mtx) {
  return(mtx %>% mutate(
    HP1 = paste("http://purl.obolibrary.org/obo/", str_replace(X1, ':', '_'), sep=''),
    HP2 = paste("http://purl.obolibrary.org/obo/", str_replace(X2, ':', '_'), sep='')
  ))
}
calculate_diversity_values <- function(mtx, counts) {
  apply(mtx, 1, function(x) {
    s <- 1-as.numeric(x[3])
    if(as.numeric(x[3]) == 0) { # special case that there's no information at all in this comparison
      s <- 0  
    }
    if(identical(x[[1]], x[[2]])) { #special case that they are the same phenotype
      s <- 0 
    }
    return(counts[[x[5]]] * counts[[x[6]]] * (s))
  })
}
calculate_abundance_values <- function(mtx, counts) {
  apply(mtx, 1, function(x) {
    return(counts[[x[5]]] * counts[[x[6]]])
  })
}


sim <- load_sim_matrix_file("/home/slater/miesim/dphens/sim_no_na.lst")
dphens <- load_phenotype_profile_file("/home/slater/miesim/dphens/dphens_merged.txt")
counts <- create_profile_counts_hash(dphens)

get_simpsons_index(counts)

sim <- add_full_iris(sim)

sim$xijw <- calculate_diversity_values(sim, counts) 
sim$xij <- calculate_abundance_values(sim, counts)

#semantic diversity before the special case ...

# average taxonomic diversity 
sum(sim$xijw) / ((total * (total-1)) / 2)
# 0.8273148

# average taxonomic distinctness
sum(sim$xijw) / sum(sim$xij)
#  0.8294299

summary(sim$xijw)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0       22       79     1742      318 38626216 

# so then we introduce the special case and now more reasonable

sum(sim$xijw) / ((total * (total-1)) / 2)
0.7474826

sum(sim$xijw) / sum(sim$xij)
0.7455765

summary(sim$xijw)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0       17       67     1570      284 16584916

# now we can look perhaps outside of two standard deviations?

ggplot(sim, aes(x=xijw)) +
  geom_histogram(binwidth=10, colour="black", fill="white")
boxplot(sim$xijw)

# now we can see that there are some million or so that are essentially contributing nothing
sim[sim$xijw==0,]

# So we're going to build 
use = data.frame(iri = character(), val = integer(), stringsAsFactors = FALSE)
for(x in ls(counts)) {
  use <- rbind(use, data.frame(iri=x, val=sum(sim[sim$HP1==x | sim$HP2 == x,]$xijw)))
}

s1 <- sim %>% group_by(HP1) %>% dplyr::summarise(h1c = sum(xijw))
s2 <- sim %>% group_by(HP2) %>% dplyr::summarise(h2c = sum(xijw))
s2$HP1 <- s2$HP2
s2$HP2 <- NULL
s3 <- left_join(s1, s2, by = "HP1")
s3 <- s3 %>% dplyr::mutate(hc = h1c+h2c, c = counts[[HP1]])

# so, we have a lots of big numbers, but mostly very small ones, and a LOT of upper outliers
boxplot(s3$hc, outline=F)

# so, how many zeros do we have?
zerobois <- s3[s3$hc==0,]
zerobois$c <- apply(zerobois, 1, function(x) {
  return(counts[[x[1]]])
})

sum(zerobois$c)
sum(zerobois$c) / total * 100

# so we can remove these from our phenotype definitions entirely
# we can removthis is kind of a trivial result - because we could hav removed these simply by looking at what would have w = 0! 
s4 <- anti_join(s3, zerobois, by="HP1")

# the intution is that we want to keep these sort of quite dominant phenotypes, they introduce a lot of splitting factors, and this is good
# to what extent can we remove our non-diverse phenotypes?
# so let's have a look at the distribution of things below the mean

# so what are high value?
mean(s4$hc)
s5 <- s4[s4$hc < mean(s4$hc),]

nrow(s5)
#5695

boxplot(s5$hc)

#high values have a high abundance and a low ic-weighted similarity (meaning they are maximally informative for telling diseases apart)
# low values have low abundance and high ic-weighted similarity (to other phenotypes)

s5$c <- apply(s5, 1, function(x) {
  return(counts[[x[1]]])
})
op <- max(s5$hc) / 100 * 5
sum(s5[s5$hc<op,]$c) / total * 100

s6 <- anti_join(s4, s5[s5$hc<op,], by="HP1")

op <- max(s5$hc) / 100 * 10
s7 <- anti_join(s4, s5[s5$hc<op,], by="HP1")

op <- max(s5$hc) / 100 * 7.5
s8 <- anti_join(s4, s5[s5$hc<op,], by="HP1")

op <- max(s5$hc) / 100 * 2.5
s9 <- anti_join(s4, s5[s5$hc<op,], by="HP1")

write_tsv(s4, "/home/slater/divsim/remove_zeros.tsv")
write_tsv(s6, "/home/slater/divsim/remove_bottom5p_belowmean.tsv")
write_tsv(s7, "/home/slater/divsim/remove_bottom10p_belowmean.tsv")
write_tsv(s8, "/home/slater/divsim/remove_bottom7.5p_belowmean.tsv")
write_tsv(s9, "/home/slater/divsim/remove_bottom2.5p_belowmean.tsv")

baseline <- load_sim_matrix_file("/home/slater/divsim/similarity/annotations.txt_dphens_baseline.txt_hp.owl.lst")
baseline_roc <- get_roc(baseline)

nozero <- load_sim_matrix_file("/home/slater/divsim/similarity/annotations.txt_dphens_merged_nozero.lst_hp.owl.lst")
nozero_roc <- get_roc(nozero)

n5pbm <- load_sim_matrix_file("/home/slater/divsim/similarity/annotations.txt_dphens_merged_no5p_belowmean.lst_hp.owl.lst")
n5pbm_roc <- get_roc(n5pbm)

n10pbm <- load_sim_matrix_file("/home/slater/divsim/similarity/annotations.txt_dphens_merged_no10p_belowmean.lst_hp.owl.lst")
n10pbm_roc <- get_roc(n10pbm)

n7p5pbm <- load_sim_matrix_file("/home/slater/divsim/similarity/annotations.txt_dphens_merged_no7.5p_belowmean.lst_hp.owl.lst")
n7p5_roc <- get_roc(n7p5pbm)

n2p5pbm <- load_sim_matrix_file("/home/slater/divsim/similarity/annotations.txt_dphens_merged_no2.5p_belowmean.lst_hp.owl.lst")
n2p5pbm_roc <- get_roc(n2p5pbm)

with_training <- load_sim_matrix_file("/home/slater/divsim/similarity/annotations.txt_dphens_with_patient_training.txt_hp.owl.lst")
with_training_roc <- get_roc(with_training)

ggroc(list(
    `LitProfileBasline (AUC=0.8734)`=baseline_roc,
    `NoZeroSemanticDiversity (AUC=0.8907)`=nz_oroc,
    `7.5SemanticDiversityCutoff (AUC=0.8998)`=n7p5_roc,
     `WithPatientTraining (AUC=0.986)`=with_training_roc
           )
      , legacy.axes = T) + labs(color = "Setting")
annotations.txt_dphens_merged_no5p_belowmean.lst_hp.owl.lst

