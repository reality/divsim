library(stringr)

# useful for sim$X3 if you're not using a normalised measure of IC or w/e
rangey <- function(x){(x-min(x))/(max(x)-min(x))}

sim <- read_csv("/home/slater/miesim/dphens/sim_no_na.lst", 
                  col_names = FALSE)

sim$X4 <- as.factor(sim$X4)
sim$X3 <- rangey(sim$X3)

ggplot(sim, aes(x=X3)) +
  geom_histogram(binwidth=.05, colour="black", fill="white")

dphens <- read_tsv("/home/slater/miesim/dphens/dphens_merged.txt", 
                  col_names = FALSE)

# Cumbersome
counts <- new.env(hash = TRUE) 
total <- 0
for(a in dphens$X2) {
  for(i in as.list(strsplit(a, ",")[[1]])) {
    if(is.null(counts[[i]])) { 
      counts[[i]] <- 0 
    }
    total <- total + 1
    counts[[i]] = counts[[i]] + 1
  }
}
a <- NULL

length(ls(counts))

weights = data.frame(iri = character(), count = integer(), stringsAsFactors = FALSE)

c = 0
for(x in ls(counts)) {
  n <- ((counts[[x]] / total) ^ 2)
  c <- c + n
  weights <- rbind(weights, data.frame(iri=x, count=n))
}

# Well, this is the simpsons diversity...
1-c

sim <- sim %>% mutate(
  HP1 = paste("http://purl.obolibrary.org/obo/", str_replace(X1, ':', '_'), sep=''),
  HP2 = paste("http://purl.obolibrary.org/obo/", str_replace(X2, ':', '_'), sep='')
)

# xi * xj * wij
sim$xijw <- apply(sim, 1, function(x) {
  s <- 1-as.numeric(x[3])
  if(as.numeric(x[3]) == 0) { # special case that there's no information at all in this comparison
    s <- 0  
  }
  if(identical(x[[1]], x[[2]])) { #special case that they are the same phenotype
    s <- 0 
  }
  return(counts[[x[5]]] * counts[[x[6]]] * (s))
  
})
# xi * xj
sim$xij <- apply(sim, 1, function(x) {
  return(counts[[x[5]]] * counts[[x[6]]])
})

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

