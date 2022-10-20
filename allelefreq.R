

## 1. the likelihood function -----
likFunction <- function (p, epison, df){
  llik = 0
  for (i in 1:dim(df)[2]){
    rinfo = as.character(df[1,i])
    rinfo = unlist(strsplit(rinfo, split="[/,]"))
    r1 = as.numeric(rinfo[2]) ## switch. The second number is derived.
    r2 = as.numeric(rinfo[1])
    #print (r1)
    t = r1+r2
    lp = p^2 * dbinom(r1,size=t,prob=1-epison) + 2*p*(1-p) * dbinom(r1,size=t,prob=0.5) + (1-p)^2 * dbinom(r1,size=t,prob=epison)
    llik = llik + log(lp)
  }
  return(llik)
}


## 2. the allele frequency function ----



alleleFreq <- function (rdata, epison){
  
  read = rdata
  read = as.data.frame(read)
  opt = optimize(likFunction, epison=epison, df=read, interval=c(0, 1), maximum=T)
  estimate = round(opt$maximum,2)
  
  
  ## start to get the confidence interval
  
  grid = seq(0,1,by=0.001)
  
  gr = c(0,0)
  read = rdata
  
  read = as.data.frame(read)
  
  for (g in grid){
    
    lout = likFunction(p=g, epison=0.01, df=read)
    lout = c(g, lout)
    gr = rbind(gr,lout)
    
  }
  gr = gr[2:dim(gr)[1],]
  gr.row = dim(gr)[1]
  max.index = which.max(gr[,2])
  p.hat = gr[max.index,1]
  
  l.index = which.min(abs(gr[1:max.index,2]-gr[max.index,2]+1.92))  ## 1.92-log likelihoods lower bound index
  u.index = which.min(abs(gr[max.index:gr.row,2]-gr[max.index,2]+1.92))  ## 1.92-log likelihoods upper bound index
  u.index = u.index + max.index - 1
  ind.estimate = c(round(estimate,2), round(gr[l.index,1],2), round(gr[u.index,1],2))
  
  names(ind.estimate) = c('estimate', 'lower', 'upper')
  ind.estimate
  
}


## 3. run with an example----

n_anc_read_ind1 = 9  ## ancestral read count of individual 1
n_der_read_ind1 = 1  ## derived read count of individual 1

n_anc_read_ind2 = 15 ## ancestral read count of individual 2
n_der_read_ind2 = 14 ## derived read count of individual 2

input.data = c(paste(n_anc_read_ind1,n_der_read_ind1, sep='/'), paste(n_anc_read_ind2,n_der_read_ind2, sep='/'))

input.data = as.data.frame(t(input.data))

error.rate = 0.001

print (alleleFreq(rdata = input.data, epison = error.rate))


