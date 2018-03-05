sparsePCARef <- function(O, K, t) {

num_components = K;

#RUN PCA
pcs = prcomp(scale(t(O)));

coeff = pcs$rotation
score = pcs$x

x = score[,1:K]%*%t(coeff[,1:K]);
An = scale(t(O),center=T,scale=F)
Bn = scale(x,center=T,scale=F)
An = t(t(An)*(1/sqrt(apply(An^2,2,sum))))
Bn = t(t(Bn)*(1/sqrt(apply(Bn^2,2,sum))))


# Find the distance of each site from its low rank approximation.
distances = apply((An-Bn)^2,2,sum)^0.5 ;
dsort = sort(distances,index.return=T);
ranked_list = dsort$ix

#run PCA on selected sites
sites = ranked_list[1:t];
pcs = prcomp(scale(t(O[sites,])));
score = pcs$x


SPCs = t(score[,1:num_components])

return(SPCs)
}

sparsePCANew <- function(O, K, t) {

num_components = K;
step_size = -5
nsites = nrow(O)
ksites = 1:nsites 
Ot = O
#start with a big step size and get smaller
for (iter in seq(nsites-1,t,step_size)) {
#cat('running iter ',iter,'\n')
#RUN PCA on remaining sites
pcs = prcomp(scale(t(Ot)))

coeff = pcs$rotation
score = pcs$x

x = score[,1:K]%*%t(coeff[,1:K]);
An = scale(t(Ot),center=T,scale=F)
Bn = scale(x,center=T,scale=F)
An = t(t(An)*(1/sqrt(apply(An^2,2,sum))))
Bn = t(t(Bn)*(1/sqrt(apply(Bn^2,2,sum))))

# Find the distance of each site from its low rank approximation.
distances = apply((An-Bn)^2,2,sum)^0.5 ;
dsort = sort(distances,index.return=T);
ranked_list = dsort$ix
ksites = ranked_list[1:iter]
Op = Ot
Ot = Ot[ksites,]
}
Ot = Op[ranked_list[1:t],]
pcs = prcomp(scale(t(Ot)));
score = pcs$x


SPCs = t(score[,1:num_components])
return(SPCs)
}


#test 
N = 1*10^2 #number of individuals
M = 10^3 #number of sites
t = 5*10^1 #number of significant sites
K = 2 #number of cell types
sig = 1 #sd difference of cell types at each significant site

regCor = refCor = desCor = 0
for (it in 1:100) {
#assume a different mean at all t columns for all K cell types
means = matrix(rnorm(K*t,0,sig),nrow=K,ncol=t)

#cell type proporitions of each individual
mixture = matrix(runif(N*K),nrow=N,ncol=K)
csum = apply(mixture,1,sum)
mixp = mixture/csum 

#sample means of signficant cols
smeans = mixp%*%means

signal = matrix(0,nrow=N,ncol=t)
for(i in 1:N) {
      signal[i,] = rnorm(t,smeans[i,],1)
}
noise = matrix(rnorm(N*(M-t)),nrow=N,ncol=(M-t))

data = cbind(signal,noise)

#RUN PCA
pcs = prcomp(scale(data));
coeff = pcs$rotation #loadings
score = pcs$x #pcs   . . . scale(data) = score%*%t(coeff)
regular_pcs = pcs$x

#refactor original
refactor_pcs = t(sparsePCARef(t(data),K,t))

#refactor original
descent_pcs = t(sparsePCANew(t(data),K,t))

#measure correlation of pcs to true proportion of cell types
regCor[it] = cor(regular_pcs[,1],smeans[,1])^2
refCor[it] = cor(refactor_pcs[,1],smeans[,1])^2
desCor[it] = cor(descent_pcs[,1],smeans[,1])^2

cat(mean(regCor),mean(refCor),mean(desCor),'\n')
}