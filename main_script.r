f=function(x) {
  result=sqrt(1+x[1])*x[2]*cos(pi*x[2])*sin(2*pi*(x[1]-x[2]))
  #result=exp(.2*(x[1]-.5)^2)*sin(2*pi*x[2])+exp(-x[2]^2)*cos(2*pi*x[1])
  #result=sin(pi*(x[1]-x[2]))*exp(-abs(x[1])/(abs(x[2])+.01));
  return(result)
}


f0=function(t) {
  L=as.matrix(sapply(t,l));
  result=apply(L,2,f);
  return(result);
}

J=15;
cheb.rmse=cheb.max=rep(0,J-4);
leg.rmse=leg.max=rep(0,J-4);
reg.rmse=reg.max=rep(0,J-4);
trig.rmse=trig.max=rep(0,J-4)

for( m in 5:J) {

d=params[2:3];
k=function(X,Y) {
  n=dim(X)[1]; m=dim(Y)[1];
  if (n==1 & m==1) {
    return(exp(-sum(d*abs(X-Y)^1.9)));
  }
  x.1=rep(X[,1],times=1,each=m);
  x.2=rep(X[,2],times=1,each=m);
  X.new=cbind(x.1,x.2);
  y.1=rep(Y[,1],times=n,each=1);
  y.2=rep(Y[,2],times=n,each=1);
  Y.new=cbind(y.1,y.2);
  dist=abs(X.new-Y.new)^1.9;
  dist[,1]=d[1]*dist[,1]; dist[,2]=d[2]*dist[,2];
  dist.scale=apply(dist,1,sum);
  vals=exp(-dist.scale);
  result=matrix(vals,nrow=n,ncol=m,byrow=T);
  return(params[1]*result);
}

#library(orthopolynom)

cheb.base=lapply(chebyshev.c.polynomials(m,normalized=T),as.function);
var.change=function(f,scale,f.arg) {
  function(x) scale*f(f.arg(x));
}

cheb=lapply(cheb.base,function(x) var.change(x,2^(3/2),function(y) 2*y));
quad=gausschb(4*m); quad[[1]]=seq(-1,1,length.out=length(quad[[2]])); n.nodes=length(quad[[1]]);
quad=gausslg(-1,1,3*m); n.nodes=length(quad[[1]]);
#trig=trig.poly.basis(-1,1,ceiling(m/2));
leg=lapply(legendre.polynomials(m,normalized=T),as.function); #list of legendre polynomials

xg=list(); xg[["1"]]=cbind(quad[[1]],-1); xg[["2"]]=cbind(1,quad[[1]]); #goes left to right and down to up
           xg[["3"]]=cbind(quad[[1]],1); xg[["4"]]=cbind(-1,quad[[1]]);

basis=leg;

ip=function(f,g) {
  ###assume f and g are already in vector form at points quad[[1]]
  value=sum(quad[[2]]*f*g);
  return(value);
}

#diff types of covariances needed: 0 is same side, 1 is adjacent side, 2 is opposite side
#if the kernel is stationary, these can be translated to get the correct combinations
#K.bas[[0]] in reality is not needed bc when applying e-vectors you just get back the e-val
start=Sys.time()
K.bas=list(); 
for (i in 1:4) {
  for (j in 1:4) {
    name=paste(i,j,sep="");
    K.bas[[name]]=k(xg[[paste(i)]],xg[[paste(j)]]);
    if (i==j) {
      K.bas[[name]]=K.bas[[name]]+1e-10*diag(dim(K.bas[[name]])[1]);
    }
  }
}



L=length(basis);
phi=matrix(0,ncol=L,nrow=length(quad[[1]]));
for (i in 1:L) {
  phi[,i]=sapply(quad[[1]],basis[[i]]);
}
K.phi=K.bas[['11']]%*%diag(quad[[2]])%*%phi;
K.int=t(phi)%*%diag(quad[[2]])%*%K.phi;
E=eigen(K.int);

###contains eigenvectors of integral operator evaluated at quad[[1]]
V=t(E$vectors)%*%t(phi);



Cov=matrix(0,nrow=4*length(basis),ncol=4*length(basis));
sparse=0;
for (i in 1:4) {
  for (j in 1:4) {
    name=paste(i,j,sep="");
    if( abs(i-j)==2 & sparse==1) {
      print('independence');
    } else{
      Cov[((i-1)*length(basis)+1):(i*length(basis)),((j-1)*length(basis)+1):(j*length(basis))]=K.1=V%*%diag(quad[[2]])%*%K.bas[[name]]%*%diag(quad[[2]])%*%t(V);
    }
  }
}

E.covsp=eigen(Cov.sp,symmetric=T);
E.cov=eigen(Cov+1e-12*diag(4*length(basis)),symmetric=T);

t=seq(0,8,.1);
test=t(sapply(t,l)); n.test=dim(test)[1];

K.f=matrix(0,nrow=1,ncol=4*length(basis));

K.e=matrix(0,nrow=n.test,ncol=4*length(basis));
for (i in 1:4) {
  k.test=k(test,xg[[paste(i)]]);
  f.test=matrix(apply(xg[[paste(i)]],1,f),nrow=1);
  K.e[,((i-1)*length(basis)+1):(i*length(basis))]=k.test%*%diag(quad[[2]])%*%t(V);
  K.f[1,((i-1)*length(basis)+1):(i*length(basis))]=f.test%*%diag(quad[[2]])%*%t(V);
}
stop=Sys.time()

reg.bayes=seq(0,1,length.out=((2*m+1)))

reg.vals=t(sapply(reg.bayes,l));
E.reg=eigen(k(reg.vals,reg.vals),symmetric=T);
k.reg.inv=E.reg$vectors%*%diag(1/abs(E.reg$values))%*%t(E.reg$vectors);
plot(t,abs(k(test,reg.vals)%*%k.reg.inv%*%as.matrix(f0(reg.bayes))-apply(test,1,f)),col='blue',type='l');

Cov.inv=E.cov$vectors%*%diag(1/abs(E.cov$values))%*%t(E.cov$vectors);

lines(t,abs(K.e%*%Cov.inv%*%t(K.f)-f0(t)),type='l');

K.post=function(X,Y) {
  n.x=dim(X)[1]; n.y=dim(Y)[1];
  K.x=matrix(0,nrow=n.x,ncol=4*length(basis));
  K.y=matrix(0,nrow=n.y,ncol=4*length(basis));
  for (i in 1:4) {
    k.x=k(X,xg[[paste(i)]]);
    k.y=k(Y,xg[[paste(i)]]);
    K.x[,((i-1)*length(basis)+1):(i*length(basis))]=k.x%*%diag(quad[[2]])%*%t(V);
    K.y[,((i-1)*length(basis)+1):(i*length(basis))]=k.y%*%diag(quad[[2]])%*%t(V);
  }
  return(k(X,Y)-K.x%*%Cov.inv%*%t(K.y));
}

mean.post=function(X,g) {
  n.x=dim(X);
  K.g=matrix(0,nrow=1,ncol=4*length(basis));
  K.e=matrix(0,nrow=n.x,ncol=4*length(basis));
  for (i in 1:4) {
    k.x=k(X,xg[[paste(i)]]);
    g.x=matrix(apply(xg[[paste(i)]],1,g),nrow=1);
    K.e[,((i-1)*length(basis)+1):(i*length(basis))]=k.x%*%diag(quad[[2]])%*%t(V);
    K.g[1,((i-1)*length(basis)+1):(i*length(basis))]=g.x%*%diag(quad[[2]])%*%t(V);
  }
  return(K.e%*%Cov.inv%*%t(K.g));
}

N=10;




#D=2*randomLHS(n=N,k=2)-1;
reg.bayes=seq(0,8,length.out=(4*length(basis)))
n.reg=length(reg.bayes);
reg.vals=t(sapply(reg.bayes,l));
D.reg=rbind(D,reg.vals);

x1=seq(-.99,.99,length.out=20); x2=x1;
test=cbind(rep(x1,length(x2)),rep(x2,1,each=length(x1)));
t=seq(0,8,.1);
test=t(sapply(t,l)); n.test=dim(test)[1]; test=rbind(.98*test);

y=apply(D,1,f)+rnorm(N,mean=0,sd=0);
y.reg=apply(D.reg,1,f)+rnorm((N+n.reg),mean=0,sd=0);

start=Sys.time();
mu.post=mean.post(test,f);

post.mean=mu.post+K.post(test,D)%*%solve(K.post(D,D)+1e-12*diag(dim(D)[1]))%*%(y-mean.post(D,f));
stop=Sys.time();
print(stop-start);

start=Sys.time();
E.reg=eigen(k(D.reg,D.reg))#+1e-10*diag(N+n.reg));
K.reg.inv=E.reg$vectors%*%diag(1/abs(E.reg$values))%*%t(E.reg$vectors);
post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg
stop=Sys.time();
stop-start;

f.test=apply(test,1,f)
z=matrix(f.test,nrow=length(x1),ncol=length(x2),byrow=F);
plot_ly(x=~x1,y=~x2,z=~z,type='surface')

plot(post.mean,f.test);
plot(post.mean.reg,f.test)

summary(lm(f.test~post.mean-1));
summary(lm(f.test~post.mean.reg-1));

leg.max[m-4]=max(abs(f.test-post.mean));
reg.max[m-4]=max(abs(f.test-post.mean.reg));

leg.rmse[m-4]=sqrt(sum((f.test-post.mean)^2)/dim(test)[1])
reg.rmse[m-4]=sqrt(sum((f.test-post.mean.reg)^2)/dim(test)[1])


d=params[2:3];
k=function(X,Y) {
  n=dim(X)[1]; m=dim(Y)[1];
  if (n==1 & m==1) {
    return(exp(-sum(d*abs(X-Y)^1.9)));
  }
  x.1=rep(X[,1],times=1,each=m);
  x.2=rep(X[,2],times=1,each=m);
  X.new=cbind(x.1,x.2);
  y.1=rep(Y[,1],times=n,each=1);
  y.2=rep(Y[,2],times=n,each=1);
  Y.new=cbind(y.1,y.2);
  dist=abs(X.new-Y.new)^1.9;
  dist[,1]=d[1]*dist[,1]; dist[,2]=d[2]*dist[,2];
  dist.scale=apply(dist,1,sum);
  vals=exp(-dist.scale);
  result=matrix(vals,nrow=n,ncol=m,byrow=T);
  return(params[1]*result);
}

#library(orthopolynom)

cheb.base=lapply(chebyshev.c.polynomials(m,normalized=T),as.function);
var.change=function(f,scale,f.arg) {
  function(x) scale*f(f.arg(x));
}

cheb=lapply(cheb.base,function(x) var.change(x,2^(3/2),function(y) 2*y));
quad=gausschb(4*m); quad[[1]]=sort(quad[[1]]); n.nodes=length(quad[[1]]);
#quad=gausslg(-1,1,3*m); n.nodes=length(quad[[1]]);
#leg=lapply(legendre.polynomials(m,normalized=T),as.function); #list of legendre polynomials

xg=list(); xg[["1"]]=cbind(quad[[1]],-1); xg[["2"]]=cbind(1,quad[[1]]); #goes left to right and down to up
xg[["3"]]=cbind(quad[[1]],1); xg[["4"]]=cbind(-1,quad[[1]]);

basis=cheb;

ip=function(f,g) {
  ###assume f and g are already in vector form at points quad[[1]]
  value=sum(quad[[2]]*f*g);
  return(value);
}

#diff types of covariances needed: 0 is same side, 1 is adjacent side, 2 is opposite side
#if the kernel is stationary, these can be translated to get the correct combinations
#K.bas[[0]] in reality is not needed bc when applying e-vectors you just get back the e-val
start=Sys.time()
K.bas=list(); 
for (i in 1:4) {
  for (j in 1:4) {
    name=paste(i,j,sep="");
    K.bas[[name]]=k(xg[[paste(i)]],xg[[paste(j)]]);
    if (i==j) {
      K.bas[[name]]=K.bas[[name]]+1e-10*diag(dim(K.bas[[name]])[1]);
    }
  }
}



L=length(basis);
phi=matrix(0,ncol=L,nrow=length(quad[[1]]));
for (i in 1:L) {
  phi[,i]=sapply(quad[[1]],basis[[i]]);
}
K.phi=K.bas[['11']]%*%diag(quad[[2]])%*%phi;
K.int=t(phi)%*%diag(quad[[2]])%*%K.phi;
E=eigen(K.int);

###contains eigenvectors of integral operator evaluated at quad[[1]]
V=t(E$vectors)%*%t(phi);



Cov=matrix(0,nrow=(4*length(basis)),ncol=4*length(basis));
sparse=0;
for (i in 1:4) {
  for (j in 1:4) {
    name=paste(i,j,sep="");
    if( abs(i-j)==2 & sparse==1) {
      print('independence');
    } else{
      Cov[((i-1)*length(basis)+1):(i*length(basis)),((j-1)*length(basis)+1):(j*length(basis))]=K.1=V%*%diag(quad[[2]])%*%K.bas[[name]]%*%diag(quad[[2]])%*%t(V);
    }
  }
}

E.covsp=eigen(Cov.sp,symmetric=T);
E.cov=eigen(Cov+1e-12*diag(4*length(basis)),symmetric=T);

t=seq(0,8,.1);
test=t(sapply(t,l)); n.test=dim(test)[1];

K.f=matrix(0,nrow=1,ncol=4*length(basis));

K.e=matrix(0,nrow=n.test,ncol=4*length(basis));
for (i in 1:4) {
  k.test=k(test,xg[[paste(i)]]);
  f.test=matrix(apply(xg[[paste(i)]],1,f),nrow=1);
  K.e[,((i-1)*length(basis)+1):(i*length(basis))]=k.test%*%diag(quad[[2]])%*%t(V);
  K.f[1,((i-1)*length(basis)+1):(i*length(basis))]=f.test%*%diag(quad[[2]])%*%t(V);
}
stop=Sys.time()

reg.bayes=seq(0,8,length.out=(4*(2*m+1)))

reg.vals=t(sapply(reg.bayes,l));
E.reg=eigen(k(reg.vals,reg.vals),symmetric=T);
k.reg.inv=E.reg$vectors%*%diag(1/abs(E.reg$values))%*%t(E.reg$vectors);
plot(t,abs(k(test,reg.vals)%*%k.reg.inv%*%as.matrix(f0(reg.bayes))-apply(test,1,f)),col='blue',type='l');

Cov.inv=E.cov$vectors%*%diag(1/abs(E.cov$values))%*%t(E.cov$vectors);

lines(t,abs(K.e%*%Cov.inv%*%t(K.f)-f0(t)),type='l');

K.post=function(X,Y) {
  n.x=dim(X)[1]; n.y=dim(Y)[1];
  K.x=matrix(0,nrow=n.x,ncol=4*length(basis));
  K.y=matrix(0,nrow=n.y,ncol=4*length(basis));
  for (i in 1:4) {
    k.x=k(X,xg[[paste(i)]]);
    k.y=k(Y,xg[[paste(i)]]);
    K.x[,((i-1)*length(basis)+1):(i*length(basis))]=k.x%*%diag(quad[[2]])%*%t(V);
    K.y[,((i-1)*length(basis)+1):(i*length(basis))]=k.y%*%diag(quad[[2]])%*%t(V);
  }
  return(k(X,Y)-K.x%*%Cov.inv%*%t(K.y));
}

mean.post=function(X,g) {
  n.x=dim(X);
  K.g=matrix(0,nrow=1,ncol=4*length(basis));
  K.e=matrix(0,nrow=n.x,ncol=4*length(basis));
  for (i in 1:4) {
    k.x=k(X,xg[[paste(i)]]);
    g.x=matrix(apply(xg[[paste(i)]],1,g),nrow=1);
    K.e[,((i-1)*length(basis)+1):(i*length(basis))]=k.x%*%diag(quad[[2]])%*%t(V);
    K.g[1,((i-1)*length(basis)+1):(i*length(basis))]=g.x%*%diag(quad[[2]])%*%t(V);
  }
  return(K.e%*%Cov.inv%*%t(K.g));
}

N=10;




#D=2*randomLHS(n=N,k=2)-1;
reg.bayes=seq(0,8,length.out=(4*length(basis)))
n.reg=length(reg.bayes);
reg.vals=t(sapply(reg.bayes,l));
D.reg=rbind(D,reg.vals);

x1=seq(-.99,.99,length.out=20); x2=x1;
test=cbind(rep(x1,length(x2)),rep(x2,1,each=length(x1)));
t=seq(0,8,.1);
test=t(sapply(t,l)); n.test=dim(test)[1]; test=rbind(.98*test);

y=apply(D,1,f)+rnorm(N,mean=0,sd=0);
y.reg=apply(D.reg,1,f)+rnorm((N+n.reg),mean=0,sd=0);

start=Sys.time();
mu.post=mean.post(test,f);

post.mean=mu.post+K.post(test,D)%*%solve(K.post(D,D)+1e-12*diag(dim(D)[1]))%*%(y-mean.post(D,f));
stop=Sys.time();
print(stop-start);

start=Sys.time();
E.reg=eigen(k(D.reg,D.reg))#+1e-10*diag(N+n.reg));
K.reg.inv=E.reg$vectors%*%diag(1/abs(E.reg$values))%*%t(E.reg$vectors);
post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg
stop=Sys.time();
stop-start;

f.test=apply(test,1,f)
#z=matrix(f.test,nrow=length(x1),ncol=length(x2),byrow=F);
#plot_ly(x=~x1,y=~x2,z=~z,type='surface')

plot(post.mean,f.test);
plot(post.mean.reg,f.test)

summary(lm(f.test~post.mean-1));
summary(lm(f.test~post.mean.reg-1));

cheb.max[m-4]=max(abs(f.test-post.mean));
reg.max[m-4]=max(abs(f.test-post.mean.reg));

cheb.rmse[m-4]=sqrt(sum((f.test-post.mean)^2)/dim(test)[1])
reg.rmse[m-4]=sqrt(sum((f.test-post.mean.reg)^2)/dim(test)[1])


cheb.base=lapply(chebyshev.c.polynomials(m,normalized=T),as.function);
var.change=function(f,scale,f.arg) {
  function(x) scale*f(f.arg(x));
}

cheb=lapply(cheb.base,function(x) var.change(x,2^(3/2),function(y) 2*y));
quad=gausschb(4*m); quad[[1]]=seq(-1,1,length.out=length(quad[[1]])); n.nodes=length(quad[[1]]);
trig=trig.poly.basis(-1,1,ceiling(m/2));
#quad=gausslg(-1,1,3*m); n.nodes=length(quad[[1]]);
#leg=lapply(legendre.polynomials(m,normalized=T),as.function); #list of legendre polynomials

xg=list(); xg[["1"]]=cbind(quad[[1]],-1); xg[["2"]]=cbind(1,quad[[1]]); #goes left to right and down to up
xg[["3"]]=cbind(quad[[1]],1); xg[["4"]]=cbind(-1,quad[[1]]);

basis=trig;

ip=function(f,g) {
  ###assume f and g are already in vector form at points quad[[1]]
  value=sum(quad[[2]]*f*g);
  return(value);
}

#diff types of covariances needed: 0 is same side, 1 is adjacent side, 2 is opposite side
#if the kernel is stationary, these can be translated to get the correct combinations
#K.bas[[0]] in reality is not needed bc when applying e-vectors you just get back the e-val
start=Sys.time()
K.bas=list(); 
for (i in 1:4) {
  for (j in 1:4) {
    name=paste(i,j,sep="");
    K.bas[[name]]=k(xg[[paste(i)]],xg[[paste(j)]]);
    if (i==j) {
      K.bas[[name]]=K.bas[[name]]+1e-10*diag(dim(K.bas[[name]])[1]);
    }
  }
}



L=length(basis);
phi=matrix(0,ncol=L,nrow=length(quad[[1]]));
for (i in 1:L) {
  phi[,i]=sapply(quad[[1]],basis[[i]]);
}
K.phi=K.bas[['11']]%*%diag(quad[[2]])%*%phi;
K.int=t(phi)%*%diag(quad[[2]])%*%K.phi;
E=eigen(K.int);

###contains eigenvectors of integral operator evaluated at quad[[1]]
V=t(E$vectors)%*%t(phi);



Cov=matrix(0,nrow=(4*length(basis)),ncol=4*length(basis));
sparse=0;
for (i in 1:4) {
  for (j in 1:4) {
    name=paste(i,j,sep="");
    if( abs(i-j)==2 & sparse==1) {
      print('independence');
    } else{
      Cov[((i-1)*length(basis)+1):(i*length(basis)),((j-1)*length(basis)+1):(j*length(basis))]=K.1=V%*%diag(quad[[2]])%*%K.bas[[name]]%*%diag(quad[[2]])%*%t(V);
    }
  }
}

E.covsp=eigen(Cov.sp,symmetric=T);
E.cov=eigen(Cov+1e-12*diag(4*length(basis)),symmetric=T);

t=seq(0,8,.1);
test=t(sapply(t,l)); n.test=dim(test)[1];

K.f=matrix(0,nrow=1,ncol=4*length(basis));

K.e=matrix(0,nrow=n.test,ncol=4*length(basis));
for (i in 1:4) {
  k.test=k(test,xg[[paste(i)]]);
  f.test=matrix(apply(xg[[paste(i)]],1,f),nrow=1);
  K.e[,((i-1)*length(basis)+1):(i*length(basis))]=k.test%*%diag(quad[[2]])%*%t(V);
  K.f[1,((i-1)*length(basis)+1):(i*length(basis))]=f.test%*%diag(quad[[2]])%*%t(V);
}
stop=Sys.time()

reg.bayes=seq(0,8,length.out=(4*(2*m+1)))

reg.vals=t(sapply(reg.bayes,l));
E.reg=eigen(k(reg.vals,reg.vals),symmetric=T);
k.reg.inv=E.reg$vectors%*%diag(1/abs(E.reg$values))%*%t(E.reg$vectors);
plot(t,abs(k(test,reg.vals)%*%k.reg.inv%*%as.matrix(f0(reg.bayes))-apply(test,1,f)),col='blue',type='l');

Cov.inv=E.cov$vectors%*%diag(1/abs(E.cov$values))%*%t(E.cov$vectors);

lines(t,abs(K.e%*%Cov.inv%*%t(K.f)-f0(t)),type='l');

K.post=function(X,Y) {
  n.x=dim(X)[1]; n.y=dim(Y)[1];
  K.x=matrix(0,nrow=n.x,ncol=4*length(basis));
  K.y=matrix(0,nrow=n.y,ncol=4*length(basis));
  for (i in 1:4) {
    k.x=k(X,xg[[paste(i)]]);
    k.y=k(Y,xg[[paste(i)]]);
    K.x[,((i-1)*length(basis)+1):(i*length(basis))]=k.x%*%diag(quad[[2]])%*%t(V);
    K.y[,((i-1)*length(basis)+1):(i*length(basis))]=k.y%*%diag(quad[[2]])%*%t(V);
  }
  return(k(X,Y)-K.x%*%Cov.inv%*%t(K.y));
}

mean.post=function(X,g) {
  n.x=dim(X);
  K.g=matrix(0,nrow=1,ncol=4*length(basis));
  K.e=matrix(0,nrow=n.x,ncol=4*length(basis));
  for (i in 1:4) {
    k.x=k(X,xg[[paste(i)]]);
    g.x=matrix(apply(xg[[paste(i)]],1,g),nrow=1);
    K.e[,((i-1)*length(basis)+1):(i*length(basis))]=k.x%*%diag(quad[[2]])%*%t(V);
    K.g[1,((i-1)*length(basis)+1):(i*length(basis))]=g.x%*%diag(quad[[2]])%*%t(V);
  }
  return(K.e%*%Cov.inv%*%t(K.g));
}

N=10;




#D=2*randomLHS(n=N,k=2)-1;
reg.bayes=seq(0,8,length.out=(4*length(basis)))
n.reg=length(reg.bayes);
reg.vals=t(sapply(reg.bayes,l));
D.reg=rbind(D,reg.vals);

x1=seq(-.99,.99,length.out=20); x2=x1;
test=cbind(rep(x1,length(x2)),rep(x2,1,each=length(x1)));
t=seq(0,8,.1);
test=t(sapply(t,l)); n.test=dim(test)[1]; test=rbind(.98*test);

y=apply(D,1,f)+rnorm(N,mean=0,sd=0);
y.reg=apply(D.reg,1,f)+rnorm((N+n.reg),mean=0,sd=0);

start=Sys.time();
mu.post=mean.post(test,f);

post.mean=mu.post+K.post(test,D)%*%solve(K.post(D,D)+1e-12*diag(dim(D)[1]))%*%(y-mean.post(D,f));
stop=Sys.time();
print(stop-start);

start=Sys.time();
E.reg=eigen(k(D.reg,D.reg))#+1e-10*diag(N+n.reg));
K.reg.inv=E.reg$vectors%*%diag(1/abs(E.reg$values))%*%t(E.reg$vectors);
post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg
stop=Sys.time();
stop-start;

f.test=apply(test,1,f)
#z=matrix(f.test,nrow=length(x1),ncol=length(x2),byrow=F);
#plot_ly(x=~x1,y=~x2,z=~z,type='surface')

plot(post.mean,f.test);
plot(post.mean.reg,f.test)

summary(lm(f.test~post.mean-1));
summary(lm(f.test~post.mean.reg-1));

trig.max[m-4]=max(abs(f.test-post.mean));
reg.max[m-4]=max(abs(f.test-post.mean.reg));

trig.rmse[m-4]=sqrt(sum((f.test-post.mean)^2)/dim(test)[1])
reg.rmse[m-4]=sqrt(sum((f.test-post.mean.reg)^2)/dim(test)[1])
}

par(mfrow=c(1,2));
D.reg=D; y.reg=y;
E.reg=eigen(k(D.reg,D.reg))#+1e-10*diag(N+n.reg));
K.reg.inv=E.reg$vectors%*%diag(1/abs(E.reg$values))%*%t(E.reg$vectors);
post.mean.reg=k(test,D.reg)%*%K.reg.inv%*%y.reg;
basic.max=max(abs(f.test-post.mean.reg));
basic.rmse=sqrt(sum((f.test-post.mean.reg)^2)/dim(test)[1])

#plot(cheb.rmse,type='l',ylim=c(0,basic.rmse)*1.2);
#lines(reg.rmse,type='l',col='blue')
#lines(leg.rmse,type='l',col='green')
#lines(-5:35,rep(basic.rmse,length(-5:35)),lwd=2);

par(mfrow=c(1,1), mar=c(4, 4, 1, 1))
plot(5:J, cheb.max,type='l',ylim=c(0,basic.max)*1,ylab="Maximum Error",xlab='m');
lines(5:J, reg.max,type='l',col='blue')
lines(5:J, leg.max,type='l',col='green')
lines(5:J, trig.max,type='l',col='red')
lines(-5:35,rep(basic.max,length(-5:35)),lwd=2);
legend('right', legend=c('Chebyshev', 'Kernel', 'Legendre', 'Trig', 'None'), col=c('black', 'blue', 'green', 'red', 'black'), 
       lty=1, cex=1, lwd=c(rep(1, 4), 2))



