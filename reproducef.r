library(PolynomF)

#x2 = seq(-1, 1, length.out=6);
#a = runif(4, -1, 1);

f1 = function(x) {
  result = rep(0, length(x));
  for (i in 1:4) {
    result = result + a[i] * exp(-abs(x - x2[i]) ^ 1.9);
  }
  return(result);
}

poly = poly_calc(x2, f1(x2))
f2 = function(x) {return(predict(poly,x))};

test = matrix(seq(-1 , 1, .001), ncol=1);
supf1 = max(abs(f1(test)));
supf2 = max(abs(f2(test)));

f = f2;
supf = supf2;

J = seq(6, 51, 5);
max.new=rep(0, length(J));
max.bayes = rep(0, length(J));

for (j in 1:length(J)) {


m = J[j];

cheb.base=lapply(chebyshev.c.polynomials(m,normalized=T),as.function);
var.change=function(f,scale,f.arg) {
  function(x) scale*f(f.arg(x));
}


cheb=lapply(cheb.base,function(x) var.change(x,2^(3/2),function(y) 2*y));
quad=gausschb(m); quad[[1]]=seq(-1,1,length.out=length(quad[[2]])); n.nodes=length(quad[[1]]);
quad=gausslg(-1,1, 2 * m); n.nodes=length(quad[[1]]);
#trig=trig.poly.basis(-1,1,ceiling(m/2));
leg=lapply(legendre.polynomials(m,normalized=T),as.function); #list of legendre polynomials


xg=matrix(quad[[1]], ncol=1);


basis=leg


k=function(X,Y) {
  n=dim(X)[1]; m=dim(Y)[1];
  X = cbind(X, 1); Y = cbind(Y, 1);
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


ip=function(f,g) {
  ###assume f and g are already in vector form at points quad[[1]]
  value=sum(quad[[2]]*f*g);
  return(value);
}


K.bas = k(xg, xg);


L=length(basis);
phi=matrix(0,ncol=L,nrow=length(quad[[1]]));
for (i in 1:L) {
  phi[,i]=sapply(quad[[1]],basis[[i]]);
}
K.phi=K.bas%*%diag(quad[[2]])%*%phi;
K.int=t(phi)%*%diag(quad[[2]])%*%K.phi;
E=eigen(K.int);

###contains eigenvectors of integral operator evaluated at quad[[1]]
V=t(E$vectors)%*%t(phi);

Cov = K.1 = V %*% diag(quad[[2]]) %*% K.bas %*% diag(quad[[2]]) %*% t(V);


E.cov = eigen(Cov + 0 * diag(L), symmetric=T);

k.test=k(test,xg);
f.test=matrix(apply(xg,1,f),nrow=1);

K.e = k.test%*%diag(quad[[2]])%*%t(V);
K.f = f.test%*%diag(quad[[2]])%*%t(V);

stop=Sys.time()

reg.bayes=seq(-1, 1, length.out = (m));

reg.vals=matrix(reg.bayes, ncol=1);
E.reg=eigen(k(reg.vals,reg.vals),symmetric=T);
k.reg.inv=E.reg$vectors%*%diag(1/abs(E.reg$values))%*%t(E.reg$vectors);
#plot(test,abs(k(test,reg.vals)%*%k.reg.inv%*%as.matrix(f(reg.bayes))-apply(test,1,f)),col='blue',type='l');

Cov.inv=E.cov$vectors%*%diag(1/abs(E.cov$values))%*%t(E.cov$vectors);

error_new = abs(K.e%*%Cov.inv%*%t(K.f)-f(test));
error_bayes = abs(k(test,reg.vals)%*%k.reg.inv%*%as.matrix(f(reg.bayes))-apply(test,1,f))
#plot(test,error_new, ylim = c(0, max(c(error_new,error_bayes)) * 1.1), type='b');
#lines(test,error_bayes,col='blue',type='b');


max.new[j] = max(error_new);
max.bayes[j] = max(error_bayes);



}


yrange = c(min(log(max.bayes / supf), log(max.new / supf)) * 1.1, max(log(max.bayes / supf), log(max.new / supf)) * .9);
ylabel = "Log Maximum Error";
xlabel = "Log Number of Basis Functions";

plot(log(J), log(max.bayes/supf), ylim=yrange, xlab=xlabel, ylab=ylabel, type='b', col='blue', main=expression(paste(f[2])));
lines(log(J), log(max.new/supf), type='b', col = 'black')

legend('bottomleft', legend=c("Interpolation Method", "RR Method"), col=c('blue', 'black'), lty=1, cex=.7)






