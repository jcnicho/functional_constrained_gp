####link function parameterization from [0,8] to the boundary of [-1,1]^2

# l=function(t) {
#   if (t <2 & t>=0) {
#     return(c(t-1,-1));
#   } else if(t <4 & t>=2) {
#     return(c(1,-3+t));
#   } else if (t<6 & t>=4) {
#     return(c(5-t,1));
#   } else {
#     return(c(-1,7-t));
#   }
# }
# 
# 
# l.inv=function(x) {
#   if (x[2]==-1) {
#     return(x[1]+1);
#   } else if (x[1]==1) {
#     return(3+x[2]);
#   } else if (x[2]==1) {
#     return(5-x[1]);
#   } else {
#     return(7-x[2]);
#   }
# }


l=function(t) {
  return(c(t,t));
}

l.inv=function(x) {
  return(x[1]);
}






