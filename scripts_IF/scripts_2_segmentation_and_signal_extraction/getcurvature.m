%%
function ng=getcurvature(n,gp)

%%
ns=[n(1+gp:-1:1+1),n,n(end-1:-1:end-gp)];
theta1=atan2(ns(1+gp:end-gp)-ns(1:end-2*gp),gp);
theta2=atan2(ns(1+2*gp:end)-ns(1+gp:end-gp),gp);
ng_a=acos(cos(theta2-theta1));
ng_v=ns(1:end-2*gp)+ns(1+2*gp:end)-2*ns(1+gp:end-gp);
ng=ng_a.*(ng_v>0)-ng_a.*(ng_v<0);