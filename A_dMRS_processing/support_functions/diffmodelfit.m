function res=diffmodelfit(y,v,b,m,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Saad Jbabdi 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% results=diffmodelfit(data,bvecs,bvals,model,[numfib=2])
% data is NxT 
% bvecs is 3xT
% bvals is 1xT
% model one of 'DTI', 'DKI', 'PVM1', 'PVM2', 'PVM10' 'PVM20'
% numfib is ignored if model=DTI


if(strcmpi(m,'DTI'))
    res=fit_tensor(y,v,b);
elseif(strcmpi(m,'DKI'))
    res=fit_tensor_kurt(y,v,b);
elseif(strcmpi(m,'PVM1'))
    if(nargin<5);n=1;end
    res=fit_pvm1(y,v,b,n);
elseif(strcmpi(m,'PVM2'))
    if(nargin<5);n=1;end
    res=fit_pvm2(y,v,b,n);
elseif(strcmpi(m,'PVM10'))
    if(nargin<5);n=1;end
    res=fit_pvm10(y,v,b,n);
elseif(strcmpi(m,'PVM20'))
    if(nargin<5);n=1;end
    res=fit_pvm20(y,v,b,n);
end


function res=fit_tensor(y,v,b)
v=v';b=b';
A = [v(:,1).*v(:,1).*b ...
    2*v(:,1).*v(:,2).*b ...
    2*v(:,1).*v(:,3).*b ...
    v(:,2).*v(:,2).*b ...
    2*v(:,2).*v(:,3).*b ...
    v(:,3).*v(:,3).*b ...
    ones(length(b),1)];
r = -pinv(A)*log(y'); % r is 7xN
res.r=r;
res.pred = exp(-A*r)';
res.S0=exp(-r(7,:))';
res.MD=zeros(size(y,1),1);
res.L1=zeros(size(y,1),1);
res.L2=zeros(size(y,1),1);
res.L3=zeros(size(y,1),1);
res.FA=zeros(size(y,1),1);
res.V1=zeros(size(y,1),3);
res.V2=zeros(size(y,1),3);
res.V3=zeros(size(y,1),3);
for i=1:size(y,1)
tens = [r(1,i) r(2,i) r(3,i); ...
        r(2,i) r(4,i) r(5,i); ...
        r(3,i) r(5,i) r(6,i)];
[Vd,Dd] = eig(tens);
Dd = sort(diag(Dd),'descend');
res.MD(i) = sum(Dd)/3;
res.L1(i) = Dd(1);
res.L2(i) = Dd(2);
res.L3(i) = Dd(3);
res.RD(i) = .5*(res.L2(i)+res.L3(i));
res.FA(i) = sqrt(1.5*((res.L1(i)-res.MD(i))^2+...
                      (res.L2(i)-res.MD(i))^2+...
                      (res.L3(i)-res.MD(i))^2)/...
                      (res.L1(i)^2+res.L2(i)^2+res.L3(i)^2));
res.V1(i,:) = Vd(:,3)';
res.V2(i,:) = Vd(:,2)';
res.V3(i,:) = Vd(:,1)';
end
function res=fit_tensor_kurt(y,v,b)
v=v';b=b';
M = [v(:,1).*v(:,1).*b ...
    2*v(:,1).*v(:,2).*b ...
    2*v(:,1).*v(:,3).*b ...
    v(:,2).*v(:,2).*b ...
    2*v(:,2).*v(:,3).*b ...
    v(:,3).*v(:,3).*b ...
    ones(length(b),1)];
x=-b.^2/6;
x=x-M*pinv(M)*x;

A = [v(:,1).*v(:,1).*b ...
    2*v(:,1).*v(:,2).*b ...
    2*v(:,1).*v(:,3).*b ...
    v(:,2).*v(:,2).*b ...
    2*v(:,2).*v(:,3).*b ...
    v(:,3).*v(:,3).*b ...
    x ...
    ones(length(b),1)];
r = -pinv(A)*log(y'); % r is 8xN
res.r=r;
res.pred = exp(-A*r)';
res.K=r(7,:);
res.S0=exp(-r(8,:))';
res.MD=zeros(size(y,1),1);
res.L1=zeros(size(y,1),1);
res.L2=zeros(size(y,1),1);
res.L3=zeros(size(y,1),1);
res.FA=zeros(size(y,1),1);
res.V1=zeros(size(y,1),3);
res.V2=zeros(size(y,1),3);
res.V3=zeros(size(y,1),3);
for i=1:size(y,1)
tens = [r(1,i) r(2,i) r(3,i); ...
        r(2,i) r(4,i) r(5,i); ...
        r(3,i) r(5,i) r(6,i)];
[Vd,Dd] = eig(tens);
Dd = sort(diag(Dd),'descend');
res.MD(i) = sum(Dd)/3;
res.L1(i) = Dd(1);
res.L2(i) = Dd(2);
res.L3(i) = Dd(3);
res.FA(i) = sqrt(1.5*((res.L1(i)-res.MD(i))^2+...
                      (res.L2(i)-res.MD(i))^2+...
                      (res.L3(i)-res.MD(i))^2)/...
                      (res.L1(i)^2+res.L2(i)^2+res.L3(i)^2));
res.V1(i,:) = Vd(:,3)';
res.V2(i,:) = Vd(:,2)';
res.V3(i,:) = Vd(:,1)';
end

function res=fit_pvm1(y,v,b,n)
dti=fit_tensor(y,v,b);
res.S0=zeros(size(y,1),1);
res.d=zeros(size(y,1),1);
res.th=zeros(size(y,1),n);
res.ph=zeros(size(y,1),n);
res.f=zeros(size(y,1),n);
res.pred=zeros(size(y));
for i=1:size(y,1)
    x0=[dti.S0(i);2*dti.MD(i);zeros(2*n,1);dti.FA(i)*ones(n,1)];
    c=[v(1,:)';v(2,:)';v(3,:)';b(:);y(i,:)'];
    lb=[0;0;-inf*ones(n,1);-inf*ones(n,1);zeros(n,1)];
    ub=[inf;inf;inf*ones(n,1);inf*ones(n,1);ones(n,1)];
    x=lsqnonlin(@(x) nlfun1(x,c),x0,lb,ub);
    res.S0(i)=x(1);
    res.d(i)=x(2);    
    res.th(i,:)=x(3:3+n-1);
    res.ph(i,:)=x(3+n:3+2*n-1);
    res.f(i,:)=x(3+2*n:end);
    [~,j]=sort(res.f(i,:),'descend');
    res.th(i,:)=res.th(i,j);
    res.ph(i,:)=res.ph(i,j);
    res.f(i,:)=res.f(i,j);
    res.pred(i,:)=pvm_forwardmodel(res.S0(i),res.d(i),res.f(i,:)',res.th(i,:)',res.ph(i,:)',b,v,0,0);
end
function res=fit_pvm10(y,v,b,n)
dti=fit_tensor(y,v,b);
res.S0=zeros(size(y,1),1);
res.d=zeros(size(y,1),1);
res.f0=zeros(size(y,1),1);
res.th=zeros(size(y,1),n);
res.ph=zeros(size(y,1),n);
res.f=zeros(size(y,1),n);
res.pred=zeros(size(y));
for i=1:size(y,1)
    x0=[dti.S0(i);2*dti.MD(i);0.01;zeros(2*n,1);dti.FA(i)*ones(n,1)];
    c=[v(1,:)';v(2,:)';v(3,:)';b(:);y(i,:)'];
    lb=[0;0;0;-inf*ones(n,1);-inf*ones(n,1);zeros(n,1)];
    ub=[inf;inf;1;inf*ones(n,1);inf*ones(n,1);ones(n,1)];
    x=lsqnonlin(@(x) nlfun10(x,c),x0,lb,ub);
    res.S0(i)=x(1);
    res.d(i)=x(2);    
    res.f0(i)=x(3);
    res.th(i,:)=x(4:4+n-1);
    res.ph(i,:)=x(4+n:4+2*n-1);
    res.f(i,:)=x(4+2*n:end);
    [~,j]=sort(res.f(i,:),'descend');
    res.th(i,:)=res.th(i,j);
    res.ph(i,:)=res.ph(i,j);
    res.f(i,:)=res.f(i,j);
    res.pred(i,:)=pvm_forwardmodel(res.S0(i),res.d(i),res.f(i,:)',res.th(i,:)',res.ph(i,:)',b,v,res.f0(i),0);
end

function res=fit_pvm2(y,v,b,n)
dti=fit_tensor(y,v,b);
res.S0=zeros(size(y,1),1);
res.d=zeros(size(y,1),1);
res.dstd=zeros(size(y,1),1);
res.th=zeros(size(y,1),n);
res.ph=zeros(size(y,1),n);
res.f=zeros(size(y,1),n);
res.pred=zeros(size(y));
for i=1:size(y,1)
    x0=[dti.S0(i);2*dti.MD(i);2*dti.MD(i)/5;zeros(2*n,1);dti.FA(i)*ones(n,1)];
    c=[v(1,:)';v(2,:)';v(3,:)';b(:);y(i,:)'];
    lb=[0;0;0;-inf*ones(n,1);-inf*ones(n,1);zeros(n,1)];
    ub=[inf;inf;inf;inf*ones(n,1);inf*ones(n,1);ones(n,1)];
    x=lsqnonlin(@(x) nlfun2(x,c),x0,lb,ub);
    res.S0(i)=x(1);
    res.d(i)=x(2);
    res.dstd(i)=x(3);
    res.th(i,:)=x(4:4+n-1);
    res.ph(i,:)=x(4+n:4+2*n-1);
    res.f(i,:)=x(4+2*n:end);
    [~,j]=sort(res.f(i,:),'descend');
    res.th(i,:)=res.th(i,j);
    res.ph(i,:)=res.ph(i,j);
    res.f(i,:)=res.f(i,j);
    res.pred(i,:)=pvm_forwardmodel(res.S0(i),res.d(i),res.f(i,:)',res.th(i,:)',res.ph(i,:)',b,v,0,res.dstd(i));
end
function res=fit_pvm20(y,v,b,n)
dti=fit_tensor(y,v,b);
res.S0=zeros(size(y,1),1);
res.d=zeros(size(y,1),1);
res.dstd=zeros(size(y,1),1);
res.f0=zeros(size(y,1),1);
res.th=zeros(size(y,1),n);
res.ph=zeros(size(y,1),n);
res.f=zeros(size(y,1),n);
res.pred=zeros(size(y));
for i=1:size(y,1)
    x0=[dti.S0(i);2*dti.MD(i);2*dti.MD(i)/5;0.01;zeros(2*n,1);dti.FA(i)*ones(n,1)];
    c=[v(1,:)';v(2,:)';v(3,:)';b(:);y(i,:)'];
    lb=[0;0;0;0;-inf*ones(n,1);-inf*ones(n,1);zeros(n,1)];
    ub=[inf;inf;inf;1;inf*ones(n,1);inf*ones(n,1);ones(n,1)];
    x=lsqnonlin(@(x) nlfun20(x,c),x0,lb,ub);
    res.S0(i)=x(1);
    res.d(i)=x(2);
    res.dstd(i)=x(3);
    res.f0(i)=x(4);
    res.th(i,:)=x(5:5+n-1);
    res.ph(i,:)=x(5+n:5+2*n-1);
    res.f(i,:)=x(5+2*n:end);
    [~,j]=sort(res.f(i,:),'descend');
    res.th(i,:)=res.th(i,j);
    res.ph(i,:)=res.ph(i,j);
    res.f(i,:)=res.f(i,j);
    res.pred(i,:)=pvm_forwardmodel(res.S0(i),res.d(i),res.f(i,:)',res.th(i,:)',res.ph(i,:)',b,v,res.f0(i),res.dstd(i));
end
function e=nlfun1(x,c)
nt=size(c,1)/5;
v=reshape(c(1:3*nt),nt,3)';
b=c(3*nt+1:4*nt)';
y=c(4*nt+1:end);
nf=(length(x)-2)/3;
s0=x(1);d=x(2);
th=x(3:3+nf-1);
ph=x(3+nf:3+2*nf-1);
f=x(3+2*nf:end);
s=pvm_forwardmodel(s0,d,f,th,ph,b,v,0,0);
e=s'-y;
function e=nlfun10(x,c)
nt=size(c,1)/5;
v=reshape(c(1:3*nt),nt,3)';
b=c(3*nt+1:4*nt)';
y=c(4*nt+1:end);
nf=(length(x)-3)/3;
s0=x(1);d=x(2);f0=x(3);
th=x(4:4+nf-1);
ph=x(4+nf:4+2*nf-1);
f=x(4+2*nf:end);
s=pvm_forwardmodel(s0,d,f,th,ph,b,v,f0,0);
e=s'-y;
function e=nlfun2(x,c)
nt=size(c,1)/5;
v=reshape(c(1:3*nt),nt,3)';
b=c(3*nt+1:4*nt)';
y=c(4*nt+1:end);
nf=(length(x)-3)/3;
s0=x(1);d=x(2);dstd=x(3);
th=x(4:4+nf-1);
ph=x(4+nf:4+2*nf-1);
f=x(4+2*nf:end);
s=pvm_forwardmodel(s0,d,f,th,ph,b,v,0,dstd);
e=s'-y;
function e=nlfun20(x,c)
nt=size(c,1)/5;
v=reshape(c(1:3*nt),nt,3)';
b=c(3*nt+1:4*nt)';
y=c(4*nt+1:end);
nf=(length(x)-4)/3;
s0=x(1);d=x(2);dstd=x(3);f0=x(4);
th=x(5:5+nf-1);
ph=x(5+nf:5+2*nf-1);
f=x(5+2*nf:end);
s=pvm_forwardmodel(s0,d,f,th,ph,b,v,f0,dstd);
e=s'-y;
function s=pvm_forwardmodel(s0,d,f,th,ph,b,v,f0,dstd)
% if(nargin==3);s=pvm_forwardmodel(s0.s0,s0.d,s0.f,s0.th,s0.ph,d,f,s0.f0,s0.dstd);end
[beta,alpha]=cart2sph(v(1,:)',v(2,:)',v(3,:)');
alpha=pi/2-alpha;
angtmp = zeros(length(f),length(alpha));
for i=1:length(f)
 angtmp(i,:)=cos(ph(i)-beta').*sin(alpha').*sin(th(i)) + cos(alpha').*cos(th(i));
end
angtmp = abs(angtmp.^2);
bg=repmat(b,length(f),1).*angtmp;
fs=repmat(f,1,length(b));
if(dstd<1e-5)
    s=s0*(f0+sum( fs.*exp(-bg*d) ,1)+(1-sum(f)-f0)*exp(-b*d));
else
    alph=d^2/dstd^2;
    beta=d/dstd^2;
    s=s0*(f0+sum(fs.*exp(alph*log(beta./(beta+bg))),1) +(1-sum(f)-f0)*exp(alph*log(beta./(beta+b))));
end




