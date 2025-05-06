function LMAPSPressureVelocityNSGlobal % Nevenly distributed points (efficient version)
clear all;clc;
% Local MPS method; modified Helmholtz Eq.
% ni-----# of interior points
% nb-----# of boundary points
% c------shape paramete
% scale with maximum radius
%generate evenly distribute grid points on a unit square
%tic
X0=0; X1=1; Y0=0; Y1=1;  nn=50;  dx=1/nn; c=3;%20*dx; 
dt=2.0e-3;%0.0001;% 
Re=1000;
t0=clock;
dy=dx;
n=(X1*nn+1)*(Y1*nn+1) ; nb=2*X1*nn+2*Y1*nn; ni=n-nb;
[X,Y]=meshgrid(X0:dx:X1,Y0:dy:Y1);
x1=reshape(X,[],1);  y1=reshape(Y,[],1);
intnode=zeros(ni,2); bdpt=zeros(nb,2);
j=0;k=0;
%generate the interior and boundary nodes evenly
for i=1:n
    if(x1(i)==X0 || x1(i)==X1 || y1(i)==Y0 || y1(i)==Y1)
       j=j+1;
       bdpt(j,1)=x1(i); bdpt(j,2)=y1(i);
    else 
       k=k+1;
       intnode(k,1)=x1(i); intnode(k,2)=y1(i);
    end
end
x=[intnode(:,1); bdpt(:,1)]; y=[intnode(:,2); bdpt(:,2)];


uk=zeros(n,1);
vk=zeros(n,1);
%pk=zeros(n,1);
sort=zeros(n,1);
for i=1:ni   
  sort(i)=1;  
end 
for i=ni+1:n
  if(y(i)==Y1)
    sort(i)=0;
  else
    sort(i)=2;
  end
end
   ioid=find(sort(:)==0);
   bouid=find(sort(:)==2);
% Normal vector of boundary points
bnorm=zeros(n,2);
for i=1:n
 if(sort(i)==2 || sort(i)==0)   
   id2=find(((x(:)-x(i)).^2+(y(:)-y(i)).^2)<2*dx^2+1e-6);
   x1=mean(x(id2))-x(i); y1=mean(y(id2))-y(i);
   bnorm(i,:)=-[x1 y1]/norm([x1 y1]);
 end
end
 bnorm=fix(bnorm*1e5)/1e5; 
% quiver(x(:),y(:),bnorm(:,1),bnorm(:,2))
% plot(x(bouid(:)),y(bouid(:)),'.r');
%produce the boundary condition of velocity
uk(bouid)=0;
vk(bouid)=0;
uk(ioid)=1;
vk(ioid)=0; 
%interpolation for 2-order, x-direction and y-direcion  derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A_p, A_v, A_x, A_y, A_f]=globalmatrix(x,y,n,c,sort,bnorm,dx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% spy(A_f)
% pause
for k=1:15000
du=A_v*uk;dv=A_v*vk;
%produce the boundary condition of velocity    
du(bouid)=0; %??????????????????
dv(bouid)=0; %??????????????????
%du(ioid)=0; 
%dv(ioid)=0; %??????????????????
duxy=uk.*(A_x*uk)+vk.*(A_y*uk);
dvxy=uk.*(A_x*vk)+vk.*(A_y*vk); 
duxy(bouid)=0; 
dvxy(bouid)=0;


um=uk+dt*((1/Re)*du-duxy);
vm=vk+dt*((1/Re)*dv-dvxy);
%um=uk+dt*((1/Re)*du'-(uk'.*(A_x*uk'))'-(vk'.*(A_y*uk'))');
%vm=vk+dt*((1/Re)*dv'-(uk'.*(A_x*vk'))'-(vk'.*(A_y*vk'))');
%produce the boundary condition of velocity    




% %first approach handling the forcing term: The new foring term involving
% %the bd conditions
% dum=A_x*um+A_y*vm; %for LMPAS
% dum(bouid)=((um(bouid)-uk(bouid)).*bnorm(bouid,1) ...
%            + (vm(bouid)-vk(bouid)).*bnorm(bouid,2));
% dum(ioid)=0;  %????????????????????????????
% dum=dum+A_f*dum;%for LHMAPS
% 


% %second approach handling the forcing term:The new foring term not involving
% %the bd conditions
% dum=A_x*um+A_y*vm; %for LMPAS
% dum(ioid)=0;  %????????????????????????????
% dum(bouid)=0;
% dum=dum+A_f*dum;%for LHMAPS
% dum(bouid)=((um(bouid)-uk(bouid)).*bnorm(bouid,1) ...
%            + (vm(bouid)-vk(bouid)).*bnorm(bouid,2));
% %The results shows that appraoch 1 and 2 are the same when c=3.


%third approach handling the forcing term: The new foring term not involving
%the bd conditions. It is actually the same as the second approach
dum=A_x*um+A_y*vm; %for LMPAS
dum(ioid)=0;  %????????????????????????????
dum(bouid)=0;
G=0*dum;G(bouid)=((um(bouid)-uk(bouid)).*bnorm(bouid,1) ...
           + (vm(bouid)-vk(bouid)).*bnorm(bouid,2));
dum=dum+A_f*dum+G;%for LHMAPS




pk = (1 / dt) * (A_p \ dum);  


% produce the values of velocity
uk=um-dt*(A_x*pk);
vk=vm-dt*(A_y*pk);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%produce the boundary condition of velocity
uk(bouid)=0;
vk(bouid)=0;
uk(ioid)=1;
vk(ioid)=0;


if(mod(k,10)==1)
   cla 
[sx,sy]=meshgrid(X0:dx:X1,Y0:dy:Y1);
svx=griddata(x(:),y(:),uk(:),sx,sy);
svy=griddata(x(:),y(:),vk(:),sx,sy); 
sp=griddata(x(:),y(:),pk,sx,sy);
svx(isnan(svx))=0;
svy(isnan(svy))=0;
h=streamslice(sx,sy,svx,svy,5);
set(h,'color','k') 
set(gca,'linewidth',10);
%hold on
%plot(x(sort==2),y(sort==2),'.','color','r');
axis([0 X1 0 Y1]);
xlabel('x');ylabel('y');zlabel('Solution');
set(get(gca,'XLabel'),'FontSize',30);
set(get(gca,'YLabel'),'FontSize',30);
set(get(gca,'ZLabel'),'FontSize',30);
set(gca,'linewidth',1);
set(gca,'FontSize',30);
%colormap(gray);
colormap([0 0 0]);
title(['t=' num2str(k)])
pause(0.02) 
end 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%
cpu=etime(clock,t0);
fprintf('cputime= %6.2f\n',cpu);
[sx,sy]=meshgrid(X0:dx:X1,Y0:dy:Y1);
svx=griddata(x(:),y(:),uk(:),sx,sy);
svy=griddata(x(:),y(:),vk(:),sx,sy); 
sp=griddata(x(:),y(:),pk,sx,sy);
svx(isnan(svx))=0;
svy(isnan(svy))=0;
%SPY(GXY)
figure
h=streamslice(sx,sy,svx,svy,2);
set(h,'color','k')
hold on
plot(x(sort==2),y(sort==2),'.','color','r');
axis([0 X1 0 Y1]);
figure
scatter3(x(:),y(:),pk(:),'.','k');
figure
quiver(x,y,uk,vk);
%mesh(sx, sy, sp);
%surf(sx, sy, sp);
%contour(sx, sy, sp);
%quiver(x,y,uk,vk);
figure
for i=1:n
    if(x(i)==0.5)
plot(y(i),uk(i),'k pentagram','MarkerSize',12);
    end
hold on
end
legend('LMAPS,Re=100','Location','North')
xlabel('Y');ylabel('U');
set(get(gca,'XLabel'),'FontSize',13);
set(get(gca,'YLabel'),'FontSize',13);
set(gca,'linewidth',1);
set(gca,'FontSize',13);


figure
for i=1:n
    if(y(i)==1.0)
plot(x(i),vk(i),'k pentagram','MarkerSize',12);
    end
hold on
end
legend('LMAPS,Re=100','Location','North')
xlabel('X');ylabel('V');
set(get(gca,'XLabel'),'FontSize',13);
set(get(gca,'YLabel'),'FontSize',13);
set(gca,'linewidth',1);
set(gca,'FontSize',13);
return


%%Produce the global LMAPS sparse matrix
function [A_p, A_v, A_x, A_y, A_f]=globalmatrix(x,y,n,c,sort,bnorm,dx)
%z=[x', y'];
%tree=kdtree_build(z); %search
A_p=sparse(n,n);
A_f=sparse(n,n);%LHMPAS
A_v=sparse(n,n);
A_x=sparse(n,n);
A_y=sparse(n,n);


for i=1:n
%[id, D]=kdtree_k_nearest_neighbors(tree,z(i,:),ns); 
 if(sort(i)==1)
  R0=dx;
  dis2=(x(:)-x(i)).^2+(y(:)-y(i)).^2;
  id=find(dis2<=2*R0.^2+1e-6);
  dis2=dis2(id);
  ns=length(id);
  mqr=sqrt(dis2*c^2+1);%sqrt(dis2+c^2);%




  f_p=mqr;
  f_v=mqr;
  f_x=1/3*(1./mqr+mqr-1./(mqr+mqr.^2)).*(x(i)-x(id));%1/3*(c^2./mqr+mqr-c^3./(c*mqr+mqr.^2)).*(x(i)-x(id));%
  f_y=1/3*(1./mqr+mqr-1./(mqr+mqr.^2)).*(y(i)-y(id));%1/3*(c^2./mqr+mqr-c^3./(c*mqr+mqr.^2)).*(y(i)-y(id));%




  idH=find(dis2>0 & dis2<=R0.^2+1e-6);
  dis2H=dis2(idH);
  m=length(idH);
  mqrH=2*c^2./sqrt(dis2H*c^2+1)-dis2H*c^4./sqrt(dis2H*c^2+1).^3;
  f_pH=[mqr; mqrH];


elseif(sort(i)==0)
  R0=dx;
  dis2=(x(:)-x(i)).^2+(y(:)-y(i)).^2;
  id=find(dis2<=2*R0.^2+1e-6);
  dis2=dis2(id);
  ns=length(id);
  mqr=sqrt(dis2*c^2+1);%sqrt(dis2+c^2);%
  f_p=mqr.*(dis2/9+4/(9*c^2))-log(mqr+1)/(3*c^2);%(4*c^2+dis2).*mqr/9-c^3*log(c+mqr)/3;%
  f_v=f_p;
  f_x=1/3*(1./mqr+mqr-1./(mqr+mqr.^2)).*(x(i)-x(id));%1/3*(c^2./mqr+mqr-c^3./(c*mqr+mqr.^2)).*(x(i)-x(id));%
  f_y=1/3*(1./mqr+mqr-1./(mqr+mqr.^2)).*(y(i)-y(id));%1/3*(c^2./mqr+mqr-c^3./(c*mqr+mqr.^2)).*(y(i)-y(id));%


    % idH=find(dis2<=R0.^2+1e-6);
  % dis2H=dis2(idH);
  % m=length(idH);
  % mqrH=2*c^2./sqrt(1+(dis2H*c).^2)-dis2H*c^4./sqrt(1+(dis2H*c).^2).^3;
  % f_pH=[mqr; mqrH];


 else  %(sort(i)==2)
  R0=3*dx;
  dis2=(x(:)-x(i)).^2+(y(:)-y(i)).^2;
  id=find(dis2<=R0.^2+1e-6);
  dis2=dis2(id);
  ns=length(id);
  mqr=sqrt(dis2*c^2+1);%sqrt(dis2+c^2);%
     
   if(bnorm(i,1)==0)
     bdx=0;
   else
%    bdx=1/bnorm(i,1);  
     bdx=bnorm(i,1);
   end
   if(bnorm(i,2)==0)
     bdy=0;
   else
 %    bdy=1/bnorm(i,2);
     bdy=bnorm(i,2);
   end 
   % why not just use bdx=bnorm(i,1); and bdy=bnorm(i,2);


   f_x=1/3*(1./mqr+mqr-1./(mqr+mqr.^2)).*(x(i)-x(id));%1/3*(c^2./mqr+mqr-c^3./(c*mqr+mqr.^2)).*(x(i)-x(id));%
   f_y=1/3*(1./mqr+mqr-1./(mqr+mqr.^2)).*(y(i)-y(id));%1/3*(c^2./mqr+mqr-c^3./(c*mqr+mqr.^2)).*(y(i)-y(id));%
   f_p=bdx.*f_x+bdy.*f_y;
   f_v=mqr.*(dis2/9+4/(9*c^2))-log(mqr+1)/(3*c^2);% (4*c^2+dis2).*mqr/9-c^3*log(c+mqr)/3;   %
 end %end of if


    DM=zeros(ns,ns);
    for j=1:ns
     DM(j,:)=(x(id(j))-x(id(1:ns))).^2+(y(id(j))-y(id(1:ns))).^2;
    end


      
    mq=sqrt(DM*c^2+1);%sqrt(DM+c^2);%
    % Construct standard RBF matrix


     A=mq.*(DM./9+4/(9*c^2))-log(mq+1)/(3*c^2);%(4*c^2+DM).*mq/9-c^3*log(c+mq)/3;%
     
     if sort(i)==1
       A12=mq(:,idH);%sqrt(DM(:,idH)*c^2+1);
       A22=2*c^2./mq(idH,idH)-DM(idH,idH)*c^4./mq(idH,idH).^3;
       AH=[A A12; A12' A22];
       c_pH=AH\f_pH;
       A_p(i,id)=c_pH(1:ns);
       A_f(i,id(idH))=c_pH(ns+1:end);
     else
       c_p=A\f_p;   
       A_p(i,id)=c_p;
     end
     %%% Replace A_p with compute_weights(x, y, coor, int_ind, c, ns)
     
 
   c_v=A\f_v;
   c_x=A\f_x;
   c_y=A\f_y;


   % 
   % for k=1:ns %sparse matrix storage for interior points
   %   %A_p(i,id(k))=c_p(k);
   %   A_v(i,id(k))=c_v(k);
   %   A_x(i,id(k))=c_x(k);
   %   A_y(i,id(k))=c_y(k);
   % end 
   A_v(i,id)=c_v;
   A_x(i,id)=c_x;
   A_y(i,id)=c_y;
  
end %end of for i=1:n


return
