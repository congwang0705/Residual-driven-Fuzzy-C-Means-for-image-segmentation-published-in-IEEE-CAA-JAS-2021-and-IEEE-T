function [E,W,V,U,t,obj,J,diff] = RFCM(IMG,Noi,C,error,m,lambda,lambda1,lambda2)
[height,width] = size(IMG);
 N = height*width;
 data = reshape(IMG,1,N);
 data_C = repmat(data,[C,1]);
 t=0;
 U =rand(C,N);
 U = U./repmat(sum(U),[C,1]);
 E = zeros(size(data_C));
 a = 1/(1+2.^0.5);
 b = 1/(1+1);
 H = [a,b,a;b,0,b;a,b,a]; 
 H1 =[a,b,a;b,1,b;a,b,a]; 
 lambda = lambda*sum(sum(H1));
 W=ones(C,N);
 obj=[];J_prev = inf; J = [];
 while t<100000
     data1 = reshape((data_C-E).',[height,width,C]);
     data_nei = imfilter(data1,H,'replicate');
     data_neighbour = reshape(data_nei,[N,C]);
     data_neighbour = data_neighbour.';     
     V = sum(U.^m.*((data_C-E)+data_neighbour),2)./sum((1+sum(sum(H)))*(U.^m),2);
     Xi_Vj = (data_C-E-repmat(V,[1,N])).^2;
     XI_VJ = reshape(Xi_Vj.',[height,width,C]);
     g = imfilter(XI_VJ,H,'replicate');
     G = reshape(g,[N,C]);
     G = G.';  
     U_NEW = (1./(Xi_Vj+G)).^(1/(m-1))./repmat(sum((1./(Xi_Vj+G)).^(1/(m-1))),[C,1]); 
     temp1 = reshape((U_NEW.^m).',[height,width,C]);
     g = imfilter(temp1,H1,'replicate');
     U_now = reshape(g,[N,C]);
     U_now = U_now.';
     M = sum(U_now.*(data_C-repmat(V,[1,N])));   
     E = soft_threshold(repmat(M,[C,1]),lambda/2);
     E = E./(repmat(sum(U_now),[C,1])+2*lambda1*W.^2);
     W=exp(-1*lambda2*E.^2);
%%%%%%%%%    
     diff(t+1)=norm(E(1,:)'-Noi,2);
%%%%%%%%%          
     X1=data-E(1,:);X1=X1';
     dist = sum(X1.*X1, 2)*ones(1, C) + (sum(V.*V, 2)*ones(1, N))'-2*X1*V';
     J_cur = sum(sum((U_NEW'.^2).*dist))+lambda1*(norm(W(1,:).*E(1,:),2))^2;
     J = [J J_cur];
     obj(t+1)=norm(U_NEW-U,2);
     if obj(t+1)/N<error
         break
     else
         U = U_NEW;
     end
     t=t+1;
     
     display(sprintf('#iteration: %03d, objective function: %f', t, J_cur));
     J_prev = J_cur;
 end