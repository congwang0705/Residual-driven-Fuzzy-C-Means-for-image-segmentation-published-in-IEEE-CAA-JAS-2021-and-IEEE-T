function [E,V,U] = RFCM_color(IMG,C,error,m,lambda,lambda1,lambda2)
[height,width,~] = size(IMG);
II1 = IMG(:,:,1);
II2 = IMG(:,:,2);
II3 = IMG(:,:,3);

 N = height*width;
 data1 = reshape(II1,1,N);
 data_C1 = repmat(data1,[C,1]);
 data2 = reshape(II2,1,N);
 data_C2 = repmat(data2,[C,1]);
 data3 = reshape(II3,1,N);
 data_C3 = repmat(data3,[C,1]);
 t=0;
 U =rand(C,N);
 U = U./repmat(sum(U),[C,1]);
 E1 = zeros(size(data_C1));
 E2 = zeros(size(data_C2));
 E3 = zeros(size(data_C3));
 a = 1/(1+2.^0.5);
 b = 1/(1+1);
 H = [a,b,a;b,0,b;a,b,a]; 
 H1 = [a,b,a;b,1,b;a,b,a]; 
 lambda(1) = lambda(1)*sum(sum(H1));
 lambda(2) = lambda(2)*sum(sum(H1));
 lambda(3) = lambda(3)*sum(sum(H1));
 W1=ones(C,N);W2=ones(C,N);W3=ones(C,N);
 while t<100
     %Compute V1
     data1 = reshape((data_C1-E1).',[height,width,C]);
     data_nei1 = imfilter(data1,H,'replicate');
     data_neighbour1 = reshape(data_nei1,[N,C]);
     data_neighbour1 = data_neighbour1.';     
     V1 = sum(U.^m.*((data_C1-E1)+data_neighbour1),2)./sum((1+sum(sum(H)))*(U.^m),2);
     %Compute V2
     data2 = reshape((data_C2-E2).',[height,width,C]);
     data_nei2 = imfilter(data2,H,'replicate');
     data_neighbour2 = reshape(data_nei2,[N,C]);
     data_neighbour2 = data_neighbour2.';     
     V2 = sum(U.^m.*((data_C2-E2)+data_neighbour2),2)./sum((1+sum(sum(H)))*(U.^m),2);
     %Compute V3
     data3 = reshape((data_C3-E3).',[height,width,C]);
     data_nei3 = imfilter(data3,H,'replicate');
     data_neighbour3 = reshape(data_nei3,[N,C]);
     data_neighbour3 = data_neighbour3.';     
     V3 = sum(U.^m.*((data_C3-E3)+data_neighbour3),2)./sum((1+sum(sum(H)))*(U.^m),2);     
     %Compute G
     Xi_Vj = (data_C1-E1-repmat(V1,[1,N])).^2+(data_C2-E2-repmat(V2,[1,N])).^2+(data_C3-E3-repmat(V3,[1,N])).^2;
     Xi_Vj2 = (data_C1-E1-repmat(V1,[1,N])).^2+(data_C2-E2-repmat(V2,[1,N])).^2+(data_C3-E3-repmat(V3,[1,N])).^2;
     XI_VJ = reshape(Xi_Vj2.',[height,width,C]);
     g = imfilter(XI_VJ,H,'replicate');
     G = reshape(g,[N,C]);
     G = G.';   
     U_NEW = (1./(Xi_Vj+G)).^(1/(m-1))./repmat(sum((1./(Xi_Vj+G)).^(1/(m-1))),[C,1]); 
     % Fixed U and V, Update E 
     temp1 = reshape((U_NEW.^m).',[height,width,C]);
     g = imfilter(temp1,H1,'replicate');
     U_now = reshape(g,[N,C]);
     U_now = U_now.';
     M1 = sum(U_now.*(data_C1-repmat(V1,[1,N])));   
      E1 = soft_threshold(repmat(M1,[C,1]),lambda(1)/2);
      E1 = E1./(repmat(sum(U_now),[C,1])+2*lambda1(1)*W1.^2);
      W1=exp(-1*lambda2(1)*E1.^2);
      % Fixed U and V, Update E 
      M2 = sum(U_now.*(data_C2-repmat(V2,[1,N])));   
      E2 = soft_threshold(repmat(M2,[C,1]),lambda(2)/2);
      E2 = E2./(repmat(sum(U_now),[C,1])+2*lambda1(2)*W2.^2);
      W2=exp(-1*lambda2(2)*E2.^2);
      % Fixed U and V, Update E 
      M3 = sum(U_now.*(data_C3-repmat(V3,[1,N])));   
      E3 = soft_threshold(repmat(M3,[C,1]),lambda(3)/2);
      E3 = E3./(repmat(sum(U_now),[C,1])+2*lambda1(3)*W3.^2);
      W3=exp(-1*lambda2(3)*E3.^2);

     if max(max((abs(U_NEW-U))))<error
         break
     else
         U = U_NEW;
     end
     t=t+1;
 end
 E = [E1;E2;E3];
 V = [V1,V2,V3];
end