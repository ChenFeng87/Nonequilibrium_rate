function x = euler_simD_r(IM,AB,S,n,K, x0, save_step, h, d)
%IM: Interaction matrix.
%AB:the sthength of activation and inhibution.
%S:Threshold matrix.
%n:Hill function matrix.
time_n=save_step+1;
[num,m]=size(x0);
x=zeros(time_n,num);   %%the process
x(1,:)=x0;

for os=2:time_n
    temp=x(os-1,:);
    H=zeros(num,m,num);
    for i=1:num
        for j=1:num
            if IM(i,j)==1
                H(i,:,j)=Ha(temp(:,i),S(i,j),AB(i,j),n(i,j));
            end
            if IM(i,j)==-1
                H(i,:,j)=Hr(temp(:,i),S(i,j),AB(i,j),n(i,j));
            end
        end
    end
    for i=1:num
        temp(:,i)=temp(:,i) + h * (sum(H(:,:,i))-K(i)*temp(:,i)) + sqrt(2 * d) * sqrt(h) * randn(1);
%         F(i,:)=sum(H(:,:,i))-K(i)*x(i,:);
        if temp(:,i)<0
            temp(:,i)=0;
        end
    end
    x(os,:) = temp;
end

end
    
% A=x(1,:);
% B=x(2,:);
% C=x(3,:);
% a=0.7;b=0.7;S=0.5;n=4;k=1;
% F(1,:)=Ha(A,S,a,n)+Hr(B,S,b,n)-k*A;
% F(2,:)=Hr(C,S,b,n)+Hr(A,S,b,n)-k*B;
% F(3,:)=Ha(A,S,a,n)-k*C;
function H=Hr(X,S,b,n)
H=b*S.^n./(X.^n+S.^n);
end
function H=Ha(X,S,a,n)
H=a*X.^n./(S.^n+X.^n);
end