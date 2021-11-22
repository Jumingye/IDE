function t_low_determined=Fc(A_,Ave_a,uu,p1,p2,sourcePic,A_channel,m,n,S)
iter=1;
for t_low=0.01:0.01:1
a=log(t_low)*A_;
b=log(t_low).*Ave_a-A_*uu*p1+log(t_low).*A_*p2-log(t_low)*A_;
c=p2*log(t_low)*Ave_a-log(t_low)*A_*p2;
t=abs((-b-sqrt(b.^2-4*a.*c))./(2*a));
ILCC=(p1./(t+p2))./(p1/(t_low+p2));
t_refine=t;
J(:,:,1) = double((sourcePic(:,:,1) -A_channel(1).*(1-t_refine))./t_refine./A_channel(1)./ILCC);
J(:,:,2) = double((sourcePic(:,:,2) -A_channel(2).*(1-t_refine))./t_refine./A_channel(2)./ILCC);
J(:,:,3) = double((sourcePic(:,:,3) -A_channel(3).*(1-t_refine))./t_refine./A_channel(3)./ILCC);%figure,imshow(J)
J=min(max(J,0),1);
SA(iter)=(sum(J(:)==0)+sum(J(:)==1))/m/n/3;
% SA(iter)=(sum(J(:)==0))/m/n/3;
iter=iter+1;
end
X=0.01:0.01:1;
figure,plot(X,SA)
POS2=find(abs(SA)==min(abs(SA)));
POS=find(abs(SA(1:POS2)-S)==min(abs(SA(1:POS2)-S)));
t_low_determined=POS(1)*0.01;
end