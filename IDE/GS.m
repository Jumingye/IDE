function SA=GS(t_low,A_,Ave_a,uu,p1,p2,sourcePic,A_channel,m,n)
a=log(t_low)*A_;
b=log(t_low).*Ave_a-A_*uu*p1+log(t_low).*A_*p2-log(t_low)*A_;
c=p2*log(t_low)*Ave_a-log(t_low)*A_*p2;
t=abs((-b-sqrt(b.^2-4*a.*c))./(2*a));
ILCC=(p1./(t+p2))./(p1/(t_low+p2));
J(:,:,1) = double((sourcePic(:,:,1) -A_channel(1).*(1-t))./t./A_channel(1)./ILCC);
J(:,:,2) = double((sourcePic(:,:,2) -A_channel(2).*(1-t))./t./A_channel(2)./ILCC);
J(:,:,3) = double((sourcePic(:,:,3) -A_channel(3).*(1-t))./t./A_channel(3)./ILCC);%figure,imshow(J)
J=min(max(J,0),1);
SA=(sum(J(:)==0)+sum(J(:)==1))/m/n/3;
% SA=(sum(J(:)==0))/m/n/3;