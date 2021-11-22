%  This code is made by Mingye Ju
%
%  Mingye Ju:  
%  Nanjing University of Posts and Telecommunications, Nanjing, China.
%  Email: jumingye@njupt.edu.cn
%
%  WORK SETTING:
%  This code has been compiled and checked on MATLAB R2016b.
%  Please note that the code does not use downsampling operator to accelerate processing. 
%
%  For more detials, please refer to the paper:
%  Mingye Ju, Can Ding, Wenqi Ren, Yi Yang, Dengyin Zhang, Y.jay Guo. 
%  IDE: Image Dehazing and Exposure Using an Enhanced Atmospheric Scattering Model, 
%  IEEE TIP 2021.

clear all;close all 
sourcePic=double(imread('./test_images/realworld_5.png'))/256;
[m,n,o]=size(sourcePic);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%initialization
Graymean=0.5;p1=-0.397;p2=0.07774;Patch_size_GWA=75;T=0.02;%T can be adjusted to 0.01~0.05.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%airlight estimation
A_channel=Airlight(sourcePic,'our',10);
A_=mean(A_channel);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Impose GWA on ASM 0.05
Kernel=ones(Patch_size_GWA,Patch_size_GWA);%DE=sourcePic(:,:,3);
I_mean=sourcePic(:,:,1)+sourcePic(:,:,2)+sourcePic(:,:,3);
I_locolmean=imfilter(I_mean, Kernel,'symmetric')/sum(sum(Kernel))/3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Find minimum point 
a=0.001;b=1;Theta_error1=0.001;Theta_error2=0.001;N=1;
r=0.618;a1=b-r*(b-a);a2=a+r*(b-a);stepNum=0;
I_locolmean_downsample=imresize(I_locolmean,N);
sourcePic_downsample=imresize(sourcePic,N);
[m_downsample,n_downsample]=size(I_locolmean_downsample);
while abs(b-a)>Theta_error1
    stepNum=stepNum+1;
    f1=GS(a1,A_,I_locolmean_downsample,Graymean,p1,p2,sourcePic_downsample,A_channel,m_downsample,n_downsample);
    f2=GS(a2,A_,I_locolmean_downsample,Graymean,p1,p2,sourcePic_downsample,A_channel,m_downsample,n_downsample);
    if f1>f2
       a=a1;
       f1=f2;
       a1=a2;
       a2=a+r*(b-a);     
    else
       b=a2;
       a2=a1;
       f2=f1;
       a1=b-r*(b-a);       
    end 
   t_low_MIN=(a1+a2)/2;
end
Z=GS(t_low_MIN,A_,I_locolmean_downsample,Graymean,p1,p2,sourcePic_downsample,A_channel,m_downsample,n_downsample);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1-D Search
Bound_low=0.01;Bound_high=t_low_MIN;
interval=Bound_high-Bound_low;StepNUM2=1;
while interval>Theta_error2
    cen=Bound_high/2+Bound_low/2;
    if (GS(Bound_low,A_,I_locolmean_downsample,Graymean,p1,p2,sourcePic_downsample,A_channel,m_downsample,n)-Z-T)*(GS(cen,A_,I_locolmean_downsample,Graymean,p1,p2,sourcePic_downsample,A_channel,m_downsample,n_downsample)-Z-T)<0
    Bound_high=cen;
    elseif (GS(Bound_low,A_,I_locolmean_downsample,Graymean,p1,p2,sourcePic_downsample,A_channel,m_downsample,n_downsample)-Z-T)*(GS(cen,A_,I_locolmean_downsample,Graymean,p1,p2,sourcePic_downsample,A_channel,m_downsample,n_downsample)-Z-T)>0
    Bound_low=cen;
    else
        Bound_high=cen;Bound_low=cen;
    end
    interval=interval/2;StepNUM2=StepNUM2+1;
end
t_low_determined=Bound_high/2+Bound_low/2+0.02;
%t_low_determined=Fc(A_,I_locolmean,Graymean,p1,p2,sourcePic,A_channel,m,n,T). % You can see curve between t_min and saturation by executing this line.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Recover
A=log(t_low_determined)*A_;
B=log(t_low_determined).*I_locolmean-A_*Graymean*p1+log(t_low_determined).*A_*p2-log(t_low_determined)*A_;
C=p2*log(t_low_determined)*I_locolmean-log(t_low_determined)*A_*p2;
t=abs((-B-sqrt(B.^2-4*A.*C))./(2*A));
t_refine=guidedfilter(sourcePic(:,:,3),t,50, 0.02); 
% t_refine=guidedfilter(sourcePic(:,:,3),t_refine,30, 0.03);
t_refine=1.0*max(t_refine,t_low_determined);
ILCC=1*(p1./(t_refine+p2))./(p1/(t_low_determined+p2));
J(:,:,1) = ((sourcePic(:,:,1) -A_channel(1).*(1-t_refine))./t_refine./A_channel(1)./ILCC);
J(:,:,2) = ((sourcePic(:,:,2) -A_channel(2).*(1-t_refine))./t_refine./A_channel(2)./ILCC);
J(:,:,3) = ((sourcePic(:,:,3) -A_channel(3).*(1-t_refine))./t_refine./A_channel(3)./ILCC);
figure,
subplot(121),imshow(sourcePic);
subplot(122),imshow(J)




