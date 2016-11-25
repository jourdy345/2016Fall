% clear; 
% clc;


p = 4;
mlag = 24;

Data = xlsread('watson.xls',1,'A1:C167'); 

Y = centre(Data);
k = size(Y,2);
n0 = 500;
n1 = 2000;
n = n0 + n1;

nu = 10;
R0 = invpd(0.2*eye(k))/nu; 

pkk = p*k*k;
b_ = zeros(pkk,1);
var_ = 1*eye(pkk);

Spec.b_ = b_;
Spec.var_ = var_;
Spec.p = p;
Spec.nu = nu;
Spec.R0 = R0;
Spec.Y = Y;
Spec.mlag = mlag;

% [ImpulseRespm, MHm]= recursive_VAR(n0, n1, Spec);

Y_before80 = Y(1:80,:);
Y_after80  = Y(81:end,:);
Spec_before80 = Spec;
Spec_before80.Y = Y_before80;
Spec_after80 = Spec;
Spec_after80.Y = Y_after80;

[ImpulseRespmbefore80,MHmbefore80] = recursive_VAR(n0,n1,Spec_before80);
% [ImpulseRespmafter80,MHmafter80] = recursive_VAR(n0,n1,Spec_after80);

   
   