clc, clear all, close all

%%

F = [1 0.5;
    0.5 1];

cvx_begin %sdp
    cvx_solver sedumi;
    variable Xsdp(2,2) symmetric;
    expression Zsdp;
    Zsdp = trace(F*Xsdp);

    minimize(Zsdp)
    subject to
    
    trace([1 0; 0 0]*Xsdp) <= 16
    trace([1 0; 0 0]*Xsdp) >= 4
    
%     trace([0 0; 1 0]*Xsdp) <= 10
%     -trace([0 0; 1 0]*Xsdp) <= 10

    trace([0 1; 0 0]*Xsdp) <= 10
    -trace([0 1; 0 0]*Xsdp) <= 10
    
    % Feeder Head Voltage
    Xsdp == semidefinite(2);
    
cvx_end;

%%

F = blkdiag(F,0);

cvx_begin %sdp
    cvx_solver sedumi;
    variable Xsdp(3,3) symmetric;
    expression Zsdp;
    Zsdp = trace(F*Xsdp);

    minimize(Zsdp)
    subject to
    
    trace([0 0 1; 0 0 0; 0 0 0]*Xsdp) <= 4
    trace([0 0 1; 0 0 0; 0 0 0]*Xsdp) >= 2
    
%     trace([0 0; 1 0]*Xsdp) <= 10
%     -trace([0 0; 1 0]*Xsdp) <= 10

    trace([0 1 0; 0 0 0; 0 0 0]*Xsdp) <= 10
    -trace([0 1 0; 0 0 0; 0 0 0]*Xsdp) <= 10
    
    trace([0 0 0; 0 0 0; 0 0 1]*Xsdp) == 1
    
    % Feeder Head Voltage
    Xsdp == semidefinite(3);
    
cvx_end;

[X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(@myfun,[1;1],[],[],[],[],[],[],@mycon)