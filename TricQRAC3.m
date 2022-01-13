% This is the code for simulating the sequential random access codes
%   Ref: arXiv:1905.06726
%        PhysRevA.95.052345(2017)
% First, we simulate/check the original 3-1 quantum random access codes:
E = [1 0; 0 1]; X = [0 1; 1 0]; Y = [0 -1i; 1i 0]; Z = [1 0; 0 -1];
H = [1; 0]; V = [0; 1]; D = [1; 1]/sqrt(2); A = [1; -1]/sqrt(2);
L = [1; 1i]/sqrt(2); R = [1; -1i]/sqrt(2);
PY{3} = H*H'; PY{2} = D*D'; PY{1} = L*L';
PZ{3} = V*V'; PZ{2} = A*A'; PZ{1} = R*R';
xx = [0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1];
% consider that the measurements by Bob are unsharp, with the measurement strength \eta
et = 0:0.001:1;
for k = 1:length(et)
    PO{1} = X; PO{2} = Y; PO{3} = Z;
    eta = et(k);
    eta = eta^0.615;
    % next five lines of code defining the unsharp measurement performed by BOB
    for i = 1:3
        PPY{i} = PY{i}*eta+(1-eta)*E/2;
        PYY{i} = sqrt((1+eta)/2)*PY{i}+sqrt((1-eta)/2)*PZ{i};
        PZZ{i} = sqrt((1-eta)/2)*PY{i}+sqrt((1+eta)/2)*PZ{i};
    end

    seta = et(k);
    for i = 1:3
        SPPY{i} = PY{i}*seta+(1-seta)*E/2;
        SPYY{i} = sqrt((1+seta)/2)*PY{i}+sqrt((1-seta)/2)*PZ{i};
        SPZZ{i} = sqrt((1-seta)/2)*PY{i}+sqrt((1+seta)/2)*PZ{i};
    end

    % Here we compute the successful probability for all eight input states
    for i = 1:8
    % the state generation according the input 3-bit string x = x0/x1/x2
        x0 = xx(i,1); x1 = xx(i,2); x2 = xx(i,3);
        cs = sqrt(1/2+1/2/sqrt(3)*(-1)^x2); ss = sqrt(1/2-1/2/sqrt(3)*(-1)^x2); 
        phi = pi/4*(1+4*x0+2*(mod(x0+x1,2)));
        psi = [cs; ss*exp(1i*phi)];
        rhoaa = psi*psi';
        % the m-th bit reveal by Bob's weak measurement.
        for y = 1:3
            pp(1) = trace(SPYY{y}*rhoaa*(SPYY{y})'); %
            pp(2) = trace(SPZZ{y}*rhoaa*(SPZZ{y})');
            P_Y(i,y) = pp(xx(i,y)+1);
            Bl(i,y) = trace(psi*psi'*PO{y})*sqrt(3);
        end
        % the next five lines generate the states out from Bob's device
        rhobb = E*0;
        for y = 1:3
            rhobb = rhobb+SPYY{y}*rhoaa*(SPYY{y})'+ SPZZ{y}*rhoaa*(SPZZ{y})';
        end
        rhobb = rhobb/3;
        % measurement by Charlie    
        for z = 1:3
            pp(1) = trace(PYY{z}*rhobb*(PYY{z})'); %
            pp(2) = trace(PZZ{z}*rhobb*(PZZ{z})');
            P_Z(i,z) = pp(xx(i,z)+1);
        end
        rhocc = E*0;
        for y = 1:3
            rhocc = rhocc+PYY{y}*rhobb*(PYY{y})'+PZZ{y}*rhobb*(PZZ{y})';
        end        
        rhocc = rhocc/3;    
        % measurement by Dave 
        for z = 1:3
            pp(1) = trace(PY{z}*rhocc*(PY{z})'); %
            pp(2) = trace(PZ{z}*rhocc*(PZ{z})');
            P_W(i,z) = pp(xx(i,z)+1);
        end
    end
    Psus(k) = abs(sum(sum(P_Y))/24);
    Pzus(k) = abs(sum(sum(P_Z))/24);
    Pwus(k) = abs(sum(sum(P_W))/24);
end
Psus;  % Bob's average success probabilities
Pzus;  % Charlie's average success probabilities
Pwus;  % Dave's average success probabilities
figure(1); 
plot(et,Psus,'r');hold on; 
plot(et,Pzus,'b-.');hold on;
plot(et,Pwus,'k:','LineWidth',2);hold on;  
plot(et,et*0+2/3,'--');
xlabel('¦Ç_1');
ylabel('Probability');
legend('Bob','Charlie','Dave','0.667');

