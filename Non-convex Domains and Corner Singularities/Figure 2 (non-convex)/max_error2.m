function [ max_err,C ] = max_error2( Col,tvec, M )
% this code plots figure 3 from fornberg's version for different types of
% collocations points

L=length(tvec);
z  = [0,1,1+2i,1i];    % The four corners 

max_err=zeros(L,1);
C=max_err;

for jj=1:L
    tic
    
    z(4)=tvec(jj)+1i*(1-tvec(jj));


    if strcmp(Col.collocation_type,'halton') % use Halton nodes
        
        % Create the list of kn complex k-points (column vector)
        % Choose here Halton points within the radius R of the origin 
        H   = haltonset(2);  % 2-D node set
        kr  = 2*H(1:2*Col.K,:)-1;
        kt  = Col.R*complex(kr(:,1),kr(:,2)); % Convert to complex
        k   = kt(find(abs(kt)<Col.R,Col.K));     % This completes the k-point calculation
              
        
        % Create the matrices that form the system matrix
        [RD1,RN1,SD1,SN1] = AB(z(1),z(2),k,M);          % Side 1
        [RD2,RN2,SD2,SN2] = AB(z(2),z(3),k,M);          % Side 2
        [RD3,RN3,SD3,SN3] = AB(z(3),z(4),k,M);          % Side 3
        [RD4,RN4,SD4,SN4] = AB(z(4),z(1),k,M);          % Side 4
        [RD5,RN5,SD5,SN5] = AB(z(2),z(4),k,M);          % Side 5

        % Compute the Legendre coefficients that should satisfy Test Problem 1:
        % f(z) = real(exp(1+2i+z)) = exp(1+x)cos(2+y); again one side at a time
        a = 1+2i;   b = 1;
        [qd1,qn1] = BV(z(1),z(2),a,b,M); qd1 = real(qd1); qn1 = real(qn1);
        [qd2,qn2] = BV(z(2),z(3),a,b,M); qd2 = real(qd2); qn2 = real(qn2);
        [qd3,qn3] = BV(z(3),z(4),a,b,M); qd3 = real(qd3); qn3 = real(qn3);
        [qd4,qn4] = BV(z(4),z(1),a,b,M); qd4 = real(qd4); qn4 = real(qn4);
        [qd5,qn5] = BV(z(2),z(4),a,b,M); qd5 = real(qd5); qn5 = real(qn5);

        ZZZ = zeros(Col.K,M);
        % Form the linear system that should be satisfied 
        AA =  [ RN1, ZZZ, ZZZ, RN4, RD5, RN5; ... % Coefficient matrix
                SN1, ZZZ, ZZZ, SN4, SD5, SN5;...
                ZZZ, RN2, RD3, ZZZ, -RD5, -RN5;...
                ZZZ, SN2, SD3, ZZZ, -SD5, -SN5];
        BB = -[ RD1, ZZZ, ZZZ, RD4; ...
                SD1, ZZZ, ZZZ, SD4;...
                ZZZ, RD2, RN3, ZZZ;...
                ZZZ, SD2, SN3, ZZZ];
        rA =  [ qn1; qn2; qd3; qn4; qd5; qn5];
        rB =  [ qd1; qd2; qn3; qd4];

        sc = 1./sum(abs(AA),2);

        % Scale the system; make each row of the AA matrix have 1-norm equal to one
        AA  = bsxfun(@times,sc,AA);
        BB  = bsxfun(@times,sc,BB);

        % Verify that the resudial becomes small
        res = AA*rA-BB*rB; max(abs(res));

        RH = BB*rB;
        sol = AA\RH;
        nc = 4*M;
        max_err(jj) = (max(abs(sol(1:nc)-rA(1:nc))));
        C(jj)=cond(AA);

    else
        k1=[];
        k2=k1;
        for ii=1:3
            if ii==1
                hj1=z(2)-z(1);
                hj2=z(3)-z(2);               
            elseif ii==2
                hj1=z(4)-z(2);
                hj2=z(4)-z(3);
            else
                hj1=z(1)-z(4);
                hj2=z(2)-z(4);
            end
                
            hj1=hj1/2;
            hj2=hj2/2;
            lvec=linspace(0,1,round(Col.K/3))*(Col.R2-Col.R1)+Col.R1;
            k1=[k1,-lvec/hj1];
            k2=[k2,-lvec/hj2];
        end
        
        
        
        % Create the matrices that form the system matrix
        [RD1,RN1,~,~] = AB(z(1),z(2),k1,M);          % Side 1
        [RD2,RN2,~,~] = AB(z(2),z(3),k2,M);          % Side 2
        [RD3,RN3,~,~] = AB(z(3),z(4),k2,M);          % Side 3
        [RD4,RN4,~,~] = AB(z(4),z(1),k1,M);          % Side 4
        [RD5,RN5,~,~] = AB(z(2),z(4),k1,M);          % Side 5
        [RD52,RN52,SD52,SN52] = AB(z(4),z(2),k2,M);          % Side 5

        % Compute the Legendre coefficients that should satisfy Test Problem 1:
        % f(z) = real(exp(1+2i+z)) = exp(1+x)cos(2+y); again one side at a time
        a = 1+2i;   b = 1;
        [qd1,qn1] = BV(z(1),z(2),a,b,M); qd1 = real(qd1); qn1 = real(qn1);
        [qd2,qn2] = BV(z(2),z(3),a,b,M); qd2 = real(qd2); qn2 = real(qn2);
        [qd3,qn3] = BV(z(3),z(4),a,b,M); qd3 = real(qd3); qn3 = real(qn3);
        [qd4,qn4] = BV(z(4),z(1),a,b,M); qd4 = real(qd4); qn4 = real(qn4);
        [qd5,qn5] = BV(z(2),z(4),a,b,M); qd5 = real(qd5); qn5 = real(qn5);

        ZZZ = zeros(3*round(Col.K/3),M);
        
        MATCH = (-1).^(0:(M-1));
        MATCH = repmat(MATCH,3*round(Col.K/3),1);
        
        % Form the linear system that should be satisfied 
        AA =  [ RN1, ZZZ, ZZZ, RN4, RD5, RN5; ... % Coefficient matrix
                conj(RN1), ZZZ, ZZZ, conj(RN4), conj(RD5), conj(RN5);...
                ZZZ, RN2, RD3, ZZZ, MATCH.*RD52, -MATCH.*RN52;...
                ZZZ, conj(RN2), conj(RD3), ZZZ, MATCH.*conj(RD52), -MATCH.*conj(RN52)];
        BB = -[ (RD1), ZZZ, ZZZ, (RD4); ...
                conj(RD1), ZZZ, ZZZ, conj(RD4);...
                ZZZ, RD2, RN3, ZZZ;...
                ZZZ, conj(RD2), conj(RN3), ZZZ];
        rA =  [ qn1; qn2; qd3; qn4; qd5; qn5];
        rB =  [ qd1; qd2; qn3; qd4];

        sc = 1./sum(abs(AA),2);

        % Scale the system; make each row of the AA matrix have 1-norm equal to one
        AA  = bsxfun(@times,sc,AA);
        BB  = bsxfun(@times,sc,BB);

        % Verify that the resudial becomes small
        res = AA*rA-BB*rB; max(abs(res));

        RH = BB*rB;
        sol = AA\RH;
        nc = 4*M;
        max_err(jj) = (max(abs(sol(1:nc)-rA(1:nc))));
        C(jj)=cond(AA);
    end
    toc

end

end
