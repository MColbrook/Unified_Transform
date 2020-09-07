function [ max_err,C] = max_error1( Col,tvec,M )
% this code plots figure 3 from fornberg's version for different types of
% collocations points

L=length(tvec);
z  = [0,1,1+2i,1i];    % The four corners 

max_err=zeros(L,1);
C=max_err;
n=4;

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
    else
        k=[];
        for ii=1:n
            if ii<n
                hj=z(ii+1)-z(ii);
                mj=z(ii+1)+z(ii);
            else
                hj=z(1)-z(4);
                mj=z(1)+z(4);
            end
            hj=hj/2;
            mj=mj/2;
            lvec=linspace(0,1,round(Col.K/4))*(Col.R2-Col.R1)+Col.R1;
            k=[k,-lvec/hj];
        end
    end

    % Create the matrices that form the system matrix

    [RD1,RN1,SD1,SN1] = AB(z(1),z(2),k,M);          % Side 1
    [RD2,RN2,SD2,SN2] = AB(z(2),z(3),k,M);          % Side 2
    [RD3,RN3,SD3,SN3] = AB(z(3),z(4),k,M);          % Side 3
    [RD4,RN4,SD4,SN4] = AB(z(4),z(1),k,M);          % Side 4
    
    % Compute the Legendre coefficients that should satisfy Test Problem 1:
    % f(z) = real(exp(1+2i+z)) = exp(1+x)cos(2+y); again one side at a time
    a = 1+2i;   b = 1;
    [qd1,qn1] = BV(z(1),z(2),a,b,M); qd1 = real(qd1); qn1 = real(qn1);
    [qd2,qn2] = BV(z(2),z(3),a,b,M); qd2 = real(qd2); qn2 = real(qn2);
    [qd3,qn3] = BV(z(3),z(4),a,b,M); qd3 = real(qd3); qn3 = real(qn3);
    [qd4,qn4] = BV(z(4),z(1),a,b,M); qd4 = real(qd4); qn4 = real(qn4);
    
    % Form the linear system that should be satisfied
    if strcmp(Col.collocation_type,'halton')
        AA =  [ RN1, RN2, RD3, RN4; ... % Coefficient matrix
                SN1, SN2, SD3, SN4];
        BB = -[ RD1, RD2, RN3, RD4; ...
                SD1, SD2, SN3, SD4];
    else
        AA =  [ RN1, RN2, RD3, RN4];
        AA=[AA;conj(AA)];
        BB = -[ RD1, RD2, RN3, RD4];
        BB=[BB;conj(BB)];
    end
    rA =  [ qn1; qn2; qd3; qn4];
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
    toc
end

end

