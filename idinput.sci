function[u]=idinput(N,types,band,levels)
// generates random binary input signal
// 
// Calling Seqence
// u=idinput(N);
// u=idinput([n,nu])
// u=idinput([n,nu,m])
// u=idinput(__,type)
// u=idinput(__,type,band)
// u=idinput(__,type,band,levels)

// 
// Parameters
// N       : no of binary signals 
// [n,nu]  : nu channel signal of length n.
// [n,nu,m]: nu channel signal of length nxm and period m.
// type    : "rbs" ,"prbs"
// band    : frequency band of the signal.
// levels  : binary output values
//
// Description   
// u = idinput(n) returns a single-channel random binary signal u of length N. By default it’s values are either -1 or 1.
// u = idinput([n,nu]) returns an nu-channel random binary input signal, where each channel signal has length n,the 
// signals in each channel are independent from each other.
// u = idinput([n,nu,m]) returns an Nu-channel periodic random binary input signal with specified period and number of 
// periods. Each input channel signal is of length m*n.
// u=idinput(___,type) specifies the type of input to be generated:
        // 'rbs' — Random binary signal.
        // 'prbs'—pseudorandom binary signal. 
// u =idinput(___,type,band) specifies the frequency band of the signal.
// u =idinput(___,type,band,levels) specifies the level of the binary generated signal
//
//    
// Examples
// u=idinput([20,2,2],'rbs',[0 0.3],[-1 2]);
// 
// Author
// Ayush Kumar
    
    
    [lhs,rhs]=argn(0)
    if  rhs<4 then                    // checking the input given with function 
        levels=[-1 1]                 // arguments and allocating default value if required
    end
    if  rhs<3 then
        band=[0 0.5]
    end
    if  rhs<2 then
        types='rbs'
    end
    if  size(N,2)==3 then
        P=N(1);
        nu=N(2);
        M=N(3);
    elseif  size(N,2)==2 then
        P=N(1);
        nu=N(2);
        M=1;
    elseif  size(N,2)==1 then
        P=N;
        nu=1;
        M=1;
    else
        error('erroininputargument');
    end

    if levels(2)<levels(1) then           // error checking
        error('levels (1)should be less then levels(2)');
        end
    if type(types)~=10 then
        error('types can be rbs only input valid argument');
        end
    if(band(2)>1) then
        error('band should be in b/w 0-1');
    end
    
    function[vec]=randnu(NN,nuu)   // nuu channel random normally distributed signal generator
    [lhs,rhs]=argn(0)
    if  rhs<2 then                
        nuu=1                      
    end
    rand('seed',getdate('s'))
    for j=1:1:nuu
        for i= 1:NN
            vec(i,j)=rand(NN,'normal')
        end
    end
    endfunction


    if convstr(types)=='rbs' then
        u=randnu(5*P,nu);
        band=band/2;
        if ~and(band==[0 0.5]) then
            if(band(1)==0) then
                [hz]=iir(8,'lp','butt',[band(2) band(1)],[0 0]);  // 8th order butterwoth filter
                num=hz(2);
                den=hz(3);
                for i=1:1:nu
                    y(:,i)=filter(num,den,u(:,i));
                end
            elseif(band(2)==0.5) then
                [hz]=iir(8,'hp','butt',[band(1) band(2)],[0 0]);  
                num=hz(2);
                den=hz(3);
                for i=1:1:nu
                    y(:,i)=filter(num,den,u(:,i));
                end
            else
                [hz]=iir(8,'bp','butt',band,[0 0]);  
                num=hz(2);
                den=hz(3);
                for i=1:1:nu
                    y(:,i)=filter(num,den,u(:,i));
                end
            end
            
            
            u = sign(y(2*P+1:3*P,:)); // taking the middle terms
            u = (levels(2)-levels(1))*(u+1)/2+levels(1);// adjusting binary values according to levels
        else
            u = sign(u(2*P+1:3*P,:));
            u = (levels(2)-levels(1))*(u+1)/2+levels(1);
            
        end
    elseif convstr(types)=='prbs' then
        clockP = floor(1/band(2));
        possP = 2.^(3:18)-1;
        P1 = max(possP(P/clockP-possP>=0));
        if isempty(P1)
            P1 = 7;
        end
        n = find(P1==possP)+2;
   
        if (clockP*P1~=P)
            if M>1
                warning(msprintf(gettext("%s: period of prbs signal changed to %d and length change to %d"),"idinput",clockP*P1,P1*M*clockP));
                P = P1*clockP;
            else
                n = min(n+1,18);
                P1 = 2^n -1;
                warning(msprintf(gettext("%s :generated signal is first %d values of a sequence of length %d"),"idinput",P,clockP*P1));
            end
        end
   
        P1 = 2^n-1;
   
        if n<3 || n>18,
            error(msprintf(gettext('Bad conditioning')));
        end
        fi = -ones(n,1);
        if n==3
            ind = [1,3];
        elseif n==4
            ind = [1,4];
        elseif n==5
            ind = [2,5];
        elseif n==6
            ind = [1,6];
        elseif n==7
            ind = [1,7];
        elseif n==8
            ind = [1,2,7,8];
        elseif n==9
            ind = [4,9];
        elseif n==10
            ind = [3,10];
        elseif n==11
            ind = [9,11];
        elseif n==12
            ind = [6,8,11,12];
        elseif n==13
            ind = [9,10,12,13];
        elseif n==14
            ind = [4,8,13,14];
        elseif n==15
            ind = [14,15];
        elseif n==16
            ind = [4,13,15,16];
        elseif n==17
            ind = [14,17];
        elseif n==18
            ind = [11,18];
        end
        for t = 1:clockP:P1*clockP
            u(t:t+clockP-1,1) = ones(clockP,1)*fi(n); // LL%% multivariable !!
            fi = [prod(fi(ind));fi(1:n-1,1)];
        end
        u = (levels(2)-levels(1))*(u+1)/2+levels(1);
   
        u = u(1:P,1);
        if nu >1
            u1 = [u;u];
            shift = floor(P/nu);
            for ku = 2:nu
                u = [u,u1(shift*(ku-1)+1:P+shift*(ku-1))];
            end
        end
    else
        error(msprintf(gettext('type can be rbs or prbs only')));
    end
    if M>1 then              // generating periodic input if no. of periods>1
        uu = u;
        for i = 2:M
            u = [uu;u];
        end
    end
endfunction

    

    
    
