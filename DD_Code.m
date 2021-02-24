function [WE,NS]=rows_unique(Om,Op) 
%find unique rows of Om and common uniquencess of Om and Op
%Om and Op are diagonals
WE=cell(3,1);
NS=cell(3,1);

%Unique of Om
[~,I,J]=unique(Om);
WE{1}=numel(I);
WE{2}=I;
WE{3}=J;

%Unique rows of [Om Op]
[~,I,J]=unique([Om Op],'rows');
NS{1}=numel(I);
NS{2}=I;
NS{3}=J;
end

function Idk=Itensor(dk,Eta,Om,Op,WE,NS) 
    ns=NS{2};
    we=WE{2};
%EQS 14 15 Makri or Eq. 12 Lovett dk here same as in TEMPO
    if dk==0
    eta_dk=Eta(1);
    Idk=exp(-(real(eta_dk)*Om+1i*imag(eta_dk)*Op).*Om);
    Idk=Idk(ns);
    else
    eta_dk=Eta(dk+1);
    Idk=exp(-(real(eta_dk)*Om+1i*imag(eta_dk)*Op)*Om.');
%    Idk=Idk(ns,we).';  %in TEMPO they take transpose
    end
end


function tab=btensor(dk,freeprop,Eta,Om,Op,WE,NS) 
%converts rank-2 I_dk tensor into a 4-leg b tensor Eq.(22) taking account of degeneracy
%%need to see the use first to make sure indexing is OK
        ndim2=numel(Om);   %ndim2=ndim^2
        N0=NS{1};
        W0=WE{1};
        if dk==1
%           if dk=1 then multiply in free propagator - note we also include I_0 here instead of b_0 like in Methods section    
            v0=Itensor(0,Eta,Om,Op,WE,NS);
            V0=Itensor(1,Eta,Om,Op,WE,NS);
            for k=1:ndim2
                V0(:,k)=V0(:,k).*v0;  %This for loop is not needed in Python because of broadcasting
            end
            iffac=V0.*freeprop^2;
%           initialise 4-leg tensor dimensions based on degeneracy
%           for dk=1 we can only use the degeneracy/partial summing technique on South and East legs (see mpsmpo_class.py)
            tab=zeros(NS{0},ndim2,ndim2,WE{0});
%           loop through assigning elements of 4-leg from elements the 2-leg
            for i1=1:ndim2
                for a1=1:ndim2
                    tab(NS{3}(i1),i1,a1,WE{3}(a1))=iffac(i1,a1);
                end
            end
       
        else
%           start by constructing an array out of the unique elements in itab
            tab=(Itensor(dk,Eta,Om,Op,WE,NS)).';  %this will correspond to python
%           combine the 2 legs(axes) into vector and create matrix with this vector as diagonal
            tab=diag(tab(:));   %16^2x16^2=4^8
%           reshape the matrix into a 4-leg
            tab=reshape(tab,W0,N0,W0,N0);
            tab=permute(tab,[2 4 1 3]);
        if dk==point || (dk==dkmax && dkmax>0)
 %          if at an end site then sum over external leg and replace with 1d dummy leg
            tab=sum(tab,4);   %?
        end    
        end        

 en


main script-work in progress:

% op=A; must be diagonal; 
% op=A; must be diagonal;
sigz=[1 0;0 -1];
sigx=[0 1;1 0];
da=0.1;
sa=0.05;
alpha=0.05;
H0=((kron(sigx,eye(2))+(1-sa)*kron(eye(2),sigx)))/2;
A=(((1-da)*kron(sigz,eye(2))+kron(eye(2),sigz)))/2;
ndim=4;
omegac=10;
dt=0.01;
kmax=200;
Om=diag(kron(A,eye(ndim)) - kron(eye(ndim),A.')); 
Op=diag(kron(A,eye(ndim)) + kron(eye(ndim),A.')); 

[WE,NS]=rows_unique(Om,Op);

ham=-1i*(kron(H0,eye(ndim)) - kron(eye(ndim),H0.')); 
%         note dt/2 due to symmetrized trotter splitting
%         also note the transpose .T -- we are using the nonstanard convention of
%         propagating row vectors (instead of column), multiplying matrices from the right
freeprop=expm(ham*dt/2).';

% % % % % freeprop=expm(ham*dt/2);   
% % % % %note that ham acts on vectorized state as ordinary matrix multiplication,
% % % % %e.g.g left to right, in contrasts to what they did in Python. Thus there is
% % % % %no transposition of ham here. Use this later in anticipation of
% MATLAB

%discretize eta(t)
T=(0:kmax)'*dt;
BCF=0.5*alpha*omegac^2./(1+1i*omegac*T).^2;
Eta=zeros(kmax+1,1);
I0=(2:kmax+1)';T0=T(I0);
Eta(I0)=-0.5*alpha*log((1+1i*omegac*T0).^2./((1+1i*omegac*(T0-dt)).*...
    (1+1i*omegac*(T0+dt))));
Eta(1)=-0.5*alpha*(1i*omegac*dt-log(1+1i*omegac*dt));


%   these are allready integrated eta arrays, using Ohmic BCF

%   preparing tempo system
%   propagate initial state, multiply in I_0 and insert into
%   mps to create initial rank-1 ADT, as in Eq.(17)
%   Note only propagating init state dt/2 due to symmetric trotter splitting
%   and using expand dims to turn 1-leg init state into 3-leg tensor with 2 1d dummy indices

I_0=Itensor(0,Eta,Om,Op,WE,NS);
