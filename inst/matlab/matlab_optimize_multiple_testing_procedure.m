function matlab_optimize_multiple_testing_procedure()

load A11.mat
A1   = tmp;
number_A1_files = importdata('number_A1_files.txt');
for i = 2:number_A1_files
    eval(['load A1',num2str(i),'.mat;']);
    A1 = [A1;tmp];


end

load A2.mat
A2  = sparse(A2(:,1),A2(:,2),A2(:,3));

load A3.mat


[~,col] = size(A2);
%[r4,tmp] = size(A4);

%if(tmp<col)
%    A4(r4,col)=0;
%end

load alphaValue.mat
load iteration.mat

[tmp,~] = size(A1);
a1      = alphaValue*ones(tmp,1);

load a21.mat
[tmp,~] = size(A2);
a2      = [ones(a21,1);zeros(tmp-a21,1)];

%[tmp,~] = size(A4);
%a4      = zeros(tmp,1);



load a3.mat
power = a3;

% AA      = [A1;A4;-A3];

AA      = [A1;-A3];


[tmpA1,~]   = size(A1);
%[tmpA4,~]   = size(A4);
%a3s     = tmpA1 + tmpA4 + 1;
% aa      = [a1;a4;-power];

aa = [a1;-power];

load a4status.mat

if (a4status == 1)

   load A4.mat
   A4  = sparse(A4(:,1),A4(:,2),A4(:,3));
   [~,col] = size(A2);
   [r4,tmp] = size(A4);

   if(tmp<col)
        A4(r4,col)=0;
   end

   AA  = [AA;A4];
   a4  = zeros(r4,1);

   aa  = [aa;a4];
end

cc = -(A3(3,:)+A3(6,:));

%load cc.mat
tmp     = length(cc);
lb      = zeros(tmp,1);
ub      = ones(tmp,1);


%options = cplexoptimset('lpmethod',0);

%options.barrier.crossover = -1;
%options.barrier.convergetol = 1e-10;

%options = cplexoptimset('Algorithm','primal');

%[z,val,status,output,dual] = cplexlp(cc,AA,aa,A2,a2,lb,ub,'options',options);

%dual = [dual.ineqlin;dual.eqlin];
%eval(['save sln1M',num2str(power),'.mat',' output z val status dual'])

status = 1;
if (status==1)||(status==5)
    load A5.mat
    A5   = sparse(A5(:,1),A5(:,2),A5(:,3));
    [~,col] = size(A2);
    [r5,tmp] = size(A5);

    if(tmp<col)
        A5(r5,col)=0;
    end
    aa       = [aa;zeros(r5,1)];
    AAn       = [AA;A5];
    [z,val,status,output,dual] = linprog(cc,AAn,aa,A2,a2,lb,ub);
    if(iteration<5)
    dual = [dual.ineqlin;dual.eqlin];
    end

    eval(['save sln2M',num2str(iteration),'.mat',' output z val status dual'])

end



