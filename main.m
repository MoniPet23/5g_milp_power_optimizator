close all;
clear all; %#ok<CLALL>
clc;

%% Initiate the variables
CimaxMACRO = 100;
CimaxSMALL = 100;
Cimax = 100;
NTXan = 8;
NTXbh = 8;
p0ANMACRO = 130;
p0ANSMALL = 6.8;
deltapMACRO = 4.7;
deltapSMALL = 4.0;
p0BH = 3.9;
deltapBH = 100000;
pmaxBH = 0.0631; %1.9954;
M1 = 1000;
M2 = 1000;
M3 = 1000;
M4 = 1000;

%% retrieve data from file
% all_input = fileread('input_file_N=13_time=02_Snap01.dat');
%all_input = fileread('provascen1_new.dat');
% all_input = fileread('provascen1_fortest.dat');
% all_input = fileread('provascen1.dat');
% all_input = fileread('2610_test.dat');
%all_input = fileread('test_con_2_user.dat');
all_input = fileread('test3user.dat');
% all_input = fileread('test_without1vn.dat');
% all_input = fileread('new_scenario.dat');
%all_input = fileread('scenario_2.dat');

curly_separated_content = regexp(all_input,'{(.*?)}','match');
curly_separated_content = regexprep(curly_separated_content,'{', '');
curly_separated_content = regexprep(curly_separated_content,'}', '');
curly_separated_content = regexprep(curly_separated_content,'<', '');
curly_separated_content = regexprep(curly_separated_content,'>,', '');

purged_content = {};
for raw_content_line = curly_separated_content
   group = [];
   separated = splitlines(strtrim(raw_content_line));
   for line = 1:1:length(separated)
       group = [group; str2double(strsplit(separated{line},','))];
   end
   purged_content{end+1} = group;
   
end

ANfinal = purged_content{1};
BHfinal = purged_content{2};
Ufinal = purged_content{3};
Dfinal = purged_content{4};

NumUEs_raw = regexp(all_input,'NumUEs(.*?);','match');
NumUEs = remove_val_identifiers('NumUEs', NumUEs_raw);

NumVirtualNodes_raw = regexp(all_input,'NumVirtualNodes(.*?);','match');
VirtNode = remove_val_identifiers('NumVirtualNodes', NumVirtualNodes_raw);

Macro_raw = regexp(all_input,'NumMacro(.*?);','match');
NumMacro = remove_val_identifiers('NumMacro', Macro_raw);

Small_raw = regexp(all_input,'NumSmall(.*?);','match');
NumSmall = remove_val_identifiers('NumSmall', Small_raw);

numbkpt_raw = regexp(all_input,'numbreak(.*?);','match');
numbkpt = remove_val_identifiers('numbreak', numbkpt_raw);

breakpt_raw = regexp(all_input,'breakpt(.*?);','match');
breakpt = remove_val_identifiers('breakpt', breakpt_raw);

slope_raw = regexp(all_input,'coeff1(.*?);','match');
slope = remove_val_identifiers('coeff1', slope_raw);

%% Compute the Binary constraints

% dimensions AN or BH
NumBs = NumMacro + NumSmall; 
NumBs_BH = NumMacro + NumSmall + VirtNode; 

% constraints active links
xu = binvar(NumBs,NumUEs,'full'); % connessione AN e UE
xv = binvar(NumBs_BH,NumBs_BH,NumUEs,'full'); % matrice dei link sul BH

for u = 1 : 1 : NumUEs
    for i = 1 : 1 : NumBs_BH
        xv(i,i,u) = 0;
    end
end

% constraint BSi or BH(i,j)
sAN = binvar(NumBs,1,'full');
sBH = binvar(NumBs_BH,NumBs_BH,'full');

%switch ON/OFF constraint
yAN = binvar(NumBs,1,'full');
yBH = binvar(NumBs_BH,NumBs_BH,'full');

% Constrainst
C = [];

%% matrice capacità
Ciufinal = zeros(NumBs, NumUEs);

indicebs_in = VirtNode + 1; % inizializzare a 2  
indicebs_end = NumBs + VirtNode; %18 numero delle base station
indiceue_in = indicebs_end + 1; %19 inizializzazione user
indiceue_end = indicebs_end + NumUEs; %31 valore finale degli user
numrig = 1; %inizializzo il numero delle righe
dimANfinal = size(ANfinal);

 for indicebs = indicebs_in : 1 : indicebs_end
   for indexmatrix = 1 : 1 : dimANfinal(1)
     if ANfinal(indexmatrix,1) == indicebs      
        numcol=1; %inizializzo il primo valore delle colonne
        for indnumue = indiceue_in:1:indiceue_end %per far scorrere i valori di A2
          if ANfinal(indexmatrix,2) == indnumue %confronto i valori dentro la matrice con i valori che puo assumere gli user 
            Ciufinal(numrig,numcol) = ANfinal(indexmatrix,3);  %(eq 2)
          end
            numcol = numcol + 1; %scorro le colonne per salvare i risultati           
        end
     end
   end
    numrig = numrig + 1; %scorro le righe per salvare i risultati
 end
 
 % capacity constraint
 product = xu .* Ciufinal; 
 toconstr = sum(product');
 C = [toconstr(:) <= Cimax];

%% power constraint
pmaxANCimax = zeros(NumBs, 1);
numrig = 1; %inizializzo il numero delle righe

for indicebs = indicebs_in : 1 : indicebs_end
  for i = 1 : 1 : dimANfinal(1)
   if ANfinal(i,1) == indicebs
          pmaxANCimax(numrig,1) = ANfinal(i,4);
   end
  end
  numrig = numrig + 1; %scorro le righe per salvare i risultati
end

PAN_out = pmaxANCimax .* toconstr';

PAN_macro = NTXan * (sAN * p0ANMACRO + deltapMACRO * PAN_out);
PAN_small = NTXan * (sAN * p0ANSMALL + deltapSMALL * PAN_out);
PAN = PAN_macro + PAN_small;
%PAN = PAN_small;
PAN=[0; PAN];
 
%% potenza BH
Du = Dfinal(:,3);
BW = BHfinal(:,5);
BW_final = zeros(NumBs_BH,NumBs_BH);
numrig = 1; %inizializzo il numero delle righe
dimBHfinal=size(BHfinal);

for indicebs = 1 : 1 : indicebs_end
  for indexmatrix = 1 : 1 : dimBHfinal(1)
   if BHfinal(indexmatrix,1) == indicebs
    numcol = 1; %inizializzo il primo valore delle colonne
    for numbs = 1 : 1 : indicebs_end %per far scorrere i valori di A2
          if BHfinal(indexmatrix,2) == numbs %confronto i valori dentro la matrice con i valori che puo assumere gli user  
            BW_final(numrig,numcol) = BHfinal(indexmatrix,5);
          end
       numcol = numcol + 1; %scorro le colonne per salvare i risultati           
    end
   end
  end
  numrig = numrig + 1; %scorro le righe per salvare i risultati
end 

load_partial = sdpvar(NumBs_BH,NumBs_BH,NumUEs);

for i = 1 : 1 : NumUEs
    load_partial(:,:, i) = xv(:,:,i) * Du(i);
end
loadBH = sum(load_partial,3) ./ BW_final;

dimBw_BH = size(BW_final);
for i = 1 : 1 : dimBw_BH(1)
    for j = 1 : 1 : dimBw_BH(2)
        if BW_final(i,j) <= 0
            loadBH(i,j) = 0;
        end
    end
end
    
alpha = BHfinal(:,4);
numrig = 1; %inizializzo il numero delle righe
alpha_final = zeros(NumBs_BH,NumBs_BH);

for indicebs = 1 : 1 : indicebs_end
  for indexmatrix = 1 : 1 : dimBHfinal(1)
   if BHfinal(indexmatrix,1) == indicebs
    numcol = 1; %inizializzo il primo valore delle colonne
    for INDnumbs = 1 : 1 : indicebs_end %per far scorrere i valori di A2
          if BHfinal(indexmatrix,2) == INDnumbs %confronto i valori dentro la matrice con i valori che puo assumere gli user  
            alpha_final(numrig,numcol) = BHfinal(indexmatrix,4);
          end
       numcol = numcol + 1; %scorro le colonne per salvare i risultati           
    end
   end
  end
  numrig = numrig + 1; %scorro le righe per salvare i risultati
end

pBH_out = (2 .^ loadBH - 1) .* alpha_final;
%pBH_out(1,:)=0;

C = [C; 0 <= pBH_out(:) <= pmaxBH];

PBH = sum(NTXbh * (sBH * p0BH + deltapBH * pBH_out),2);
%PBH(1) = 0;

%% on/off constraints
% access network
cAN1 = (sum(xu, 2) + M1 * yAN);
cAN2 = sAN + M1 * yAN;
cAN3 = sum(xu, 2) - M2 * (1-yAN);
C = [C; cAN1(:) >= 1; cAN2(:) >= 1; cAN3(:) <= 0];

% backhaul network
cBH1 = (sum(xv, 3) + M3 * yBH);
cBH2 = sBH + M3 * yBH;
cBH3 = sum(xv, 3) - M4 * (1-yBH);
C = [C; cBH1(:) >= 1; cBH2(:) >= 1; cBH3(:) <= 0];

%% path conservation

% matrice access azzera link non esistenti
conlink = zeros(NumBs,NumUEs);
 numrig = 1;
for indicebs = indicebs_in : 1 : indicebs_end
   for indexmatrix = 1 : 1 : dimANfinal(1)
     if ANfinal(indexmatrix,1)  == indicebs      
        numcol = 1; %inizializzo il primo valore delle colonne
        for indnumue = indiceue_in : 1 : indiceue_end %per far scorrere i valori di A2
          if ANfinal(indexmatrix,2) == indnumue %confronto i valori dentro la matrice con i valori che puo assumere gli user 
             conlink(numrig,numcol) = 1;          
          end
            numcol = numcol + 1; %scorro le colonne per salvare i risultati           
        end
     end
   end
    numrig = numrig + 1; %scorro le righe per salvare i risultati
end

for i = 1 : 1 : NumBs
     for j = 1 : 1 : NumUEs
         if conlink(i,j) == 0
             xu(i,j) = 0;
         end
     end
end
 

% azzera link BH non esistenti
sizeBHfinal = size(BHfinal);
conlink_BH = zeros(NumBs_BH, NumBs_BH);

for i= 1 : 1 : sizeBHfinal(1)
    k = BHfinal(i,1);
    j = BHfinal(i,2);
    conlink_BH(k,j) = 1;
end

for i = 1 : 1 : NumBs_BH
    for j = 1 : 1 : NumBs_BH
        if conlink_BH(i,j) == 0
            xv(i,j,:) = 0;
        end
    end
end

xu_single = zeros(1,NumUEs);

xu_0 = [xu_single; xu];


% eq 18
%rowsXu_0 = sum(xu_0, 2);
for u = 1 : 1 : NumUEs
    rowsXv = sum(xv(:,:, u), 2);
    for i = 1 : 1 : NumBs_BH
        C = [C; rowsXv(i) + xu_0(i,u) <= 1];
    end
end


% p2 = sdpvar(NumBs_BH, NumUEs);
% somma = sum(xv,2);
% for i=1:1:NumBs_BH
%     for j=1:1:NumUEs
%         p2(i,j) = somma(1, i,j) + xu_0(i,j);
%     end
% end

% eq 19
p3 = sum(xu, 1);

% eq 17

for u = 1 : 1 : NumUEs
    
    rowsXv = sum(xv(:,:, u), 2);
    columnsXv = sum(xv(:,:, u), 1)';
    
    C = [C; (rowsXv(1)-columnsXv(1)) == 1];
    
    for i = 2 : 1 : NumBs_BH
       C = [C; (rowsXv(i)-columnsXv(i)+xu_0(i,u)) == 0]; 
    end
end

% C = [C; p2(:) <= 1; p3(:) == 1;];
C = [C; p3(:) == 1];

% funzione obiettivo
Objective = sum(PAN + PBH);

ops = sdpsettings('solver','bmibnb','bmibnb.upper','fmincon');

solution = optimize(C,Objective,ops)
