clc; close all; clear;
d1 = fopen('source_JPARK2019_vfinal.txt', 'r');
f1 = '%1s';

%-------------------------------------
%  총 Source와 Count 구하기
%-------------------------------------

TEST_H = [];
reading = fscanf(d1,f1);
length_ALL = length(reading);
Source_ALL = [];
Count_ALL = [];
for i = 1:length_ALL
    length2 = length(Source_ALL);
    same = 0;
    for j = 1:length2
        if reading(i) == Source_ALL(j)
           same = j;
        end
     end

        if same > 0
           Count_ALL(same) = Count_ALL(same) + 1;
        else
           Source_ALL = [Source_ALL, reading(i)];
           Count_ALL = [Count_ALL, 1];
        end
end

%------------------------------------------
% t이전의 source에 부여할수있는 코드길이 계산
%------------------------------------------

R = 1;
while 2^R < length(Source_ALL)
    R = R + 1;
end

%----------------------------------
% Learning 시작
%----------------------------------
for x = R : 7
    LearningPhase_array = length_ALL / 1000 :length_ALL / 1000 : 4*length_ALL/50;
    for idx_learning = 1:length(LearningPhase_array)
        length1 = LearningPhase_array(idx_learning);
        Source = [];
        Count = [];
            for i = 1:length1
                length2 = length(Source);
                same = 0;
                for j = 1:length2
                    if reading(i) == Source(j)
                        same = j;
                    end
                end
                if same > 0
                   Count(same) = Count(same) + 1;
                else
                   Source = [Source, reading(i)];
                   Count = [Count, 1];
                end
            end
        Prob = Count/length1;
        Count_f = Count;

        %-------------------------------------
        %  Huffman encoding에서 code length 구하기
        %-------------------------------------

        Prob_ex = Prob;
        [SS] = C_L(Source, Prob_ex);
        SS_s = sort(SS, 'descend');
        SS_f = SS;

        %-------------------------------------
        %  code length이용해서 최적화
        %-------------------------------------

        DIM = length(Count_ALL) - length(SS_f);
        if DIM >0
            for i = 1 : DIM
                Count_f = [Count_f, 0];
                SS_f = [SS_f, x];
            end
        end
        Count_H = Count_ALL - Count_f;
        HUFF = 0;
        for i = 1:length(Count_H)
            HUFF = HUFF + (Count_H(i) * SS_f(i));
        end
        TEST = x*length1 + HUFF;
        TEST_H = [TEST_H, TEST];
    end
    Avg_H = TEST_H / length_ALL;
    plot(LearningPhase_array, Avg_H((x-R)*80+1:(x-R+1)*80))
    title('Plots')
    hold on
end

TT = [];
Good = [];
Source_code_t = [];
Length_code_t = [];

for x = R:7
    y = (x-R)*80+1;
    [good, index] = min(TEST_H(y:y+79));
    t = index*length_ALL/100;
    TT = [TT, t];
    Good = [Good, good];
    
    %--------------------------------------------
    %  t에 관해서 코드길이 다시 계산
    %--------------------------------------------
    
    Source_t = [];
    Count_t = [];   
    for idx = 1:t
        if not(ismember(reading(idx), Source_t))
           Source_t = [Source_t, reading(idx)];
        end

    end
    for idx = 1:length(Source_t)
        Count_t(idx) = count(reading(1:t), (Source_t(idx)));
    end

    Prob_t = Count_t / TT(x-R+1);
    [Source_code, Prob_code] = S(Source_t, Prob_t);
    Length = C_L(Source_code, Prob_code);
    Count_code = Prob_code * TT(x-R+1);
    length_sc = length(Source_code);
    DIM_t = length(Count_ALL) - length_sc;
    
    if DIM_t >0
       for i = 1 : DIM_t
           Count_code = [Count_code, 0];
           Length = [Length, x];
           Source_code = [Source_code, Source_ALL(length_sc+i)];
       end
    end
    Source_code_t = [Source_code_t, Source_code];
    Length_code_t = [Length_code_t, Length];
end
fclose(d1);

%----------------------------------------------------
% t이전코드의 길이에 관해 변하는 t와 각각의 허프만코드길이 출력
%----------------------------------------------------

HUFFMAN = {};
length3 = length(Source_ALL);
for i = 1 : 7-R+1
    H_L = Length_code_t((i-1)*length3+1 : i*length3);
    for j = 1: length3
        HUFFMAN{j, 2*i} = H_L(j);
        HUFFMAN{j, 2*i-1} = Source_code_t(length3*(i-1)+j);
    end
end

for i = R:7
    txt5 = [];
    txt1 = ['ASCII코드 길이가', string(i), '일때, t는 ', TT(i-R+1)];
    for j = 1:7-R+1
        txt4 = [];
        for k = 1:length(Source_ALL)
            txt2 = HUFFMAN{k, 2*j-1};
            txt3 = HUFFMAN{k, 2*j};
            txt4 = [txt4, string(txt2), '의 길이는', string(txt3)];
        end
        txt5 = [txt5, txt4];
    end
    txt5((i-R)*3*length(Source_ALL)+1: (i-R+1)*3*length(Source_ALL));
end

%--------Function : Sorting Source&Prob--------
function [L1, L2] = S(S, P)
    Source_s = [];
    Prob_s = [];
    length_s = length(S);
    for k = 1 : length_s
        [M, I] = min(P);
        Source_s = [Source_s, S(I)];
        Prob_s = [Prob_s, M];
        P(I) = 1;
    end
    L1 = Source_s;
    L2 = Prob_s;
end

%-------Function : Conmbinations--------
function [y1, y2, y3] = C(S, P)
    S_comb = [S(1), S(2)];
    P_comb = P(1) + P(2);
    S_h = S;
    P_h = P;
    S_h(1:2) = [];
    P_h(1:2) = [];
    S_h = [S_h, S(1)];
    P_h = [P_h, P_comb];
    y1 = S_comb;
    y2 = S_h;
    y3 = P_h;
end

%----------Function : Find index-------
function [N] = F(SSS,S)
    N = [];
    SSS1 = unique(SSS);
    for i = 1:length(SSS1)
        for j = 1:length(S)
            if SSS1(i) == S(j)
                N = [N, j];
            end
        end
    end
end

%-------Function : code length------
function [SS] = C_L(Source, Prob)
length2 = length(Source);
SS = [];
    for i = 1:length(Source)
        SS = [SS,0];
    end
    SSS = [];
    [Source_0, Prob_0] = S(Source, Prob);
    for i = 1 : length2 - 1
        [COMB, Source_0, Prob_0] = C(Source_0, Prob_0);
        [Source_0, Prob_0] = S(Source_0, Prob_0);
        SSS = [SSS, COMB];
        N = F(SSS, Source);
        for j = 1:length(N)
            k = N(j);
            SS(k) = SS(k) + 1;
        end
    end    
     L = length(SSS);
     length3 = L / 2;
     for m = 2 : length3
         M = 2*m;
         k = (ismember(SSS(M-1), SSS(1:M-2))) + (ismember(SSS(M), SSS(1:M-2)));
         if k == 0
            Nm = F(SSS(1:M-2),Source);
            for i = 1:length(Nm)
                v = Nm(i);
                SS(v) = SS(v) - 1;
            end
         end
     end
end
