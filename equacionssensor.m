function [output]= equacionssensor (glucose,n,total_time,IncT)
    %% 22/9/16 -FUNCIONA BÉ!
    % change the letters of A B C for mc ma mb and add the fabrication of
    % protein proinsulin and protein insulin
    % finally also model the normal transcription and translation of insulin
    % from and glucose inducible promoter.
    Ga=zeros(1,total_time);
    Gb=zeros(1,total_time);
    Gc=zeros(1,total_time);
    ma=zeros(1,total_time);
    mb=zeros(1,total_time);
    mc=zeros(1,total_time);
    mab=zeros(1,total_time);
    mca=zeros(1,total_time);


    Ga(1)=10;
    Gb(1)=10;
    Gc(1)=10;
    ma(1)=0;
    mb(1)=0;
    mc(1)=0;
    mab(1)=0;
    mca(1)=0;


    Ka(n) =(900 + (7000 / (1 + 200 * (glucose)^ 2))) * 0.001;
    Kb= 3.5;
    delta_a = 0.01;
    delta_b = 0.01;
    K1 = 0.14;
    K_1 = 0.1;
    delta_ab = 0.01;
    Kc = 0.5;
    delta_c = 0.001;
    K2 = 0.14;
    K_2 = 0.1;
    delta_ca = 0.001;
    %K3=0.5;%100; % de mc a ppi proteina.pro-insulina
    delta_ppi=0.0000001;% degradacio ppi proteina.pro-insulina
    K4=0.5; %de ppi a pi; proteina.insulina
    delta_pi =0.0000001; % degracaio pi proteina.insulin

    % euler integration

    ma(n) = ma(n-1) +(Ka(n)* Ga(n-1) -delta_a * ma(n-1) - K1 * ma(n-1) * mb(n-1) + K_1 * mab(n-1) - K2 * mc(n-1) * ma(n-1) + K_2 * mca(n-1)) * IncT;
    mb(n)= mb(n-1) +(Kb * Gb(n-1) - delta_b * mb(n-1) - K1 * ma(n-1) * mb(n-1) + K_1 * mab(n-1))*IncT;
    mab(n)= mab(n-1) + (K1 * ma(n-1) * mb(n-1) - K_1 * mab(n-1) - delta_ab * mab(n-1))*IncT;
    mc(n)= mc(n-1) +(Kc * Gc(n-1) - delta_c * mc(n-1) - K2 * ma(n-1) * mc(n-1) + K_2 * mca(n-1)  )*IncT;
    mca(n) = mca(n-1) + (K2 * ma(n-1) * mc(n-1) - K_2 * mca(n-1) - delta_ca * mca(n-1))*IncT;

    output=mc(n);
end
