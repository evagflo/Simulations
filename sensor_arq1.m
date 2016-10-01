function [mc]= equacionssensor (glucose)
%% 22/9/16 -FUNCIONA BÉ!
% change the letters of A B C for mc ma mb and add the fabrication of
% protein proinsulin and protein insulin
% finally also model the normal transcription and translation of insulin
% from and glucose inducible promoter.


%% glucose vector definition: each position corresponds to a different glucose concentration
L_glc = 1;

%% variables definition
Ga = 10*ones(L_glc,1);
Gb = 10*ones(L_glc,1);
Gc =10*ones(L_glc,1);
ma =zeros(L_glc,1);
ma_New=zeros(L_glc,1);
mb = zeros(L_glc,1);
mb_New=zeros(L_glc,1);
mab = zeros(L_glc,1);
mab_New=zeros(L_glc,1);
mc = zeros(L_glc,1);
mc_New=zeros(L_glc,1);
mca = zeros(L_glc,1);
mca_New=zeros(L_glc,1);

ppi=zeros(L_glc,1);
ppi_New=zeros(L_glc,1);
pi=zeros(L_glc,1);
pi_New=zeros(L_glc,1);



Incma=10*ones(L_glc,1);
Incmb=10*ones(L_glc,1);
Incmc=10*ones(L_glc,1);
temps=0;
cont=0;
hill=[];
IncT = 0.00001;
    for m= 1:L_glc  
        cont =0;
        while (cont <= 500000)
            % parameters definition
            f_G =(900 + (7000 / (1 + 200 * (glucose(m))^ 2))) * 0.001;
            Ka(m) = f_G;
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
            K3=0.5;%100; % de mc a ppi proteina.pro-insulina
            delta_ppi=0.0000001;% degradacio ppi proteina.pro-insulina
            K4=0.5; %de ppi a pi; proteina.insulina
            delta_pi =0.0000001; % degracaio pi proteina.insulin
            

            % euler integration
            ma_New(m) = ma(m) +(Ka(m) * Ga(m) -delta_a * ma(m) - K1 * ma(m) * mb(m) + K_1 * mab(m) - K2 * mc(m) * ma(m) + K_2 * mca(m)) * IncT;
            mb_New(m)= mb(m) +(Kb * Gb(m) - delta_b * mb(m) - K1 * ma(m) * mb(m) + K_1 * mab(m))*IncT;
            mab_New (m)= mab(m) + (K1 * ma(m) * mb(m) - K_1 * mab(m) - delta_ab * mab(m))*IncT;
            mc_New(m)= mc(m) +(Kc * Gc(m) - delta_c * mc(m) - K2 * ma(m) * mc(m) + K_2 * mca(m) - K3 *mc(m) + K3 * mc(m) )*IncT;
            mca_New (m) = mca(m) + (K2 * ma(m) * mc(m) - K_2 * mca(m) - delta_ca * mca(m))*IncT;
            ppi_New (m)= ppi(m) +( K3 * mc(m) - delta_ppi * ppi(m) - K4 * ppi(m))*IncT;
            pi_New (m)= pi(m) +(K4 *ppi(m) - delta_pi * pi(m))*IncT;
            
            % difference between steps calculation -> necessary to find the
            % steady state
            Incmc(m) = abs(mc(m) - mc_New(m));
            Incma(m) = abs(ma(m) - ma_New(m));
            Incmb(m) = abs(mb(m) - mb_New(m));
            Incppi(m) = abs(ppi(m) - ppi_New(m));
            Incpi(m) = abs(pi(m) - pi_New(m));
            % variables actualization
            ma(m)= ma_New(m);
            mb(m)=mb_New(m);
            mc(m)=mc_New(m);
            mab(m)=mab_New(m);
            mca(m)=mca_New(m);
            ppi(m)=ppi_New(m);
            pi(m)=pi_New(m);
            
            %time calculation
            temps = temps + IncT;
            cont = cont + 1;

        end
%     
    end
% figure(1)
% subplot(2,1,1)
% plot(glucose_levels,ma)
% xlabel('glucosa')
% ylabel('A mRNA')
% subplot(2,2,2)
% plot(glucose_levels,mb)
% xlabel('glucosa')
% ylabel('B mRNA')
% subplot(2,2,3)
% plot(glucose_levels,mc)%
% xlabel('glucosa')
% ylabel('C mRNA')
% subplot(2,2,4)
% plot(glucose_levels,Ka)%
% xlabel('glucose concentration(%)')
% ylabel('Ka')
% 
%    
% figure(2)
% subplot(2,2,1)
% plot(glucose_levels,mc)
% xlabel('glucosa')
% ylabel('mRNA')
% subplot(2,2,2)
% plot(glucose_levels,ppi)
% xlabel('glucosa')
% ylabel('pro-insulin protein')
% subplot(2,2,3)
% plot(glucose_levels,pi)%
% xlabel('glucosa')
% ylabel('insulin protein')
% subplot(2,2,4)
% plot(glucose_levels,Ka)%
% xlabel('glucosa')
% ylabel('Ka')
% %
% 
% figure(3)
% subplot(1,2,1)
% plot(glucose_levels,Ka)%
% xlabel('glucose concentration(%)')
% ylabel('Ka')
% title('Natural biosensor')
% subplot(1,2,2)
% plot(glucose_levels,mc)%
% xlabel('glucose concentration(%)')
% ylabel('C mRNA')
% title('Synthetic biosensor')

figure()
plot(glucose, hill,'r')
hold on
plot(glucose,mc,'b')

