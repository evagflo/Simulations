clear all
close all
clc
%% ALGORTIME EVOLUTIU 26/9 -> seems work good, but re-pasar dema mati bans de seguir!!!
%inicialitzo les variables que es calcularan a dintre del loop, encara que
%la funcio diabetisubject.m ja tingui definides les dimensions segons el
%temps
TotalPatients=2; % 100
numiter=2; %trobar la manera que no sigui jo qui el fixo: depenent de la fitness optima
g1=[];
half_index=zeros(numiter,1);
half_fitness=zeros(numiter,1);

N=zeros(numiter,TotalPatients);
N_ordenat=zeros(numiter,TotalPatients);

fitness_rep=zeros(numiter,TotalPatients);
fitness=zeros(numiter,TotalPatients);
fitness_sort=zeros(numiter,TotalPatients);
fitness_probability=zeros(numiter,TotalPatients);
fitness_cumsum=zeros(numiter,TotalPatients);

number=zeros(numiter,1);
pos=zeros(numiter,1);


%normalsubject.m execute el codi del pacient virtual sa, entra com a valor
%la glucosa basal(definit al paper) i retorna la glucosa al llarg del temps
%i el periode de temps.
[g(1,:)] = normalsubject(91.76);
%per un nombre definit de pacients virtuals, i, calculo la glucosa al llard
%del temps i el temps. despres faig per cada i, pacient virtual diabetic,
%la diferencia quadràtica punt a punt entre el pacient sa i el pacient
%virtual, després en calculo el sumatori.
mu=0.3;
sigma=0.1;


N(1,:)=100*rand(1,TotalPatients);% inicialitzar EL Numero 100 veure q....
for i=2:numiter %i= numero d'iteracions -> ¡CANVIAR PER WHILE!
    ind=1;
    for np=1:TotalPatients %np= pacient virtual
        
        [g1(np,:)] = diabeticimplant(N(i-1,np));
        diferencia(np,:)= (g(1,:) - g1(np,:)).^2;
        fitness(i,np)= sum(diferencia(np,:));
    end   
        
        
        % ordenar les fitness de major a menor : comanda sort(vector)
        % [A,B]= sort (C)
        % A contidrà els valors de C ordenats de gran a petit i B contindrà
        % l'index on prèviament es trobaven els valors ara ordenats
        [fitness_sort(i,:), index_fitness(i,:)]= sort(fitness(i,:),'descend');  
        
        % reordenar el vector que guarda les característiques de cada pacient segons l'ordre de la fitness.
        N_ordenat(i,:) = N(i-1, index_fitness(i,:));
        half_fitness(i)= median(fitness_sort(i,:));%mean(fitness_sort(i,:));
        %%% el problema esta aqui: cm trobar el indec que ocupa el valor
        %%% mitja  de la fitness
        half_index(i)= find(fitness_sort(i,:)>=half_fitness(i),1,'last');      
        
        %replico ara només de la posició 1 a la posició half_index 
        N(i,half_index(i):end)=N(i-1,half_index(i):end);
        fitness_rep(i,half_index(i):end)=fitness(i,half_index(i):end);
        
        % escalar fitness 1-> pitjor fitnees 0-> millor fitness divideixo cada valor de fitness de cada pacient pel total  
        suma(i)=sum(fitness_rep(i,half_index(i):end));
        fitness_probability(i,:)=fitness_rep(i,:)./suma(i);
        %replico els llindars anteriors i introdueixo mutacions segons una normal distirbution
        %l'altre meitat: de half_index+1:end ha de ser proporcional a la
        %resta de individus que ja tinc -> mirarho
        fitness_cumsum(i,:)= cumsum(fitness_probability(i,:));
        
        
    % half index-> aqui hi ha un problema
        while (ind <= half_index(i))
            number(i,ind)= rand(1);% 
            pos(i,ind)= find(number(i,ind)<fitness_cumsum(i,:),1,'first')-1;
            nova_pos(ind)= pos(i,ind); 
            N(i,ind)= N(i,nova_pos(ind)) + normrnd(mu,sigma);
            ind= ind + 1;
        end
        
end

