function [g,time] = diabeticsubject(gb)

    % GIM Model Eva
    % bona 11/3/16
    %forward-euler
    clear all
    %close all

    %% parameters: set the parameters of table 1 and steady state condictions
    % given variables
    % given variables
    gb=180; % mg/dl
    ib= 0; % pmol/L
    egpb= 2.40; % mg/kg/min
    g_tar=130; % mg/dL
    % time

    dt=0.01;
    t=0:dt:(24*60);
    total_time= length(t)-1;

    ing_glc= zeros(1,total_time);
    ing_glc(1)=0;
    ing_glc((8*60)/dt:(8*60)/dt+(1/dt))=45000;
    ing_glc((12*60)/dt:(12*60)/dt+(1/dt))=70000;
    ing_glc((20*60)/dt:(20*60)/dt+(1/dt))=70000;
    ing_glc((24*60)/dt)=0;



    % glucose kinetics
    vg = 1.88; %distribution volume of glucose (dl/kg)
    k1= 0.065; % rate parameter gp-> gt (1/min)
    k2=0.079; % rate parameter gt-> gp (1/min)
    fcns=1;% glucose uptake by brain an eithrocytes

    gpb= gb * vg; %glucose mass in plasma (mg/kg)
    gtb=(fcns - egpb + k1*gpb)/k2; %glucose mass in tissue (mg/kg)

    % insulin kinetics
    vl=0.05; %distrubution volume of insulin( l(kg)
    m1=0.190; %rate parameter (1/min)
    ipb=ib*vl; % insulin mass in plasma ( pmol/l)
    heb=0.6; % hepatic extraction of insulin
    m2=0.484;%(( sb/ipb ) - ( m4/ (1-heb))) * (1-heb)/(heb);%0.484; %rate parameter (1/min)
    m4=2/3*m2*heb;% (2/5)* (sb/ipb)*(1-heb);%0.194;%rate parameter (1/min)
    m5=0.0304;%rate parameter (1/min)
    m3b=(heb*m1)/(1-heb);%rate parameter (1/min)
    ilb=ipb*(m2/(m1+m3b)); % insulin mass in liver (pmol/l)
    sb=ipb*m4+ilb*m3b; % insulin secretion (pmol/kg/min)
    m6=m5*sb + heb ;

    %4. insulin secretion
    KK=2.30; %pancreatic responsivity to the gl rate of change (pmol/l per mg/dl))
    alpha=0.050;%delay between glc signal and insulin secretion (1/min)
    beta=0.11; % pancreatic responsivity to glc (pmol/kg/min per mg/dl)
    gamma=0.5; % transfer rate constant between portal vein and liver (1/min)
    db=sb;
    yb=0;
    spob=sb;
    h=gb;% threshold level of glc above b-cell initiate to produce new insulun ( mg/dl) 

    %1. endogenous glucose productions
    kp2=0.0021;% liver glc efectiveness (1/min)
    kp3=0.009; % amplitude of insulin action on the liver ( mg/kg/min per pmol/l)
    kp4=0.0618; %amplitude of portal insulin action on the liver (mg/kg/min /(pmol/kg))
    ki=0.0079; % rate for delay between insulin signal and insullin action (1/min)
    kp1=egpb+kp2*gpb+kp3*ib;%2.70;%extrapolated egp at zero glc and ins (mg/kg/min)%egpb + kp2*gpb + kp3 *ib + kp4 *ipob;%
    i1b=ib;
    idb=ib; %delayed  insulin signal (pmol/kg)

    %2.glucose rate of appearence
    kmax=0.0558;
    kmin=0.0080;
    kabs=0.057; % rate constant of gastric emptying (1/min)
    kgri=0.0558; %rate of grinding(1/min)
    bw=78;% body weight (mg)
    f=0.90; % fracion of intestinal absortion which appears in plasma 
    a=0.00013;
    b=0.82;
    c=0.00236;
    d=0.010;% , dd is the impulse function, D is the amount of ingested glucose,
    %dosekempt= 90000;
    qstob=0; %amount of glucose in the stomach(mg)
    qsto1b=0;%amount of solid-glucose in the stomach(mg)
    qsto2b=0;%amount of liquid-glucose in the stomach(mg)
    qgutb=0; %glucose mass in the intestine (mg)
    rab=0; %glucose rate of appearence in plasma (mg/kg/min)

    %3. glucose utilization
    fcns=1;% glucose uptake by brain an eithrocytes
    vmx=0.047; % m.m vel max
    kmo=255.59; % m.m km
    p2u=0.0331; % insulin action in the periphereal glc utilization(1/min)
    vmo=((egpb-fcns)*(kmo+gtb))/(gtb);%2.50;%
    ub=egpb; % glucose utilization (mg/kg/min)
    xb=0; %insulin in the intersticial fluid (pmol/l)
    uiib=fcns;% insulin-independent glc utilization (mg/kg/min)
    uidb=(vmo * gtb)/(kmo + gtb);% insulin-dependent glc utilization (mg/kg/min)

    %5. glucose renal excretion
    ke1=0.0005; % glomerular filtration rate (1/min)
    ke2=339;  % renal threshold of glucose (mg/kg)
    eb=0; % renal excretion (mg/kg/min)
    %subcutaneous insulin kinetics
    kd=0.0164;
    ka1=0.0018;
    ka2=0.0182;
    ksc=0.03;%
    g_mb=gb;
    Ts=10;%min
    Kp=0.032;
    TD=66;
    Tl=450;
    isc1ss= (ipb/(kd+ka1))*(m2+m4 - (m1*m2)/(m1+m3b));
    isc2ss= (kd/ka2)*isc1ss;
    %% variables: declare all the variables as vectors
    time=zeros(1,total_time);
    %glucose system
    gp=zeros(1,total_time);
    gt=zeros(1,total_time);
    g=zeros(1,total_time);
    g_der=zeros(1,total_time);%%
        % endogenous glucose production
    egp=zeros(1,total_time);
    i1=zeros(1,total_time);
    id=zeros(1,total_time);
    ipo=zeros(1,total_time);
        % glucose rate of appearence
    qsto=zeros(1,total_time);
    qsto1=zeros(1,total_time);
    qsto2=zeros(1,total_time);
    qgut=zeros(1,total_time);
    ra=zeros(1,total_time);
    %ing_glc=zeros(1,total_time);
    kempt=zeros(1,total_time);
    dosekempt=zeros(1,total_time);
        %glucose utilization
    uii=zeros(1,total_time);
    uid=zeros(1,total_time);
    vm=zeros(1,total_time);
    km=zeros(1,total_time);
    x=zeros(1,total_time);
    u=zeros(1,total_time);
        %glucose renal excretion
    e=zeros(1,total_time);
    % insulin system=
    il=zeros(1,total_time);
    ip=zeros(1,total_time);
    i=zeros(1,total_time);
    n=zeros(1,total_time);
    m3=zeros(1,total_time);
    he=zeros(1,total_time);
        % insulin secretion
    s1=zeros(1,total_time);
    s2=zeros(1,total_time);
    s=zeros(1,total_time);
    g_m=zeros(1,total_time);
    y=zeros(1,total_time); 
    g_s=zeros(1,total_time);

    %% initialize the variables using the steady state condition
    %glucose system
    gp(1)=gpb;
    gt(1)=gtb;
    g(1)=gb;
    g_der(1)=0;
        % endogenous glucose production
    egp(1)=egpb;
    i1(1)=ib;
    id(1)=ib;

        % glucose rate of appearence
    qsto(1)=0;
    qsto1(1)=0;%qsto1b;
    qsto2(1)=0;%qsto2b;
    qgut(1)=0;%qgutb;
    ra(1)=0;%rab;

        %glucose utilization
    uii(1)=uiib;
    uid(1)=uidb;
    vm(1)=vmo;
    x(1)=xb;
    u(1)=ub;
        %glucose renal excretion
    e(1)=eb;
    % insulin system=
    il(1)=ilb;
    ip(1)=ipb;
    i(1)=ib;
    m3(1)=m3b;
    he(1)=heb;
      % insulin secretion
    g_m(1)=g_mb;
    g_s(1)=gb;
    isc1(1)=isc1ss;
    isc2(1)=isc2ss;
    %% set the model equations and iterate / integrate (euler)
    cont=1;
    for n= 2:total_time
    %variables

        %endogenous glucose production
        egp(n)= max(0,kp1 - kp2 * gp(n-1) - kp3 * id(n-1) - kp4 * ipo(n-1));
        g(n)= gp(n-1)/vg;
        i(n)=ip(n-1)/vl;
        %rate of glucose appearence - correctt
        qsto(n)= qsto1(n-1) + qsto2(n-1);
        qsto1(n)= qsto1(n-1) + dt * ( ing_glc(n) - kgri* qsto1(n-1));
        if ing_glc(n) ~= 0 && cont ~= 0
            cont= n;
        end
            dosekempt= qsto(cont) ;
            aa(n)=5/2*(1-b)*dosekempt;
            bb(n)=5/2*d*dosekempt;
            xx(n)=aa(n)*(qsto(n) - b * dosekempt);%inf*0 = nan
            yy(n)=bb(n)*(qsto(n) - c * dosekempt); 
            kempt(n)= kmin + ((kmax-kmin)/2) * (tanh(xx(n)) - tanh(yy(n)) + 2);
            qsto2(n)= qsto2(n-1) + dt*(kgri * qsto1(n-1) - kempt(n) * qsto2(n-1));
            qgut(n)=qgut(n-1) + dt*(kempt(n) * qsto2(n-1) - kabs * qgut(n-1));
        ra(n)= max(0,f*kabs*qgut(n-1)/bw);%be





        %glucose utilization
        vm(n)= vmo + vmx * x(n-1);
        uii(n)= fcns;
        uid(n)= (vm(n)*gt(n-1))/(kmo + gt(n-1));
        u(n)= uii(n) + uid(n);
        % insulin secretion
        % glucose renal excretion
        if gp(n-1) > ke2
            e(n) =  ke1 * ( gp(n-1) - ke2);
        else
            e(n)= 0;
        end
        %subcutaneus insulin kinetics

        ri(n)= (ka1 * s1(n-1) + ka2 *s2(n-1));

        pro(n)= Kp *( (g_s(n-1)) - g_tar);
        int(n)= (Kp/Tl) * ( (g_s(n-1)) - g_tar);
        der(n)= Kp*TD * ( -(1/Ts) * g_s(n-1) + (1/Ts)*g(n)); 
        pid(n)= (pro(n) + int(n) + der(n)); 

        he(n)= -m5 * s(n) + m6;
        m3(n)=(he(n)* m1)/(1- he(n));

        %deriv
    %     % glucose utilization
        x(n)= x(n-1) +dt*floor (( p2u * ((ip(n-1)/vl) - ib) - p2u * x(n-1)));

        %glucose endogenous production
        i1(n)= i1(n-1) +dt*floor ( (ki * (ip(n)/vl)- ki * i1(n-1)));
        id(n)=id(n-1) + dt*floor (( ki* i1(n-1) - ki* id(n-1)));

        % subcunatenus insulin secretion
       if pid(n)>0
           s1(n)= s1(n-1) + dt* (pid(n) -(ka1+kd)*s1(n-1));
       else
           s1(n)=s1(n-1) + dt* ( -(ka1+kd)*s1(n-1));
       end
        s2(n)= s2(n-1) + dt*( kd*s1(n-1) - ka2*s2(n-1));

        %glucose subsytem
        gp(n)= gp(n-1) + dt* (egp(n) + ra(n) - uii(n) - e(n) - k1*gp(n-1) + k2*gt(n-1));
        gt(n)= gt(n-1) + dt*( k1*gp(n-1) - k2*gt(n-1) - uid(n));

        g_m(n)= g_m(n-1) + dt*(- ksc * g_m(n-1) + ksc*(gp(n-1)/vg));

        g_s(n)= g_s(n-1) + dt*( (-(1/Ts)*g_s(n-1)+ (1/Ts)*(gp(n-1)/vg)));
        % insulin subsystem
        il(n)= il(n-1) + dt*(- (m1+m3(n)) * il(n-1) + m2*ip(n-1) );
        ip(n)= ip(n-1) + dt*(- (m2+m4)*ip(n-1) + m1*il(n-1)+ ka1*s1(n-1) + ka2*s2(n-1) );

        time(n)= time(n-1)+ dt; 

    end

    %% plots

%     figure()
%     subplot(3,2,1)
%     plot(time, g)
%     hold on
%     plot( [time(1) time(end)],[gb gb])
%     title('plasma glucose')
%     xlabel('time (h)')
%     ylabel('mg/dl')
% 
% 
%     subplot(3,2,2)
%     plot(time, i)
%     hold on
%     plot( [time(1) time(end)],[ib ib])
%     title('plasma insulin ')
%     xlabel('time (h)')
%     ylabel('pmol/L')
% 
%     subplot(3,2,3)
%     plot(time, egp)
%     xlabel('time (h)')
%     title('glucose endogenous production')
%     ylabel('mg/kg*min')
% 
%     subplot(3,2,5)
%     plot(time, u)
%     xlabel('time (h)')
%     title('glucose utilization')
%     ylabel('mg/kg*min')
% 
% 
%     subplot(3,2,4)
%     plot(time, ra)
%     xlabel('time (h)')
%     title('glucose rate of appearence')
%     ylabel('mg/kg*min')
% 
%     subplot(3,2,6)
%     plot(time, ri)
%     xlabel('time (h)')
%     title('insulin rate of appearence')
%     ylabel('pmol/kg/min')

end