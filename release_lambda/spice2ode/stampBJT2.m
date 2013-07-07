function core=stampBJT(core,t)
    core.isNonlinear=1;

    Bf=100;% forward beta
    Br=4; % reverse beta
    Va=10;	    % Early voltage
    Is = 1e-15;% saturation current
    tau=2*10e-11;
    vt = 0.02585126075 ; % 1/40  % thermal voltage
    Vtf= 1;
    Cj = 10 ^ (-14) ;
    Vj=0.8 ;
    fc=0.5;
    mj = 0.5 ;
    beta=100;
    fgain=beta/(beta+1);
    gmin=0;

    
    
    
    
    
%     Cje0=2.6734e-11;
%     Vje=0.67373;
%     Mje= 0.32339;
%     Mjs=0;
%     Cjbe = 10e-11;
%     Cjbc = Cjbe;
%     Cjcs = 2*10e-14;
%     Cjs0 = 0 ;
%     Vjs= 0.50000;
%     Fc=0.95000;
%     Cjc0=Cj;
%     Xcjc= 0.94518 ;
%     Vjc= 0.47372 ;
%     Mjc= 0.34700 ;
%     Ikf=0.402925734377646;
%     Ikr=0;
%     Ise = 8.01e-15;
%     Tr = 380.00e-9;
%     Tf= 9.693e-10;
%     Vaf=100;
%     Var=100;
%     Nf	= 1.00531 ;
%     Ne  = 1.31919 ;
%     Nr  = 1.00688 ; 
%     Nc  = 1.10024 ;
%     Isc = 5E-9; 
%     
%     
%     Xtf=1;
%     Itf= 0.50000;

    Is=6.734e-15 ;
    Xti=3 ;
    Eg=1.11 ;
    Vaf=74.03 ; 
    Bf=416.4 ;
    Ne=1.259 ;
    Ise=6.734e-15; 
    Ikf=66.78e-3 ;
    Xtb=1.5 ;
    Br=0.7371 ; 
    Nc=2 ;
    Isc=0 ; 
    Ikr=0 ;
    Rc=1 ;
    Cjc=3.638e-12;
    Mjc=0.3085; 
    Vjc=0.75; 
    Fc=0.5 ;
    Cje=4.493e-13;
    Mje=0.2593 ; 
    Vje=0.75 ;
    Tr=239.5e-9;
    Tf=301.2e-12;
    Itf=0.4 ;
    Vtf=4 ;
    Xtf=2 ;
    Rb=10 ;
    Nf=1 ;
    Nr=1 ;
    Var=0 ;
    Rbm=0 ;
    Irb=0 ;
    Xcjc=1 ;
    Cjs=0 ;
    Vjs=0.75 ; 
    Mjs=0 ;


    vdcoef = 1/vt;
    vcrit = vt*log( vt/(sqrt(2)*Is) );
    rgain = .5;

    c=core.Element(t).Node1 ;
    b=core.Element(t).Node2 ;
    e=core.Element(t).Node3 ;
    if (strncmpi(core.Element(t).Type, 'npn',4)==1) polarity=1 ;
    else polarity=-1; end
       
    if(c ~= 0) vc = core.vNodeset( c  );
    else vc = 0; end
    if(b ~= 0) vb = core.vNodeset( b  );
    else vb = 0; end
    if(e ~= 0) ve = core.vNodeset( e  );
    else ve = 0; end
    
    disp(sprintf('TTTttTTTTTTTT Transistor %d\n', t));
    vbe=vb-ve
    vce=vc-ve
    vbc=vb-vc
    
    
    
    
    
    
%   %interpret zero as infinity for these model parameters
% 	if(Ikf>0), Ikf=1/Ikf; else Ikf=0; end
% 	if(Ikr>0), Ikr=1/Ikr; else Ikr=0; end
%     if(Vaf>0), Vaf=1/Vaf; else Vaf=0; end
%     if(Var>0), Var=1/Var; else Var=0; end
% 
% 
%  
%   
%   
%     %base-emitter diodes
%     if (vbe < - 10 * vt * Nf) %gtiny = Ube < - 10 * Ut * Nf ? (Is + Ise) : 0;
%         gtiny = (Is + Ise) ;
%     else
%         gtiny = 0 ;
%     end
%     
%     If = pnCurrent (vbe, Is, vt * Nf);
%     Ibei = If / Bf;
%     gif = pnConductance (vbe, Is, vt * Nf);
%     gbei = gif / Bf;
%     Iben = pnCurrent (vbe, Ise, vt * Ne);
%     gben = pnConductance (vbe, Ise, vt * Ne);
%     Ibe = Ibei + Iben + gtiny * vbe;
%     gbe = gbei + gben + gtiny;
% 
%     %base-collector diodes
%     if (vbc < - 10 * vt * Nr ) %gtiny = vbc < - 10 * vt * Nr ? (Is + Isc) : 0;
%         gtiny = (Is + Isc) ;
%     else
%         gtiny = 0 ;
%     end
%     
%     Ir = pnCurrent (vbc, Is, vt * Nr);
%     Ibci = Ir / Br;
%     gir = pnConductance (vbc, Is, vt * Nr);
%     gbci = gir / Br;
%     Ibcn = pnCurrent (vbc, Isc, vt * Nc);
%     gbcn = pnConductance (vbc, Isc, vt * Nc);
%     Ibc = Ibci + Ibcn + gtiny * vbc;
%     gbc = gbci + gbcn + gtiny;
% 
%     %compute base charge quantities
%     Q1 = 1 / (1 - vbc * Vaf - vbe * Var);
%     Q2 = If * Ikf + Ir * Ikr;
%     SArg = 1.0 + 4.0 * Q2;
%     if( SArg>0 )
%         Sq=sqrt (SArg);
%     else
%         Sq= 1;
%     end
%     
%     Qb = Q1 * (1 + Sq) / 2;
%     dQbdUbe = Q1 * (Qb * Var + gif * Ikf / Sq);
%     dQbdUbc = Q1 * (Qb * Vaf + gir * Ikr / Sq);
% 
%     %compute transfer current
%     It = (If - Ir) / Qb;
% 
%     %compute forward and backward transconductance
%     gmf = (gif - If * dQbdUbe / Qb) / Qb;
%     gmr = (gir - Ir * dQbdUbc / Qb) / Qb;
% 
%     %compute old SPICE values
%     go = (gir + It * dQbdUbc) / Qb;
%     gm = (gif - It * dQbdUbe) / Qb - go;
% 
%    
%    disp(sprintf('gm=%f, go=%f, gmf=%f, gmr=%f, Qb=%f, dQbdUbe=%f \n',  gm, go, gmf, gmr, Qb, dQbdUbe));
%     
    
  
    
    if(core.mode==1) %AC small-signal model
        w=core.freq;    
        %I need Cbe, Cbci, CCs, gmf, gmr, 
        
%         
%         %computing charges & capacitances
%         % depletion capacitance of base-emitter diode
%         Cbe = pnCapacitance (vbe, Cje0, Vje, Mje, Fc);
%         Qbe = pnCharge (vbe, Cje0, Vje, Mje, Fc);
% 
%         %diffusion capacitance of base-emitter diode
%   
%         er = 2 * exp (vbc * Vtf);
%         Tff = Tf * (1 + Xtf * sqrt (If / (If + Itf)) * er);
%         dTffdUbe = Tf * Xtf * 2 * gif * If * Itf / (If + Itf)^3 * er;
%         Cbe = Cbe + (If * dTffdUbe + Tff * (gif - If / Qb * dQbdUbe)) / Qb;
%         Qbe = Qbe + If * Tff / Qb;
% 
%         %depletion and diffusion capacitance of base-collector diode
%         Cbc = pnCapacitance (vbc, Cjc0, Vjc, Mjc, Fc);
%         Cbci = Cbc * Xcjc + Tr * gir;
%         Qbc = pnCharge (vbc, Cjc0, Vjc, Mjc, Fc);
%         Qbci = Xcjc * Qbc + Tr * Ir;
%         Cbcx = Cbc * (1 - Xcjc);
%         Qbcx = Qbc * (1 - Xcjc);
% 
%         %depletion capacitance of collector-substrate diode
%         Ccs = pnCapacitance2 (vc, Cjs0, Vjs, Mjs);
%         Qcs = pnCharge2 (vc, Cjs0, Vjs, Mjs);
%         
%         
%         gbc=0;
%         Cbci=0;
%         %these are ac-specific
%         Ybe = gbe + 1i* 2.0 * pi * w * Cbe ;
%         Ybc = gbc + 1i* 2.0 * pi * w * Cbci;
%         Ycs =       1i* 2.0 * pi * w * Ccs;
%         
%         disp(sprintf('frequency=%f', w));
%         disp(sprintf('gbe=%f, gbc=%f, Cbe=%f, Cbci=%f, Ccs=%f, gmf=%f \n',  gbe, gbc, Cbe, Cbci, Ccs, gmf));
% %        pause
%         core = stampY(core, b, b, Ybc + Ybe);
%         core = stampY(core, b, c, -Ybc);
%         core = stampY(core, b, e, -Ybe);
%         core = stampY(core, c, b, -Ybc + gmf - gmr);
%         core = stampY(core, c, c, Ybc + gmr + Ycs);
%         core = stampY(core, c, e, -gmf);
%         core = stampY(core, e, b, -Ybe - gmf + gmr);
%         core = stampY(core, e, c, -gmr);
%         core = stampY(core, e, e, Ybe + gmf);
%         
%         %Ibf = (Is/bF)*(exp(vbe/vt)-1);   %Eq. 13
%         %Ibr = (Is/bR)*(exp(vbc/vt)-1);   %Eq. 14
%         %Icc = bF*Ibf-bR*Ibr; %*(1+vce/va)
%         %gPiF = Ibf/vt ;                  %Eq. 15, g_ce
%         %gPiR = Ibr/vt ;                  %Eq. 16, g_ec
%         %gmf  = bF*gPiF ;                  % g_ee
%         %gmR  = bR*gPiR ;                  % g_cc
%         %go  = Icc/Va;
%         
%         %rPi = bF*vt/Icc ;
%         %ro  = abs(Va)/Icc ;
% 
%         %disp('rPi')
%         %rPi
%         %disp('gmf')
%         %gmf
%         %pause
%         %core = stampConductance(core, b, e, 1/rPi);
%         %core = stampConductance(core, c, e, 1/ro);
%         %core = stampConductance(core, b, e, 1/rPi);
%         %core = StampVCCS(core, gmf, b, e, c, e);
%         
%         %Cbe=1e-9;
%         %Cbci=1e-9;
%         %Ccs=1e-9;
%         %compute admittance matrix entries
%         %Ybe = 1i * 2.0 * pi * w * Cbe ; %+ gbe ;
%         %Ybc = 1i * 2.0 * pi * w * Cbci;%+gbc;
%         %Ycs = 0.0 + 1i * 2.0 * pi * w * Ccs;    
%         %if(g>1e12) g=1e12;end
%         %core=stampConductance(core, b, e, Ybe); 
%         %core=stampConductance(core, b, c, Ybc); 
%         %core=stampConductance(core, c, 0, Ycs); 
%         
    else %DC and Transient Analysis, large-signal model

         disp('-------------------');
        %power-consumption: (v0-v2)*ib + (v1-v2)*ic
        gmin = 0;
        if (core.dcIteration > 10),
            % if we have trouble converging, put a conductance in parallel with all P-N junctions.
            % Gradually increase the conductance value for each iteration.
            gmin = exp(-9*log(10)*(1-core.dcIteration/3000.));
            if (gmin > .1)
                gmin = .1;
            end
        end


        %vbc = pnp*limitStep(pnp*vbc, pnp*lastvbc);
        %vbe = pnp*limitStep(pnp*vbe, pnp*lastvbe);
        %lastvbc = vbc;
        %lastvbe = vbe;



        pcoef = vdcoef*polarity;
        expbc = exp(vbc*pcoef);

        expbe = exp(vbe*pcoef);
        if (expbe < 1)
            expbe = 1;
        end

        ie = polarity*Is*(-(expbe-1)+rgain*(expbc-1));
        ic = polarity*Is*(fgain*(expbe-1)-(expbc-1));
        ib = -(ie+ic);

        gee = -Is*vdcoef*expbe;
       
        gec = rgain*Is*vdcoef*expbc;
        gce = -gee*fgain;
        gcc = -gec*(1/rgain);

        disp(sprintf('voltages: %f %f %f',  vb, vc, ve));
        disp(sprintf('param: beta=%f vcrit=%f fgain=%f', beta, vcrit, fgain ));
        disp(sprintf('T: %f %f ',  vbc, vbe));
        disp(sprintf('gain %f ',  ic/ib));
        disp(sprintf('T %f %f %f %f ', vbc, vbe, ie, ic));

        disp(sprintf('gee= %f ',  gee));
        disp(sprintf('gec= %f ',  gec));
        disp(sprintf('gce= %f ',  gce));
        disp(sprintf('gcc= %f ',  gcc));
        disp(sprintf('gce+gcc= %f ',  gce+gcc));
        disp(sprintf('gee+gec= %f ',  gee+gec));
        
        %pause

        % stamps from page 302 of Pillage.  Node 0 is the base,
        % node 1 the collector, node 2 the emitter.  Also stamp
        % minimum conductance (gmin) between b,e and b,c
        core = stampY(core, b, b, -gee-gec-gce-gcc + gmin*2);
        core = stampY(core, b, c, gec+gcc - gmin);
        core = stampY(core, b, e, gee+gce - gmin);
        core = stampY(core, c, b, gce+gcc - gmin);
        core = stampY(core, c, c, -gcc + gmin);
        core = stampY(core, c, e, -gce);
        core = stampY(core, e, b, gee+gec - gmin);
        core = stampY(core, e, c, -gec);
        core = stampY(core, e, e, -gee + gmin);

        % we are solving for v(k+1), not delta v, so we use formula
        % 10.5.13, multiplying J by v(k)
        core = stampJ(core, b, -ib - (gec+gcc)*vbc - (gee+gce)*vbe);
        core = stampJ(core, c, -ic + gce*vbe + gcc*vbc);
        core = stampJ(core, e, -ie + gee*vbe + gec*vbc);
    end

end