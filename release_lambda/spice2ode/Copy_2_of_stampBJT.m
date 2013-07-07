function core=stampBJT(core,t)
    bF=100;% forward beta
    bR=4; % reverse beta
    Va=10;	    % Early voltage
    %tau = 2 * 10^(-11) ;
    leakage = 1e-13;% saturation current
	vt = 0.02585126075 ; % 1/40  % thermal voltage
    Vtf= 1;
    Cj = 10 ^ (-14) ;
    %Vj = 0.8 ;
    %fc = 0.5 ;
    mj = 0.5 ;
    beta=100;
    fgain=beta/(beta+1);
    gmin=0;
    Is = 10^(-15);  % saturation current
% .MODEL t2n2222a NPN
% + ISS= 0. XTF= 1. NS = 1.00000
% + CJS= 0. VJS= 0.50000 PTF= 0.
% + MJS= 0. EG = 1.10000 AF = 1.
% + ITF= 0.50000 VTF= 1.00000
% + BR = 40.00000 IS = 1.6339e-14 VAF= 103.40529
% + VAR= 17.77498 IKF= 1.00000
% + NE = 1.31919 IKR= 1.00000 ISC= 3.6856e-13
% + NC = 1.10024 IRB= 4.3646e-05 NF = 1.00531
% + NR = 1.00688 RBM= 1.0000e-02 RB = 71.82988
% + RC = 0.42753 RE = 3.0503e-03 MJE= 0.32339
% + MJC= 0.34700 VJE= 0.67373 VJC= 0.47372
% + TF = 9.693e-10 TR = 380.00e-9 CJE= 2.6734e-11
% + CJC= 1.4040e-11 FC = 0.95000 XCJC= 0.94518
    Cje0=2.6734e-11;
    Vje=0.67373;
    Mje= 0.32339;
    Mjs=0;
    Cjbe = 10e-11;
    Cjbc = Cjbe;
    Cjcs = 2*10e-14;
    Cjs0 = 0 ;
    Vjs= 0.50000;
    Fc=0.95000;
    Cjc0=Cj;
    Xcjc= 0.94518 ;
    Vjc= 0.47372 ;
    Mjc= 0.34700 ;
    
    Tr = 380.00e-9;
    Tf= 9.693e-10;
    VAf=100;
    VAr=100;
    
    vdcoef = 1/vt;
    rgain = .5;
    vcrit = vt*log( vt/(sqrt(2)*leakage) );

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
    vbe=vb-ve;
    vce=vc-ve;
    vbc=vb-vc;
    
    core.isNonlinear=1;
    
    if(core.mode==1) %AC small-signal model
       
        w=core.freq;       
        
        %interpret zero as infinity for that model parameter
        if( Vtf > 0 )
            Vtf=1.0 / Vtf ;
        else
            Vtf=0;
        end

        %depletion capacitance of base-emitter diode
        Cbe = pnCapacitance (vbe, Cje0, Vje, Mje, Fc);
        Qbe = pnCharge (vbe, Cje0, Vje, Mje, Fc);

 
        %depletion and diffusion capacitance of base-collector diode
        Cbci = pnCapacitance (vbc, Cjc0 * Xcjc, Vjc, Mjc, Fc) ; %+ Tr * gir;
        Qbci = pnCharge (vbc, Cjc0 * Xcjc, Vjc, Mjc, Fc) ;%+ Tr * Ir;

        %depletion and diffusion capacitance of external base-collector capacitor
        Cbcx = pnCapacitance (vce, Cjc0 * (1 - Xcjc), Vjc, Mjc, Fc);
        Qbcx = pnCharge (vce, Cjc0 * (1 - Xcjc), Vjc, Mjc, Fc);

        %depletion capacitance of collector-substrate diode
        Ccs = pnCapacitance2 (vc, Cjs0, Vjs, Mjs);
        Qcs = pnCharge2 (vc, Cjs0, Vjs, Mjs);


        Ibf = (Is/bF)*(exp(vbe/vt)-1);   %Eq. 13
        Ibr = (Is/bR)*(exp(vbc/vt)-1);   %Eq. 14
        Icc = bF*Ibf-bR*Ibr; %*(1+vce/va)
        gPiF = Ibf/vt ;                  %Eq. 15, g_ce
        gPiR = Ibr/vt ;                  %Eq. 16, g_ec
        gmf  = bF*gPiF ;                  % g_ee
        gmR  = bR*gPiR ;                  % g_cc
        go  = Icc/Va;
       
        %gbe  = gbei + gben;
        %gbc  = gbci + gbcn;
        %gmfr = gm;
        

        Cbe=1e-9;
        Cbci=1e-9;
        Ccs=1e-9;
        %compute admittance matrix entries
        Ybe = 1i * 2.0 * pi * w * Cbe ; %+ gbe ;
        Ybc = 1i * 2.0 * pi * w * Cbci;%+gbc;
        Ycs = 0.0 + 1i * 2.0 * pi * w * Ccs;
        disp(sprintf('capacitances= %f %f %f , %f %f %f', Cbe, Cbci, Ccs, Ybe, Ybc, Ycs ));

        
        core = stampY(core, b, b, Ybc + Ybe);
        core = stampY(core, b, c, -Ybc);
        core = stampY(core, b, e, -Ybe);
        core = stampY(core, c, b, -Ybc + gmf);
        core = stampY(core, c, c, Ybc + Ycs + go);
        core = stampY(core, c, e, -gmf - go);
        core = stampY(core, e, b, -Ybe - gmf);
        core = stampY(core, e, c, -go);
        core = stampY(core, e, e, Ybe + gmf + go);

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

        ie = polarity*leakage*(-(expbe-1)+rgain*(expbc-1));
        ic = polarity*leakage*(fgain*(expbe-1)-(expbc-1));
        ib = -(ie+ic);

        gee = -leakage*vdcoef*expbe;
        gec = rgain*leakage*vdcoef*expbc;
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