function core=stampBJT(core,i)

    %bF=100;% forward beta
    %bR=4; % reverse beta
    %Va=10;	    % Early voltage
    %tau = 2 * 10^(-11) ;
    %Is = 10^(-15);  % saturation current
    %Vt = 0.02585126075 ;    % 1/40  % thermal voltage
    %Cj = 10 ^ (-14) ;
    %Vj = 0.8 ;
    %fc = 0.5 ;
    %mj = 0.5 ;

    
	leakage = 1e-15;%1e-13;
	vt = 0.02585126075 ;
	vdcoef = 1/vt;
	rgain = .5;
    vcrit = vt*log( vt/(sqrt(2)*leakage) );
    lastvbc=0;
    lastvbe=0;
    beta=100;
    fgain=beta/(beta+1);
    gmin=0;

    Cjbe = 10e-11;
    Cjbc = Cjbe;
    Cjcs = 2*10e-14;
    fc=0.5;
    mj=0.33;
    VAf=100;
    VAr=100;
    tf=2*10e-11;
    %v0=0;
    %v1=-lastvbe;
    %v2=-lastvbc;

    %power-consumption: (v0-v2)*ib + (v1-v2)*ic
    c=core.Element(i).Node1 ;
    b=core.Element(i).Node2 ;
    e=core.Element(i).Node3 ;
    if (strncmpi(core.Element(i).Type, 'npn',4)==1) polarity=1 ;
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

    if(core.mode==1)
        w=core.freq;
        Cjbe=1e-12;
        c1= j*2*pi*w*Cjbe ;
        c2= j*2*pi*w*Cjbe ;
        c3= j*2*pi*w*Cjbe ;
        core = stampY(core, b, c, c1);
        core = stampY(core, b, e, c2);
        core = stampY(core, c, 0, c3);
    end

    core.isNonlinear=1;
end