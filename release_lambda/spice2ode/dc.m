function core  = dc( config, core )
    diff=1;
    iteration=0;
    
    % Main loop for non-linear DC Analysis using Newton-Raphson Method
    while (diff>config.DCErrorThreshold)
        %do stamping and stuff
        core = stamp(config, core);
        disp('DC iterations ...');
%         core.vNodeset = [+1.00000000000e+00
% 1.00056805948e+00
% -3.06337886154e-01
% +1.43292813284e+01		
% -1.42345556261e+01
% +1.42350559196e+01
% -1.36756936314e+01
% -1.34854441661e+01
% -1.43378483480e+01
% -1.49906487031e+01
% -1.49906487031e+01
% +3.46831056923e-01
% +3.47115086665e-01
% -1.41625109432e+01
% -1.49277097252e+01
% -4.85651138006e-01
% +1.72212380071e+00
% +1.04297243622e+00
% +2.79035671847e-01
% -1.49999999999e+01
% +9.97381318928e-01
% +1.00056805948e+00
% +1.00375736380e+00
% -1.49055366526e+01
% +1.50000000000e+01		
% -1.50000000000e+01	];
        
        %first approach:
        core.operatingPoint = core.Y\core.J ; % this is faster than dcPointVoltage = inv(Y)*J ;
%            disp('Y');
%         core.Y
%         disp('J');
%         core.J
%         disp('X');
%         core.operatingPoint
%         
%         pause
    
        %   [L,U] = lu(Y);
        %   y = L\J;
        %   x = U\y;    %internally, matlab uses LU factorization in Y\J.
        %Second Approach: Using LU Factorization and Guassian Elimination

%         core.operatingPoint = [
%         +1.00000000000e+00
% 1.00056805948e+00
% -3.06337886154e-01
% +1.43292813284e+01		
% -1.42345556261e+01
% +1.42350559196e+01
% -1.36756936314e+01
% -1.34854441661e+01
% -1.43378483480e+01
% -1.49906487031e+01
% -1.49906487031e+01
% +3.46831056923e-01
% +3.47115086665e-01
% -1.41625109432e+01
% -1.49277097252e+01
% -4.85651138006e-01
% +1.72212380071e+00
% +1.04297243622e+00
% +2.79035671847e-01
% -1.49999999999e+01
% +9.97381318928e-01
% +1.00056805948e+00
% +1.00375736380e+00
% -1.49055366526e+01
% +1.50000000000e+01		
% -1.50000000000e+01
% -2.33876888433e-03
% +2.33886280095e-03
% -9.49541493439e-08
% -9.39166219288e-08    
%         ];
        %% Printing out some data:
        if( diff == 1 ),
            [s t] = size(core.operatingPoint);
%             for i=1:s
%                 disp( sprintf('[DC] %s = ', core.V{i} ) )
%                 core.operatingPoint(i)
%             end
        end
        
        d = core.vNodeset - core.operatingPoint(1:core.numNode) ;
        core.vNodeset = core.operatingPoint(1:core.numNode) ;
        core.operatingPoint
        diff = sqrt(d'*d);
        iteration = iteration + 1; 
        core.dcIteration = iteration;
        disp( sprintf('iteration %d, difference=%f', iteration, diff))
        if(core.isNonlinear==0) diff=0; end
    end
    
     disp( sprintf('converged after %d iteration, final error=%f', iteration, diff))

end

