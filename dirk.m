clear all
% Solve ODE: y'(t) = f(y,t)

f 	   = @(y,t) y+t;
fprime = @(y,t) 1;
exact  = @(t) 1*exp(t)-t-1;
y0 	   = 0;

% f 	   = @(y,t) y*t;
% fprime = @(y,t) t; 
% exact  = @(t) exp(t.^2/2);
% y0 	   = 1;

% f 	   = @(y,t) t.^2;
% fprime = @(y,t) 0; 
% exact  = @(t) 1 + (t.^3)/3;
% y0 	   = 1;

% f 	   = @(y,t) t.^3;
% fprime = @(y,t) 0; 
% exact  = @(t) 1 + (t.^4)/4;
% y0 	   = 1;

% f 	   = @(y,t) t.^4;
% fprime = @(y,t) 0; 
% exact  = @(t) 1 + (t.^5)/5;
% y0 	   = 1;

% f      = @(y,t) -10000*(y - sin(t + pi/4)) + cos(t + pi/4);
% fprime = @(y,t) -10000;
% exact  = @(t) (sin(t) + cos(t))/sqrt(2);
% y0     = sin(pi/4);

NT     = 5;
T      = 1;
tol    = 1e-15;

dT = ones(1,NT)*0.1;
for i = 2:length(dT)
	dT(i) = dT(i-1)/2;
end

lines = ['-ko';'-bo';'-ro'];

%%% Loop over methods
for method = 1:3 

	switch method
		case 1
			% --- Alexander DIRK(2,2) Butcher tableau
			alpha = 0.5*(2-sqrt(2));
			A = [   alpha,   0  ;
				  1-alpha, alpha];
			c = [   alpha,   1  ];
			b = [ 1-alpha, alpha];
			% ---
		case 2
			% --- Cash DIRK(3,3) Butcher tableau
			alpha = 0.43586652150845899941601945;
			tau2  = (alpha^2 - 3*alpha/2 + 1/3)/(alpha^2 - 2*alpha + 0.5);
			b1    = (0.5*tau2  - 1/6)/( (tau2 - alpha)*(1 - alpha) );
			b2    = (0.5*alpha - 1/6)/( (alpha - tau2)*(1 - tau2 ) );
			A = [    alpha    ,   0  ,   0  ;
				  tau2 - alpha, alpha,   0  ;
				      b1      ,  b2  , alpha];
			c = [    alpha    , tau2 ,   1  ];
			b = [     b1      ,  b2  , alpha];
			% ---
		case 3
			% --- Hairer & Wanner DIRK(5,4) Butcher tableau
			A = [  1/4   ,      0   ,     0  ,     0 ,   0;
				   1/2   ,    1/4   ,     0  ,     0 ,   0;
				  17/50  , -  1/25  ,   1/4  ,     0 ,   0;
				 371/1360, -137/2720,  15/544,   1/4 ,   0;
				  25/24  , - 49/48  , 125/16 , -85/12, 1/4];
			c = [  1/4   ,    3/4   ,  11/20 ,   1/2 ,   1];
			b = [ 25/24  , - 49/48  , 125/16 , -85/12, 1/4];
			% ---
		otherwise
	end

	%%% Loop over time step sizes (NT)
	for N = 1:length(dT)   

		dt 		   = dT(N);
		nStages    = length(A);
		nTimesteps = T/dt;

		% time
		t 		   = 0;			
		tn 		   = t;

		% Solution vector
		y = zeros(nTimesteps+1,1); 	
		y(1) = y0;

		% Intermediate solution vector
		prevsol = zeros(nStages+1,1); 

		%%%% Loop over number timesteps
		for i = 2:(nTimesteps+1)   

			% Save time
			told = t;

			fprintf("\nTIMESTEP: %d\t TIME = %.3f\n",(i-1),(told+dt));

			prevsol(1) = y(i-1);		% solution from previous timestep
			yi   	   = prevsol(1);	% initial guess 

			%%% Loop over DIRK stages
			for s = 1:nStages

				fprintf("\tStage %d\n",s);

				% Get DIRK coeffs
				coeffs = zeros(s+1,1);
				coeffs(1) = 1;
				coeffs(2) = 1;
				if (s > 1)
					for q = 3:length(coeffs)
						coeffs(q) = A(s,q-2);
					end
				end
				for q = 1:length(coeffs)
					coeffs(q) = coeffs(q)/A(s,s);
				end

				c_j = c(s);
				t = told + c_j*dt;

				% Newton
				deltaY = 2*tol;

				fprintf("\t|  #Newton | resDelta \n");
				n = 1;

				while abs(deltaY) > tol
					
					% Assembly "matrix"
					RHS = - coeffs(1)*yi/dt + f(yi,t); 
					for k = 1:s
						RHS = RHS + coeffs(k+1)*prevsol(k)/dt;
					end

					% Assembly residual
					LHS = coeffs(1)*1/dt - fprime(yi,t);

					% Solve
					deltaY = RHS/LHS;

					% Update
					yi = yi + deltaY;

					fprintf("\t|\t%d  | %.15f \n",n,abs(deltaY));
					n = n + 1;

				end


				% Store intermediate solution
				prevsol(s+1) = coeffs(1)*yi;
				for k = 1:s
					prevsol(s+1) = prevsol(s+1) - coeffs(k+1)*prevsol(k);
				end
			end

			% Final step
			y(i) = yi;
			% Not needed (stiffly-accurate)
			% y(i) = prevsol(1);
			% for k = 1:nStages
			% 	y(i) = y(i) + b(k)*prevsol(k+1);
			% end

			% Increase time
			t = told + dt;
			tn(i) = t;

		end

		err(N) = abs(y(end) - exact(t));

	end

	% exactSol = exact(tn);
	% plot(tn,y,'-ko',tn,exactSol,'-r')
	% legend('DIRK','exact')
	% plot(tn,abs(y-exactSol'))

	loglog(dT,err,lines(method,:))
	hold on

	for i = 2:length(err)
		log(err(i-1)/err(i))/log(2)
	end
end
set(gca,'FontSize',20);
legend('DIRK(2,2)','DIRK(3,3)','DIRK(5,4)');
xlim([min(dT)*0.5,max(dT)*2]);