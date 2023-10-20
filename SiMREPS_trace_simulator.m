% SiMREPS Trace Simulator
%
% Takes input intensity states and transition probabilities, total
% fluorescence counts, frames, and noise level, and simulates SiMREPS traces
% using a kinetic Monte Carlo Approach.
%
% Alex Johnson-Buck, 3/19/2014
%

% *******************Parameters to set*************************

noise_level = 0.15; % Noise level (currently a constant multiple of totcounts, and independent of channel)
totcounts = 1000; % Total intensity counts for fluorophore in bound state
frames = 1200; % Number of frames
N = 100; % Number of traces to simulate
exposure_time = 0.5; % Time corresponding to each movie frame

framerate = 1./exposure_time;

plotting = 1; % Plot individual traces?  Set = 1 for yes.

%Define bound state lifetimes, in seconds
tau_on = 10;
tau_off = 50;

% *******************End of Parameters to Set*************************

noise = noise_level*totcounts;

tau_on = tau_on./exposure_time; % Convert lifetimes from seconds to frames
tau_off = tau_off./exposure_time;

%Define intensity states
states(1,1) = 1;
states(2,1) = 0;

% Transition probabilities
% State 1
P11 = exp(-1/tau_on);
P12 = 1-P11;
% State 2
P22 = exp(-1/tau_off);
P21 = 1-P22;

P = [P11 P12; P21 P22];
traces = zeros(N,frames); % Initialize trace matrix

for i = 1:size(P,1)
    P(i,:) = P(i,:)./sum(P(i,:));
end

for m = 1:N

    bound_state_true = zeros(frames,1);

    k = unidrnd(2); % Start in random state k
    bound_state_true(1,1) = states(k,1);

    decider = random('uniform',0,1,frames); % Random variable for determining transition points

    for frame = 2:frames
        if decider(frame) <= P(k,1)
            k = 1;
        else
            k = 2;
        end
        bound_state_true(frame,1) = states(k,1);
    end

    Ival = bound_state_true.*totcounts + normrnd(0,noise,frames,1);

    if plotting == 1
      if m == 1
        disp('Press Q in figure window to quit, or any other key to advance to the next simulated trace.');
      end
        f1 = figure(1);
        plot(Ival,'k-','LineWidth',1);
        xlabel('Frame');
        ylabel('Intensity');
        title(strcat('Simulated Trace:',num2str(m)));
        set(gca,'LineWidth',2,'FontSize',20);
        waitforbuttonpress;
        key1 = get(f1,'CurrentCharacter');
        if strcmp(key1,'q')==1 || strcmp(key1,'Q')==1
          disp('Quitting plotting and generating remaining traces...');
          close(f1);
          plotting=0;
        end
    end

    traces(m,:)=Ival;

    if mod(m,25)==0
        disp(strcat('Done generating trace...',num2str(m)));
    end

end
