% function slgkmodel()
clear global; clear functions;

b=zeros(20,20);
cr=0;
cs=zeros(20,1);
err=zeros(20,1);
hidden=0;
i=0;
ii=0;
in=zeros(20,1);
ipre=zeros(20,1);
is=zeros(20,1);
itotal=0;
j=0;
manypre=0;
novel=zeros(20,1);
novtotal=0;
or_fv=0;
pp=zeros(2000+1,20,20,20);
pre=zeros(20,1);
pwin=zeros(20,1);
r=zeros(20,20);
s=0;
sal=zeros(20,1);
sps=0;
spu=0;
sts=0;
stu=0;
t=0;
teacher=zeros(20,1);
tempi=zeros(20,1);
time_fv=0;
tipe=0;
tr=zeros(20,1);
trh=zeros(20,1);
trial=zeros(2000+1,1);
us_intensity=0;
v=zeros(20,20);
win=zeros(20,1);
x=zeros(20,1);
xp=zeros(20,1);
z=zeros(20,1);
z1=zeros(20,1);
zz=zeros(20,1);
zzz=zeros(20,1);

 
 
k=[.2,.1,.005,.02,.005,1.,2.,.4,.995,.995,.7,.5].';
time_units = 50; % # of time units
stim = 4; % # of stimuli including CX
ty = 3; % Types of trials
n = 60; % # of trials
trial_types_order = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; % Order of trial types
pn = 1; % Print from trial #
test_time = 39; % Test CR at time
hidden = 0; % Hidden
output_trials = zeros(n,1);
pre1_test_time = zeros(n,1); % CR values per trial at test time
output_time_units = zeros(time_units, n);
cs1_time_units = zeros(time_units, n); % CS(1) per time unit
pre1_keys_time_units(1:time_units, 1:n) =  {0}; % pre1 determined keys from us intensity values
%output_time_units_single = zeros(time_units * n, 1);
pre1_determined = 0;
pre1_determined_key = 0;
us_intensity_values(1:2000, 1:20) = {0}; % eg. us_intensity_values{1,5} = us_intensity_order{1};
 
% stim * ty many values should be given for:
% manypre_values, sts_values, sps_values, sal_values
%% # of presentations of stimulus
manypre_values = [1 1 1 1 1 1; 1 1 1 1 1 1; 1 1 1 1 1 1];
%% CS starts
sts_values = [1 30 10 20 20 20; 1 30 10 1 20 20; 1 20 20 20 20 20];
%% CS stops
sps_values = [2000 40 20 40 40 40; 2000 40 20 40 40 40; 2000 40 40 40 40 40];
%% CS salience
sal_values = [.1 0 0 0 0 0; .1 0 0 0 0 0; ; .1 0 1 0 0 0];
 
% ty many values should be given for:
% stu_values, spu_values, us_intensity_values
%% US starts
stu_values = [20 20 20];
%% US stops
spu_values = [40 40 40];
%% US intensity
% eg. us_intensity_order{1}(1) us_intensity_order{1}(2) ...
% [USnoact dUS]  dUS = Action Outcome = (USact - USnoact)*4
us_intensity_order = {[1 1]; [0 0] ; [0 0]};




itotal = stim + 1;
for tipe = 1: ty
    for s = 1: stim
        manypre = manypre_values(tipe, s);
        for time_fv = 1: time_units
            for ii = 1: manypre
                pp(time_fv+1,tipe,s,ii) = 0;
            end
            ii = manypre+1;
        end
        time_fv = time_units+1;
        
        for ii = 1: manypre
            sts = sts_values(tipe, s);
            win(ii) = fix(sts);
            sps = sps_values(tipe, s);
            pwin(ii) = fix(sps);
            sal(s) = sal_values(tipe, s);
            
            for time_fv = 1: time_units
                if(time_fv >= win(ii) && time_fv <= pwin(ii))
                    pp(time_fv+1,tipe,s,ii) = sal(s);
                end
            end
            time_fv = time_units+1;
        end
        ii = manypre+1;
    end
    
    s = stim+1;   
    stu = stu_values(tipe);
    spu = spu_values(tipe);
   
    for time_fv = 1: time_units
        us_intensity_values{time_fv+1,tipe} = 0;
        
        if(time_fv >= stu)
            us_intensity_values{time_fv+1,tipe} = us_intensity_order{tipe};
        end
        
        if(time_fv > spu)
            us_intensity_values{time_fv+1,tipe} = 0;
        end
    end

    time_fv = time_units+1;
end
tipe = ty+1;

if(pn == 0)
    pn = 1;
end

for t = 1: n
    i = trial_types_order(t); % Order of trial types
    trial(t+1) = fix(i);
end
t = n+1;

% j's of hidden go from itotal+1 to itotal+1+hidden
% connections of inputs to total to hidden
for j = itotal + 1: itotal + 1 + hidden
    for i = 1: itotal
        r(i,j) =(((1.-rand).*1.)-.5)./2;
    end
    i = fix(itotal+1);
end

j = itotal + 1 + hidden + 1;

for i = 1: itotal + hidden
    for j = 1: itotal
        v(i,j) = 0;
    end
    j = itotal + 1;
end
i = itotal + hidden + 1;

for t = 1: n
    for time_fv = 1: time_units
        us_intensity = us_intensity_values{time_fv+1,trial(t+1)};
        if (time_fv == stu)
            if (pre1_test_time(t,1) <= 0.01)% change between 0 and 0.001
                if (rand() > 0.2)
                   pre1_determined = us_intensity(1);
                   pre1_determined_key = 1;
                else
                   pre1_determined = us_intensity(2);
                   pre1_determined_key = 2;
                end
            else 
                if (rand() > 0.8)
                   pre1_determined = us_intensity(1);
                   pre1_determined_key = 1;
                else
                   pre1_determined = us_intensity(2);
                   pre1_determined_key = 2;
                end
            end
        end
        
        cs(1) = 0;
        if(time_fv >= stu)
            cs(1) = pre1_determined;
            pre1_keys_time_units{time_fv, t} = [num2str(trial(t+1), '%d') ',' num2str(pre1_determined_key, '%d')];
        end
        if(time_fv > spu)
            cs(1) = 0;
            pre1_keys_time_units{time_fv, t} = 0;
        end
        cs1_time_units(time_fv, t) = cs(1);
        
        for i = 1: stim
            cs(i+1) = 0;
            for ii = 1: manypre
                cs(i+1) = cs(i+1) + pp(time_fv+1,trial(t+1),i,ii);
            end
            ii = manypre+1;
        end
        i = stim+1;
        
        % Assign teaching signals of the output units
        for j = 1: itotal
            teacher(j) = cs(j);
        end
        j = itotal+1;
        
        % T and X
        for i = 1: itotal
            % Traces for direct connections
            tr(i) = tr(i) + k(1).*(cs(i)-tr(i));
            in(i) = tr(i) + k(8).*pre(i);
            tempi(i) = in(i).*z(i);
            if(tempi(i) < 0.)
                tempi(i) = 0.;
            end
            
            % Traces for hidden connections
            trh(i) = trh(i) + k(1).*cs(i).*(1-trh(i)) - 0.1.*k(1).*trh(i);
            % X for direct connections
            xii(i) = k(7).*(in(i).*k(2)+tempi(i));
            x(i) = abs(xii(i)); 
        end
        i = itotal+1;
        
       % x(1) = 0;
        % Output of hidden
        for j = itotal + 1: itotal + 1 + hidden
            xp(j) = 0;
            x(j) = 0;
        end
        j = fix(itotal + 1 + hidden+1);
        
        for j = itotal + 1: itotal + 1 + hidden
            for i = 1: itotal
                xp(j) = xp(j) + trh(i).*r(i,j);
            end
            i = itotal+1;
        end
        j = itotal + 1 + hidden + 1;
        
        % Changing attention to the hidden units
        for j = itotal + 1: itotal + 1 + hidden
            zz(j) = zz(j) + .00001.*novel(1).*(1-zz(j));
        end
        j = itotal + 1 + hidden + 1;
        
        for j = itotal + 1: itotal + 1 + hidden
            zzz(j) = zz(j).^8./(zz(j).^8+.1.^8);
            x(j) =(5.*xp(j).*zzz(j)).^2;
            if(x(j) <= 0)
                x(j) = 0;
            end
        end
        j = itotal + 1 + hidden + 1;
        
        % Aggregate Predictions and CR
        for i = 1: itotal + hidden
            for j = 1: itotal
                b(i,j) = x(i).*v(i,j);
            end
            j = itotal+1;
        end
        i = itotal + hidden + 1;
        
        for j = 1: itotal
            pre(j) = 0.;
            for i = 1: itotal + hidden
                pre(j) = pre(j) + b(i,j);
            end
            i = itotal + hidden + 1;
            
            if(pre(j) < 0)
                pre(j) = 0;
            end
        end
        j = itotal+1;
        
        cr = k(6).*(pre(1).^2./(pre(1).^2+.15.^2));
        if(cr < 0)
            cr = 0.;
        end
        
        % novelty
        novtotal = 0.;
        
        % integrator of CS/US and Prediction of CS/US
        for j = 1: itotal
            if(j == 1)
                is(1) =(k(9).*is(1)+cs(1).*(1-is(1)));
                ipre(1) =(k(10).*ipre(1)+pre(1).*(1-ipre(1)));
            else
                is(j) =(k(9).*is(j)+cs(j).*(1-is(j)));
                ipre(j) =(k(10).*ipre(j)+pre(j).*(1-ipre(j)));
            end
            
            % novelty is the difference between the two integrators
            novel(j) = abs(is(j)-ipre(j));
            % Sum of the individual novelty
            novtotal = novtotal + novel(j);
        end
        j = itotal+1;
        
        or_fv =(novtotal.^2./(novtotal.^2+.75.^2));
        cr = cr .* (1 - or_fv);
        % error signal
        for j = 1: itotal
            err(j) =(teacher(j)-pre(j));
        end
        j = itotal+1;
        
        % Calculate Weight Change
        % V
        for i = 1: itotal + hidden
            for j = 1: itotal
                v(i,j) = v(i,j) + .6.*k(3).*x(i).*err(j).*(1.0-abs(v(i,j)));
            end
            j = itotal+1;
        end
        i = itotal + hidden + 1;
        
        % Z
        for i = 1: itotal
            z(i) = z(i) +(in(i).*((k(4).*or_fv.*(1.-abs(z(i)))-k(5).*(1.+z(i)))));
            if(z1(i) < -1.)
                z1(i) = -1.;
            end
        end
        i = itotal+1;
        
        % Output results
        if(time_fv == test_time && t >= pn)
            pre1_test_time(t,1) = pre(1); % DONT change this
            output_trials0(t,1) = err(1);  % change this
            output_trials(t,1) = x(1);  % change this
            output_trials2(t,1) = v(3,1);  % change this
            output_trials3(t,1) = v(1,3);  % change this
            output_trials4(t,1) = b(1);  % change this
            output_trials5(t,1) = novtotal;  % change this
            output_trials6(t,1) = cr;  % change this
        end
      
        if(t >= pn)
            output_time_units0(time_fv, t) = err(1); % change this
            output_time_units(time_fv, t) = x(1); % change this
            output_time_units2(time_fv, t) = v(3,1); % change this
            output_time_units3(time_fv, t) = v(1,3); % change this
            output_time_units4(time_fv, t) = b(1); % change this
            output_time_units5(time_fv, t) = novtotal; % change this
            output_time_units6(time_fv, t) = cr; % change this
        end
    end
    time_fv = time_units+1;
end

output_time_units_single0 = output_time_units0(:);
output_time_units_single = output_time_units(:);
output_time_units_single2 = output_time_units2(:);
output_time_units_single3 = output_time_units3(:);
output_time_units_single4 = output_time_units4(:);
output_time_units_single5 = output_time_units5(:);
output_time_units_single6 = output_time_units6(:);
disp(output_trials);
disp(output_time_units);

plot(output_trials);
%plot(output_time_units)
%plot(output_time_units_single)
% clear all;
% end
