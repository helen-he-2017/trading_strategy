% Simulate trading strategy of Option Trading by Euan Sinclair, Chapter 10

%% Simulate Geometric Brownian Motion of Stock Price
N = 1e3;
r = 1;
alpha = 0.1;
T = 1;
npaths = 1e3;        % Number of simulations

t = T*(0:1:N)/N;     % Time vector
y0 = ones(npaths,1); % Vector of initial conditions, must match number of paths
opts = sdeset('RandSeed',0,'SDEType','Ito'); % Set seed
y = sde_gbm(r,alpha,t,y0,opts);

figure;
plot(t,y,'b',t,y0*exp(r*t),'r--');
xlabel('t');
ylabel('y(t)');

%% Strategy Input
%  The initial price is 99/101.
%  We always need to make a price that is no more than 2 ticks wide.
%  All trades are in single units (this is not essential, but makes analysis
% easier).
%  After a buy, we lower both the bid and the offer by one.
%  After a sale, we raise both the bid and the offer by one.


n_t = 1000;
n_sim = 10000;

bid = zeros(n_t,1);
ask = zeros(n_t,1);
inventory = zeros(n_t,1);
inventory_end = zeros(n_sim,1);
bid(1) = 99;
ask(1) = 101;
t= 1:n_t;
profit = zeros(n_sim,1);

test = 0;
for j = 1:n_sim
    for i = 2:n_t
        test = randn();
        if test >= 0
            bid(i) = bid(i-1) + 1;
            ask(i) = ask(i-1) + 1;
            inventory(i) = inventory(i-1)-1;
            profit(j) = profit(j) + ask(i-1);
        else
            bid(i) = bid(i-1) - 1;
            ask(i) = ask(i-1) - 1;
            inventory(i) = inventory(i-1)+1;
            profit(j) = profit(j) - bid(i-1);
        end
    end
    inventory_end(j) = inventory(end);
    profit(j) = profit(j)+inventory(end)*(bid(end)+ask(end))/2;
end


mid = (bid+ask)/2;
disp('profit');
disp(mean(profit));
disp('inventory');
disp(mean(inventory_end));


figure
plot(t,mid);


hist(profit,20);

figure
plot(inventory_end,profit,'.');

%% Strategy Input - Information-based Market Making
%  We start with an idea of the instrument?s value, S, and also some uncertainty
% associated with this value, ?S. That is, we have our best estimate,
% but we also do not claim that it is perfect.
%  We quote bid and ask prices around this value, Sb and Sa.
%  We do a trade at one of these prices. Let?s assume we buy at Sb. This
% gives us information (someone thinks our price is too high), and of
% course inventory (we now have a long position).
%  Now we need to update our estimate of value.
clc
n_sim = 10000;
n_t = 100;
sigma_S = zeros(n_t,1);
sigma_S(1) = 5;
mid = zeros(n_t,1);
bid = zeros(n_t,1);
ask = zeros(n_t,1);
inventory = zeros(n_t,1);
inventory_end = zeros(n_sim,1);
bid(1) = 99;
ask(1) = 101;
mid(1) = (bid(1)+ask(1))/2;
t= 1:n_t;
profit = zeros(n_sim,1);
T = zeros(n_t,1);
Tmax = 20;
sigma = 20;
k = zeros(n_t,1);
sigma_e = zeros(n_t,1);

for j = 1:n_sim
    for i = 2:n_t
        T(i) = sigma*randn();
%         if i == 2
%             Tmax = 50;
%         else
%             Tmax = max(abs(T));
%         end
        
        
        if T >= 0
            inventory(i) = inventory(i-1) - T(i);
            if T(i) > Tmax
                sigma_e(i) = sigma_S(i-1)* (Tmax/abs(T(i)) - 1);
                k(i) = sigma_S(i-1)^2/(sigma_S(i-1)^2+sigma_e(i)^2);
                mid(i) = mid(i-1) + k(i)*(ask(i-1)-mid(i-1));
                sigma_S(i) = sqrt(1-k(i))*sigma_S(i-1);
                bid(i) = mid(i) - sigma_S(i);
                ask(i) = mid(i) + sigma_S(i);
            else
                sigma_e(i) = sigma_e(i-1);
                sigma_S(i) = sigma_S(i-1);
                k(i) = k(i-1);
                bid(i) = bid(i-1);
                ask(i) = ask(i-1);
                mid(i) = mid(i-1);
               
            end
        profit(j) = profit(j) + ask(i-1)*T(i);
        else
            inventory(i) = inventory(i-1) - T(i);
            if abs(T(i)) > Tmax
                sigma_e(i) = sigma_S(i-1)* (abs(Tmax)/abs(T(i)) - 1);
                k(i) = sigma_S(i-1)^2/(sigma_S(i-1)^2+sigma_e(i)^2);
                mid(i) = mid(i-1) + k(i)*(mid(i-1)-bid(i-1));
                sigma_S(i) = sqrt(1-k(i))*sigma_S(i-1);
                bid(i) = mid(i) - sigma_S(i);
                ask(i) = mid(i) + sigma_S(i);
            else
                sigma_e(i) = sigma_e(i-1);
                sigma_S(i) = sigma_S(i-1);
                k(i) = k(i-1);
                bid(i) = bid(i-1);
                ask(i) = ask(i-1);
                mid(i) = mid(i-1);
            end
        profit(j) = profit(j) - bid(i-1)*T(i);

        end

    end
    inventory_end(j) = inventory(end);
    profit(j) = profit(j)+inventory(end)*T(end);
end

disp('profit');
disp(mean(profit));
disp('inventory');
disp(mean(inventory_end));


figure
plot(t,mid);


hist(profit,20);

figure
plot(inventory_end,profit,'.');
