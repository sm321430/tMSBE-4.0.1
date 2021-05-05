%% Test script for visualizing simple exponential recovery

setupConstants;
setupPlot;

tau_recov=20;
t_rt=44.2;
l_short=t_rt/8;
l_long=t_rt-l_short;
l_equal=t_rt/2;
num_rt=30;
t=linspace(0,t_rt*num_rt,num_rt*100);
hole_depth=0.1;

unequal=ones(1,num_rt*100);
for j=1:num_rt
   unequal=unequal-burnHole(t,(j-1)*t_rt+l_long,hole_depth,tau_recov); 
   unequal=unequal-burnHole(t,(j-1)*t_rt+l_long+l_short,hole_depth,tau_recov);
end

equal=ones(1,num_rt*100);
for j=1:num_rt
    equal=equal-burnHole(t,(j-1)*t_rt+l_equal,hole_depth,tau_recov);
    equal=equal-burnHole(t,j*t_rt,hole_depth,tau_recov); 
end

unequal_ctr=0;
equal_ctr=0;
unequal_drops=zeros(1,num_rt*2-1);
equal_drops=zeros(1,num_rt*2-1);
%Capture all drops in density
for k=1:num_rt*100-1
    if unequal(k+1)<unequal(k)
       unequal_ctr=unequal_ctr+1;
       unequal_drops(unequal_ctr)=unequal(k);
    end
    if equal(k+1)<equal(k)
       equal_ctr=equal_ctr+1;
       equal_drops(equal_ctr)=equal(k);
    end
end

figure
hold on
plot(t,unequal,'DisplayName','Unequal');
plot(t,equal,'DisplayName','Equal');
hold off
grid on
legend show
xlim([0,4*t_rt]);
xlabel('t [ps]');
ylabel('Density [a.u.]');

figure
hold on
plot(unequal_drops,'DisplayName','Unequal');
plot(equal_drops,'DisplayName','Equal');
hold off
grid on
legend show
xlabel('Pulse Hit');
ylabel('Density [a.u.]');

unequal_mean=mean(unequal_drops(40:49))
equal_mean=mean(equal_drops(40:49))

function hole=burnHole(t,t_burn,depth,recov)
    hole=depth.*heaviside(t-t_burn).*exp(-(t-t_burn)/recov);
end