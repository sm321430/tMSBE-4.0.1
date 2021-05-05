load test
dp=p_qw-p_fm1_qw;
dp_abs=abs(p_qw)-abs(p_fm1_qw);
dp_abs_rel=abs(dp_abs(1:10:end,1:28))./abs(p_qw(1:10:end,1:28));
dp_rel=abs(dp(1:10:end,1:28))./abs(p_qw(1:10:end,1:28));

dI=Inv_qw-Inv_00_qw;
dI_rel=abs(dI(1:10:end,1:28))./abs(Inv_qw(1:10:end,1:28));

figure('Name', 'Micro pol. p_fm1 for MockLinear')
contourf(k_QW6(1:28),t(1:10:end)/ps,log(abs(p_fm1_qw(1:10:end,1:28))));
grid on;
xlabel('k [1/a_0]');
ylabel('t [ps]');
colorbar;

figure('Name', 'Micro pol. for linear')
contourf(k_QW6(1:28),t(1:10:end)/ps,log(abs(p_fm1_qw(1:10:end,1:28))));
grid on;
xlabel('k [1/a_0]');
ylabel('t [ps]');
colorbar;

figure('Name', 'Inversion for MockLinear')
contourf(k_QW6(1:28),t(1:10:end)/ps,log(abs(Inv_00_qw(1:10:end,1:28))));
grid on;
xlabel('k [1/a_0]');
ylabel('t [ps]');
colorbar;

figure('Name', 'Inversion total error')
contourf(k_QW6(1:28),t(1:10:end)/ps,log(abs(dI(1:10:end,1:28))));
grid on;
xlabel('k [1/a_0]');
ylabel('t [ps]');
colorbar;

figure('Name', 'Inversion relative error')
contourf(k_QW6(1:28),t(1:10:end)/ps,log(abs(dI(1:10:end,1:28))));
grid on;
xlabel('k [1/a_0]');
ylabel('t [ps]');
colorbar;

figure('Name','Error in micro pol. total')
contourf(k_QW6(1:28),t(1:10:end)/ps,log(abs(dp(1:10:end,1:28))));
xlabel('k [1/a_0]');
ylabel('t [ps]');
grid on;
colorbar;

figure('Name','Error in micro pol. magnitude')
contourf(k_QW6(1:28),t(1:10:end)/ps,log(abs(dp_abs(1:10:end,1:28))));
xlabel('k [1/a_0]');
ylabel('t [ps]');
grid on;
colorbar;

figure('Name','Relative error in micro pol.')
contourf(k_QW6(1:28),t(1:10:end)/ps,log(dp_rel));
xlabel('k [1/a_0]');
ylabel('t [ps]');
grid on;
colorbar;

figure('Name','Relative error in micro pol. mag')
contourf(k_QW6(1:28),t(1:10:end)/ps,log(dp_abs_rel));
xlabel('k [1/a_0]');
ylabel('t [ps]');
grid on;
colorbar;

figure('Name','Relative error in micro pol. overlay')
hold on
plot(t(1:10:end)/ps,log(dp_rel(:,3)));
plot(t(1:10:end)/ps,log(dp_rel(:,7)));
plot(t(1:10:end)/ps,log(dp_rel(:,14)));
plot(t(1:10:end)/ps,log(dp_rel(:,21)));
lgd=legend(['k=',num2str(k_QW6(3))],['k=',num2str(k_QW6(7))],['k=',num2str(k_QW6(14))],...
    ['k=',num2str(k_QW6(21))]);
xlabel('t [ps]');
ylabel('log(|\Delta p_k|/|p_{k,lin}|)');