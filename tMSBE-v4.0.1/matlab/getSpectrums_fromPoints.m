function [w_qw, Y1_qw, Y2_qw] = getSpectrums_fromPoints(t_cur,sizeIl,sizeIr,points,E_cav)
    setupConstants

    x0 = t_cur-points(1);
    indp = x0>=-sizeIl/ps;
    indm = x0<=sizeIr/ps;
    ind3 = indp==indm;
    t0_qw = t_cur(ind3)*ps;
    y0_qw = E_cav(ind3);

    x0 = t_cur-points(2);
    indp = x0>=-sizeIl/ps;
    indm = x0<=sizeIr/ps;
    ind3 = indp==indm;
    t1_qw = t_cur(ind3)*ps;
    y1_qw = E_cav(ind3);
    
    y0_qw=transpose(y0_qw);
    y1_qw=transpose(y1_qw);

    
    [w_qw,Y1_qw,Y2_qw] = getSpectrums(y0_qw,y1_qw,t0_qw(2)-t0_qw(1),t1_qw(2)-t1_qw(1));

    %     plot(hbar*w_qw/e,abs(Y1_qw),hbar*w_qw/e,abs(Y2_qw))
%     grid on
%     xlim([1.17,1.23])
%     xlabel('Energy [eV]');
%     ylabel('Spectrum [a.u.]');
%     legend(['t-t_0=',num2str(points(1),2),'ps'],['t-t_0=',num2str(points(2),2),'ps']);
end