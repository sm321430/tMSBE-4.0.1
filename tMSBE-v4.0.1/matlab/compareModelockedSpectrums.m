function compareModelockedSpectrums(savedData,saveKey)
    setupPlot
    setupConstants
    tmp_fig=figure(101);
    w0 = 1.2037; %Central frequency in eV
    set(tmp_fig,'Name','Pulse spectrum comparison');
    hold on
    
    max_spectrum=0;
    for j=1:length(savedData)
       loadedData(j)=load(char(savedData(j).key));
       if max(abs(loadedData(j).pulse_spectrum1_integrated))>max_spectrum
        max_spectrum=max(abs(loadedData(j).pulse_spectrum1_integrated));
       end
    end
    
    for j=1:length(savedData)
       plot(hbar*loadedData(j).w_qw/e,abs(loadedData(j).pulse_spectrum1_integrated)/max_spectrum,...
           'DisplayName',char(savedData(j).legend));
    end
    plot(w0*[1,1],[0,1],'k--','DisplayName','\omega_0');
    ylim([0,1]);
    xlim([1.17,1.24]);
    set(gca,'YTick', [0,0.25,0.5,0.75,1]);
    xlabel('Energy [eV]');
    ylabel('Spectrum [a.u.]');
    legend('show');
    grid on
    saveas(tmp_fig,[saveKey,'pulseSpectrumComparison.png']); 
end
