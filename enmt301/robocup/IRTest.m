clearvars -except IRLog2
t = (0:0.01:(width(IRLog2)/100)-1/100);
vals = table2array(IRLog2);
vals = double(vals);
plot(t, vals)
ylabel('Arduino ADC output')
xlabel('Time (s)')