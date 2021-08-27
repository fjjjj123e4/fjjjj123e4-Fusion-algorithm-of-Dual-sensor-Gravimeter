clear all
Ts = 60;
t = 0.05:0.05:3600;

dr = 1e-17 * t.^4 - 1e-14 * t.^3 - 8e-12 * t.^2 - 0.0000012 * t;
acc_ESA = 9.8 + 0.000001 * randn(1, 20 * 3600) - dr + 0.0001 * randn(1, 20 * 3600);
acc_RIMU = 9.8 + sin(2 * pi / Ts * t) * 0.000323 + 0.0001 * randn(1, 20 * 3600) -...
    cos(2 * pi / Ts * t) * 0.00042 + 0.0001 * randn(1, 20 * 3600) + 0.1 * sin(2 * pi / Ts * t);
acc_ESA = smooth(acc_ESA,600);
for ii = 1:60
   acc_RIMU2(ii) = mean(acc_RIMU((ii-1)*1200+1:ii*1200));
end
c = [8e-17; -5e-14; -2e-12; -0.0000012];
p1 = eye(4);
p2 = 1;
for jj = 1:length(acc_ESA)
    H =[(jj/20)^4, (jj/20)^3, (jj/20)^2, jj/20];
    if mod(jj, 1200) == 0
        [c, p1] = fusion_KF(c, p1, acc_RIMU2(jj/1200), acc_ESA(jj), jj/1200);
    end
    drift = H * c;
    acc_new(jj) = acc_ESA(jj) + drift;
end

p = 0;
acc_new1(1) = acc_new(1);
for jj = 2:length(acc_new)
    [acc_new1(jj), p] = noise_KF(acc_new1(jj-1), p, acc_new(jj));
end


%% plot 

Fig = figure(...
    'Units',           'pixels',...
    'Name',            'simu_static_ESA',...
    'NumberTitle',     'off',...
    'IntegerHandle',   'off');

AxesH = axes(...
    'Parent',  Fig,...
    'Xlim',    [0 3600],...
    'XGrid',   'on',...
    'YGrid',   'on',...
    'Visible', 'on');


realtrace = line(AxesH, t, acc_ESA, 'linewidth',2);
precisetrace = line(AxesH, [0, 3600], [9.8, 9.8],...
                    'color','r',...
                    'linewidth',2,...
                    'linestyle','--');
xlabel(AxesH, 't/s');
ylabel(AxesH, 'g/(m/s^2)');
legend(AxesH,[realtrace,precisetrace],...
    'ESA Output','True Value',...
    'location','northwest');
% edit the saving path and format of image
saveas(Fig,['E:\E盘\论文相关\IEEE_Sensor_Journal\simu_code\simu\image\' Fig.Name],'pdf')


