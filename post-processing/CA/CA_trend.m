%% Contact Angle Estimation from LAMMPS trajectory
%
% M. Provenzano – Politecnico di Torino, Italy - Aug 2024
%
% Loads wet_*.mat (from CA_data_parser.m) and plots contact angle versus
% trajectory frame with uncertainty bounds.

load wet_0.mat

kk = 1;
CA_vec = zeros(1000,1);
CA_dev_vec = zeros(1000,1);

for ii = 1:size(tot,1)
    CA_vec(kk) = tot(ii).CA_tan;
    CA_dev_vec(kk) = tot(ii).CA_dev;
    kk = kk+1;
end

load wet_1.mat

kk = kk - 1

for ii = 1:size(tot,1)
    CA_vec(kk) = tot(ii).CA_tan;
    CA_dev_vec(kk) = tot(ii).CA_dev;
    kk = kk+1;
end

load wet_2.mat

kk = kk - 1

for ii = 1:size(tot,1)
    CA_vec(kk) = tot(ii).CA_tan;
    CA_dev_vec(kk) = tot(ii).CA_dev;
    kk = kk+1;
end

CA_vec(kk:end) = [];
CA_dev_vec(kk:end) = [];

num_el = round(length(CA_vec)*0.5)+1;
ref_value = mean(CA_vec(num_el:end));
figure
plot([0;length(CA_vec)],[ref_value;ref_value],'-r')
hold all

ni = 16;
k = length(CA_vec)-num_el+1;
denom = ni*k-1;
num = 0;
for ii = num_el:length(CA_vec)
    num = num+CA_dev_vec(ii)^2*(ni-1)+ni*(CA_vec(ii)-ref_value)^2;
end
error_bar = sqrt(num/denom);
if error_bar > 20
    warning('ContactAngle:highUncertainty', 'Contact angle uncertainty is high: %.1f°.\nAverage contact angle may be unreliable.\n',error_bar)
end

for ii = 1:length(CA_vec)
    plot([ii,ii],[CA_vec(ii)+CA_dev_vec(ii),CA_vec(ii)-CA_dev_vec(ii)],'-c')
    hold all
end

plot(1:length(CA_vec),CA_vec,'.-b','MarkerSize',5)
grid on
hold all
legend(['Errorbar: ',num2str(error_bar),'°'])
title('Wettability')
ylabel('Angle (°)')

T = table(CA_vec, CA_dev_vec, 'VariableNames', { 'CA', 'std'} );
%writetable(T, 'wet.txt')

