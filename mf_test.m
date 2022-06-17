save_folder = 'test\';

kStep = 0.001;
k = -2:kStep:(2-kStep);
v = 0.2*rand(size(k));

s1 = 2.^(-(k.^2)./2);
s1 = 0.5*sin(k*4);
s2 = 2.^(-(k.^2)./4);
sig1 = s1 + 0.5*(randn(length(s1),1));
sig2 = s2 + 0.5*(randn(length(s2),1)); %+ 1j*randn(length(s2),1));

signal = sig1; %declaring signal
template = s1; %declaring template
length_template = length(template); %length of template
length_signal = length(signal); %length of signal
length_convolution = length_signal + length_template - 1; %length of convolution

mf_output = zeros(2*length_signal, 1);

subplot(3,1,1)
plot(sig1);

subplot(3,1,2)
plot(s1);

for I = 1:length_signal %zeroing mf_output array
    mf_output(I) = 0;
end

for I = 1:length_signal %iteration loop for each index of input signal
    for J = 1:length_template %iteration loop for each index of the template
        mf_output(I + J) = mf_output(I + J) + (signal(I) * template(J));
    end
end

subplot(3,1,3)
plot(mf_output);
saveas(gcf,strcat('figures\',save_folder,'mf_output.png'))
writematrix(mf_output, strcat('data\',save_folder,'mf_output.csv'))