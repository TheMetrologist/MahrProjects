clear all
load scale_error.mat
%loaded variables with data:
    %error1
    %error2
    %error3
    %error4
    %travel

%variables
period_scale = 4;   %period of scale in um
nbins = 40;     %number of bins for histogram
step_scale = 0.001; %distance in um to step scale
nmeas = 20; %number of diameters to sample
sigma = .120; %sigma of error of ball in um
% nball = 100;  %number of balls to sample
% step = period_scale / nball;
step = 0.508;    %distance to step ball, um

%fit error data to spline
xx = 0 : step_scale : 4;  %define x data for last 3 curves
    e1 = spline(travel(1:length(error1)), error1, xx);      %fit spline to first curve
    e2 = spline(travel(1:length(error2)), error2, xx);      %fit spline to each curve
    e3 = spline(travel(1:length(error3)), error3, xx);
    e4 = spline(travel(1:length(error4)), error4, xx);

%correct for drift so ends of error splines line up
drift1 = (e1(end) - e1(1)) / length(e1);   %incremental drift
drift2 = (e2(end) - e2(1)) / length(e2);   %total drift
drift3 = (e3(end) - e3(1)) / length(e3);   %total drift
drift4 = (e4(end) - e4(1)) / length(e4);   %total drift

n = 1:length(e1);
e1 = e1 - drift1*n;
e2 = e2 - drift2*n;
e3 = e3 - drift3*n;
e4 = e4 - drift4*n;
eavg = (e1 + e2 + e3 + e4) / 4;  %average splines together

%make 3 pitch long error curve
totalerror = [eavg eavg(1:end-1) eavg(1:end-1)];
totaltravel = -4:step_scale:8;

% figure(1)
%     plot(xx1, e1, xx, e2, xx, e3, xx, e4)
%     title('Individual Spline Error Curves')
%     xlabel('Scale Displacement, um'); ylabel('Error, um')
%     legend('err1','err2','err3','err4')
% 
% figure(2)
%     subplot(4, 1, 1)
%         hist(e1, nbins)
%     subplot(4, 1, 2)
%         hist(e2, nbins)
%     subplot(4, 1, 3)
%         hist(e3, nbins)
%     subplot(4, 1, 4)
%         hist(e4, nbins)

figure(1)   %plot spline curve
    plot(totaltravel, totalerror)
    title('Average of Spline Error Curves')
    xlabel('Scale Displacement, um'); ylabel('Error, um')
figure(2)   %plot histogram
    [freq, bins] = hist(eavg, nbins); bar(bins,freq)
    title('Error Curve Distribution')

%ball data
step_ball = (0 : step : period_scale); %step_ball(end) = [];%define step of ball, throw last point, will be the same as first
[~, nball] = size(step_ball);
[~, step_ball] = meshgrid(linspace(0,1,nmeas), step_ball');
ball = zeros(nball, nmeas); result_ball = zeros(nball,nmeas);
%generate normally distributed ball data individually
for i = 1 : nball   %need loops because interp1 expects vectors, no matrices
    ball(i,:) = step_ball(i,:) + sigma * randn(1, nmeas);    %generate normal ball data
    result_ball(i,:) = interp1(totaltravel, totalerror, ball(i,:), 'spline');   %calculate spline interpolation of ball positions on error curve
    devresult(i,:) = sum(result_ball(i,:) / nmeas);
end

figure(3)
    plot(step_ball(:,1), devresult, 'b-*'); hold on
    plot(xx, eavg, 'r-')
    hold off