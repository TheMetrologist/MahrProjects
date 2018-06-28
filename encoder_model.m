%authored by M. Braine
%scale and encoder parameters
clear all
period_scale = 4;                                                               %define period/pitch of scale, slit spacing, um
x_step = 0.05;                                                                  %how far to step scale and take point
travel = 16;                                                                    %distance to displace scale in um
disp(' ');

%ask for encoder errors. no input is all 0% error
error_encoder = input('Enter the errors for each sensor in array [e1 e2 e3 ...] format: ')/100;%ask for sensor errors
n_encoders = length(error_encoder);                                             %count the number of sensors

%encoder spacing and locations
spacing_encoder = period_scale / n_encoders;                                    %define encoder spacing, assume evenly spaced
x_encoder = zeros(1, n_encoders); x_encoder(1) = period_scale / (2*n_encoders); %define first sensor position on encoder
i = 2;
while i <= n_encoders                                                           %calculate positions of each encoder in um within one scale period
    x_encoder(i) = x_encoder(i-1) + spacing_encoder;
    i=i+1;
end %end while

%preallocate large matrices for iterative calculations
intensity = zeros((travel/period_scale)*(period_scale/x_step + 1), n_encoders); %generate matrix for intensity values
[n, trash] = size(intensity);
phase = linspace(0, 360*(travel/period_scale), n)';
[phase, x_encoder] = meshgrid(phase, x_encoder);                                %generate matrices for phase, x_encoder
[error_encoder, phase] = meshgrid(error_encoder, phase(1,:));
solution = zeros(324,3);
x_encoder = x_encoder';

%calculate intensity with simulated encoder error
intensity = sind(x_encoder * (360/period_scale) + phase) + error_encoder .* sind(x_encoder * (360/period_scale) + phase);
%intensity = sind(x_encoder * (360/period_scale) + phase).^2 + error_encoder .* sind(x_encoder * (360/period_scale) + phase).^2 - 0.5;

%sin fit
f = fittype('sin1');                                                            %define the type of fit: a1 * sin(b1 * x + c1)
%f = fittype(@(a1,b1,c1,x) a1*sin(b1*x+c1).^2);
    options = fitoptions('Method','NonLinearLeastSquares');%options for fit

for i = 1 : n                                                                   %fit data to defined fit type
    abc = fit(x_encoder(i,:)',intensity(i,:)',f,options);
    solution(i,1) = abc.a1;                                                     %store amplitude
    solution(i,2) = abc.b1;                                                     %store frequency
    solution(i,3) = abc.c1;                                                     %store phase shift
end
solution(:,3) = solution(:,3)'*180/pi;                                          %convert phase solution from radians to degrees

%fit solves for each phase of 0 thru 360. add additional 360 degrees to each pitch for linear line
index = sign(solution(:,3));                                                    %find sign of phase
for i = 2:n
    if index(i) < index(i-1)                                                    %if phase transition from 360 -> 0, correct
        solution(i:end,3) = solution(i:end,3) + 360;
%        solution_phase(i:end) = solution_phase(i:end)+90;
    end
end

%convert displacement and scale error from degrees to um for plotting
displacement = phase(:,1)/360 * period_scale;
error_scale = solution(:,3)/360 * period_scale;
figure(1)
    plot(displacement,error_scale - displacement,'-*')
    title('Simulated Error Curve'); xlabel('Scale Displacement, um'); ylabel('Error, um')

%Plot 0-deg phase shift intensities
figure(2); subplot(4,1,1); title('Encoder Position and Simulated Intensity, 0-deg. Phase')
    plot(x_encoder(1,:), intensity(1,:), 'o'); hold on
    xlabel('Encoder Position, um'); ylabel('Intensity')

    %plot fit sine curve
    xx = 0 : 0.01 : period_scale;   %x data for plotting fit sin(x)
    yy = solution(1,1) * sind(solution(1,2)*(xx*180/pi) + solution(1,3)); %y data for fit sin(x)
    plot(xx, yy, 'r-')  %plot fit curve
    axis([0 period_scale min(intensity(1,:))-0.5 max(intensity(1,:))+0.5])

%     xx = 0 : 0.01 : period_scale;   %data for fitting spline
%     yy = spline(x_encoder(1,:), intensity(1,:), xx);    %fit spline
%     plot(xx, yy, 'r-')  %plot spline
clear yy

%plot 90-deg phase shift
subplot(4,1,2); title('Encoder Position and Simulated Intensity, 90-deg. Phase')
    plot(x_encoder(21,:), intensity(21,:), 'o'); hold on
    xlabel('Encoder Position, um'); ylabel('Intensity')

    %plot fit sine curve
    yy = solution(21,1) * sind(solution(21,2)*(xx*180/pi) + solution(21,3)); %y data for fit sin(x)
    plot(xx, yy, 'r-')  %plot fit curve
    axis([0 period_scale min(intensity(21,:))-0.5 max(intensity(21,:))+0.5])
clear yy

%plot 180-deg phase shift
subplot(4,1,3); title('Encoder Position and Simulated Intensity, 180-deg. Phase')
    plot(x_encoder(42,:), intensity(42,:), 'o'); hold on
    xlabel('Encoder Position, um'); ylabel('Intensity')

    %plot fit sine curve
    yy = solution(42,1) * sind(solution(42,2)*(xx*180/pi) + solution(42,3)); %y data for fit sin(x)
    plot(xx, yy, 'r-')  %plot fit curve
    axis([0 period_scale min(intensity(42,:))-0.5 max(intensity(42,:))+0.5])
clear yy

%plot 270-deg phase shift
subplot(4,1,4); title('Encoder Position and Simulated Intensity, 270-deg. Phase')
    plot(x_encoder(62,:), intensity(62,:), 'o'); hold on
    xlabel('Encoder Position, um'); ylabel('Intensity')

    %plot fit sine curve
    yy = solution(62,1) * sind(solution(62,2)*(xx*180/pi) + solution(62,3)); %y data for fit sin(x)
    plot(xx, yy, 'r-')  %plot fit curve
    axis([0 period_scale min(intensity(62,:))-0.5 max(intensity(62,:))+0.5])

figure(3); title('Error Distribution')
    hist(error_scale - displacement, 15)

disp(' ')