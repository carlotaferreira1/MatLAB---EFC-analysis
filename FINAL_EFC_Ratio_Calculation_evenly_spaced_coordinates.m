% Load the XY boundary coordinates
data = readmatrix(['/Users/carlotaferreira/Documents/HEK coordinates max 4.csv']);
x = data(:, 1);  % X-coordinates
y = data(:, 2);  % Y-coordinates

% Ensure the shape is closed
x = [x; x(1)];
y = [y; y(1)];

% Step 1: Compute the distances between consecutive points
distances = sqrt(diff(x).^2 + diff(y).^2);  % Euclidean distance

% Step 2: Compute the cumulative arc length
arc_length = [0; cumsum(distances)];

% Step 3: Define the desired number of evenly spaced points
num_interpolated_points = length(x);  % Automatically set to the number of points in x

% Step 4: Create a new set of evenly spaced arc lengths
new_arc_length = linspace(0, arc_length(end), num_interpolated_points);

% Step 5: Interpolate the X and Y coordinates based on the new arc lengths
x_interp = interp1(arc_length, x, new_arc_length, 'linear');
y_interp = interp1(arc_length, y, new_arc_length, 'linear');

% Number of points
n_points = num_interpolated_points;  % Use interpolated points
t = linspace(0, 2*pi, n_points)';  % Uniformly parameterize the shape boundary

% Perform Discrete Fourier Transform (DFT) on the X and Y coordinates
cX = fft(x_interp) / n_points;  % Fourier coefficients for X (normalized)
cY = fft(y_interp) / n_points;  % Fourier coefficients for Y (normalized)

% Step 3: Extract the first 20 harmonics
num_harmonics = 20;  % Number of harmonics to use
harmonics = 1:num_harmonics;

% Initialize major and minor axes storage
major_axes = zeros(num_harmonics, 1);
minor_axes = zeros(num_harmonics, 1);

% Step 4: Calculate major and minor axes for each harmonic
for n = harmonics
    % Calculate the amplitude of each harmonic (magnitude of the Fourier coefficients)
    a_n = abs(cX(n + 1));  % X-component amplitude (cosine term)
    b_n = abs(cY(n + 1));  % Y-component amplitude (sine term)
    
    % Major and minor axes are proportional to these amplitudes
    major_axes(n) = a_n;  % Treat as "major axis"
    minor_axes(n) = b_n;  % Treat as "minor axis"
end

% Step 5: Calculate the EFC ratio
% EFC numerator: sum of major and minor axes for the first harmonic
EFC_numerator = major_axes(1) + minor_axes(1);

% EFC denominator: sum of major and minor axes for the remaining harmonics
EFC_denominator = sum(major_axes(2:end) + minor_axes(2:end));

% EFC ratio calculation
EFC_ratio = EFC_numerator / EFC_denominator;

% Display the results
disp(['EFC Ratio: ', num2str(EFC_ratio)]);

% Step 6: Reconstruct the shape using the first 'num_harmonics'
x_reconstructed = zeros(size(t));
y_reconstructed = zeros(size(t));

for n = harmonics
    % Add contributions of each harmonic to the reconstruction
    x_reconstructed = x_reconstructed + real(cX(n + 1) * exp(1i * n * t));
    y_reconstructed = y_reconstructed + real(cY(n + 1) * exp(1i * n * t));
end

% Step 7: Plot the original and reconstructed shapes
figure;

% Original shape
subplot(1, 2, 1);
plot(x, y, 'b-', 'LineWidth', 1.5);
title('Original Shape');
axis equal;
xlabel('X'); ylabel('Y');

% Reconstructed shape
subplot(1, 2, 2);
plot(x_reconstructed, y_reconstructed, 'r-', 'LineWidth', 1.5);
title(['Reconstructed Shape (', num2str(num_harmonics), ' Harmonics)']);
axis equal;
xlabel('X'); ylabel('Y');

% Display major and minor axes for debugging
%disp('Major Axes (Harmonics):');
%disp(major_axes);

%disp('Minor Axes (Harmonics):');
%disp(minor_axes);
