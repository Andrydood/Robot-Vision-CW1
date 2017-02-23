function [ updated_state_estimates,updated_error_covariance,output_latitude,output_longitude,output_height,output_velocity   ] = gnss_kalman( pseudo_ranges, pseudo_range_rates,state_estimate,error_covariance,epoch)

addpath('Libraries/');

deg_to_rad = 0.01745329252; % Degrees to radians conversion factor
rad_to_deg = 1/deg_to_rad; % Radians to degrees conversion factor
c = 299792458; % Speed of light in m/s
omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s
Omega_ie = Skew_symmetric([0,0,omega_ie]);
R_0 = 6378137; %WGS84 Equatorial radius in meters
e = 0.0818191908425; %WGS84 eccentricity


%Compute the transition matrix
tao = 0.5;

I = eye(3,3);

transition_matrix = [I, tao*I, zeros(3,1), zeros(3,1);
              zeros(3,3), I, zeros(3,1), zeros(3,1);
              zeros(1,3), zeros(1,3), 1, tao;
              zeros(1,3), zeros(1,3), 0, 1];

%Compute the system noise covariance matrix 
S_a = 10; %acceleration power spectral density          
S_c = 1; %clock phase   
S_f = 1; %clock frequency   

noise_covariance = [zeros(3,3), zeros(3,3), zeros(3,1), zeros(3,1);
             zeros(3,3), S_a*tao*I,  zeros(3,1), zeros(3,1);
             zeros(1,3), zeros(1,3), S_c*tao, 0;
              zeros(1,3), zeros(1,3), 0, S_f*tao];
          

%Propagate the state estimates
state_estimate_new = transition_matrix*state_estimate;

user_position = state_estimate_new(1:3);
user_velocity = state_estimate_new(4:6);
clock_offset_estimate = state_estimate_new(7);
clock_drift_estimate = state_estimate_new(8);

%Propagate the error covariance matrix
error_covariance_new = transition_matrix * error_covariance * transition_matrix' + noise_covariance;

%Compute the Cartesian ECEF positions of the satellites at time 0

satellites = pseudo_ranges(1,2:end);
sat_num = size(satellites,2);

sat_position = zeros(3,sat_num);
sat_velocity = zeros(3,sat_num);

sat_pseudo_range = pseudo_ranges(2+epoch,2:end);
sat_pseudo_range_rate = pseudo_range_rates(2+epoch,2:end);

time = pseudo_ranges(2+epoch,1);

for i = 1:sat_num
    
    [sat_position(:,i),sat_velocity(:,i)] = Satellite_position_and_velocity(time,satellites(i));

end

predicted_range = zeros(1, sat_num);
predicted_range_rate = zeros(1, sat_num);

%Predict the ranges and range rates from the approximate user position
%Set sagnac compensation as identity
for sat = 1: sat_num
    
    predicted_range(sat) = sqrt([eye(3)*sat_position(:,sat)-user_position]'*[eye(3)*sat_position(:,sat)-user_position]);

end

%Calculate sagnac compensation and calculate predicted range and predicted
%range rate
line_of_sight = zeros(3, sat_num);

for sat = 1: sat_num
    
    sagnac_compensation = (omega_ie/ c) .*[0 predicted_range(sat) 0; ...
                                            -predicted_range(sat) 0 0; ...
                                            0 0 0] +eye(3);
        
    predicted_range(sat) = sqrt([sagnac_compensation*sat_position(:,sat)-user_position]'*[sagnac_compensation*sat_position(:,sat)-user_position]);
    
    pos = (sagnac_compensation*sat_position(:,sat)-user_position)./   predicted_range(sat) ;
    
    predicted_range_rate(sat) = pos'*(sagnac_compensation*(sat_velocity(:,sat)+Omega_ie*sat_position(:,sat))-(user_velocity+Omega_ie*user_position));

    line_of_sight(:,sat) =  (sat_position(:,sat)-user_position)./   predicted_range(sat) ;
    
end

%Compute the measurement matrix 
measurements = zeros(sat_num*2,8);

for sat = 1:sat_num
     measurements(sat,:) = [-line_of_sight(1,sat), -line_of_sight(2,sat), -line_of_sight(3,sat), 0, 0, 0, 1, 0];
end

for sat = 1:sat_num
     measurements((sat+sat_num),:) = [0,0,0,-line_of_sight(1,sat), -line_of_sight(2,sat), -line_of_sight(3,sat), 0,1];
end

%Compute measurement noise covariance matrix
pseudo_sd = 5;
pseudo_rate_sd = 0.05;

measurement_noise = eye(sat_num*2,sat_num*2);

measurement_noise(1:sat_num,1:sat_num) = measurement_noise(1:sat_num,1:sat_num).*pseudo_sd^2;
measurement_noise(sat_num+1:sat_num*2,sat_num+1:sat_num*2) = measurement_noise(sat_num+1:sat_num*2,sat_num+1:sat_num*2).*pseudo_rate_sd^2;

%Compute kalman gain matrix

kalman_gain = (error_covariance_new*measurements')*inv(measurements*error_covariance_new*measurements'+measurement_noise);

%Compute measurement innovation vector

measurement_innovation = zeros(sat_num*2,1);

for sat = 1:sat_num

    measurement_innovation(sat) = [sat_pseudo_range(sat)-predicted_range(sat)-clock_offset_estimate];
    
end

for sat = 1:sat_num

    measurement_innovation(sat+sat_num) = [sat_pseudo_range_rate(sat)-predicted_range_rate(sat)-clock_drift_estimate];
    
end

%Update the state estimates and error covariance matrix

updated_state_estimates = state_estimate_new + kalman_gain*measurement_innovation;

updated_error_covariance = (eye(8,8)-kalman_gain*measurements)*error_covariance_new;

%Convert solution to longitude, latitude and height
[output_latitude,output_longitude,output_height,output_velocity] = pv_ECEF_to_NED(updated_state_estimates(1:3),updated_state_estimates(4:6));
output_latitude = output_latitude/0.01745329252;
output_longitude = output_longitude/0.01745329252;
end

