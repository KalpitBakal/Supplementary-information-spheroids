
% Geometry file should look like
%A file with coordinates in the following form: Output from the ImageJ 
% PointPicker plugin
%
%point     x     y slice color
%    0   844   425     1     0
%    1   962   412     1     1
%    2   856   572     1     2
%    3   962   574     1     3
%    4   810   500     1     4
%    5  1010   489     1     5
%    6   962   574     1     6
%    0   760   436     3     0
%    1   893   421     3     1
%    2   766   566     3     2
%    3   904   570     3     3
%    4   729   502     3     4
%    5   942   490     3     5
%    6   962   574     1     6

% The pressure file should look like (pressure in mm H2O):
% Frame	Actual pressure
% 1	    0.0
% 2	    12.87
% 3	    38.62
% 4	    57.93
% 5	    83.67
% 6	    102.98
% 7	    122.29
% 8	    141.6
% 9	    160.9
% 10	186.65


clear
[FileName,PathName]=uigetfile('*.txt','Select the Pressure file (mm H2O)..');
cd ([PathName,'\']);
p_data_file=([PathName,'\',FileName]);
[FileName,PathName]=uigetfile('*.txt','Select the ball geometry data file..');
cd ([PathName,'\']);
ball_data_file=([PathName,'\',FileName]);

stress_balance_pc(ball_data_file,p_data_file,PathName,['stress_balance']);

