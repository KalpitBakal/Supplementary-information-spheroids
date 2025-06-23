function [e_G_, p_wall_min_p_, V_strain_,p_av_, n_images,K_step_av,K_step_stdev,G_step_av,G_step_stdev,K_av,K_stdev,G_av,G_stdev,V_average, p_average]=stress_balance_pc(ball_data_file,p_data_file,output_path,output_name)


cd (output_path);

%ball_data_file: full path to ball data file
%pressure_file: full path to pressure data file
%output_path: path were output is written
%

%
% This code is for analyzing the geometry of soft particles deformed in
% glass capillaries by an applied pressure. ("Micromechanical squeezing tests")
% It converts the manually (by using ImageJ) determined positions to a
% variety of geometrical features of the particle at different pressures
% and times.
%
%
%  Version where we balance the work done by the pressure between two
%  measurement points with the elastic energy.
% Hans Wyss, Apr. 2009

%Edited by Kalpit Jagdish Bakal, August 2022




%INPUT:
%*******
%A file with coordinates in the following form:
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
%
%
%point: 6 points for each image
% convention: the images have the tip of the capillary pointing towards the
% left side of the image. With this convention the indices are as follows:
%           0: reference point (in case the field of view is changing
%           around.
%           1: upper left corner of the gel ball (triple boundary
%            particle/glass/water)
%           2: upper right corner
%           3: lower left corner
%           4: lower right corner
%           5: leftmost point of the particle
%           6: rightmost point of the particle
%
%x: x-coordinate of the point
%y: y-coordinate of the point
%slice: image# in a stack of images being analyzed. (this is the image# in
%alphabetical order inside a folder when images were imported from a folder using the
%command "import image series" in ImageJ.
%



%********** INPUTS **************
%********************************

%pixels_per_mu=3.2*1e6
%pixels_per_mu=1.5*1e6
%pixels_per_mu=3.2*1e6
pixels_per_mu=1;

K=1;

%********************************
%********************************

%Reading ball geometry data file:

ball_data = importdata(ball_data_file);
ball_data;
colheaders = genvarname(ball_data.colheaders);

point=NaN;
x=NaN;
y=NaN;
slice=NaN;
color=NaN;
point=ball_data.data(:,1);
x=ball_data.data(:,2);
y=ball_data.data(:,3);
slice=ball_data.data(:,4);
% color=ball_data.data(:,5);


p_data = importdata(p_data_file);
p_data;
colheaders = genvarname(p_data.colheaders);
frame=NaN;
height=NaN;
frame=p_data.data(:,1);
height=p_data.data(:,2);
% disp(height)



% main analysis:
%----------------
n_points=numel(x);
n_images=floor(n_points/7);


% Fit top line:
%****************
x_coords=0; y_coords=0;
i=1;
for n=1:n_images
    x0=x(7*n-6);% (x-coordinate of reference point)
    y0=y(7*n-6);
    x_coords(i)=x(7*n-5)-x0;y_coords(i)=y(7*n-5)-y0;%(point 1, according to manual tracking convention)
    i=i+1;
    x_coords(i)=x(7*n-4)-x0;y_coords(i)=y(7*n-4)-y0;%(point 2,  "  
    i=i+1;
end
p_top=polyfit(x_coords,y_coords,1);
a_top_from_fit=p_top(2);
b_top_from_fit=p_top(1);


% Fit bottom line:
%****************
x_coords=0; y_coords=0;
i=1;
for n=1:n_images
    x0=x(7*n-6);% (x-coordinate of reference point)
    y0=y(7*n-6);
    x_coords(i)=x(7*n-3)-x0;y_coords(i)=y(7*n-3)-y0;%(point 1, according to manual tracking convention)
    i=i+1;
    x_coords(i)=x(7*n-2)-x0;y_coords(i)=y(7*n-2)-y0;%(point 2,  "  
    i=i+1;
end
p_bot=polyfit(x_coords,y_coords,1);
a_bottom_from_fit=p_bot(2);
b_bottom_from_fit=p_bot(1);

alpha=abs(b_top_from_fit-b_bottom_from_fit)/2;


% For each image:
for n=1:n_images
    %get positions relative to reference point:
    x0=x(7*n-6);% (x-coordinate of reference point)
    y0=y(7*n-6);% ..
    p1.x=x(7*n-5)-x0;p1.y=y(7*n-5)-y0;%(point 1, according to manual tracking convention)
    p2.x=x(7*n-4)-x0;p2.y=y(7*n-4)-y0;%(point 2,  "                                     )
    p3.x=x(7*n-3)-x0;p3.y=y(7*n-3)-y0;% ..
    p4.x=x(7*n-2)-x0;p4.y=y(7*n-2)-y0;
    p5.x=x(7*n-1)-x0;p5.y=y(7*n-1)-y0;
    p6.x=x(7*n)-x0;p6.y=y(7*n)-y0;

    if n==1
        %top line (connection between points p1 and p2)
           a_top=a_top_from_fit;
           b_top=b_top_from_fit;
        %bottom line (connection between points p3 and p4)
           a_bott=a_bottom_from_fit;
           b_bott=b_bottom_from_fit;
        %middle line: line between top and bottom lines (through the center of
        %the capillary)
        b_mid=(b_top+b_bott)/2;
        a_mid=(a_top+a_bott)/2;
        %intersection between top and bottom lines:
        p_intsect.x=(a_bott-a_top)/(b_top-b_bott);
        p_intsect.y=a_top+p_intsect.x*b_top;
    end
    %project manual points onto top, bottom and middle lines:
    [p1.x,p1.y]=projection(a_top,b_top, p1.x,p1.y);
    [p2.x,p2.y]=projection(a_top,b_top, p2.x,p2.y);
    [p3.x,p3.y]=projection(a_bott,b_bott, p3.x,p3.y);
    [p4.x,p4.y]=projection(a_bott,b_bott, p4.x,p4.y);
    [p5.x,p5.y]=projection(a_mid,b_mid, p5.x,p5.y);
    [p6.x,p6.y]=projection(a_mid,b_mid, p6.x,p6.y);
    
    
    %normalized vectors parallel and perpendicular to the middle line;
    c_temp=(1+b_mid^2)^0.5;
    n_par.x=1/c_temp;
    n_par.y=b_mid/c_temp;
    n_perp.x=-b_mid/c_temp;
    n_perp.y=1/c_temp;
    %Projections onto middle line:
    [p0mid.x,p0mid.y]=projection(a_mid,b_mid,0,0);
    [p1mid.x,p1mid.y]=projection(a_mid,b_mid,p1.x,p1.y);
    [p2mid.x,p2mid.y]=projection(a_mid,b_mid,p2.x,p2.y);
    [p3mid.x,p3mid.y]=projection(a_mid,b_mid,p3.x,p3.y);
    [p4mid.x,p4mid.y]=projection(a_mid,b_mid,p4.x,p4.y);
    [p5mid.x,p5mid.y]=projection(a_mid,b_mid,p5.x,p5.y);
    [p6mid.x,p6mid.y]=projection(a_mid,b_mid,p6.x,p6.y);
    %Distance from tip:
    d1=dist_2D(p1mid.x,p1mid.y,p0mid.x,p0mid.y);
    d2=dist_2D(p2mid.x,p2mid.y,p0mid.x,p0mid.y);
    d3=dist_2D(p3mid.x,p3mid.y,p0mid.x,p0mid.y);
    d4=dist_2D(p4mid.x,p4mid.y,p0mid.x,p0mid.y);
    d5=dist_2D(p5mid.x,p5mid.y,p0mid.x,p0mid.y);
    d6=dist_2D(p6mid.x,p6mid.y,p0mid.x,p0mid.y);
    %Contact band (the contact area between the particle and the glass)
    %******************************************************************
    d_band_front=(d1+d3)/2; % distance from the front of the band to the tip (measured along the middle line)
    d_band_back=(d2+d4)/2;  % distance from the front of the band to the tip (measured along the middle line)
    L_band=d_band_back-d_band_front; %length of the band (measured along the middle line)
    W_band=(dist_2D(p1.x,p1.y,p2.x,p2.y)+dist_2D(p3.x,p3.y,p4.x,p4.y))/2; % length of the band (measured along the side)
    R_band_front=(dist_2D(p1.x,p1.y,p1mid.x,p1mid.y)+dist_2D(p3.x,p3.y,p3mid.x,p3mid.y))/2; %radius of the circle described by the front of the band.
    R_band_back=(dist_2D(p2.x,p2.y,p2mid.x,p2mid.y)+dist_2D(p4.x,p4.y,p4mid.x,p4mid.y))/2;  %radius of the circle described by the back of the band.
    c_temp=(R_band_back-R_band_front)/W_band;
    A_band=2*pi*(R_band_front*W_band+c_temp/2*W_band^2);% <=== band area: total contact area between the particle and the wall
    Length=dist_2D(p5mid.x,p5mid.y,p6mid.x,p6mid.y); % length of the particle measured along the middle line (distance between projections of points 5 and 6 onto middle line)
    
    %Conversion from applied pressure to wall pressure:
    p=height(n)*10;
    p_(n)=p;
    F_p=p*(pi*R_band_back^2);
    F_wall=F_p/sin(alpha);
    p_wall=F_wall/A_band;
    p_wall_(n)=p_wall;
    p_wall_min_p=p_wall-p;
    A_ratio=pi*R_band_back^2/A_band;
    p_wall_min_p_(n)=p_wall_min_p/2;
    p_av_(n)=(2*p_wall_(n)+p_(n))/3;
    
    %Volumes and Surfaces:
    h_cap_front=(d1+d3)/2-d5;
    h_cap_back=d6-(d2+d4)/2;
    a_cap_front=(dist_2D(p5mid.x,p5mid.y,p1.x,p1.y)+dist_2D(p5mid.x,p5mid.y,p3.x,p3.y))/2;
    a_cap_back=(dist_2D(p6mid.x,p6mid.y,p2.x,p2.y)+dist_2D(p6mid.x,p6mid.y,p4.x,p4.y))/2;
    V_cap_front=pi*h_cap_front*(3*a_cap_front^2+h_cap_front^2)/6;
    V_cap_back=pi*h_cap_back*(3*a_cap_back^2+h_cap_back^2)/6;
    L_tot=L_band/(1-(R_band_front/R_band_back));   %Total length along the middle line from the back of the band to the intersection of the top and the middle lines.
    V_band=pi/3*(R_band_back^2*L_tot-R_band_front^2*(L_tot-L_band));
    V=V_band+V_cap_back+V_cap_front;
    V_(n)=V;
    log_V_(n)=log(V_(n));
    minlog_V_(n)=-log(V_(n));
    log_R_band_(n)=log((R_band_front+R_band_back)/2);
    minlog_R_band_(n)=-log_R_band_(n);
    log_L_band_(n)=log(L_band);
    log_Length_(n)=log(Length);
    log_L_tot_(n)=log(L_tot);
    R_band=(R_band_front+R_band_back)/2;
    A=pi*R_band^2;
    
    R_band_(n)=(R_band_front+R_band_back)/2;
    Length_(n)=Length;
    
    %strains:
    if n==1
        R_0=(3*V/4/pi)^(1/3);
        R_band_strain=0;
        R_band0=R_band;
        A_0=pi*R_0^2;
        A_strain=(A_0-A)/A_0;  %"Area strain": relative change in cross-section area.
                               % Comparison point: cross-section area of
                               % sphere corresponding to the volume
                               % measured at the lowest pressure.
        minlog_A_rel=-log(A/A_0);
        minlog_R_rel=-log(R_band/R_0);
        
        L_band_strain=0;
        L_band0=L_band;
        Length_strain=0;
        Length0=Length;
        V_strain=0;
        V0=V;
        % epsilon_1=0; %strain due to pure compression, with compressive modulus K. Set to zero at the lowest pressure
        % V_strain1=0; %Volumetric strain due to pure compression, with compressive modulus K
        % V_strain2=0; %Volumetric strain after subtracting strain due to pure compression
        % Length_strain2=0; %Length strain after subtracting strain due to pure compression
        % R_band_strain2=0; %Radial strain after subtracting strain due to pure compression
        e_r=0;
        e_z=0;
    else
        
       
        R_band_strain=(R_band0-R_band)/R_band0;
        A_strain=(A_0-A)/A_0;  %"Area strain": relative change in cross-section area.
                               % Comparison point: cross-section area of
                               % sphere corresponding to the volume
                               % measured at the lowest pressure.
        minlog_A_rel=-log(A/A_0);
        minlog_R_rel=-log(R_band/R_0);
        
        L_band_strain=(L_band0-L_band)/L_band0;
        Length_strain=(Length0-Length)/Length0;
        V_strain=abs((V0-V)/V0);
        epsilon_1=1-exp(-(p-p_(1))/(6*K));
        V_strain1=1-exp(-(p-p_(1))/(2*K));
        V_strain2=(V_strain-V_strain1)/(1-V_strain1);
        Length_strain2=(Length_strain-epsilon_1)/(1-epsilon_1);
        %R_band_strain2=(R_band_strain-epsilon_1)/(1-epsilon_1);
        
        e_r=R_band_strain;
        e_z=L_band_strain; 
        %e_G_(n)=2*(e_r - e_z);
        e_G_(n)=abs((e_r-e_z));
        %By balancing the elastic stress, we get:
        K_(n)= (2*p_wall_(n) - p_(n))/3 / (2*e_r + e_z);
        G_(n)= (p_wall_(n) - p_(n))/2/(e_r - e_z);
        % AGain balancing the elastic stress, but evaluating each pressure step separately, we get:
        delta_p_wall_(n)=p_wall_(n)-p_wall_(1);
        delta_p_(n)=p_(n)-p_(1);
        G_step_(n)= (delta_p_wall_(n) - delta_p_(n))/2/(e_r - e_z );
        K_step_(n)= (2* delta_p_wall_(n) + delta_p_(n))/3/V_strain;
        
    end
    
    L_band_strain_(n)=L_band_strain;
    R_band_strain_(n)=R_band_strain;
    
    minlog_A_rel_(n)=minlog_A_rel;
    minlog_R_rel_(n)=minlog_R_rel;
    A_(n)=A;
    A_strain_(n)=A_strain;
   
    
    alpha_(n)=alpha;
    % Length_strain_(n)=Length_strain;
    V_strain_(n)=V_strain;
    R_shape=R_band*(V0/V)^(1/3);
    % R_shape_(n)=R_shape;
    minlog_R_shape=-log(R_shape);


    

    %Output:
    if n==1
        fid = fopen(([output_path,output_name,'.txt']), 'wt');
        outputstring=['n',     '\t','L_band','\t','R_band','\t',         'R_band_strain','\t',       'L_band_strain','\t','p [Pa]','\t',        'p_wall [Pa]','\t', '(2 p_wall + p) /3','\t',          '(p_wall - p) / 2',    '\t',      'V_strain',    '\t', '(e_r - e_z)',   '\t',       '(V_strain - e_r)',            '\t',   'V[pix^3]','\t',' -log(V)',        '\t',' -log(R_band)',     '\t',     ' -log(L_band)','\t',      ' -log(Length)',    '\t',   ' -log(L_tot)','\t',      'p_wall-p',          '\t','-','\t',                     'K_step','\t',                     'G_step'];
        fprintf(fid,outputstring);
        outputstring=['\n',num2str(n),'\t',num2str(L_band),'\t',num2str(R_band),'\t',num2str(R_band_strain),'\t',num2str(L_band_strain),'\t',num2str(p),'\t',num2str(p_wall),'\t',num2str(((2*p_wall+p)/3)),'\t',  num2str((p_wall-p)/2),     '\t', num2str(V_strain),    '\t', num2str(e_r-e_z),   '\t', num2str(V_strain-e_r),       '\t',  num2str(V),'\t',num2str(-log(V)), '\t',num2str(-log(R_band)), '\t',num2str(-log(L_band)), '\t',num2str(-log(Length)), '\t',num2str(-log(L_tot)), '\t',num2str(p_wall_min_p), '\t','-','\t',          '','\t',          ''];
        fprintf(fid,outputstring);
    else
        outputstring=['\n',num2str(n),'\t',num2str(L_band),'\t',num2str(R_band),'\t',num2str(R_band_strain),'\t',num2str(L_band_strain),'\t',num2str(p),'\t',num2str(p_wall),'\t',num2str(((2*p_wall+p)/3)),'\t',  num2str((p_wall-p)/2),     '\t', num2str(V_strain),    '\t', num2str(e_r-e_z),   '\t', num2str(V_strain-e_r),       '\t',  num2str(V),'\t',num2str(-log(V)), '\t',num2str(-log(R_band)), '\t',num2str(-log(L_band)), '\t',num2str(-log(Length)), '\t',num2str(-log(L_tot)), '\t',num2str(p_wall_min_p), '\t','-','\t',          num2str(K_step_(n)),'\t',          num2str(G_step_(n))];
        fprintf(fid,outputstring);
    end

    
end
 
   K_step_av=mean(K_step_);
   K_step_stdev=std(K_step_);
   G_step_av=mean(G_step_);
   G_step_stdev=std(G_step_);
   
   p=polyfit((V_strain_(1:n_images)),p_av_(1:n_images),1);
   
   K_av=p(1);
   K_stdev=NaN;

   p=polyfit((e_G_(1:n_images)),p_wall_min_p_(1:n_images),1);

   G_av=p(1);
   G_stdev=NaN;
   V_average=mean(V_);
   p_average=mean(p_);
    
%Figure 1: p_wall_min_p_ vs. 2(eps_r-eps_z)
%******************************************
%(a) Fits:
p=polyfit((e_G_(1:n_images)),p_wall_min_p_(1:n_images),1);
a1=p(2);
b1=p(1);
%(b) Figure:
figure(1);
hold on;
 clf reset;
plot((e_G_(1:n_images)),p_wall_min_p_(1:n_images),'bo',(e_G_(1:n_images)),a1+b1*(e_G_(1:n_images)),'b--','MarkerSize',10,'LineWidth',1.5);


axis square; axis on;
h=legend(' ',['slope: ',num2str(round(b1)),' Pa']);
h.NumColumns=2;
set(h,'Box','off');
set(gca,'Box','on');
set(gca,'FontSize',14);
set(gca,'LineWidth',2);
set(gca,'FontName','Times');
%xlabel('$$dr/r - dL/L$$','fontsize',24, 'interpreter','latex');
%ylabel('$$p_\mathrm{wall} - p$$ [Pa]','fontsize',24, 'interpreter','latex');



curr_fig=figure(1);

Image = getframe(curr_fig);
saveas(h,[output_path,'\',output_name,'_G_plot.fig'],'fig'); 
saveas(h,[output_path,'\',output_name,'_G_plot.pdf'],'pdf'); 

%print -depsc p_wall_min_p__vs_eps_r_min_eps_z;
%print -dtiff -r200 p_vs_minlog_A_rel.tiff
%********************************
hold off;


%Figure 2: p_av vs. V_strain_
%******************************************
%(a) Fits:
p=polyfit((V_strain_(1:n_images)),p_av_(1:n_images),1);
a1=p(2);
b1=p(1);
%(b) Figure:
figure(2);
hold on;
 clf reset;
plot((V_strain_(1:n_images)),p_av_(1:n_images),'bo',(V_strain_(1:n_images)),a1+b1*(V_strain_(1:n_images)),'b--','MarkerSize',10,'LineWidth',1.5);


axis square; axis on;
h=legend(' ',['slope: ',num2str(round(b1)),' Pa']);
h.NumColumns=2;
set(h,'Box','off');
set(gca,'Box','on');
set(gca,'FontSize',14);
set(gca,'LineWidth',2);
set(gca,'FontName','Times');
%xlabel('$$dr/r - dL/L$$','fontsize',24, 'interpreter','latex');
%ylabel('$$p_\mathrm{wall} - p$$ [Pa]','fontsize',24, 'interpreter','latex');



curr_fig=figure(2);

Image = getframe(curr_fig);
saveas(h,[output_path,'\',output_name,'_K_plot.fig'],'fig'); 
saveas(h,[output_path,'\',output_name,'_K_plot.pdf'],'pdf'); 
%print -depsc p_wall_min_p__vs_eps_r_min_eps_z;
%print -dtiff -r200 p_vs_minlog_A_rel.tiff
%********************************
hold off;


end

