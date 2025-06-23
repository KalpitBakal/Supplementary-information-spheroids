function [x,y]=projection(a,b,p_x,p_y)
    c_temp=(1+b^2)^0.5;
    n_perp.x=-b/c_temp;
    n_perp.y=1/c_temp;
    x=(p_y-a-n_perp.y/n_perp.x*p_x)/(b-n_perp.y/n_perp.x);
    y=a+b*x;
end