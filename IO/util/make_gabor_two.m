function G = make_gabor_two(width, sigma, scale, theta)

[x,y] = meshgrid(-width:width, -width:width);
%x = x*sigma;

u = (x*cos(theta)-y*sin(theta))*sigma*scale;
v = (x*sin(theta)+y*cos(theta))*scale;

d = u.^2+v.^2;

G = exp(-.5*d) .* cos(2*pi*u/width+pi/2);
G = G/sum(abs(G(:)));

end
