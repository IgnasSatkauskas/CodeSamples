function y = GaussLegendre10(f,a,b)

%10 pnt Gauss-Legendre quadrature to approximate the integral of f on [a,b]

%WEIGHTS AND NODES

%1	0.2955242247147529	-0.1488743389816312
%2	0.2955242247147529	0.1488743389816312
%3	0.2692667193099963	-0.4333953941292472
%4	0.2692667193099963	0.4333953941292472
%5	0.2190863625159820	-0.6794095682990244
%6	0.2190863625159820	0.6794095682990244
%7	0.1494513491505806	-0.8650633666889845
%8	0.1494513491505806	0.8650633666889845
%9	0.0666713443086881	-0.9739065285171717
%10	0.0666713443086881	0.9739065285171717

weights = [	0.2955242247147529;
            0.2955242247147529;	
            0.2692667193099963;	
            0.2692667193099963;	
            0.2190863625159820;	
            0.2190863625159820;	
            0.1494513491505806;	
            0.1494513491505806;
            0.0666713443086881;
            0.0666713443086881];

        
nodes =   [	-0.1488743389816312;
            0.1488743389816312;
            -0.4333953941292472;
            0.4333953941292472;
            -0.6794095682990244;
            0.6794095682990244;
            -0.8650633666889845;
            0.8650633666889845;
            -0.9739065285171717;
            0.9739065285171717];

 %x = .5*(b-a)*nodes + .5*(b+a);

 y = .5* (b-a) * weights.' * col( f( .5*(b-a)*nodes + .5*(b+a) ) );


end