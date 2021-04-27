function x = pulse(t) 

f = 5.0E+0; %[Hz]
w = 2*pi*f; 
peak = 5e2; %[V]
pulse = @(t) peak*sin(w*t).^2;
x = pulse(t);
end