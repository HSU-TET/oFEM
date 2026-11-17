function pulse = pulseGen(endtime, timestep,peak,f,c)
t = 0:timestep:endtime;
x=0; 

pulse = peak*sin(2*pi*f(x/c -t));

end