
Wp = 5.2;
Wm = 2.2;
Tp = 20;
Tm = 100;

 s1 = 20;
 s2 = 40;
 s3 = 60;
 s4 = 80;
 s5 = 90;
 s6 = 100;
 
 t1 = -100:0;
 t2 = 0:100;
 
 figure; hold;
 plot([-100 100], [0 0], 'k');
 plot([0 0], [-6 6], 'k');
 
 eVal = exp((s1-100)/25);
 wm = ((-Wp - Wm) * eVal) + Wm;
 tm = ((Tp - Tm) * eVal) + Tm;
 plot(t1, Wp * exp(t1/Tp) - eVal);
 plot(t2, -wm * exp(-t2/tm) - eVal);
 
 eVal = exp((s2-100)/25);
 wm = ((-Wp - Wm) * eVal) + Wm;
 tm = ((Tp - Tm) * eVal) + Tm;
 plot(t1, Wp * exp(t1/Tp) - eVal, 'c');
 plot(t2, -wm * exp(-t2/tm) - eVal, 'c');
 
 eVal = exp((s3-100)/25);
 wm = ((-Wp - Wm) * eVal) + Wm;
 tm = ((Tp - Tm) * eVal) + Tm;
 plot(t1, Wp * exp(t1/Tp) - eVal, 'g');
 plot(t2, -wm * exp(-t2/tm) - eVal, 'g');
 
 eVal = exp((s4-100)/25);
 wm = ((-Wp - Wm) * eVal) + Wm;
 tm = ((Tp - Tm) * eVal) + Tm;
 plot(t1, Wp * exp(t1/Tp) - eVal, 'y');
 plot(t2, -wm * exp(-t2/tm) - eVal, 'y');
 
 eVal = exp((s5-100)/25);
 wm = ((-Wp - Wm) * eVal) + Wm;
 tm = ((Tp - Tm) * eVal) + Tm;
 plot(t1, Wp * exp(t1/Tp) - eVal, 'r');
 plot(t2, -wm * exp(-t2/tm) - eVal, 'r');
 
  eVal = exp((s6-100)/25);
 wm = ((-Wp - Wm) * eVal) + Wm;
 tm = ((Tp - Tm) * eVal) + Tm;
 plot(t1, Wp * exp(t1/Tp) - eVal, 'm');
 plot(t2, -wm * exp(-t2/tm) - eVal, 'm');